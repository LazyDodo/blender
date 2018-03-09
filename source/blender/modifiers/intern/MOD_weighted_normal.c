/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_weighted_normal.c
 *  \ingroup modifiers
 */

#include "limits.h"

#include "MEM_guardedalloc.h"

#include "DNA_mesh_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_deform.h"
#include "BKE_mesh.h"

#include "BLI_math.h"
#include "BLI_stack.h"

#include "bmesh_class.h"

#include "MOD_modifiertypes.h"
#include "MOD_util.h"

#define INDEX_UNSET INT_MIN
#define INDEX_INVALID -1
#define IS_EDGE_SHARP(_e2l) (ELEM((_e2l)[1], INDEX_UNSET, INDEX_INVALID))

#define CLNORS_VALID_VEC_LEN (1e-4f)

typedef struct ModePair {
	float val;  /* Contains mode based value (face area / corner angle). */
	int index;  /* Index value per poly or per loop. */
} ModePair;

/* Sorting function used in modifier, sorts in decreasing order. */
static int modepair_cmp_by_val_inverse(const void *p1, const void *p2)
{
	ModePair *r1 = (ModePair *)p1;
	ModePair *r2 = (ModePair *)p2;

	if (r1->val < r2->val)
		return 1;
	else if (r1->val > r2->val)
		return -1;

	return 0;
}

/* Sorts by index in increasing order. */
static int modepair_cmp_by_index(const void *p1, const void *p2)
{
	ModePair *r1 = (ModePair *)p1;
	ModePair *r2 = (ModePair *)p2;

	if (r1->index > r2->index)
		return 1;
	else if (r1->index < r2->index)
		return -1;

	return 0;
}

#define NUM_CACHED_INVERSE_POWERS_OF_WEIGHT 128

typedef struct WeightedNormalData {
	const int numVerts;
	const int numEdges;
	const int numLoops;
	const int numPolys;

	MVert *mvert;
	MEdge *medge;

	MLoop *mloop;
	short (*clnors)[2];

	MPoly *mpoly;
	float (*polynors)[3];
	int *poly_strength;

	MDeformVert *dvert;
	const int defgrp_index;
	const bool use_invert_vgroup;

	const float weight;
	short mode;

	/* Lower-level, internal processing data. */
	float cached_inverse_powers_of_weight[NUM_CACHED_INVERSE_POWERS_OF_WEIGHT];

	ModePair *mode_pair;

	int *loop_to_poly;

	float (*vert_normals)[3];

	int *vert_loops_count;  /* Count number of loops using this vertex so far. */
	float *curr_vert_val;  /* Current max val for this vertex. */
	int *curr_vert_strength;  /* Current max strength encountered for this vertex. */
} WeightedNormalData;


static bool check_strength(WeightedNormalData *wn_data, const int mv_index, const int mp_index)
{
	BLI_assert (wn_data->poly_strength != NULL);
	BLI_assert (wn_data->curr_vert_strength != NULL);

	const int mp_strength = wn_data->poly_strength[mp_index];

	int *curr_vert_strength = wn_data->curr_vert_strength;
	float *curr_vert_val = wn_data->curr_vert_val;
	int *vert_loops_count = wn_data->vert_loops_count;
	float (*vert_normals)[3] = wn_data->vert_normals;

	if ((mp_strength == FACE_STRENGTH_STRONG && curr_vert_strength[mv_index] != FACE_STRENGTH_STRONG) ||
	    (mp_strength == FACE_STRENGTH_MEDIUM && curr_vert_strength[mv_index] == FACE_STRENGTH_WEAK))
	{
		curr_vert_strength[mv_index] = mp_strength;
		curr_vert_val[mv_index] = 0.0f;
		vert_loops_count[mv_index] = 0;
		zero_v3(vert_normals[mv_index]);
	}

	return mp_strength == curr_vert_strength[mv_index];
}

static void apply_weights_sharp_loops(
        WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data,
        int *loop_index, int size, float(*loop_normal)[3])
{
	ModePair *mode_pair = wn_data->mode_pair;

	int *loop_to_poly = wn_data->loop_to_poly;
	float(*polynors)[3] = wn_data->polynors;
	int *poly_strength = wn_data->poly_strength;
	const int weight = wn_data->weight;

	for (int i = 0; i < size - 1; i++) {
		for (int j = 0; j < size - i - 1; j++) {
			if (wnmd->mode == MOD_WEIGHTEDNORMAL_MODE_FACE
				&& mode_pair[loop_to_poly[loop_index[j]]].val < mode_pair[loop_to_poly[loop_index[j + 1]]].val) {
				int temp = loop_index[j];
				loop_index[j] = loop_index[j + 1];
				loop_index[j + 1] = temp;
			}
			else if ((wnmd->mode == MOD_WEIGHTEDNORMAL_MODE_ANGLE || wnmd->mode == MOD_WEIGHTEDNORMAL_MODE_FACE_ANGLE)
				&& mode_pair[loop_index[j]].val < mode_pair[loop_index[j + 1]].val) {
				int temp = loop_index[j];
				loop_index[j] = loop_index[j + 1];
				loop_index[j + 1] = temp;
			}
		}
	}

	zero_v3(wn_data->vert_normals[0]);
	wn_data->vert_loops_count[0] = 0;
	wn_data->curr_vert_val[0] = 0.0f;
	wn_data->curr_vert_strength[0] = FACE_STRENGTH_WEAK;
	const bool has_face_influence = (wnmd->flag & MOD_WEIGHTEDNORMAL_FACE_INFLUENCE) != 0 && poly_strength != NULL;

	float *vert_normals = wn_data->vert_normals[0];
	int *vert_loops_count = &wn_data->vert_loops_count[0];
	float *curr_vert_val = &wn_data->curr_vert_val[0];

	for (int i = 0; i < size; i++) {
		int j, mp_index;
		bool do_loop = true;

		if (wnmd->mode == MOD_WEIGHTEDNORMAL_MODE_FACE) {
			j = loop_to_poly[loop_index[i]];
			mp_index = mode_pair[j].index;
		}
		else {
			j = loop_index[i];
			mp_index = loop_to_poly[j];
		}

		if (has_face_influence && poly_strength) {
			do_loop = check_strength(wn_data, 0, mp_index);
		}
		if (do_loop) {
			const float curr_val = mode_pair[j].val;

			float *cached_inverse_powers_of_weight = wn_data->cached_inverse_powers_of_weight;

			if (*curr_vert_val == 0.0f) {
				*curr_vert_val = curr_val;
			}
			if (!compare_ff(*curr_vert_val, mode_pair[j].val, wnmd->thresh)) {
				(*vert_loops_count)++;
				*curr_vert_val = curr_val;
			}

			/* Exponentially divided weight for each normal (since a few values will be used by most vertices, we cache those). */
			const int vl_count = *vert_loops_count;
			if (vl_count < NUM_CACHED_INVERSE_POWERS_OF_WEIGHT && cached_inverse_powers_of_weight[vl_count] == 0.0f) {
				cached_inverse_powers_of_weight[vl_count] = 1.0f / powf(weight, vl_count);
			}
			const float inverted_n_weight = vl_count < NUM_CACHED_INVERSE_POWERS_OF_WEIGHT ?
			                                    cached_inverse_powers_of_weight[vl_count] : 1.0f / powf(weight, vl_count);

			madd_v3_v3fl(vert_normals, polynors[mp_index], curr_val * inverted_n_weight);
		}
	}
	if (normalize_v3(vert_normals) < CLNORS_VALID_VEC_LEN) {
		zero_v3(vert_normals);
	}

	for (int i = 0; i < size; i++) {
		copy_v3_v3(loop_normal[loop_index[i]], vert_normals);
	}
}

/* Modified version of loop_split_worker_do which sets custom_normals without considering smoothness of faces or
 * loop normal space array.
 * Used only to work on sharp edges. */
static void loop_split_worker(
        WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data,
        MLoop *ml_curr, MLoop *ml_prev,
        int ml_curr_index, int ml_prev_index, int *e2l_prev, int mp_index,
        float (*loop_normals)[3], int (*edge_to_loops)[2])
{
	MEdge *medge = wn_data->medge;
	MLoop *mloop = wn_data->mloop;
	MPoly *mpoly = wn_data->mpoly;

	int *loop_to_poly = wn_data->loop_to_poly;
	float (*polynors)[3] = wn_data->polynors;

	if (e2l_prev) {
		int *e2lfan_curr = e2l_prev;
		const MLoop *mlfan_curr = ml_prev;
		int mlfan_curr_index = ml_prev_index;
		int mlfan_vert_index = ml_curr_index;
		int mpfan_curr_index = mp_index;

		BLI_Stack *loop_index = BLI_stack_new(sizeof(int), __func__);

		while (true) {
			const unsigned int mv_pivot_index = ml_curr->v;
			const MEdge *me_curr = &medge[mlfan_curr->e];
			const MEdge *me_org = &medge[ml_curr->e];

			BLI_stack_push(loop_index, &mlfan_vert_index);

			if (IS_EDGE_SHARP(e2lfan_curr) || (me_curr == me_org)) {
				break;
			}

			BKE_mesh_loop_manifold_fan_around_vert_next(
			            mloop, mpoly, loop_to_poly, e2lfan_curr, mv_pivot_index,
			            &mlfan_curr, &mlfan_curr_index, &mlfan_vert_index, &mpfan_curr_index);

			e2lfan_curr = edge_to_loops[mlfan_curr->e];
		}

		int *index = MEM_malloc_arrayN((size_t)BLI_stack_count(loop_index), sizeof(*index), __func__);
		int cur = 0;
		while (!BLI_stack_is_empty(loop_index)) {
			BLI_stack_pop(loop_index, &index[cur]);
			cur++;
		}
		apply_weights_sharp_loops(wnmd, wn_data, index, cur, loop_normals);
		MEM_freeN(index);
		BLI_stack_free(loop_index);
	}
	else {
		copy_v3_v3(loop_normals[ml_curr_index], polynors[loop_to_poly[ml_curr_index]]);
	}
}

static void aggregate_vertex_normal(
        WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data,
        const int mv_index, const int mp_index,
        const float curr_val,
        const bool use_face_influence)
{
	float (*polynors)[3] = wn_data->polynors;

	MDeformVert *dvert = wn_data->dvert;
	const int defgrp_index = wn_data->defgrp_index;
	const bool use_invert_vgroup = wn_data->use_invert_vgroup;

	const float weight = wn_data->weight;

	float (*vert_normals)[3] = wn_data->vert_normals;
	int *vert_loops_count = wn_data->vert_loops_count;
	float *curr_vert_val = wn_data->curr_vert_val;

	float *cached_inverse_powers_of_weight = wn_data->cached_inverse_powers_of_weight;

	const bool has_vgroup = dvert != NULL;
	const bool vert_of_group = has_vgroup && defvert_find_index(&dvert[mv_index], defgrp_index) != NULL;

	if (has_vgroup && ((vert_of_group && use_invert_vgroup) || (!vert_of_group && !use_invert_vgroup))) {
		return;
	}

	if (use_face_influence && !check_strength(wn_data, mv_index, mp_index)) {
		return;
	}

	/* If cur_val is 0 init it to present value. */
	if (curr_vert_val[mv_index] == 0.0f) {
		curr_vert_val[mv_index] = curr_val;
	}
	if (!compare_ff(curr_vert_val[mv_index], curr_val, wnmd->thresh)) {
		/* curr_vert_val and present value differ more than threshold, update. */
		vert_loops_count[mv_index]++;
		curr_vert_val[mv_index] = curr_val;
	}

	/* Exponentially divided weight for each normal (since a few values will be used by most vertices, we cache those). */
	const int vl_count = vert_loops_count[mv_index];
	if (vl_count < NUM_CACHED_INVERSE_POWERS_OF_WEIGHT && cached_inverse_powers_of_weight[vl_count] == 0.0f) {
		cached_inverse_powers_of_weight[vl_count] = 1.0f / powf(weight, vl_count);
	}
	const float inverted_n_weight = vl_count < NUM_CACHED_INVERSE_POWERS_OF_WEIGHT ?
	                                    cached_inverse_powers_of_weight[vl_count] : 1.0f / powf(weight, vl_count);

	madd_v3_v3fl(vert_normals[mv_index], polynors[mp_index], curr_val * inverted_n_weight);
}

static void apply_weights_vertex_normal(WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data)
{
	const int numVerts = wn_data->numVerts;
	const int numEdges = wn_data->numEdges;
	const int numLoops = wn_data->numLoops;
	const int numPolys = wn_data->numPolys;

	MVert *mvert = wn_data->mvert;
	MEdge *medge = wn_data->medge;

	MLoop *mloop = wn_data->mloop;
	short (*clnors)[2] = wn_data->clnors;
	int *loop_to_poly = wn_data->loop_to_poly;

	MPoly *mpoly = wn_data->mpoly;
	float (*polynors)[3] = wn_data->polynors;
	int *poly_strength = wn_data->poly_strength;

	MDeformVert *dvert = wn_data->dvert;

	const short mode = wn_data->mode;
	ModePair *mode_pair = wn_data->mode_pair;

	float (*vert_normals)[3] = MEM_calloc_arrayN((size_t)numVerts, sizeof(*vert_normals), __func__);
	int *vert_loops_count = MEM_calloc_arrayN((size_t)numVerts, sizeof(*vert_loops_count), __func__);
	float *curr_vert_val = MEM_calloc_arrayN((size_t)numVerts, sizeof(*curr_vert_val), __func__);
	int *curr_vert_strength = NULL;

	const bool keep_sharp = (wnmd->flag & MOD_WEIGHTEDNORMAL_KEEP_SHARP) != 0;
	const bool use_face_influence = (wnmd->flag & MOD_WEIGHTEDNORMAL_FACE_INFLUENCE) != 0 && poly_strength != NULL;
	const bool has_vgroup = dvert != NULL;

	if (use_face_influence) {
		curr_vert_strength = MEM_malloc_arrayN((size_t)numVerts, sizeof(*curr_vert_strength), __func__);
		for (int i = 0; i < numVerts; i++) {
			curr_vert_strength[i] = FACE_STRENGTH_WEAK;
		}
	}

	wn_data->vert_normals = vert_normals;
	wn_data->vert_loops_count = vert_loops_count;
	wn_data->curr_vert_val = curr_vert_val;
	wn_data->curr_vert_strength = curr_vert_strength;

	switch (mode) {
		case MOD_WEIGHTEDNORMAL_MODE_FACE:
			for (int i = 0; i < numPolys; i++) {
				const int mp_index = mode_pair[i].index;
				const float mp_val = mode_pair[i].val;

				int ml_index = mpoly[mp_index].loopstart;
				const int ml_index_end = ml_index + mpoly[mp_index].totloop;
				for (; ml_index < ml_index_end; ml_index++) {
					const int mv_index = mloop[ml_index].v;

					aggregate_vertex_normal(wnmd, wn_data, mv_index, mp_index, mp_val, use_face_influence);
				}
			}
			break;
		case MOD_WEIGHTEDNORMAL_MODE_ANGLE:
		case MOD_WEIGHTEDNORMAL_MODE_FACE_ANGLE:
			BLI_assert(loop_to_poly != NULL);

			for (int i = 0; i < numLoops; i++) {
				const int ml_index = mode_pair[i].index;
				const float ml_val = mode_pair[i].val;

				const int mp_index = loop_to_poly[ml_index];
				const int mv_index = mloop[ml_index].v;

				aggregate_vertex_normal(wnmd, wn_data, mv_index, mp_index, ml_val, use_face_influence);
			}
			break;
		default:
			BLI_assert(0);
	}

	for (int mv_index = 0; mv_index < numVerts; mv_index++) {
		if (normalize_v3(vert_normals[mv_index]) < CLNORS_VALID_VEC_LEN) {
			zero_v3(vert_normals[mv_index]);
		}
	}

	if (!keep_sharp && !has_vgroup) {
		BKE_mesh_normals_loop_custom_from_vertices_set(mvert, vert_normals, numVerts, medge, numEdges,
		                                               mloop, numLoops, mpoly, polynors, numPolys, clnors);
	}
	else {
		float (*loop_normal)[3] = MEM_calloc_arrayN((size_t)numLoops, sizeof(*loop_normal), "__func__");
		int *loop_to_poly_mem = NULL;

		/* We need loop to poly mapping at this stage, but conviniently BKE_mesh_normals_loop_split
		 * will generate it for us if we don't have it yet. */
		if (loop_to_poly == NULL) {
			loop_to_poly_mem = MEM_malloc_arrayN((size_t)numLoops, sizeof(*loop_to_poly_mem), __func__);
			loop_to_poly = loop_to_poly_mem;
		}

		BKE_mesh_normals_loop_split(mvert, numVerts, medge, numEdges, mloop, loop_normal, numLoops, mpoly, polynors,
		                            numPolys, true, (float)M_PI, NULL, clnors, loop_to_poly);

		for (int mp_index = 0; mp_index < numPolys; mp_index++) {
			const int ml_index = mpoly[mp_index].loopstart;
			const int ml_index_end = ml_index + mpoly[mp_index].totloop;

			for (int i = ml_index; i < ml_index_end; i++) {
				const int mv_index = mloop[i].v;
				if (!is_zero_v3(vert_normals[mv_index])) {
					copy_v3_v3(loop_normal[i], vert_normals[mv_index]);
				}
			}
		}

		if (keep_sharp) {
			int (*edge_to_loops)[2] = MEM_calloc_arrayN((size_t)numEdges, sizeof(*edge_to_loops), __func__);

			if (wnmd->mode == MOD_WEIGHTEDNORMAL_MODE_FACE) {
				qsort(mode_pair, numPolys, sizeof(*mode_pair), modepair_cmp_by_index);
			}
			else {
				qsort(mode_pair, numLoops, sizeof(*mode_pair), modepair_cmp_by_index);
			}
			MPoly *mp;
			int mp_index;
			for (mp = mpoly, mp_index = 0; mp_index < numPolys; mp++, mp_index++) {
				int ml_curr_index = mp->loopstart;
				const int ml_last_index = (ml_curr_index + mp->totloop) - 1;

				MLoop *ml_curr = &mloop[ml_curr_index];

				for (; ml_curr_index <= ml_last_index; ml_curr++, ml_curr_index++) {
					int *e2l = edge_to_loops[ml_curr->e];

					if ((e2l[0] | e2l[1]) == 0) {
						e2l[0] = ml_curr_index;
						/* Not considering smoothness of faces, UNSET if first loop encountered on this edge. */
						e2l[1] = INDEX_UNSET;
					}
					else if (e2l[1] == INDEX_UNSET) {
						if ((medge[ml_curr->e].flag & ME_SHARP) || ml_curr->v == mloop[e2l[0]].v) {
							e2l[1] = INDEX_INVALID;
						}
						else {
							e2l[1] = ml_curr_index;
						}
					}
					else if (!IS_EDGE_SHARP(e2l)) {
						e2l[1] = INDEX_INVALID;
					}
				}
			}

			for (mp = mpoly, mp_index = 0; mp_index < numPolys; mp++, mp_index++) {
				const int ml_last_index = (mp->loopstart + mp->totloop) - 1;
				int ml_curr_index = mp->loopstart;
				int ml_prev_index = ml_last_index;

				MLoop *ml_curr = &mloop[ml_curr_index];
				MLoop *ml_prev = &mloop[ml_prev_index];

				for (; ml_curr_index <= ml_last_index; ml_curr++, ml_curr_index++) {
					int *e2l_curr = edge_to_loops[ml_curr->e];
					int *e2l_prev = edge_to_loops[ml_prev->e];

					if (IS_EDGE_SHARP(e2l_curr)) {
						if (IS_EDGE_SHARP(e2l_curr) && IS_EDGE_SHARP(e2l_prev)) {
							loop_split_worker(wnmd, wn_data, ml_curr, ml_prev, ml_curr_index, -1, NULL,
							                  mp_index, loop_normal, edge_to_loops);
						}
						else {
							loop_split_worker(wnmd, wn_data, ml_curr, ml_prev, ml_curr_index, ml_prev_index, e2l_prev,
							                  mp_index, loop_normal, edge_to_loops);
						}
					}
					ml_prev = ml_curr;
					ml_prev_index = ml_curr_index;
				}
			}
			MEM_freeN(edge_to_loops);
		}
		BKE_mesh_normals_loop_custom_set(mvert, numVerts, medge, numEdges,
		                                 mloop, loop_normal, numLoops, mpoly, polynors, numPolys, clnors);

		MEM_freeN(loop_normal);
	}
}

static void wn_face_area(WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data)
{
	const int numPolys = wn_data->numPolys;

	MVert *mvert = wn_data->mvert;
	MLoop *mloop = wn_data->mloop;
	MPoly *mpoly = wn_data->mpoly;

	MPoly *mp;
	int mp_index;

	ModePair *face_area = MEM_malloc_arrayN((size_t)numPolys, sizeof(*face_area), __func__);

	ModePair *f_area = face_area;
	for (mp_index = 0, mp = mpoly; mp_index < numPolys; mp_index++, mp++, f_area++) {
		f_area->val = BKE_mesh_calc_poly_area(mp, &mloop[mp->loopstart], mvert);
		f_area->index = mp_index;
	}

	qsort(face_area, numPolys, sizeof(*face_area), modepair_cmp_by_val_inverse);

	wn_data->mode_pair = face_area;
	apply_weights_vertex_normal(wnmd, wn_data);
}

static void wn_corner_angle(WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data)
{
	const int numLoops = wn_data->numLoops;
	const int numPolys = wn_data->numPolys;

	MVert *mvert = wn_data->mvert;
	MLoop *mloop = wn_data->mloop;
	MPoly *mpoly = wn_data->mpoly;

	MPoly *mp;
	int mp_index;

	int *loop_to_poly = MEM_malloc_arrayN((size_t)numLoops, sizeof(*loop_to_poly), __func__);

	ModePair *corner_angle = MEM_malloc_arrayN((size_t)numLoops, sizeof(*corner_angle), __func__);

	for (mp_index = 0, mp = mpoly; mp_index < numPolys; mp_index++, mp++) {
		MLoop *ml_start = &mloop[mp->loopstart];

		float *index_angle = MEM_malloc_arrayN((size_t)mp->totloop, sizeof(*index_angle), __func__);
		BKE_mesh_calc_poly_angles(mp, ml_start, mvert, index_angle);

		ModePair *c_angl = &corner_angle[mp->loopstart];
		float *angl = index_angle;
		for (int ml_index = mp->loopstart; ml_index < mp->loopstart + mp->totloop; ml_index++, c_angl++, angl++) {
			c_angl->val = (float)M_PI - *angl;
			c_angl->index = ml_index;

			loop_to_poly[ml_index] = mp_index;
		}
		MEM_freeN(index_angle);
	}

	qsort(corner_angle, numLoops, sizeof(*corner_angle), modepair_cmp_by_val_inverse);

	wn_data->loop_to_poly = loop_to_poly;
	wn_data->mode_pair = corner_angle;
	apply_weights_vertex_normal(wnmd, wn_data);
}

static void wn_face_with_angle(WeightedNormalModifierData *wnmd, WeightedNormalData *wn_data)
{
	const int numLoops = wn_data->numLoops;
	const int numPolys = wn_data->numPolys;

	MVert *mvert = wn_data->mvert;
	MLoop *mloop = wn_data->mloop;
	MPoly *mpoly = wn_data->mpoly;

	MPoly *mp;
	int mp_index;

	int *loop_to_poly = MEM_malloc_arrayN((size_t)numLoops, sizeof(*loop_to_poly), __func__);

	ModePair *combined = MEM_malloc_arrayN((size_t)numLoops, sizeof(*combined), __func__);

	for (mp_index = 0, mp = mpoly; mp_index < numPolys; mp_index++, mp++) {
		MLoop *ml_start = &mloop[mp->loopstart];

		float face_area = BKE_mesh_calc_poly_area(mp, ml_start, mvert);
		float *index_angle = MEM_malloc_arrayN((size_t)mp->totloop, sizeof(*index_angle), __func__);
		BKE_mesh_calc_poly_angles(mp, ml_start, mvert, index_angle);

		ModePair *cmbnd = &combined[mp->loopstart];
		float *angl = index_angle;
		for (int ml_index = mp->loopstart; ml_index < mp->loopstart + mp->totloop; ml_index++, cmbnd++, angl++) {
			/* In this case val is product of corner angle and face area. */
			cmbnd->val = ((float)M_PI - *angl) * face_area;
			cmbnd->index = ml_index;

			loop_to_poly[ml_index] = mp_index;
		}
		MEM_freeN(index_angle);
	}

	qsort(combined, numLoops, sizeof(*combined), modepair_cmp_by_val_inverse);

	wn_data->loop_to_poly = loop_to_poly;
	wn_data->mode_pair = combined;
	apply_weights_vertex_normal(wnmd, wn_data);
}

static DerivedMesh *applyModifier(ModifierData *md, Object *ob, DerivedMesh *dm, ModifierApplyFlag UNUSED(flag))
{
	WeightedNormalModifierData *wnmd = (WeightedNormalModifierData *)md;

	Mesh *me = ob->data;

	if (!(me->flag & ME_AUTOSMOOTH)) {
		modifier_setError((ModifierData *)wnmd, "Enable 'Auto Smooth' option in mesh settings");
		return dm;
	}

	const int numVerts = dm->getNumVerts(dm);
	const int numEdges = dm->getNumEdges(dm);
	const int numLoops = dm->getNumLoops(dm);
	const int numPolys = dm->getNumPolys(dm);

	MPoly *mpoly = dm->getPolyArray(dm);
	MVert *mvert = dm->getVertArray(dm);
	MEdge *medge = dm->getEdgeArray(dm);
	MLoop *mloop = dm->getLoopArray(dm);

	bool free_polynors = false;

	/* Add some comment here about why this is needed? */
	float weight = ((float)wnmd->weight) / 50.0f;
	if (wnmd->weight == 100) {
		weight = (float)SHRT_MAX;
	}
	else if (wnmd->weight == 1) {
		weight = 1 / (float)SHRT_MAX;
	}
	else if ((weight - 1) * 25 > 1) {
		weight = (weight - 1) * 25;
	}

	float (*polynors)[3] = dm->getPolyDataArray(dm, CD_NORMAL);
	if (!polynors) {
		polynors = MEM_malloc_arrayN((size_t)numPolys, sizeof(*polynors), __func__);
		BKE_mesh_calc_normals_poly(mvert, NULL, numVerts, mloop, mpoly, numLoops, numPolys, polynors, false);
		free_polynors = true;
	}

	short (*clnors)[2];
	clnors = CustomData_duplicate_referenced_layer(&dm->loopData, CD_CUSTOMLOOPNORMAL, numLoops);
	if (!clnors) {
		DM_add_loop_layer(dm, CD_CUSTOMLOOPNORMAL, CD_CALLOC, NULL);
		clnors = dm->getLoopDataArray(dm, CD_CUSTOMLOOPNORMAL);
	}

	MDeformVert *dvert;
	int defgrp_index;
	modifier_get_vgroup(ob, dm, wnmd->defgrp_name, &dvert, &defgrp_index);

	WeightedNormalData wn_data = {
		.numVerts = numVerts,
		.numEdges = numEdges,
		.numLoops = numLoops,
		.numPolys = numPolys,

		.mvert = mvert,
		.medge = medge,

		.mloop = mloop,
		.clnors = clnors,

		.mpoly = mpoly,
		.polynors = polynors,
		.poly_strength = CustomData_get_layer_named(&dm->polyData, CD_PROP_INT, MOD_WEIGHTEDNORMALS_FACEWEIGHT_CDLAYER_ID),

		.dvert = dvert,
		.defgrp_index = defgrp_index,
		.use_invert_vgroup = (wnmd->flag & MOD_WEIGHTEDNORMAL_INVERT_VGROUP) != 0,

		.weight = weight,
		.mode = wnmd->mode,
	};

	switch (wnmd->mode) {
		case MOD_WEIGHTEDNORMAL_MODE_FACE:
			wn_face_area(wnmd, &wn_data);
			break;
		case MOD_WEIGHTEDNORMAL_MODE_ANGLE:
			wn_corner_angle(wnmd, &wn_data);
			break;
		case MOD_WEIGHTEDNORMAL_MODE_FACE_ANGLE:
			wn_face_with_angle(wnmd, &wn_data);
			break;
	}

	if (free_polynors) {
		MEM_freeN(polynors);
	}

	MEM_SAFE_FREE(wn_data.loop_to_poly);
	MEM_SAFE_FREE(wn_data.mode_pair);
	MEM_SAFE_FREE(wn_data.curr_vert_strength);
	MEM_SAFE_FREE(wn_data.curr_vert_val);
	MEM_SAFE_FREE(wn_data.vert_loops_count);
	MEM_SAFE_FREE(wn_data.vert_normals);

	return dm;
}

static void copyData(ModifierData *md, ModifierData *target)
{
	modifier_copyData_generic(md, target);
}

static void initData(ModifierData *md)
{
	WeightedNormalModifierData *wnmd = (WeightedNormalModifierData *)md;
	wnmd->mode = MOD_WEIGHTEDNORMAL_MODE_FACE;
	wnmd->weight = 50;
	wnmd->thresh = 1e-2f;
	wnmd->flag = 0;
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *md)
{
	WeightedNormalModifierData *wnmd = (WeightedNormalModifierData *)md;
	CustomDataMask dataMask = CD_CUSTOMLOOPNORMAL;

	if (wnmd->defgrp_name[0]) {
		dataMask |= CD_MASK_MDEFORMVERT;
	}

	if (wnmd->flag & MOD_WEIGHTEDNORMAL_FACE_INFLUENCE) {
		dataMask |= CD_MASK_PROP_INT;
	}

	return dataMask;
}

static bool dependsOnNormals(ModifierData *UNUSED(md))
{
	return true;
}

ModifierTypeInfo modifierType_WeightedNormal = {
	/* name */              "Weighted Normal",
	/* structName */        "WeightedNormalModifierData",
	/* structSize */        sizeof(WeightedNormalModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsMapping |
	                        eModifierTypeFlag_SupportsEditmode |
	                        eModifierTypeFlag_EnableInEditmode,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  requiredDataMask,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepgraph */    NULL,
	/* updateDepsgraph */   NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  dependsOnNormals,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
