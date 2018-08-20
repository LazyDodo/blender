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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Copyright (C) 2014, 2018 by Martin Felke.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/fracture.c
 *  \ingroup blenkernel
 */

#include <stdio.h>
#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_customdata.h"
#include "BKE_deform.h"
#include "BKE_DerivedMesh.h"
#include "BKE_fracture.h"
#include "BKE_fracture_util.h"
#include "BKE_global.h"
#include "BKE_material.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_particle.h"
#include "BKE_pointcache.h"
#include "BKE_rigidbody.h"

#include "BLI_edgehash.h"
#include "BLI_kdtree.h"
#include "BLI_listbase.h"
#include "BLI_math_vector.h"
#include "BLI_mempool.h"
#include "BLI_path_util.h"
#include "BLI_rand.h"
#include "BLI_string.h"
#include "BLI_sort.h"
#include "BLI_task.h"
#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_fracture_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_group_types.h"
#include "DNA_material_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_particle_types.h"

#include "DNA_object_types.h"
#include "DEG_depsgraph_query.h"

#include "bmesh.h"

#include "RBI_api.h"

/* debug timing */
#define USE_DEBUG_TIMER

#ifdef USE_DEBUG_TIMER
#include "PIL_time.h"
#endif

#ifdef WITH_VORO
#include "../../../../extern/voro++/src/c_interface.hh"
#endif

/* prototypes */
static MeshIsland *parse_cell(cell c);
static void parse_cell_verts(cell c, MVert *mvert, int totvert);
static void parse_cell_polys(cell c, MPoly *mpoly, int totpoly);
static void parse_cell_loops(cell c, MLoop *mloop, MPoly *mpoly, int totpoly);
static void parse_cell_neighbors(cell c, int *neighbors, int totpoly);
static void do_island_index_map(FractureModifierData *fmd, Object *ob);
static void do_rigidbody(Main* bmain, Scene* scene, MeshIsland* mi, Object* ob, short rb_type, int i);


static void fracture_meshisland_add(FractureModifierData *fmd, MeshIsland *mi)
{
	MVert *mv;
	int i;

	mul_m4_v3(fmd->shared->splinter_matrix, mi->centroid);
	for (i = 0, mv = mi->mesh->mvert; i < mi->mesh->totvert; i++, mv++ )
	{
		mul_m4_v3(fmd->shared->splinter_matrix, mv->co);
	}

	BLI_addtail(&fmd->shared->mesh_islands, mi);
}


static int mesh_sortsize(const void *s1, const void *s2, void* UNUSED(context))
{
	Mesh **me1 = (Mesh **)s1;
	Mesh **me2 = (Mesh **)s2;

	float size1[3], size2[3], loc[3];
	float val_a,  val_b;

	if ((*me1 == NULL) || (*me2 == NULL)) {
		return -1;
	}

	BKE_fracture_mesh_boundbox_calc(*me1, loc, size1);
	BKE_fracture_mesh_boundbox_calc(*me2, loc, size2);

	//squared diameter
	val_a = size1[0]*size1[0] + size1[1]*size1[1] + size1[2]*size1[2];
	val_b = size2[0]*size2[0] + size2[1]*size2[1] + size2[2]*size2[2];

	/* sort */
	if      (val_a < val_b) return -1;
	else if (val_a > val_b) return 1;
	return 0;
}

#if 0
/* copied from mesh_evaluate.c */
/**
 * Calculate the volume and volume-weighted centroid of the volume formed by the polygon and the origin.
 * Results will be negative if the origin is "outside" the polygon
 * (+ve normal side), but the polygon may be non-planar with no effect.
 *
 * Method from:
 * - http://forums.cgsociety.org/archive/index.php?t-756235.html
 * - http://www.globalspec.com/reference/52702/203279/4-8-the-centroid-of-a-tetrahedron
 *
 * \note volume is 6x actual volume, and centroid is 4x actual volume-weighted centroid
 * (so division can be done once at the end)
 * \note results will have bias if polygon is non-planar.
 */
static float mesh_calc_poly_volume_and_weighted_centroid(
		const MPoly *mpoly, const MLoop *loopstart, const MVert *mvarray,
		float r_cent[3])
{
	const float *v_pivot, *v_step1;
	float total_volume = 0.0f;

	zero_v3(r_cent);

	v_pivot = mvarray[loopstart[0].v].co;
	v_step1 = mvarray[loopstart[1].v].co;

	for (int i = 2; i < mpoly->totloop; i++) {
		const float *v_step2 = mvarray[loopstart[i].v].co;

		/* Calculate the 6x volume of the tetrahedron formed by the 3 vertices
		 * of the triangle and the origin as the fourth vertex */
		float v_cross[3];
		cross_v3_v3v3(v_cross, v_pivot, v_step1);
		const float tetra_volume = dot_v3v3 (v_cross, v_step2);
		total_volume += tetra_volume;

		/* Calculate the centroid of the tetrahedron formed by the 3 vertices
		 * of the triangle and the origin as the fourth vertex.
		 * The centroid is simply the average of the 4 vertices.
		 *
		 * Note that the vector is 4x the actual centroid so the division can be done once at the end. */
		for (uint j = 0; j < 3; j++) {
			r_cent[j] += tetra_volume * (v_pivot[j] + v_step1[j] + v_step2[j]);
		}

		v_step1 = v_step2;
	}

	return total_volume;
}

#endif

/* note, results won't be correct if polygon is non-planar */
/* copied from mesh_evaluate.c */
static float mesh_calc_poly_planar_area_centroid(
		const MPoly *mpoly, const MLoop *loopstart, const MVert *mvarray,
		float r_cent[3])
{
	int i;
	float tri_area;
	float total_area = 0.0f;
	float v1[3], v2[3], v3[3], normal[3], tri_cent[3];

	BKE_mesh_calc_poly_normal(mpoly, loopstart, mvarray, normal);
	copy_v3_v3(v1, mvarray[loopstart[0].v].co);
	copy_v3_v3(v2, mvarray[loopstart[1].v].co);
	zero_v3(r_cent);

	for (i = 2; i < mpoly->totloop; i++) {
		copy_v3_v3(v3, mvarray[loopstart[i].v].co);

		tri_area = area_tri_signed_v3(v1, v2, v3, normal);
		total_area += tri_area;

		mid_v3_v3v3v3(tri_cent, v1, v2, v3);
		madd_v3_v3fl(r_cent, tri_cent, tri_area);

		copy_v3_v3(v2, v3);
	}

	mul_v3_fl(r_cent, 1.0f / total_area);

	return total_area;
}

// old method, keep for now in case new has different results
/* modified from BKE_mesh_center_centroid */
bool BKE_fracture_mesh_center_centroid_area(Mesh *shard, float cent[3])
{
	int i = shard->totpoly;
	MPoly *mpoly;
	float poly_area;
	float total_area = 0.0f;
	float poly_cent[3];

	zero_v3(cent);

	/* calculate a weighted average of polygon centroids */
	for (mpoly = shard->mpoly; i--; mpoly++) {
		BKE_mesh_calc_poly_center(mpoly, shard->mloop + mpoly->loopstart, shard->mvert, poly_cent);
//		poly_area = BKE_mesh_calc_poly_area(mpoly, shard->mloop + mpoly->loopstart, shard->mvert);
		poly_area = mesh_calc_poly_planar_area_centroid(mpoly, shard->mloop + mpoly->loopstart, shard->mvert,
														poly_cent);
		madd_v3_v3fl(cent, poly_cent, poly_area);
		total_area += poly_area;
	}
	/* otherwise we get NAN for 0 polys */
	if (shard->totpoly) {
		mul_v3_fl(cent, 1.0f / total_area);
	}

	/* zero area faces cause this, fallback to median */
	if (UNLIKELY(!is_finite_v3(cent))) {
		return BKE_mesh_center_median(shard, cent);
	}

	return (shard->totpoly != 0);
}

static void calculate_fast_bisect(FractureModifierData *fmd, Mesh* me, BisectContext *ctx)
{
	float factor = 1 - fmd->orthogonality_factor;
	float vec[3];
	float loc[3], size[3];
	int max_axis;

	BKE_fracture_mesh_boundbox_calc(me, loc, size);

	//make a random vector (interpret as cutter plane)
	vec[0] = BLI_thread_frand(0) * 2 - 1;
	vec[1] = BLI_thread_frand(0) * 2 - 1;
	vec[2] = BLI_thread_frand(0) * 2 - 1;

	//multiply two minor dimensions with a factor to emphasize the max dimension
	max_axis = axis_dominant_v3_single(size);
	switch (max_axis) {
		case 0:
			vec[1] *= factor;
			vec[2] *= factor;
			break;
		case 1:
			vec[0] *= factor;
			vec[2] *= factor;
			break;
		case 2:
			vec[0] *= factor;
			vec[1] *= factor;
			break;
	}

	copy_v3_v3(ctx->normal, vec);
	BKE_fracture_mesh_center_centroid_area(me, ctx->centroid);
}

static void calculate_fractal(FractureModifierData* fmd, Mesh* me, BooleanContext *ctx)
{
	float factor = 1 - fmd->orthogonality_factor;
	float radius;
	float size[3];
	float quat[4];
	float loc[3], vec[3];
	float one[3] = {1.0f, 1.0f, 1.0f};
	float matrix[4][4];
	int max_axis;

	BKE_fracture_mesh_boundbox_calc(me, loc, size);
	radius = sqrt(size[0]*size[0] + size[1]*size[1] + size[2]*size[2]) * 1.5f;

	vec[0] = BLI_thread_frand(0) * 2 - 1;
	vec[1] = BLI_thread_frand(0) * 2 - 1;
	vec[2] = BLI_thread_frand(0) * 2 - 1;

	//multiply two minor dimensions with a factor to emphasize the max dimension
	max_axis = axis_dominant_v3_single(size);
	switch (max_axis) {
		case 0:
			vec[1] *= factor;
			vec[2] *= factor;
			break;
		case 1:
			vec[0] *= factor;
			vec[2] *= factor;
			break;
		case 2:
			vec[0] *= factor;
			vec[1] *= factor;
			break;
	}

	//printf("(%f %f %f) (%f %f %f) \n", size[0], size[1], size[2], eul[0], eul[1], eul[2]);*/
	//loc_eul_size_to_mat4(matrix, loc, vec, one);
	vec_to_quat(quat, vec, OB_POSZ, OB_POSX);
	loc_quat_size_to_mat4(matrix, loc, quat, one);

	ctx->cutter_plane_radius = radius;
	copy_m4_m4(ctx->cutter_plane_matrix, matrix);
}

static bool needs_process(FractureModifierData *fmd, MeshIsland* mi)
{
	//check against tree ? TODO
	//here comes the logic from prepare cells, compare locations against tree etc
	int n, j;
	float max = 0.0f;
	KDTreeNearest nearest;

	if (!fmd->shared->last_islands || !fmd->shared->last_island_tree)
	{
		return true;
	}

	/* check how far the farthest neighbor centroid is */
	for (n = 0; n < mi->neighbor_count; n++)
	{
		int index = mi->neighbors[n];
		if (index > -1 && index < fmd->shared->last_expected_islands)
		{
			MeshIsland *mii = fmd->shared->last_islands[index];
			float dist = len_squared_v3v3(mi->raw_centroid, mii->raw_centroid);
			if (dist > max)
			{
				max = dist;
			}
		}
	}

	/* check in a distance around us (all neighbors are sure to be included) whether our neighborhood changed*/
	j = BLI_kdtree_find_nearest(fmd->shared->last_island_tree, mi->raw_centroid, &nearest);
	if (j > -1 && j < fmd->shared->last_expected_islands)
	{
		float epsilon = 0.00001;
		MeshIsland *mii = fmd->shared->last_islands[j];
		if (mii != NULL && nearest.dist < max)
		{
			if (nearest.dist < epsilon) {
				if ((fabsf(mi->raw_volume - mii->raw_volume) < epsilon))
				{
					return false;
				}
			}
		}
	}
	return true;
}

static void prepare_boolean(FractureModifierData* fmd, Object* ob, BooleanContext* ctx)
{
	ctx->use_fractal = false;
	ctx->operation = 0;
	ctx->thresh = fmd->boolean_double_threshold;
	BLI_strncpy(ctx->uv_layer, fmd->uvlayer_name, 64);
	ctx->inner_material_index = BKE_object_material_slot_find_index(ob, fmd->inner_material) - 1;
	/*if no inner material has been found, just pick the first one */
	if (ctx->inner_material_index < 0)
		ctx->inner_material_index = 0;
}


static void prepare_boolean_fractal(FractureModifierData* fmd, Object* ob, Mesh* me, BooleanContext* ctx)
{
	prepare_boolean(fmd, ob, ctx);
	calculate_fractal(fmd, me, ctx);
	ctx->use_fractal = true;
	ctx->fractal_amount = fmd->fractal_amount;
	ctx->num_cuts = fmd->fractal_cuts;
	ctx->num_iterations = fmd->fractal_iterations;
	ctx->use_smooth_inner = fmd->use_smooth;
}


static void prepare_bisect(FractureModifierData *fmd, Object* ob, BisectContext* ctx)
{
	ctx->clear_inner = false;
	ctx->clear_outer = true;
	ctx->do_fast_bisect = false;
	ctx->geometry_limitation_tree = fmd->shared->last_island_tree;
	ctx->inner_material_index = BKE_object_material_slot_find_index(ob, fmd->inner_material) - 1;

	/*if no inner material has been found, just pick the first one */
	if (ctx->inner_material_index < 0)
		ctx->inner_material_index = 0;

	ctx->use_fill = false;
	BLI_strncpy(ctx->uv_layer, fmd->uvlayer_name, 64);
	unit_m4(ctx->obmat); // hmm, why necessary ?
}

static void prepare_bisect_fill(FractureModifierData *fmd, Object* ob, BisectContext* ctx)
{
	prepare_bisect(fmd, ob, ctx);
	ctx->use_fill = true;
}

static void prepare_fast_bisect(FractureModifierData *fmd, Object *ob, Mesh* me, BisectContext* ctx)
{
	prepare_bisect(fmd, ob, ctx);
	calculate_fast_bisect(fmd, me, ctx);
	ctx->do_fast_bisect = true;
	ctx->use_fill = false;
}

static void prepare_fast_bisect_fill(FractureModifierData *fmd, Object* ob, Mesh* me, BisectContext *ctx)
{
	prepare_fast_bisect(fmd, ob, me, ctx);
	ctx->use_fill = true;
}

static Mesh* get_mesh(Mesh** meshes, int index, Mesh* mesh)
{
	if (index == 0) {
		return mesh;
	}
	else {
		return meshes[index];
	}
}

static void process_cells(FractureModifierData* fmd, MeshIsland* mi, Main* bmain, Object* ob, Scene* scene, cell *c, int count)
{
	int i, j = 0;
	BisectContext bictx = {};
	BooleanContext boctx = {};
	MeshIsland** islands = NULL;
	KDTree *tree = NULL;
	Mesh** temp_meshs = NULL;
	Mesh* me = NULL, *mesh = mi->mesh;
	int count_new = count+1;
	float frame = BKE_scene_frame_get(scene);

	/*global preparations */
	islands = MEM_callocN(sizeof(MeshIsland*) * count, "islands");
	tree = BLI_kdtree_new(count);
	temp_meshs = MEM_callocN(sizeof(Mesh*) * count_new, "temp_meshs");

	mi->endframe = frame;

	/*for each cell...*/
//#pragma omp parallel for
	for (i = 0; i < count; i++)
	{
		/* parse to raw meshisland*/
		MeshIsland *mi = parse_cell(c[i]);

		BLI_kdtree_insert(tree, count, mi->centroid);
		islands[i] = mi;

		/* check whether it needs to be processed */
		if (needs_process(fmd, mi))
		{
			/* meshB is for "halving" algorithms like fractal and bisectfast/bisectfastfill*/
			Mesh *meshA = NULL, *meshB = NULL;

			/* process according to algorithm */
			switch (fmd->frac_algorithm) {
				case MOD_FRACTURE_BOOLEAN:
					prepare_boolean(fmd, ob, &boctx);
					meshA = BKE_fracture_mesh_boolean(mesh, mi->mesh, ob, &boctx);
					break;

				case MOD_FRACTURE_BOOLEAN_FRACTAL:
					me = get_mesh(temp_meshs, i, mesh);
					if (me) {
						prepare_boolean_fractal(fmd, ob, me, &boctx);
						BKE_fracture_mesh_boolean_fractal(me, &meshA, &meshB, ob, &boctx);
					}
					break;

				case MOD_FRACTURE_BISECT:
					prepare_bisect(fmd, ob, &bictx);
					meshA = BKE_fracture_mesh_bisect(mesh, mi, &bictx);
					break;

				case MOD_FRACTURE_BISECT_FILL:
					prepare_bisect_fill(fmd, ob, &bictx);
					meshA = BKE_fracture_mesh_bisect(mesh, mi, &bictx);
					break;

				case MOD_FRACTURE_BISECT_FAST:
					me = get_mesh(temp_meshs, i, mesh);
					prepare_fast_bisect(fmd, ob, me, &bictx);
					BKE_fracture_mesh_bisect_fast(me, &meshA, &meshB, &bictx);
					break;

				case MOD_FRACTURE_BISECT_FAST_FILL:
					me = get_mesh(temp_meshs, i, mesh);
					prepare_fast_bisect_fill(fmd, ob, me, &bictx);
					BKE_fracture_mesh_bisect_fast(me, &meshA, &meshB, &bictx);
					break;
			}

			/* if successful, create processed meshisland in FM */
			if (temp_meshs[i]) {
				BKE_fracture_mesh_free(temp_meshs[i]);
				temp_meshs[i] = NULL;
			}
			if (temp_meshs[i+1]) {
				BKE_fracture_mesh_free(temp_meshs[i+1]);
				temp_meshs[i+1] = NULL;
			}

			if (meshA != me) {
				temp_meshs[i] = meshA;
			}

			if (meshB != me) {
				temp_meshs[i+1] = meshB;
			}

			/*sort meshs by size*/
			if(fmd->frac_algorithm == MOD_FRACTURE_BISECT_FAST ||
			   fmd->frac_algorithm == MOD_FRACTURE_BISECT_FAST_FILL ||
			   fmd->frac_algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL)
			{
				BLI_qsort_r(temp_meshs, i+2, sizeof(Mesh *), mesh_sortsize, NULL);
			}
		}
	}

	if (fmd->split_islands)
	{
		int diff = 1;
		for (i = 0; i < count+1; i++)
		{
			if (temp_meshs[i]) {
				BKE_fracture_split_islands(fmd, ob, temp_meshs[i], &temp_meshs, &count_new );
				BKE_fracture_mesh_free(temp_meshs[i]);
				temp_meshs[i] = NULL;
			}

			diff = count_new - (count+1);

			if (diff > 1) {
				if (temp_meshs[i+diff-1]) {
					BKE_fracture_split_islands(fmd, ob, temp_meshs[i+diff-1], &temp_meshs, &count_new );
					BKE_fracture_mesh_free(temp_meshs[i+diff-1]);
					temp_meshs[i+diff-1] = NULL;
				}
			}
		}
	}

	for (i = 0; i < count_new; i++)
	{
		if (temp_meshs[i])
		{
			if (temp_meshs[i]->totvert > 0)
			{	/* skip invalid cells, e.g. those which are eliminated by bisect */
				MeshIsland *result = BKE_fracture_mesh_island_create(temp_meshs[i], bmain, scene, ob, frame);
				fracture_meshisland_add(fmd, result);
				result->id = j;
				j++;
			}
			else {
				BKE_fracture_mesh_free(temp_meshs[i]);
			}
		}
	}

	BLI_kdtree_balance(tree);

	MEM_freeN(temp_meshs);

	/* swap old last islands and tree against new for next run */
	if (fmd->shared->last_island_tree)
	{
		BLI_kdtree_free(fmd->shared->last_island_tree);
	}

	fmd->shared->last_island_tree = tree;

	if (fmd->shared->last_islands)
	{
		int k = 0;
		for (k = 0; k < fmd->shared->last_expected_islands; k++)
		{
			BKE_fracture_mesh_island_free(fmd->shared->last_islands[k], scene);
		}

		MEM_freeN(fmd->shared->last_islands);
	}

	fmd->shared->last_islands = islands;
	fmd->shared->last_expected_islands = count;
}

static MeshIsland *parse_cell(cell c)
{
	MeshIsland *mi = MEM_callocN(sizeof(MeshIsland), "mi_cell");
	Mesh* me = BKE_mesh_new_nomain(c.totvert, 0, 0, c.totloop, c.totpoly);

	int totpoly = 0, totloop = 0, totvert = 0;
	float centr[3];

	mi->mesh = me;

	totvert = c.totvert;
	if (totvert > 0) {
		parse_cell_verts(c, me->mvert, totvert);
	}

	totpoly = c.totpoly;
	if (totpoly > 0) {
		parse_cell_polys(c, me->mpoly, totpoly);
	}

	totloop = c.totloop;
	if (totloop > 0) {
		parse_cell_loops(c, me->mloop, me->mpoly, totpoly);
	}

	if (totpoly > 0) {
		mi->neighbors = MEM_callocN(sizeof(int) * totpoly, __func__);
		mi->neighbor_count = totpoly;
		parse_cell_neighbors(c, mi->neighbors, totpoly);
	}

	copy_v3_v3(centr, c.centroid);
	copy_v3_v3(mi->centroid, centr);
	copy_v3_v3(mi->raw_centroid, centr);
	mi->raw_volume = c.volume;

	BKE_mesh_calc_edges(me, true, true);
	BKE_mesh_calc_normals(me);

	return mi;
}

static void parse_cell_verts(cell c, MVert *mvert, int totvert)
{
	int i;

	for (i = 0; i < totvert; i++) {
		float *co = mvert[i].co;
		copy_v3_v3(co, c.verts[i]);
	}
}

static void parse_cell_polys(cell c, MPoly *mpoly, int totpoly)
{
	int i;
	int totloop = 0;

	for (i = 0; i < totpoly; ++i) {
		int numloop;

		numloop = c.poly_totvert[i];

		mpoly[i].loopstart = totloop;
		mpoly[i].totloop = numloop;

		totloop += numloop;
	}
}

static void parse_cell_loops(cell c, MLoop *mloop, MPoly *mpoly, int totpoly)
{
	int i, k;

	for (i = 0; i < totpoly; ++i) {
		int loopstart = mpoly[i].loopstart;
		int numloop = mpoly[i].totloop;

		for (k = 0; k < numloop; ++k) {
			int index;

			index = c.poly_indices[i][k];

			/* note: invert vertex order here,
			 * otherwise normals are pointing inward
			 */
			mloop[loopstart + (numloop - 1) - k].v = index;
		}
	}
}

static void parse_cell_neighbors(cell c, int *neighbors, int totpoly)
{
	int i;

	for (i = 0; i < totpoly; i++) {
		int n;
		n = c.neighbors[i];
		neighbors[i] = n;
	}
}

static void stroke_to_faces(FractureModifierData *fmd, BMesh** bm, bGPDstroke *gps, int inner_material_index)
{
	BMVert *lastv1 = NULL;
	BMVert *lastv2 = NULL;
	int p = 0;
	float thresh = (float)fmd->grease_decimate / 100.0f;
	float half[3] = {0, 0, 1};

	for (p = 0; p < gps->totpoints; p++) {

		if ((BLI_thread_frand(0) < thresh) || (p == 0) || (p == gps->totpoints-1)) {
			BMVert *v1, *v2;
			float point[3] = {0, 0, 0};

			point[0] = gps->points[p].x;
			point[1] = gps->points[p].y;
			point[2] = gps->points[p].z;

			v1 = BM_vert_create(*bm, point, NULL, 0);

			if (lastv1)
			{
				BMFace* f;
				float nvec[3] = {0.0f, 0.0f, 0.0f}, co1[3], co2[3];

				/*also "extrude" this along the normal, no...use global axises instead*/
				if (fmd->cutter_axis == MOD_FRACTURE_CUTTER_X)
				{
					nvec[0] = 1.0f;
					nvec[1] = 0.0f;
					nvec[2] = 0.0f;
				}

				if (fmd->cutter_axis == MOD_FRACTURE_CUTTER_Y)
				{
					nvec[0] = 0.0f;
					nvec[1] = 1.0f;
					nvec[2] = 0.0f;
				}

				if (fmd->cutter_axis == MOD_FRACTURE_CUTTER_Z)
				{
					nvec[0] = 0.0f;
					nvec[1] = 0.0f;
					nvec[2] = 1.0f;
				}

				mul_v3_fl(nvec, fmd->grease_offset);
				mul_v3_v3fl(half, nvec, 0.5f);

				add_v3_v3v3(co1, v1->co, nvec);
				v2 = BM_vert_create(*bm, co1, NULL, 0);

				if (!lastv2)
				{
					add_v3_v3v3(co2, lastv1->co, nvec);
					lastv2 = BM_vert_create(*bm, co2, NULL, 0);
				}

				f = BM_face_create_quad_tri(*bm, lastv1, v1, v2, lastv2, NULL, 0);
				f->mat_nr = inner_material_index;
				lastv2 = v2;
			}

			lastv1 = v1;
		}
	}

	{
		/* move the stroke mesh a bit out, half of offset */
		BMIter iter;
		BMVert *v;

		BM_ITER_MESH(v, &iter, *bm, BM_VERTS_OF_MESH)
		{
			sub_v3_v3(v->co, half);
		}

		BM_mesh_elem_hflag_enable_all(*bm, BM_FACE | BM_EDGE | BM_VERT, BM_ELEM_SELECT, false);
		BMO_op_callf(*bm, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
					 "remove_doubles verts=%av dist=%f", BM_VERTS_OF_MESH, 0.01, false);
	}
}

static void fracture_shard_vgroup_add(MeshIsland *s, Object *ob, const char* name)
{
	int index = 0, i = 0;
	MDeformVert *dvert;
	if (!defgroup_find_name(ob, name)) {
		BKE_defgroup_new(ob, name);
	}
	index = defgroup_name_index(ob, name);
	dvert = CustomData_get_layer(&s->mesh->vdata, CD_MDEFORMVERT);
	if (dvert == NULL) {
		dvert = CustomData_add_layer(&s->mesh->vdata, CD_MDEFORMVERT, CD_CALLOC, NULL, s->mesh->totvert);
	}
	for (i = 0; i < s->mesh->totvert; i++) {
		MDeformVert* dv = dvert + i;
		defvert_add_index_notest(dv, index, 1.0f);
	}
}

static void fracture_shard_material_add(MeshIsland* s, Object *ob, short mat_ofs) {

	/* only use material offsets if we have 3 or more materials; since FM material handling still is a bit odd  */
	/* todo, use index for inner material too, then this hack here isnt necessary any more */

	const short mat_nr_max = ob->totcol > 1 ? ob->totcol - 1 : 0;
	mat_ofs = mat_nr_max ? mat_ofs : 0;

	if (mat_ofs) {
		MPoly *mp;
		int i = 0;

		for (i = 0; i < s->mesh->totpoly; i++)
		{
			mp = s->mesh->mpoly + i;
			mp->mat_nr = mat_ofs;

			//hrm first material and second (inner) should be untouched...
			CLAMP(mp->mat_nr, 0, mat_nr_max);
		}
	}
}

#if 0 //TODO FIX
static void do_intersect(FractureModifierData *fmd, Object* ob, MeshIsland *t, short inner_mat_index,
						 bool is_zero, float mat[4][4], int **shard_counts, int* count,
						 int k, Mesh **dm_parent, bool keep_other_shard, float thresh)
{
	/*just keep appending items at the end here */
	MPoly *mpoly, *mp;
	int totpoly;
	MeshIsland *parent = *dm_parent;
	MeshIsland *s = NULL, *s2 = NULL;
	int shards = 0, j = 0;

	if (is_zero == false && *dm_parent == NULL) {
		parent = BLI_findlink(&fmd->shared->frac_mesh->shard_map, k);
		*dm_parent = BKE_fracture_shard_to_mesh(parent, true);
	}

	mpoly = (*dm_parent)->mpoly;
	totpoly = (*dm_parent)->totpoly;

	for (j = 0, mp = mpoly; j < totpoly; j++, mp++) {
		mp->flag &= ~ME_FACE_SEL;
	}

	if (keep_other_shard)
	{
		s = BKE_fracture_shard_boolean(ob, *dm_parent, t, inner_mat_index, 0, 0.0f, &s2, NULL, 0.0f, false, 0, fmd->uvlayer_name, thresh);
	}
	else
	{
		s = BKE_fracture_shard_boolean(ob, *dm_parent, t, inner_mat_index, 0, 0.0f, NULL, NULL, 0.0f, false, 0, fmd->uvlayer_name, thresh);
	}

	//printf("Fractured: %d\n", k);

	if (ELEM(fmd->keep_cutter_shards, MOD_FRACTURE_KEEP_BOTH, MOD_FRACTURE_KEEP_INTERSECT)) {
		if (s != NULL) {
			fracture_shard_vgroup_add(s, ob, "Intersect");
			fracture_shard_material_add(s, ob, fmd->mat_ofs_intersect);
			fracture_shard_add(fmd->shared->frac_mesh, s, mat);
			shards++;
			s = NULL;
		}
	}

	if (ELEM(fmd->keep_cutter_shards, MOD_FRACTURE_KEEP_BOTH, MOD_FRACTURE_KEEP_DIFFERENCE)) {
		if (s2 != NULL) {
			fracture_shard_vgroup_add(s2, ob, "Difference");
			fracture_shard_material_add(s2, ob, fmd->mat_ofs_difference);
			fracture_shard_add(fmd->shared->frac_mesh, s2, mat);
			shards++;
			s2 = NULL;
		}
	}

	//if ((is_zero && ob->derivedFinal == NULL) || !is_zero) {
	{
		if (is_zero) {
			*count = 0;
		}

		BKE_mesh_free(*dm_parent);
		*dm_parent = NULL;
	}

	if (is_zero) {
		shards = 0;
	}

	(*shard_counts)[k] = shards;
	//printf("k, shards: %d %d \n", k, shards);
	//shards = 0;
}
#endif


#if 0  //TODO FIX !!!!
static void intersect_shards_by_dm(FractureModifierData *fmd, Mesh *d, Object *ob, Object *ob2, short inner_mat_index,
								   float mat[4][4], bool keep_other_shard, float thresh)
{
	Shard *t = NULL;
	int i = 0, count = 0, k = 0;
	float imat[4][4];
	int* shard_counts = NULL;
	bool is_zero = false;
	MVert *mv;
	Mesh *dm_parent = NULL;

	t = BKE_fracture_shard_create(d->mvert, d->mpoly, d->mloop, d->medge,
								  d->totvert, d->totpoly, d->totloop, d->totedge, true);
	BKE_fracture_custom_data_mesh_to_shard(t, d);


	invert_m4_m4(imat, ob->obmat);
	for (i = 0, mv = t->mvert; i < t->totvert; mv++, i++){
		if (ob2)
			mul_m4_v3(ob2->obmat, mv->co);
		mul_m4_v3(imat, mv->co);
	}

	count = fmd->shared->frac_mesh->shard_count;

	/*TODO, pass modifier mesh here !!! */
	if (count == 0 && keep_other_shard) {

		if (ob->runtime.mesh_eval != NULL) {
			dm_parent = BKE_fracture_mesh_copy(ob->runtime.mesh_eval, ob);
		}

		if (dm_parent == NULL) {
			dm_parent = BKE_fracture_mesh_copy(ob->data, ob);
		}

		count = 1;
		is_zero = true;
	}

	shard_counts = MEM_mallocN(sizeof(int) * count, "shard_counts");

	for (k = 0; k < count; k++) {
		do_intersect(fmd, ob, t, inner_mat_index, is_zero, mat, &shard_counts, &count, k, &dm_parent, keep_other_shard, thresh);
	}

	for (k = 0; k < count; k++)
	{
		int cnt = shard_counts[k];

		if (cnt > 0)
		{
			if (keep_other_shard)
			{
				/*clean up old entries here to avoid unnecessary shards*/
				Shard *first = fmd->shared->frac_mesh->shard_map.first;
				BLI_remlink_safe(&fmd->shared->frac_mesh->shard_map,first);
				BKE_fracture_shard_free(first, true);
				first = NULL;
			}

			/* keep asynchronous by intent, to keep track of original shard count */
			fmd->shared->frac_mesh->shard_count--;
		}
	}

	MEM_freeN(shard_counts);
	shard_counts = NULL;

	BKE_fracture_shard_free(t, true);
}
#endif

void BKE_fracture_shard_by_greasepencil(FractureModifierData *fmd, Object *obj, short inner_material_index, float mat[4][4])
{
	bGPDlayer *gpl;
	bGPDframe *gpf;
	bGPDstroke *gps;

	//reset_shards(fmd);

	if ((obj->gpd) && (obj->gpd->layers.first)) {

		float imat[4][4];
		invert_m4_m4(imat, mat);
		for (gpl = obj->gpd->layers.first; gpl; gpl = gpl->next) {
			for (gpf = gpl->frames.first; gpf; gpf = gpf->next) {
				for (gps = gpf->strokes.first; gps; gps = gps->next) {
					BMesh *bm = BM_mesh_create(&bm_mesh_allocsize_default, &((struct BMeshCreateParams){.use_toolflags = true,}));
					Mesh *dm;

					/*create stroke mesh */
					stroke_to_faces(fmd, &bm, gps, inner_material_index);
					dm = BKE_fracture_bmesh_to_mesh(bm);
#if 0
					{
						/*create debug mesh*/
						Object* o;
						o = BKE_object_add(G.main, fmd->modifier.scene, OB_MESH, "DUMMY");
						BM_mesh_bm_to_me(bm, o->data, (&(struct BMeshToMeshParams){0}));
					}
#endif

					BM_mesh_free(bm);

					/*do intersection*/
					//TODO FIX
					//intersect_shards_by_dm(fmd, dm, obj, NULL, inner_material_index, mat, true, fmd->boolean_double_threshold);

					BKE_fracture_mesh_free(dm);
				}
			}
		}
	}
}

#if 0 //TODO FIX
void BKE_fracture_shard_by_planes(FractureModifierData *fmd, Object *obj, short inner_material_index, float mat[4][4])
{
	if (fmd->cutter_group != NULL && obj->type == OB_MESH)
	{
		CollectionObject* go;
		float imat[4][4];

		//reset_shards(fmd);
		invert_m4_m4(imat, obj->obmat);

		for (go = fmd->cutter_group->gobject.first; go; go = go->next)
		{
			Object* ob = go->ob;

			printf("Cutting with %s ...\n", ob->id.name);
			/*simple case....one cutter object per object*/
			if (ob->type == OB_MESH) {

				FractureModifierData *fmd2 = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);
				if (fmd2 && BLI_listbase_count(&fmd2->shared->mesh_islands) > 0)
				{
					MeshIsland* mi = NULL;
					int j = 0;
					for (mi = fmd2->shared->mesh_islands.first; mi; mi = mi->next)
					{
						Mesh *dm;
						MVert *mv, *v = NULL;

						//TODO FIX
					//	dm = BKE_fracture_mesh_copy(mi->physics_mesh, ob);
						mv = dm->mvert;
						int totvert = dm->totvert;
						int i = 0;

						//printf("Cutting with %s, island %d...\n", ob->id.name, j);
						for (i = 0, v = mv; i < totvert; i++, v++)
						{
							add_v3_v3(v->co, mi->centroid);
						}

					//	TODO FIX
					//	intersect_shards_by_dm(fmd, dm, obj, ob, inner_material_index, mat, false, fmd->boolean_double_threshold);

						BKE_mesh_free(dm);
						dm = NULL;
						j++;
					}

					/*now delete first shards, those are the old ones */
					while (fmd->shared->frac_mesh->shard_count > 0)
					{
						/*clean up old entries here to avoid unnecessary shards*/
						Shard *first = fmd->shared->frac_mesh->shard_map.first;
						BLI_remlink_safe(&fmd->shared->frac_mesh->shard_map,first);
						BKE_fracture_shard_free(first, true);
						first = NULL;
						fmd->shared->frac_mesh->shard_count--;
					}

					/* re-synchronize counts, was possibly different before */
					fmd->shared->frac_mesh->shard_count = BLI_listbase_count(&fmd->shared->frac_mesh->shard_map);
				}
				else
				{

					Mesh *d = ob->runtime.mesh_eval;
					if (d == NULL) {
						d = ob->data;
					}

					intersect_shards_by_dm(fmd, d, obj, ob, inner_material_index, mat, true, fmd->boolean_double_threshold);
				}
			}
		}
	}
}
#endif

void BKE_fracture_shard_by_points(FractureModifierData *fmd, FracPointCloud *pointcloud, Object *ob, MeshIsland* mi,
								 Scene *scene, Main *bmain)
{
	int n_size = 8;

	float min[3], max[3];
	float theta = 0.001f; /* TODO, container enlargement, because boundbox exact container and boolean might create artifacts */
	int p;

	container *voro_container;
	particle_order *voro_particle_order;
	cell *voro_cells;

#ifdef USE_DEBUG_TIMER
	double time_start;
#endif

	printf("Fracturing with %d points...\n", pointcloud->totpoints);
	/* calculate bounding box with theta margin */
	INIT_MINMAX(min, max);
	BKE_mesh_minmax(mi->mesh, min, max);

	add_v3_fl(min, -theta);
	add_v3_fl(max, theta);

	mul_m4_v3(fmd->shared->splinter_matrix, min);
	mul_m4_v3(fmd->shared->splinter_matrix, max);

	if (fmd->point_source & MOD_FRACTURE_GRID)
	{
		float off[3] =  {0, 0, 0};
		float eps[3] = {theta, theta, theta};

		for (p = 0; p < pointcloud->totpoints; p++)
		{
			//find "max" offset (where atleast 1 axis is > 0)
			if (pointcloud->points[p].offset[0] > 0 ||
				pointcloud->points[p].offset[1] > 0 ||
				pointcloud->points[p].offset[2] > 0)
			{
				copy_v3_v3(off, pointcloud->points[p].offset);
				break;
			}
		}

		if (off[0] > 0 || off[1] > 0 || off[2] > 0)
		{
			sub_v3_v3(min, off);
		}

		sub_v3_v3(min, eps);

		//special treatment for grid pointsource... with offsets
		voro_container = container_new(min[0], max[0], min[1], max[1], min[2], max[2],
										n_size, n_size, n_size, false, false, false,
										pointcloud->totpoints);
	}
	else {
		voro_container = container_new(min[0], max[0], min[1], max[1], min[2], max[2],
									   n_size, n_size, n_size, false, false, false,
									   pointcloud->totpoints);
	}

	voro_particle_order = particle_order_new();
	for (p = 0; p < pointcloud->totpoints; p++) {
		float *co = pointcloud->points[p].co;
		container_put(voro_container, voro_particle_order, p, co[0], co[1], co[2]);
	}

#ifdef USE_DEBUG_TIMER
	time_start = PIL_check_seconds_timer();
#endif

	/* we expect as many raw cells as we have particles */
	voro_cells = cells_new(pointcloud->totpoints);

	/*Compute directly...*/
	container_compute_cells(voro_container, voro_particle_order, voro_cells);

	/*Apply offsets (if any, grid only */
	if (fmd->point_source & MOD_FRACTURE_GRID)
	{
		int v = 0;
		float fact[3] = {1 + fmd->grid_spacing[0], 1 + fmd->grid_spacing[1], 1 + fmd->grid_spacing[2]};
		for (p = 0; p < pointcloud->totpoints; p++)
		{
			//adjust centroid and...
			float off[3], cent[3];

			copy_v3_v3(off, pointcloud->points[p].offset);
			add_v3_v3(voro_cells[p].centroid, off);
			copy_v3_v3(cent, voro_cells[p].centroid);
			mul_v3_v3(cent, fact);

			//vertex coordinates
			for (v = 0; v < voro_cells[p].totvert; v++)
			{
				add_v3_v3(voro_cells[p].verts[v], off);

				//print_v3("Vert", voro_cells[p].verts[v]);
				add_v3_v3(voro_cells[p].verts[v], cent);
				sub_v3_v3(voro_cells[p].verts[v], voro_cells[p].centroid);
			}
		}
	}

	/*Evaluate result*/
	process_cells(fmd, mi, bmain, ob, scene, voro_cells, pointcloud->totpoints);


	/*Free structs in C++ area of memory */
	cells_free(voro_cells, pointcloud->totpoints);
	particle_order_free(voro_particle_order);
	container_free(voro_container);

#ifdef USE_DEBUG_TIMER
	printf("Fracture done, %g\n", PIL_check_seconds_timer() - time_start);
#endif

}

static void do_marking(FractureModifierData *fmd, Mesh *result)
{
	MEdge *medge = result->medge;
	MPoly *mpoly = result->mpoly, *mp = NULL;
	MLoop *mloop = result->mloop;
	MVert *mvert = result->mvert;
	int totpoly = result->totpoly;
	int i = 0;
	for (i = 0, mp = mpoly; i < totpoly; i++, mp++)
	{
		if (mp->flag & ME_FACE_SEL)
		{
			int j = 0;
			for (j = 0; j < mp->totloop; j++)
			{
				MLoop ml;
				ml = mloop[mp->loopstart + j];
				medge[ml.e].flag |= ME_SHARP;
				medge[ml.e].crease = fmd->inner_crease * 255.0f;
				mvert[ml.v].flag |= ME_VERT_TMP_TAG;
			}

			if (fmd->use_smooth)
				mp->flag |= ME_SMOOTH;
		}
		else
		{
			/*remove verts from unselected faces again*/
			int j = 0;
			for (j = 0; j < mp->totloop; j++)
			{
				MLoop ml;
				ml = mloop[mp->loopstart + j];
				mvert[ml.v].flag &= ~ME_VERT_TMP_TAG;
			}
		}
	}
}

void fracture_copy_customdata(CustomData* src, CustomData* dst,CustomDataMask mask, int src_ofs, int dst_ofs,
                              int copyelem, int totelem)
{
	CustomDataLayer *layer;
	int i;
	for (i = 0; i < src->totlayer; i++)
	{
		layer = src->layers + i;
		if (mask & CD_TYPE_AS_MASK(layer->type))
		{
			if (!CustomData_has_layer(dst, layer->type))
			{
				CustomData_add_layer(dst, layer->type, CD_CALLOC, NULL, totelem);
			}

			CustomData_copy_data_layer(src, dst, i, CustomData_get_layer_index(dst, layer->type), src_ofs, dst_ofs, copyelem);
		}
	}
}


Mesh* BKE_fracture_assemble_mesh_from_islands(FractureModifierData* fmd, ListBase *islands, Object* ob, float ctime)
{
	float imat[4][4];
	MeshIsland *mi;
	Mesh *mesh = NULL;
	int vertstart, polystart, loopstart, edgestart, num_verts, num_polys, num_loops, num_edges;
	vertstart = polystart = loopstart = edgestart = num_verts = num_polys = num_loops = num_edges = 0;
	const CustomDataMask CD_MASK_ISLAND =
		    CD_MASK_MDEFORMVERT |
		    CD_MASK_PROP_FLT | CD_MASK_PROP_INT | CD_MASK_PROP_STR | CD_MASK_MDISPS |
		    CD_MASK_MLOOPUV | CD_MASK_MLOOPCOL |
		    CD_MASK_RECAST | CD_MASK_PAINT_MASK |
		    CD_MASK_GRID_PAINT_MASK | CD_MASK_MVERT_SKIN | CD_MASK_FREESTYLE_EDGE | CD_MASK_FREESTYLE_FACE |
		    CD_MASK_CUSTOMLOOPNORMAL | CD_MASK_FACEMAP;

	for (mi = islands->first; mi; mi = mi->next)
	{
		RigidBodyOb *rbo = mi->rigidbody;

		if (BKE_fracture_meshisland_check_frame(fmd, mi, (int)ctime)) {
 			continue;
		}

		num_verts += mi->mesh->totvert;
		num_polys += mi->mesh->totpoly;
		num_loops += mi->mesh->totloop;
		num_edges += mi->mesh->totedge;
	}

	mesh = BKE_mesh_new_nomain(num_verts, num_edges, 0, num_loops, num_polys);
	mi = islands->first;

	if (fmd->shared->vert_index_map) {
		BLI_ghash_free(fmd->shared->vert_index_map, NULL, NULL);
		fmd->shared->vert_index_map = NULL;
	}

	fmd->shared->vert_index_map = BLI_ghash_int_new("vert_index_map");

	invert_m4_m4(imat, ob->obmat);
	for (mi = islands->first; mi; mi = mi->next)
	{
		MVert *mv;
		MPoly *mp;
		MLoop *ml;
		MEdge *me;
		int i, v;
		float iquat[4], irot[4], quat[4], size[3];
		RigidBodyOb *rbo = mi->rigidbody;

		if (BKE_fracture_meshisland_check_frame(fmd, mi, (int)ctime)) {
			continue;
		}

		memcpy(mesh->mvert + vertstart, mi->mesh->mvert, mi->mesh->totvert * sizeof(MVert));

		invert_qt_qt(irot, mi->rot);

		mat4_to_quat(quat, ob->obmat);
		invert_qt_qt(iquat, quat);
		mat4_to_size(size, ob->obmat);

		/*transform meshisland meshes, perform calculation here */
		for (v = 0, mv = mesh->mvert + vertstart; v < mi->mesh->totvert; v++, mv++)
		{
			float fno[3], centr[3];

			if (mi->rigidbody)
			{
				if (fmd->fix_normals) {
					/*ignore global quaternion rotation here */
					normal_short_to_float_v3(fno, mi->mesh->mvert[v].no);
					mul_qt_v3(mi->rigidbody->orn, fno);
					mul_qt_v3(iquat, fno);
					normal_float_to_short_v3(mv->no, fno);
				}

				mul_v3_v3(mv->co, size);
				mul_qt_v3(mi->rigidbody->orn, mv->co);
				copy_v3_v3(centr, mi->centroid);
				mul_v3_v3(centr, size);
				mul_qt_v3(mi->rigidbody->orn, centr);
				sub_v3_v3(mv->co, centr);
				add_v3_v3(mv->co, mi->rigidbody->pos);
				mul_m4_v3(imat, mv->co);
			}

			BLI_ghash_insert(fmd->shared->vert_index_map, SET_INT_IN_POINTER(vertstart + v), SET_INT_IN_POINTER(mi->id));
		}

		memcpy(mesh->mpoly + polystart, mi->mesh->mpoly, mi->mesh->totpoly * sizeof(MPoly));

		for (i = 0, mp = mesh->mpoly + polystart; i < mi->mesh->totpoly; ++i, ++mp) {
			/* adjust loopstart index */
			mp->loopstart += loopstart;

			/* material index lookup and correction, avoid having the same material in different slots */
			//index = GET_INT_FROM_POINTER(BLI_ghash_lookup(mat_index_map, SET_INT_IN_POINTER(mp->mat_nr + fmd->shared->matstart)));
			//mp->mat_nr = index-1;
		}

		memcpy(mesh->mloop + loopstart, mi->mesh->mloop, mi->mesh->totloop * sizeof(MLoop));

		for (i = 0, ml = mesh->mloop + loopstart; i < mi->mesh->totloop; ++i, ++ml) {
			/* adjust vertex index */
			ml->v += vertstart;
			ml->e += edgestart;
		}

		memcpy(mesh->medge + edgestart, mi->mesh->medge, mi->mesh->totedge * sizeof(MEdge));

		for (i = 0, me = mesh->medge + edgestart; i < mi->mesh->totedge; ++i, ++me) {
			/* adjust vertex indices */
			me->v1 += vertstart;
			me->v2 += vertstart;
		}

		fracture_copy_customdata(&mi->mesh->vdata, &mesh->vdata, CD_MASK_ISLAND, 0, vertstart, mi->mesh->totvert, num_verts);
		fracture_copy_customdata(&mi->mesh->edata, &mesh->edata, CD_MASK_ISLAND, 0, edgestart, mi->mesh->totedge, num_edges);
		fracture_copy_customdata(&mi->mesh->ldata, &mesh->ldata, CD_MASK_ISLAND, 0, loopstart, mi->mesh->totloop, num_loops);
		fracture_copy_customdata(&mi->mesh->pdata, &mesh->pdata, CD_MASK_ISLAND, 0, polystart, mi->mesh->totpoly, num_polys);

		vertstart += mi->mesh->totvert;
		polystart += mi->mesh->totpoly;
		loopstart += mi->mesh->totloop;
		edgestart += mi->mesh->totedge;
	}

	do_marking(fmd, mesh);

	BKE_mesh_calc_normals(mesh);
	return mesh;
}

#if 0
void BKE_match_vertex_coords(MeshIsland* mi, MeshIsland *par, Object *ob, int frame, bool is_parent, bool shards_to_islands)
{
	float loc[3] = {0.0f, 0.0f, 0.0f};
	float rot[4] = {1.0f, 0.0f, 0.0f, 0.0f};
	int j = 0;

	float centr[3] = {0.0f, 0.0f, 0.0f};

	float mat[4][4];
	float quat[4] = {1.0f, 0.0f, 0.0f, 0.0f};
	float qrot[4] = {1.0f, 0.0f, 0.0f, 0.0f};
	float iquat[4] =  {1.0f, 0.0f, 0.0f, 0.0f};

	int val = shards_to_islands ? -1 : 0;

	invert_m4_m4(mat, ob->obmat);

	mi->locs[0] = loc[0] = par->locs[3*frame];
	mi->locs[1] = loc[1] = par->locs[3*frame+1];
	mi->locs[2] = loc[2] = par->locs[3*frame+2];

	mi->rots[0] = rot[0] = par->rots[4*frame];
	mi->rots[1] = rot[1] = par->rots[4*frame+1];
	mi->rots[2] = rot[2] = par->rots[4*frame+2];
	mi->rots[3] = rot[3] = par->rots[4*frame+3];

	mul_m4_v3(mat, loc);
	mat4_to_quat(quat, ob->obmat);
	invert_qt_qt(iquat, quat);

	if (par->id == val)
	{
		invert_qt_qt(qrot, par->rot);
		mul_qt_qtqt(qrot, rot, qrot);
		mul_qt_qtqt(qrot, iquat, qrot);
	}
	else
	{
		mul_qt_qtqt(qrot, rot, par->rot);
		mul_qt_qtqt(qrot, iquat, qrot);
	}

	if (is_parent)
	{
		copy_v3_v3(centr, mi->centroid);
		mul_qt_v3(qrot, centr);
		add_v3_v3(centr, loc);
	}
	else
	{
		copy_v3_v3(centr, loc);
	}

	for (j = 0; j < mi->vertex_count; j++)
	{
		float co[3];

		//first add vert to centroid, then rotate
		copy_v3_v3(co, mi->vertices_cached[j]->co);

		sub_v3_v3(co, mi->centroid);
		mul_qt_v3(qrot, co);
		add_v3_v3(co, centr);

		copy_v3_v3(mi->vertices_cached[j]->co, co);

		mi->vertco[3*j]   = co[0];
		mi->vertco[3*j+1] = co[1];
		mi->vertco[3*j+2] = co[2];
	}

	{
		Mesh *dm = mi->physics_mesh;
		MVert* mv, *mvert = dm->mvert;
		int numVert = dm->totvert;

		for (j = 0, mv = mvert; j < numVert; j++, mv++)
		{
			//also rotate physicsmesh (shouldnt be necessary,
			//but lets do it to check whether its correct then)
			mul_qt_v3(qrot, mv->co);
		}
	}

	//init rigidbody properly ?
	copy_v3_v3(mi->centroid, centr);
	copy_qt_qt(mi->rot, qrot);
}
#endif

/* flush a hflag to from verts to edges/faces */
void BKE_bm_mesh_hflag_flush_vert(BMesh *bm, const char hflag)
{
	BMEdge *e;
	BMLoop *l_iter;
	BMLoop *l_first;
	BMFace *f;

	BMIter eiter;
	BMIter fiter;

	int ok;

	BM_ITER_MESH (e, &eiter, bm, BM_EDGES_OF_MESH) {
		if (BM_elem_flag_test(e->v1, hflag) &&
			BM_elem_flag_test(e->v2, hflag))
		{
			BM_elem_flag_enable(e, hflag);
		}
		else {
			BM_elem_flag_disable(e, hflag);
		}
	}
	BM_ITER_MESH (f, &fiter, bm, BM_FACES_OF_MESH) {
		ok = true;
		l_iter = l_first = BM_FACE_FIRST_LOOP(f);
		do {
			if (!BM_elem_flag_test(l_iter->v, hflag)) {
				ok = false;
				break;
			}
		} while ((l_iter = l_iter->next) != l_first);

		BM_elem_flag_set(f, hflag, ok);
	}
}

void BKE_update_velocity_layer(FractureModifierData *fmd, Mesh *dm)
{
	float *velX=NULL, *velY=NULL, *velZ = NULL;
	RigidBodyOb *rbo = NULL;
	int i = 0;
	MeshIsland *mi;
	int totvert;

	if (!dm)
		return;

	totvert = dm->totvert;

	velX = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "velX");
	velY = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "velY");
	velZ = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "velZ");

	if (!velX)
		velX = CustomData_add_layer_named(&dm->vdata, CD_PROP_FLT, CD_CALLOC, NULL, totvert, "velX");

	if (!velY)
		velY = CustomData_add_layer_named(&dm->vdata, CD_PROP_FLT, CD_CALLOC, NULL, totvert, "velY");

	if (!velZ)
		velZ = CustomData_add_layer_named(&dm->vdata, CD_PROP_FLT, CD_CALLOC, NULL, totvert, "velZ");

	for (i = 0; i < totvert; i++)
	{
		mi = BLI_ghash_lookup(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(i));
		rbo = mi->rigidbody;

		velX[i] = rbo->lin_vel[0] + rbo->ang_vel[0];
		velY[i] = rbo->lin_vel[1] + rbo->ang_vel[1];
		velZ[i] = rbo->lin_vel[2] + rbo->ang_vel[2];
	}
}

static void fracture_anim_bind_activate(MeshIsland *mi, AnimBind *bind)
{
	bind->v = -1;
	bind->v1 = -1;
	bind->v2 = -1;
	bind->mi = -1;
	zero_v3(bind->no);
	zero_v3(bind->offset);
	unit_qt(bind->quat);

	if (mi->rigidbody->type == RBO_TYPE_ACTIVE)
	{
		RigidBodyOb* rbo = mi->rigidbody;

		rbo->flag &= ~RBO_FLAG_KINEMATIC;
		rbo->flag &= ~RBO_FLAG_KINEMATIC_BOUND;
		rbo->flag |= RBO_FLAG_NEEDS_VALIDATE;

		if (rbo->shared->physics_object)
		{
			RB_body_set_mass(rbo->shared->physics_object, rbo->mass);
			RB_body_set_kinematic_state(rbo->shared->physics_object, false);
			RB_body_activate(rbo->shared->physics_object);
		}
	}
}

void BKE_fracture_animated_loc_rot(FractureModifierData *fmd, Object *ob, bool do_bind, Depsgraph *depsgraph)
{
	//to be called after rigidbodies have been actually created... from MOD_fracture.c
	//rotation is optional, remesher + particlesystem can provide it
	float *quatX, *quatY, *quatZ, *quatW;
	MVert *mvert = NULL;
	MPoly *mpoly = NULL;
	MLoop *mloop = NULL;
	MeshIsland *mi;
	Mesh *dm = NULL;
	int totvert, count = 0, i = 0, *orig_index = NULL, totpoly, items;
	KDTree *tree = NULL;
	float anim_imat[4][4], imat[4][4];
	Object *ob_eval;
	bool mesh_free = false;

	if (!fmd->anim_mesh_ob)
		return;

	if (fmd->anim_mesh_ob == ob)
		return;

	ob_eval = DEG_get_evaluated_object(depsgraph, fmd->anim_mesh_ob);
	dm = BKE_modifier_get_evaluated_mesh_from_evaluated_object(ob_eval, &mesh_free);

	if (!dm)
		return;

	totvert = dm->totvert;
	totpoly = dm->totpoly;

	invert_m4_m4(anim_imat, fmd->anim_mesh_ob->obmat);
	invert_m4_m4(imat, ob->obmat);

	if (do_bind) {

		items = totpoly > 0 ? totpoly : totvert;

		count = BLI_listbase_count(&fmd->shared->mesh_islands);
		tree = BLI_kdtree_new(items);

		fmd->shared->anim_bind_len = count;
		if (fmd->shared->anim_bind) {
			MEM_freeN(fmd->shared->anim_bind);
		}

		//TODO, possibly a weak solution, but do we really want to store the length too ?
		fmd->shared->anim_bind = MEM_mallocN(sizeof(AnimBind) * fmd->shared->anim_bind_len, "anim_bind");
		for (i = 0; i < fmd->shared->anim_bind_len; i++)
		{
			fmd->shared->anim_bind[i].mi = -1;
			fmd->shared->anim_bind[i].v = -1;
			fmd->shared->anim_bind[i].v1 = -1;
			fmd->shared->anim_bind[i].v2 = -1;
			zero_v3(fmd->shared->anim_bind[i].offset);
			zero_v3(fmd->shared->anim_bind[i].no);
			unit_qt(fmd->shared->anim_bind[i].quat);
		}
	}

	i = 0;
	mvert = dm->mvert;
	mpoly = dm->mpoly;
	mloop = dm->mloop;
	if (do_bind)
	{
		if (totpoly > 0)
		{
			//poly based bind
			for (i = 0; i < totpoly; i++)
			{
				float co[3];
				copy_v3_v3(co, mvert[mloop[mpoly[i].loopstart].v].co);
				BLI_kdtree_insert(tree, i, co);
			}
		}
		else
		{
			//no faces -> vertex based bind
			for (i = 0; i < totvert; i++)
			{
				float co[3];
				copy_v3_v3(co, mvert[i].co);
				BLI_kdtree_insert(tree, i, co);
			}
		}

		BLI_kdtree_balance(tree);
	}


	quatX = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "quatX");
	quatY = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "quatY");
	quatZ = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "quatZ");
	quatW = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "quatW");
	orig_index = CustomData_get_layer(&dm->vdata, CD_ORIGINDEX);

	//check vertexcount and islandcount, TODO for splitshards... there it might differ, ignore then for now
	//later do interpolation ? propagate to islands somehow then, not yet now...
	//also maybe skip for dynamic for now, since this is totally different

	//bind loop, bind the verts to the shards
	if (do_bind)
	{
		i = 0;
		for (mi = fmd->shared->mesh_islands.first; mi; mi = mi->next, i++)
		{
			KDTreeNearest n;
			float co[3], diff[3] = {0, 0, 0}, f_no[3];

			copy_v3_v3(co, mi->rigidbody->pos);
			mul_m4_v3(anim_imat, co);
			BLI_kdtree_find_nearest(tree, co, &n);

			if (n.dist <= fmd->anim_bind_limit || fmd->anim_bind_limit == 0)
			{
				if (totpoly > 0)
				{
					int v1, v2, v3;
					float limit = 0.0001f;
					MPoly *mp  = mpoly + n.index;
					MLoop *ml = mloop + mp->loopstart;
					BKE_mesh_calc_poly_normal(mp, ml, mvert, f_no);

					v1 = ml->v;
					v2 = (ml + 1)->v;
					v3 = (ml + 2)->v;

					if (mpoly[n.index].totloop < 3)
					{
						printf("Degenerate face, skipping\n");
						fracture_anim_bind_activate(mi, &fmd->shared->anim_bind[i]);
						continue;
					}

					if (compare_v3v3(mvert[v1].co, mvert[v2].co, limit) ||
						compare_v3v3(mvert[v1].co, mvert[v3].co, limit) ||
						compare_v3v3(mvert[v2].co, mvert[v3].co, limit))
					{
						printf("Very close coordinates, skipping %d %d %d in %d\n", v1, v2, v3, mi->id);
						fracture_anim_bind_activate(mi, &fmd->shared->anim_bind[i]);
						continue;
					}

					fmd->shared->anim_bind[i].v = v1;
					fmd->shared->anim_bind[i].v1 = v2;
					fmd->shared->anim_bind[i].v2 = v3;
					fmd->shared->anim_bind[i].poly = n.index;
					mi->rigidbody->flag |= RBO_FLAG_KINEMATIC_BOUND;
				}
				else {
					fmd->shared->anim_bind[i].v = n.index;
					fmd->shared->anim_bind[i].v1 = -1;
					fmd->shared->anim_bind[i].v2 = -1;
					mi->rigidbody->flag |= RBO_FLAG_KINEMATIC_BOUND;
				}

				if (fmd->shared->anim_bind[i].v != -1)
				{
					fmd->shared->anim_bind[i].mi = i;
					sub_v3_v3v3(diff, n.co, co);

					copy_v3_v3(fmd->shared->anim_bind[i].offset, diff);

					if ((fmd->shared->anim_bind[i].v1 == -1 || fmd->shared->anim_bind[i].v2 == -1)) {
						//fallback if not enough verts around
						normal_short_to_float_v3(fmd->shared->anim_bind[i].no, mvert[n.index].no);
						normalize_v3(fmd->shared->anim_bind[i].no);
					}
					else if (totpoly > 0) {
						tri_to_quat_ex(fmd->shared->anim_bind[i].quat,
								mvert[fmd->shared->anim_bind[i].v].co,
								mvert[fmd->shared->anim_bind[i].v1].co,
								mvert[fmd->shared->anim_bind[i].v2].co, f_no);
						copy_v3_v3(fmd->shared->anim_bind[i].no, f_no);
					}
				}
			}
			else
			{
				fracture_anim_bind_activate(mi, &fmd->shared->anim_bind[i]);
			}
		}

		if (tree)
			BLI_kdtree_free(tree);
	}
	else
	{
		mi = fmd->shared->mesh_islands.first;
		for (i = 0; i < fmd->shared->anim_bind_len; i++, mi = mi->next)
		{
			float co[3];
			int index = -1;
			int vindex = -1;

			index = fmd->shared->anim_bind[i].mi;
			vindex = fmd->shared->anim_bind[i].v;

			if (index == -1 || vindex == -1)
			{
				fracture_anim_bind_activate(mi, &fmd->shared->anim_bind[i]);
				continue;
			}

			//only let kinematic rbs do this, active ones are being taken care of by bullet
			if (mi && mi->rigidbody && (mi->rigidbody->flag & RBO_FLAG_KINEMATIC))
			{
				//the 4 rot layers *should* be aligned, caller needs to ensure !
				bool quats = quatX && quatY && quatZ && quatW;
				float quat[4], vec[3], no[3], off[3];
				int v = fmd->shared->anim_bind[i].v;
				unit_qt(quat);

				if (v >= totvert) {
					fracture_anim_bind_activate(mi, &fmd->shared->anim_bind[i]);
					continue;
				}

				if ((orig_index && orig_index[v] != v && (fmd->shared->anim_bind[i].v1 == -1 || fmd->shared->anim_bind[i].v2 == -1)))
				{
					fracture_anim_bind_activate(mi, &fmd->shared->anim_bind[i]);
					continue;
				}

				copy_v3_v3(co, mvert[v].co);
				copy_v3_v3(off, fmd->shared->anim_bind[i].offset);

				if (fmd->anim_mesh_rot)
				{
					if (quats)
					{
						quat[0] = quatX[v];
						quat[1] = quatY[v];
						quat[2] = quatZ[v];
						quat[3] = quatW[v];
					}
					else
					{
						copy_v3_v3(vec, fmd->shared->anim_bind[i].no);
						if (fmd->shared->anim_bind[i].v1 == -1 || fmd->shared->anim_bind[i].v2 == -1) {
							//fallback if not enough verts around;
							normal_short_to_float_v3(no, mvert[v].no);
							normalize_v3(no);
							rotation_between_vecs_to_quat(quat, vec, no);
						}
						else {
							float rot[4], iquat[4], fno[3];
							MPoly *mp = mpoly + fmd->shared->anim_bind[i].poly;
							MLoop *ml = mloop + mp->loopstart;
							BKE_mesh_calc_poly_normal(mp, ml, mvert, fno);

							tri_to_quat_ex(rot, mvert[fmd->shared->anim_bind[i].v].co,
											  mvert[fmd->shared->anim_bind[i].v1].co,
											  mvert[fmd->shared->anim_bind[i].v2].co,
											  fno);
							invert_qt_qt(iquat, fmd->shared->anim_bind[i].quat);
							mul_qt_qtqt(quat, rot, iquat);
						}
					}
				}

				mul_qt_v3(quat, off);
				sub_v3_v3(co, off);
				mul_m4_v3(fmd->anim_mesh_ob->obmat, co);

				copy_v3_v3(mi->rigidbody->pos, co);

				if (fmd->anim_mesh_rot)
				{
					float ob_quat[4];
					mat4_to_quat(ob_quat, ob->obmat);
					mul_qt_qtqt(quat, ob_quat, quat);
					copy_qt_qt(mi->rigidbody->orn, quat);
				}

				mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
				if (mi->rigidbody->shared->physics_object)
				{
					RB_body_set_loc_rot(mi->rigidbody->shared->physics_object, mi->rigidbody->pos, mi->rigidbody->orn);
				}
			}
		}
	}

	if (mesh_free) {
		BKE_fracture_mesh_free(dm);
		return;
	}
}

bool BKE_fracture_meshisland_check_frame(FractureModifierData *fmd, MeshIsland* mi, int frame)
{
	return ((frame < mi->startframe && mi->startframe > fmd->shared->last_cache_start) ||
		   (frame >= mi->endframe && mi->endframe < fmd->shared->last_cache_end));
}

void BKE_fracture_meshisland_check_realloc_cache(FractureModifierData *fmd, RigidBodyWorld *rbw, MeshIsland* mi, int frame)
{
	int endframe = rbw->shared->pointcache->endframe;
	int startframe = rbw->shared->pointcache->startframe;

	//only grow...
	if (endframe != fmd->shared->last_cache_end)
	{
		//keep invalid shards untouched
		if (!BKE_fracture_meshisland_check_frame(fmd, mi, frame) && mi->endframe == fmd->shared->last_cache_end)
		{
			mi->endframe = endframe;
			frame = mi->endframe - frame + 1;

			mi->locs = MEM_reallocN(mi->locs, sizeof(float*) * 3 * frame);
			mi->rots = MEM_reallocN(mi->rots, sizeof(float*) * 4 * frame);
			mi->vels = MEM_reallocN(mi->vels, sizeof(float*) * 3 * frame);
			mi->aves = MEM_reallocN(mi->aves, sizeof(float*) * 3 * frame);
		}
	}

	fmd->shared->last_cache_end = endframe;
	fmd->shared->last_cache_start = startframe;
}

void BKE_fracture_meshislands_free(FractureModifierData* fmd, Scene* scene)
{
	MeshIsland *mi;

	while (fmd->shared->mesh_islands.first) {
		mi = fmd->shared->mesh_islands.first;
		BLI_remlink_safe(&fmd->shared->mesh_islands, mi);
		BKE_fracture_mesh_island_free(mi, scene);
		mi = NULL;
	}

	fmd->shared->mesh_islands.first = NULL;
	fmd->shared->mesh_islands.last = NULL;
}

#if 0
void BKE_fracture_dynamic_free(FractureModifierData *fmd, Scene *scene)
{
	/* in dynamic mode we have to get rid of the entire Meshisland sequence */
	/* either at manual refresh or when removing the modifier */
	MeshIslandSequence *msq;

	while (fmd->shared->meshIsland_sequence.first) {
		msq = fmd->shared->meshIsland_sequence.first;
		BLI_remlink(&fmd->shared->meshIsland_sequence, msq);
		BKE_fracture_meshislands_free(fmd, &msq->meshIslands, scene);
		MEM_freeN(msq);
		msq = NULL;
	}

	fmd->shared->meshIsland_sequence.first = NULL;
	fmd->shared->meshIsland_sequence.last = NULL;

	fmd->shared->mesh_islands.first = NULL;
	fmd->shared->mesh_islands.last = NULL;

	fmd->shared->current_mi_entry = NULL;
}
#endif

void BKE_fracture_modifier_free(FractureModifierData *fmd, Scene *scene)
{
	BKE_fracture_constraints_free(fmd, scene);
	BKE_fracture_meshislands_free(fmd, scene);

	if (fmd->shared->material_index_map)
	{
		BLI_ghash_free(fmd->shared->material_index_map, NULL, NULL);
		fmd->shared->material_index_map = NULL;
		fmd->shared->matstart = 1;
	}

	if (fmd->shared->defgrp_index_map)
	{
		BLI_ghash_free(fmd->shared->defgrp_index_map, NULL, NULL);
		fmd->shared->defgrp_index_map = NULL;
		fmd->shared->defstart = 0;
	}

	if (fmd->shared->vertex_island_map) {
		BLI_ghash_free(fmd->shared->vertex_island_map, NULL, NULL);
		fmd->shared->vertex_island_map = NULL;
	}

	if (fmd->shared->nor_tree != NULL) {
		BLI_kdtree_free(fmd->shared->nor_tree);
		fmd->shared->nor_tree = NULL;
	}

	if (fmd->shared->face_pairs != NULL) {
		BLI_ghash_free(fmd->shared->face_pairs, NULL, NULL);
		fmd->shared->face_pairs = NULL;
	}

	BKE_fracture_shared_verts_free(&fmd->shared->shared_verts);

	if (fmd->shared->vert_index_map != NULL) {
		BLI_ghash_free(fmd->shared->vert_index_map, NULL, NULL);
		fmd->shared->vert_index_map = NULL;
	}

	if (fmd->shared->anim_bind)
		MEM_freeN(fmd->shared->anim_bind);

	if (fmd->shared->last_island_tree)
	{
		BLI_kdtree_free(fmd->shared->last_island_tree);
		fmd->shared->last_island_tree = NULL;
	}

	if (fmd->shared->last_islands)
	{
		int j;
		for (j = 0; j < fmd->shared->last_expected_islands; j++)
		{
			BKE_fracture_mesh_island_free(fmd->shared->last_islands[j], scene);
		}

		MEM_freeN(fmd->shared->last_islands);
		fmd->shared->last_islands = NULL;

		fmd->shared->last_expected_islands = 0;
	}
}

//XXXX TODO same applies for autohide prep and normals fixing, latter could be a separate operator or so, called from refresh op
static KDTree *build_nor_tree(Mesh *dm)
{
	int i = 0, totvert = dm->totvert;
	KDTree *tree = BLI_kdtree_new(totvert);
	MVert *mv, *mvert = dm->mvert;

	for (i = 0, mv = mvert; i < totvert; i++, mv++) {
		BLI_kdtree_insert(tree, i, mv->co);
	}

	BLI_kdtree_balance(tree);

	return tree;
}


//XXX TODO plan is to get rid of this, since we have a packing mechanism now, but wrap that functionality to useful op (via C, python api is optional)
static int getGroupObjects(Collection *gr, Object ***obs, int g_exist)
{
	int ctr = g_exist;
	CollectionObject *go;
	if (gr == NULL) return ctr;

	for (go = gr->gobject.first; go; go = go->next) {

		*obs = MEM_reallocN(*obs, sizeof(Object *) * (ctr + 1));
		(*obs)[ctr] = go->ob;
		ctr++;
	}

	return ctr;
}

#if 0
static Mesh* get_object_dm(Object* o)
{
	Mesh *dm_ob = NULL;

	/*ensure o->derivedFinal*/
	FractureModifierData* fmd2 = (FractureModifierData*) modifiers_findByType(o, eModifierType_Fracture);
	if (fmd2)
	{
		dm_ob = fmd2->shared->visible_mesh_cached;
	}
	else
	{
		dm_ob = o->runtime.mesh_eval;
	}

	return dm_ob;
}
#endif

#if 0

static void adjustPolys(MPoly **mpoly, Mesh *dm_ob, GHash *mat_index_map, short matstart, int loopstart)
{
	MPoly *mp;
	int j;

	for (j = 0, mp = *mpoly; j < dm_ob->totpoly; ++j, ++mp) {
		short index = 0;
		/* adjust loopstart index */
		mp->loopstart += loopstart;

		/* material index lookup and correction, avoid having the same material in different slots */
		index = GET_INT_FROM_POINTER(BLI_ghash_lookup(mat_index_map, SET_INT_IN_POINTER(mp->mat_nr + matstart)));
		mp->mat_nr = index-1;
	}
}

static void adjustLoops(MLoop **mloop, Mesh *dm_ob, int vertstart, int loopstart, Mesh *result)
{
	MLoop *ml;
	int j;

	for (j = 0, ml = *mloop; j < dm_ob->totloop; ++j, ++ml) {
		/* adjust vertex index */
		if (CustomData_has_layer(&dm_ob->ldata, CD_MLOOPUV))
		{
			MLoopUV *mluv = CustomData_get(&dm_ob->ldata, j, CD_MLOOPUV);
			if (mluv)
				CustomData_set(&result->ldata, loopstart + j, CD_MLOOPUV, mluv);
		}
		ml->v += vertstart;
	}
}

static void adjustVerts(MVert **mvert, FractureModifierData *fmd, Object *o, Mesh* dm_ob, int vertstart, int i, Mesh* result)
{
	MVert *mv;
	int v;

	for (v = 0, mv = *mvert; v < dm_ob->totvert; v++, mv++)
	{
		if (CustomData_has_layer(&dm_ob->vdata, CD_MDEFORMVERT))
		{
			MDeformVert *mdv = CustomData_get(&dm_ob->vdata, v, CD_MDEFORMVERT);
			if (mdv)
				CustomData_set(&result->vdata, vertstart + v, CD_MDEFORMVERT, mdv);
		}
		mul_m4_v3(o->obmat, mv->co);
		BLI_ghash_insert(fmd->shared->vert_index_map, SET_INT_IN_POINTER(vertstart + v), SET_INT_IN_POINTER(i));
	}
}
#endif

#if 0
static void collect_derivedmeshes(Main* bmain, FractureModifierData* fmd, Object *ob, MVert** mvert, MLoop** mloop, MPoly **mpoly,
								  Mesh* result, GHash** mat_index_map)
{
	int vertstart = 0, polystart = 0, loopstart = 0;
	int matstart = 1;
	MVert *mverts = *mvert;
	MLoop *mloops = *mloop;
	MPoly *mpolys = *mpoly;

	MVert *mv;
	MLoop *ml;
	MPoly *mp;

	CollectionObject* go;
	int totcol;
	int i = 0;

	for (go = fmd->dm_group->gobject.first; go; go = go->next)
	{
		Mesh* dm_ob = NULL;
		Object *o = go->ob;

		dm_ob = get_object_dm(o);
		if (dm_ob == NULL)
		{   /* avoid crash atleast...*/
			return;
		}

		totcol = BKE_fracture_collect_materials(bmain, o, ob, matstart, mat_index_map);

		mv = mverts + vertstart;
		memcpy(mv, dm_ob->mvert, dm_ob->totvert * sizeof(MVert));
		adjustVerts(&mv, fmd, o, dm_ob, vertstart, i, result);

		mp = mpolys + polystart;
		memcpy(mp, dm_ob->mpoly, dm_ob->totpoly * sizeof(MPoly));
		adjustPolys(&mp, dm_ob, *mat_index_map, matstart, loopstart);

		ml = mloops + loopstart;
		memcpy(ml, dm_ob->mloop, dm_ob->totloop * sizeof(MLoop));
		adjustLoops(&ml, dm_ob, vertstart, loopstart, result);

		vertstart += dm_ob->totvert;
		polystart += dm_ob->totpoly;
		loopstart += dm_ob->totloop;
		matstart += totcol;
		i++;
	}
}
#endif
#if 0
static void count_dm_contents(FractureModifierData *fmd, int *num_verts, int *num_loops, int *num_polys)
{
	CollectionObject* go;

	for (go = fmd->dm_group->gobject.first; go; go = go->next)
	{
		Mesh* dm_ob = NULL;
		Object *o = go->ob;

		/*ensure o->derivedFinal*/
		FractureModifierData* fmd2 = (FractureModifierData*) modifiers_findByType(o, eModifierType_Fracture);
		if (fmd2)
		{
			dm_ob = fmd2->shared->visible_mesh_cached;
		}
		else
		{
			dm_ob = o->runtime.mesh_eval;
		}

		if (dm_ob == NULL) continue;

		(*num_verts) += dm_ob->totvert;
		(*num_polys) += dm_ob->totpoly;
		(*num_loops) += dm_ob->totloop;
	}
}


static Mesh *get_group_dm(Main* bmain, FractureModifierData *fmd, Mesh *dm, Object* ob, bool do_refresh)
{
	/* combine derived meshes from group objects into 1, trigger submodifiers if ob->derivedFinal is empty */
	int num_verts = 0, num_polys = 0, num_loops = 0;
	Mesh *result;
	MVert *mverts;
	MPoly *mpolys;
	MLoop *mloops;

	GHash *mat_index_map = NULL;

	if (fmd->dm_group && do_refresh && !fmd->use_constraint_group)
	{
		mat_index_map = BLI_ghash_int_new("mat_index_map");
		if (fmd->shared->vert_index_map != NULL) {
			BLI_ghash_free(fmd->shared->vert_index_map, NULL, NULL);
			fmd->shared->vert_index_map = NULL;
		}

		fmd->shared->vert_index_map = BLI_ghash_int_new("vert_index_map");

		count_dm_contents(fmd, &num_verts, &num_loops, &num_polys);
		if (num_verts == 0)
		{
			return dm;
		}

		result = BKE_mesh_new_nomain(num_verts, 0, 0, num_loops, num_polys);
		mverts = result->mvert;
		mloops = result->mloop;
		mpolys = result->mpoly;

		CustomData_add_layer(&result->vdata, CD_MDEFORMVERT, CD_CALLOC, NULL, num_verts);
		CustomData_add_layer(&result->ldata, CD_MLOOPUV, CD_CALLOC, NULL, num_loops);

		collect_derivedmeshes(bmain, fmd, ob, &mverts, &mloops, &mpolys, result, &mat_index_map);
		BKE_mesh_calc_edges(result, true, true);
		BKE_mesh_calc_normals(result);

		BLI_ghash_free(mat_index_map, NULL, NULL);
		mat_index_map = NULL;
		return result;
	}

	return dm;
}
#endif

static bool in_bbox(float p[3], float min[3], float max[3])
{
	return (p[0] > min[0]) && (p[0] < max[0]) && (p[1] > min[1]) && (p[1] < max[1]) && (p[2] > min[2]) && (p[2] < max[2]);
}

static void points_from_verts(Object **ob, int totobj, FracPointCloud *points, float mat[4][4], float thresh,
							  FractureModifierData *emd, Object *obj, MeshIsland *mi)
{
	int v, o, pt = points->totpoints;
	float co[3];

	for (o = 0; o < totobj; o++) {
		if (ob[o]->type == OB_MESH) {
			/* works for mesh objects only, curves, surfaces, texts have no verts */
			float imat[4][4];
			Mesh *d;
			MVert *vert;

			if (ob[o] == obj) {
				/* same object, use given derivedmesh */
				d = mi->mesh;
			}
			else {
				d = ob[o]->runtime.mesh_eval;
			}

			invert_m4_m4(imat, mat);
			vert = d->mvert;

			for (v = 0; v < d->totvert; v++) {
				if (BLI_thread_frand(0) < thresh) {

					copy_v3_v3(co, vert[v].co);

					//XXXX TODO, own verts transform seems to have a bug here as well
					if (emd->point_source & MOD_FRACTURE_EXTRA_VERTS) {
						mul_m4_v3(ob[o]->obmat, co);
						mul_m4_v3(imat, co);
					}

					if (mi->id > 0)
					{
						float min[3], max[3], cent[3];
						if (mi)
						{
							//arrange_shard(emd, mi->id, false, cent);

							//TODO FIX fill / cache mi->min / max
							add_v3_v3v3(min, mi->min, cent);
							add_v3_v3v3(max, mi->max, cent);
							if (in_bbox(co, min, max))
							{
								points->points = MEM_reallocN(points->points, (pt + 1) * sizeof(FracPoint));
								copy_v3_v3(points->points[pt].co, co);
								zero_v3(points->points[pt].offset);
								pt++;
							}
						}
					}
					else
					{
						points->points = MEM_reallocN(points->points, (pt + 1) * sizeof(FracPoint));
						copy_v3_v3(points->points[pt].co, co);
						zero_v3(points->points[pt].offset);
						pt++;
					}
				}
			}
		}
	}

	points->totpoints = pt;
}

static void points_from_particles(Object **ob, int totobj, Scene *scene, FracPointCloud *points, float mat[4][4],
								  float thresh, FractureModifierData *fmd, MeshIsland *mi)
{
	int o, p, pt = points->totpoints;
	ParticleSystemModifierData *psmd;
	ParticleData *pa;
	ParticleSimulationData sim = {NULL};
	ParticleKey birth;
	ModifierData *mod;

	for (o = 0; o < totobj; o++) {
		for (mod = ob[o]->modifiers.first; mod; mod = mod->next) {
			if (mod->type == eModifierType_ParticleSystem) {
				float imat[4][4];
				psmd = (ParticleSystemModifierData *)mod;
				sim.scene = scene;
				sim.ob = ob[o];
				sim.psys = psmd->psys;
				sim.psmd = psmd;
				invert_m4_m4(imat, mat);

				for (p = 0, pa = psmd->psys->particles; p < psmd->psys->totpart; p++, pa++) {
					/* XXX was previously there to choose a particle with a certain state */
					bool particle_unborn = pa->alive == PARS_UNBORN;
					bool particle_alive = pa->alive == PARS_ALIVE;
					bool particle_dead = pa->alive == PARS_DEAD;
					bool particle_mask = particle_unborn || particle_alive || particle_dead;

					if ((BLI_thread_frand(0) < thresh) && particle_mask) {
						float co[3];

						/* birth coordinates are not sufficient in case we did pre-simulate the particles, so they are not
						 * aligned with the emitter any more BUT as the particle cache is messy and shows initially wrong
						 * positions "sabotaging" fracture, default use case is using birth coordinates, let user decide... */
						if (fmd->use_particle_birth_coordinates && fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED)
						{
							psys_get_birth_coords(&sim, pa, &birth, 0, 0);
						}
						else {
							psys_get_particle_state(&sim, p, &birth, 1);
						}

						copy_v3_v3(co, birth.co);
						mul_m4_v3(imat, co);

						if (mi->id > 0)
						{
							float min[3], max[3], cent[3];
//							arrange_shard(fmd, mi->id, false, cent);
							add_v3_v3v3(min, mi->min, cent);
							add_v3_v3v3(max, mi->max, cent);
							if (in_bbox(co, min, max))
							{
								points->points = MEM_reallocN(points->points, (pt + 1) * sizeof(FracPoint));
								copy_v3_v3(points->points[pt].co, co);
								zero_v3(points->points[pt].offset);
								pt++;
							}
						}
						else
						{
							points->points = MEM_reallocN(points->points, (pt + 1) * sizeof(FracPoint));
							copy_v3_v3(points->points[pt].co, co);
							zero_v3(points->points[pt].offset);
							pt++;
						}
					}
				}
			}
		}
	}

	points->totpoints = pt;
}

//XXX TODO maybe remove, hardly used
static void points_from_greasepencil(Object **ob, int totobj, FracPointCloud *points, float mat[4][4], float thresh)
{
	bGPDlayer *gpl;
	bGPDframe *gpf;
	bGPDstroke *gps;
	int pt = points->totpoints, p, o;

	for (o = 0; o < totobj; o++) {
		if ((ob[o]->gpd) && (ob[o]->gpd->layers.first)) {
			float imat[4][4];
			invert_m4_m4(imat, mat);
			for (gpl = ob[o]->gpd->layers.first; gpl; gpl = gpl->next) {
				for (gpf = gpl->frames.first; gpf; gpf = gpf->next) {
					for (gps = gpf->strokes.first; gps; gps = gps->next) {
						for (p = 0; p < gps->totpoints; p++) {
							if (BLI_thread_frand(0) < thresh) {
								float point[3] = {0, 0, 0};
								points->points = MEM_reallocN(points->points, (pt + 1) * sizeof(FracPoint));

								point[0] = gps->points[p].x;
								point[1] = gps->points[p].y;
								point[2] = gps->points[p].z;

								mul_m4_v3(imat, point);

								copy_v3_v3(points->points[pt].co, point);
								zero_v3(points->points[pt].offset);
								pt++;
							}
						}
					}
				}
			}
		}
	}

	points->totpoints = pt;
}

FracPointCloud BKE_fracture_points_get(Depsgraph *depsgraph, FractureModifierData *emd, Object *ob, MeshIsland* mi)
{
	FracPointCloud points;
	Scene* scene = DEG_get_evaluated_scene(depsgraph); //do we need the eval one ?

	/* global settings, for first fracture only, or global secondary and so on fracture, apply to entire fracmesh */
	int totgroup = 0;
	Object **go = MEM_mallocN(sizeof(Object *), "groupobjects");
	float thresh = (float)emd->percentage / 100.0f;
	float min[3], max[3];
	int i;

	points.points = MEM_mallocN(sizeof(FracPoint), "points");
	points.totpoints = 0;

	if (emd->point_source & (MOD_FRACTURE_EXTRA_PARTICLES | MOD_FRACTURE_EXTRA_VERTS)) {
		if (((emd->point_source & MOD_FRACTURE_OWN_PARTICLES) && (emd->point_source & MOD_FRACTURE_EXTRA_PARTICLES)) ||
			((emd->point_source & MOD_FRACTURE_OWN_VERTS) && (emd->point_source & MOD_FRACTURE_EXTRA_VERTS)) ||
			((emd->point_source & MOD_FRACTURE_GREASEPENCIL) && (emd->point_source & MOD_FRACTURE_EXTRA_PARTICLES)) ||
			((emd->point_source & MOD_FRACTURE_GREASEPENCIL) && (emd->point_source & MOD_FRACTURE_EXTRA_VERTS)))
		{
			go = MEM_reallocN(go, sizeof(Object *) * (totgroup + 1));
			go[totgroup] = ob;
			totgroup++;
		}

		totgroup = getGroupObjects(emd->extra_group, &go, totgroup);
	}
	else {
		totgroup = 1;
		go[0] = ob;
	}

	if (emd->point_source & (MOD_FRACTURE_OWN_PARTICLES | MOD_FRACTURE_EXTRA_PARTICLES)) {
		points_from_particles(go, totgroup, scene, &points, ob->obmat, thresh, emd, mi);
	}

	if (emd->point_source & (MOD_FRACTURE_OWN_VERTS | MOD_FRACTURE_EXTRA_VERTS)) {
		points_from_verts(go, totgroup, &points, ob->obmat, thresh, emd, mi, ob);
	}

	if (emd->point_source & MOD_FRACTURE_GREASEPENCIL && !emd->use_greasepencil_edges) {
		points_from_greasepencil(go, totgroup, &points, ob->obmat, thresh);
	}


	/* local settings, apply per shard!!! Or globally too first. */
	if (emd->point_source & MOD_FRACTURE_UNIFORM)
	{
		float cent[3], bmin[3], bmax[3];
		int count = emd->shard_count;

		INIT_MINMAX(min, max);
		//for limit impact we need entire container always, because we need to determine secondary impacts on the shards at their original pos

		BKE_mesh_minmax(mi->mesh, min, max);
		copy_v3_v3(mi->min, min);
		copy_v3_v3(mi->max, max);

		//arrange shards according to their original centroid (parent centroid sum) position in shard-space (else they are centered at 0, 0, 0)
//		arrange_shard(emd, mi->id, false, cent);
		add_v3_v3v3(bmax, max, cent);
		add_v3_v3v3(bmin, min, cent);

#if 0
		//first impact only, so shard has id 0
		if (emd->fracture_mode == MOD_FRACTURE_DYNAMIC) {
			//shrink pointcloud container around impact point, to a size

			copy_v3_v3(max, bmax);
			copy_v3_v3(min, bmin);

			if (mi && mi->impact_size[0] > 0.0f && emd->limit_impact) {
				float size[3], nmin[3], nmax[3], loc[3], tmin[3], tmax[3], rloc[3] = {0,0,0}, quat[4] = {1,0,0,0};
				//MeshIslandSequence *msq = emd->shared->current_mi_entry->prev ? emd->shared->current_mi_entry->prev :
				//                                                                emd->shared->current_mi_entry;
				MeshIsland *mi = NULL;
				RigidBodyOb *rbo = NULL;

				mat4_to_quat(quat, ob->obmat);
				invert_qt(quat);

				if (msq) {
					if (mi->parent) {
						//TODO; this should be the same meshisland later, it should be solved via pointer redirect
						//mi = find_meshisland(&msq->meshIslands, mi->id);
						mi = mi->parent; // ????? TODO FIX
					}

					if (mi) {
						rbo = mi->rigidbody;
						copy_v3_v3(rloc, rbo->pos);
						mul_qt_qtqt(quat, rbo->orn, quat);
					}
				}

				print_v3("Impact Loc\n", mi->impact_loc);
				print_v3("Impact Size\n", mi->impact_size);

				copy_v3_v3(loc, mi->impact_loc);

				sub_v3_v3(loc, rloc);
				mul_qt_v3(quat, loc);
				add_v3_v3(loc, mi->centroid);

				copy_v3_v3(tmax, mi->max);
				copy_v3_v3(tmin, mi->min);

				mul_v3_v3fl(size, mi->impact_size, 0.75f);
				sub_v3_v3v3(nmin, loc, size);
				add_v3_v3v3(nmax, loc, size);

				//clamp
				if (tmin[0] > nmin[0]) {
					nmin[0] = tmin[0];
				}

				if (tmin[1] > nmin[1]) {
					nmin[1] = tmin[1];
				}

				if (tmin[2] > nmin[2]) {
					nmin[2] = tmin[2];
				}

				if (tmax[0] < nmax[0]) {
					nmax[0] = tmax[0];
				}

				if (tmax[1] < nmax[1]) {
					nmax[1] = tmax[1];
				}

				if (tmax[2] < nmax[2]) {
					nmax[2] = tmax[2];
				}

				copy_v3_v3(max, nmax);
				copy_v3_v3(min, nmin);

				/*print_v3("SMIN:", s->min);
				print_v3("SMAX:", s->max);
				print_v3("CENTR", s->centroid);
				print_v3("POS", cent);
				print_v3("MIN:", min);
				print_v3("MAX:", max);*/
			}
		}
#endif

		//omg, vary the seed here
		if (emd->split_islands && emd->fracture_mode == MOD_FRACTURE_DYNAMIC) {
			BLI_thread_srandom(0, mi->id);
		}
		else
		{
			BLI_thread_srandom(0, emd->point_seed);
		}
		for (i = 0; i < count; ++i) {
			if (BLI_thread_frand(0) < thresh) {
				float co[3];
				co[0] = min[0] + (max[0] - min[0]) * BLI_thread_frand(0);
				co[1] = min[1] + (max[1] - min[1]) * BLI_thread_frand(0);
				co[2] = min[2] + (max[2] - min[2]) * BLI_thread_frand(0);

				if (mi->id > 0 && emd->cutter_group == NULL)
				{
					if (in_bbox(co, bmin, bmax))
					{
						points.points = MEM_reallocN(points.points, sizeof(FracPoint) * (points.totpoints + 1));
						copy_v3_v3(points.points[points.totpoints].co, co);
						zero_v3(points.points[points.totpoints].offset);
						points.totpoints++;
					}
				}
				else
				{
					points.points = MEM_reallocN(points.points, sizeof(FracPoint) * (points.totpoints + 1));
					copy_v3_v3(points.points[points.totpoints].co, co);
					zero_v3(points.points[points.totpoints].offset);
					points.totpoints++;
				}
			}
		}
	}

	if (emd->point_source & MOD_FRACTURE_GRID)
	{
		float cent[3], bmin[3], bmax[3];
		int x, y, z, k = 0;

		if (emd->grid_resolution[0] < 1)
		{	//sanity check
			emd->grid_resolution[0] = 1;
		}

		if (emd->grid_resolution[1] < 1)
		{	//sanity check
			emd->grid_resolution[1] = 1;
		}

		if (emd->grid_resolution[2] < 1)
		{	//sanity check
			emd->grid_resolution[2] = 1;
		}

		//just draw over bbox
		INIT_MINMAX(min, max);
		//for limit impact we need entire container always, because we need to determine secondary impacts on the shards at their original pos
		//if (mi->id == 0) //TODO FIX....
		//	return points; //id 0 should be entire mesh

		BKE_mesh_minmax(mi->mesh, min, max);
		copy_v3_v3(mi->min, min);
		copy_v3_v3(mi->max, max);

		//arrange shards according to their original centroid (parent centroid sum) position in shard-space (else they are centered at 0, 0, 0)
//		arrange_shard(emd, mi->id, false, cent);
		add_v3_v3v3(bmax, max, cent);
		add_v3_v3v3(bmin, min, cent);

		//centroid of grid cells is
		for (z = 0; z < emd->grid_resolution[2]; z++) {
			for (y = 0; y < emd->grid_resolution[1]; y++) {
				for (x = 0; x < emd->grid_resolution[0]; x++) {
					float co[3], off[3] =  {0, 0, 0};
					co[0] = min[0] + (max[0] - min[0]) * ((float)x + 0.5f)/(float)emd->grid_resolution[0];
					co[1] = min[1] + (max[1] - min[1]) * ((float)y + 0.5f)/(float)emd->grid_resolution[1];
					co[2] = min[2] + (max[2] - min[2]) * ((float)z + 0.5f)/(float)emd->grid_resolution[2];


					//alternating offset for bricks
					if (((x % 2) == 1) && (emd->grid_offset[0] == 0))
					{
						off[1] = emd->grid_offset[1] * ((max[1] - min[1]) / (float)emd->grid_resolution[1]);
						off[2] = emd->grid_offset[2] * ((max[2] - min[2]) / (float)emd->grid_resolution[2]);
					}

					if (((y % 2) == 1) && (emd->grid_offset[1] == 0))
					{
						off[0] = emd->grid_offset[0] * ((max[0] - min[0]) / (float)emd->grid_resolution[0]);
						off[2] = emd->grid_offset[2] * ((max[2] - min[2]) / (float)emd->grid_resolution[2]);
					}

					if (((z % 2) == 1) && (emd->grid_offset[2] == 0))
					{
						off[0] = emd->grid_offset[0] * ((max[0] - min[0]) / (float)emd->grid_resolution[0]);
						off[1] = emd->grid_offset[1] * ((max[1] - min[1]) / (float)emd->grid_resolution[1]);
					}

					//print_v3("offset", off);

					if (mi->id > 0 && emd->cutter_group == NULL)
					{
						if (in_bbox(co, bmin, bmax))
						{
							points.points = MEM_reallocN(points.points, sizeof(FracPoint) * (points.totpoints + 1));
							copy_v3_v3(points.points[points.totpoints].co, co);
							copy_v3_v3(points.points[points.totpoints].offset, off);
							points.totpoints++;
						}
					}
					else
					{
						points.points = MEM_reallocN(points.points, sizeof(FracPoint) * (points.totpoints + 1));
						copy_v3_v3(points.points[points.totpoints].co, co);
						copy_v3_v3(points.points[points.totpoints].offset, off);
						points.totpoints++;
					}

					k++;
				}
			}
		}
	}

	MEM_freeN(go);
	return points;
}

static Material* find_material(Main* bmain, const char* name)
{
	ID* mat;

	for (mat = bmain->mat.first; mat; mat = mat->next)
	{
		char *cmp = BLI_strdupcat("MA", name);
		if (strcmp(cmp, mat->name) == 0)
		{
			MEM_freeN(cmp);
			cmp = NULL;
			return (Material*)mat;
		}
		else
		{
			MEM_freeN(cmp);
			cmp = NULL;
		}
	}

	return BKE_material_add(bmain, name);
}

//splinter handling is a case for BKE too
static MeshIsland* do_splinters(FractureModifierData *fmd, FracPointCloud points, float(*mat)[4][4], MeshIsland *mi)
{
	float imat[4][4];

	/*need to add island / shard centroid...*/
	unit_m4(*mat);

	/*splinters... just global axises and a length, for rotation rotate the object */
	if (fmd->splinter_axis & MOD_FRACTURE_SPLINTER_X)
	{
		(*mat)[0][0] *= fmd->splinter_length;
	}
	if (fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Y)
	{
		(*mat)[1][1] *= fmd->splinter_length;
	}
	if (fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Z)
	{
		(*mat)[2][2] *= fmd->splinter_length;
	}

	if ((fmd->splinter_axis & MOD_FRACTURE_SPLINTER_X) ||
		(fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Y) ||
		(fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Z))
	{
		int i = 0, num_verts = 0;
		MVert* mvert = NULL, *mv;
		mvert = mi->mesh->mvert;
		num_verts =  mi->mesh->totvert;
		invert_m4_m4(imat, *mat);

		for (i = 0; i < points.totpoints; i++)
		{
			mul_m4_v3(imat, points.points[i].co);
		}

		for (i = 0, mv = mvert; i < num_verts; i++, mv++)
		{
			mul_m4_v3(imat, mv->co);
		}
	}

	return mi;
}

//so is material handling too, XXX TODO move to BKE
static short do_materials(Main* bmain, FractureModifierData *fmd, Object* obj)
{
	short mat_index = 0;
	//Object *obj = DEG_get_original_object(ob);

	if (fmd->inner_material) {
		/* assign inner material as secondary mat to ob if not there already */
		mat_index = BKE_object_material_slot_find_index(obj, fmd->inner_material);
		if (mat_index == 0) {
			BKE_object_material_slot_add(bmain, obj);
			assign_material(bmain, obj, fmd->inner_material, obj->totcol, BKE_MAT_ASSIGN_OBDATA);
		}

		/* get index again */
		mat_index = BKE_object_material_slot_find_index(obj, fmd->inner_material);
	}
	else
	{
		/* autogenerate materials */
		char name[MAX_ID_NAME];

		short* totmat = give_totcolp(obj);

		BLI_strncpy(name, obj->id.name + 2, strlen(obj->id.name));
		if (*totmat == 0)
		{
			/*create both materials if no material is present*/
			Material* mat_inner;
			char *matname = BLI_strdupcat(name, "_Outer");
			Material* mat_outer = find_material(bmain, matname);
			BKE_object_material_slot_add(bmain, obj);
			assign_material(bmain, obj, mat_outer, obj->totcol, BKE_MAT_ASSIGN_OBDATA);

			MEM_freeN(matname);
			matname = NULL;
			matname = BLI_strdupcat(name, "_Inner");
			mat_inner = find_material(bmain, matname);
			BKE_object_material_slot_add(bmain, obj);
			assign_material(bmain, obj, mat_inner, obj->totcol, BKE_MAT_ASSIGN_OBDATA);

			MEM_freeN(matname);
			matname = NULL;

			fmd->inner_material = mat_inner;
			mat_index = 1;
		}
		else if (*totmat > 0)
		{
			/* append inner material to the stack if materials are present */
			char* matname = BLI_strdupcat(name, "_Inner");
			Material* mat_inner = find_material(bmain, matname);
			BKE_object_material_slot_add(bmain, obj);
			assign_material(bmain, obj, mat_inner, obj->totcol, BKE_MAT_ASSIGN_OBDATA);
			MEM_freeN(matname);
			matname = NULL;

			fmd->inner_material = mat_inner;
			mat_index = *totmat;
		}
	}

	return mat_index;
}

static void cleanup_splinters(FractureModifierData *fmd, float mat[4][4], MeshIsland *mi)
{
	if ((fmd->splinter_axis & MOD_FRACTURE_SPLINTER_X) ||
		(fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Y) ||
		(fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Z))
	{
		int i = 0, num_verts = 0;
		MVert* mvert = NULL, *mv;

		mvert = mi->mesh->mvert;
		num_verts = mi->mesh->totvert;

		for (i = 0, mv = mvert; i < num_verts; i++, mv++)
		{
			mul_m4_v3(mat, mv->co);
		}
	}
}

#if 0
static void arrange_shard(FractureModifierData *fmd, MeshIsland* mi, bool do_verts, float cent[3])
{
	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		zero_v3(cent);
		bool found = false;
		MeshIslandSequence *seq = fmd->shared->current_mi_entry;
		while (seq->prev && mi)
		{
			MeshIsland *par = mi->parent;
			if (par)
			{
				add_v3_v3(cent, par->centroid);
				found = true;
			}

			seq = seq->prev;
		}

		if (found && do_verts)
		{
			add_v3_v3(mi->centroid, cent);
			add_v3_v3(mi->min, cent);
			add_v3_v3(mi->max, cent);

			int i = 0;
			for (i = 0; i < mi->mesh->totvert; i++)
			{
				add_v3_v3(mi->mesh->mvert[i].co, cent);
			}
		}
	}
}
#endif

#if 0
static void cleanup_arrange_shard(FractureModifierData *fmd, MeshIsland *sh, float cent[3])
{
	int i = 0;
	sub_v3_v3(sh->centroid, cent);

	for (i = 0; i < sh->mesh->totvert; i++)
	{
		sub_v3_v3(sh->mesh->mvert[i].co, cent);
	}
}
#endif


void BKE_fracture_points(FractureModifierData *fmd, Object* obj, MeshIsland *mi, Depsgraph* depsgraph, Main *bmain,
                         Scene *scene)
{
	/* dummy point cloud, random */
	FracPointCloud points;

	points = BKE_fracture_points_get(depsgraph, fmd, obj, mi);

	if (points.totpoints > 0 || fmd->use_greasepencil_edges) {
		short mat_index = 0;
		float mat[4][4];
		float cent[3];

//		arrange_shard(fmd, mi, true, cent);

		/*splinters... just global axises and a length, for rotation rotate the object */
		mi = do_splinters(fmd, points, &fmd->shared->splinter_matrix, mi);


		mat_index = 1; //do_materials(bmain, fmd, obj);
		mat_index = mat_index > 0 ? mat_index - 1 : mat_index;

		if (points.totpoints > 0) {
			BKE_fracture_shard_by_points(fmd, &points, obj, mi, scene, bmain);
		}

		/*TODO, limit this to settings shards !*/
		if (fmd->point_source & MOD_FRACTURE_GREASEPENCIL && fmd->use_greasepencil_edges) {
			BKE_fracture_shard_by_greasepencil(fmd, obj, mat_index, mat);
		}

		if (mi->mesh)
		{
//			cleanup_arrange_shard(fmd, mi, cent);

			/* here we REALLY need to fracture so deactivate the shards to islands flag and activate afterwards */
			//fmd->shards_to_islands = false;
			//BKE_fracture_assemble_mesh_from_islands(fmd, &fmd->shared->mesh_islands, obj);
			//fmd->shards_to_islands = temp;

			cleanup_splinters(fmd, fmd->shared->splinter_matrix, mi);
		}

#if 0
		/* The original mesh island must not be in the collection any more, so unlink and free */
		BLI_remlink(&fmd->shared->mesh_islands, mi);
		BKE_fracture_mesh_island_free(mi, scene);
		mi = NULL;
#endif

		if (!fmd->auto_execute && fmd->execute_threaded) {
			fmd->shared->reset_shards = false;
		}
	}
	MEM_freeN(points.points);
}

//this is the main fracture function, outsource to BKE, so op or rb system can call it
void BKE_fracture_do(FractureModifierData *fmd, MeshIsland *mi, Object *obj, Depsgraph *depsgraph, Main* bmain, Scene *scene)
{
	MeshIsland *mii = NULL;

	/* no pointsource means re-use existing mesh islands*/
	if (fmd->point_source == 0) {
		int count = 1, i, j = 0;
		int frame = (int)(BKE_scene_frame_get(scene)); //TODO ensure original scene !!!
		Mesh** temp_meshs = MEM_callocN(sizeof(Mesh*) * count, "temp_islands no pointsource");
		BKE_fracture_split_islands(fmd, obj, mi->mesh, &temp_meshs, &count);

		/* The original mesh island must not be in the collection any more, so unlink and free */
		BLI_remlink(&fmd->shared->mesh_islands, mi);
		BKE_fracture_mesh_island_free(mi, scene);
		mi = NULL;

		for (i = 0; i < count; i++)
		{
			if (temp_meshs[i])
			{
				if (temp_meshs[i]->totvert > 0)
				{	/* skip invalid cells, e.g. those which are eliminated by bisect */
					MeshIsland *result = BKE_fracture_mesh_island_create(temp_meshs[i], bmain, scene, obj, frame);
					fracture_meshisland_add(fmd, result);
					result->id = j;
					j++;
				}
				else {
					BKE_fracture_mesh_free(temp_meshs[i]);
				}
			}
		}

		MEM_freeN(temp_meshs);

		return;
	}

	if (fmd->cutter_group != NULL) {
		//attempt to combine fracture by cutter group with regular fracture
		float mat[4][4];
		bool reset = fmd->shared->reset_shards;

		unit_m4(mat);
		//mat_index = do_materials(fmd, obj);
		//mat_index = mat_index > 0 ? mat_index - 1 : mat_index;

//		BKE_fracture_shard_by_planes(fmd, obj, mat_index, mat);

		fmd->shared->reset_shards = false;

		for (mii = fmd->shared->mesh_islands.first; mii; mii = mii->next)
		{
			printf("Fracturing MeshIsland %d\n", mii->id);
			BKE_fracture_points(fmd, obj, mii, depsgraph, bmain, scene);
		}

		fmd->shared->reset_shards = reset;
	}
	else {
		BKE_fracture_points(fmd, obj, mi, depsgraph, bmain, scene);
	}
}

//XXXX TODO, is BB really useds still ? aint there exact volume calc now ?
/* mi->bb, its for volume fraction calculation.... */
static float bbox_vol(BoundBox *bb)
{
	float x[3], y[3], z[3];

	sub_v3_v3v3(x, bb->vec[4], bb->vec[0]);
	sub_v3_v3v3(y, bb->vec[3], bb->vec[0]);
	sub_v3_v3v3(z, bb->vec[1], bb->vec[0]);

	return len_v3(x) * len_v3(y) * len_v3(z);
}

static void bbox_dim(BoundBox *bb, float dim[3])
{
	float x[3], y[3], z[3];

	sub_v3_v3v3(x, bb->vec[4], bb->vec[0]);
	sub_v3_v3v3(y, bb->vec[3], bb->vec[0]);
	sub_v3_v3v3(z, bb->vec[1], bb->vec[0]);

	dim[0] = len_v3(x);
	dim[1] = len_v3(y);
	dim[2] = len_v3(z);
}

//This still necessary ? if yes move to fracture.c for now
static int BM_calc_center_centroid(BMesh *bm, float cent[3], int tagged)
{
	BMFace *f;
	BMIter iter;
	float face_area;
	float total_area = 0.0f;
	float face_cent[3];

	zero_v3(cent);

	/* calculate a weighted average of face centroids */
	BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
		if (BM_elem_flag_test(f, BM_ELEM_TAG) || !tagged) {
			BM_face_calc_center_mean(f, face_cent);
			face_area = BM_face_calc_area(f);

			madd_v3_v3fl(cent, face_cent, face_area);
			total_area += face_area;
		}
	}
	/* otherwise we get NAN for 0 polys */
	if (bm->totface) {
		mul_v3_fl(cent, 1.0f / total_area);
	}
	else if (bm->totvert == 1) {
		copy_v3_v3(cent, BM_vert_at_index_find(bm, 0)->co);
	}

	return (bm->totface != 0);
}

#if 0 //TODO FIX...
static short do_vert_index_map(FractureModifierData *fmd, MeshIsland *mi, MeshIsland *par)
{
	short rb_type = mi->ground_weight > 0.01f ?  RBO_TYPE_PASSIVE : (par && par->rigidbody ? par->rigidbody->type : RBO_TYPE_ACTIVE);

	if (fmd->shared->vert_index_map && fmd->dm_group && fmd->cluster_count == 0 && mi->vertex_indices)
	{
		CollectionObject* go = NULL;
		/* autocreate clusters out of former objects, if we dont override */
		mi->particle_index = GET_INT_FROM_POINTER(BLI_ghash_lookup(fmd->shared->vert_index_map,
		                                          SET_INT_IN_POINTER(mi->vertex_indices[0])));

		/*look up whether original object is active or passive */
		go = BLI_findlink(&fmd->dm_group->gobject, mi->particle_index);
		if (go && go->ob && go->ob->rigidbody_object) {
			rb_type = go->ob->rigidbody_object->type;
		}
	}

	return rb_type;
}
#endif

#if 0 //TODO FIX
static void do_fix_normals(FractureModifierData *fmd, MeshIsland *mi)
{
	/* copy fixed normals to physicsmesh too, for convert to objects */
	if (fmd->fix_normals) {
		MVert *verts, *mv;
		int j = 0, totvert = 0;
		totvert = mi->mesh->totvert;
		verts = mi->mesh->mvert;
		for (mv = verts, j = 0; j < totvert; mv++, j++) {
			short no[3];
			no[0] = mi->vertno[j * 3];
			no[1] = mi->vertno[j * 3 + 1];
			no[2] = mi->vertno[j * 3 + 2];

			copy_v3_v3_short(mv->no, no);
		}
	}
}
#endif

#if 0
static float do_setup_meshisland(FractureModifierData *fmd, Object *ob, int totvert, float centroid[3],
								 BMVert **verts, float *vertco, short *vertno, BMesh **bm_new, Mesh *orig_dm,
								 int id, Scene *scene, Main *bmain)
{
	MeshIsland *mi;
	Mesh *dm;
	float min[3], max[3], vol = 0;
	int i = 0;
	short rb_type = RBO_TYPE_ACTIVE;
	//struct BMeshToMeshParams bmt = {.calc_object_remap = 0};

	mi = MEM_callocN(sizeof(MeshIsland), "meshIsland");
	unit_qt(mi->rot);

#if 0
	if (fmd->fracture_mode != MOD_FRACTURE_DYNAMIC)
	{
		mi->locs = NULL;
		mi->rots = NULL;
		mi->frame_count = 0;

		if (scene->rigidbody_world)
		{
			mi->start_frame = scene->rigidbody_world->shared->pointcache->startframe;
		}
		else
		{
			mi->start_frame = 1;
		}
	}
	else
#endif
	{
		/* in dynamic case preallocate cache here */
		int start = scene->rigidbody_world->shared->pointcache->startframe;
		int end = 10; //fmd->modifier.scene->rigidbody_world->pointcache->endframe;

		if (fmd->shared->current_mi_entry) {
			MeshIslandSequence *prev = fmd->shared->current_mi_entry->prev;
			if (prev)
			{
				start = prev->frame;
			}
		}

		end = start + 10;

		mi->frame_count = end - start + 1;
		mi->start_frame = start;
		mi->locs = MEM_mallocN(sizeof(float)*3* mi->frame_count, "mi->locs");
		mi->rots = MEM_mallocN(sizeof(float)*4* mi->frame_count, "mi->rots");
	}

	mi->particle_index = -1;
	mi->thresh_weight = 0.0f;
	mi->ground_weight = 0.0f;
	mi->id = id;
	BLI_snprintf(mi->name, 64, "%d", mi->id);

	BM_mesh_normals_update(*bm_new);
	dm = BKE_fracture_bmesh_to_mesh(*bm_new);
	BM_mesh_free(*bm_new);
	*bm_new = NULL;

	mi->mesh = dm;

	//do_fix_normals(fmd, mi); TODO FIX...

	copy_v3_v3(mi->centroid, centroid);
	//mi->bb = BKE_boundbox_alloc_unit();
	//BKE_boundbox_init_from_minmax(mi->bb, min, max);
	mi->participating_constraints = NULL;
	mi->participating_constraint_count = 0;

/*	vol = bbox_vol(mi->bb);
	if (vol > fmd->max_vol) {
		fmd->max_vol = vol;
	}*/

//	mi->vertices_cached = NULL;

#if 0
	TODO FIX
	rb_type = do_vert_index_map(fmd, mi, NULL);
	i = BLI_listbase_count(&fmd->shared->meshIslands);
	do_rigidbody(bmain, scene, mi, DEG_get_original_object(ob), orig_dm, rb_type, i);
#endif

	mi->start_frame = scene->rigidbody_world->shared->pointcache->startframe;
	BLI_addtail(&fmd->shared->mesh_islands, mi);

	return vol;
}

#endif

static void mesh_separate_tagged(FractureModifierData *fmd, Object *ob, BMVert** v_tag, int v_count,
								  BMesh *bm_work, Mesh*** temp_islands, int* count)
{
	BMesh *bm_new;
	BMesh *bm_old = bm_work;
	Mesh *me = NULL;
	struct BMeshCreateParams bmc = {};
	bmc.use_toolflags = true;

	bm_new = BM_mesh_create(&bm_mesh_allocsize_default, &bmc);
	BM_mesh_elem_toolflags_ensure(bm_new);  /* needed for 'duplicate' bmo */

	CustomData_copy(&bm_old->vdata, &bm_new->vdata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_old->edata, &bm_new->edata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_old->ldata, &bm_new->ldata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_old->pdata, &bm_new->pdata, CD_MASK_BMESH, CD_CALLOC, 0);

	CustomData_bmesh_init_pool(&bm_new->vdata, bm_mesh_allocsize_default.totvert, BM_VERT);
	CustomData_bmesh_init_pool(&bm_new->edata, bm_mesh_allocsize_default.totedge, BM_EDGE);
	CustomData_bmesh_init_pool(&bm_new->ldata, bm_mesh_allocsize_default.totloop, BM_LOOP);
	CustomData_bmesh_init_pool(&bm_new->pdata, bm_mesh_allocsize_default.totface, BM_FACE);

	BMO_op_callf(bm_old, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
				 "duplicate geom=%hvef dest=%p", BM_ELEM_TAG, bm_new);

	//BM_calc_center_centroid(bm_new, centroid, false);
	BM_mesh_elem_index_ensure(bm_new, BM_VERT | BM_EDGE | BM_FACE);

#if 0
	BM_ITER_MESH (v, &iter, bm_new, BM_VERTS_OF_MESH) {
		/* eliminate centroid in vertex coords */
		sub_v3_v3(v->co, centroid);
	}
#endif

	me = BKE_fracture_bmesh_to_mesh(bm_new);
	*temp_islands = MEM_reallocN(*temp_islands, sizeof(Mesh*) * ((*count)+1));
	(*temp_islands)[(*count)] = me;
	(*count)++;

	BM_mesh_free(bm_new);

	/* deselect loose data - this used to get deleted,
	 * we could de-select edges and verts only, but this turns out to be less complicated
	 * since de-selecting all skips selection flushing logic */
	BM_mesh_elem_hflag_disable_all(bm_old, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_TAG, false);
}

static void handle_vert(BMVert* vert, BMVert** orig_work,
						BMVert*** v_tag, int *tot, int *tag_counter)
{
	/* treat the specified vert and put it into the tagged array, also store its coordinates and normals
	 * for usage in meshislands later on */

	//short no[3];
	//short vno[3];

	if (*v_tag == NULL)
		*v_tag = MEM_callocN(sizeof(BMVert *), "v_tag");

	//if (*startco == NULL)
	//	*startco = MEM_callocN(sizeof(float), "mesh_separate_loose->startco");

	//if (*startno == NULL)
	//	*startno = MEM_callocN(sizeof(short), "mesh_separate_loose->startno");

	BM_elem_flag_enable(vert, BM_ELEM_TAG);
	BM_elem_flag_enable(vert, BM_ELEM_INTERNAL_TAG);
	*v_tag = MEM_reallocN(*v_tag, sizeof(BMVert *) * ((*tag_counter) + 1));
	(*v_tag)[(*tag_counter)] = orig_work[vert->head.index];

	/**startco = MEM_reallocN(*startco, ((*tag_counter) + 1) * 3 * sizeof(float));
	(*startco)[3 * (*tag_counter)] = vert->co[0];
	(*startco)[3 * (*tag_counter) + 1] = vert->co[1];
	(*startco)[3 * (*tag_counter) + 2] = vert->co[2];

	*startno = MEM_reallocN(*startno, ((*tag_counter) + 1) * 3 * sizeof(short));*/

	/*normal_float_to_short_v3(vno, vert->no);
	normal_float_to_short_v3(no, vert->no);
	if (fmd->fix_normals)
		BKE_fracture_normal_find(dm, fmd->shared->nor_tree, vert->co, vno, no, fmd->nor_range);
	(*startno)[3 * (*tag_counter)] = no[0];
	(*startno)[3 * (*tag_counter) + 1] = no[1];
	(*startno)[3 * (*tag_counter) + 2] = no[2];*/

	(*tot)++;
	(*tag_counter)++;
}

static void mesh_separate_loose_partition(FractureModifierData *fmd, Object *ob, BMesh *bm_work, BMVert **orig_work,
										  Mesh *** temp_islands, int* count)
{
	int i, tag_counter = 0;
	BMEdge *e;
	BMVert *v_seed = NULL, **v_tag = NULL;
	BMWalker walker;
	int tot = 0;
	BMesh *bm_old = bm_work;
	int max_iter = bm_old->totvert;
	BMIter iter;


	/* Clear all selected vertices */
	BM_mesh_elem_hflag_disable_all(bm_old, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_INTERNAL_TAG | BM_ELEM_TAG, false);


	/* A "while (true)" loop should work here as each iteration should
	 * select and remove at least one vertex and when all vertices
	 * are selected the loop will break out. But guard against bad
	 * behavior by limiting iterations to the number of vertices in the
	 * original mesh.*/
	for (i = 0; i < max_iter; i++) {
		tag_counter = 0;

		BM_ITER_MESH (v_seed, &iter, bm_old, BM_VERTS_OF_MESH) {
			/* Hrm need to look at earlier verts to for unused ones.*/
			if (!BM_elem_flag_test(v_seed, BM_ELEM_TAG) && !BM_elem_flag_test(v_seed, BM_ELEM_INTERNAL_TAG)) {
				break;
			}
		}

		/* No vertices available, can't do anything */
		if (v_seed == NULL) {
			break;
		}
		/* Select the seed explicitly, in case it has no edges */
		if (!BM_elem_flag_test(v_seed, BM_ELEM_TAG) && !BM_elem_flag_test(v_seed, BM_ELEM_INTERNAL_TAG)) {
			handle_vert(v_seed, orig_work, &v_tag, &tot, &tag_counter);
		}

		/* Walk from the single vertex, selecting everything connected
		 * to it */
		BMW_init(&walker, bm_old, BMW_VERT_SHELL,
				 BMW_MASK_NOP, BMW_MASK_NOP, BMW_MASK_NOP,
				 BMW_FLAG_NOP,
				 BMW_NIL_LAY);

		e = BMW_begin(&walker, v_seed);
		for (; e; e = BMW_step(&walker)) {
			if (!BM_elem_flag_test(e->v1, BM_ELEM_TAG) && !BM_elem_flag_test(e->v1, BM_ELEM_INTERNAL_TAG)) {
				handle_vert(e->v1, orig_work, &v_tag, &tot, &tag_counter);
			}
			if (!BM_elem_flag_test(e->v2, BM_ELEM_TAG) && !BM_elem_flag_test(e->v2, BM_ELEM_INTERNAL_TAG)) {
				handle_vert(e->v2, orig_work, &v_tag, &tot, &tag_counter);
			}
		}
		BMW_end(&walker);

		/* Flush the selection to get edge/face selections matching
		 * the vertex selection */
		BKE_bm_mesh_hflag_flush_vert(bm_old, BM_ELEM_TAG);

		/* Move selection into a separate object */
		mesh_separate_tagged(fmd, ob, v_tag, tag_counter, bm_old, temp_islands, count);

		MEM_freeN(v_tag);
		v_tag = NULL;

		if (tot >= bm_old->totvert) {
			break;
		}
	}
}

/* inlined select_linked functionality here, because not easy to reach without modifications */
static void select_linked(BMesh **bm_in)
{
	BMIter iter;
	BMVert *v;
	BMEdge *e;
	BMWalker walker;
	BMesh *bm_work = *bm_in;


	BM_ITER_MESH (v, &iter, bm_work, BM_VERTS_OF_MESH) {
		if (BM_elem_flag_test(v, BM_ELEM_SELECT)) {
			BM_elem_flag_enable(v, BM_ELEM_TAG);
		}
		else {
			BM_elem_flag_disable(v, BM_ELEM_TAG);
		}
	}

	BMW_init(&walker, bm_work, BMW_VERT_SHELL,
			 BMW_MASK_NOP, BMW_MASK_NOP, BMW_MASK_NOP,
			 BMW_FLAG_TEST_HIDDEN,
			 BMW_NIL_LAY);

	BM_ITER_MESH (v, &iter, bm_work, BM_VERTS_OF_MESH) {
		if (BM_elem_flag_test(v, BM_ELEM_TAG)) {
			for (e = BMW_begin(&walker, v); e; e = BMW_step(&walker)) {
				BM_edge_select_set(bm_work, e, true);
			}
		}
	}
	BMW_end(&walker);

	BM_mesh_select_flush(bm_work);
}

static void mesh_separate_selected(BMesh **bm_work, BMesh **bm_out, BMVert **orig_work, BMVert ***orig_out1, BMVert ***orig_out2)
{
	BMesh *bm_old = *bm_work;
	BMesh *bm_new = *bm_out;
	BMVert *v, **orig_new = *orig_out1, **orig_mod = *orig_out2;
	BMIter iter;
	int new_index = 0, mod_index = 0;

	BM_mesh_elem_hflag_disable_all(bm_old, BM_FACE | BM_EDGE | BM_VERT, BM_ELEM_TAG, false);
	/* sel -> tag */
	BM_mesh_elem_hflag_enable_test(bm_old, BM_FACE | BM_EDGE | BM_VERT, BM_ELEM_TAG, true, false, BM_ELEM_SELECT);

	BM_mesh_elem_toolflags_ensure(bm_new);  /* needed for 'duplicate' bmo */

	CustomData_copy(&bm_old->vdata, &bm_new->vdata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_old->edata, &bm_new->edata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_old->ldata, &bm_new->ldata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_old->pdata, &bm_new->pdata, CD_MASK_BMESH, CD_CALLOC, 0);

	CustomData_bmesh_init_pool(&bm_new->vdata, bm_mesh_allocsize_default.totvert, BM_VERT);
	CustomData_bmesh_init_pool(&bm_new->edata, bm_mesh_allocsize_default.totedge, BM_EDGE);
	CustomData_bmesh_init_pool(&bm_new->ldata, bm_mesh_allocsize_default.totloop, BM_LOOP);
	CustomData_bmesh_init_pool(&bm_new->pdata, bm_mesh_allocsize_default.totface, BM_FACE);

	BMO_op_callf(bm_old, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
				 "duplicate geom=%hvef dest=%p", BM_ELEM_TAG, bm_new);

	/* lets hope the order of elements in new mesh is the same as it was in old mesh */
	BM_ITER_MESH (v, &iter, bm_old, BM_VERTS_OF_MESH) {
		if (BM_elem_flag_test(v, BM_ELEM_TAG)) {
			orig_new[new_index] = orig_work[v->head.index];
			new_index++;
		}
		else {
			orig_mod[mod_index] = orig_work[v->head.index];
			mod_index++;
		}
	}

	new_index = 0;
	BM_ITER_MESH (v, &iter, bm_new, BM_VERTS_OF_MESH) {
		v->head.index = new_index;
		new_index++;
	}

	BMO_op_callf(bm_old, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
				 "delete geom=%hvef context=%i", BM_ELEM_TAG, DEL_FACES);

	/* deselect loose data - this used to get deleted,
	 * we could de-select edges and verts only, but this turns out to be less complicated
	 * since de-selecting all skips selection flushing logic */
	BM_mesh_elem_hflag_disable_all(bm_old, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_SELECT | BM_ELEM_TAG, false);

	BM_mesh_normals_update(bm_new);
}

static void halve(FractureModifierData *rmd, Object *ob, int minsize, BMesh **bm_work, BMVert ***orig_work, bool separated,
				  Mesh*** temp_islands, int *count)
{

	int half;
	int i = 0, new_count = 0;
	BMIter iter;
	BMVert **orig_old = *orig_work, **orig_new, **orig_mod;
	BMVert *v;
	BMesh *bm_old = *bm_work;
	BMesh *bm_new = NULL;
	separated = false;

	bm_new = BM_mesh_create(&bm_mesh_allocsize_default, &((struct BMeshCreateParams){.use_toolflags = true,}));

	BM_mesh_elem_hflag_disable_all(bm_old, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_SELECT | BM_ELEM_TAG, false);

	half = bm_old->totvert / 2;
	BM_ITER_MESH (v, &iter, bm_old, BM_VERTS_OF_MESH) {
		if (i >= half) {
			break;
		}
		BM_elem_select_set(bm_old, (BMElem *)v, true);
		i++;
	}

	BKE_bm_mesh_hflag_flush_vert(bm_old, BM_ELEM_SELECT);
	select_linked(&bm_old);

	new_count = bm_old->totvertsel;
	printf("Halving...%d => %d %d\n", bm_old->totvert, new_count, bm_old->totvert - new_count);

	orig_new = MEM_callocN(sizeof(BMVert *) * new_count, "orig_new");
	orig_mod = MEM_callocN(sizeof(BMVert *) * bm_old->totvert - new_count, "orig_mod");
	mesh_separate_selected(&bm_old, &bm_new, orig_old, &orig_new, &orig_mod);

	//printf("Old New: %d %d\n", bm_old->totvert, bm_new->totvert);
	if ((bm_old->totvert <= minsize && bm_old->totvert > 0) || (bm_new->totvert == 0)) {
		mesh_separate_loose_partition(rmd, ob, bm_old, orig_mod, temp_islands, count);
		separated = true;
	}

	if ((bm_new->totvert <= minsize && bm_new->totvert > 0) || (bm_old->totvert == 0)) {
		mesh_separate_loose_partition(rmd, ob, bm_new, orig_new, temp_islands, count);
		separated = true;
	}

	if ((bm_old->totvert > minsize && bm_new->totvert > 0) || (bm_new->totvert == 0 && !separated)) {
		halve(rmd, ob, minsize, &bm_old, &orig_mod, separated, temp_islands, count);
	}

	if ((bm_new->totvert > minsize && bm_old->totvert > 0) || (bm_old->totvert == 0 && !separated)) {
		halve(rmd, ob, minsize, &bm_new, &orig_new, separated, temp_islands, count);
	}


	MEM_freeN(orig_mod);
	MEM_freeN(orig_new);
	BM_mesh_free(bm_new);
	bm_new = NULL;
}

static void mesh_separate_loose(BMesh* bm_work, FractureModifierData *rmd, Object *ob, Mesh*** temp_islands, int* count)
{
	int minsize = 500;
	BMVert *vert, **orig_start;
	BMIter iter;

	//TODO FIX work on copy or regenerate bmesh
	BM_mesh_elem_hflag_disable_all(bm_work, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_SELECT | BM_ELEM_TAG, false);

	orig_start = MEM_callocN(sizeof(BMVert *) * bm_work->totvert, "orig_start");
	/* associate new verts with old verts, here indexes should match still */
	BM_ITER_MESH (vert, &iter, bm_work, BM_VERTS_OF_MESH)
	{
		orig_start[vert->head.index] = vert;
	}

	BM_mesh_elem_index_ensure(bm_work, BM_VERT);
	BM_mesh_elem_table_ensure(bm_work, BM_VERT);

	halve(rmd, ob, minsize, &bm_work, &orig_start, false, temp_islands, count);

	MEM_freeN(orig_start);
	orig_start = NULL;
	BM_mesh_free(bm_work);
	bm_work = NULL;
}

#if 0
//TODO FIX
void BKE_fracture_fill_vgroup(FractureModifierData *rmd, Mesh *dm, MDeformVert *dvert, Object *ob, Mesh *old_cached)
{
	/* use fallback over inner material (no more, now directly via tagged verts) */
	if (rmd->inner_defgrp_name[0]) {
		int ind = 0, mat_index = BKE_object_material_slot_find_index(ob, rmd->inner_material);
		bool fallback = false, dynamic = false;
		MPoly *mp = dm->mpoly, *p;
		MLoop *ml = dm->mloop;
		MVert *mv = dm->mvert;
		int totpoly = dm->totpoly;
		int totvert = dm->totvert;
		const int inner_defgrp_index = defgroup_name_index(ob, rmd->inner_defgrp_name);
		MDeformVert *old_dvert = NULL;
		int old_totvert = 0;

		dynamic = rmd->fracture_mode == MOD_FRACTURE_DYNAMIC;
		fallback = rmd->frac_algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL;
		dvert = dm->dvert;

		if (dvert == NULL)
		{
			dvert = CustomData_add_layer(&dm->vdata, CD_MDEFORMVERT, CD_CALLOC,
									 NULL, totvert);
		}

		if (old_cached) {
			old_dvert = old_cached->dvert;
			old_totvert = old_cached->totvert;
		}

		for (ind = 0, p = mp; ind < totpoly; ind++, p++) {
			int j;
			for (j = 0; j < p->totloop; j++) {
				MLoop *l;
				MVert *v;
				MDeformVert *dv;
				int l_index = p->loopstart + j;
				l = ml + l_index;
				v = mv + l->v;
				dv = dvert + l->v;

				if (dv->dw == NULL)
				{
					if ((v->flag & ME_VERT_TMP_TAG) && !fallback) {
						defvert_add_index_notest(dv, inner_defgrp_index, 1.0f);
					}
					else if ((p->mat_nr == mat_index-1) && fallback) {
						defvert_add_index_notest(dv, inner_defgrp_index, 1.0f);
					}
					else {
						defvert_add_index_notest(dv, inner_defgrp_index, 0.0f);
					}
				}
				else
				{
					MDeformWeight *dw;
					int w;

					for (dw = dv->dw, w = 0; w < dv->totweight; dw++, w++)
					{
						if (dw->def_nr == inner_defgrp_index) {
							if ((v->flag & ME_VERT_TMP_TAG) && !fallback) {
								dw->weight = 1.0f;
							}
							else if ((p->mat_nr == mat_index-1) && fallback) {
								dw->weight = 1.0f;
							}
							else {
								dw->weight = 0.0f;
							}
						}
					}
				}
			}
		}

		if (dynamic) {
			ssq = rmd->shared->current_shard_entry->prev;
		}

		if (old_cached && ssq) {
			Shard *s;
			int last_id = -1;
			int offset = 0;
			int vertstart = 0;

			for (s = rmd->shared->frac_mesh->shard_map.first; s; s = s->next) {
				MDeformVert *old_dv, *dv;
				int i = 0;

				if (s->shard_id != last_id + 1) {
					Shard *t = BKE_fracture_shard_find(&ssq->frac_mesh->shard_map, last_id + 1);
					if (t) {
						offset += t->totvert;
						printf("Shard offset %d %d\n", t->shard_id, offset);
					}
				}

				for (i = 0; i < s->totvert; i++) {
					if ((vertstart + i + offset) < old_totvert)
					{
						old_dv = old_dvert + vertstart + i + offset;
						dv = dvert + vertstart + i;
						if (old_dv->dw && old_dv->dw->def_nr == inner_defgrp_index) {
							if (dv->dw && dv->dw->def_nr == inner_defgrp_index) {
								dv->dw->weight = old_dv->dw->weight;
							}
						}
					}
				}

				last_id = s->shard_id;
				vertstart += s->totvert;
			}
		}
	}
}
#endif

/* inline face center calc here */
void BKE_fracture_face_calc_center_mean(Mesh *dm, MPoly *mp, float r_cent[3])
{
	MLoop *ml = NULL;
	MLoop *mloop = dm->mloop;
	MVert *mvert = dm->mvert;
	int i = 0;

	zero_v3(r_cent);

	for (i = mp->loopstart; i < mp->loopstart + mp->totloop; i++) {
		MVert *mv = NULL;
		ml = mloop + i;
		mv = mvert + ml->v;

		add_v3_v3(r_cent, mv->co);

	}

	mul_v3_fl(r_cent, 1.0f / (float) mp->totloop);
}

#if 0
TODO FIX
static void do_verts_weights(FractureModifierData *fmd, MeshIsland *mi, int vertstart,
							 int thresh_defgrp_index, int ground_defgrp_index)
{
	MVert *mverts;
	int k;
	MDeformVert *dvert = fmd->shared->dm->dvert;

	mi->vertices_cached = MEM_mallocN(sizeof(MVert *) * s->totvert, "vert_cache");
	mverts = fmd->shared->visible_mesh_cached->mvert;

	mi->vertex_indices = MEM_mallocN(sizeof(int) * mi->vertex_count, "mi->vertex_indices");

	for (k = 0; k < s->totvert; k++) {
		mi->vertices_cached[k] = mverts + vertstart + k;
		mi->vertex_indices[k] = vertstart + k;
		/* sum up vertexweights and divide by vertcount to get islandweight*/
		if (dvert && fmd->thresh_defgrp_name[0]) {
			float vweight = defvert_find_weight(dvert + vertstart + k, thresh_defgrp_index);
			mi->thresh_weight += vweight;
		}

		if (dvert && fmd->ground_defgrp_name[0]) {
			float gweight = defvert_find_weight(dvert + vertstart + k, ground_defgrp_index);
			mi->ground_weight += gweight;
		}
	}

	if (mi->vertex_count > 0) {
		mi->thresh_weight /= mi->vertex_count;
		mi->ground_weight /= mi->vertex_count;
	}
}

#endif

#if 0 //TODO FIX (for dynamic ?)
static void do_handle_parent_mi(FractureModifierData *fmd, MeshIsland *mi, MeshIsland *par, Object* ob, int frame,
								bool is_parent, Scene* scene)
{
	frame -= par->start_frame;
	BKE_match_vertex_coords(mi, par, ob, frame, is_parent, fmd->shards_to_islands);
	if (!is_parent && fmd->is_dynamic_external) {
		//keep the damn names...
		BLI_snprintf(mi->name, sizeof(par->name), "%s", par->name);
	}

	BKE_rigidbody_remove_shard(scene, par);
	scene->rigidbody_world->flag |= RBW_FLAG_OBJECT_CHANGED;
	par->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
}
#endif

#if 0
static bool contains(float loc[3], float size[3], float point[3])
{
	if ((fabsf(loc[0] - point[0]) < size[0]) &&
		(fabsf(loc[1] - point[1]) < size[1]) &&
		(fabsf(loc[2] - point[2]) < size[2]))
	{
		return true;
	}

	return false;
}
#endif

#if 0
void set_rigidbody_type(FractureModifierData *fmd, Shard *s, MeshIsland *mi)
{
	//how far is impact location away from this shard, if beyond a bbox, keep passive
	if (fmd->current_shard_entry)
	{
		ShardSequence *prev_shards = fmd->current_shard_entry->prev;

		if (prev_shards && (prev_shards->prev == NULL)) //only affect primary fracture
		{
			Shard *par_shard = BKE_shard_by_id(prev_shards->frac_mesh, s->parent_id, NULL);
			if (par_shard)
			{
				float impact_loc[3], impact_size[3];
				copy_v3_v3(impact_loc, par_shard->impact_loc);
				copy_v3_v3(impact_size, par_shard->impact_size);

				if (contains(impact_loc, impact_size, s->centroid))
				{
					mi->rigidbody->flag &= ~RBO_FLAG_KINEMATIC;
				}
				else
				{
					mi->rigidbody->flag |= RBO_FLAG_KINEMATIC;
				}

				mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
			}
		}
	}
}
#endif

#if 0
//TODO FIX
static void do_island_from_shard(FractureModifierData *fmd, Object *ob, Shard* s, Mesh *orig_dm,
								 int i, int thresh_defgrp_index, int ground_defgrp_index, int vertstart,
								 Scene *scene, Main *bmain)
{
	MeshIsland *mi;
	MeshIsland *par = NULL;
	bool is_parent = false;
	short rb_type = RBO_TYPE_ACTIVE;
	//float dummyloc[3], rot[4];

	if (s->totvert == 0) {
		return;
	}

	fmd->shared->frac_mesh->progress_counter++;

	mi = MEM_callocN(sizeof(MeshIsland), "meshIsland");
	BLI_addtail(&fmd->shared->mesh_islands, mi);
	unit_qt(mi->rot);

	if (fmd->fracture_mode != MOD_FRACTURE_DYNAMIC)
	{
		mi->locs = NULL;
		mi->rots = NULL;
		mi->frame_count = 0;

		if (scene->rigidbody_world)
		{
			mi->start_frame = scene->rigidbody_world->shared->pointcache->startframe;
		}
		else
		{
			mi->start_frame = 1;
		}
	}
	else
	{
		/* in dynamic case preallocate cache here */
		int start = 1;
		int end = 10;

		if (scene->rigidbody_world)
		{
			start = scene->rigidbody_world->shared->pointcache->startframe;
		}

		if (fmd->shared->current_mi_entry) {
			MeshIslandSequence *prev = fmd->shared->current_mi_entry->prev;
			if (prev)
			{
				start = prev->frame + 1;
			}
		}

		end = start + 10;

		mi->frame_count = end - start + 1;
		mi->start_frame = start;
		mi->locs = MEM_mallocN(sizeof(float)*3* mi->frame_count, "mi->locs");
		mi->rots = MEM_mallocN(sizeof(float)*4* mi->frame_count, "mi->rots");
	}

	mi->participating_constraints = NULL;
	mi->participating_constraint_count = 0;
	mi->thresh_weight = 0.0f;
	mi->ground_weight = 0.0f;
	mi->vertex_count = s->totvert;

	do_verts_weights(fmd, s, mi, vertstart, thresh_defgrp_index, ground_defgrp_index);

	/*copy fixed normals to physics mesh too (needed for convert to objects)*/

	BKE_fracture_physics_normals_fix(fmd, s, mi, i, orig_dm);

	BKE_shard_calc_minmax(s);
	copy_v3_v3(mi->centroid, s->centroid);

	//mat4_to_loc_quat(dummyloc, rot, ob->obmat);
	//copy_qt_qt(mi->rot, rot);
	//unit_qt(mi->rot);
	mi->id = s->shard_id;
	BLI_snprintf(mi->name, 64, "%d", mi->id);

	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		/*take care of previous transformation, if any*/
		MeshIslandSequence *prev = NULL;
		int val = fmd->shards_to_islands ? -1 : 0;

		if (fmd->shared->current_mi_entry) {
			prev = fmd->shared->current_mi_entry->prev;
		}

		/*also take over the UNFRACTURED last shards transformation !!! */
		if (s->parent_id == val)
		{
			//float quat[4];

			//TODO, scale ?
			//mat4_to_quat(quat, ob->obmat);

			mi->locs[0] = mi->centroid[0];
			mi->locs[1] = mi->centroid[1];
			mi->locs[2] = mi->centroid[2];

			mi->rots[0] = mi->rot[0];
			mi->rots[1] = mi->rot[1];
			mi->rots[2] = mi->rot[2];
			mi->rots[3] = mi->rot[3];
		}

		if (prev)
		{
			int frame = prev->frame;

			par = find_meshisland(&prev->meshIslands, s->parent_id);
			if (par)
			{
				is_parent = true;
				do_handle_parent_mi(fmd, mi, par, ob, frame, is_parent, scene);
			}
			else
			{
				par = find_meshisland(&prev->meshIslands, s->shard_id);
				if (par)
				{
					is_parent = false;
					do_handle_parent_mi(fmd, mi, par, ob, frame, is_parent, scene);
				}
			}
		}
	}

	mi->bb = BKE_boundbox_alloc_unit();
	BKE_boundbox_init_from_minmax(mi->bb, s->min, s->max);

	mi->particle_index = -1;
	mi->neighbor_ids = s->neighbor_ids;
	mi->neighbor_count = s->neighbor_count;

	rb_type = do_vert_index_map(fmd, mi, par);
	do_rigidbody(bmain, scene, mi, DEG_get_original_object(ob), orig_dm, rb_type, i);

	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		if (par != NULL)
		{
			int val = fmd->shards_to_islands ? -1 : 0;
			copy_v3_v3(mi->rigidbody->lin_vel, par->rigidbody->lin_vel);
			copy_v3_v3(mi->rigidbody->ang_vel, par->rigidbody->ang_vel);
			mi->rigidbody->flag = par->rigidbody->flag;

			//keep 1st level shards kinematic if parent is triggered
			if ((par->rigidbody->flag & RBO_FLAG_USE_KINEMATIC_DEACTIVATION) && fmd->limit_impact && !fmd->is_dynamic_external) {

				ShardSequence *prev_shards = fmd->shared->current_shard_entry ? fmd->shared->current_shard_entry->prev : NULL;
				Shard *par_shard = prev_shards ? BKE_fracture_shard_find(&prev_shards->frac_mesh->shard_map, s->parent_id) : NULL;

				if (!par_shard) {
					par_shard = prev_shards ? BKE_fracture_shard_find(&prev_shards->frac_mesh->shard_map, s->shard_id) : NULL;
				}

				if (par_shard) {
					float size[3];
					copy_v3_v3(size, par_shard->impact_size);
					mul_v3_fl(size, 2.0f);

					if (contains(par_shard->impact_loc, size, mi->rigidbody->pos)) {
						mi->rigidbody->flag &= ~RBO_FLAG_KINEMATIC;
						mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
					}
					else if (par->id == val) {
						mi->rigidbody->flag |= RBO_FLAG_KINEMATIC;
						mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
					}
				}
				else if (par->id > val) {
					mi->rigidbody->flag &= ~RBO_FLAG_KINEMATIC;
					mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
				}
			}
		}

		mi->rigidbody->meshisland_index = mi->id;
	}
}

MDeformVert* BKE_fracture_shards_to_islands(FractureModifierData* fmd, Object* ob, Mesh *orig_dm, Scene *scene)
{
	/* can be created without shards even, when using fracturemethod = NONE (re-using islands)*/
	Shard *s;
	int i = 0, vertstart = 0;

	MDeformVert *ivert = NULL;
	ListBase shardlist;
	const int thresh_defgrp_index = defgroup_name_index(ob, fmd->thresh_defgrp_name);
	const int ground_defgrp_index = defgroup_name_index(ob, fmd->ground_defgrp_name);

	/*XXX should rename this... this marks the fracture case, to distinguish from halving case */
	fmd->valid_mesh = true;

	if (fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED)
	{
		/* exchange cached mesh after fracture, XXX looks like double code */
		if (fmd->shared->visible_mesh_cached) {
			BKE_mesh_free(fmd->shared->visible_mesh_cached);
			fmd->shared->visible_mesh_cached = NULL;
		}

		fmd->shared->visible_mesh_cached = BKE_fracture_mesh_copy(fmd->shared->dm, ob);

		/* to write to a vgroup (inner vgroup) use the copied cached mesh */
		ivert = fmd->shared->visible_mesh_cached->dvert;

		if (ivert == NULL) {    /* add, if not there */
			int totvert = fmd->shared->visible_mesh_cached->totvert;
			ivert = CustomData_add_layer(&fmd->shared->visible_mesh_cached->vdata, CD_MDEFORMVERT, CD_CALLOC,
										 NULL, totvert);
		}
	}
	else
	{
		fmd->shared->visible_mesh_cached = BKE_fracture_mesh_copy(fmd->shared->dm, ob);
	}

	shardlist = fmd->shared->frac_mesh->shard_map;

	for (s = shardlist.first; s; s = s->next) {

		// XXX need some solution for G.main...
		do_island_from_shard(fmd, ob, s, orig_dm, i, thresh_defgrp_index, ground_defgrp_index, vertstart,
							 scene, G.main);
		vertstart += s->totvert;
		i++;
	}

	return ivert;
}

Mesh *BKE_fracture_result_mesh(FractureModifierData* fmd, Mesh *dm, Object* ob, bool exploOK, Scene* scene)
{
	if ((fmd->shared->visible_mesh_cached != NULL) && exploOK) {
		Mesh *dm_final;

		MDeformVert *dvert = fmd->shared->visible_mesh_cached->dvert;

		//fade out weights in dynamic mode
		if (dvert && (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)) {
			int i;
			int defgrp = defgroup_name_index(ob, fmd->inner_defgrp_name);
			for (i = 0; i < fmd->shared->visible_mesh_cached->totvert; i++) {
				if (dvert[i].dw && dvert[i].dw->def_nr == defgrp && dvert[i].dw->weight >= 0.0f) {
					dvert[i].dw->weight -= 0.1f;
				}
			}
		}

		if (fmd->autohide_dist > 0 || fmd->automerge_dist > 0 || fmd->use_centroids || fmd->use_vertices)
		{
			//printf("Autohide \n");
			dm_final = BKE_fracture_autohide_do(fmd, fmd->shared->visible_mesh_cached, ob, scene);
		}
		else {
			dm_final = BKE_fracture_mesh_copy(fmd->shared->visible_mesh_cached, ob);
			if (!fmd->fix_normals) {
			   BKE_mesh_calc_normals(dm_final);
			}
		}

		return dm_final;
	}
	else {
		if (fmd->shared->visible_mesh == NULL && fmd->shared->visible_mesh_cached == NULL) {
			/* oops, something went definitely wrong... */
			fmd->shared->refresh = true;
			BKE_fracture_free(fmd, fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED, true, scene);
			fmd->shared->visible_mesh_cached = NULL;
			fmd->shared->refresh = false;
		}
	}

	return dm;
}
#endif

#if 0
static void do_post_island_creation(FractureModifierData *fmd, Object *ob, Mesh *dm, Scene* scene, Depsgraph *depsgraph)
{
	double start;

	if (((fmd->shared->visible_mesh != NULL && fmd->shared->refresh && (!fmd->valid_mesh)) || (fmd->shared->visible_mesh_cached == NULL))
		&& (fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED))
	{
		start = PIL_check_seconds_timer();
		/*post process ... convert to DerivedMesh only at refresh times, saves permanent conversion during execution */
		if (fmd->shared->visible_mesh_cached != NULL) {
			BKE_mesh_free(fmd->shared->visible_mesh_cached);
			fmd->shared->visible_mesh_cached = NULL;
		}

		fmd->shared->visible_mesh_cached = createCache(fmd, ob, dm);
		printf("Building cached DerivedMesh done, %g\n", PIL_check_seconds_timer() - start);
	}
	else
	{
		/* fallback, this branch is executed when the modifier data has been loaded via readfile.c,
		 * although this might not be directly visible due to complex logic */

		MDeformVert* dvert = NULL;
		if (fmd->shared->visible_mesh_cached) {
			dvert = fmd->shared->visible_mesh_cached->dvert;
			BKE_fracture_fill_vgroup(fmd, fmd->shared->visible_mesh_cached, dvert, ob, NULL);
		}
	}

	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && fmd->shared->refresh == true)
	{
		fmd->shared->current_mi_entry->is_new = false;
	}

	if (fmd->shared->refresh)
	{
		fmd->shared->refresh = false;
		BKE_rigidbody_update_ob_array(scene->rigidbody_world, false);

		DEG_id_tag_update(&ob->id, OB_RECALC_OB);
	}

	fmd->shared->refresh_constraints = true;
	fmd->shared->refresh_autohide = true;
	fmd->distortion_cached = false;
}
#endif

void BKE_fracture_split_islands(FractureModifierData *fmd, Object* ob, Mesh *me, Mesh*** temp_islands, int* count)
{
	double start;
	BMesh *bm_work;

	bm_work = BKE_fracture_mesh_to_bmesh(me);
	start = PIL_check_seconds_timer();
	//printf("Steps: %d \n", fmd->frac_mesh->progress_counter);
	mesh_separate_loose(bm_work, fmd, ob, temp_islands, count);
	printf("Splitting to islands done, %g \n"/*  Steps: %d \n"*/, PIL_check_seconds_timer() - start);//, fmd->frac_mesh->progress_counter);
}

#if 0
static void do_refresh(FractureModifierData *fmd, Object *ob, Mesh* dm, Mesh *orig_dm, Mesh *old_cached, Scene *scene)
{
	double start = 0.0;
	MDeformVert *ivert = NULL;

	copy_m4_m4(fmd->origmat, ob->obmat);

	{
		/* refracture, convert the fracture shards to new meshislands here *
		 * shards = fracture datastructure
		 * meshisland = simulation datastructure */
		if (fmd->shared->frac_mesh && fmd->shared->frac_mesh->shard_count > 0 && fmd->shared->dm && fmd->shared->dm->totvert > 0 &&
			(!fmd->shards_to_islands || fmd->fracture_mode == MOD_FRACTURE_DYNAMIC))
		{
			if (fmd->fix_normals)
			{
				start = PIL_check_seconds_timer();
			}

			ivert = BKE_fracture_shards_to_islands(fmd, ob, orig_dm, scene);

			if (fmd->fix_normals) {
				printf("Fixing normals done, %g\n", PIL_check_seconds_timer() - start);
			}

			BKE_fracture_fill_vgroup(fmd, fmd->shared->visible_mesh_cached, ivert, ob, old_cached);
		}
		else if (fmd->fracture_mode != MOD_FRACTURE_DYNAMIC){
			if (fmd->shared->visible_mesh != NULL) {
				BM_mesh_free(fmd->shared->visible_mesh);
				fmd->shared->visible_mesh = NULL;
			}
			BKE_fracture_do_halving(fmd, ob, dm, orig_dm, false, -1, scene);
			fmd->valid_mesh = false;
		}
	}

	printf("Islands: %d\n", BLI_listbase_count(&fmd->shared->mesh_islands));

	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		/* Grrr, due to stupid design of mine (listbase as value in struct instead of pointer)
		 * we have to synchronize the lists here again */

		/* need to ensure(!) old pointers keep valid, else the whole meshisland concept is broken */
		fmd->shared->current_mi_entry->visible_dm = fmd->shared->visible_mesh_cached;
		fmd->shared->current_mi_entry->meshIslands = fmd->shared->mesh_islands;
	}
}
#endif

static void do_island_index_map(FractureModifierData *fmd, Object* ob)
{
	MeshIsland *mi;

	if (fmd->shared->vertex_island_map) {
		BLI_ghash_free(fmd->shared->vertex_island_map, NULL, NULL);
	}

	fmd->shared->vertex_island_map = BLI_ghash_ptr_new("island_index_map");

	if (fmd->dm_group && fmd->use_constraint_group)
	{
		CollectionObject *go;
		int i,j = 0;
		for (go = fmd->dm_group->gobject.first; go; go = go->next)
		{
			if (go->ob != ob)
			{
				FractureModifierData *fmdi = (FractureModifierData*)modifiers_findByType(go->ob, eModifierType_Fracture);
				if (fmdi)
				{
					for (mi = fmdi->shared->mesh_islands.first; mi; mi = mi->next){
						for (i = 0; i < mi->mesh->totvert; i++)
						{
							if (!BLI_ghash_haskey(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(i + j)))
							{
								BLI_ghash_insert(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(i + j), mi);
							}
						}
						j += mi->mesh->totvert;
					}
				}
			}
		}
	}
	else {
		int i,j = 0;
		for (mi = fmd->shared->mesh_islands.first; mi; mi = mi->next){
			for (i = 0; i < mi->mesh->totvert; i++)
			{
				if (!BLI_ghash_haskey(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(i+j)))
				{
					BLI_ghash_insert(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(i+j), mi);
				}
			}
			j += mi->mesh->totvert;
		}
	}
}

#if 0
Mesh *BKE_fracture_mesh_from_packdata(FractureModifierData *fmd, Mesh *derivedData)
{
	Mesh *dm = NULL;

	/* keep old way of using dynamic external working as well, without interfering with packing */
	if (fmd->shared->pack_storage.first && !fmd->is_dynamic_external)
	{
		dm = BKE_fracture_assemble_mesh_from_islands(fmd, );
	}
	else {
		dm = derivedData;
	}

	return dm;
}
#endif

Mesh* BKE_fracture_mesh_copy(Mesh* source, Object* ob)
{
	Mesh* me = BKE_mesh_new_nomain(source->totvert,
								   source->totedge,
								   0,
								   source->totloop,
								   source->totpoly);

	BKE_mesh_nomain_to_mesh(source, me, ob, CD_MASK_MESH, false);

	return me;
}

Mesh* BKE_fracture_bmesh_to_mesh(BMesh* bm)
{
	struct BMeshToMeshParams bmt = {.calc_object_remap = 0};

	Mesh* me = BKE_mesh_new_nomain(bm->totvert,
								   bm->totedge,
								   0,
								   bm->totloop,
								   bm->totface);


	BM_mesh_bm_to_me(NULL, bm, me, &bmt);

	return me;
}

BMesh* BKE_fracture_mesh_to_bmesh(Mesh* me)
{
	struct BMeshFromMeshParams bmf = {.calc_face_normal = true, };
	struct BMeshCreateParams bmc = {.use_toolflags = true, };
	BMesh* bm;

	const BMAllocTemplate allocsize = BMALLOC_TEMPLATE_FROM_ME(me);
	bm = BM_mesh_create(&allocsize, &bmc);

	BM_mesh_bm_from_me(bm, me, &bmf);

	BM_mesh_elem_toolflags_ensure(bm);
	BM_mesh_elem_index_ensure(bm, BM_VERT | BM_EDGE | BM_FACE);
	BM_mesh_elem_table_ensure(bm, BM_VERT | BM_EDGE | BM_FACE);

	return bm;
}

void BKE_fracture_external_constraints_setup(FractureModifierData *fmd, Scene *scene, Object *ob)
{
	if (fmd->shared->refresh_constraints || fmd->shared->refresh) {

		Object *pob = NULL;
		FractureModifierData *pfmd = NULL;

		fmd->shared->refresh_constraints = false;

		//force-refresh other FMs if they have "us" in our group (shouldnt be 1000s, so should be ok performance wise)
		if (scene && scene->rigidbody_world) {
			CollectionObject *go;
			for (go = scene->rigidbody_world->group->gobject.first; go; go = go->next)
			{
				FractureModifierData *fmdi = (FractureModifierData*)modifiers_findByType(go->ob, eModifierType_Fracture);
				if (fmdi && fmdi->dm_group && fmdi->use_constraint_group)
				{
					CollectionObject *go2;
					for (go2 = fmdi->dm_group->gobject.first; go2; go2 = go2->next)
					{
						if (go2->ob == ob)
						{
							pfmd = fmdi;
							pob = go->ob;

							fmdi->shared->refresh_constraints = true;
							BKE_fracture_constraints_free(fmdi,scene);
							break;
						}
					}
				}
			}
		}

		do_island_index_map(fmd, ob);
		BKE_fracture_constraints_refresh(fmd, ob, scene);

		if (fmd->dm_group && fmd->use_constraint_group)
		{	//disable the carrier object, it would interfere (it should have 1 island only)
			MeshIsland *mi = fmd->shared->mesh_islands.first;
			mi->rigidbody->flag |= RBO_FLAG_KINEMATIC;
			mi->rigidbody->flag |= RBO_FLAG_IS_GHOST;

			ob->rigidbody_object->flag |= RBO_FLAG_KINEMATIC;
			ob->rigidbody_object->flag |= RBO_FLAG_IS_GHOST;
		}

		if (pfmd && pob) {
			double start = PIL_check_seconds_timer();
			do_island_index_map(pfmd, pob);
			BKE_fracture_constraints_refresh(pfmd, pob, scene);
			pfmd->shared->refresh_constraints = false;
			printf("Rebuilding external constraints done, %g\n", PIL_check_seconds_timer() - start);
		}
	}
}

void BKE_fracture_mesh_boundbox_calc(Mesh *me, float r_loc[], float r_size[])
{
	//just calculate loc and size, but throw bbox away again
	BKE_mesh_boundbox_calc(me, r_loc, r_size);
	if (me->bb)
	{
		MEM_freeN(me->bb);
		me->bb = NULL;
	}
}

void BKE_fracture_mesh_free(Mesh *me)
{
	BKE_mesh_free(me);
	MEM_freeN(me);
	me = NULL;
}
