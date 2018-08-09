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
 * Copyright (C) 2014 by Martin Felke.
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


#include "MEM_guardedalloc.h"

#include "BLI_callbacks.h"
#include "BLI_edgehash.h"
#include "BLI_ghash.h"
#include "BLI_kdtree.h"
#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_math_matrix.h"
#include "BLI_math_vector.h"
#include "BLI_rand.h"
#include "BLI_utildefines.h"
#include "BLI_string.h"
#include "BLI_threads.h"


#include "BKE_cdderivedmesh.h"
#include "BKE_deform.h"
#include "BKE_fracture.h"
#include "BKE_global.h"
#include "BKE_collection.h"
#include "BKE_library.h"
#include "BKE_library_query.h"
#include "BKE_main.h"
#include "BKE_material.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_particle.h"
#include "BKE_pointcache.h"
#include "BKE_rigidbody.h"
#include "BKE_scene.h"
#include "BKE_mesh.h"
#include "BKE_curve.h"
#include "BKE_multires.h"

#include "bmesh.h"

#include "DNA_fracture_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_group_types.h"
#include "DNA_listBase.h"
#include "DNA_material_types.h"
#include "DNA_object_types.h"
#include "DNA_particle_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"
#include "DNA_curve_types.h"

#include "../../rigidbody/RBI_api.h"
#include "PIL_time.h"
#include "../../bmesh/tools/bmesh_decimate.h" /* decimate_dissolve function */
#include "limits.h"

/* prototypes */
static Shard *parse_cell(cell c);
static void parse_cell_verts(cell c, MVert *mvert, int totvert);
static void parse_cell_polys(cell c, MPoly *mpoly, int totpoly, int *r_totloop);
static void parse_cell_loops(cell c, MLoop *mloop, int totloop, MPoly *mpoly, int totpoly);
static void parse_cell_neighbors(cell c, int *neighbors, int totpoly);
static void arrange_shard(FractureModifierData *fmd, ShardID id, bool do_verts, float cent[3]);
static void fracture_shard_add(FracMesh *fm, Shard *s, float mat[4][4]);
static Shard* fracture_initial_shard_create(Mesh *dm);
static void fracture_shard_add(FracMesh *fm, Shard *s, float mat[4][4])
{
	MVert *mv;
	int i = 0;

	for (i = 0, mv = s->mvert; i < s->totvert; i++, mv++ )
	{
		mul_m4_v3(mat, mv->co);
	}

	mul_m4_v3(mat, s->centroid);

	BLI_addtail(&fm->shard_map, s);
	s->shard_id = fm->shard_count;
	fm->shard_count++;
}

static BMesh* fracture_shard_to_bmesh(Shard *s)
{
	Mesh *dm_parent;
	BMesh *bm_parent;
	BMIter iter;
	BMFace *f;

	dm_parent = BKE_fracture_shard_to_mesh(s, true);
	bm_parent = BKE_fracture_mesh_to_bmesh(dm_parent);

	BM_mesh_elem_table_ensure(bm_parent, BM_VERT | BM_FACE);

	BM_ITER_MESH (f, &iter, bm_parent, BM_FACES_OF_MESH)
	{
		BM_elem_flag_disable(f, BM_ELEM_SELECT);
	}

	BKE_mesh_free(dm_parent);
	dm_parent = NULL;

	return bm_parent;
}

static void fracture_shard_boundbox(Shard *s, float r_loc[3], float r_size[3])
{
	float min[3], max[3];
	float mloc[3], msize[3];

	if (!r_loc) r_loc = mloc;
	if (!r_size) r_size = msize;

	if (!BKE_shard_calc_minmax(s)) {
		min[0] = min[1] = min[2] = -1.0f;
		max[0] = max[1] = max[2] = 1.0f;
	}

	copy_v3_v3(max, s->max);
	copy_v3_v3(min, s->min);

	mid_v3_v3v3(r_loc, min, max);

	r_size[0] = (max[0] - min[0]) / 2.0f;
	r_size[1] = (max[1] - min[1]) / 2.0f;
	r_size[2] = (max[2] - min[2]) / 2.0f;
}

static int shard_sortsize(const void *s1, const void *s2, void* UNUSED(context))
{
	Shard **sh1 = (Shard **)s1;
	Shard **sh2 = (Shard **)s2;

	float size1[3], size2[3], loc[3];
	float val_a,  val_b;

	if ((*sh1 == NULL) || (*sh2 == NULL)) {
		return -1;
	}

	fracture_shard_boundbox(*sh1, loc, size1);
	fracture_shard_boundbox(*sh2, loc, size2);

	//squared diameter
	val_a = size1[0]*size1[0] + size1[1]*size1[1] + size1[2]*size1[2];
	val_b = size2[0]*size2[0] + size2[1]*size2[1] + size2[2]*size2[2];

	/* sort descending */
	if      (val_a < val_b) return 1;
	else if (val_a > val_b) return -1;
	return 0;
}

static void* check_add_layer(CustomData *src, CustomData *dst, int type, int totelem, const char* name)
{
	void *layer =  CustomData_get_layer_named(dst, type, name);

	if (!layer) {
		void* orig = NULL;

		if (src) {
			orig = CustomData_get_layer_named(src, type, name);
		}

		if (orig) {
			return CustomData_add_layer_named(dst, type, CD_DUPLICATE, orig, totelem, name);
		}
		else{
			return CustomData_add_layer_named(dst, type, CD_CALLOC, NULL, totelem, name);
		}
	}
	else {
		return layer;
	}
}

void BKE_fracture_custom_data_mesh_to_shard(Shard *s, Mesh *dm)
{
	CustomData_reset(&s->vertData);
	CustomData_reset(&s->loopData);
	CustomData_reset(&s->polyData);
	CustomData_reset(&s->edgeData);

	CustomData_copy(&dm->vdata, &s->vertData, CD_MASK_MDEFORMVERT, CD_DUPLICATE, s->totvert);
	CustomData_copy(&dm->ldata, &s->loopData, CD_MASK_MLOOPUV, CD_DUPLICATE, s->totloop);
	CustomData_copy(&dm->edata, &s->edgeData, CD_MASK_CREASE | CD_MASK_BWEIGHT | CD_MASK_MEDGE, CD_DUPLICATE, s->totedge);

	//add velocity vertex layers...
	check_add_layer(&dm->vdata, &s->vertData, CD_PROP_FLT, s->totvert, "velX");
	check_add_layer(&dm->vdata, &s->vertData, CD_PROP_FLT, s->totvert, "velY");
	check_add_layer(&dm->vdata, &s->vertData, CD_PROP_FLT, s->totvert, "velZ");
}

/* modified from BKE_mesh_center_median */
bool BKE_fracture_shard_center_median(Shard *shard, float cent[3])
{
	int i = shard->totvert;
	MVert *mvert;
	zero_v3(cent);
	for (mvert = shard->mvert; i--; mvert++) {
		add_v3_v3(cent, mvert->co);
	}
	/* otherwise we get NAN for 0 verts */
	if (shard->totvert) {
		mul_v3_fl(cent, 1.0f / (float)shard->totvert);
	}

	return (shard->totvert != 0);
}

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

/* modified from BKE_mesh_center_centroid */
bool BKE_fracture_shard_center_centroid(Shard *shard, float r_cent[3])
{
	int i = shard->totpoly;
	MPoly *mpoly;
	float poly_volume;
	float total_volume = 0.0f;
	float poly_cent[3];

	zero_v3(r_cent);

	/* calculate a weighted average of polyhedron centroids */
	for (mpoly = shard->mpoly; i--; mpoly++) {
		poly_volume = mesh_calc_poly_volume_and_weighted_centroid(mpoly, shard->mloop + mpoly->loopstart, shard->mvert, poly_cent);

		/* poly_cent is already volume-weighted, so no need to multiply by the volume */
		add_v3_v3(r_cent, poly_cent);
		total_volume += poly_volume;
	}
	/* otherwise we get NAN for 0 polys */
	if (total_volume != 0.0f) {
		/* multipy by 0.25 to get the correct centroid */
		/* no need to divide volume by 6 as the centroid is weighted by 6x the volume, so it all cancels out */
		mul_v3_fl(r_cent, 0.25f / total_volume);
	}

	/* this can happen for non-manifold objects, first fallback to old area based method, then fallback to median there */
	if (!is_finite_v3(r_cent) || total_volume < 0.000001f) {
		return BKE_fracture_shard_center_centroid_area(shard, r_cent);
	}

	copy_v3_v3(shard->centroid, r_cent);
	return (shard->totpoly != 0);
}

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
bool BKE_fracture_shard_center_centroid_area(Shard *shard, float cent[3])
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
		return BKE_fracture_shard_center_median(shard, cent);
	}
	copy_v3_v3(shard->centroid, cent);

	return (shard->totpoly != 0);
}

void BKE_fracture_shard_free(Shard *s, bool doCustomData)
{
	if (s->mvert) {
		MEM_freeN(s->mvert);
	}
	if (s->mloop) {
		MEM_freeN(s->mloop);
	}
	if (s->mpoly) {
		MEM_freeN(s->mpoly);
	}
	if (s->medge) {
		MEM_freeN(s->medge);
	}
	if (s->neighbor_ids) {
		MEM_freeN(s->neighbor_ids);
	}
	if (s->cluster_colors) {
		MEM_freeN(s->cluster_colors);
	}

	if (doCustomData) {
		CustomData_free(&s->vertData, s->totvert);
		CustomData_free(&s->loopData, s->totloop);
		CustomData_free(&s->polyData, s->totpoly);
		CustomData_free(&s->edgeData, s->totedge);
	}

	MEM_freeN(s);
}

float BKE_shard_calc_minmax(Shard *shard)
{
	float min[3], max[3], diff[3];
	int i;

	INIT_MINMAX(min, max);
	for (i = 0; i < shard->totvert; i++) {
		minmax_v3v3_v3(min, max, shard->mvert[i].co);
	}

	copy_v3_v3(shard->min, min);
	copy_v3_v3(shard->max, max);

	sub_v3_v3v3(diff, max, min);
	return len_v3(diff);
}

static Shard* fracture_initial_shard_create(Mesh *dm)
{
	/* create temporary shard covering the entire mesh */
	Shard *s = BKE_fracture_shard_create(dm->mvert, dm->mpoly, dm->mloop, dm->medge,
										 dm->totvert, dm->totpoly, dm->totloop, dm->totedge, true);
	BKE_fracture_custom_data_mesh_to_shard(s, dm);
	s->flag = SHARD_INTACT;
	s->shard_id = -2;
	return s;
}

/*access shard directly by index / id*/
Shard *BKE_shard_by_id(FracMesh *mesh, ShardID id, Mesh *dm) {
	if ((id >= 0)) {
		Shard *s = mesh->shard_map.first;
		while (s)
		{
			if (s->shard_id == id)
			{
				return s;
			}
			s = s->next;
		}

		return NULL;
	}
	else if (id == -1 && dm != NULL)
	{
		/* create temporary shard covering the entire mesh */
		return fracture_initial_shard_create(dm);
	}

	return NULL;
}

bool BKE_get_shard_minmax(FracMesh *mesh, ShardID id, float min_r[3], float max_r[3], Mesh *dm)
{
	Shard *shard = BKE_shard_by_id(mesh, id, dm);
	if (shard != NULL) {
		BKE_shard_calc_minmax(shard);
		copy_v3_v3(min_r, shard->min);
		copy_v3_v3(max_r, shard->max);

		if (shard->shard_id == -2)
		{
			BKE_fracture_shard_free(shard, true);
		}

		return true;
	}
	return false;
}

Shard *BKE_fracture_shard_create(MVert *mvert, MPoly *mpoly, MLoop *mloop, MEdge* medge,  int totvert, int totpoly,
								 int totloop, int totedge, bool copy)
{
	Shard *shard = MEM_mallocN(sizeof(Shard), __func__);
	shard->totvert = totvert;
	shard->totpoly = totpoly;
	shard->totloop = totloop;
	shard->totedge = totedge;
	shard->cluster_colors = NULL;
	shard->neighbor_ids = NULL;
	shard->neighbor_count = 0;

	if (copy) {
		shard->mvert = MEM_mallocN(sizeof(MVert) * totvert, "shard vertices");
		shard->mpoly = MEM_mallocN(sizeof(MPoly) * totpoly, "shard polys");
		shard->mloop = MEM_mallocN(sizeof(MLoop) * totloop, "shard loops");
		shard->medge = MEM_mallocN(sizeof(MEdge) * totedge, "shard edges");
		memcpy(shard->mvert, mvert, sizeof(MVert) * totvert);
		memcpy(shard->mpoly, mpoly, sizeof(MPoly) * totpoly);
		memcpy(shard->mloop, mloop, sizeof(MLoop) * totloop);
		memcpy(shard->medge, medge, sizeof(MEdge) * totedge);
	}
	else {
		shard->mvert = mvert;
		shard->mpoly = mpoly;
		shard->mloop = mloop;
		shard->medge = medge;
	}

	shard->shard_id = -1;
	shard->setting_id = -1;
	shard->parent_id = -1;

	shard->flag = SHARD_INTACT;
	BKE_shard_calc_minmax(shard);

	BKE_fracture_shard_center_centroid_area(shard, shard->centroid);
	copy_v3_v3(shard->raw_centroid, shard->centroid);
	zero_v3(shard->impact_loc);
	shard->impact_size[0] = 1.0f;
	shard->impact_size[1] = 1.0f;
	shard->impact_size[2] = 1.0f;

	return shard;
}

FracMesh *BKE_fracture_fracmesh_create(void)
{
	FracMesh *fmesh;

	fmesh = MEM_mallocN(sizeof(FracMesh), __func__);
	fmesh->shard_map.first = NULL;
	fmesh->shard_map.last = NULL;
	fmesh->shard_count = 0;
	fmesh->cancel = 0;
	fmesh->running = 0;
	fmesh->progress_counter = 0;
	fmesh->last_shards = NULL;
	fmesh->last_shard_tree = NULL;
	fmesh->last_expected_shards = 0;

	return fmesh;
}

static void handle_fast_bisect(FracMesh *fm, int expected_shards, int algorithm, BMesh** bm_parent, float obmat[4][4],
							   float centroid[3], short inner_material_index, int parent_id, Shard **tempshards, Shard ***tempresults,
							   char uv_layer[64], cell* cells, float fac, Shard* parent)
{
	int i = 0, index = 0;
	float factor = 1 - fac;

	float dim[3];
	sub_v3_v3v3(dim, parent->max, parent->min);

	for (i = 0; i < expected_shards; i++) {
		Shard *s = NULL;
		Shard *s2 = NULL;
		Shard *t;
		float vec[3];
		int max_axis;

		if (fm->cancel == 1) {
			break;
		}

		printf("Processing shard: %d\n", i);
		t = tempshards[i];

		if (t != NULL) {
			t->parent_id = parent_id;
			t->flag = SHARD_INTACT;
		}

		if (t == NULL || t->totvert == 0 || t->totloop == 0 || t->totpoly == 0) {
			/* invalid shard, stop parsing*/
			continue;
		}

		//index = (int)(BLI_frand() * (t->totpoly - 1));

		if (index > (t->totpoly - 1)){
			index = 0;
		}

		//make a random vector (interpret as cutter plane)
		vec[0] = BLI_thread_frand(0) * 2 - 1;
		vec[1] = BLI_thread_frand(0) * 2 - 1;
		vec[2] = BLI_thread_frand(0) * 2 - 1;

		//multiply two minor dimensions with a factor to emphasize the max dimension
		max_axis = axis_dominant_v3_single(dim);
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

		printf("Bisecting cell %d...\n", i);
		printf("Bisecting cell %d...\n", i + 1);

		s = BKE_fracture_shard_bisect(*bm_parent, t, obmat, algorithm == MOD_FRACTURE_BISECT_FAST_FILL,
									  false, true, index, centroid, inner_material_index, uv_layer, NULL, vec);
		s2 = BKE_fracture_shard_bisect(*bm_parent, t, obmat, algorithm == MOD_FRACTURE_BISECT_FAST_FILL,
									   true, false, index, centroid, inner_material_index, uv_layer, NULL, vec);

		index++;

		if (s == NULL || s2 == NULL) {
			printf("Shard missed....\n");
			continue;
		}

		if (s != NULL && s2 != NULL && tempresults != NULL) {
			int j = 0;

			fm->progress_counter++;

			s->parent_id = parent_id;
			s->flag = SHARD_INTACT;

			s2->parent_id = parent_id;
			s2->flag = SHARD_INTACT;

			if (*bm_parent != NULL) {
				BM_mesh_free(*bm_parent);
				*bm_parent = NULL;
			}

			(*tempresults)[i] = s;
			(*tempresults)[i + 1] = s2;

			//BLI_qsort_r(*tempresults, i + 1, sizeof(Shard *), shard_sortdist, &(cells[i]));
			BLI_qsort_r(*tempresults, i + 1, sizeof(Shard *), shard_sortsize, &(cells[i]));

			while ((*tempresults)[j] == NULL && j < (i + 1)) {
				/* ignore invalid shards */
				j++;
			}

			/* continue splitting if not all expected shards exist yet */
			if ((i + 2) < expected_shards) {
				*bm_parent = fracture_shard_to_bmesh((*tempresults)[j]);
				copy_v3_v3(centroid, (*tempresults)[j]->centroid);
				sub_v3_v3v3(dim, (*tempresults)[j]->max, (*tempresults)[j]->min);

				BKE_fracture_shard_free((*tempresults)[j], true);
				(*tempresults)[j] = NULL;
			}
			i++;
		}
	}
}

static void handle_boolean_fractal(Shard* t, int expected_shards, Mesh* dm_parent, Object *obj, short inner_material_index,
								   int num_cuts, float fractal, int num_levels, bool smooth,int parent_id, int* i, Shard ***tempresults,
								   Mesh **dm_p, char uv_layer[64], int thresh, float fac)
{
	/* physics shard and fractalized shard, so we need to booleanize twice */
	/* and we need both halves, so twice again */
	Shard *s2 = NULL;
	Shard *s = NULL;
	int index = 0;
	int max_retries = 3;
	float factor = 1 - fac;

	/*continue with "halves", randomly*/
	if ((*i) == 0) {
		*dm_p = BKE_fracture_mesh_copy(dm_parent, obj);
	}

	while (s == NULL || s2 == NULL) {

		float radius;
		float size[3];
		float quat[4];
		float loc[3], vec[3];
		float min[3], max[3];
		float one[3] = {1.0f, 1.0f, 1.0f};
		float matrix[4][4];
		int max_axis;

		/*make a plane as cutter*/
//		shard_boundbox(p, loc, size);
		INIT_MINMAX(min, max);
		BKE_mesh_minmax(*dm_p, min, max);

		mid_v3_v3v3(loc, min, max);
		size[0] = (max[0] - min[0]) / 2.0f;
		size[1] = (max[1] - min[1]) / 2.0f;
		size[2] = (max[2] - min[2]) / 2.0f;

		radius = sqrt(size[0]*size[0] + size[1]*size[1] + size[2]*size[2]);

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

		/*visual shards next, fractalized cuts */
		s = BKE_fracture_shard_boolean(obj, *dm_p, t, inner_material_index, num_cuts,fractal, &s2, matrix, radius, smooth,
									   num_levels, uv_layer, thresh);

		if (index < max_retries)
		{
			printf("Retrying...%d\n", index);
			index++;
		}
		else if (s == NULL || s2 == NULL)
		{
			(*i)++;
			break;
		}
	}

	if ((s != NULL) && (s2 != NULL)) {
		int j = 0;

		s->parent_id = parent_id;
		s->flag = SHARD_INTACT;
		(*tempresults)[(*i)+1] = s;

		s2->parent_id = parent_id;
		s2->flag = SHARD_INTACT;
		(*tempresults)[*i] = s2;

		BLI_qsort_r(*tempresults, (*i) + 1, sizeof(Shard *), shard_sortsize, i);
		while ((*tempresults)[j] == NULL && j < ((*i) + 1)) {
			/* ignore invalid shards */
			j++;
		}

		/* continue splitting if not all expected shards exist yet */
		if (((*i) + 2) < expected_shards) {

			Shard *p = (*tempresults)[j];

			if (*dm_p != NULL) {
				BKE_mesh_free(*dm_p);
				*dm_p = NULL;
			}

			if (p != NULL) {
				*dm_p = BKE_fracture_shard_to_mesh(p, true);
				BKE_fracture_shard_free((*tempresults)[j], true);
				(*tempresults)[j] = NULL;
			}
		}
		(*i)++; //XXX remember to "double" the shard amount....
	}
}

static bool handle_boolean_bisect(FracMesh *fm, Object *obj, int expected_shards, int algorithm, int parent_id, Shard **tempshards,
								  Mesh *dm_parent, BMesh* bm_parent, float obmat[4][4], short inner_material_index, int num_cuts,
								  int num_levels, float fractal, int *i, bool smooth, Shard*** tempresults, Mesh **dm_p, char uv_layer[64],
								  KDTree *preselect_tree, int thresh, float fac)
{
	Shard *s = NULL, *t = NULL;
	if (fm->cancel == 1)
		return true;

	t = tempshards[*i];

	if (t != NULL) {
		t->parent_id = parent_id;
		t->flag = SHARD_INTACT;
	}

	if (t == NULL || t->totvert == 0 || t->totloop == 0 || t->totpoly == 0) {
		/* invalid shard, stop parsing */
		return true;
	}

	printf("Processing shard: %d\n", *i);

	/* XXX TODO, need object for material as well, or atleast a material index... */
	if (algorithm == MOD_FRACTURE_BOOLEAN) {
		s = BKE_fracture_shard_boolean(obj, dm_parent, t, inner_material_index, 0, 0.0f, NULL, NULL, 0.0f, false, 0, uv_layer, thresh);
	}
	else if (algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL) {
		handle_boolean_fractal(t, expected_shards, dm_parent, obj, inner_material_index, num_cuts, fractal,
							   num_levels, smooth, parent_id, i, tempresults, dm_p, uv_layer, thresh, fac);
	}
	else if (algorithm == MOD_FRACTURE_BISECT || algorithm == MOD_FRACTURE_BISECT_FILL) {
		float co[3] = {0, 0, 0}, quat[4] =  {1, 0, 0, 0};
		printf("Bisecting cell %d...\n", *i);
		s = BKE_fracture_shard_bisect(bm_parent, t, obmat, algorithm == MOD_FRACTURE_BISECT_FILL, false, true, -1, co, inner_material_index, uv_layer,
									  preselect_tree, quat);
	}
	else {
		/* do not fracture case */
		s = t;
	}

	if ((s != NULL) && (algorithm != MOD_FRACTURE_BOOLEAN_FRACTAL)) {
		s->parent_id = parent_id;
		s->flag = SHARD_INTACT;

		(*tempresults)[*i] = s;
	}

	fm->progress_counter++;
	return false;
}

static void do_prepare_cells(FracMesh *fm, cell *cells, int expected_shards, int algorithm, Shard *p, float (*centroid)[3],
							 Mesh **dm_parent, BMesh** bm_parent, Shard ***tempshards, Shard ***tempresults, int override_count,
							 FractureModifierData *fmd)
{
	int i;
	Shard *s = NULL;
	int *skipmap = MEM_callocN(sizeof(int) * expected_shards, "skipmap");
	int *deletemap = MEM_callocN(sizeof(int) * fm->shard_count, "deletemap");

	if ((algorithm == MOD_FRACTURE_BOOLEAN) || (algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL)) {
		MPoly *mpoly, *mp;
		int totpoly, po;

		*dm_parent = BKE_fracture_shard_to_mesh(p, true);
		mpoly = (*dm_parent)->mpoly;
		totpoly = (*dm_parent)->totpoly;
		for (po = 0, mp = mpoly; po < totpoly; po++, mp++) {
			mp->flag &= ~ME_FACE_SEL;
		}
	}
	else if (algorithm == MOD_FRACTURE_BISECT || algorithm == MOD_FRACTURE_BISECT_FILL ||
			 algorithm == MOD_FRACTURE_BISECT_FAST || algorithm == MOD_FRACTURE_BISECT_FAST_FILL)
	{
		*bm_parent = fracture_shard_to_bmesh(p);
		copy_v3_v3(*centroid, p->centroid);
	}

	if (algorithm == MOD_FRACTURE_BISECT_FAST ||
		algorithm == MOD_FRACTURE_BISECT_FAST_FILL ||
		algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL)
	{
		copy_vn_i(deletemap, fm->shard_count, 1);
	}

	if (fm->last_shard_tree)
	{
		copy_vn_i(deletemap, fm->shard_count, 1);
		copy_vn_i(skipmap, expected_shards, 1);

		for (i = 0; i < expected_shards; i++)
		{
			KDTreeNearest n;
			int l, j;
			float max = 0;
			for (l = 0; l < cells[i].totpoly; l++)
			{
				int index = cells[i].neighbors[l];
				if (index > -1)
				{
					float dist = len_squared_v3v3(cells[index].centroid, cells[i].centroid);
					if (dist > max)
					{
						max = dist;
					}
				}
			}

			j = BLI_kdtree_find_nearest(fm->last_shard_tree, cells[i].centroid, &n);
			if (j > -1)
			{
				float epsilon = 0.00001;
				Shard *t = fm->last_shards[j];
				if (t != NULL && n.dist < max)
				{
					if (n.dist < epsilon) {
						if ((fabsf(cells[i].volume - t->raw_volume) < epsilon))
						{
							deletemap[j] = false;
						}
						else
						{
							skipmap[i] = false;
						}
					}
					else
					{
						skipmap[i] = false;
					}
				}
				else
				{
					skipmap[i] = false;
				}
			}
			else {
				skipmap[i] = true;
			}
		}
	}

	//BLI_lock_thread(LOCK_CUSTOM1);
	//skipping /deletion pass


	for (i = 0; i < expected_shards; i++)
	{
		if (fm->cancel == 1) {
			break;
		}

		if (skipmap[i])
		{
			printf("Skipping shard: %d\n", i);
			(*tempshards)[i] = NULL;
		}
		else
		{
			printf("Parsing shard: %d\n", i);
			s = parse_cell(cells[i]);
			(*tempshards)[i] = s;
		}

		(*tempresults)[i] = NULL;

		fm->progress_counter++;
	}

	for (i = 0; i < fm->shard_count; i++)
	{
		if (deletemap[i] && fm->last_shards)
		{
			Shard *t = fm->last_shards[i];

			if (!t)
				continue;

			//seems this override count was a totally wrong thought, just passing -1 or 0 here... hmm
			if ((override_count == -1) || ((override_count > 0) && (i < override_count+1)))
			{
				printf("Deleting shard: %d %d %d\n", i, t->shard_id, t->setting_id);
				BLI_remlink_safe(&fm->shard_map, t);
				BKE_fracture_shard_free(t, false);
				fm->last_shards[i] = NULL;
			}
			else
			{
				printf("NOT Deleting shard: %d %d %d\n", i, t->shard_id, t->setting_id);
			}
		}
	}

	if (override_count > -1) {
		printf("Deleting island shards!\n");
		while (fmd->shared->islandShards.first) {
			Shard *sh = fmd->shared->islandShards.first;
			if (sh) {
				BLI_remlink_safe(&fmd->shared->islandShards, sh);
				BKE_fracture_shard_free(sh, false);
			}
		}
	}
	//BLI_unlock_thread(LOCK_CUSTOM1);

	fm->last_expected_shards = expected_shards;

	MEM_freeN(skipmap);
	MEM_freeN(deletemap);
}


/* parse the voro++ cell data */
//static ThreadMutex prep_lock = BLI_MUTEX_INITIALIZER;
static void parse_cells(cell *cells, int expected_shards, ShardID parent_id, FracMesh *fm, int algorithm, Object *obj, Mesh *dm,
						short inner_material_index, float mat[4][4], int num_cuts, float fractal, bool smooth, int num_levels, int mode,
						bool reset, char uv_layer[64], bool threaded, float thresh, int override_count, float factor)
{
	/*Parse voronoi raw data*/
	int i = 0, j = 0, count = 0;
	Shard *p = BKE_shard_by_id(fm, parent_id, dm); // *t;
	float obmat[4][4]; /* use unit matrix for now */
	float centroid[3], pcentroid[3] = {0,0,0};
	BMesh *bm_parent = NULL;
	Mesh *dm_parent = NULL;
	Mesh *dm_p = NULL;
	Shard **tempshards;
	Shard **tempresults;
	bool do_tree = (algorithm != MOD_FRACTURE_BISECT_FAST &&
					algorithm != MOD_FRACTURE_BISECT_FAST_FILL &&
					algorithm != MOD_FRACTURE_BOOLEAN_FRACTAL);
	FractureModifierData *fmd = (FractureModifierData*) modifiers_findByType(obj, eModifierType_Fracture);

	if (p == NULL || reset)
	{
		if (fm->last_shard_tree)
		{
			BLI_kdtree_free(fm->last_shard_tree);
			fm->last_shard_tree = NULL;
		}

		if (fm->last_shards)
		{
			MEM_freeN(fm->last_shards);
			fm->last_shards = NULL;
		}

		if (!p)
			return;
	}

	if (mode == MOD_FRACTURE_PREFRACTURED && reset && !threaded)
	{
		while (fm->shard_map.first)
		{
			Shard *t = fm->shard_map.first;
			BLI_remlink_safe(&fm->shard_map, t);
			printf("Resetting shard: %d\n", t->shard_id);
			BKE_fracture_shard_free(t, true);
		}
	}

	if (mode == MOD_FRACTURE_PREFRACTURED && !reset)
	{
		//rebuild tree
		if (!fm->last_shard_tree && mode == MOD_FRACTURE_PREFRACTURED)
		{
			Shard *t;
			int ti = 0;
			count = BLI_listbase_count(&fm->shard_map);
			fm->shard_count = count;
			if (do_tree)
			{
				fm->last_shard_tree = BLI_kdtree_new(fm->shard_count);
			}

			fm->last_shards = MEM_callocN(sizeof(Shard*) * fm->shard_count, "last_shards");

			//fill tree from current shardmap
			for (t = fm->shard_map.first; t; t = t->next)
			{
				t->flag &=~ (SHARD_SKIP | SHARD_DELETE);

				if (do_tree)
				{
					BLI_kdtree_insert(fm->last_shard_tree, ti, t->raw_centroid);
				}

				fm->last_shards[ti] = t;
				ti++;
			}

			if (do_tree)
			{
				BLI_kdtree_balance(fm->last_shard_tree);
			}

			p->flag |= SHARD_DELETE;
		}
	}
	else
	{
		fm->last_shard_tree = NULL;
		fm->last_shards = NULL;
	}

	tempshards = MEM_callocN(sizeof(Shard *) * expected_shards, "tempshards");
	tempresults = MEM_callocN(sizeof(Shard *) * expected_shards, "tempresults");

	p->flag = 0;
	p->flag |= SHARD_FRACTURED;

	if (mode == MOD_FRACTURE_DYNAMIC)
	{
		copy_v3_v3(pcentroid, p->centroid);
		parent_id = p->shard_id;
		//remove parent shard from map as well
		BLI_remlink(&fm->shard_map, p);
		fm->shard_count--;
		p->shard_id = -2;
	}

	unit_m4(obmat);
	do_prepare_cells(fm, cells, expected_shards, algorithm, p, &centroid, &dm_parent, &bm_parent, &tempshards, &tempresults, override_count, fmd);

	if (fm->last_shard_tree)
	{
		BLI_kdtree_free(fm->last_shard_tree);
		fm->last_shard_tree = NULL;
	}

	if (fm->last_shards)
	{
		MEM_freeN(fm->last_shards);
		fm->last_shards = NULL;
	}

	if (algorithm != MOD_FRACTURE_BISECT_FAST && algorithm != MOD_FRACTURE_BISECT_FAST_FILL) {
		int totvert = p->totvert;
		MVert *mvert = p->mvert;

		KDTree *preselect_tree = BLI_kdtree_new(totvert);
		for (i = 0; i < totvert; i++) {
			BLI_kdtree_insert(preselect_tree, i, mvert[i].co);
		}

		BLI_kdtree_balance(preselect_tree);

		if (algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL)
		{
			//attempt to have some variance atleast here too
			BLI_thread_srandom(0, fmd->point_seed);
		}

		if ((algorithm == MOD_FRACTURE_BOOLEAN) && !threaded)
		{
			#pragma omp parallel for
			for (i = 0; i < expected_shards; i++)	{
				handle_boolean_bisect(fm, obj, expected_shards, algorithm, parent_id, tempshards, dm_parent,
										bm_parent, obmat, inner_material_index, num_cuts, num_levels, fractal,
										&i, smooth, &tempresults, &dm_p, uv_layer, preselect_tree, thresh,
										fmd->orthogonality_factor);
			}
		}
		else {
			for (i = 0; i < expected_shards; i++)	{
				handle_boolean_bisect(fm, obj, expected_shards, algorithm, parent_id, tempshards, dm_parent,
										bm_parent, obmat, inner_material_index, num_cuts, num_levels, fractal,
										&i, smooth, &tempresults, &dm_p, uv_layer, preselect_tree, thresh,
										fmd->orthogonality_factor);
			}
		}

		BLI_kdtree_free(preselect_tree);
	}
	else {

		if (expected_shards == 1)
		{
			/* do not fracture case */
			tempresults[0] = p;
			p->shard_id = -1;
		}
		else
		{
			handle_fast_bisect(fm, expected_shards, algorithm, &bm_parent, obmat, centroid, inner_material_index, parent_id,
							   tempshards, &tempresults, uv_layer, cells, factor, p);
		}
	}

	if (bm_parent != NULL) {
		BM_mesh_free(bm_parent);
		bm_parent = NULL;
	}

	if (dm_parent != NULL) {
		BKE_mesh_free(dm_parent);
		dm_parent = NULL;
	}

	/*only used with fractal, and is doubly freed in case of 1 shard (doubled) */
	if (dm_p != NULL && expected_shards > 2) {
		BKE_mesh_free(dm_p);
		dm_p = NULL;
	}

	if (p)
	{
		BLI_remlink_safe(&fm->shard_map, p);
		BKE_fracture_shard_free(p, true);
		p = NULL;
	}

	fm->shard_count = 0; /* may be not matching with expected shards, so reset... did increment this for
						  *progressbar only */

	//keep empty ids... need to catch this later
	if (mode == MOD_FRACTURE_DYNAMIC)
	{
		j = 1;

		if (fm->shard_map.last)
		{
			j += ((Shard*)(fm->shard_map.last))->shard_id;
		}
	}
	else
	{
		j = 1;
	}

	for (i = 0; i < expected_shards; i++) {
		Shard *s = tempresults[i];
		Shard *t = tempshards[i];

		if (s != NULL) {
			fracture_shard_add(fm, s, mat);
			s->shard_id += j;
			s->parent_id = parent_id;
			//printf("ADDED: %d %d %d\n", i, j, s->shard_id);
			if (parent_id > -1)
			{
				int si = 0;
				MVert *v;

				sub_v3_v3(s->centroid, pcentroid);
				for (si = 0, v = s->mvert; si < s->totvert; si++, v++)
				{
					sub_v3_v3(v->co, pcentroid);
				}
			}
		}

		if (t != NULL && t != s) {
			BKE_fracture_shard_free(t, false);
		}
	}

	if (fm->shard_count == 0)
	{
		//might happen if all has been skipped, but this distracts the halving method (thinks shardmap is empty)
		//so better correct this here
		fm->shard_count = BLI_listbase_count(&fm->shard_map);
	}

	MEM_freeN(tempshards);
	MEM_freeN(tempresults);
}

static Shard *parse_cell(cell c)
{
	Shard *s;
	MVert *mvert = NULL;
	MPoly *mpoly = NULL;
	MLoop *mloop = NULL;
	int *neighbors = NULL;
	int totpoly = 0, totloop = 0, totvert = 0;
	float centr[3];

	totvert = c.totvert;
	if (totvert > 0) {
		mvert = MEM_callocN(sizeof(MVert) * totvert, __func__);
		parse_cell_verts(c, mvert, totvert);
	}

	totpoly = c.totpoly;
	if (totpoly > 0) {
		mpoly = MEM_callocN(sizeof(MPoly) * totpoly, __func__);
		parse_cell_polys(c, mpoly, totpoly, &totloop);
	}
	else
		totloop = 0;

	if (totloop > 0) {
		mloop = MEM_callocN(sizeof(MLoop) * totloop, __func__);
		parse_cell_loops(c, mloop, totloop, mpoly, totpoly);
	}

	if (totpoly > 0) {
		neighbors = MEM_callocN(sizeof(int) * totpoly, __func__);
		parse_cell_neighbors(c, neighbors, totpoly);
	}

	copy_v3_v3(centr, c.centroid);

	s = BKE_fracture_shard_create(mvert, mpoly, mloop, NULL, totvert, totpoly, totloop, 0, false);

	//s->flag &= ~(SHARD_SKIP | SHARD_DELETE);
	s->neighbor_ids = neighbors;
	s->neighbor_count = totpoly;
	copy_v3_v3(s->centroid, centr);
	copy_v3_v3(s->raw_centroid, centr);
	s->raw_volume = c.volume;

	return s;
}

static void parse_cell_verts(cell c, MVert *mvert, int totvert)
{
	int i;

	for (i = 0; i < totvert; i++) {
		float *co = mvert[i].co;
		copy_v3_v3(co, c.verts[i]);
	}
}

static void parse_cell_polys(cell c, MPoly *mpoly, int totpoly, int *r_totloop)
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

	*r_totloop = totloop;
}

static void parse_cell_loops(cell c, MLoop *mloop, int UNUSED(totloop), MPoly *mpoly, int totpoly)
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

static void fracture_shard_vgroup_add(Shard *s, Object *ob, const char* name)
{
	int index = 0, i = 0;
	MDeformVert *dvert;
	if (!defgroup_find_name(ob, name)) {
		BKE_defgroup_new(ob, name);
	}
	index = defgroup_name_index(ob, name);
	dvert = CustomData_get_layer(&s->vertData, CD_MDEFORMVERT);
	if (dvert == NULL) {
		dvert = CustomData_add_layer(&s->vertData, CD_MDEFORMVERT, CD_CALLOC, NULL, s->totvert);
	}
	for (i = 0; i < s->totvert; i++) {
		MDeformVert* dv = dvert + i;
		defvert_add_index_notest(dv, index, 1.0f);
	}
}

static void fracture_shard_material_add(Shard* s, Object *ob, short mat_ofs) {

	/* only use material offsets if we have 3 or more materials; since FM material handling still is a bit odd  */
	/* todo, use index for inner material too, then this hack here isnt necessary any more */

	const short mat_nr_max = ob->totcol > 1 ? ob->totcol - 1 : 0;
	mat_ofs = mat_nr_max ? mat_ofs : 0;

	if (mat_ofs) {
		MPoly *mp;
		int i = 0;

		for (i = 0; i < s->totpoly; i++)
		{
			mp = s->mpoly + i;
			mp->mat_nr = mat_ofs;

			//hrm first material and second (inner) should be untouched...
			CLAMP(mp->mat_nr, 0, mat_nr_max);
		}
	}
}

static void do_intersect(FractureModifierData *fmd, Object* ob, Shard *t, short inner_mat_index,
						 bool is_zero, float mat[4][4], int **shard_counts, int* count,
						 int k, Mesh **dm_parent, bool keep_other_shard, float thresh)
{
	/*just keep appending items at the end here */
	MPoly *mpoly, *mp;
	int totpoly;
	Shard *parent = NULL;
	Shard *s = NULL, *s2 = NULL;
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

static void reset_shards(FractureModifierData *fmd)
{
	if (fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED && fmd->shared->reset_shards)
	{
		FracMesh *fm = fmd->shared->frac_mesh;
		while (fm && fm->shard_map.first)
		{
			Shard *t = fm->shard_map.first;
			BLI_remlink_safe(&fm->shard_map, t);
			printf("Resetting shard (Greasepencil/Cutter Plane): %d\n", t->shard_id);
			BKE_fracture_shard_free(t, true);
		}
		fm->shard_count = 0;
		/* do not reset again afterwards, in case we have multiple point sources */
		if (!fmd->execute_threaded) {
			fmd->shared->reset_shards = false;
		}
	}
}

void BKE_fracture_shard_by_greasepencil(FractureModifierData *fmd, Object *obj, short inner_material_index, float mat[4][4])
{
	bGPDlayer *gpl;
	bGPDframe *gpf;
	bGPDstroke *gps;

	reset_shards(fmd);

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
					intersect_shards_by_dm(fmd, dm, obj, NULL, inner_material_index, mat, true, fmd->boolean_double_threshold);

					BKE_mesh_free(dm);
					dm = NULL;
				}
			}
		}
	}
}

void BKE_fracture_shard_by_planes(FractureModifierData *fmd, Object *obj, short inner_material_index, float mat[4][4])
{
	if (fmd->cutter_group != NULL && obj->type == OB_MESH)
	{
		CollectionObject* go;
		float imat[4][4];

		reset_shards(fmd);
		invert_m4_m4(imat, obj->obmat);

		for (go = fmd->cutter_group->gobject.first; go; go = go->next)
		{
			Object* ob = go->ob;

			printf("Cutting with %s ...\n", ob->id.name);
			/*simple case....one cutter object per object*/
			if (ob->type == OB_MESH) {

				FractureModifierData *fmd2 = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);
				if (fmd2 && BLI_listbase_count(&fmd2->shared->meshIslands) > 0)
				{
					MeshIsland* mi = NULL;
					int j = 0;
					for (mi = fmd2->shared->meshIslands.first; mi; mi = mi->next)
					{
						Mesh *dm;
						MVert *mv, *v = NULL;

						dm = BKE_fracture_mesh_copy(mi->physics_mesh, ob);
						mv = dm->mvert;
						int totvert = dm->totvert;
						int i = 0;

						//printf("Cutting with %s, island %d...\n", ob->id.name, j);
						for (i = 0, v = mv; i < totvert; i++, v++)
						{
							add_v3_v3(v->co, mi->centroid);
						}

						intersect_shards_by_dm(fmd, dm, obj, ob, inner_material_index, mat, false, fmd->boolean_double_threshold);

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

typedef struct FractureData {
	cell *voro_cells;
	FracMesh *fmesh;
	ShardID id;
	int totpoints;
	int algorithm;
	Object *obj;
	Mesh *dm;
	short inner_material_index;
	float mat[4][4];
	int num_cuts;
	float fractal;
	bool smooth;
	int num_levels;
	int mode;
	bool reset;
	int active_setting;
	int num_settings;
	char uv_layer[64];
	float thresh;
	int override_count;
	float factor;
} FractureData;


static void compute_fracture(TaskPool *UNUSED(pool), void *taskdata, int UNUSED(threadid))
{
	FractureData *fd = (FractureData*)taskdata;

	/*Evaluate result*/
	if (fd->totpoints > 0) {
		parse_cells(fd->voro_cells, fd->totpoints, fd->id, fd->fmesh, fd->algorithm, fd->obj, fd->dm, fd->inner_material_index, fd->mat,
					fd->num_cuts, fd->fractal, fd->smooth, fd->num_levels,fd->mode, fd->reset, fd->uv_layer,
					true, fd->thresh, fd->override_count, fd->factor);
	}
}


static FractureData segment_cells(cell *voro_cells, int startcell, int totcells, FracMesh *fmesh, ShardID id,
					int algorithm, Object *obj, Mesh *dm, short
					inner_material_index, float mat[4][4], int num_cuts, float fractal, bool smooth, int num_levels, int mode,
					bool reset, char uv_layer[64], float thresh, int override_count,
					float factor)
{
	FractureData fd;
	fd.fmesh = fmesh;
	fd.id = id;
	fd.totpoints = totcells;
	fd.algorithm = algorithm;
	fd.obj = obj;
	fd.dm = dm;
	fd.inner_material_index = inner_material_index;
	copy_m4_m4(fd.mat, mat);
	fd.num_cuts = num_cuts;
	fd.fractal = fractal;
	fd.smooth = smooth;
	fd.num_levels = num_levels;
	fd.mode = mode;
	fd.reset = reset;
	strncpy(fd.uv_layer, uv_layer, 64);
	fd.thresh = thresh;
	fd.override_count = override_count;
	fd.factor = factor;

	//cell start pointer, only take fd.totpoints cells out
	fd.voro_cells = voro_cells + startcell;

	return fd;
}

void BKE_fracture_shard_by_points(FractureModifierData *fmd, ShardID id, FracPointCloud *pointcloud, Object *obj, Mesh *dm, short
								  inner_material_index, float mat[4][4], bool reset, int override_count)
{
	int n_size = 8;

	Shard *shard;

	float min[3], max[3];
	float theta = 0.001f; /* TODO, container enlargement, because boundbox exact container and boolean might create artifacts */
	int p, i = 0, num = 0, totcell = 0, remainder_start = 0;

	container *voro_container;
	particle_order *voro_particle_order;
	cell *voro_cells;

	TaskScheduler *scheduler = BLI_task_scheduler_get();
	TaskPool *pool;
	FractureData *fdata;

#ifdef USE_DEBUG_TIMER
	double time_start;
#endif

	shard = BKE_shard_by_id(fmd->shared->frac_mesh, id, dm);
	if (!shard || (shard->flag & SHARD_FRACTURED && (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC))) {
		int val = fmd->shards_to_islands ? -1 : 0;
		if (id == val)
		{
			//fallback to entire mesh
			shard = BKE_shard_by_id(fmd->shared->frac_mesh, -1 , dm);
		}
		else
		{
			return;
		}
	}

	printf("Fracturing with %d points...\n", pointcloud->totpoints);
	/* calculate bounding box with theta margin */
	copy_v3_v3(min, shard->min);
	copy_v3_v3(max, shard->max);

	if (shard->shard_id == -2) {
		BKE_fracture_shard_free(shard, true);
	}

	add_v3_fl(min, -theta);
	add_v3_fl(max, theta);

	mul_m4_v3(mat, min);
	mul_m4_v3(mat, max);


	if (fmd->point_source & MOD_FRACTURE_GRID)
	{
		float off[3] =  {0, 0, 0};
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

	/*Disable for fast bisect/fill, dynamic and mousebased for now -> errors and crashes */
	if (fmd->fracture_mode != MOD_FRACTURE_DYNAMIC && reset == true &&
		fmd->frac_algorithm != MOD_FRACTURE_BISECT_FAST && fmd->frac_algorithm != MOD_FRACTURE_BISECT_FAST_FILL && fmd->execute_threaded == true) {
		/*segment cells, give each thread a chunk to work on */
		pool = BLI_task_pool_create(scheduler, NULL);
		num = BLI_task_scheduler_num_threads(scheduler);
		fdata = MEM_callocN(sizeof(FractureData) * (num + 1), "threaded fracture data");
		totcell = pointcloud->totpoints / num;
		for (i = 0; i < num; i++) {
			//give each task a segment of the shards...
			int startcell = i * totcell;
			FractureData fd = segment_cells(voro_cells, startcell, totcell, fmd->shared->frac_mesh, id, fmd->frac_algorithm,
			                                obj, dm, inner_material_index,
											mat, fmd->fractal_cuts, fmd->fractal_amount, fmd->use_smooth, fmd->fractal_iterations,
											fmd->fracture_mode, reset, fmd->uvlayer_name,
											fmd->boolean_double_threshold, override_count, fmd->orthogonality_factor);
			fdata[i] = fd;
		}

		//remainder...
		remainder_start = fdata[0].totpoints * num;
		if (remainder_start < pointcloud->totpoints) {
			int remainder = pointcloud->totpoints - remainder_start;
			int startcell = remainder_start;
			printf("REMAINDER %d %d\n", startcell, remainder);
			fdata[num] = segment_cells(voro_cells, startcell, totcell, fmd->shared->frac_mesh, id, fmd->frac_algorithm, obj, dm, inner_material_index,
									   mat, fmd->fractal_cuts, fmd->fractal_amount, fmd->use_smooth, fmd->fractal_iterations,
									   fmd->fracture_mode, reset, fmd->uvlayer_name,
									   fmd->boolean_double_threshold, override_count, fmd->orthogonality_factor);
		}

		for (i = 0; i < num+1; i++) {
			BLI_task_pool_push(pool, compute_fracture, &fdata[i], false, TASK_PRIORITY_HIGH);
		}

		BLI_task_pool_work_and_wait(pool);
		BLI_task_pool_free(pool);
		MEM_freeN(fdata);
	}
	else {
		/*Evaluate result*/
		parse_cells(voro_cells, pointcloud->totpoints, id, fmd->shared->frac_mesh, fmd->frac_algorithm, obj,
		            dm, inner_material_index, mat,
					fmd->fractal_cuts, fmd->fractal_amount, fmd->use_smooth, fmd->fractal_iterations, fmd->fracture_mode, reset,
					fmd->uvlayer_name, false, fmd->boolean_double_threshold, override_count, fmd->orthogonality_factor);
	}

	/*Free structs in C++ area of memory */
	cells_free(voro_cells, pointcloud->totpoints);
	particle_order_free(voro_particle_order);
	container_free(voro_container);

#ifdef USE_DEBUG_TIMER
	printf("Fracture done, %g\n", PIL_check_seconds_timer() - time_start);
#endif

}

void BKE_fracture_fracmesh_free(FracMesh *fm, bool doCustomData)
{
	if (fm == NULL) {
		return;
	}

	while (fm->shard_map.first) {
		Shard* s = (Shard*)fm->shard_map.first;
		BLI_remlink(&fm->shard_map, s);
		BKE_fracture_shard_free(s, doCustomData);
	}

	if (fm->last_shard_tree)
	{
		BLI_kdtree_free(fm->last_shard_tree);
		fm->last_shard_tree = NULL;
	}

	if (fm->last_shards)
	{
		MEM_freeN(fm->last_shards);
		fm->last_shards = NULL;
	}
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

static Mesh* do_create(FractureModifierData *fmd, int num_verts, int num_loops, int num_polys, int num_edges,
							  bool doCustomData, bool use_packed)
{
	ListBase *shardlist;
	Shard *shard;

	int vertstart, polystart, loopstart, edgestart;

	MVert *mverts;
	MPoly *mpolys;
	MLoop *mloops;
	MEdge *medges;

	Mesh *result = BKE_mesh_new_nomain(num_verts, num_edges, 0, num_loops, num_polys);
	mverts = result->mvert;
	mloops = result->mloop;
	mpolys = result->mpoly;

	if (num_edges > 0) {
		medges = result->medge;
	}

	vertstart = polystart = loopstart = edgestart = 0;
	if (use_packed)
	{
		shardlist = &fmd->shared->pack_storage;
	}
	else if (fmd->shards_to_islands) {
		shardlist = &fmd->shared->islandShards;
	}
	else {
		shardlist = &fmd->shared->frac_mesh->shard_map;
	}

	for (shard = shardlist->first; shard; shard = shard->next)
	{
		MPoly *mp;
		MLoop *ml;
		MEdge *me;
		int i;

		memcpy(mverts + vertstart, shard->mvert, shard->totvert * sizeof(MVert));
		memcpy(mpolys + polystart, shard->mpoly, shard->totpoly * sizeof(MPoly));

		for (i = 0, mp = mpolys + polystart; i < shard->totpoly; ++i, ++mp) {
			/* adjust loopstart index */
			mp->loopstart += loopstart;
		}

		memcpy(mloops + loopstart, shard->mloop, shard->totloop * sizeof(MLoop));

		for (i = 0, ml = mloops + loopstart; i < shard->totloop; ++i, ++ml) {
			/* adjust vertex index */
			ml->v += vertstart;
			ml->e += edgestart;
		}

		if (num_edges > 0) {
			memcpy(medges + edgestart, shard->medge, shard->totedge * sizeof(MEdge));

			for (i = 0, me = medges + edgestart; i < shard->totedge; ++i, ++me) {
				/* adjust vertex indices */
				me->v1 += vertstart;
				me->v2 += vertstart;
			}
		}

		if (doCustomData) {
			BKE_fracture_collect_layers(shard, result, vertstart, polystart, loopstart, edgestart);
		}

		vertstart += shard->totvert;
		polystart += shard->totpoly;
		loopstart += shard->totloop;
		edgestart += shard->totedge;
	}

	return result;
}


static Mesh *fracture_create_mesh(FractureModifierData *fmd, bool doCustomData, bool use_packed)
{
	Shard *s;
	int num_verts, num_polys, num_loops, num_edges;
	Mesh *result;

	num_verts = num_polys = num_loops = num_edges = 0;

	if (use_packed)
	{
		for (s = fmd->shared->pack_storage.first; s; s = s->next) {
			num_verts += s->totvert;
			num_polys += s->totpoly;
			num_loops += s->totloop;
			num_edges += s->totedge;
		}
	}
	else if (fmd->shards_to_islands) {
		for (s = fmd->shared->islandShards.first; s; s = s->next) {
			num_verts += s->totvert;
			num_polys += s->totpoly;
			num_loops += s->totloop;
			num_edges += s->totedge;
		}
	}
	else {

		if (!fmd->shared->frac_mesh)
			return NULL;

		for (s = fmd->shared->frac_mesh->shard_map.first; s; s = s->next) {
			num_verts += s->totvert;
			num_polys += s->totpoly;
			num_loops += s->totloop;
			num_edges += s->totedge;
		}
	}

	result = do_create(fmd, num_verts, num_loops, num_polys, num_edges, doCustomData, use_packed);

	if (num_edges == 0) {
		CustomData_free(&result->edata, 0);
		BKE_mesh_calc_edges(result, true, true);
	}

	if (fmd->fracture_mode != MOD_FRACTURE_EXTERNAL)
	{
		do_marking(fmd, result);
	}

	BKE_mesh_calc_normals(result);
	return result;
}

Mesh* BKE_fracture_assemble_mesh_from_shards(FractureModifierData *fmd, bool doCustomData, bool use_packed)
{
	Mesh *dm_final = NULL;

	if (fmd->shared->dm && !use_packed) {
		BKE_mesh_free(fmd->shared->dm);
		fmd->shared->dm = NULL;
	}

	dm_final = fracture_create_mesh(fmd, doCustomData, use_packed);

	if (!use_packed) {
		fmd->shared->dm = dm_final;
	}

	return dm_final;
}

static void fracture_customdata_layers_copy(CustomData* dest, CustomData *src, int type, int count)
{
	int layer;
	for (layer = 0; layer < src->totlayer; layer++)
	{
		if (src->layers[layer].type == type)
		{
			CustomData_add_layer(dest, type, CD_DUPLICATE, src->layers[layer].data, count);
		}
	}
}

Mesh *BKE_fracture_shard_to_mesh(Shard *s, bool doCustomData)
{
	Mesh *dm;
	MVert *mverts;
	MLoop *mloops;
	MPoly *mpolys;

	dm = BKE_mesh_new_nomain(s->totvert, s->totedge, 0, s->totloop, s->totpoly);

	mverts = dm->mvert;
	mloops = dm->mloop;
	mpolys = dm->mpoly;

	memcpy(mverts, s->mvert, s->totvert * sizeof(MVert));
	memcpy(mloops, s->mloop, s->totloop * sizeof(MLoop));
	memcpy(mpolys, s->mpoly, s->totpoly * sizeof(MPoly));

	if (s->totedge > 0) {
		MEdge *medges = dm->medge;
		memcpy(medges, s->medge, s->totedge * sizeof(MEdge));
	}
	else {
		CustomData_free(&dm->edata, 0);
		BKE_mesh_calc_edges(dm, true, true);
	}

	BKE_mesh_calc_normals(dm);

	if (doCustomData) {
		BKE_fracture_collect_layers(s, dm, 0, 0, 0, 0);
	}

	return dm;
}



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

void BKE_update_velocity_layer(FractureModifierData *fmd, MeshIsland *mi)
{
	Mesh *dm = fmd->shared->visible_mesh_cached;
	float *velX=NULL, *velY=NULL, *velZ = NULL;
	RigidBodyOb *rbo = mi->rigidbody;
	Shard *s, *t = NULL;
	void *pX, *pY, *pZ, *spX = NULL, *spY = NULL, *spZ = NULL;
	float *sX=NULL, *sY=NULL, *sZ=NULL;
	int i = 0;
	ListBase *lb;

	if (!dm)
		return;

	//XXX TODO deal with split shards to islands etc, here take only "real" shards for now
	if (fmd->shards_to_islands) {
		lb = &fmd->shared->islandShards;
	}
	else {
		lb = &fmd->shared->frac_mesh->shard_map;
	}

	for (s = lb->first; s; s = s->next)
	{
		if (s->shard_id == mi->id)
		{
			t = s;
			break;
		}
	}

	pX = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "velX");
	pY = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "velY");
	pZ = CustomData_get_layer_named(&dm->vdata, CD_PROP_FLT, "velZ");

	if (!pX ||!pY || !pZ)
		return;

	velX = (float*)pX;
	velY = (float*)pY;
	velZ = (float*)pZ;

	//XXX how to represent this in mblur ?
	//zero_v3(rbo->ang_vel);

	if (t)
	{
		spX = check_add_layer(NULL, &t->vertData, CD_PROP_FLT, t->totvert, "velX");
		spY = check_add_layer(NULL, &t->vertData, CD_PROP_FLT, t->totvert, "velY");
		spZ = check_add_layer(NULL, &t->vertData, CD_PROP_FLT, t->totvert, "velZ");
	}

	for (i = 0; i < mi->vertex_count; i++)
	{
		if (spX && spY && spZ)
		{
			sX = (float*)spX;
			sY = (float*)spY;
			sZ = (float*)spZ;

			sX[i] = rbo->lin_vel[0] + rbo->ang_vel[0];
			sY[i] = rbo->lin_vel[1] + rbo->ang_vel[1];
			sZ[i] = rbo->lin_vel[2] + rbo->ang_vel[2];
		}

		velX[mi->vertex_indices[i]] = rbo->lin_vel[0] + rbo->ang_vel[0];
		velY[mi->vertex_indices[i]] = rbo->lin_vel[1] + rbo->ang_vel[1];
		velZ[mi->vertex_indices[i]] = rbo->lin_vel[2] + rbo->ang_vel[2];
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

		count = BLI_listbase_count(&fmd->shared->meshIslands);
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
		for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next, i++)
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
		mi = fmd->shared->meshIslands.first;
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
		BKE_mesh_free(dm);
		return;
	}
}



static void cleanup_arrange_shard(FractureModifierData *fmd, Shard *s, float cent[]);
static MeshIsland* find_meshisland(ListBase* meshIslands, int id);
static void do_island_index_map(FractureModifierData *fmd, Object *ob);

FracMesh* BKE_fracture_fracmesh_copy(FracMesh* fm)
{
	FracMesh *fmesh;
	Shard* s, *t;
	int i = 0;

	fmesh = MEM_mallocN(sizeof(FracMesh), __func__);
	fmesh->shard_map.first = NULL;
	fmesh->shard_map.last = NULL;

	for (s = fm->shard_map.first; s; s = s->next)
	{
		t = BKE_fracture_shard_copy(s);
		BLI_addtail(&fmesh->shard_map, t);
		i++;
	}

	fmesh->shard_count = fm->shard_count;
	fmesh->cancel = 0;
	fmesh->running = 0;
	fmesh->progress_counter = 0;
	fmesh->last_shard_tree = NULL;
	fmesh->last_shards = NULL;

	return fmesh;
}


void BKE_fracture_meshislands_free(FractureModifierData* fmd, ListBase* meshIslands, bool do_free_rigidbody,
                                   Scene* scene)
{
	MeshIsland *mi;

	while (meshIslands->first) {
		mi = meshIslands->first;
		BLI_remlink_safe(meshIslands, mi);
		BKE_fracture_mesh_island_free(fmd, mi, do_free_rigidbody, scene);
		mi = NULL;
	}

	meshIslands->first = NULL;
	meshIslands->last = NULL;
}

void BKE_fracture_simulation_free(FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, Scene *scene)
{
	/* what happens with this in dynamic fracture ? worst case, we need a sequence for this too*/
	if (fmd->shards_to_islands) {
		while (fmd->shared->islandShards.first) {
			Shard *s = fmd->shared->islandShards.first;
			BLI_remlink(&fmd->shared->islandShards, s);
			BKE_fracture_shard_free(s, true);
			s = NULL;
		}

		fmd->shared->islandShards.first = NULL;
		fmd->shared->islandShards.last = NULL;
	}

	/* when freeing meshislands, we MUST get rid of constraints before too !!!! */
	BKE_fracture_constraints_free(fmd, scene);

	if (!do_free_seq) {

		BKE_fracture_meshislands_free(fmd, &fmd->shared->meshIslands, do_free_rigidbody, scene);
		fmd->shared->meshIslands.first = NULL;
		fmd->shared->meshIslands.last = NULL;
	}
	else
	{
		if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
		{
			/* in dynamic mode we have to get rid of the entire Meshisland sequence */
			MeshIslandSequence *msq;
			ShardSequence *ssq;

			while (fmd->shared->meshIsland_sequence.first) {
				msq = fmd->shared->meshIsland_sequence.first;
				BLI_remlink(&fmd->shared->meshIsland_sequence, msq);
				BKE_fracture_meshislands_free(fmd, &msq->meshIslands, do_free_rigidbody, scene);
				MEM_freeN(msq);
				msq = NULL;
			}

			fmd->shared->meshIsland_sequence.first = NULL;
			fmd->shared->meshIsland_sequence.last = NULL;

			fmd->shared->meshIslands.first = NULL;
			fmd->shared->meshIslands.last = NULL;

			fmd->shared->current_mi_entry = NULL;

			while (fmd->shared->shard_sequence.first)
			{
				ssq = fmd->shared->shard_sequence.first;
				BLI_remlink(&fmd->shared->shard_sequence, ssq);
				BKE_fracture_fracmesh_free(ssq->frac_mesh, true);
				MEM_freeN(ssq->frac_mesh);
				MEM_freeN(ssq);
			}

			fmd->shared->shard_sequence.first = NULL;
			fmd->shared->shard_sequence.last = NULL;
			fmd->shared->current_shard_entry = NULL;
			fmd->shared->frac_mesh = NULL;
		}
	}

	if (!fmd->explo_shared && fmd->shared->visible_mesh != NULL) {
		BM_mesh_free(fmd->shared->visible_mesh);
		fmd->shared->visible_mesh = NULL;
	}
}

static void free_shards(FractureModifierData *fmd)
{
	Shard *s;

	if (fmd->shared->frac_mesh) {

		if (fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED ||
			fmd->fracture_mode == MOD_FRACTURE_EXTERNAL)
		{
			BKE_fracture_fracmesh_free(fmd->shared->frac_mesh, true);
			MEM_freeN(fmd->shared->frac_mesh);
			fmd->shared->frac_mesh = NULL;
		}
		else
		{
			/* free entire shard sequence here */
			while(fmd->shared->shard_sequence.first)
			{
				ShardSequence* ssq = (ShardSequence*)fmd->shared->shard_sequence.first;
				BLI_remlink(&fmd->shared->shard_sequence, ssq);
				BKE_fracture_fracmesh_free(ssq->frac_mesh, true);
				MEM_freeN(ssq->frac_mesh);
				MEM_freeN(ssq);
			}
			fmd->shared->frac_mesh = NULL;
			fmd->shared->shard_sequence.first = NULL;
			fmd->shared->shard_sequence.last = NULL;

			fmd->shared->current_shard_entry = NULL;
		}
	}

	while (fmd->shared->pack_storage.first) {
		s = fmd->shared->pack_storage.first;
		BLI_remlink(&fmd->shared->pack_storage, s);
		BKE_fracture_shard_free(s, true);
	}
}

void BKE_fracture_modifier_free(FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, Scene *scene)
{
	BKE_fracture_simulation_free(fmd, do_free_seq, (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC) && do_free_rigidbody, scene);

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

	//called on deleting modifier, object or quitting blender...
	//why was this necessary again ?!
	if (fmd->shared->dm) {
		BKE_mesh_free(fmd->shared->dm);
		fmd->shared->dm = NULL;
	}

	if (fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED || fmd->fracture_mode == MOD_FRACTURE_EXTERNAL)
	{
		if (fmd->shared->visible_mesh_cached) {
			BKE_mesh_free(fmd->shared->visible_mesh_cached);
			fmd->shared->visible_mesh_cached = NULL;
		}
	}

	free_shards(fmd);

	if (fmd->shared->vert_index_map != NULL) {
		BLI_ghash_free(fmd->shared->vert_index_map, NULL, NULL);
		fmd->shared->vert_index_map = NULL;
	}

	/*needs to be freed in any case here ?*/
	if (fmd->shared->visible_mesh != NULL) {
		BM_mesh_free(fmd->shared->visible_mesh);
		fmd->shared->visible_mesh = NULL;
	}

	if (fmd->shared->anim_bind)
		MEM_freeN(fmd->shared->anim_bind);

}

void BKE_fracture_free(FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, Scene *scene)
{
	//for prefractured and external case

	if ((!fmd->shared->refresh && !fmd->shared->refresh_constraints)) {
		/* free entire modifier or when job has been cancelled */
		BKE_fracture_modifier_free(fmd, do_free_seq, do_free_rigidbody, scene);

		if (fmd->shared->visible_mesh_cached && !fmd->shards_to_islands)
		{
			/* free visible_mesh_cached in any case ?!*/
			BKE_mesh_free(fmd->shared->visible_mesh_cached);
			fmd->shared->visible_mesh_cached = NULL;
		}
	}
	else if (!fmd->shared->refresh_constraints) {
		/* refreshing all simulation data only, no refracture */
		BKE_fracture_simulation_free(fmd, false, do_free_rigidbody, scene);
	}
	else if (fmd->shared->refresh_constraints && !fmd->is_dynamic_external) {
		/* refresh constraints only */
		BKE_fracture_constraints_free(fmd, scene);
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

static bool in_bbox(float p[3], float min[3], float max[3])
{
	return (p[0] > min[0]) && (p[0] < max[0]) && (p[1] > min[1]) && (p[1] < max[1]) && (p[2] > min[2]) && (p[2] < max[2]);
}

//XXX MOve to BKE_Fracture.h / fracture.c, prefracture stuff should be a function called from op, or in dynamic case from Rigidbody system callback
static void points_from_verts(Object **ob, int totobj, FracPointCloud *points, float mat[4][4], float thresh,
							  FractureModifierData *emd, Mesh *dm, Object *obj, ShardID id)
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
				d = dm;
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

					if (id > 0)
					{
						Shard *sh;
						float min[3], max[3], cent[3];
						arrange_shard(emd, id, false, cent);
						sh = BKE_fracture_shard_find(&emd->shared->frac_mesh->shard_map, id);
						if (sh)
						{
							add_v3_v3v3(min, sh->min, cent);
							add_v3_v3v3(max, sh->max, cent);
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
								  float thresh, FractureModifierData *fmd, ShardID id)
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

						if (id > 0)
						{
							Shard *sh;
							float min[3], max[3], cent[3];
							arrange_shard(fmd, id, false, cent);
							sh = BKE_fracture_shard_find(&fmd->shared->frac_mesh->shard_map, id);
							if (sh)
							{
								add_v3_v3v3(min, sh->min, cent);
								add_v3_v3v3(max, sh->max, cent);
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

FracPointCloud BKE_fracture_points_get(Depsgraph *depsgraph, FractureModifierData *emd, Object *ob, Mesh *fracmesh, ShardID id)
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
		points_from_particles(go, totgroup, scene, &points, ob->obmat, thresh, emd, id);
	}

	if (emd->point_source & (MOD_FRACTURE_OWN_VERTS | MOD_FRACTURE_EXTRA_VERTS)) {
		points_from_verts(go, totgroup, &points, ob->obmat, thresh, emd, fracmesh, ob, id);
	}

	if (emd->point_source & MOD_FRACTURE_GREASEPENCIL && !emd->use_greasepencil_edges) {
		points_from_greasepencil(go, totgroup, &points, ob->obmat, thresh);
	}


	/* local settings, apply per shard!!! Or globally too first. */
	if (emd->point_source & MOD_FRACTURE_UNIFORM)
	{
		float cent[3], bmin[3], bmax[3];
		int count = emd->shard_count;
		Shard *s = NULL;

		INIT_MINMAX(min, max);
		//for limit impact we need entire container always, because we need to determine secondary impacts on the shards at their original pos
		if (!BKE_get_shard_minmax(emd->shared->frac_mesh, id, min, max, fracmesh))
			return points; //id 0 should be entire mesh

		//arrange shards according to their original centroid (parent centroid sum) position in shard-space (else they are centered at 0, 0, 0)
		arrange_shard(emd, id, false, cent);
		add_v3_v3v3(bmax, max, cent);
		add_v3_v3v3(bmin, min, cent);

		//first impact only, so shard has id 0
		if (emd->fracture_mode == MOD_FRACTURE_DYNAMIC) {
			//shrink pointcloud container around impact point, to a size
			s = BKE_shard_by_id(emd->shared->frac_mesh, id, fracmesh);

			copy_v3_v3(max, bmax);
			copy_v3_v3(min, bmin);

			if (s != NULL && s->impact_size[0] > 0.0f && emd->limit_impact) {
				float size[3], nmin[3], nmax[3], loc[3], tmin[3], tmax[3], rloc[3] = {0,0,0}, quat[4] = {1,0,0,0};
				MeshIslandSequence *msq = emd->shared->current_mi_entry->prev ? emd->shared->current_mi_entry->prev :
				                                                                emd->shared->current_mi_entry;
				MeshIsland *mi = NULL;
				RigidBodyOb *rbo = NULL;

				mat4_to_quat(quat, ob->obmat);
				invert_qt(quat);

				if (msq) {
					mi = find_meshisland(&msq->meshIslands, s->parent_id);
					if (!mi) {
						mi = find_meshisland(&msq->meshIslands, id);
					}

					if (mi) {
						rbo = mi->rigidbody;
						copy_v3_v3(rloc, rbo->pos);
						mul_qt_qtqt(quat, rbo->orn, quat);
					}
				}

				print_v3("Impact Loc\n", s->impact_loc);
				print_v3("Impact Size\n", s->impact_size);

				copy_v3_v3(loc, s->impact_loc);

				sub_v3_v3(loc, rloc);
				mul_qt_v3(quat, loc);
				add_v3_v3(loc, s->centroid);

				copy_v3_v3(tmax, s->max);
				copy_v3_v3(tmin, s->min);

				mul_v3_v3fl(size, s->impact_size, 0.75f);
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

		if (emd->frac_algorithm == MOD_FRACTURE_BISECT_FAST || emd->frac_algorithm == MOD_FRACTURE_BISECT_FAST_FILL ||
			emd->frac_algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL) {
			/* XXX need double amount of shards, because we create 2 islands at each cut... so this matches the input count */
			if ((count > 1) || emd->frac_algorithm == MOD_FRACTURE_BOOLEAN_FRACTAL) {
				count--;
				count *= 2;
			}
		}

		//omg, vary the seed here
		if (emd->shards_to_islands && emd->fracture_mode == MOD_FRACTURE_DYNAMIC) {
			BLI_thread_srandom(0, id);
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

				if (id > 0 && emd->cutter_group == NULL)
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
		if (!BKE_get_shard_minmax(emd->shared->frac_mesh, id, min, max, fracmesh))
			return points; //id 0 should be entire mesh

		//arrange shards according to their original centroid (parent centroid sum) position in shard-space (else they are centered at 0, 0, 0)
		arrange_shard(emd, id, false, cent);
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

					if (id > 0 && emd->cutter_group == NULL)
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

static Material* find_material(const char* name)
{
	ID* mat;

	for (mat = G.main->mat.first; mat; mat = mat->next)
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

	return BKE_material_add(G.main, name);
}

//splinter handling is a case for BKE too
static Shard* do_splinters(FractureModifierData *fmd, FracPointCloud points, float(*mat)[4][4], ShardID id, Mesh *dm)
{
	float imat[4][4];

	/*need to add island / shard centroid...*/
	Shard *s = BKE_shard_by_id(fmd->shared->frac_mesh, id, NULL);

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


		if (s) {
			mvert = s->mvert;
			num_verts = s->totvert;
		}
		else
		{
			mvert = dm->mvert;
			num_verts = dm->totvert;
		}

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

	return s;
}

//so is material handling too, XXX TODO move to BKE
static short do_materials(Main* bmain, FractureModifierData *fmd, Object* obj)
{
	short mat_index = 0;

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
			Material* mat_outer = find_material(matname);
			BKE_object_material_slot_add(bmain, obj);
			assign_material(bmain, obj, mat_outer, obj->totcol, BKE_MAT_ASSIGN_OBDATA);

			MEM_freeN(matname);
			matname = NULL;
			matname = BLI_strdupcat(name, "_Inner");
			mat_inner = find_material(matname);
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
			Material* mat_inner = find_material(matname);
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

static void cleanup_splinters(FractureModifierData *fmd, float mat[4][4], Shard *s, Mesh *dm)
{
	if ((fmd->splinter_axis & MOD_FRACTURE_SPLINTER_X) ||
		(fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Y) ||
		(fmd->splinter_axis & MOD_FRACTURE_SPLINTER_Z))
	{
		int i = 0, num_verts = 0;
		MVert* mvert = NULL, *mv;

		if (s) {
			mvert = s->mvert;
			num_verts = s->totvert;
		}
		else
		{
			mvert = dm->mvert;
			num_verts = dm->totvert;
		}

		for (i = 0, mv = mvert; i < num_verts; i++, mv++)
		{
			mul_m4_v3(mat, mv->co);
		}
	}
}

Shard* BKE_fracture_shard_find(ListBase *shards, ShardID id)
{
	Shard *s = shards->first;
	while (s)
	{
		if (s->shard_id == id)
		{
			return s;
		}
		s = s->next;
	}

	return NULL;
}

static void arrange_shard(FractureModifierData *fmd, ShardID id, bool do_verts, float cent[3])
{
	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		zero_v3(cent);
		bool found = false;
		Shard *sh = BKE_fracture_shard_find(&fmd->shared->frac_mesh->shard_map, id);
		ShardSequence *seq = fmd->shared->current_shard_entry;
		while (seq->prev && sh)
		{
			Shard *par = BKE_fracture_shard_find(&seq->prev->frac_mesh->shard_map, sh->parent_id);
			if (par)
			{
				add_v3_v3(cent, par->centroid);
				found = true;
			}

			seq = seq->prev;
		}

		if (found && do_verts)
		{
			add_v3_v3(sh->centroid, cent);
			add_v3_v3(sh->min, cent);
			add_v3_v3(sh->max, cent);

			int i = 0;
			for (i = 0; i < sh->totvert; i++)
			{
				add_v3_v3(sh->mvert[i].co, cent);
			}
		}
	}
}

static void cleanup_arrange_shard(FractureModifierData *fmd, Shard* sh, float cent[3])
{
	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && sh)
	{
		sub_v3_v3(sh->centroid, cent);
		sub_v3_v3(sh->min, cent);
		sub_v3_v3(sh->max, cent);

		int i = 0;
		for (i = 0; i < sh->totvert; i++)
		{
			sub_v3_v3(sh->mvert[i].co, cent);
		}
	}
}


void BKE_fracture_points(FractureModifierData *fmd, Object* obj, Mesh *dm, ShardID id,
							   int override_count, Depsgraph* depsgraph, Main *bmain)
{
	/* dummy point cloud, random */
	FracPointCloud points;

	points = BKE_fracture_points_get(depsgraph, fmd, obj, dm, id);

	if (points.totpoints > 0 || fmd->use_greasepencil_edges) {
		bool temp = fmd->shards_to_islands;
		short mat_index = 0;
		float mat[4][4];
		Shard *s = NULL;
		float cent[3];

		arrange_shard(fmd, id, true, cent);

		/*splinters... just global axises and a length, for rotation rotate the object */
		s = do_splinters(fmd, points, &mat, id, dm);


		mat_index = do_materials(bmain, fmd, obj);
		mat_index = mat_index > 0 ? mat_index - 1 : mat_index;

		if (points.totpoints > 0) {
			BKE_fracture_shard_by_points(fmd, id, &points, obj, dm, mat_index, mat, fmd->shared->reset_shards, override_count);
		}

		/*TODO, limit this to settings shards !*/
		if (fmd->point_source & MOD_FRACTURE_GREASEPENCIL && fmd->use_greasepencil_edges) {
			BKE_fracture_shard_by_greasepencil(fmd, obj, mat_index, mat);
		}

		cleanup_arrange_shard(fmd, s, cent);

		/* here we REALLY need to fracture so deactivate the shards to islands flag and activate afterwards */
		fmd->shards_to_islands = false;
		BKE_fracture_assemble_mesh_from_shards(fmd, true, false);
		fmd->shards_to_islands = temp;

		cleanup_splinters(fmd, mat, s, dm);

		if (!fmd->auto_execute && fmd->execute_threaded) {
			fmd->shared->reset_shards = false;
		}
	}
	MEM_freeN(points.points);
}

//this is the main fracture function, outsource to BKE, so op or rb system can call it
void BKE_fracture_do(FractureModifierData *fmd, ShardID id, Object *obj, Mesh *dm, Depsgraph *depsgraph, Main* bmain)
{
	short mat_index = 0;
	ShardID* ids = NULL;

	if (fmd->cutter_group != NULL) {
		//attempt to combine fracture by cutter group with regular fracture
		float mat[4][4];
		Shard* s = NULL;
		int count = 0, i = 0;
		bool reset = fmd->shared->reset_shards;

		unit_m4(mat);
		//mat_index = do_materials(fmd, obj);
		//mat_index = mat_index > 0 ? mat_index - 1 : mat_index;

		BKE_fracture_shard_by_planes(fmd, obj, mat_index, mat);

		ids = (ShardID*)MEM_callocN(sizeof(ShardID), "iDs");
		for (s = fmd->shared->frac_mesh->shard_map.first; s; s = s->next)
		{
			printf("Adding Shard ID: %d %d\n", count, s->shard_id);
			if (count > 0) {
				ids = MEM_reallocN_id(ids, sizeof(ShardID) * (count+1), "iDs");
			}

			ids[count] = s->shard_id;
			count++;
		}

		fmd->shared->reset_shards = false;

		for (i = 0; i < count; i++)
		{
			printf("Fracturing Shard ID: %d %d\n", i, ids[i]);
			BKE_fracture_points(fmd, obj, dm, ids[i], 0, depsgraph, bmain);
		}

		fmd->shared->reset_shards = reset;

		MEM_freeN(ids);
	}
	else {
		BKE_fracture_points(fmd, obj, dm, id, -1, depsgraph, bmain);
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

//XXX BKE
static int do_shard_to_island(FractureModifierData *fmd, BMesh* bm_new, ShardID par_id, float centroid[3])
{
	Mesh *dmtemp;
	Shard *s;

	if ((fmd->shards_to_islands || (fmd->shared->frac_mesh && fmd->shared->frac_mesh->shard_count < 2)) && (!fmd->dm_group)) {
		/* store temporary shards for each island */
		int id = 0;

		dmtemp = BKE_fracture_bmesh_to_mesh(bm_new);
		s = BKE_fracture_shard_create(dmtemp->mvert, dmtemp->mpoly,
									  dmtemp->mloop, dmtemp->medge,
									  dmtemp->totvert, dmtemp->totpoly,
									  dmtemp->totloop, dmtemp->totedge,
									  true);
		BKE_fracture_custom_data_mesh_to_shard(s, dmtemp);

		/*for dynamic mode, store this in the main shardmap instead of separately */
		if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC) {

			id = BLI_listbase_count(&fmd->shared->frac_mesh->shard_map);
			s->shard_id = id;
			s->parent_id = par_id;
			s->flag = SHARD_INTACT;
			BLI_addtail(&fmd->shared->frac_mesh->shard_map, s);
			fmd->shared->frac_mesh->shard_count = id + 1;
		}
		else {
			id = BLI_listbase_count(&fmd->shared->islandShards);
			s->shard_id = id;
			s->parent_id = -1;
			BLI_addtail(&fmd->shared->islandShards, s);
		}

		copy_v3_v3(centroid, s->centroid);

		BKE_mesh_free(dmtemp);
		dmtemp = NULL;

		return id;
	}

	return -1;
}

//bke
static void do_rigidbody(Main* bmain, Scene* scene, MeshIsland* mi, Object* ob, Mesh *orig_dm, short rb_type, int i)
{
	mi->rigidbody = NULL;
	mi->rigidbody = BKE_rigidbody_create_shard(bmain, scene, ob, NULL, mi);
	mi->rigidbody->type = rb_type;
	mi->rigidbody->meshisland_index = i;
	BKE_rigidbody_calc_shard_mass(ob, mi, orig_dm);
}

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

static void do_fix_normals(FractureModifierData *fmd, MeshIsland *mi)
{
	/* copy fixed normals to physicsmesh too, for convert to objects */
	if (fmd->fix_normals) {
		MVert *verts, *mv;
		int j = 0, totvert = 0;
		totvert = mi->vertex_count;
		verts = mi->physics_mesh->mvert;
		for (mv = verts, j = 0; j < totvert; mv++, j++) {
			short no[3];
			no[0] = mi->vertno[j * 3];
			no[1] = mi->vertno[j * 3 + 1];
			no[2] = mi->vertno[j * 3 + 2];

			copy_v3_v3_short(mv->no, no);
		}
	}
}

static float do_setup_meshisland(FractureModifierData *fmd, Object *ob, int totvert, float centroid[3],
								 BMVert **verts, float *vertco, short *vertno, BMesh **bm_new, Mesh *orig_dm,
								 int id, Scene *scene, Main *bmain)
{
	MeshIsland *mi;
	Mesh *dm;
	float min[3], max[3], vol = 0;
	int i = 0;
	short rb_type = RBO_TYPE_ACTIVE;
	struct BMeshToMeshParams bmt = {.calc_object_remap = 0};

	mi = MEM_callocN(sizeof(MeshIsland), "meshIsland");
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
	mi->vertices = verts; /*those are temporary only !!! */
	mi->vertco = MEM_mallocN(sizeof(float) * 3 * totvert, "mi->vertco");
	memcpy(mi->vertco, vertco, 3 * totvert * sizeof(float));

	mi->vertno = MEM_mallocN(sizeof(short) * 3 * totvert, "mi->vertco");
	memcpy(mi->vertno, vertno, 3 * totvert * sizeof(short));
	zero_v3(mi->start_co);

	BM_mesh_normals_update(*bm_new);
	dm = BKE_fracture_bmesh_to_mesh(*bm_new);
	BM_mesh_free(*bm_new);
	*bm_new = NULL;

	mi->physics_mesh = dm;
	mi->vertex_count = totvert;

	mi->vertex_indices = MEM_mallocN(sizeof(int) * mi->vertex_count, "mi->vertex_indices");
	for (i = 0; i < mi->vertex_count; i++) {
		mi->vertex_indices[i] = mi->vertices[i]->head.index;
	}

	do_fix_normals(fmd, mi);

	copy_v3_v3(mi->centroid, centroid);
	mi->bb = BKE_boundbox_alloc_unit();
	BKE_boundbox_init_from_minmax(mi->bb, min, max);
	mi->participating_constraints = NULL;
	mi->participating_constraint_count = 0;

	vol = bbox_vol(mi->bb);
	if (vol > fmd->max_vol) {
		fmd->max_vol = vol;
	}

	mi->vertices_cached = NULL;

	rb_type = do_vert_index_map(fmd, mi, NULL);
	i = BLI_listbase_count(&fmd->shared->meshIslands);
	do_rigidbody(bmain, scene, mi, DEG_get_original_object(ob), orig_dm, rb_type, i);

	mi->start_frame = scene->rigidbody_world->shared->pointcache->startframe;
	BLI_addtail(&fmd->shared->meshIslands, mi);

	return vol;
}

static float mesh_separate_tagged(FractureModifierData *fmd, Object *ob, BMVert **v_tag, int v_count,
								  float *startco, BMesh *bm_work, short *startno, Mesh *orig_dm, ShardID par_id,
								  Scene *scene)
{
	BMesh *bm_new;
	BMesh *bm_old = bm_work;
	float centroid[3];
	float vol;
	int id;

	BMVert *v;
	BMIter iter;

	bm_new = BM_mesh_create(&bm_mesh_allocsize_default, &((struct BMeshCreateParams){.use_toolflags = true,}));
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

	BM_calc_center_centroid(bm_new, centroid, false);
	BM_mesh_elem_index_ensure(bm_new, BM_VERT | BM_EDGE | BM_FACE);

	//overwrite centroid with shard centroid here if we have a valid shard
	id = do_shard_to_island(fmd, bm_new, par_id, centroid);

	BM_ITER_MESH (v, &iter, bm_new, BM_VERTS_OF_MESH) {
		/* eliminate centroid in vertex coords */
		sub_v3_v3(v->co, centroid);
	}

	//xxxx need some solution for G.main
	vol = do_setup_meshisland(fmd, ob, v_count, centroid, v_tag, startco, startno, &bm_new, orig_dm, id, scene, G.main);

	/* deselect loose data - this used to get deleted,
	 * we could de-select edges and verts only, but this turns out to be less complicated
	 * since de-selecting all skips selection flushing logic */
	BM_mesh_elem_hflag_disable_all(bm_old, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_TAG, false);

	return vol;
}

static void handle_vert(FractureModifierData *fmd, Mesh *dm, BMVert* vert, BMVert** orig_work,
						float **startco, short **startno, BMVert*** v_tag, int *tot, int *tag_counter)
{
	/* treat the specified vert and put it into the tagged array, also store its coordinates and normals
	 * for usage in meshislands later on */

	short no[3];
	short vno[3];

	if (*v_tag == NULL)
		*v_tag = MEM_callocN(sizeof(BMVert *), "v_tag");

	if (*startco == NULL)
		*startco = MEM_callocN(sizeof(float), "mesh_separate_loose->startco");

	if (*startno == NULL)
		*startno = MEM_callocN(sizeof(short), "mesh_separate_loose->startno");

	BM_elem_flag_enable(vert, BM_ELEM_TAG);
	BM_elem_flag_enable(vert, BM_ELEM_INTERNAL_TAG);
	*v_tag = MEM_reallocN(*v_tag, sizeof(BMVert *) * ((*tag_counter) + 1));
	(*v_tag)[(*tag_counter)] = orig_work[vert->head.index];

	*startco = MEM_reallocN(*startco, ((*tag_counter) + 1) * 3 * sizeof(float));
	(*startco)[3 * (*tag_counter)] = vert->co[0];
	(*startco)[3 * (*tag_counter) + 1] = vert->co[1];
	(*startco)[3 * (*tag_counter) + 2] = vert->co[2];

	*startno = MEM_reallocN(*startno, ((*tag_counter) + 1) * 3 * sizeof(short));

	normal_float_to_short_v3(vno, vert->no);
	normal_float_to_short_v3(no, vert->no);
	if (fmd->fix_normals)
		BKE_fracture_normal_find(dm, fmd->shared->nor_tree, vert->co, vno, no, fmd->nor_range);
	(*startno)[3 * (*tag_counter)] = no[0];
	(*startno)[3 * (*tag_counter) + 1] = no[1];
	(*startno)[3 * (*tag_counter) + 2] = no[2];

	(*tot)++;
	(*tag_counter)++;
}

static void mesh_separate_loose_partition(FractureModifierData *fmd, Object *ob, BMesh *bm_work, BMVert **orig_work,
										  Mesh *dm, ShardID id, Scene* scene)
{
	int i, tag_counter = 0;
	BMEdge *e;
	BMVert *v_seed = NULL, **v_tag = NULL;
	BMWalker walker;
	int tot = 0;
	BMesh *bm_old = bm_work;
	int max_iter = bm_old->totvert;
	BMIter iter;
	float *startco = NULL;
	short *startno = NULL;

	if (max_iter > 0 && fmd->shared->frac_mesh) {
		fmd->shared->frac_mesh->progress_counter++;
	}

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
			handle_vert(fmd, dm, v_seed, orig_work, &startco, &startno, &v_tag, &tot, &tag_counter);
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
				handle_vert(fmd, dm, e->v1, orig_work, &startco, &startno, &v_tag, &tot, &tag_counter);
			}
			if (!BM_elem_flag_test(e->v2, BM_ELEM_TAG) && !BM_elem_flag_test(e->v2, BM_ELEM_INTERNAL_TAG)) {
				handle_vert(fmd, dm, e->v2, orig_work, &startco, &startno, &v_tag, &tot, &tag_counter);
			}
		}
		BMW_end(&walker);

		/* Flush the selection to get edge/face selections matching
		 * the vertex selection */
		BKE_bm_mesh_hflag_flush_vert(bm_old, BM_ELEM_TAG);

		/* Move selection into a separate object */
		mesh_separate_tagged(fmd, ob, v_tag, tag_counter, startco, bm_old, startno, dm, id, scene);

		MEM_freeN(v_tag);
		v_tag = NULL;

		MEM_freeN(startco);
		startco = NULL;

		MEM_freeN(startno);
		startno = NULL;

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
				  Mesh *dm, ShardID id, Scene* scene)
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
		mesh_separate_loose_partition(rmd, ob, bm_old, orig_mod, dm, id, scene);
		separated = true;
	}

	if ((bm_new->totvert <= minsize && bm_new->totvert > 0) || (bm_old->totvert == 0)) {
		mesh_separate_loose_partition(rmd, ob, bm_new, orig_new, dm, id, scene);
		separated = true;
	}

	if ((bm_old->totvert > minsize && bm_new->totvert > 0) || (bm_new->totvert == 0 && !separated)) {
		halve(rmd, ob, minsize, &bm_old, &orig_mod, separated, dm, id, scene);
	}

	if ((bm_new->totvert > minsize && bm_old->totvert > 0) || (bm_old->totvert == 0 && !separated)) {
		halve(rmd, ob, minsize, &bm_new, &orig_new, separated, dm, id, scene);
	}


	MEM_freeN(orig_mod);
	MEM_freeN(orig_new);
	BM_mesh_free(bm_new);
	bm_new = NULL;
}

static void mesh_separate_loose(FractureModifierData *rmd, Object *ob, Mesh *dm, ShardID id, Scene* scene)
{
	int minsize = 500;
	BMesh *bm_work;
	BMVert *vert, **orig_start;
	BMIter iter;

	BM_mesh_elem_hflag_disable_all(rmd->shared->visible_mesh, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_SELECT | BM_ELEM_TAG, false);
	bm_work = BM_mesh_copy(rmd->shared->visible_mesh);

	orig_start = MEM_callocN(sizeof(BMVert *) * rmd->shared->visible_mesh->totvert, "orig_start");
	/* associate new verts with old verts, here indexes should match still */
	BM_ITER_MESH (vert, &iter, rmd->shared->visible_mesh, BM_VERTS_OF_MESH)
	{
		orig_start[vert->head.index] = vert;
	}

	BM_mesh_elem_index_ensure(bm_work, BM_VERT);
	BM_mesh_elem_table_ensure(bm_work, BM_VERT);

	/* free old islandshards first, if any */
	while (rmd->shared->islandShards.first && rmd->fracture_mode != MOD_FRACTURE_DYNAMIC) {
		Shard *s = rmd->shared->islandShards.first;
		BLI_remlink(&rmd->shared->islandShards, s);
		BKE_fracture_shard_free(s, true);
		s = NULL;
	}

	rmd->shared->islandShards.first = NULL;
	rmd->shared->islandShards.last = NULL;

	halve(rmd, ob, minsize, &bm_work, &orig_start, false, dm, id, scene);

	MEM_freeN(orig_start);
	orig_start = NULL;
	BM_mesh_free(bm_work);
	bm_work = NULL;

}





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
		ShardSequence *ssq = NULL;

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

static void do_cache_regular(FractureModifierData* fmd, MeshIsland *mi, int thresh_defgrp_index,
							 int ground_defgrp_index, MVert** verts, MDeformVert** dvert, int *vertstart)
{
	int i;

	for (i = 0; i < mi->vertex_count; i++) {
		mi->vertices_cached[i] = (*verts) + (*vertstart) + i;

		/* sum up vertexweights and divide by vertcount to get islandweight*/
		if (*dvert && ((*dvert) + (*vertstart) + i)->dw && fmd->thresh_defgrp_name[0]) {
			float vweight = defvert_find_weight((*dvert) + (*vertstart) + i, thresh_defgrp_index);
			mi->thresh_weight += vweight;
		}

		if (*dvert && ((*dvert) + (*vertstart) + i)->dw && fmd->ground_defgrp_name[0]) {
			float gweight = defvert_find_weight((*dvert) + (*vertstart) + i, ground_defgrp_index);
			mi->ground_weight += gweight;
		}

		if (mi->vertno != NULL && fmd->fix_normals) {
			short sno[3];
			sno[0] = mi->vertno[i * 3];
			sno[1] = mi->vertno[i * 3 + 1];
			sno[2] = mi->vertno[i * 3 + 2];
			copy_v3_v3_short(mi->vertices_cached[i]->no, sno);
		}
	}

	(*vertstart) += mi->vertex_count;
}

static void do_cache_split_islands(FractureModifierData* fmd, MeshIsland *mi, int thresh_defgrp_index,
								   int ground_defgrp_index, MVert** verts, MDeformVert** dvert)
{
	int i;

	for (i = 0; i < mi->vertex_count; i++) {

		int index = mi->vertex_indices[i];
		if (index >= 0 && index <= fmd->shared->visible_mesh->totvert) {
			mi->vertices_cached[i] = (*verts) + index;
		}
		else {
			mi->vertices_cached[i] = NULL;
		}

		if (*dvert && ((*dvert) + index)->dw && fmd->thresh_defgrp_name[0]) {
			float vweight = defvert_find_weight((*dvert) + index, thresh_defgrp_index);
			mi->thresh_weight += vweight;
		}

		if (*dvert && ((*dvert) + index)->dw && fmd->ground_defgrp_name[0]) {
			float gweight = defvert_find_weight((*dvert) + index, ground_defgrp_index);
			mi->ground_weight += gweight;
		}

		if (mi->vertno != NULL && fmd->fix_normals) {
			short sno[3];
			sno[0] = mi->vertno[i * 3];
			sno[1] = mi->vertno[i * 3 + 1];
			sno[2] = mi->vertno[i * 3 + 2];
			copy_v3_v3_short(mi->vertices_cached[i]->no, sno);
		}
	}
}

static Mesh *createCache(FractureModifierData *fmd, Object *ob, Mesh *origdm)
{
	MeshIsland *mi;
	Mesh *dm;
	MVert *verts;
	MDeformVert *dvert = NULL;
	int vertstart = 0;
	const int thresh_defgrp_index = defgroup_name_index(ob, fmd->thresh_defgrp_name);
	const int ground_defgrp_index = defgroup_name_index(ob, fmd->ground_defgrp_name);
	bool orig_chosen = false;

	/*regular fracture case */
	if (fmd->shared->dm && !fmd->shards_to_islands && (fmd->shared->dm->mvert > 0)) {
		dm = BKE_fracture_mesh_copy(fmd->shared->dm, ob);
	}
	/* split to islands or halving case (fast bisect e.g.) */
	else if (fmd->shared->visible_mesh && (fmd->shared->visible_mesh->totface > 0) &&
	         BLI_listbase_count(&fmd->shared->meshIslands) > 1)
	{
		dm = BKE_fracture_bmesh_to_mesh(fmd->shared->visible_mesh);
	}
	else if (origdm != NULL) {
		dm = BKE_fracture_mesh_copy(origdm, ob);
		orig_chosen = true;
	}
	else {
		return NULL;
	}

	BKE_mesh_tessface_ensure(dm);
	BKE_mesh_calc_normals(dm);
	//DM_update_tessface_data(dm);

	verts = dm->mvert;

	if (dvert == NULL)
		dvert = dm->dvert;

	/* we reach this code when we fracture without "split shards to islands", but NOT when we load such a file...
	 * readfile.c has separate code for dealing with this XXX WHY ? there were problems with the mesh...*/
	for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
		if (mi->vertices_cached) {
			MEM_freeN(mi->vertices_cached);
			mi->vertices_cached = NULL;
		}

		if (fmd->thresh_defgrp_name[0]) {
			mi->thresh_weight = 0;
		}

		mi->vertices_cached = MEM_mallocN(sizeof(MVert *) * mi->vertex_count, "mi->vertices_cached");
		if (fmd->shared->dm != NULL && !fmd->shards_to_islands && !orig_chosen && fmd->shared->visible_mesh == NULL) {
			do_cache_regular(fmd, mi, thresh_defgrp_index, ground_defgrp_index, &verts, &dvert, &vertstart);
		}
		else {  /* halving case... */
			do_cache_split_islands(fmd, mi, thresh_defgrp_index, ground_defgrp_index, &verts, &dvert);
		}

		if (mi->vertex_count > 0) {
			mi->thresh_weight /= mi->vertex_count;
			mi->ground_weight /= mi->vertex_count;
		}

		/*disable for dm_group, cannot paint onto this mesh at all */
		if (mi->rigidbody != NULL && fmd->dm_group == NULL && !fmd->is_dynamic_external) {
			mi->rigidbody->type = mi->ground_weight > 0.01f ? RBO_TYPE_PASSIVE : RBO_TYPE_ACTIVE;
		}

		/* use fallback over inner material*/
		BKE_fracture_fill_vgroup(fmd, dm, dvert, ob, NULL);
	}

	return dm;
}

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

static void do_verts_weights(FractureModifierData *fmd, Shard *s, MeshIsland *mi, int vertstart,
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

#define OUT(name, id, co) printf("%s : %d -> (%.2f, %.2f, %.2f) \n", (name), (id), (co)[0], (co)[1], (co)[2]);
#define OUT4(name,id, co) printf("%s : %d -> (%.2f, %.2f, %.2f, %.2f) \n", (name), (id), (co)[0], (co)[1], (co)[2], (co)[3]);



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

static MeshIsland* find_meshisland(ListBase* meshIslands, int id)
{
	MeshIsland* mi = meshIslands->first;
	while (mi)
	{
		if (mi->id == id)
		{
			return mi;
		}

		mi = mi->next;
	}

	return NULL;
}



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
	BLI_addtail(&fmd->shared->meshIslands, mi);
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

	BKE_fracture_physics_mesh_normals_fix(fmd, s, mi, i, orig_dm);

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
	fmd->explo_shared = true;

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
			dm_final = BKE_fracture_autohide_do(fmd, dm, ob, scene);
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

static void do_post_island_creation(FractureModifierData *fmd, Object *ob, Mesh *dm, Scene* scene, Depsgraph *depsgraph)
{
	double start;

	if (((fmd->shared->visible_mesh != NULL && fmd->shared->refresh && (!fmd->explo_shared)) || (fmd->shared->visible_mesh_cached == NULL))
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
	}

	fmd->shared->refresh_constraints = true;
	fmd->shared->refresh_autohide = true;
	fmd->distortion_cached = false;
}




#if 0
/*XXX should never happen */
static void do_clear(FractureModifierData* fmd)
{
	MeshIsland *mi;
	/* nullify invalid data */
	for (mi = fmd->meshIslands.first; mi; mi = mi->next) {
		mi->vertco = NULL;
		mi->vertex_count = 0;
		mi->vertices = NULL;
		if (mi->vertices_cached)
		{
			MEM_freeN(mi->vertices_cached);
			mi->vertices_cached = NULL;
		}
	}

	if (fmd->visible_mesh_cached) {
		fmd->visible_mesh_cached->needsFree = 1;
		fmd->visible_mesh_cached->release(fmd->visible_mesh_cached);
		fmd->visible_mesh_cached = NULL;
	}
}
#endif

void BKE_fracture_do_halving(FractureModifierData *fmd, Object* ob, Mesh *dm, Mesh *orig_dm, bool is_prehalving, ShardID id, Scene* scene)
{
	double start;

	if (fmd->shared->dm && fmd->shards_to_islands && !is_prehalving) {
		fmd->shared->visible_mesh = BKE_fracture_mesh_to_bmesh(fmd->shared->dm);
	}
	else {
		/* split to meshislands now */
		/* ensures indexes automatically*/
		fmd->shared->visible_mesh = BKE_fracture_mesh_to_bmesh(dm);
	}

	start = PIL_check_seconds_timer();
	//printf("Steps: %d \n", fmd->frac_mesh->progress_counter);
	mesh_separate_loose(fmd, ob, orig_dm, id, scene);
	printf("Splitting to islands done, %g \n"/*  Steps: %d \n"*/, PIL_check_seconds_timer() - start);//, fmd->frac_mesh->progress_counter);
}

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
			fmd->explo_shared = false;
		}
	}

	printf("Islands: %d\n", BLI_listbase_count(&fmd->shared->meshIslands));

	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		/* Grrr, due to stupid design of mine (listbase as value in struct instead of pointer)
		 * we have to synchronize the lists here again */

		/* need to ensure(!) old pointers keep valid, else the whole meshisland concept is broken */
		fmd->shared->current_mi_entry->visible_dm = fmd->shared->visible_mesh_cached;
		fmd->shared->current_mi_entry->meshIslands = fmd->shared->meshIslands;
	}
}

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
		int j = 0;
		for (go = fmd->dm_group->gobject.first; go; go = go->next)
		{
			if (go->ob != ob)
			{
				FractureModifierData *fmdi = (FractureModifierData*)modifiers_findByType(go->ob, eModifierType_Fracture);
				if (fmdi)
				{
					int k = 0;
					for (mi = fmdi->shared->meshIslands.first; mi; mi = mi->next){
						if (mi->vertex_indices != NULL)
						{	/* might not existing yet for older files ! */
							int i = 0;
							for (i = 0; i < mi->vertex_count; i++)
							{
								if (!BLI_ghash_haskey(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(mi->vertex_indices[i] + j)))
								{
									BLI_ghash_insert(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(mi->vertex_indices[i] + j), mi);
								}
							}
						}
						k += mi->vertex_count;
					}

					j += k;
					k = 0;
				}
			}
		}
	}
	else {
		for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next){
			int i = 0;
			if (mi->vertex_indices != NULL)
			{	/* might not existing yet for older files ! */
				for (i = 0; i < mi->vertex_count; i++)
				{
					if (!BLI_ghash_haskey(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(mi->vertex_indices[i])))
					{
						BLI_ghash_insert(fmd->shared->vertex_island_map, SET_INT_IN_POINTER(mi->vertex_indices[i]), mi);
					}
				}
			}
		}
	}
}

Mesh *BKE_fracture_prefractured_do(FractureModifierData *fmd, Object *ob, Mesh *dm, Mesh *orig_dm,
							   char names [][66], int count, Scene* scene, Depsgraph *depsgraph)
{
	bool exploOK = false; /* doFracture */

	if ((fmd->shared->refresh) || (fmd->shared->refresh_constraints))
	{
		/* if we changed the fracture parameters */
		BKE_fracture_free(fmd, true, true, scene);

		/* 2 cases, we can have a visible mesh or a cached visible mesh, the latter primarily when loading blend from file or using halving */
		/* free cached mesh in case of "normal refracture here if we have a visible mesh, does that mean REfracture ?*/
		if (fmd->shared->visible_mesh != NULL && !fmd->shards_to_islands && fmd->shared->frac_mesh &&
			fmd->shared->frac_mesh->shard_count > 0 && fmd->shared->refresh)
		{
			if (fmd->shared->visible_mesh_cached) {
			   BKE_mesh_free(fmd->shared->visible_mesh_cached);
			   fmd->shared->visible_mesh_cached = NULL;
			}
		}

		if (fmd->shared->refresh) {
			do_refresh(fmd, ob, dm, orig_dm, NULL, scene);

			if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && names)
			{
				MeshIsland *mi;
				int i = 0;
				for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
					if (i < count) {
						BLI_snprintf(mi->name, sizeof(names[i]), "%s", names[i]);
					}
					i++;
				}
			}
		}

		do_post_island_creation(fmd, ob, dm, scene, depsgraph);
	}

	if (fmd->shared->refresh_autohide) {
		BKE_fracture_autohide_refresh(fmd, ob);

		if (fmd->autohide_dist > 0 && !fmd->distortion_cached) {
			BKE_fracture_automerge_refresh(fmd);
		}
	}

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
			MeshIsland *mi = fmd->shared->meshIslands.first;
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

	/*XXX better rename this, it checks whether we have a valid fractured mesh */
	exploOK = !fmd->explo_shared || (fmd->explo_shared && fmd->shared->dm && fmd->shared->frac_mesh);

#if 0
	if ((!exploOK) || (fmd->visible_mesh == NULL && fmd->visible_mesh_cached == NULL)) {
		do_clear(fmd);
	}
#endif

	return BKE_fracture_result_mesh(fmd, dm, ob, exploOK, scene);
}

void BKE_fracture_initialize(FractureModifierData *fmd, Object *ob, Mesh *dm, Depsgraph *depsgraph)
{
	/*TODO_1 refresh, move to BKE and just call from operator for prefractured case*/
	if (fmd->shared->refresh)
	{
		printf("ADD NEW 1: %s \n", ob->id.name);
		if (fmd->shared->reset_shards)
		{
			if (fmd->shared->dm != NULL) {
				BKE_mesh_free(fmd->shared->dm);
				fmd->shared->dm = NULL;
			}
		}

		if (fmd->shared->frac_mesh != NULL) {
			BKE_fracture_fracmesh_free(fmd->shared->frac_mesh, true);
			fmd->shared->frac_mesh = NULL;
		}

		/* here we just create the fracmesh, in dynamic case we add the first sequence entry as well */
		if (fmd->shared->frac_mesh == NULL) {
			fmd->shared->frac_mesh = BKE_fracture_fracmesh_create();
		}

		/*normal trees and autohide should work in dynamic too, in theory, but disable for now */
		/* build normaltree from origdm */
		if (fmd->shared->nor_tree != NULL) {
			BLI_kdtree_free(fmd->shared->nor_tree);
			fmd->shared->nor_tree = NULL;
		}

		fmd->shared->nor_tree = build_nor_tree(dm);
		if (fmd->shared->face_pairs != NULL) {
			BLI_ghash_free(fmd->shared->face_pairs, NULL, NULL);
			fmd->shared->face_pairs = NULL;
		}

		fmd->shared->face_pairs = BLI_ghash_int_new("face_pairs");

		/*HERE we must know which shard(s) to fracture... hmm shards... we should "merge" states which happen in the same frame automatically !*/
		BKE_fracture_do(fmd, -1, ob, dm, depsgraph, G.main);
	}
}

Shard* BKE_fracture_shard_copy(Shard *s)
{
	Shard *t = BKE_fracture_shard_create(s->mvert, s->mpoly, s->mloop, s->medge, s->totvert, s->totpoly, s->totloop, s->totedge, true);

	CustomData_reset(&t->vertData);
	CustomData_reset(&t->loopData);
	CustomData_reset(&t->polyData);
	CustomData_reset(&t->edgeData);

	CustomData_copy(&s->vertData, &t->vertData, CD_MASK_MDEFORMVERT, CD_DUPLICATE, s->totvert);
	CustomData_copy(&s->loopData, &t->loopData, CD_MASK_MLOOPUV, CD_DUPLICATE, s->totloop);
	CustomData_copy(&s->edgeData, &t->edgeData, CD_MASK_CREASE | CD_MASK_BWEIGHT | CD_MASK_MEDGE, CD_DUPLICATE, s->totedge);

	t->neighbor_count = t->neighbor_count;
	t->neighbor_ids = MEM_mallocN(sizeof(int) * s->neighbor_count, __func__);
	memcpy(t->neighbor_ids, s->neighbor_ids, sizeof(int) * s->neighbor_count);
	copy_v3_v3(t->centroid, s->centroid);
	copy_v3_v3(t->raw_centroid, s->raw_centroid);
	t->raw_volume = s->raw_volume;
	t->shard_id = s->shard_id;
	t->parent_id = s->parent_id;
	copy_v3_v3(t->impact_loc, s->impact_loc);
	copy_v3_v3(t->impact_size, s->impact_size);

	return t;
}

Mesh *BKE_fracture_mesh_from_packdata(FractureModifierData *fmd, Mesh *derivedData)
{
	Mesh *dm = NULL;

	/* keep old way of using dynamic external working as well, without interfering with packing */
	if (fmd->shared->pack_storage.first && !fmd->is_dynamic_external)
	{
		dm = BKE_fracture_assemble_mesh_from_shards(fmd, true, true);
	}
	else {
		dm = derivedData;
	}

	return dm;
}

static void fracture_refresh_all_constraints(Scene* scene, Object* ob, FractureModifierData* fmd)
{
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
						BKE_fracture_constraints_free(fmdi, scene);
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
		MeshIsland *mi = fmd->shared->meshIslands.first;
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
	struct BMeshFromMeshParams bmf = {.calc_face_normal = true};
	struct BMeshCreateParams bmc = {.use_toolflags = true};
	BMesh* bm;

	bm = BM_mesh_create(&bm_mesh_allocsize_default, &bmc);
	BM_mesh_bm_from_me(bm, me, &bmf);

	return bm;
}
