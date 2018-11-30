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

/** \file blender/blenkernel/intern/fracture_util.c
 *  \ingroup blenkernel
 *  \brief CSG operations
 */

#include "BKE_mesh.h"
#include "BKE_deform.h"
#include "BKE_editmesh.h"
#include "BKE_fracture.h"
#include "BKE_fracture_util.h"
#include "BKE_material.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_boolean.h"

#include "BLI_alloca.h"
#include "BLI_boxpack_2d.h"
#include "BLI_convexhull_2d.h"
#include "BLI_ghash.h"
#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_sys_types.h"
#include "BLI_kdtree.h"

#include "DNA_fracture_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_material_types.h"
#include "DNA_modifier_types.h"

#include "MEM_guardedalloc.h"

#include "bmesh.h"
#include "bmesh_tools.h"

/*prototypes*/
void uv_bbox(float uv[][2], int num_uv, float minv[2], float maxv[2]);
void uv_translate(float uv[][2], int num_uv, float trans[2]);
void uv_scale(float uv[][2], int num_uv, float scale);
void uv_transform(float uv[][2], int num_uv, float mat[2][2]);
void uv_unwrap_raw_geometry(Mesh *dm, char uv_layer[], bool do_boxpack);
static void do_set_inner_material(Mesh *shard, short inner_material_index, Object *ob);

/* UV Helpers */
void uv_bbox(float uv[][2], int num_uv, float minv[2], float maxv[2])
{
	int v;
	INIT_MINMAX2(minv, maxv);

	for (v = 0; v < num_uv; v++) {
		minmax_v2v2_v2(minv, maxv, uv[v]);
	}
}

void uv_translate(float uv[][2], int num_uv, float trans[2])
{
	int v;
	for (v = 0; v < num_uv; v++) {
		uv[v][0] += trans[0];
		uv[v][1] += trans[1];
	}
}

void uv_scale(float uv[][2], int num_uv, float scale)
{
	int v;
	for (v = 0; v < num_uv; v++) {
		uv[v][0] *= scale;
		uv[v][1] *= scale;
	}
}

void uv_transform(float uv[][2], int num_uv, float mat[2][2])
{
	int v;
	for (v = 0; v < num_uv; v++) {
		mul_m2v2(mat, uv[v]);
	}
}

static void do_clean_uv(Mesh *dm, char uv_layer[64])
{
	MLoopUV* mluv = CustomData_get_layer_named(&dm->ldata, CD_MLOOPUV, uv_layer);
	int i, totpoly = dm->totpoly;
	MPoly *mp, *mpoly = dm->mpoly;

	if (mluv)
	{
		for (i = 0, mp = mpoly; i < totpoly; i++, mp++)
		{
			if (mp->mat_nr != 1)
			{	//clean up (set uv coords to zero) all except inner faces (material based)
				int j;
				for (j = mp->loopstart; j < mp->loopstart + mp->totloop; j++)
				{
					mluv[j].uv[0] = 0.0f;
					mluv[j].uv[1] = 0.0f;
				}
			}
		}
	}
}

static void do_unwrap(MPoly *mp, MVert *mvert, MLoop* mloop, int i, MLoopUV **mluv, BoxPack **boxpack)
{
	MLoop *ml;
	int j = 0;
	float (*verts)[3] = MEM_mallocN(sizeof(float[3]) * mp->totloop, "unwrap_shard_dm verts");
	float nor[3];
	float mat[3][3];
	float (*uv)[2] = MEM_mallocN(sizeof(float[2]) * mp->totloop, "unwrap_shard_dm_uv");
	BoxPack *box;
	float uvbbox[2][2];
	float angle;

	/* uv unwrap cells, so inner faces get a uv map */
	for (j = 0; j < mp->totloop; j++) {
		ml = mloop + mp->loopstart + j;
		copy_v3_v3(verts[j], (mvert + ml->v)->co);
	}

	normal_poly_v3(nor, (const float (*)[3])verts, mp->totloop);
	normalize_v3(nor);
	axis_dominant_v3_to_m3(mat, nor);

	for (j = 0; j < mp->totloop; j++) {
		mul_v2_m3v3(uv[j], mat, verts[j]);
	}

	/* rotate uvs for better packing */
	angle = BLI_convexhull_aabb_fit_points_2d((const float (*)[2])uv, mp->totloop);

	if (angle != 0.0f) {
		float matt[2][2];
		angle_to_mat2(matt, angle);
		uv_transform((float (*)[2])uv, mp->totloop, matt);
	}

	/* prepare box packing... one poly is a box */
	box = (*boxpack) + i;
	uv_bbox((float (*)[2])uv, mp->totloop, uvbbox[0], uvbbox[1]);

	uvbbox[0][0] = -uvbbox[0][0];
	uvbbox[0][1] = -uvbbox[0][1];

	uv_translate((float (*)[2])uv, mp->totloop, uvbbox[0]);

	box->w = uvbbox[1][0] + uvbbox[0][0];
	box->h = uvbbox[1][1] + uvbbox[0][1];
	box->index = i;

	/* copy coords back */
	for (j = 0; j < mp->totloop; j++) {
		copy_v2_v2((*mluv)[j + mp->loopstart].uv, uv[j]);
		(*mluv)[j + mp->loopstart].flag = 0;
	}

	MEM_freeN(uv);
	MEM_freeN(verts);
}

void uv_unwrap_raw_geometry(Mesh *dm, char uv_layer[64], bool do_boxpack)
{
	MVert *mvert;
	MLoop *mloop;
	MPoly *mpoly, *mp;
	int totpoly, i = 0;
	MLoopUV *mluv = MEM_callocN(sizeof(MLoopUV) * dm->totloop, "mluv");
	BoxPack *boxpack = MEM_mallocN(sizeof(BoxPack) * dm->totpoly, "boxpack");
	float scale, tot_width, tot_height;

	/* set inner material on child shard */
	mvert = dm->mvert;
	mpoly = dm->mpoly;
	mloop = dm->mloop;
	totpoly = dm->totpoly;
	for (i = 0, mp = mpoly; i < totpoly; i++, mp++) {
		do_unwrap(mp, mvert, mloop, i, &mluv, &boxpack);
	}

	if (do_boxpack)
	{
		/* do box packing and match uvs according to it */
		BLI_box_pack_2d(boxpack, totpoly, &tot_width, &tot_height);

		if (tot_height > tot_width)
			scale = 1.0f / tot_height;
		else
			scale = 1.0f / tot_width;

		for (i = 0, mp = mpoly; i < totpoly; i++, mp++) {
			float trans[2];
			BoxPack *box;
			int j;

			box = boxpack + i;
			trans[0] = box->x;
			trans[1] = box->y;

			for (j = 0; j < mp->totloop; j++)
			{
				uv_translate((float (*)[2])mluv[j + mp->loopstart].uv, 1, trans);
				uv_scale((float (*)[2])mluv[j + mp->loopstart].uv, 1, scale);
			}
		}
	}

	MEM_freeN(boxpack);

	CustomData_add_layer_named(&dm->ldata, CD_MLOOPUV, CD_ASSIGN, mluv, dm->totloop, uv_layer);
}

static bool check_non_manifold(Mesh* dm)
{
	BMesh *bm;
	BMVert* v;
	BMIter iter;
	BMEdge *e;

	/*check for watertightness*/
	bm = BKE_fracture_mesh_to_bmesh(dm);

	if (bm->totvert == 0) {
		BM_mesh_free(bm);
		printf("Empty mesh...\n");
		return true;
	}

	BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
		if (!BM_vert_is_manifold(v)) {
			BM_mesh_free(bm);
			printf("Mesh not manifold...\n");
			return true;
		}
	}

	BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
		if (BM_edge_is_wire(e) ||
			BM_edge_is_boundary(e) ||
			(BM_edge_is_manifold(e) && !BM_edge_is_contiguous(e)) ||
			BM_edge_face_count(e) > 2)
		{
			/* check we never select perfect edge (in test above) */
			BLI_assert(!(BM_edge_is_manifold(e) && BM_edge_is_contiguous(e)));
			BM_mesh_free(bm);
			printf("Mesh not watertight...\n");
			return true;
		}
	}

	BM_mesh_free(bm);
	return false;
}

static bool compare_size(Mesh *result, Mesh *check)
{
	float min[3], max[3];
	float size[3];
	float v1, v2;
	double eps = 0.000001;

	INIT_MINMAX(min, max);
	BKE_mesh_minmax(result, min, max);
	sub_v3_v3v3(size, max, min);

	v1 = size[0] * size[1] * size[2];

	INIT_MINMAX(min, max);
	BKE_mesh_minmax(check, min, max);
	sub_v3_v3v3(size, max, min);

	v2 = size[0] * size[1] * size[2];

	if (v2 > (v1 + eps))
	{
		printf("Size mismatch !\n");
	}

	return v2 <= (v1 + eps);
}

static Mesh* do_fractal(BooleanContext *ctx)
{
	BMFace* f;
	BMIter iter;
	BMesh *bm;
	int i;
	Mesh *ret = NULL;

	/*create a grid plane */
	bm = BM_mesh_create(&bm_mesh_allocsize_default,  &((struct BMeshCreateParams){.use_toolflags = true,}));
	BMO_op_callf(bm, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			"create_grid x_segments=%i y_segments=%i size=%f matrix=%m4",
			1, 1, ctx->cutter_plane_radius, ctx->cutter_plane_matrix);

	/*subdivide the plane fractally*/
	for (i = 0; i < ctx->num_iterations; i++)
	{
		BMO_op_callf(bm,(BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
					 "subdivide_edges edges=ae "
					 "smooth=%f smooth_falloff=%i use_smooth_even=%b "
					 "fractal=%f along_normal=%f "
					 "cuts=%i "
					 "quad_corner_type=%i "
					 "use_single_edge=%b use_grid_fill=%b "
					 "use_only_quads=%b "
					 "seed=%i",
					 0.0f, SUBD_FALLOFF_ROOT, false,
					 ctx->fractal_amount, 1.0f,
					 ctx->num_cuts,
					 SUBD_CORNER_INNERVERT,
					 false, true,
					 true,
					 0);
	}

	//TODO extrude the plane to give thickness... solidify op
	BMO_op_callf(bm, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			"solidify geom=af thickness=%f offset=%f", ctx->cutter_plane_radius * 2, -1.0f);

	BMO_op_callf(bm, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			"triangulate faces=af quad_method=%i ngon_method=%i", MOD_TRIANGULATE_QUAD_BEAUTY, MOD_TRIANGULATE_NGON_BEAUTY);

	BMO_op_callf(bm, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			"recalc_face_normals faces=af");

	BM_ITER_MESH(f, &iter, bm, BM_FACES_OF_MESH)
	{
		if (ctx->use_smooth_inner)
		{
			BM_elem_flag_enable(f, BM_ELEM_SMOOTH);
		}

		if (ctx->inner_material_index > 0)
		{
			f->mat_nr = ctx->inner_material_index;
		}
	}

	/*convert back*/
	ret = BKE_fracture_bmesh_to_mesh(bm);
	BM_mesh_free(bm);

	return ret;
}

static void do_set_inner_material(Mesh *shard, short inner_material_index, Object *ob)
{
	MPoly *mpoly = shard->mpoly, *mp;
	int totpoly = shard->totpoly, i = 0;
	MDeformVert *dvert;
	int totvert = shard->totvert;
	FractureModifierData *fmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);

	/* set inner material and inner group on raw shard to let boolean transfer it over to the target geometry
	 * in the correct places */
	for (i = 0, mp = mpoly; i < totpoly; i++, mp++) {
		if (inner_material_index > 0) {
			mp->mat_nr = inner_material_index;
		}
		mp->flag |= ME_FACE_SEL;
	}

	if (fmd && fmd->inner_defgrp_name[0]) {
		int defgrp = defgroup_name_index(ob, fmd->inner_defgrp_name);
		dvert = CustomData_add_layer(&shard->vdata, CD_MDEFORMVERT, CD_CALLOC, NULL, totvert);
		for (i = 0; i < totvert; i++) {
			defvert_add_index_notest(dvert + i, defgrp, 1.0f);
		}
	}
}

Mesh* BKE_fracture_mesh_boolean(Mesh* geometry, Mesh* shard, Object* obj, BooleanContext* ctx)
{
	Mesh* result = NULL;

#if 0
	int totvert = MIN2(shard->totvert, geometry->totvert);
	int totedge = MIN2(shard->totedge, geometry->totedge);
	int totloop = MIN2(shard->totloop, geometry->totloop);
	int totpoly = MIN2(shard->totpoly, geometry->totpoly);
#endif

	if (ctx->use_fractal == false)
	{
		uv_unwrap_raw_geometry(shard, ctx->uv_layer, true);
		do_set_inner_material(shard, ctx->inner_material_index, obj);
	}

#if 0
	BKE_fracture_copy_customdata(&geometry->vdata, &shard->vdata, CD_MASK_ISLAND, 0, 0, totvert, totvert);
	BKE_fracture_copy_customdata(&geometry->edata, &shard->edata, CD_MASK_ISLAND, 0, 0, totedge, totedge);
	BKE_fracture_copy_customdata(&geometry->ldata, &shard->ldata, CD_MASK_ISLAND, 0, 0, totloop, totloop);
	BKE_fracture_copy_customdata(&geometry->pdata, &shard->pdata, CD_MASK_ISLAND, 0, 0, totpoly, totpoly);
#endif

	result = BKE_boolean_operation(geometry, obj, shard, obj, ctx->operation, ctx->thresh, NULL);
	/*0 == intersection, 2 == difference*/

#if 0
	if (ctx->use_fractal == false)
	{
		if (!result || check_non_manifold(result) || !compare_size(geometry, result))
		{	/* watertightness check */
			if (result) {
				BKE_fracture_mesh_free(result);
			}

			return NULL;
		}
	}
#endif

	return result;
}

void BKE_fracture_mesh_boolean_fractal(Mesh* geometry, Mesh **outputA, Mesh** outputB, Object* obj, BooleanContext* ctx)
{
	Mesh* cutter = do_fractal(ctx);
	/* dont do fancy boxpacking of uv here, the uv of the fractal surface takes ages to process else */
	uv_unwrap_raw_geometry(cutter, ctx->uv_layer, false);
	do_set_inner_material(cutter, ctx->inner_material_index, obj);

	/* first intersect, then difference with the cutter on inverted positions*/
	ctx->use_fractal = true; //just make sure
	ctx->operation = 0;
	*outputA = BKE_fracture_mesh_boolean(geometry, cutter, obj, ctx);

	ctx->operation = 2;
	*outputB = BKE_fracture_mesh_boolean(geometry, cutter, obj, ctx);

	//BKE_mesh_calc_normals(*outputA);
	//BKE_mesh_calc_normals(*outputB);

	BKE_fracture_mesh_free(cutter);
}

static void do_fill(float plane_no[3], BMOperator bmop, BMesh* bm_parent, BisectContext* ctx)
{
	float normal_fill[3];
	BMOperator bmop_fill;
	//BMOperator bmop_attr;

	normalize_v3_v3(normal_fill, plane_no);
	if (ctx->clear_outer == true && ctx->clear_inner == false) {
		negate_v3(normal_fill);
	}
/*
	if (inner_mat_index == 0) { // dont use inner material here
		BMO_op_initf(
			bm_parent, &bmop_fill, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			"triangle_fill edges=%S normal=%v use_dissolve=%b use_beauty=%b",
			&bmop, "geom_cut.out", normal_fill, true, true);
		BMO_op_exec(bm_parent, &bmop_fill);

		BMO_op_initf(bm_parent, &bmop_attr, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
					 "face_attribute_fill faces=%S use_normals=%b use_data=%b",
					 &bmop_fill, "geom.out", false, true);
		BMO_op_exec(bm_parent, &bmop_attr);
		BMO_op_finish(bm_parent, &bmop_attr);

		BMO_slot_buffer_hflag_enable(bm_parent, bmop_fill.slots_out, "geom.out", BM_FACE, BM_ELEM_TAG | BM_ELEM_SELECT, true);
	}
	else { */
		/* use edgenet fill with inner material */
		BMO_op_initf(
			bm_parent, &bmop_fill, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			"edgenet_fill edges=%S mat_nr=%i use_smooth=%b sides=%i",
			&bmop, "geom_cut.out", ctx->inner_material_index, ctx->use_smooth_inner, 2);
		BMO_op_exec(bm_parent, &bmop_fill);
	//}

	/*gah, edgenetfill now overwrites the set materials internally, so attempt to re-set them */
	BMO_slot_buffer_hflag_enable(bm_parent, bmop_fill.slots_out, "faces.out", BM_FACE, BM_ELEM_TAG, true);
	{
		BMFace* f;
		BMIter iter;
		BM_ITER_MESH(f, &iter, bm_parent, BM_FACES_OF_MESH) {
			if (BM_elem_flag_test(f, BM_ELEM_TAG)) {
				f->mat_nr = ctx->inner_material_index;
				if (!ctx->use_smooth_inner) {
					BM_elem_flag_disable(f, BM_ELEM_SMOOTH);
				}
			}
		}
	}

	BMO_op_finish(bm_parent, &bmop_fill);
}

static void bisect_op(BMesh* bm_geometry, float plane_co[3], float plane_no[3], BisectContext* ctx)
{
	BMOperator bmop;
	float thresh = 0.00001f;

	BM_mesh_elem_hflag_enable_all(bm_geometry, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_TAG, false);
	BMO_op_initf(bm_geometry, &bmop, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
				 "bisect_plane geom=%hvef dist=%f plane_co=%v plane_no=%v use_snap_center=%b clear_inner=%b clear_outer=%b",
				 BM_ELEM_TAG, thresh, plane_co, plane_no, false, ctx->clear_inner, ctx->clear_outer);
	BMO_op_exec(bm_geometry, &bmop);

	BM_mesh_elem_hflag_disable_all(bm_geometry, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_TAG, false);

	if (ctx->use_fill) {
		do_fill(plane_no, bmop, bm_geometry, ctx);
	}

	BMO_slot_buffer_hflag_enable(bm_geometry, bmop.slots_out, "geom_cut.out", BM_VERT | BM_EDGE, BM_ELEM_TAG, true);

	BMO_op_finish(bm_geometry, &bmop);
}

static void do_bisect(BMesh* bm_geometry, BMesh* bm_raw_shard, BisectContext* ctx)
{
	BMIter iter;
	BMFace *f;

	float plane_co[3];
	float plane_no[3];
	float imat[4][4];

	invert_m4_m4(imat, ctx->obmat);

	if (ctx->do_fast_bisect)
	{
		/* with fast bisect, only take one given face into account */
		copy_v3_v3(plane_co, ctx->centroid);
		copy_v3_v3(plane_no, ctx->normal);

		bisect_op(bm_geometry, plane_co, plane_no, ctx);
	}
	else {

		BM_ITER_MESH(f, &iter, bm_raw_shard, BM_FACES_OF_MESH)
		{
			/* with bisect, consider all faces (slower, but nicer looking) */
			copy_v3_v3(plane_co, f->l_first->v->co);
			copy_v3_v3(plane_no, f->no);

			mul_m4_v3(imat, plane_co);
			mul_mat3_m4_v3(imat, plane_no);

			bisect_op(bm_geometry, plane_co, plane_no, ctx);
		}
	}
}

static BMesh *limit_geometry(Mesh* geometry, Shard *shard, KDTree *preselect_tree)
{
	int i = 0, r = 0;
	float max_dist = 0;
	KDTreeNearest* n = NULL;
	BMesh *bm_orig = BKE_fracture_mesh_to_bmesh(geometry);
	BMesh *bm_new = BM_mesh_create(&bm_mesh_allocsize_default, &((struct BMeshCreateParams){.use_toolflags = true,}));
	BMIter iter;
	BMFace *f;
#define MY_TAG (1 << 6)

	BM_mesh_elem_toolflags_ensure(bm_new);  /* needed for 'duplicate' bmo */

	CustomData_copy(&bm_orig->vdata, &bm_new->vdata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_orig->edata, &bm_new->edata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_orig->ldata, &bm_new->ldata, CD_MASK_BMESH, CD_CALLOC, 0);
	CustomData_copy(&bm_orig->pdata, &bm_new->pdata, CD_MASK_BMESH, CD_CALLOC, 0);

	CustomData_bmesh_init_pool(&bm_new->vdata, bm_mesh_allocsize_default.totvert, BM_VERT);
	CustomData_bmesh_init_pool(&bm_new->edata, bm_mesh_allocsize_default.totedge, BM_EDGE);
	CustomData_bmesh_init_pool(&bm_new->ldata, bm_mesh_allocsize_default.totloop, BM_LOOP);
	CustomData_bmesh_init_pool(&bm_new->pdata, bm_mesh_allocsize_default.totface, BM_FACE);

	for (i = 0; i < shard->mesh->totvert; i++)
	{
		float dist = len_squared_v3v3(shard->raw_centroid, shard->mesh->mvert[i].co);
		if (dist > max_dist) {
			max_dist = dist;
		}
	}

	BM_mesh_elem_hflag_disable_all(bm_orig, BM_VERT | BM_EDGE | BM_FACE, MY_TAG, false);

	//do a range search first in case we have many verts as in dense geometry
	r = BLI_kdtree_range_search(preselect_tree, shard->raw_centroid, &n, sqrt(max_dist)*1.5f);

	//if we have sparse geometry, just return all
	if (r < shard->mesh->totvert) {

		int j = 0, s = 0;
		KDTreeNearest *n2 = MEM_callocN(sizeof(KDTreeNearest) * shard->mesh->totvert, "n2 kdtreenearest");
		s = BLI_kdtree_find_nearest_n(preselect_tree, shard->raw_centroid, n2, shard->mesh->totvert);
		for (j = 0; j < s; j++)
		{
			BMVert* v;
			int index = n2[j].index;
			if (index < bm_orig->totvert)
			{
				v = BM_vert_at_index(bm_orig, index);
				BM_elem_flag_enable(v, MY_TAG);
			}
		}
		if (n2) {
			MEM_freeN(n2);
			n2 = NULL;
		}
	}
	else {
		BMVert *v = NULL;
		for (i = 0; i < r; i++) {
			int index = n[i].index;
			if (index < bm_orig->totvert)
			{
				v = BM_vert_at_index(bm_orig, index);
				BM_elem_flag_enable(v, MY_TAG);
			}
		}
	}

	BM_ITER_MESH(f, &iter, bm_orig, BM_FACES_OF_MESH) {
		//select bigger squarish faces and long skinny ones
		float area = BM_face_calc_area(f);
		BMIter eiter;
		BMEdge *e;
		BM_ITER_ELEM(e, &eiter, f, BM_EDGES_OF_FACE) {
			if ((area > 0.75f) || (BM_edge_calc_length_squared(e) > max_dist * 3.0f) ||
				BM_elem_flag_test(e->v1, MY_TAG) || BM_elem_flag_test(e->v2, MY_TAG))
			{
				BM_elem_flag_enable(e->v1, MY_TAG);
				BM_elem_flag_enable(e->v2, MY_TAG);
			}
		}
	}

	/* Flush the selection to get edge/face selections matching
	 * the vertex selection */
	BKE_bm_mesh_hflag_flush_vert(bm_orig, MY_TAG);

	BMO_op_callf(bm_orig, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
				 "duplicate geom=%hvef dest=%p", MY_TAG, bm_new);
#undef MY_TAG
	BM_mesh_elem_index_ensure(bm_new, BM_VERT | BM_EDGE | BM_FACE);

	return bm_new;
}

Mesh* BKE_fracture_mesh_bisect(Mesh* geometry, Shard* raw_shard, BisectContext *ctx)
{
	BMesh *bm_raw_shard;
	BMesh *bm_geometry;
	Mesh* result = NULL;

	if (!ctx->do_fast_bisect) {

		uv_unwrap_raw_geometry(raw_shard->mesh, ctx->uv_layer, true);
		bm_raw_shard = BKE_fracture_mesh_to_bmesh(raw_shard->mesh);

		/* if we detected what shards did change and built a limitation tree, we use it here to chop away unneeded
		 * geometry from the original geometry. So less geometry has to be processed and the bisect should be faster */
		if (ctx->geometry_limitation_tree != NULL) {
			bm_geometry = limit_geometry(geometry, raw_shard, ctx->geometry_limitation_tree);
		}
		else {
			/* without limitation tree, just use the full geometry */
			bm_geometry = BKE_fracture_mesh_to_bmesh(geometry);
		}
	}
	else {
		/* limitation tree makes no sense with fast bisect, and raw shard also not */
		bm_geometry = BKE_fracture_mesh_to_bmesh(geometry);
		bm_raw_shard = NULL;
	}

	if (bm_geometry)
	{
		//BM_mesh_elem_hflag_enable_all(bm_parent, BM_VERT | BM_EDGE | BM_FACE, BM_ELEM_TAG, false);
		do_bisect(bm_geometry, bm_raw_shard, ctx);
		result = BKE_fracture_bmesh_to_mesh(bm_geometry);

		BM_mesh_free(bm_geometry);
	}

	if (bm_raw_shard)
	{
		BM_mesh_free(bm_raw_shard);
	}

	return result;
}

void BKE_fracture_mesh_bisect_fast(Mesh* geometry, Mesh **outputA, Mesh** outputB, BisectContext* ctx)
{
	ctx->clear_inner = false;
	ctx->clear_outer = true;
	*outputA = BKE_fracture_mesh_bisect(geometry, NULL, ctx);


	ctx->clear_inner = true;
	ctx->clear_outer = false;
	*outputB = BKE_fracture_mesh_bisect(geometry, NULL, ctx);
}
