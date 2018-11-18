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
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/hair.c
 *  \ingroup bke
 */

#include <limits.h>
#include <string.h>

#include "MEM_guardedalloc.h"

#include "BLI_kdtree.h"
#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_sort.h"
#include "BLI_string_utf8.h"
#include "BLI_string_utils.h"

#include "DNA_mesh_types.h"
#include "DNA_hair_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_animsys.h"
#include "BKE_customdata.h"
#include "BKE_global.h"
#include "BKE_hair.h"
#include "BKE_hair_iterators.h"
#include "BKE_hair_runtime.h"
#include "BKE_idcode.h"
#include "BKE_library.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_mesh_sample.h"
#include "BKE_object.h"

#include "DEG_depsgraph_query.h"

#include "BLT_translation.h"

void BKE_hair_init(HairSystem *hsys)
{
	BLI_assert(MEMCMP_STRUCT_OFS_IS_ZERO(hsys, id));

	hsys->bb = BKE_boundbox_alloc_unit();

	BKE_hair_curve_data_init(&hsys->curve_data);
	hsys->draw_settings = BKE_hair_draw_settings_new();
}

void *BKE_hair_add(Main *bmain, const char *name)
{
	HairSystem *hsys = BKE_libblock_alloc(bmain, ID_HA, name, 0);

	BKE_hair_init(hsys);

	return hsys;
}

void BKE_hair_curve_data_init(HairCurveData *data)
{
	CustomData_reset(&data->vdata);
	CustomData_reset(&data->cdata);
	CustomData_reset(&data->fdata);
}

/* Does not free the data pointer itself! */
void BKE_hair_curve_data_free(HairCurveData *data)
{
	MEM_SAFE_FREE(data->curves);
	MEM_SAFE_FREE(data->verts);
	MEM_SAFE_FREE(data->follicles);
	CustomData_free(&data->vdata, data->totverts);
	CustomData_free(&data->cdata, data->totcurves);
	CustomData_free(&data->fdata, data->totfollicles);
}

void BKE_hair_curve_data_copy(HairCurveData *dst_data, const HairCurveData *src_data, int flag)
{
	if (src_data->curves) {
		dst_data->curves = MEM_dupallocN(src_data->curves);
	}
	if (src_data->verts) {
		dst_data->verts = MEM_dupallocN(src_data->verts);
	}
	if (src_data->follicles) {
		dst_data->follicles = MEM_dupallocN(src_data->follicles);
	}

	dst_data->totcurves = src_data->totcurves;
	dst_data->totverts = src_data->totverts;
	dst_data->totfollicles = src_data->totfollicles;

	CustomDataMask mask = CD_MASK_EVERYTHING;
	const eCDAllocType alloc_type = (flag & LIB_ID_COPY_CD_REFERENCE) ? CD_REFERENCE : CD_DUPLICATE;
	CustomData_copy(&src_data->vdata, &dst_data->vdata, mask, alloc_type, dst_data->totverts);
	CustomData_copy(&src_data->cdata, &dst_data->cdata, mask, alloc_type, dst_data->totcurves);
	CustomData_copy(&src_data->fdata, &dst_data->fdata, mask, alloc_type, dst_data->totfollicles);
}

/** Free (or release) any data used by this hair system (does not free the hair system itself). */
void BKE_hair_free(HairSystem *hsys)
{
	BKE_hair_batch_cache_free(hsys);

	MEM_SAFE_FREE(hsys->bb);

	if (hsys->edithair)
	{
		EditHair *edit = hsys->edithair;

		BKE_hair_curve_data_free(&edit->curve_data);

		MEM_freeN(edit);
		hsys->edithair = NULL;
	}

	BKE_hair_curve_data_free(&hsys->curve_data);
	MEM_SAFE_FREE(hsys->draw_settings);
	MEM_SAFE_FREE(hsys->mat);

	BKE_animdata_free(&hsys->id, false);
}

/**
 * Only copy internal data of HairSystem ID from source to already allocated/initialized destination.
 * You probably never want to use that directly, use id_copy or BKE_id_copy_ex for typical needs.
 *
 * WARNING! This function will not handle ID user count!
 *
 * \param flag  Copying options (see BKE_library.h's LIB_ID_COPY_... flags for more).
 */
void BKE_hair_copy_data(Main *UNUSED(bmain), HairSystem *hsys_dst, const HairSystem *hsys_src, const int flag)
{
	hsys_dst->bb = MEM_dupallocN(hsys_src->bb);

	hsys_dst->edithair = NULL;

	BKE_hair_curve_data_copy(&hsys_dst->curve_data, &hsys_src->curve_data, flag);

	hsys_dst->mat = MEM_dupallocN(hsys_src->mat);

	if (hsys_src->draw_settings)
	{
		hsys_dst->draw_settings = BKE_hair_draw_settings_copy(hsys_src->draw_settings);
	}

	BKE_hair_runtime_reset(hsys_dst);
}

HairSystem *BKE_hair_copy(Main *bmain, const HairSystem *hsys)
{
	HairSystem *hsys_copy;
	BKE_id_copy_ex(bmain, &hsys->id, (ID **)&hsys_copy, 0, false);
	return hsys_copy;
}

/* Custom data layer functions; those assume that totXXX are set correctly. */
static void hair_ensure_cdlayers_primary(HairSystem *hsys)
{
	// if (!CustomData_get_layer(&hsys->curve_data.fdata, CD_MVERT))
	// 	CustomData_add_layer(&mesh->vdata, CD_MVERT, CD_CALLOC, NULL, mesh->totvert);
}

static void hair_ensure_cdlayers_origindex(HairSystem *hsys)
{
	// if (!CustomData_get_layer(&mesh->vdata, CD_ORIGINDEX))
	// 	CustomData_add_layer(&mesh->vdata, CD_ORIGINDEX, CD_CALLOC, NULL, mesh->totvert);
}

HairSystem *BKE_hair_new_nomain(int verts_len, int curves_len, int follicles_len)
{
	HairSystem *hsys = BKE_libblock_alloc(
	        NULL, ID_HA,
	        BKE_idcode_to_name(ID_HA),
	        LIB_ID_CREATE_NO_MAIN |
	        LIB_ID_CREATE_NO_USER_REFCOUNT |
	        LIB_ID_CREATE_NO_DEG_TAG);
	BKE_libblock_init_empty(&hsys->id);

	/* don't use CustomData_reset(...); because we dont want to touch customdata */
	copy_vn_i(hsys->curve_data.fdata.typemap, CD_NUMTYPES, -1);
	copy_vn_i(hsys->curve_data.vdata.typemap, CD_NUMTYPES, -1);
	copy_vn_i(hsys->curve_data.cdata.typemap, CD_NUMTYPES, -1);

	hsys->curve_data.totverts = verts_len;
	hsys->curve_data.totcurves = curves_len;
	hsys->curve_data.totfollicles = follicles_len;

	hair_ensure_cdlayers_primary(hsys);
	hair_ensure_cdlayers_origindex(hsys);

	return hsys;
}

static HairSystem *hair_new_nomain_from_template_ex(
        const HairSystem *hsys_src,
        int verts_len, int curves_len, int follicles_len,
        CustomDataMask mask)
{
	HairSystem *hsys_dst = BKE_id_new_nomain(ID_HA, NULL);

	hsys_dst->mat = MEM_dupallocN(hsys_src->mat);

	hsys_dst->curve_data.totverts = verts_len;
	hsys_dst->curve_data.totcurves = curves_len;
	hsys_dst->curve_data.totfollicles = follicles_len;

	CustomData_copy(&hsys_src->curve_data.fdata, &hsys_dst->curve_data.fdata, mask, CD_CALLOC, follicles_len);
	CustomData_copy(&hsys_src->curve_data.vdata, &hsys_dst->curve_data.vdata, mask, CD_CALLOC, verts_len);
	CustomData_copy(&hsys_src->curve_data.cdata, &hsys_dst->curve_data.cdata, mask, CD_CALLOC, curves_len);

	/* The destination hair system should at least have valid primary CD layers,
	 * even in cases where the source hair system does not. */
	hair_ensure_cdlayers_primary(hsys_dst);
	hair_ensure_cdlayers_origindex(hsys_dst);

	return hsys_dst;
}

HairSystem *BKE_hair_new_nomain_from_template(
        const HairSystem *hsys_src,
        int verts_len, int curves_len, int follicles_len)
{
	return hair_new_nomain_from_template_ex(
	        hsys_src,
	        verts_len, curves_len, follicles_len,
	        CD_MASK_EVERYTHING);
}

HairSystem *BKE_hair_copy_for_eval(struct HairSystem *source, bool reference)
{
	int flags = (LIB_ID_CREATE_NO_MAIN |
	             LIB_ID_CREATE_NO_USER_REFCOUNT |
	             LIB_ID_CREATE_NO_DEG_TAG |
	             LIB_ID_COPY_NO_PREVIEW);

	if (reference) {
		flags |= LIB_ID_COPY_CD_REFERENCE;
	}

	HairSystem *result;
	BKE_id_copy_ex(NULL, &source->id, (ID **)&result, flags, false);
	return result;
}

void BKE_hair_make_local(Main *bmain, HairSystem *hsys, const bool lib_local)
{
	BKE_id_make_local_generic(bmain, &hsys->id, true, lib_local);
}

bool BKE_hair_minmax(HairSystem *hsys, float min[3], float max[3])
{
	if (hsys->curve_data.totverts == 0)
	{
		return false;
	}

	HairFiberVertex *vert = hsys->curve_data.verts;
	for (int i = 0; i < hsys->curve_data.totverts; ++i) {
		minmax_v3v3_v3(min, max, vert->co);
	}
	return true;
}

BoundBox *BKE_hair_boundbox_get(Object *ob)
{
	BLI_assert(ob->type == OB_HAIR);
	HairSystem *hsys = ob->data;

	if (ob->bb)
		return ob->bb;

	if (hsys->bb == NULL || (hsys->bb->flag & BOUNDBOX_DIRTY)) {
		BKE_hair_boundbox_calc(hsys);
	}

	return hsys->bb;
}

void BKE_hair_boundbox_calc(HairSystem *hsys)
{
	if (hsys->bb == NULL)
	{
		hsys->bb = MEM_callocN(sizeof(BoundBox), "boundbox");
	}

	float min[3], max[3];
	INIT_MINMAX(min, max);
	if (!BKE_hair_minmax(hsys, min, max)) {
		min[0] = min[1] = min[2] = -1.0f;
		max[0] = max[1] = max[2] = 1.0f;
	}

	BKE_boundbox_init_from_minmax(hsys->bb, min, max);
	hsys->bb->flag &= ~BOUNDBOX_DIRTY;
}

/* Find the mesh used as the scalp surface */
struct Mesh* BKE_hair_get_scalp(
        const Depsgraph *depsgraph,
        const Object *ob,
        const HairSystem *hsys)
{
	Object *scalp_object = BKE_hair_get_scalp_object(ob, hsys);
	if (scalp_object)
	{
		return (struct Mesh *)DEG_get_evaluated_id(depsgraph, scalp_object->data);
	}

	return NULL;
}

/* Find the object used as the scalp surface */
struct Object* BKE_hair_get_scalp_object(
        const Object *ob,
        const HairSystem *hsys)
{
	/* TODO add scalp mode and optional object pointer */
	UNUSED_VARS(hsys);

	if (ob->parent && ob->parent->type == OB_MESH) {
		return ob->parent;
	}

	return NULL;
}

/* Calculate surface area of a scalp mesh */
float BKE_hair_calc_surface_area(const Mesh *scalp)
{
	BLI_assert(scalp != NULL);
	
	int numpolys = scalp->totpoly;
	MPoly *mpolys = scalp->mpoly;
	MLoop *mloops = scalp->mloop;
	MVert *mverts = scalp->mvert;

	float area = 0.0f;
	for (int i = 0; i < numpolys; ++i)
	{
		area += BKE_mesh_calc_poly_area(&mpolys[i], mloops + mpolys[i].loopstart, mverts);
	}
	return area;
}

/* Calculate a density value based on surface area and sample count */
float BKE_hair_calc_density_from_count(float area, int count)
{
	return area > 0.0f ? count / area : 0.0f;
}

/* Calculate maximum sample count based on surface area and density */
int BKE_hair_calc_max_count_from_density(float area, float density)
{
	return (int)(density * area);
}

/* Calculate a density value based on a minimum distance */
float BKE_hair_calc_density_from_min_distance(float min_distance)
{
	// max. circle packing density (sans pi factor): 1 / (2 * sqrt(3))
	static const float max_factor = 0.288675135;
	
	return min_distance > 0.0f ? max_factor / (min_distance * min_distance) : 0.0f;
}

/* Calculate a minimum distance based on density */
float BKE_hair_calc_min_distance_from_density(float density)
{
	// max. circle packing density (sans pi factor): 1 / (2 * sqrt(3))
	static const float max_factor = 0.288675135;
	
	return density > 0.0f ? sqrt(max_factor / density) : 0.0f;
}

/* Distribute hair follicles on a scalp mesh */
void BKE_hair_generate_follicles(
        HairSystem* hsys,
        struct Mesh *scalp,
        unsigned int seed,
        int count)
{
	BKE_hair_generate_follicles_ex(hsys, scalp, seed, count, NULL);
}

/* Distribute hair follicles on a scalp mesh.
 * Optional per-loop weights control follicle density on the scalp.
 */
void BKE_hair_generate_follicles_ex(
        HairSystem* hsys,
        struct Mesh *scalp,
        unsigned int seed,
        int count,
        const float *loop_weights)
{
	// Limit max_count to theoretical limit based on area
	float scalp_area = BKE_hair_calc_surface_area(scalp);
	float density = BKE_hair_calc_density_from_count(scalp_area, count);
	float min_distance = BKE_hair_calc_min_distance_from_density(density);
	
	if (hsys->curve_data.follicles)
	{
		MEM_freeN(hsys->curve_data.follicles);
	}
	hsys->curve_data.follicles = MEM_callocN(sizeof(HairFollicle) * count, "hair follicles");
	
	{
		MeshSampleGenerator *gen = BKE_mesh_sample_gen_surface_poissondisk(seed, min_distance, count, loop_weights);
		
		BKE_mesh_sample_generator_bind(gen, scalp);
		
		static const bool use_threads = false;
		hsys->curve_data.totfollicles = BKE_mesh_sample_generate_batch_ex(
		            gen,
		            &hsys->curve_data.follicles->mesh_sample,
		            sizeof(HairFollicle),
		            count,
		            use_threads);
		
		BKE_mesh_sample_free_generator(gen);
	}
	
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

/* ================================= */

void BKE_hair_fiber_curves_begin(HairSystem *hsys, int totcurves)
{
	if (totcurves != hsys->curve_data.totcurves)
	{
		hsys->curve_data.curves = MEM_reallocN(hsys->curve_data.curves, sizeof(HairFiberCurve) * totcurves);
		hsys->curve_data.totcurves = totcurves;

		BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
	}
}

void BKE_hair_set_fiber_curve(HairSystem *hsys, int index, int numverts,
                              float taper_length, float taper_thickness)
{
	BLI_assert(index <= hsys->curve_data.totcurves);
	
	HairFiberCurve *curve = &hsys->curve_data.curves[index];
	curve->numverts = numverts;
	curve->taper_length = taper_length;
	curve->taper_thickness = taper_thickness;
	
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

/* Calculate vertex start indices on all curves based on length.
 * Returns the total number of vertices.
 */
static int hair_curve_calc_vertstart(HairSystem *hsys)
{
	/* Recalculate vertex count and start offsets in curves */
	int vertstart = 0;
	for (int i = 0; i < hsys->curve_data.totcurves; ++i)
	{
		hsys->curve_data.curves[i].vertstart = vertstart;
		vertstart += hsys->curve_data.curves[i].numverts;
	}
	
	return vertstart;
}

void BKE_hair_fiber_curves_end(HairSystem *hsys)
{
	const int totverts = hair_curve_calc_vertstart(hsys);

	if (totverts != hsys->curve_data.totverts)
	{
		hsys->curve_data.verts = MEM_reallocN(hsys->curve_data.verts, sizeof(HairFiberVertex) * totverts);
		hsys->curve_data.totverts = totverts;

		BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
	}
}

void BKE_hair_set_fiber_vertex(HairSystem *hsys, int index, int flag, const float co[3])
{
	BLI_assert(index <= hsys->curve_data.totverts);
	
	HairFiberVertex *vertex = &hsys->curve_data.verts[index];
	vertex->flag = flag;
	copy_v3_v3(vertex->co, co);
	
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

void BKE_hair_set_fiber_curves(HairSystem *hsys, HairCurveData *curves)
{
	if (hsys->curve_data.curves)
	{
		MEM_freeN(hsys->curve_data.curves);
	}
	hsys->curve_data.curves = MEM_dupallocN(hsys->curve_data.curves);
	hsys->curve_data.totcurves = curves->totcurves;

	if (hsys->curve_data.verts)
	{
		MEM_freeN(hsys->curve_data.verts);
	}
	hsys->curve_data.verts = MEM_dupallocN(hsys->curve_data.verts);
	hsys->curve_data.totverts = curves->totverts;

#ifndef NDEBUG
	const int vertcount = hair_curve_calc_vertstart(hsys);
	BLI_assert(vertcount <= hsys->curve_data.totverts);
#endif

	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

void BKE_hair_clear_fiber_curves(HairSystem *hsys)
{
	if (hsys->curve_data.curves) {
		MEM_freeN(hsys->curve_data.curves);
		hsys->curve_data.curves = NULL;
	}
	hsys->curve_data.totcurves = 0;

	if (hsys->curve_data.verts) {
		MEM_freeN(hsys->curve_data.verts);
		hsys->curve_data.verts = NULL;
	}
	hsys->curve_data.totverts = 0;

	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

/* === Depsgraph evaluation === */

void BKE_hair_eval_geometry(const Depsgraph *depsgraph, HairSystem *hsys)
{
	if (G.debug & G_DEBUG_DEPSGRAPH) {
		printf("%s on %s\n", __func__, hsys->id.name);
	}

	if (hsys->bb == NULL || (hsys->bb->flag & BOUNDBOX_DIRTY)) {
		BKE_hair_boundbox_calc(hsys);
	}
}

void BKE_hair_ensure_follicle_space(const Mesh *scalp, HairCurveData *curve_data)
{
	const int num_follicles = curve_data->totfollicles;
	float (*orco)[3] = NULL;
	float (*normals)[3] = NULL;
	float (*tangents)[3] = NULL;
	CustomData *fdata = &curve_data->fdata;
	
	if (CustomData_has_layer(fdata, CD_ORCO)) {
		orco = CustomData_get_layer(fdata, CD_ORCO);
	}
	else {
		orco = CustomData_add_layer(fdata, CD_ORCO, CD_CALLOC, NULL, num_follicles);
	}

	if (CustomData_has_layer(fdata, CD_NORMAL)) {
		normals = CustomData_get_layer(fdata, CD_NORMAL);
	}
	else {
		normals = CustomData_add_layer(fdata, CD_NORMAL, CD_CALLOC, NULL, num_follicles);
	}
	
	if (CustomData_has_layer(fdata, CD_TANGENT)) {
		tangents = CustomData_get_layer(fdata, CD_TANGENT);
	}
	else {
		tangents = CustomData_add_layer(fdata, CD_TANGENT, CD_CALLOC, NULL, num_follicles);
	}

	float (*co)[3] = orco;
	float (*nor)[3] = normals;
	float (*tang)[3] = tangents;
	HairFollicle *follicle = curve_data->follicles;
	for (int i = 0; i < num_follicles; ++i) {
		BKE_mesh_sample_eval(scalp, &follicle->mesh_sample, co, nor, tang);

		++co;
		++nor;
		++tang;
		++follicle;
	}
}

static void hair_build_data(const Mesh *scalp, HairCurveData *curve_data, CustomDataMask datamask)
{
	if (datamask & (CD_MASK_ORCO | CD_MASK_NORMAL | CD_MASK_TANGENT)) {
		BKE_hair_ensure_follicle_space(scalp, curve_data);
	}
}

void BKE_hair_modifiers_calc(const Depsgraph *depsgraph, Scene *scene, Object *ob)
{
	HairSystem *hsys = ob->data;
	const Mesh *scalp = BKE_hair_get_scalp(depsgraph, ob, hsys);
	EditHair *edit = hsys->edithair;

	CustomDataMask datamask = CD_MASK_ORCO | CD_MASK_NORMAL | CD_MASK_TANGENT;

	if (edit) {
		hair_build_data(scalp, &edit->curve_data, datamask);
	}
	else {
		hair_build_data(scalp, &hsys->curve_data, datamask);
	}
}


/* === Export === */

/* Returns number of vertices in a curve after subdivision */
BLI_INLINE int hair_get_strand_subdiv_length(int orig_length, int subdiv)
{
	return ((orig_length - 1) << subdiv) + 1;
}

/* Returns total number of vertices after subdivision */
BLI_INLINE int hair_get_strand_subdiv_numverts(int numstrands, int numverts, int subdiv)
{
	return ((numverts - numstrands) << subdiv) + numstrands;
}

/* Subdivide a curve */
static int hair_curve_subdivide(const HairFiberCurve* curve, const HairFiberVertex* verts,
                                int subdiv, HairFiberVertex *r_verts)
{
	{
		/* Move vertex positions from the dense array to their initial configuration for subdivision.
		 * Also add offset to ensure the curve starts on the scalp surface.
		 */
		const int step = (1 << subdiv);
		BLI_assert(curve->numverts > 0);
		
		HairFiberVertex *dst = r_verts;
		for (int i = 0; i < curve->numverts; ++i) {
			copy_v3_v3(dst->co, verts[i].co);
			dst += step;
		}
	}
	
	/* Subdivide */
	for (int d = 0; d < subdiv; ++d) {
		const int num_edges = (curve->numverts - 1) << d;
		const int hstep = 1 << (subdiv - d - 1);
		const int step = 1 << (subdiv - d);
		
		/* Calculate edge points */
		{
			int index = 0;
			for (int k = 0; k < num_edges; ++k, index += step) {
				add_v3_v3v3(r_verts[index + hstep].co, r_verts[index].co, r_verts[index + step].co);
				mul_v3_fl(r_verts[index + hstep].co, 0.5f);
			}
		}
		
		/* Move original points */
		{
			int index = step;
			for (int k = 1; k < num_edges; ++k, index += step) {
				add_v3_v3v3(r_verts[index].co, r_verts[index - hstep].co, r_verts[index + hstep].co);
				mul_v3_fl(r_verts[index].co, 0.5f);
			}
		}
	}
	
	const int num_verts = ((curve->numverts - 1) << subdiv) + 1;
	return num_verts;
}

/* Calculate tangent and normal vector changes from one segment to the next */
static void hair_curve_transport_frame(const float co1[3], const float co2[3],
                                       float prev_tang[3], float prev_nor[3],
                                       float r_tang[3], float r_nor[3])
{
	/* segment direction */
	sub_v3_v3v3(r_tang, co2, co1);
	normalize_v3(r_tang);
	
	/* rotate the frame */
	float rot[3][3];
	rotation_between_vecs_to_mat3(rot, prev_tang, r_tang);
	mul_v3_m3v3(r_nor, rot, prev_nor);
	
	copy_v3_v3(prev_tang, r_tang);
	copy_v3_v3(prev_nor, r_nor);
}

/* Calculate tangent and normal vectors for all vertices on a curve */
static void hair_curve_calc_vectors(const HairFiberVertex* verts, int numverts,
                                    float (*r_tangents)[3], float (*r_normals)[3])
{
	BLI_assert(numverts >= 2);
	
	float prev_tang[3] = {0.0f, 0.0f, 1.0f};
	float prev_nor[3] = {1.0f, 0.0f, 0.0f};
	
	hair_curve_transport_frame(
	        verts[0].co, verts[1].co,
	        prev_tang, prev_nor,
	        r_tangents[0], r_normals[0]);
	
	for (int i = 1; i < numverts - 1; ++i)
	{
		hair_curve_transport_frame(
		        verts[i-1].co, verts[i+1].co,
		        prev_tang, prev_nor,
		        r_tangents[i], r_normals[i]);
	}
	
	hair_curve_transport_frame(
	        verts[numverts-2].co, verts[numverts-1].co,
	        prev_tang, prev_nor,
	        r_tangents[numverts-1], r_normals[numverts-1]);
}

/* Create a new export cache.
 * This can be used to construct full fiber data for rendering.
 */

HairExportCache* BKE_hair_export_cache_new(void)
{
	HairExportCache *cache = MEM_callocN(sizeof(HairExportCache), "hair export cache");
	return cache;
}

/* Returns flags for missing data parts */

static int hair_export_cache_get_required_updates(const HairExportCache *cache)
{
	int data = 0;
	if (!cache->fiber_curves)
	{
		data |= HAIR_EXPORT_FIBER_CURVES;
	}
	if (!cache->fiber_verts || !cache->fiber_normals || !cache->fiber_tangents)
	{
		data |= HAIR_EXPORT_FIBER_VERTICES;
	}
	if (!cache->follicles)
	{
		data |= HAIR_EXPORT_FOLLICLE_BINDING;
	}
	if (!cache->follicle_root_position)
	{
		data |= HAIR_EXPORT_FOLLICLE_ROOT_POSITIONS;
	}
	return data;
}

/* Include data dependencies of the given flags */

static int hair_export_cache_get_dependencies(int data)
{
	/* Ordering here is important to account for recursive dependencies */
	
	if (data & HAIR_EXPORT_FIBER_CURVES)
		data |= HAIR_EXPORT_FIBER_VERTICES | HAIR_EXPORT_FOLLICLE_BINDING;
	
	if (data & HAIR_EXPORT_FOLLICLE_BINDING)
		data |= HAIR_EXPORT_FOLLICLE_ROOT_POSITIONS;
	
	return data;
}

/* Update an existing export cache to ensure it contains the requested data.
 * Returns flags for data that has been updated.
 */

int BKE_hair_export_cache_update(HairExportCache *cache, const HairSystem *hsys,
                                 int subdiv, Mesh *scalp, int requested_data)
{
	/* Include dependencies */
	int data = hair_export_cache_get_dependencies(requested_data);
	
	int uncached = hair_export_cache_get_required_updates(cache);
	/* Invalid data should already include all dependencies */
	BLI_assert(uncached == hair_export_cache_get_dependencies(uncached));
	
	/* Only update invalidated parts */
	data &= uncached;
	
	if (data & HAIR_EXPORT_FIBER_CURVES) {
		/* Cache subdivided curves */
		const int totcurves = cache->totcurves = hsys->curve_data.totcurves;
		cache->fiber_curves = MEM_reallocN_id(cache->fiber_curves, sizeof(HairFiberCurve) * totcurves, "hair export curves");
		
		int totverts = 0;
		HairFiberCurve *curve_orig;
		HairIterator iter;
		int i;
		BKE_HAIR_ITER_CURVES_INDEX(curve_orig, &iter, &hsys->curve_data, i) {
			HairFiberCurve *curve = &cache->fiber_curves[i];
			
			memcpy(curve, curve_orig, sizeof(HairFiberCurve));
			curve->numverts = hair_get_strand_subdiv_length(curve_orig->numverts, subdiv);
			curve->vertstart = totverts;
			
			totverts  += curve->numverts;
		}
		cache->totverts = totverts;
	}
	
	if (data & HAIR_EXPORT_FIBER_VERTICES) {
		const int totverts = cache->totverts;
		cache->fiber_verts = MEM_reallocN_id(cache->fiber_verts, sizeof(HairFiberVertex) * totverts, "hair export verts");
		cache->fiber_tangents = MEM_reallocN_id(cache->fiber_tangents, sizeof(float[3]) * totverts, "hair export tangents");
		cache->fiber_normals = MEM_reallocN_id(cache->fiber_normals, sizeof(float[3]) * totverts, "hair export normals");
		
		HairFiberCurve *curve_orig;
		HairIterator iter;
		int i;
		BKE_HAIR_ITER_CURVES_INDEX(curve_orig, &iter, &hsys->curve_data, i) {
			const HairFiberVertex *verts_orig = &hsys->curve_data.verts[curve_orig->vertstart];
			const HairFiberCurve *curve = &cache->fiber_curves[i];
			HairFiberVertex *verts = &cache->fiber_verts[curve->vertstart];
			float (*tangents)[3] = &cache->fiber_tangents[curve->vertstart];
			float (*normals)[3] = &cache->fiber_normals[curve->vertstart];
			
			hair_curve_subdivide(curve_orig, verts_orig, subdiv, verts);
			
			hair_curve_calc_vectors(verts, curve->numverts, tangents, normals);
		}
	}

	if (data & HAIR_EXPORT_FOLLICLE_BINDING) {
		cache->follicles = hsys->curve_data.follicles;
		cache->totfollicles = hsys->curve_data.totfollicles;
	}

	if (data & HAIR_EXPORT_FOLLICLE_ROOT_POSITIONS) {
		cache->follicle_root_position = MEM_reallocN_id(cache->follicle_root_position, sizeof(float[3]) * cache->totfollicles, "fiber root position");
		const HairFollicle *follicle;
		HairIterator iter;
		int i;
		BKE_HAIR_ITER_FOLLICLES_INDEX(follicle, &iter, &hsys->curve_data, i) {
			/* Cache fiber root position */
			BKE_mesh_sample_eval(scalp, &follicle->mesh_sample, cache->follicle_root_position[i], NULL, NULL);
		}
	}
	else {
		MEM_SAFE_FREE(cache->follicle_root_position);
	}
	
	return data;
}

/* Free the given export cache */

void BKE_hair_export_cache_free(HairExportCache *cache)
{
	if (cache->fiber_curves)
	{
		MEM_freeN(cache->fiber_curves);
	}
	if (cache->fiber_verts)
	{
		MEM_freeN(cache->fiber_verts);
	}
	if (cache->fiber_tangents)
	{
		MEM_freeN(cache->fiber_tangents);
	}
	if (cache->fiber_normals)
	{
		MEM_freeN(cache->fiber_normals);
	}
	if (cache->follicle_root_position)
	{
		MEM_freeN(cache->follicle_root_position);
	}
	MEM_freeN(cache);
}

/* Invalidate all data in a hair export cache */

void BKE_hair_export_cache_clear(HairExportCache *cache)
{
	/* Invalidate everything */
	BKE_hair_export_cache_invalidate(cache, HAIR_EXPORT_ALL);
}

/* Invalidate part of the data in a hair export cache.
 *
 * Note some parts may get invalidated automatically based on internal dependencies.
 */

void BKE_hair_export_cache_invalidate(HairExportCache *cache, int invalidate)
{
	/* Include dependencies */
	int data = hair_export_cache_get_dependencies(invalidate);

	if (data & HAIR_EXPORT_FIBER_CURVES)
	{
		if (cache->fiber_curves)
		{
			MEM_freeN(cache->fiber_curves);
			cache->fiber_curves = 0;
		}
	}
	if (data & HAIR_EXPORT_FIBER_VERTICES)
	{
		if (cache->fiber_verts)
		{
			MEM_freeN(cache->fiber_verts);
			cache->fiber_verts = NULL;
		}
		if (cache->fiber_tangents)
		{
			MEM_freeN(cache->fiber_tangents);
			cache->fiber_tangents = NULL;
		}
		if (cache->fiber_normals)
		{
			MEM_freeN(cache->fiber_normals);
			cache->fiber_tangents = NULL;
		}
	}
	if (data & HAIR_EXPORT_FOLLICLE_BINDING)
	{
		cache->follicles = NULL;
	}
	if (data & HAIR_EXPORT_FOLLICLE_ROOT_POSITIONS)
	{
		if (cache->follicle_root_position)
		{
			MEM_freeN(cache->follicle_root_position);
			cache->follicle_root_position = NULL;
		}
	}
}
