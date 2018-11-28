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
 * The Original Code is Copyright (C) 2017 by Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Blender Foundation, Mike Erwin, Dalai Felinto
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file draw_cache_impl_strands.c
 *  \ingroup draw
 *
 * \brief Strands API for render engines
 */

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"
#include "BLI_math_vector.h"
#include "BLI_ghash.h"

#include "DNA_hair_types.h"
#include "DNA_mesh_types.h"
#include "DNA_scene_types.h"

#include "BKE_customdata.h"
#include "BKE_hair.h"
#include "BKE_hair_iterators.h"
#include "BKE_mesh_sample.h"

#include "DEG_depsgraph.h"

#include "GPU_batch.h"
#include "GPU_extensions.h"
#include "GPU_texture.h"

#include "draw_common.h"
#include "draw_cache_impl.h"  /* own include */
#include "draw_hair_private.h"

#include "DRW_render.h"

/* ---------------------------------------------------------------------- */
/* Hair GPUBatch Cache */

/* GPUBatch cache management. */

typedef struct HairBatchCache {
	ParticleHairCache hair;

	/* Control points when in edit mode. */
	ParticleHairCache edit_hair;

	GPUVertBuf *edit_follicle_pos;
	GPUBatch *edit_follicle_points;
	int edit_follicle_point_len;

	bool is_dirty;
	bool is_editmode;
} HairBatchCache;

static void hair_batch_cache_clear(HairSystem *hsys);

static bool hair_batch_cache_valid(HairSystem *hsys)
{
	HairBatchCache *cache = hsys->runtime.draw_batch_cache;

	if (cache == NULL) {
		return false;
	}

	if (cache->is_dirty) {
		return false;
	}

	if (cache->is_editmode != (hsys->edithair != NULL)) {
		return false;
	}

	return true;
}

static void hair_batch_cache_init(HairSystem *hsys)
{
	HairBatchCache *cache = hsys->runtime.draw_batch_cache;
	
	if (!cache) {
		cache = hsys->runtime.draw_batch_cache = MEM_callocN(sizeof(*cache), __func__);
	}
	else {
		memset(cache, 0, sizeof(*cache));
	}
	
	cache->is_dirty = false;
	cache->is_editmode = (hsys->edithair != NULL);
}

static HairBatchCache *hair_batch_cache_get(HairSystem *hsys)
{
	if (!hair_batch_cache_valid(hsys)) {
		hair_batch_cache_clear(hsys);
		hair_batch_cache_init(hsys);
	}
	return hsys->runtime.draw_batch_cache;
}

void DRW_hair_batch_cache_dirty(HairSystem *hsys, int mode)
{
	HairBatchCache *cache = hsys->runtime.draw_batch_cache;
	if (cache == NULL) {
		return;
	}
	switch (mode) {
		case BKE_HAIR_BATCH_DIRTY_ALL:
			cache->is_dirty = true;
			break;
		case BKE_HAIR_BATCH_DIRTY_SELECT:
			cache->is_dirty = true;
			break;
		default:
			BLI_assert(0);
	}
}

static void hair_batch_cache_clear(HairSystem *hsys)
{
	HairBatchCache *cache = hsys->runtime.draw_batch_cache;
	if (cache) {
		particle_batch_cache_clear_hair(&cache->hair);
		particle_batch_cache_clear_hair(&cache->edit_hair);

		GPU_BATCH_DISCARD_SAFE(cache->edit_follicle_points);
		GPU_VERTBUF_DISCARD_SAFE(cache->edit_follicle_pos);
	}
}

void DRW_hair_batch_cache_free(HairSystem *hsys)
{
	hair_batch_cache_clear(hsys);
	MEM_SAFE_FREE(hsys->runtime.draw_batch_cache);
}

static void hair_batch_cache_ensure_count(
        const HairCurveData *curve_data,
        ParticleHairCache *cache)
{
	cache->strands_len = curve_data->totfollicles;
	cache->elems_len = 0;
	cache->point_len = 0;

	HairIterator iter;
	HairFollicle *follicle;
	HairFiberCurve *curve;
	BKE_HAIR_ITER_FOLLICLE_CURVES(follicle, curve, &iter, curve_data) {
		/* Add one for primitive restart */
		cache->elems_len += curve->numverts + 1;
		cache->point_len += curve->numverts;
	}
}

static bool hair_get_local_space(const float (*orco)[3], const float (*normals)[3], const float (*tangents)[3], int index, float r[4][4])
{
	if (!orco || !normals || !tangents) {
		unit_m4(r);
		return false;
	}

	copy_v3_v3(r[1], tangents[index]);
	copy_v3_v3(r[2], normals[index]);
	cross_v3_v3v3(r[0], r[1], r[2]);
	copy_v3_v3(r[3], orco[index]);
	return true;
}

static void hair_batch_cache_fill_segments_proc_pos(
        const HairCurveData *curve_data,
        GPUVertBufRaw *attr_step)
{
	const HairFiberVertex *verts = BKE_hair_get_verts(curve_data);
	float (*orco)[3] = CustomData_get_layer(&curve_data->fdata, CD_ORCO);
	float (*normals)[3] = CustomData_get_layer(&curve_data->fdata, CD_NORMAL);
	float (*tangents)[3] = CustomData_get_layer(&curve_data->fdata, CD_TANGENT);
	HairIterator iter;
	HairFollicle *follicle;
	HairFiberCurve *curve;
	int curve_idx;
	BKE_HAIR_ITER_FOLLICLE_CURVES_INDEX(follicle, curve, &iter, curve_data, curve_idx) {
		if (curve->numverts < 2) {
			continue;
		}
		float hmat[4][4];
		hair_get_local_space(orco, normals, tangents, curve_idx, hmat);
		const HairFiberVertex *vert_start = &verts[curve->vertstart];
		float total_len = 0.0f;
		const float *co_prev = NULL;
		float *seg_data_first;
		for (int j = 0; j < curve->numverts; j++) {
			float *seg_data = (float *)GPU_vertbuf_raw_step(attr_step);
			mul_v3_m4v3(seg_data, hmat, vert_start[j].co);
			if (co_prev) {
				total_len += len_v3v3(co_prev, seg_data);
			}
			else {
				seg_data_first = seg_data;
			}
			seg_data[3] = total_len;
			co_prev = vert_start[j].co;
		}
		if (total_len > 0.0f) {
			/* Divide by total length to have a [0-1] number. */
			for (int j = 0; j < curve->numverts; j++, seg_data_first += 4) {
				seg_data_first[3] /= total_len;
			}
		}
	}
}

static void hair_batch_cache_ensure_procedural_pos(
        const HairCurveData *curve_data,
        ParticleHairCache *cache)
{
	if (cache->proc_point_buf != NULL) {
		return;
	}

	/* initialize vertex format */
	GPUVertFormat format = {0};
	uint pos_id = GPU_vertformat_attr_add(&format, "posTime", GPU_COMP_F32, 4, GPU_FETCH_FLOAT);

	cache->proc_point_buf = GPU_vertbuf_create_with_format(&format);
	GPU_vertbuf_data_alloc(cache->proc_point_buf, cache->point_len);

	GPUVertBufRaw pos_step;
	GPU_vertbuf_attr_get_raw_data(cache->proc_point_buf, pos_id, &pos_step);

	hair_batch_cache_fill_segments_proc_pos(curve_data, &pos_step);

	/* Create vbo immediatly to bind to texture buffer. */
	GPU_vertbuf_use(cache->proc_point_buf);

	cache->point_tex = GPU_texture_create_from_vertbuf(cache->proc_point_buf);
}

static void hair_pack_mcol(MLoopCol *mcol, ushort r_scol[4])
{
	/* Convert to linear ushort and swizzle */
	r_scol[0] = unit_float_to_ushort_clamp(BLI_color_from_srgb_table[mcol->r]);
	r_scol[1] = unit_float_to_ushort_clamp(BLI_color_from_srgb_table[mcol->g]);
	r_scol[2] = unit_float_to_ushort_clamp(BLI_color_from_srgb_table[mcol->b]);
}

static int hair_batch_cache_fill_strands_data(
        const HairCurveData *curve_data,
        GPUVertBufRaw *data_step,
        GPUVertBufRaw *uv_step, int num_uv_layers,
        GPUVertBufRaw *col_step, int num_col_layers)
{
	int curr_point = 0;
	HairIterator iter;
	HairFollicle *follicle;
	HairFiberCurve *curve;
	BKE_HAIR_ITER_FOLLICLE_CURVES(follicle, curve, &iter, curve_data) {
		if (curve->numverts < 2) {
			continue;
		}

		uint *seg_data = (uint *)GPU_vertbuf_raw_step(data_step);
		const uint numseg = curve->numverts - 1;
		*seg_data = (curr_point & 0xFFFFFF) | (numseg << 24);
		curr_point += curve->numverts;

		float (*uv)[2] = NULL;
		MCol *mcol = NULL;
		
#if 0
		particle_calculate_uvs(
				psys, psmd,
				is_simple, num_uv_layers,
				is_child ? psys->child[i].parent : i,
				is_child ? i : -1,
				mtfaces,
				*r_parent_uvs, &uv);

		particle_calculate_mcol(
				psys, psmd,
				is_simple, num_col_layers,
				is_child ? psys->child[i].parent : i,
				is_child ? i : -1,
				mcols,
				*r_parent_mcol, &mcol);
#else
		/* XXX dummy uvs and mcols, TODO */
		uv = MEM_mallocN(sizeof(*uv) * num_uv_layers, __func__);
		mcol = MEM_mallocN(sizeof(*mcol) * num_col_layers, __func__);
		for (int k = 0; k < num_uv_layers; k++) {
			zero_v3(uv[k]);
		}
		for (int k = 0; k < num_col_layers; k++) {
			mcol[k].a = 255;
			mcol[k].r = 255;
			mcol[k].g = 0;
			mcol[k].b = 255;
		}
#endif
		
//		for (int k = 0; k < num_uv_layers; k++) {
//			float *t_uv = (float *)GPU_vertbuf_raw_step(uv_step + k);
//			copy_v2_v2(t_uv, uv[k]);
//		}
//		for (int k = 0; k < num_col_layers; k++) {
//			unsigned short *scol = (unsigned short *)GPU_vertbuf_raw_step(col_step + k);
//			hair_pack_mcol(&mcol[k], scol);
//		}

		if (uv) {
			MEM_freeN(uv);
		}
		if (mcol) {
			MEM_freeN(mcol);
		}
	}
	return curr_point;
}

static void hair_batch_cache_ensure_procedural_strand_data(
        const HairCurveData *curve_data,
        ParticleHairCache *cache)
{
	int active_uv = 0;
	int active_col = 0;

#if 0 // TODO
	ParticleSystemModifierData *psmd = (ParticleSystemModifierData *)md;

	if (psmd != NULL && psmd->mesh_final != NULL) {
		if (CustomData_has_layer(&psmd->mesh_final->ldata, CD_MLOOPUV)) {
			cache->num_uv_layers = CustomData_number_of_layers(&psmd->mesh_final->ldata, CD_MLOOPUV);
			active_uv = CustomData_get_active_layer(&psmd->mesh_final->ldata, CD_MLOOPUV);
		}
		if (CustomData_has_layer(&psmd->mesh_final->ldata, CD_MLOOPCOL)) {
			cache->num_col_layers = CustomData_number_of_layers(&psmd->mesh_final->ldata, CD_MLOOPCOL);
			active_col = CustomData_get_active_layer(&psmd->mesh_final->ldata, CD_MLOOPCOL);
		}
	}
#endif

	GPUVertBufRaw data_step;
	GPUVertBufRaw uv_step[MAX_MTFACE];
	GPUVertBufRaw col_step[MAX_MCOL];

	MTFace *mtfaces[MAX_MTFACE] = {NULL};
	MCol *mcols[MAX_MCOL] = {NULL};

	GPUVertFormat format_data = {0};
	uint data_id = GPU_vertformat_attr_add(&format_data, "data", GPU_COMP_U32, 1, GPU_FETCH_INT);

	GPUVertFormat format_uv = {0};
	uint uv_id = GPU_vertformat_attr_add(&format_uv, "uv", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);

	GPUVertFormat format_col = {0};
	uint col_id = GPU_vertformat_attr_add(&format_col, "col", GPU_COMP_U16, 4, GPU_FETCH_INT_TO_FLOAT_UNIT);

	memset(cache->uv_layer_names, 0, sizeof(cache->uv_layer_names));
	memset(cache->col_layer_names, 0, sizeof(cache->col_layer_names));

	/* Strand Data */
	cache->proc_strand_buf = GPU_vertbuf_create_with_format(&format_data);
	GPU_vertbuf_data_alloc(cache->proc_strand_buf, cache->strands_len);
	GPU_vertbuf_attr_get_raw_data(cache->proc_strand_buf, data_id, &data_step);

#if 0 // TODO
	/* UV layers */
	for (int i = 0; i < cache->num_uv_layers; i++) {
		cache->proc_uv_buf[i] = GPU_vertbuf_create_with_format(&format_uv);
		GPU_vertbuf_data_alloc(cache->proc_uv_buf[i], cache->strands_len);
		GPU_vertbuf_attr_get_raw_data(cache->proc_uv_buf[i], uv_id, &uv_step[i]);

		const char *name = CustomData_get_layer_name(&psmd->mesh_final->ldata, CD_MLOOPUV, i);
		uint hash = BLI_ghashutil_strhash_p(name);
		int n = 0;
		BLI_snprintf(cache->uv_layer_names[i][n++], MAX_LAYER_NAME_LEN, "u%u", hash);
		BLI_snprintf(cache->uv_layer_names[i][n++], MAX_LAYER_NAME_LEN, "a%u", hash);

		if (i == active_uv) {
			BLI_snprintf(cache->uv_layer_names[i][n], MAX_LAYER_NAME_LEN, "u");
		}
	}
	/* Vertex colors */
	for (int i = 0; i < cache->num_col_layers; i++) {
		cache->proc_col_buf[i] = GPU_vertbuf_create_with_format(&format_col);
		GPU_vertbuf_data_alloc(cache->proc_col_buf[i], cache->strands_len);
		GPU_vertbuf_attr_get_raw_data(cache->proc_col_buf[i], col_id, &col_step[i]);

		const char *name = CustomData_get_layer_name(&psmd->mesh_final->ldata, CD_MLOOPCOL, i);
		uint hash = BLI_ghashutil_strhash_p(name);
		int n = 0;
		BLI_snprintf(cache->col_layer_names[i][n++], MAX_LAYER_NAME_LEN, "c%u", hash);

		/* We only do vcols auto name that are not overridden by uvs */
		if (CustomData_get_named_layer_index(&psmd->mesh_final->ldata, CD_MLOOPUV, name) == -1) {
			BLI_snprintf(cache->col_layer_names[i][n++], MAX_LAYER_NAME_LEN, "a%u", hash);
		}

		if (i == active_col) {
			BLI_snprintf(cache->col_layer_names[i][n], MAX_LAYER_NAME_LEN, "c");
		}
	}

	if (cache->num_uv_layers || cache->num_col_layers) {
		BKE_mesh_tessface_ensure(psmd->mesh_final);
		if (cache->num_uv_layers) {
			for (int j = 0; j < cache->num_uv_layers; j++) {
				mtfaces[j] = (MTFace *)CustomData_get_layer_n(&psmd->mesh_final->fdata, CD_MTFACE, j);
			}
		}
		if (cache->num_col_layers) {
			for (int j = 0; j < cache->num_col_layers; j++) {
				mcols[j] = (MCol *)CustomData_get_layer_n(&psmd->mesh_final->fdata, CD_MCOL, j);
			}
		}
	}
#endif

	hair_batch_cache_fill_strands_data(
	            curve_data,
	            &data_step,
	            uv_step, cache->num_uv_layers,
	            col_step, cache->num_col_layers);

	/* Create vbo immediatly to bind to texture buffer. */
	GPU_vertbuf_use(cache->proc_strand_buf);
	cache->strand_tex = GPU_texture_create_from_vertbuf(cache->proc_strand_buf);

	for (int i = 0; i < cache->num_uv_layers; i++) {
		GPU_vertbuf_use(cache->proc_uv_buf[i]);
		cache->uv_tex[i] = GPU_texture_create_from_vertbuf(cache->proc_uv_buf[i]);
	}
	for (int i = 0; i < cache->num_col_layers; i++) {
		GPU_vertbuf_use(cache->proc_col_buf[i]);
		cache->col_tex[i] = GPU_texture_create_from_vertbuf(cache->proc_col_buf[i]);
	}
}

static void hair_batch_cache_ensure_final_count(
        const HairCurveData *curve_data,
        ParticleHairFinalCache *cache,
        int subdiv,
        int thickness_res)
{
	cache->strands_len = curve_data->totfollicles;
	cache->elems_len = 0;
	cache->point_len = 0;

	HairIterator iter;
	HairFollicle *follicle;
	HairFiberCurve *curve;
	BKE_HAIR_ITER_FOLLICLE_CURVES(follicle, curve, &iter, curve_data) {
		/* +1 for primitive restart */
		cache->elems_len += (((curve->numverts - 1) << subdiv) + 2) * thickness_res;
		cache->point_len += ((curve->numverts - 1) << subdiv) + 1;
	}
}

static void hair_batch_cache_ensure_procedural_final_points(
        ParticleHairCache *cache,
        int subdiv)
{
	/* Same format as point_tex. */
	GPUVertFormat format = { 0 };
	GPU_vertformat_attr_add(&format, "pos", GPU_COMP_F32, 4, GPU_FETCH_FLOAT);

	cache->final[subdiv].proc_point_buf = GPU_vertbuf_create_with_format(&format);

	/* Create a destination buffer for the tranform feedback. Sized appropriately */
	/* Thoses are points! not line segments. */
	GPU_vertbuf_data_alloc(cache->final[subdiv].proc_point_buf, cache->final[subdiv].point_len);

	/* Create vbo immediatly to bind to texture buffer. */
	GPU_vertbuf_use(cache->final[subdiv].proc_point_buf);

	cache->final[subdiv].proc_tex = GPU_texture_create_from_vertbuf(cache->final[subdiv].proc_point_buf);
}

static int hair_batch_cache_fill_segments_indices(
        const HairCurveData *curve_data,
        const int subdiv,
        const int thickness_res,
        GPUIndexBufBuilder *elb)
{
	int curr_point = 0;
	HairIterator iter;
	HairFollicle *follicle;
	HairFiberCurve *curve;
	BKE_HAIR_ITER_FOLLICLE_CURVES(follicle, curve, &iter, curve_data) {
		if (curve->numverts < 2) {
			continue;
		}

		const int res = (((curve->numverts - 1) << subdiv) + 1) * thickness_res;
		for (int k = 0; k < res; k++) {
			GPU_indexbuf_add_generic_vert(elb, curr_point++);
		}
		GPU_indexbuf_add_primitive_restart(elb);
	}
	return curr_point;
}

static void hair_batch_cache_ensure_procedural_indices(
        const HairCurveData *curve_data,
        ParticleHairCache *cache,
        int thickness_res,
        int subdiv)
{
	BLI_assert(thickness_res <= MAX_THICKRES); /* Cylinder strip not currently supported. */

	if (cache->final[subdiv].proc_hairs[thickness_res - 1] != NULL) {
		return;
	}

	int element_count = cache->final[subdiv].elems_len;
	GPUPrimType prim_type = (thickness_res == 1) ? GPU_PRIM_LINE_STRIP : GPU_PRIM_TRI_STRIP;

	static GPUVertFormat format = { 0 };
	GPU_vertformat_clear(&format);

	/* initialize vertex format */
	GPU_vertformat_attr_add(&format, "dummy", GPU_COMP_U8, 1, GPU_FETCH_INT_TO_FLOAT_UNIT);

	GPUVertBuf *vbo = GPU_vertbuf_create_with_format(&format);
	GPU_vertbuf_data_alloc(vbo, 1);

	GPUIndexBufBuilder elb;
	GPU_indexbuf_init_ex(&elb, prim_type, element_count, element_count, true);

	hair_batch_cache_fill_segments_indices(curve_data, subdiv, thickness_res, &elb);

	cache->final[subdiv].proc_hairs[thickness_res - 1] = GPU_batch_create_ex(
	        prim_type,
	        vbo,
	        GPU_indexbuf_build(&elb),
	        GPU_BATCH_OWNS_VBO | GPU_BATCH_OWNS_INDEX);
}

/* Ensure all textures and buffers needed for GPU accelerated drawing. */
bool hair_ensure_procedural_data(
        Object *object,
        HairSystem *hsys,
        ParticleHairCache **r_hair_cache,
        int subdiv,
        int thickness_res)
{
	const DRWContextState *draw_ctx = DRW_context_state_get();
	struct Mesh *scalp = BKE_hair_get_scalp(draw_ctx->depsgraph, object, hsys);
	bool need_ft_update = false;

	HairCurveData *curve_data;
	if (hsys->edithair) {
		curve_data = &hsys->edithair->curve_data;
	}
	else {
		curve_data = &hsys->curve_data;
	}

	HairBatchCache *cache = hair_batch_cache_get(hsys);
	*r_hair_cache = &cache->hair;

	/* Refreshed on combing and simulation. */
	if (cache->hair.proc_point_buf == NULL) {
		hair_batch_cache_ensure_count(curve_data, &cache->hair);
		
		hair_batch_cache_ensure_procedural_pos(curve_data, &cache->hair);
		need_ft_update = true;
	}

	/* Refreshed if active layer or custom data changes. */
	if (cache->hair.strand_tex == NULL) {
		hair_batch_cache_ensure_procedural_strand_data(curve_data, &cache->hair);
	}

	/* Refreshed only on subdiv count change. */
	if (cache->hair.final[subdiv].proc_point_buf == NULL) {
		hair_batch_cache_ensure_final_count(curve_data, &cache->hair.final[subdiv], subdiv, thickness_res);
		
		hair_batch_cache_ensure_procedural_final_points(&cache->hair, subdiv);
		need_ft_update = true;
	}
	if (cache->hair.final[subdiv].proc_hairs[thickness_res - 1] == NULL) {
		hair_batch_cache_ensure_procedural_indices(curve_data, &cache->hair, thickness_res, subdiv);
	}

	return need_ft_update;
}

typedef struct HairScalpAttributeData
{
	uint *uv_id;
	uint *col_id;
	int num_uv_layers;
	int num_col_layers;
	int active_uv;
	int active_col;

	float (*hair_uv)[2];
	ushort (*hair_col)[4];
} HairScalpAttributeData;

static void hair_batch_cache_scalp_attributes_get(Mesh *scalp, HairScalpAttributeData *data, GPUVertFormat *format)
{
	memset(data, 0, sizeof(*data));

	if (scalp == NULL) {
		return;
	}

	if (CustomData_has_layer(&scalp->ldata, CD_MLOOPUV)) {
		data->num_uv_layers = CustomData_number_of_layers(&scalp->ldata, CD_MLOOPUV);
		data->active_uv = CustomData_get_active_layer(&scalp->ldata, CD_MLOOPUV);
	}
	if (CustomData_has_layer(&scalp->ldata, CD_MLOOPCOL)) {
		data->num_col_layers = CustomData_number_of_layers(&scalp->ldata, CD_MLOOPCOL);
		data->active_col = CustomData_get_active_layer(&scalp->ldata, CD_MLOOPCOL);
	}

	data->uv_id = MEM_mallocN(sizeof(*data->uv_id) * data->num_uv_layers, "UV attrib format");
	data->col_id = MEM_mallocN(sizeof(*data->col_id) * data->num_col_layers, "Col attrib format");

	for (int i = 0; i < data->num_uv_layers; i++) {
		const char *name = CustomData_get_layer_name(&scalp->ldata, CD_MLOOPUV, i);
		char uuid[32];

		BLI_snprintf(uuid, sizeof(uuid), "u%u", BLI_ghashutil_strhash_p(name));
		data->uv_id[i] = GPU_vertformat_attr_add(format, uuid, GPU_COMP_F32, 2, GPU_FETCH_FLOAT);

		if (i == data->active_uv) {
			GPU_vertformat_alias_add(format, "u");
		}
	}

	for (int i = 0; i < data->num_col_layers; i++) {
		const char *name = CustomData_get_layer_name(&scalp->ldata, CD_MLOOPCOL, i);
		char uuid[32];

		BLI_snprintf(uuid, sizeof(uuid), "c%u", BLI_ghashutil_strhash_p(name));
		data->col_id[i] = GPU_vertformat_attr_add(format, uuid, GPU_COMP_F32, 2, GPU_FETCH_FLOAT);

		if (i == data->active_col) {
			GPU_vertformat_alias_add(format, "c");
		}
	}

	data->hair_uv = MEM_callocN(sizeof(*data->hair_uv) * data->num_uv_layers, "Hair UVs");
	data->hair_col = MEM_callocN(sizeof(*data->hair_col) * data->num_col_layers, "Hair MCol");
}

static void hair_batch_cache_scalp_attributes_free(HairScalpAttributeData *data)
{
	MEM_SAFE_FREE(data->uv_id);
	MEM_SAFE_FREE(data->col_id);
	MEM_SAFE_FREE(data->hair_uv);
	MEM_SAFE_FREE(data->hair_col);
}

typedef struct HairAttributeID {
	uint pos;
	uint tan;
	uint ind;
} HairAttributeID;

static int hair_batch_cache_fill_segments(
        const HairCurveData *curve_data,
        const Mesh *scalp,
        const HairScalpAttributeData *scalp_attr,
        GPUIndexBufBuilder *elb,
        HairAttributeID *attr_id,
        ParticleHairCache *hair_cache)
{
	int curr_point = 0;
	const HairFiberVertex *verts = BKE_hair_get_verts(curve_data);
	float (*orco)[3] = CustomData_get_layer(&curve_data->fdata, CD_ORCO);
	float (*normals)[3] = CustomData_get_layer(&curve_data->fdata, CD_NORMAL);
	float (*tangents)[3] = CustomData_get_layer(&curve_data->fdata, CD_TANGENT);
	const HairFollicle *follicle;
	const HairFiberCurve *curve;
	HairIterator iter;
	int curve_idx;
	BKE_HAIR_ITER_FOLLICLE_CURVES_INDEX(follicle, curve, &iter, curve_data, curve_idx) {
		if (curve->numverts < 2) {
			continue;
		}

		float hmat[4][4];
		hair_get_local_space(orco, normals, tangents, curve_idx, hmat);
		for (int k = 0; k < scalp_attr->num_uv_layers; k++) {
			BKE_mesh_sample_eval_uv(scalp, &follicle->mesh_sample, k, scalp_attr->hair_uv[k]);
		}
		for (int k = 0; k < scalp_attr->num_col_layers; k++) {
			MLoopCol col;
			BKE_mesh_sample_eval_col(scalp, &follicle->mesh_sample, k, &col);
			hair_pack_mcol(&col, scalp_attr->hair_col[k]);
		}

		const HairFiberVertex *vert_start = &verts[curve->vertstart];
		for (int j = 0; j < curve->numverts; ++j) {
			const HairFiberVertex *vert = vert_start + j;
			const HairFiberVertex *vert_prev = (j > 0 ? vert_start + j - 1 : vert_start + j);
			const HairFiberVertex *vert_next = (j < curve->numverts-1 ? vert_start + j + 1 : vert_start + j);
			float co[3];
			mul_v3_m4v3(co, hmat, vert->co);
			float tangent[3];
			sub_v3_v3v3(tangent, vert_next->co, vert_prev->co);
			mul_mat3_m4_v3(hmat, tangent);

			GPU_vertbuf_attr_set(hair_cache->pos, attr_id->pos, curr_point, co);
			GPU_vertbuf_attr_set(hair_cache->pos, attr_id->tan, curr_point, tangent);
			GPU_vertbuf_attr_set(hair_cache->pos, attr_id->ind, curr_point, &j);

			for (int k = 0; k < scalp_attr->num_uv_layers; k++) {
				GPU_vertbuf_attr_set(hair_cache->pos, scalp_attr->uv_id[k], curr_point, scalp_attr->hair_uv[k]);
			}
			for (int k = 0; k < scalp_attr->num_col_layers; k++) {
				GPU_vertbuf_attr_set(hair_cache->pos, scalp_attr->col_id[k], curr_point, scalp_attr->hair_col[k]);
			}

			GPU_indexbuf_add_generic_vert(elb, curr_point);
			++curr_point;
		}

		/* Add restart primitive. */
		GPU_indexbuf_add_primitive_restart(elb);
	}

	return curr_point;
}

static void hair_batch_cache_ensure_pos_and_seg(
        const HairCurveData *curve_data,
        Mesh *scalp,
        ParticleHairCache *hair_cache)
{
	if (hair_cache->pos != NULL && hair_cache->indices != NULL) {
		return;
	}

	GPU_VERTBUF_DISCARD_SAFE(hair_cache->pos);
	GPU_INDEXBUF_DISCARD_SAFE(hair_cache->indices);

	static GPUVertFormat format = { 0 };
	HairAttributeID attr_id;
	HairScalpAttributeData scalp_attr;

	GPU_vertformat_clear(&format);

	/* initialize vertex format */
	attr_id.pos = GPU_vertformat_attr_add(&format, "pos", GPU_COMP_F32, 3, GPU_FETCH_FLOAT);
	attr_id.tan = GPU_vertformat_attr_add(&format, "nor", GPU_COMP_F32, 3, GPU_FETCH_FLOAT);
	attr_id.ind = GPU_vertformat_attr_add(&format, "ind", GPU_COMP_I32, 1, GPU_FETCH_INT);
	hair_batch_cache_scalp_attributes_get(scalp, &scalp_attr, &format);

	hair_cache->pos = GPU_vertbuf_create_with_format(&format);
	GPU_vertbuf_data_alloc(hair_cache->pos, hair_cache->point_len);

	GPUIndexBufBuilder elb;
	GPU_indexbuf_init_ex(
	        &elb,
	        GPU_PRIM_LINE_STRIP,
	        hair_cache->elems_len, hair_cache->point_len,
	        true);

	hair_batch_cache_fill_segments(curve_data, scalp, &scalp_attr, &elb, &attr_id, hair_cache);

	/* Cleanup. */
	hair_batch_cache_scalp_attributes_free(&scalp_attr);

	hair_cache->indices = GPU_indexbuf_build(&elb);
}

static void ensure_edit_follicle_points_count(
        const HairCurveData *curve_data,
        HairBatchCache *cache)
{
	if (cache->edit_follicle_pos != NULL) {
		return;
	}
	cache->edit_follicle_point_len = curve_data->totfollicles;
}

static void hair_batch_cache_ensure_edit_follicle_pos(
        const HairCurveData *curve_data,
        const Mesh *scalp,
        HairBatchCache *cache)
{
	if (cache->edit_follicle_pos != NULL) {
		return;
	}

	static GPUVertFormat format = { 0 };
	static uint pos_id, data_id;

	GPU_VERTBUF_DISCARD_SAFE(cache->edit_follicle_pos);

	if (format.attr_len == 0) {
		/* initialize vertex format */
		pos_id = GPU_vertformat_attr_add(&format, "pos", GPU_COMP_F32, 3, GPU_FETCH_FLOAT);
		data_id = GPU_vertformat_attr_add(&format, "data", GPU_COMP_U32, 1, GPU_FETCH_INT);
	}

	cache->edit_follicle_pos = GPU_vertbuf_create_with_format(&format);
	GPU_vertbuf_data_alloc(cache->edit_follicle_pos, cache->edit_follicle_point_len);

	const HairFollicle *follicle;
	HairIterator iter;
	int point_index;
	BKE_HAIR_ITER_FOLLICLES_INDEX(follicle, &iter, curve_data, point_index) {
		float loc[3];
		BKE_mesh_sample_eval(scalp, &follicle->mesh_sample, loc, NULL, NULL);

		GPU_vertbuf_attr_set(cache->edit_follicle_pos, pos_id, point_index, loc);

		GPU_vertbuf_attr_set(cache->edit_follicle_pos, data_id, point_index, &follicle->flag);
	}
}

GPUBatch *DRW_hair_batch_cache_get_edit_follicle_points(Object *ob, HairSystem *hsys)
{
	HairBatchCache *cache = hair_batch_cache_get(hsys);
	if (cache->edit_follicle_points != NULL) {
		return cache->edit_follicle_points;
	}

	Mesh *scalp = BKE_hair_get_scalp(DRW_context_state_get()->depsgraph, ob, hsys);
	const HairCurveData *curve_data;
	if (hsys->edithair) {
		curve_data = &hsys->edithair->curve_data;
	}
	else {
		curve_data = &hsys->curve_data;
	}

	ensure_edit_follicle_points_count(curve_data, cache);
	hair_batch_cache_ensure_edit_follicle_pos(curve_data, scalp, cache);
	cache->edit_follicle_points = GPU_batch_create(
	        GPU_PRIM_POINTS,
	        cache->edit_follicle_pos,
	        NULL);
	return cache->edit_follicle_points;
}

GPUBatch *DRW_hair_batch_cache_get_edit_follicle_normals(Object *ob, HairSystem *hsys)
{
	// TODO
	UNUSED_VARS(ob, hsys);
	return NULL;
}

GPUBatch *DRW_hair_batch_cache_get_edit_follicle_axes(Object *ob, HairSystem *hsys)
{
	// TODO
	UNUSED_VARS(ob, hsys);
	return NULL;
}

GPUBatch *DRW_hair_batch_cache_get_edit_strands(Object *ob, HairSystem *hsys)
{
	HairBatchCache *cache = hair_batch_cache_get(hsys);
	if (cache->edit_hair.hairs != NULL) {
		return cache->edit_hair.hairs;
	}

	Mesh *scalp = BKE_hair_get_scalp(DRW_context_state_get()->depsgraph, ob, hsys);
	const HairCurveData *curve_data;
	if (hsys->edithair) {
		curve_data = &hsys->edithair->curve_data;
	}
	else {
		curve_data = &hsys->curve_data;
	}
	hair_batch_cache_ensure_count(curve_data, &cache->edit_hair);
	hair_batch_cache_ensure_pos_and_seg(curve_data, scalp, &cache->edit_hair);
	cache->edit_hair.hairs = GPU_batch_create(
	        GPU_PRIM_LINE_STRIP,
	        cache->edit_hair.pos,
	        cache->edit_hair.indices);
	return cache->edit_hair.hairs;
}
