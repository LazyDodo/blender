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
#include "DNA_scene_types.h"

#include "BKE_hair.h"
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
/* Hair Gwn_Batch Cache */

/* Gwn_Batch cache management. */

typedef struct HairBatchCache {
	ParticleHairCache hair;

	bool is_dirty;
} HairBatchCache;

static void hair_batch_cache_clear(HairSystem *hsys);

static bool hair_batch_cache_valid(HairSystem *hsys)
{
	HairBatchCache *cache = hsys->draw_batch_cache;

	if (cache == NULL) {
		return false;
	}

	if (cache->is_dirty) {
		return false;
	}

	return true;
}

static void hair_batch_cache_init(HairSystem *hsys)
{
	HairBatchCache *cache = hsys->draw_batch_cache;
	
	if (!cache) {
		cache = hsys->draw_batch_cache = MEM_callocN(sizeof(*cache), __func__);
	}
	else {
		memset(cache, 0, sizeof(*cache));
	}
	
	cache->is_dirty = false;
}

static HairBatchCache *hair_batch_cache_get(HairSystem *hsys)
{
	// Hair follicle binding needs to be updated after changes
	BLI_assert(!(hsys->flag & HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING));
	
	if (!hair_batch_cache_valid(hsys)) {
		hair_batch_cache_clear(hsys);
		hair_batch_cache_init(hsys);
	}
	return hsys->draw_batch_cache;
}

void DRW_hair_batch_cache_dirty(HairSystem *hsys, int mode)
{
	HairBatchCache *cache = hsys->draw_batch_cache;
	if (cache == NULL) {
		return;
	}
	switch (mode) {
		case BKE_HAIR_BATCH_DIRTY_ALL:
			cache->is_dirty = true;
			break;
		default:
			BLI_assert(0);
	}
}

static void hair_batch_cache_clear(HairSystem *hsys)
{
	HairBatchCache *cache = hsys->draw_batch_cache;
	particle_batch_cache_clear_hair(&cache->hair);
}

void DRW_hair_batch_cache_free(HairSystem *hsys)
{
	hair_batch_cache_clear(hsys);
	MEM_SAFE_FREE(hsys->draw_batch_cache);
}

static void ensure_seg_pt_count(
        const HairExportCache *hair_export,
        ParticleHairCache *cache)
{
	if ((cache->pos != NULL && cache->indices != NULL) ||
	    (cache->proc_point_buf != NULL))
	{
		return;
	}

	cache->strands_count = hair_export->totfollicles;
	cache->elems_count = hair_export->totverts - hair_export->totcurves;
	cache->point_count = hair_export->totverts;
}

static void hair_batch_cache_fill_segments_proc_pos(
        const HairExportCache *hair_export,
        Gwn_VertBufRaw *attr_step)
{
	for (int i = 0; i < hair_export->totcurves; i++) {
		const HairFiberCurve *curve = &hair_export->fiber_curves[i];
		const HairFiberVertex *verts = &hair_export->fiber_verts[curve->vertstart];
		if (curve->numverts < 2) {
			continue;
		}
		float total_len = 0.0f;
		const float *co_prev = NULL;
		float *seg_data_first;
		for (int j = 0; j < curve->numverts; j++) {
			float *seg_data = (float *)GWN_vertbuf_raw_step(attr_step);
			copy_v3_v3(seg_data, verts[j].co);
			if (co_prev) {
				total_len += len_v3v3(co_prev, verts[j].co);
			}
			else {
				seg_data_first = seg_data;
			}
			seg_data[3] = total_len;
			co_prev = verts[j].co;
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
        const HairExportCache *hair_export,
        ParticleHairCache *cache)
{
	if (cache->proc_point_buf != NULL) {
		return;
	}

	/* initialize vertex format */
	Gwn_VertFormat format = {0};
	uint pos_id = GWN_vertformat_attr_add(&format, "posTime", GWN_COMP_F32, 4, GWN_FETCH_FLOAT);

	cache->proc_point_buf = GWN_vertbuf_create_with_format(&format);
	GWN_vertbuf_data_alloc(cache->proc_point_buf, cache->point_count);

	Gwn_VertBufRaw pos_step;
	GWN_vertbuf_attr_get_raw_data(cache->proc_point_buf, pos_id, &pos_step);

	hair_batch_cache_fill_segments_proc_pos(hair_export, &pos_step);

	/* Create vbo immediatly to bind to texture buffer. */
	GWN_vertbuf_use(cache->proc_point_buf);

	cache->point_tex = GPU_texture_create_from_vertbuf(cache->proc_point_buf);
}

static void hair_batch_cache_ensure_procedural_strand_data(
        const HairExportCache *hair_export,
        ParticleHairCache *cache)
{
#if 0 // TODO
	int active_uv = 0;
	int active_col = 0;

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

	Gwn_VertBufRaw data_step;
	Gwn_VertBufRaw uv_step[MAX_MTFACE];
	Gwn_VertBufRaw col_step[MAX_MCOL];

	MTFace *mtfaces[MAX_MTFACE] = {NULL};
	MCol *mcols[MAX_MCOL] = {NULL};
	float (**parent_uvs)[2] = NULL;
	MCol **parent_mcol = NULL;

	Gwn_VertFormat format_data = {0};
	uint data_id = GWN_vertformat_attr_add(&format_data, "data", GWN_COMP_U32, 1, GWN_FETCH_INT);

	Gwn_VertFormat format_uv = {0};
	uint uv_id = GWN_vertformat_attr_add(&format_uv, "uv", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);

	Gwn_VertFormat format_col = {0};
	uint col_id = GWN_vertformat_attr_add(&format_col, "col", GWN_COMP_U16, 4, GWN_FETCH_INT_TO_FLOAT_UNIT);

	memset(cache->uv_layer_names, 0, sizeof(cache->uv_layer_names));
	memset(cache->col_layer_names, 0, sizeof(cache->col_layer_names));

	/* Strand Data */
	cache->proc_strand_buf = GWN_vertbuf_create_with_format(&format_data);
	GWN_vertbuf_data_alloc(cache->proc_strand_buf, cache->strands_count);
	GWN_vertbuf_attr_get_raw_data(cache->proc_strand_buf, data_id, &data_step);

	/* UV layers */
	for (int i = 0; i < cache->num_uv_layers; i++) {
		cache->proc_uv_buf[i] = GWN_vertbuf_create_with_format(&format_uv);
		GWN_vertbuf_data_alloc(cache->proc_uv_buf[i], cache->strands_count);
		GWN_vertbuf_attr_get_raw_data(cache->proc_uv_buf[i], uv_id, &uv_step[i]);

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
		cache->proc_col_buf[i] = GWN_vertbuf_create_with_format(&format_col);
		GWN_vertbuf_data_alloc(cache->proc_col_buf[i], cache->strands_count);
		GWN_vertbuf_attr_get_raw_data(cache->proc_col_buf[i], col_id, &col_step[i]);

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

	if (edit != NULL && edit->pathcache != NULL) {
		particle_batch_cache_fill_strands_data(
		        psys, psmd, edit->pathcache, PARTICLE_SOURCE_PARENT,
		        0, edit->totcached,
		        &data_step,
		        &parent_uvs, uv_step, (MTFace **)mtfaces, cache->num_uv_layers,
		        &parent_mcol, col_step, (MCol **)mcols, cache->num_col_layers);
	}
	else {
		int curr_point = 0;
		if ((psys->pathcache != NULL) &&
		    (!psys->childcache || (psys->part->draw & PART_DRAW_PARENT)))
		{
			curr_point = particle_batch_cache_fill_strands_data(
			        psys, psmd, psys->pathcache, PARTICLE_SOURCE_PARENT,
			        0, psys->totpart,
			        &data_step,
			        &parent_uvs, uv_step, (MTFace **)mtfaces, cache->num_uv_layers,
			        &parent_mcol, col_step, (MCol **)mcols, cache->num_col_layers);
		}
		if (psys->childcache) {
			const int child_count = psys->totchild * psys->part->disp / 100;
			curr_point = particle_batch_cache_fill_strands_data(
			        psys, psmd, psys->childcache, PARTICLE_SOURCE_CHILDREN,
			        curr_point, child_count,
			        &data_step,
			        &parent_uvs, uv_step, (MTFace **)mtfaces, cache->num_uv_layers,
			        &parent_mcol, col_step, (MCol **)mcols, cache->num_col_layers);
		}
	}
	/* Cleanup. */
	if (parent_uvs != NULL) {
		/* TODO(sergey): For edit mode it should be edit->totcached. */
		for (int i = 0; i < psys->totpart; i++) {
			MEM_SAFE_FREE(parent_uvs[i]);
		}
		MEM_freeN(parent_uvs);
	}
	if (parent_mcol != NULL) {
		for (int i = 0; i < psys->totpart; i++) {
			MEM_SAFE_FREE(parent_mcol[i]);
		}
		MEM_freeN(parent_mcol);
	}
#else
	UNUSED_VARS(hair_export);
#endif

	/* Create vbo immediatly to bind to texture buffer. */
	GWN_vertbuf_use(cache->proc_strand_buf);
	cache->strand_tex = GPU_texture_create_from_vertbuf(cache->proc_strand_buf);

	for (int i = 0; i < cache->num_uv_layers; i++) {
		GWN_vertbuf_use(cache->proc_uv_buf[i]);
		cache->uv_tex[i] = GPU_texture_create_from_vertbuf(cache->proc_uv_buf[i]);
	}
	for (int i = 0; i < cache->num_col_layers; i++) {
		GWN_vertbuf_use(cache->proc_col_buf[i]);
		cache->col_tex[i] = GPU_texture_create_from_vertbuf(cache->proc_col_buf[i]);
	}
}

static void hair_batch_cache_ensure_procedural_final_points(
        ParticleHairCache *cache,
        int subdiv)
{
	/* Same format as point_tex. */
	Gwn_VertFormat format = { 0 };
	GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 4, GWN_FETCH_FLOAT);

	cache->final[subdiv].proc_buf = GWN_vertbuf_create_with_format(&format);

	/* Create a destination buffer for the tranform feedback. Sized appropriately */
	/* Thoses are points! not line segments. */
	GWN_vertbuf_data_alloc(cache->final[subdiv].proc_buf, cache->final[subdiv].strands_res * cache->strands_count);

	/* Create vbo immediatly to bind to texture buffer. */
	GWN_vertbuf_use(cache->final[subdiv].proc_buf);

	cache->final[subdiv].proc_tex = GPU_texture_create_from_vertbuf(cache->final[subdiv].proc_buf);
}

static int hair_batch_cache_fill_segments_indices(
        const HairExportCache *hair_export,
        const int res,
        Gwn_IndexBufBuilder *elb)
{
	int curr_point = 0;
	for (int i = 0; i < hair_export->totcurves; i++) {
		const HairFiberCurve *curve = &hair_export->fiber_curves[i];
		if (curve->numverts < 2) {
			continue;
		}
		for (int k = 0; k < res; k++) {
			GWN_indexbuf_add_generic_vert(elb, curr_point++);
		}
		GWN_indexbuf_add_primitive_restart(elb);
	}
	return curr_point;
}

static void hair_batch_cache_ensure_procedural_indices(
        const HairExportCache *hair_export,
        ParticleHairCache *cache,
        int thickness_res,
        int subdiv)
{
	BLI_assert(thickness_res <= MAX_THICKRES); /* Cylinder strip not currently supported. */

	if (cache->final[subdiv].proc_hairs[thickness_res - 1] != NULL) {
		return;
	}

	int verts_per_hair = cache->final[subdiv].strands_res * thickness_res;
	/* +1 for primitive restart */
	int element_count = (verts_per_hair + 1) * cache->strands_count;
	Gwn_PrimType prim_type = (thickness_res == 1) ? GWN_PRIM_LINE_STRIP : GWN_PRIM_TRI_STRIP;

	static Gwn_VertFormat format = { 0 };
	GWN_vertformat_clear(&format);

	/* initialize vertex format */
	GWN_vertformat_attr_add(&format, "dummy", GWN_COMP_U8, 1, GWN_FETCH_INT_TO_FLOAT_UNIT);

	Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);
	GWN_vertbuf_data_alloc(vbo, 1);

	Gwn_IndexBufBuilder elb;
	GWN_indexbuf_init_ex(&elb, prim_type, element_count, element_count, true);

	hair_batch_cache_fill_segments_indices(hair_export, verts_per_hair, &elb);

	cache->final[subdiv].proc_hairs[thickness_res - 1] = GWN_batch_create_ex(
	        prim_type,
	        vbo,
	        GWN_indexbuf_build(&elb),
	        GWN_BATCH_OWNS_VBO | GWN_BATCH_OWNS_INDEX);
}

/* Ensure all textures and buffers needed for GPU accelerated drawing. */
bool hair_ensure_procedural_data(
        Object *UNUSED(object),
        HairSystem *hsys,
        struct Mesh *scalp,
        ParticleHairCache **r_hair_cache,
        int subdiv,
        int thickness_res)
{
	bool need_ft_update = false;

	HairExportCache *hair_export = BKE_hair_export_cache_new();
	BKE_hair_export_cache_update(hair_export, hsys, subdiv, scalp, HAIR_EXPORT_ALL);

	HairBatchCache *cache = hair_batch_cache_get(hsys);
	*r_hair_cache = &cache->hair;

	const int hsys_subdiv = 0; // XXX TODO per-hsys or per-fiber subdiv
	cache->hair.final[subdiv].strands_res = 1 << (hsys_subdiv + subdiv);

	/* Refreshed on combing and simulation. */
	if (cache->hair.proc_point_buf == NULL) {
		ensure_seg_pt_count(hair_export, &cache->hair);
		hair_batch_cache_ensure_procedural_pos(hair_export, &cache->hair);
		need_ft_update = true;
	}

	/* Refreshed if active layer or custom data changes. */
	if ((*r_hair_cache)->strand_tex == NULL) {
		hair_batch_cache_ensure_procedural_strand_data(hair_export, &cache->hair);
	}

	/* Refreshed only on subdiv count change. */
	if ((*r_hair_cache)->final[subdiv].proc_buf == NULL) {
		hair_batch_cache_ensure_procedural_final_points(&cache->hair, subdiv);
		need_ft_update = true;
	}
	if ((*r_hair_cache)->final[subdiv].proc_hairs[thickness_res - 1] == NULL) {
		hair_batch_cache_ensure_procedural_indices(hair_export, &cache->hair, thickness_res, subdiv);
	}

	BKE_hair_export_cache_free(hair_export);

	return need_ft_update;
}
