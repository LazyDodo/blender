/*
 * Copyright 2016, Blender Foundation.
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
 * Contributor(s): Blender Institute
 *
 */

/** \file eevee_lightprobes.c
 *  \ingroup draw_engine
 */

#include "DRW_render.h"

#include "BLI_utildefines.h"
#include "BLI_string_utils.h"
#include "BLI_rand.h"

#include "DNA_world_types.h"
#include "DNA_texture_types.h"
#include "DNA_image_types.h"
#include "DNA_lightprobe_types.h"
#include "DNA_view3d_types.h"

#include "BKE_collection.h"
#include "BKE_object.h"
#include "MEM_guardedalloc.h"

#include "GPU_material.h"
#include "GPU_texture.h"
#include "GPU_glew.h"

#include "DEG_depsgraph_query.h"

#include "eevee_engine.h"
#include "eevee_private.h"

#include "ED_screen.h"

#define HAMMERSLEY_SIZE 1024

static struct {
	struct GPUShader *probe_default_sh;
	struct GPUShader *probe_default_studiolight_sh;
	struct GPUShader *probe_filter_glossy_sh;
	struct GPUShader *probe_filter_diffuse_sh;
	struct GPUShader *probe_filter_visibility_sh;
	struct GPUShader *probe_grid_fill_sh;
	struct GPUShader *probe_grid_display_sh;
	struct GPUShader *probe_planar_display_sh;
	struct GPUShader *probe_planar_downsample_sh;
	struct GPUShader *probe_cube_display_sh;

	struct GPUTexture *hammersley;
	struct GPUTexture *planar_pool_placeholder;
	struct GPUTexture *depth_placeholder;
	struct GPUTexture *depth_array_placeholder;
	struct GPUTexture *cube_face_minmaxz;

	struct Gwn_VertFormat *format_probe_display_cube;
	struct Gwn_VertFormat *format_probe_display_planar;

	EEVEE_LightCache dummy_cache;
} e_data = {NULL}; /* Engine data */

extern char datatoc_background_vert_glsl[];
extern char datatoc_default_world_frag_glsl[];
extern char datatoc_lightprobe_filter_glossy_frag_glsl[];
extern char datatoc_lightprobe_filter_diffuse_frag_glsl[];
extern char datatoc_lightprobe_filter_visibility_frag_glsl[];
extern char datatoc_lightprobe_geom_glsl[];
extern char datatoc_lightprobe_vert_glsl[];
extern char datatoc_lightprobe_planar_display_frag_glsl[];
extern char datatoc_lightprobe_planar_display_vert_glsl[];
extern char datatoc_lightprobe_planar_downsample_frag_glsl[];
extern char datatoc_lightprobe_planar_downsample_geom_glsl[];
extern char datatoc_lightprobe_planar_downsample_vert_glsl[];
extern char datatoc_lightprobe_cube_display_frag_glsl[];
extern char datatoc_lightprobe_cube_display_vert_glsl[];
extern char datatoc_lightprobe_grid_display_frag_glsl[];
extern char datatoc_lightprobe_grid_display_vert_glsl[];
extern char datatoc_lightprobe_grid_fill_frag_glsl[];
extern char datatoc_irradiance_lib_glsl[];
extern char datatoc_lightprobe_lib_glsl[];
extern char datatoc_octahedron_lib_glsl[];
extern char datatoc_bsdf_common_lib_glsl[];
extern char datatoc_common_uniforms_lib_glsl[];
extern char datatoc_common_view_lib_glsl[];
extern char datatoc_bsdf_sampling_lib_glsl[];

extern GlobalsUboStorage ts;

/* *********** FUNCTIONS *********** */

void EEVEE_lightprobes_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
	EEVEE_CommonUniformBuffer *common_data = &sldata->common_data;
	EEVEE_TextureList *txl = vedata->txl;
	EEVEE_StorageList *stl = vedata->stl;

	const DRWContextState *draw_ctx = DRW_context_state_get();
	const Scene *scene_eval = DEG_get_evaluated_scene(draw_ctx->depsgraph);

	/* XXX TODO remove only for testing */
	if (e_data.dummy_cache.grid_tx == NULL) {
#ifdef IRRADIANCE_SH_L2
		/* we need a signed format for Spherical Harmonics */
		int irradiance_format = GPU_RGBA16F;
#else
		int irradiance_format = GPU_RGBA8;
#endif
		float rgba[4] = {0.0f, 1.0f, 0.0f, 0.0f};
		e_data.dummy_cache.grid_tx = DRW_texture_create_2D_array(1, 1, 1, irradiance_format, DRW_TEX_FILTER, rgba);
		e_data.dummy_cache.cube_tx = DRW_texture_create_2D_array(1, 1, 1, GPU_RGBA8, DRW_TEX_FILTER, rgba);
		e_data.dummy_cache.cube_data = MEM_callocN(sizeof(EEVEE_LightProbe), "EEVEE Cube Data Cache");
		e_data.dummy_cache.grid_data = MEM_callocN(sizeof(EEVEE_LightGrid), "EEVEE Grid Data Cache");
		e_data.dummy_cache.cube_count = 1;
		e_data.dummy_cache.grid_count = 1;
	}

	if (!sldata->probes) {
		sldata->probes = MEM_callocN(sizeof(EEVEE_LightProbesInfo), "EEVEE_LightProbesInfo");
		sldata->probe_ubo = DRW_uniformbuffer_create(sizeof(EEVEE_LightProbe) * MAX_PROBE, NULL);
		sldata->grid_ubo = DRW_uniformbuffer_create(sizeof(EEVEE_LightGrid) * MAX_GRID, NULL);
		sldata->planar_ubo = DRW_uniformbuffer_create(sizeof(EEVEE_PlanarReflection) * MAX_PLANAR, NULL);
	}

	if (scene_eval->eevee.light_cache) {
		stl->g_data->light_cache = scene_eval->eevee.light_cache;
	}
	else {
		stl->g_data->light_cache = &e_data.dummy_cache;
	}

	common_data->prb_num_planar = 0;
	common_data->prb_num_render_cube = 1;
	common_data->prb_num_render_grid = 1;

	common_data->spec_toggle = true;
	common_data->ssr_toggle = true;
	common_data->sss_toggle = true;

	if (!txl->planar_pool) {
		/* Makes Opengl Happy : Create a placeholder texture that will never be sampled but still bound to shader. */
		txl->planar_pool = DRW_texture_create_2D_array(1, 1, 1, GPU_RGBA8, DRW_TEX_FILTER | DRW_TEX_MIPMAP, NULL);
		txl->planar_depth = DRW_texture_create_2D_array(1, 1, 1, GPU_DEPTH_COMPONENT24, 0, NULL);
	}
}

void EEVEE_lightprobes_cache_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
}

void EEVEE_lightprobes_cache_add(EEVEE_ViewLayerData *sldata, Object *ob)
{
}

void EEVEE_lightprobes_grid_data_from_object(Object *ob, EEVEE_LightGrid *egrid, int *offset)
{
	LightProbe *probe = (LightProbe *)ob->data;

	copy_v3_v3_int(egrid->resolution, &probe->grid_resolution_x);

	/* Save current offset and advance it for the next grid. */
	egrid->offset = *offset;
	*offset += egrid->resolution[0] * egrid->resolution[1] * egrid->resolution[2];

	/* Add one for level 0 */
	float fac = 1.0f / max_ff(1e-8f, probe->falloff);
	egrid->attenuation_scale = fac / max_ff(1e-8f, probe->distinf);
	egrid->attenuation_bias = fac;

	/* Update transforms */
	float cell_dim[3], half_cell_dim[3];
	cell_dim[0] = 2.0f / egrid->resolution[0];
	cell_dim[1] = 2.0f / egrid->resolution[1];
	cell_dim[2] = 2.0f / egrid->resolution[2];

	mul_v3_v3fl(half_cell_dim, cell_dim, 0.5f);

	/* Matrix converting world space to cell ranges. */
	invert_m4_m4(egrid->mat, ob->obmat);

	/* First cell. */
	copy_v3_fl(egrid->corner, -1.0f);
	add_v3_v3(egrid->corner, half_cell_dim);
	mul_m4_v3(ob->obmat, egrid->corner);

	/* Opposite neighbor cell. */
	copy_v3_fl3(egrid->increment_x, cell_dim[0], 0.0f, 0.0f);
	add_v3_v3(egrid->increment_x, half_cell_dim);
	add_v3_fl(egrid->increment_x, -1.0f);
	mul_m4_v3(ob->obmat, egrid->increment_x);
	sub_v3_v3(egrid->increment_x, egrid->corner);

	copy_v3_fl3(egrid->increment_y, 0.0f, cell_dim[1], 0.0f);
	add_v3_v3(egrid->increment_y, half_cell_dim);
	add_v3_fl(egrid->increment_y, -1.0f);
	mul_m4_v3(ob->obmat, egrid->increment_y);
	sub_v3_v3(egrid->increment_y, egrid->corner);

	copy_v3_fl3(egrid->increment_z, 0.0f, 0.0f, cell_dim[2]);
	add_v3_v3(egrid->increment_z, half_cell_dim);
	add_v3_fl(egrid->increment_z, -1.0f);
	mul_m4_v3(ob->obmat, egrid->increment_z);
	sub_v3_v3(egrid->increment_z, egrid->corner);

	/* Visibility bias */
	egrid->visibility_bias = 0.05f * probe->vis_bias;
	egrid->visibility_bleed = probe->vis_bleedbias;
	egrid->visibility_range = 1.0f + sqrtf(max_fff(len_squared_v3(egrid->increment_x),
	                                               len_squared_v3(egrid->increment_y),
	                                               len_squared_v3(egrid->increment_z)));
}

void EEVEE_lightprobes_cube_data_from_object(Object *ob, EEVEE_LightProbe *eprobe)
{
	LightProbe *probe = (LightProbe *)ob->data;

	/* Update transforms */
	copy_v3_v3(eprobe->position, ob->obmat[3]);

	/* Attenuation */
	eprobe->attenuation_type = probe->attenuation_type;
	eprobe->attenuation_fac = 1.0f / max_ff(1e-8f, probe->falloff);

	unit_m4(eprobe->attenuationmat);
	scale_m4_fl(eprobe->attenuationmat, probe->distinf);
	mul_m4_m4m4(eprobe->attenuationmat, ob->obmat, eprobe->attenuationmat);
	invert_m4(eprobe->attenuationmat);

	/* Parallax */
	unit_m4(eprobe->parallaxmat);

	if ((probe->flag & LIGHTPROBE_FLAG_CUSTOM_PARALLAX) != 0) {
		eprobe->parallax_type = probe->parallax_type;
		scale_m4_fl(eprobe->parallaxmat, probe->distpar);
	}
	else {
		eprobe->parallax_type = probe->attenuation_type;
		scale_m4_fl(eprobe->parallaxmat, probe->distinf);
	}

	mul_m4_m4m4(eprobe->parallaxmat, ob->obmat, eprobe->parallaxmat);
	invert_m4(eprobe->parallaxmat);
}

static void eevee_lightprobes_extract_from_cache(EEVEE_LightProbesInfo *pinfo, EEVEE_LightCache *lcache)
{
	/* copy the entire cache for now (up to MAX_PROBE) */
	memcpy(pinfo->probe_data, lcache->cube_data, sizeof(EEVEE_LightProbe) * min_ii(lcache->cube_count, MAX_PROBE));
	/* TODO compute the max number of grid based on sample count. */
	memcpy(pinfo->grid_data, lcache->grid_data, sizeof(EEVEE_LightGrid) * min_ii(lcache->grid_count, MAX_GRID));
}

void EEVEE_lightprobes_refresh_planar(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
}

bool EEVEE_lightprobes_all_probes_ready(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
	return true;
}

void EEVEE_lightprobes_cache_finish(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
	EEVEE_StorageList *stl = vedata->stl;

	eevee_lightprobes_extract_from_cache(sldata->probes, stl->g_data->light_cache);

	DRW_uniformbuffer_update(sldata->probe_ubo, &sldata->probes->probe_data);
	DRW_uniformbuffer_update(sldata->grid_ubo, &sldata->probes->grid_data);
}

void EEVEE_lightprobes_refresh_world(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
}

void EEVEE_lightprobes_refresh(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
}

void EEVEE_lightprobes_free(void)
{
	MEM_SAFE_FREE(e_data.format_probe_display_cube);
	MEM_SAFE_FREE(e_data.format_probe_display_planar);
	DRW_SHADER_FREE_SAFE(e_data.probe_default_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_default_studiolight_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_filter_glossy_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_filter_diffuse_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_filter_visibility_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_grid_fill_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_grid_display_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_planar_display_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_planar_downsample_sh);
	DRW_SHADER_FREE_SAFE(e_data.probe_cube_display_sh);
	DRW_TEXTURE_FREE_SAFE(e_data.hammersley);
	DRW_TEXTURE_FREE_SAFE(e_data.planar_pool_placeholder);
	DRW_TEXTURE_FREE_SAFE(e_data.depth_placeholder);
	DRW_TEXTURE_FREE_SAFE(e_data.depth_array_placeholder);

	DRW_TEXTURE_FREE_SAFE(e_data.dummy_cache.grid_tx);
	DRW_TEXTURE_FREE_SAFE(e_data.dummy_cache.cube_tx);
	MEM_SAFE_FREE(e_data.dummy_cache.cube_data);
	MEM_SAFE_FREE(e_data.dummy_cache.grid_data);
}
