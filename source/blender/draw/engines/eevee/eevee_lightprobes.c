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
static struct GPUTexture *create_hammersley_sample_texture(int samples)
{
	struct GPUTexture *tex;
	float (*texels)[2] = MEM_mallocN(sizeof(float[2]) * samples, "hammersley_tex");
	int i;

	for (i = 0; i < samples; i++) {
		double dphi;
		BLI_hammersley_1D(i, &dphi);
		float phi = (float)dphi * 2.0f * M_PI;
		texels[i][0] = cosf(phi);
		texels[i][1] = sinf(phi);
	}

	tex = DRW_texture_create_1D(samples, GPU_RG16F, DRW_TEX_WRAP, (float *)texels);
	MEM_freeN(texels);
	return tex;
}

static void lightprobe_shaders_init(void)
{
	const char *filter_defines = "#define HAMMERSLEY_SIZE " STRINGIFY(HAMMERSLEY_SIZE) "\n"
#if defined(IRRADIANCE_SH_L2)
	                             "#define IRRADIANCE_SH_L2\n"
#elif defined(IRRADIANCE_CUBEMAP)
	                             "#define IRRADIANCE_CUBEMAP\n"
#elif defined(IRRADIANCE_HL2)
	                             "#define IRRADIANCE_HL2\n"
#endif
	                             "#define NOISE_SIZE 64\n";

	char *shader_str = NULL;
	char *vert_str = NULL;

	shader_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_common_uniforms_lib_glsl,
	        datatoc_bsdf_common_lib_glsl,
	        datatoc_bsdf_sampling_lib_glsl,
	        datatoc_lightprobe_filter_glossy_frag_glsl);

	e_data.probe_filter_glossy_sh = DRW_shader_create(
	        datatoc_lightprobe_vert_glsl, datatoc_lightprobe_geom_glsl, shader_str, filter_defines);

	e_data.probe_default_sh = DRW_shader_create(
	        datatoc_background_vert_glsl, NULL, datatoc_default_world_frag_glsl, NULL);

	e_data.probe_default_studiolight_sh = DRW_shader_create(
	        datatoc_background_vert_glsl, NULL, datatoc_default_world_frag_glsl, "#define LOOKDEV\n");

	MEM_freeN(shader_str);

	shader_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_common_uniforms_lib_glsl,
	        datatoc_bsdf_common_lib_glsl,
	        datatoc_bsdf_sampling_lib_glsl,
	        datatoc_lightprobe_filter_diffuse_frag_glsl);

	e_data.probe_filter_diffuse_sh = DRW_shader_create_fullscreen(shader_str, filter_defines);

	MEM_freeN(shader_str);

	shader_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_common_uniforms_lib_glsl,
	        datatoc_bsdf_common_lib_glsl,
	        datatoc_bsdf_sampling_lib_glsl,
	        datatoc_lightprobe_filter_visibility_frag_glsl);

	e_data.probe_filter_visibility_sh = DRW_shader_create_fullscreen(shader_str, filter_defines);

	MEM_freeN(shader_str);

	shader_str = BLI_string_joinN(
	        datatoc_octahedron_lib_glsl,
	        datatoc_common_view_lib_glsl,
	        datatoc_common_uniforms_lib_glsl,
	        datatoc_bsdf_common_lib_glsl,
	        datatoc_irradiance_lib_glsl,
	        datatoc_lightprobe_lib_glsl,
	        datatoc_lightprobe_grid_display_frag_glsl);

	vert_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_lightprobe_grid_display_vert_glsl);

	e_data.probe_grid_display_sh = DRW_shader_create(vert_str, NULL, shader_str, filter_defines);

	MEM_freeN(vert_str);
	MEM_freeN(shader_str);

	e_data.probe_grid_fill_sh = DRW_shader_create_fullscreen(
	        datatoc_lightprobe_grid_fill_frag_glsl, filter_defines);

	shader_str = BLI_string_joinN(
	        datatoc_octahedron_lib_glsl,
	        datatoc_common_view_lib_glsl,
	        datatoc_common_uniforms_lib_glsl,
	        datatoc_bsdf_common_lib_glsl,
	        datatoc_lightprobe_lib_glsl,
	        datatoc_lightprobe_cube_display_frag_glsl);

	vert_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_lightprobe_cube_display_vert_glsl);

	e_data.probe_cube_display_sh = DRW_shader_create(vert_str, NULL, shader_str, NULL);

	MEM_freeN(vert_str);
	MEM_freeN(shader_str);

	vert_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_lightprobe_planar_display_vert_glsl);

	shader_str = BLI_string_joinN(
	        datatoc_common_view_lib_glsl,
	        datatoc_lightprobe_planar_display_frag_glsl);

	e_data.probe_planar_display_sh = DRW_shader_create(vert_str, NULL, shader_str, NULL);

	MEM_freeN(vert_str);
	MEM_freeN(shader_str);

	e_data.probe_planar_downsample_sh = DRW_shader_create(
	        datatoc_lightprobe_planar_downsample_vert_glsl,
	        datatoc_lightprobe_planar_downsample_geom_glsl,
	        datatoc_lightprobe_planar_downsample_frag_glsl,
	        NULL);

	e_data.hammersley = create_hammersley_sample_texture(HAMMERSLEY_SIZE);
}

void EEVEE_lightprobes_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
	EEVEE_CommonUniformBuffer *common_data = &sldata->common_data;
	EEVEE_TextureList *txl = vedata->txl;
	EEVEE_StorageList *stl = vedata->stl;

	const DRWContextState *draw_ctx = DRW_context_state_get();
	const Scene *scene_eval = DEG_get_evaluated_scene(draw_ctx->depsgraph);

	if (!e_data.probe_filter_glossy_sh) {
		lightprobe_shaders_init();
	}

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

/* Only init the passes usefull for rendering the light cache. */
void EEVEE_lightbake_cache_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata, GPUTexture *rt_color, GPUTexture *rt_depth)
{
	EEVEE_PassList *psl = vedata->psl;
	EEVEE_LightProbesInfo *pinfo = sldata->probes;

	{
		psl->probe_glossy_compute = DRW_pass_create("LightProbe Glossy Compute", DRW_STATE_WRITE_COLOR);

		DRWShadingGroup *grp = DRW_shgroup_create(e_data.probe_filter_glossy_sh, psl->probe_glossy_compute);
		DRW_shgroup_uniform_float(grp, "intensityFac", &pinfo->intensity_fac, 1);
		DRW_shgroup_uniform_float(grp, "sampleCount", &pinfo->samples_ct, 1);
		DRW_shgroup_uniform_float(grp, "invSampleCount", &pinfo->invsamples_ct, 1);
		DRW_shgroup_uniform_float(grp, "roughnessSquared", &pinfo->roughness, 1);
		DRW_shgroup_uniform_float(grp, "lodFactor", &pinfo->lodfactor, 1);
		DRW_shgroup_uniform_float(grp, "lodMax", &pinfo->lod_rt_max, 1);
		DRW_shgroup_uniform_float(grp, "texelSize", &pinfo->texel_size, 1);
		DRW_shgroup_uniform_float(grp, "paddingSize", &pinfo->padding_size, 1);
		DRW_shgroup_uniform_int(grp, "Layer", &pinfo->layer, 1);
		DRW_shgroup_uniform_texture(grp, "texHammersley", e_data.hammersley);
		// DRW_shgroup_uniform_texture(grp, "texJitter", e_data.jitter);
		DRW_shgroup_uniform_texture(grp, "probeHdr", rt_color);
		DRW_shgroup_call_add(grp, DRW_cache_fullscreen_quad_get(), NULL);
	}

	{
		psl->probe_diffuse_compute = DRW_pass_create("LightProbe Diffuse Compute", DRW_STATE_WRITE_COLOR);

		DRWShadingGroup *grp = DRW_shgroup_create(e_data.probe_filter_diffuse_sh, psl->probe_diffuse_compute);
#ifdef IRRADIANCE_SH_L2
		DRW_shgroup_uniform_int(grp, "probeSize", &pinfo->shres, 1);
#else
		DRW_shgroup_uniform_float(grp, "sampleCount", &pinfo->samples_ct, 1);
		DRW_shgroup_uniform_float(grp, "invSampleCount", &pinfo->invsamples_ct, 1);
		DRW_shgroup_uniform_float(grp, "lodFactor", &pinfo->lodfactor, 1);
		DRW_shgroup_uniform_float(grp, "lodMax", &pinfo->lod_rt_max, 1);
		DRW_shgroup_uniform_texture(grp, "texHammersley", e_data.hammersley);
#endif
		DRW_shgroup_uniform_float(grp, "intensityFac", &pinfo->intensity_fac, 1);
		DRW_shgroup_uniform_texture(grp, "probeHdr", rt_color);

		struct Gwn_Batch *geom = DRW_cache_fullscreen_quad_get();
		DRW_shgroup_call_add(grp, geom, NULL);
	}

	{
		psl->probe_visibility_compute = DRW_pass_create("LightProbe Visibility Compute", DRW_STATE_WRITE_COLOR);

		DRWShadingGroup *grp = DRW_shgroup_create(e_data.probe_filter_visibility_sh, psl->probe_visibility_compute);
		DRW_shgroup_uniform_int(grp, "outputSize", &pinfo->shres, 1);
		DRW_shgroup_uniform_float(grp, "visibilityRange", &pinfo->visibility_range, 1);
		DRW_shgroup_uniform_float(grp, "visibilityBlur", &pinfo->visibility_blur, 1);
		DRW_shgroup_uniform_float(grp, "sampleCount", &pinfo->samples_ct, 1);
		DRW_shgroup_uniform_float(grp, "invSampleCount", &pinfo->invsamples_ct, 1);
		DRW_shgroup_uniform_float(grp, "storedTexelSize", &pinfo->texel_size, 1);
		DRW_shgroup_uniform_float(grp, "nearClip", &pinfo->near_clip, 1);
		DRW_shgroup_uniform_float(grp, "farClip", &pinfo->far_clip, 1);
		DRW_shgroup_uniform_texture(grp, "texHammersley", e_data.hammersley);
		DRW_shgroup_uniform_texture(grp, "probeDepth", rt_depth);

		struct Gwn_Batch *geom = DRW_cache_fullscreen_quad_get();
		DRW_shgroup_call_add(grp, geom, NULL);
	}
}

void EEVEE_lightprobes_cache_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
	EEVEE_TextureList *txl = vedata->txl;
	EEVEE_PassList *psl = vedata->psl;
	EEVEE_StorageList *stl = vedata->stl;
	EEVEE_LightProbesInfo *pinfo = sldata->probes;
	EEVEE_LightCache *lcache = stl->g_data->light_cache;

	{
		psl->probe_background = DRW_pass_create("World Probe Background Pass", DRW_STATE_WRITE_COLOR | DRW_STATE_DEPTH_EQUAL);

		struct Gwn_Batch *geom = DRW_cache_fullscreen_quad_get();
		DRWShadingGroup *grp = NULL;

		const DRWContextState *draw_ctx = DRW_context_state_get();
		Scene *scene = draw_ctx->scene;
		World *wo = scene->world;

		float *col = ts.colorBackground;

		/* LookDev */
		EEVEE_lookdev_cache_init(vedata, &grp, e_data.probe_default_studiolight_sh, psl->probe_background, wo, pinfo);
		/* END */
		if (!grp && wo) {
			col = &wo->horr;
			bool wo_sh_compiled = true;

			if (wo->use_nodes && wo->nodetree) {
				static float error_col[3] = {1.0f, 0.0f, 1.0f};
				static float compile_col[3] = {0.5f, 0.5f, 0.5f};
				struct GPUMaterial *gpumat = EEVEE_material_world_lightprobe_get(scene, wo);

				GPUMaterialStatus status = GPU_material_status(gpumat);

				switch (status) {
					case GPU_MAT_SUCCESS:
						grp = DRW_shgroup_material_create(gpumat, psl->probe_background);
						DRW_shgroup_uniform_float(grp, "backgroundAlpha", &stl->g_data->background_alpha, 1);
						DRW_shgroup_call_add(grp, geom, NULL);
						wo_sh_compiled = true;
						break;
					case GPU_MAT_QUEUED:
						pinfo->all_materials_updated = false;
						wo_sh_compiled = false;
						/* TODO Bypass probe compilation. */
						col = compile_col;
						break;
					case GPU_MAT_FAILED:
					default:
						wo_sh_compiled = true;
						col = error_col;
						break;
				}
			}

			if (wo->update_flag != 0 || pinfo->prev_world != wo || pinfo->prev_wo_sh_compiled != wo_sh_compiled) {
				pinfo->update_world |= PROBE_UPDATE_ALL;
				pinfo->studiolight_index = 0;
				pinfo->prev_wo_sh_compiled = wo_sh_compiled;
				pinfo->prev_world = wo;
			}
			wo->update_flag = 0;
		}
		else if (pinfo->prev_world) {
			pinfo->update_world |= PROBE_UPDATE_ALL;
			pinfo->studiolight_index = 0;
			pinfo->prev_wo_sh_compiled = false;
			pinfo->prev_world = NULL;
		}

		/* Fallback if shader fails or if not using nodetree. */
		if (grp == NULL) {
			grp = DRW_shgroup_create(e_data.probe_default_sh, psl->probe_background);
			DRW_shgroup_uniform_vec3(grp, "color", col, 1);
			DRW_shgroup_uniform_float(grp, "backgroundAlpha", &stl->g_data->background_alpha, 1);
			DRW_shgroup_call_add(grp, geom, NULL);
		}
	}

	{
		DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL | DRW_STATE_CULL_BACK;
		psl->probe_display = DRW_pass_create("LightProbe Display", state);

		DRW_shgroup_instance_format(e_data.format_probe_display_cube, {
		    {"probe_id",       DRW_ATTRIB_INT, 1},
		    {"probe_location", DRW_ATTRIB_FLOAT, 3},
		    {"sphere_size",    DRW_ATTRIB_FLOAT, 1},
		});

		DRWShadingGroup *grp = DRW_shgroup_instance_create(
		        e_data.probe_cube_display_sh,
		        psl->probe_display,
		        DRW_cache_sphere_get(),
		        e_data.format_probe_display_cube);
		stl->g_data->cube_display_shgrp = grp;
		DRW_shgroup_uniform_texture_ref(grp, "probeCubes", &lcache->cube_tx);
		DRW_shgroup_uniform_block(grp, "common_block", sldata->common_ubo);

		DRW_shgroup_instance_format(e_data.format_probe_display_planar, {
		    {"probe_id", DRW_ATTRIB_INT, 1},
		    {"probe_mat", DRW_ATTRIB_FLOAT, 16},
		});

		grp = DRW_shgroup_instance_create(
		        e_data.probe_planar_display_sh,
		        psl->probe_display,
		        DRW_cache_quad_get(),
		        e_data.format_probe_display_planar);
		stl->g_data->planar_display_shgrp = grp;
		DRW_shgroup_uniform_texture_ref(grp, "probePlanars", &txl->planar_pool);
	}

	{
		psl->probe_planar_downsample_ps = DRW_pass_create("LightProbe Planar Downsample", DRW_STATE_WRITE_COLOR);

		DRWShadingGroup *grp = DRW_shgroup_create(e_data.probe_planar_downsample_sh, psl->probe_planar_downsample_ps);
		DRW_shgroup_uniform_texture_ref(grp, "source", &txl->planar_pool);
		DRW_shgroup_uniform_float(grp, "fireflyFactor", &sldata->common_data.ssr_firefly_fac, 1);
		DRW_shgroup_call_instances_add(grp, DRW_cache_fullscreen_quad_get(), NULL, (uint *)&pinfo->num_planar);
	}
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

	/* For shading, save max level of the octahedron map */
	int mipsize = GPU_texture_width(stl->g_data->light_cache->cube_tx);
	sldata->common_data.prb_lod_cube_max = (float)(floorf(log2f(mipsize)) - MIN_CUBE_LOD_LEVEL) - 1.0f;
}

void EEVEE_lightbake_render_world(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata, struct GPUFrameBuffer *face_fb[6])
{
	EEVEE_PassList *psl = vedata->psl;
	DRWMatrixState matstate;
	float (*viewmat)[4] = matstate.mat[DRW_MAT_VIEW];
	float (*viewinv)[4] = matstate.mat[DRW_MAT_VIEWINV];
	float (*persmat)[4] = matstate.mat[DRW_MAT_PERS];
	float (*persinv)[4] = matstate.mat[DRW_MAT_PERSINV];
	float (*winmat)[4] = matstate.mat[DRW_MAT_WIN];
	float (*wininv)[4] = matstate.mat[DRW_MAT_WININV];

	perspective_m4(winmat, -0.1f, 0.1f, -0.1f, 0.1f, 0.1f, 1.0f);
	invert_m4_m4(wininv, winmat);

	for (int i = 0; i < 6; ++i) {
		/* Setup custom matrices */
		copy_m4_m4(viewmat, cubefacemat[i]);
		mul_m4_m4m4(persmat, winmat, viewmat);
		invert_m4_m4(persinv, persmat);
		invert_m4_m4(viewinv, viewmat);
		DRW_viewport_matrix_override_set_all(&matstate);

		GPU_framebuffer_bind(face_fb[i]);
		/* For world probe, we don't need to clear since we render the background directly. */
		GPU_framebuffer_clear_depth(face_fb[i], 1.0f);
		DRW_draw_pass(psl->probe_background);
	}
}

/* Glossy filter rt_color to light_cache->cube_tx at index probe_idx */
void EEVEE_lightbake_filter_glossy(
        EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata,
        struct GPUTexture *rt_color, struct GPUFrameBuffer *fb,
        int probe_idx, float intensity)
{
	EEVEE_PassList *psl = vedata->psl;
	EEVEE_LightProbesInfo *pinfo = sldata->probes;
	EEVEE_LightCache *light_cache = vedata->stl->g_data->light_cache;

	float target_size = (float)GPU_texture_width(rt_color);

	/* Max lod used from the render target probe */
	pinfo->lod_rt_max = floorf(log2f(target_size)) - 2.0f;
	pinfo->intensity_fac = intensity;

	/* Start fresh */
	GPU_framebuffer_ensure_config(&fb, {
		GPU_ATTACHMENT_NONE,
		GPU_ATTACHMENT_NONE
	});

	/* 2 - Let gpu create Mipmaps for Filtered Importance Sampling. */
	/* Bind next framebuffer to be able to gen. mips for probe_rt. */
	EEVEE_downsample_cube_buffer(vedata, rt_color, (int)(pinfo->lod_rt_max));

	/* 3 - Render to probe array to the specified layer, do prefiltering. */
	int mipsize = GPU_texture_width(light_cache->cube_tx);
	const int maxlevel = (int)floorf(log2f(mipsize));
	const int min_lod_level = MIN_CUBE_LOD_LEVEL;
	for (int i = 0; i < maxlevel - min_lod_level; i++) {
		float bias = (i == 0) ? -1.0f : 1.0f;
		pinfo->texel_size = 1.0f / (float)mipsize;
		pinfo->padding_size = (float)(1 << (maxlevel - min_lod_level - 1 - i));
		/* XXX : WHY THE HECK DO WE NEED THIS ??? */
		/* padding is incorrect without this! float precision issue? */
		if (pinfo->padding_size > 32) {
			pinfo->padding_size += 5;
		}
		if (pinfo->padding_size > 16) {
			pinfo->padding_size += 4;
		}
		else if (pinfo->padding_size > 8) {
			pinfo->padding_size += 2;
		}
		else if (pinfo->padding_size > 4) {
			pinfo->padding_size += 1;
		}
		pinfo->layer = probe_idx;
		pinfo->roughness = (float)i / ((float)maxlevel - 4.0f);
		pinfo->roughness *= pinfo->roughness; /* Disney Roughness */
		pinfo->roughness *= pinfo->roughness; /* Distribute Roughness accros lod more evenly */
		CLAMP(pinfo->roughness, 1e-8f, 0.99999f); /* Avoid artifacts */

#if 1 /* Variable Sample count (fast) */
		switch (i) {
			case 0: pinfo->samples_ct = 1.0f; break;
			case 1: pinfo->samples_ct = 16.0f; break;
			case 2: pinfo->samples_ct = 32.0f; break;
			case 3: pinfo->samples_ct = 64.0f; break;
			default: pinfo->samples_ct = 128.0f; break;
		}
#else /* Constant Sample count (slow) */
		pinfo->samples_ct = 1024.0f;
#endif

		pinfo->invsamples_ct = 1.0f / pinfo->samples_ct;
		pinfo->lodfactor = bias + 0.5f * log((float)(target_size * target_size) * pinfo->invsamples_ct) / log(2);

		GPU_framebuffer_ensure_config(&fb, {
			GPU_ATTACHMENT_NONE,
			GPU_ATTACHMENT_TEXTURE_MIP(light_cache->cube_tx, i)
		});
		GPU_framebuffer_bind(fb);
		GPU_framebuffer_viewport_set(fb, 0, 0, mipsize, mipsize);
		DRW_draw_pass(psl->probe_glossy_compute);

		mipsize /= 2;
		CLAMP_MIN(mipsize, 1);
	}
}

/* Diffuse filter rt_color to light_cache->grid_tx at index grid_offset */
void EEVEE_lightbake_filter_diffuse(
        EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata,
        struct GPUTexture *rt_color, struct GPUFrameBuffer *fb,
        int grid_offset, float intensity)
{
	EEVEE_PassList *psl = vedata->psl;
	EEVEE_LightProbesInfo *pinfo = sldata->probes;
	EEVEE_LightCache *light_cache = vedata->stl->g_data->light_cache;

	float target_size = (float)GPU_texture_width(rt_color);

	pinfo->intensity_fac = intensity;

	/* find cell position on the virtual 3D texture */
	/* NOTE : Keep in sync with load_irradiance_cell() */
#if defined(IRRADIANCE_SH_L2)
	int size[2] = {3, 3};
#elif defined(IRRADIANCE_CUBEMAP)
	int size[2] = {8, 8};
	pinfo->samples_ct = 1024.0f;
#elif defined(IRRADIANCE_HL2)
	int size[2] = {3, 2};
	pinfo->samples_ct = 1024.0f;
#endif

	int cell_per_row = GPU_texture_width(light_cache->grid_tx) / size[0];
	int x = size[0] * (grid_offset % cell_per_row);
	int y = size[1] * (grid_offset / cell_per_row);

#ifndef IRRADIANCE_SH_L2
	/* Tweaking parameters to balance perf. vs precision */
	const float bias = 0.0f;
	pinfo->invsamples_ct = 1.0f / pinfo->samples_ct;
	pinfo->lodfactor = bias + 0.5f * log((float)(target_size * target_size) * pinfo->invsamples_ct) / log(2);
	pinfo->lod_rt_max = floorf(log2f(target_size)) - 2.0f;
#else
	pinfo->shres = 32; /* Less texture fetches & reduce branches */
	pinfo->lod_rt_max = 2.0f; /* Improve cache reuse */
#endif

	/* Start fresh */
	GPU_framebuffer_ensure_config(&fb, {
		GPU_ATTACHMENT_NONE,
		GPU_ATTACHMENT_NONE
	});

	/* 4 - Compute diffuse irradiance */
	EEVEE_downsample_cube_buffer(vedata, rt_color, (int)(pinfo->lod_rt_max));

	GPU_framebuffer_ensure_config(&fb, {
		GPU_ATTACHMENT_NONE,
		GPU_ATTACHMENT_TEXTURE_LAYER(light_cache->grid_tx, 0)
	});
	GPU_framebuffer_bind(fb);
	GPU_framebuffer_viewport_set(fb, x, y, size[0], size[1]);
	DRW_draw_pass(psl->probe_diffuse_compute);
}

/* Filter rt_depth to light_cache->grid_tx at index grid_offset */
void EEVEE_lightbake_filter_visibility(
        EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata,
        struct GPUTexture *UNUSED(rt_depth), struct GPUFrameBuffer *fb,
        float clipsta, float clipend,
        float vis_range, float vis_blur, int vis_size,
        int grid_offset)
{
	EEVEE_PassList *psl = vedata->psl;
	EEVEE_LightProbesInfo *pinfo = sldata->probes;
	EEVEE_LightCache *light_cache = vedata->stl->g_data->light_cache;

	pinfo->samples_ct = 512.0f; /* TODO refine */
	pinfo->invsamples_ct = 1.0f / pinfo->samples_ct;
	pinfo->shres = vis_size;
	pinfo->visibility_range = vis_range;
	pinfo->visibility_blur = vis_blur;
	pinfo->near_clip = -clipsta;
	pinfo->far_clip = -clipend;
	pinfo->texel_size = 1.0f / (float)vis_size;

	int cell_per_col = GPU_texture_height(light_cache->grid_tx) / vis_size;
	int cell_per_row = GPU_texture_width(light_cache->grid_tx) / vis_size;
	int x = vis_size * (grid_offset % cell_per_row);
	int y = vis_size * ((grid_offset / cell_per_row) % cell_per_col);
	int layer = 1 + ((grid_offset / cell_per_row) / cell_per_col);

	GPU_framebuffer_ensure_config(&fb, {
		GPU_ATTACHMENT_NONE,
		GPU_ATTACHMENT_TEXTURE_LAYER(light_cache->grid_tx, layer)
	});
	GPU_framebuffer_bind(fb);
	GPU_framebuffer_viewport_set(fb, x, y, vis_size, vis_size);
	DRW_draw_pass(psl->probe_visibility_compute);
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
