/*
 * Copyright 2017, Blender Foundation.
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
 * Contributor(s): Blender Institute, Antonio Vazquez
 *
 */

/** \file blender/draw/engines/gpencil/gpencil_depth_of_field.c
 *  \ingroup draw
 *
 * Depth of field post process effect. This code is based on Eevee DOF code
 */

#include "DRW_render.h"

#include "DNA_camera_types.h"
#include "DNA_view3d_types.h"

#include "BKE_camera.h"

#include "DEG_depsgraph_query.h"

#include "GPU_framebuffer.h"
#include "GPU_texture.h"

#include "gpencil_engine.h"

extern char datatoc_gpencil_dof_vert_glsl[];
extern char datatoc_gpencil_dof_frag_glsl[];

/* create dof shaders */
static void gpencil_create_shader_depth_of_field(GPENCIL_e_data *e_data)
{
	e_data->gpencil_dof_downsample_sh = DRW_shader_create(datatoc_gpencil_dof_vert_glsl, NULL,
	                                             datatoc_gpencil_dof_frag_glsl, "#define STEP_DOWNSAMPLE\n");
	e_data->gpencil_dof_scatter_sh = DRW_shader_create(datatoc_gpencil_dof_vert_glsl, NULL,
	                                          datatoc_gpencil_dof_frag_glsl, "#define STEP_SCATTER\n");
	e_data->gpencil_dof_resolve_sh = DRW_shader_create(datatoc_gpencil_dof_vert_glsl, NULL,
	                                          datatoc_gpencil_dof_frag_glsl, "#define STEP_RESOLVE\n");
}

/* init depth of field effect */
int GPENCIL_depth_of_field_init(DrawEngineType *draw_engine_gpencil_type, GPENCIL_e_data *e_data, GPENCIL_Data *vedata, Object *camera)
{
	GPENCIL_StorageList *stl = ((GPENCIL_Data *)vedata)->stl;
	GPENCIL_FramebufferList *fbl = ((GPENCIL_Data *)vedata)->fbl;

	const DRWContextState *draw_ctx = DRW_context_state_get();

	Scene *scene = draw_ctx->scene;
	RegionView3D *rv3d = draw_ctx->rv3d;

	if (!e_data->gpencil_dof_downsample_sh) {
		gpencil_create_shader_depth_of_field(e_data);
	}

	if (camera) {
		const float *viewport_size = DRW_viewport_size_get();
		Camera *cam = (Camera *)camera->data;

		/* Retreive Near and Far distance */
		stl->storage->dof_near_far[0] = -cam->clipsta;
		stl->storage->dof_near_far[1] = -cam->clipend;

		int buffer_size[2] = { (int)viewport_size[0] / 2, (int)viewport_size[1] / 2 };

		/* Setup buffers */
		e_data->gpencil_dof_down_near = DRW_texture_pool_query_2D(buffer_size[0], buffer_size[1], DRW_TEX_RGB_11_11_10,
														   draw_engine_gpencil_type);
		e_data->gpencil_dof_down_far =  DRW_texture_pool_query_2D(buffer_size[0], buffer_size[1], DRW_TEX_RGB_11_11_10,
														   draw_engine_gpencil_type);
		e_data->gpencil_dof_coc =       DRW_texture_pool_query_2D(buffer_size[0], buffer_size[1], DRW_TEX_RG_16,
														   draw_engine_gpencil_type);
		e_data->gpencil_dof_alpha = DRW_texture_pool_query_2D(buffer_size[0], buffer_size[1], DRW_TEX_R_16,
														   draw_engine_gpencil_type);

		GPU_framebuffer_ensure_config(&fbl->dof_down_fb, {
			GPU_ATTACHMENT_NONE,
			GPU_ATTACHMENT_TEXTURE(e_data->gpencil_dof_down_near),
			GPU_ATTACHMENT_TEXTURE(e_data->gpencil_dof_down_far),
			GPU_ATTACHMENT_TEXTURE(e_data->gpencil_dof_coc),
			GPU_ATTACHMENT_TEXTURE(e_data->gpencil_dof_alpha)
		});
			
		/* Go full 32bits for rendering and reduce the color artifacts. */
		DRWTextureFormat fb_format = DRW_state_is_image_render() ? DRW_TEX_RGBA_32 : DRW_TEX_RGBA_16;

		e_data->gpencil_dof_far_blur = DRW_texture_pool_query_2D(buffer_size[0], buffer_size[1], fb_format,
														  draw_engine_gpencil_type);
		GPU_framebuffer_ensure_config(&fbl->dof_scatter_far_fb, {
			GPU_ATTACHMENT_NONE,
			GPU_ATTACHMENT_TEXTURE(e_data->gpencil_dof_far_blur),
		});

		e_data->gpencil_dof_near_blur = DRW_texture_pool_query_2D(buffer_size[0], buffer_size[1], fb_format,
														   draw_engine_gpencil_type);
		GPU_framebuffer_ensure_config(&fbl->dof_scatter_near_fb, {
			GPU_ATTACHMENT_NONE,
			GPU_ATTACHMENT_TEXTURE(e_data->gpencil_dof_near_blur),
			});

		/* Parameters */
		/* TODO UI Options */
		float fstop = cam->gpu_dof.fstop;
		float blades = cam->gpu_dof.num_blades;
		float rotation = cam->gpu_dof.rotation;
		float ratio = 1.0f / cam->gpu_dof.ratio;
		float sensor = BKE_camera_sensor_size(cam->sensor_fit, cam->sensor_x, cam->sensor_y);
		float focus_dist = BKE_camera_object_dof_distance(camera);
		float focal_len = cam->lens;

		UNUSED_VARS(rotation, ratio);

		/* this is factor that converts to the scene scale. focal length and sensor are expressed in mm
		 * unit.scale_length is how many meters per blender unit we have. We want to convert to blender units though
		 * because the shader reads coordinates in world space, which is in blender units.
		 * Note however that focus_distance is already in blender units and shall not be scaled here (see T48157). */
		float scale = (scene->unit.system) ? scene->unit.scale_length : 1.0f;
		float scale_camera = 0.001f / scale;
		/* we want radius here for the aperture number  */
		float aperture = 0.5f * scale_camera * focal_len / fstop;
		float focal_len_scaled = scale_camera * focal_len;
		float sensor_scaled = scale_camera * sensor;

		if (rv3d != NULL) {
			sensor_scaled *= rv3d->viewcamtexcofac[0];
		}

		stl->storage->dof_params[0] = aperture * fabsf(focal_len_scaled / (focus_dist - focal_len_scaled));
		stl->storage->dof_params[1] = -focus_dist;
		stl->storage->dof_params[2] = viewport_size[0] / sensor_scaled;
		stl->storage->dof_bokeh[0] = blades;
		stl->storage->dof_bokeh[1] = rotation;
		stl->storage->dof_bokeh[2] = ratio;
		stl->storage->dof_bokeh[3] = 100.0f; /* TODO: Review this value */

		return 1;
	}

	return 0;
}

/* init cache for depth of field */
void GPENCIL_depth_of_field_cache_init(GPENCIL_e_data *e_data, GPENCIL_Data *vedata)
{
	GPENCIL_PassList *psl = ((GPENCIL_Data *)vedata)->psl;
	GPENCIL_StorageList *stl = ((GPENCIL_Data *)vedata)->stl;

	if (stl->storage->enable_dof == true) {
		/**  Depth of Field algorithm
		 *
		 * Overview :
		 * - Downsample the color buffer into 2 buffers weighted with
		 *   CoC values. Also output CoC into a texture.
		 * - Shoot quads for every pixel and expand it depending on the CoC.
		 *   Do one pass for near Dof and one pass for far Dof.
		 * - Finally composite the 2 blurred buffers with the original render.
		 **/
		DRWShadingGroup *grp;
		struct Gwn_Batch *quad = DRW_cache_fullscreen_quad_get();

		psl->dof_down = DRW_pass_create("DoF Downsample", DRW_STATE_WRITE_COLOR);

		grp = DRW_shgroup_create(e_data->gpencil_dof_downsample_sh, psl->dof_down);
		DRW_shgroup_uniform_texture(grp, "colorBuffer", e_data->input_color_tx);
		DRW_shgroup_uniform_texture(grp, "depthBuffer", e_data->input_depth_tx);
		DRW_shgroup_uniform_vec2(grp, "nearFar", stl->storage->dof_near_far, 1);
		DRW_shgroup_uniform_vec3(grp, "dofParams", stl->storage->dof_params, 1);
		DRW_shgroup_call_add(grp, quad, NULL);

		psl->dof_scatter = DRW_pass_create("DoF Scatter", DRW_STATE_WRITE_COLOR | DRW_STATE_ADDITIVE_FULL);

		/* This create an empty batch of N triangles to be positioned
		 * by the vertex shader 0.4ms against 6ms with instancing */
		const float *viewport_size = DRW_viewport_size_get();
		const int sprite_ct = ((int)viewport_size[0] / 2) * ((int)viewport_size[1] / 2); /* brackets matters */
		grp = DRW_shgroup_empty_tri_batch_create(e_data->gpencil_dof_scatter_sh, psl->dof_scatter, sprite_ct);

		DRW_shgroup_uniform_texture(grp, "colorBuffer", stl->storage->unf_source_buffer);
		DRW_shgroup_uniform_texture(grp, "cocBuffer", e_data->gpencil_dof_coc);
		DRW_shgroup_uniform_vec2(grp, "layerSelection", stl->storage->dof_layer_select, 1);
		DRW_shgroup_uniform_vec4(grp, "bokehParams", stl->storage->dof_bokeh, 1);

		psl->dof_resolve = DRW_pass_create("DoF Resolve", DRW_STATE_WRITE_COLOR);

		grp = DRW_shgroup_create(e_data->gpencil_dof_resolve_sh, psl->dof_resolve);
		DRW_shgroup_uniform_texture(grp, "colorBuffer", e_data->input_color_tx);
		DRW_shgroup_uniform_texture(grp, "nearBuffer", e_data->gpencil_dof_near_blur);
		DRW_shgroup_uniform_texture(grp, "farBuffer", e_data->gpencil_dof_far_blur);
		DRW_shgroup_uniform_texture(grp, "depthBuffer", e_data->input_depth_tx);
		DRW_shgroup_uniform_vec2(grp, "nearFar", stl->storage->dof_near_far, 1);
		DRW_shgroup_uniform_vec3(grp, "dofParams", stl->storage->dof_params, 1);
		DRW_shgroup_call_add(grp, quad, NULL);
	}
}

/* draw depth of field effect */
void GPENCIL_depth_of_field_draw(GPENCIL_e_data *e_data, GPENCIL_Data *vedata)
{
	GPENCIL_PassList *psl = ((GPENCIL_Data *)vedata)->psl;
	GPENCIL_StorageList *stl = ((GPENCIL_Data *)vedata)->stl;
	GPENCIL_FramebufferList *fbl = ((GPENCIL_Data *)vedata)->fbl;
	DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();
	bool is_render = stl->storage->is_render;

	/* Depth Of Field */
	if (stl->storage->enable_dof == true) {
		float clear_col[4] = {0.0f, 0.0f, 0.0f, 0.0f};

		/* Downsample */
		GPU_framebuffer_bind(fbl->dof_down_fb);
		DRW_draw_pass(psl->dof_down);

		/* Scatter Far */
		stl->storage->unf_source_buffer = e_data->gpencil_dof_down_far;
		copy_v2_fl2(stl->storage->dof_layer_select, 0.0f, 1.0f);
		GPU_framebuffer_bind(fbl->dof_scatter_far_fb);
		GPU_framebuffer_clear_color(fbl->dof_scatter_far_fb, clear_col);
		DRW_draw_pass(psl->dof_scatter);

		/* Scatter Near */
		stl->storage->unf_source_buffer = e_data->gpencil_dof_down_near;
		copy_v2_fl2(stl->storage->dof_layer_select, 1.0f, 0.0f);
		GPU_framebuffer_bind(fbl->dof_scatter_near_fb);
		GPU_framebuffer_clear_color(fbl->dof_scatter_near_fb, clear_col);
		DRW_draw_pass(psl->dof_scatter);

		/* Resolve to default pass */
		if ((!is_render) || (fbl->main == NULL)) {
			GPU_framebuffer_bind(dfbl->default_fb);
		}
		else {
			GPU_framebuffer_bind(fbl->main);
		}

		DRW_draw_pass(psl->dof_resolve);
	}
}
