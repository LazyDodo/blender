#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
#include "BLI_math_matrix.h"
#include "lanpr_all.h"
#include "DRW_render.h"
#include "BKE_object.h"
#include "DNA_mesh_types.h"
#include "DNA_camera_types.h"
#include "GPU_immediate.h"
#include "GPU_immediate_util.h"
#include "GPU_framebuffer.h"
#include "DNA_lanpr_types.h"
#include "GPU_draw.h"
#include "DEG_depsgraph_query.h"
#include "RE_pipeline.h"
#include "BLI_rect.h"

#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"
#include "bmesh.h"

#include <math.h>


extern char datatoc_common_fullscreen_vert_glsl[];
extern char datatoc_gpu_shader_3D_normal_smooth_color_vert_glsl[];
extern char datatoc_lanpr_snake_multichannel_fragment[];
extern char datatoc_lanpr_snake_edge_fragment[];
extern char datatoc_lanpr_image_peel_fragment[];
extern char datatoc_lanpr_line_connection_vertex[];
extern char datatoc_lanpr_line_connection_fragment[];
extern char datatoc_lanpr_line_connection_geometry[];
extern char datatoc_lanpr_software_line_width_geometry[];
extern char datatoc_lanpr_atlas_project_passthrough_vertex[];
extern char datatoc_lanpr_atlas_preview_fragment[];
extern char datatoc_lanpr_software_scale_compensate_vertex[];



LANPROneTimeInit OneTime;


static void lanpr_engine_init(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	const DRWContextState *draw_ctx = DRW_context_state_get();
	Scene *scene = DEG_get_evaluated_scene(draw_ctx->depsgraph);
	SceneLANPR *lanpr = &scene->lanpr;
	View3D *v3d = draw_ctx->v3d;
	Object *camera;
	if (v3d) {
		RegionView3D *rv3d = draw_ctx->rv3d;
		camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	}
	else {
		camera = scene->camera;
	}

	if (!OneTime.InitComplete) {
		//lanpr->depth_clamp = 0.01;
		//lanpr->depth_strength = 800;
		//lanpr->normal_clamp = 2;
		//lanpr->normal_strength = 10;
		//lanpr->line_thickness = 2;
		//lanpr->taper_left_distance = 20;
		//lanpr->taper_left_strength = 0.9;
		//lanpr->taper_right_distance = 20;
		//lanpr->taper_right_strength = 0.9;

		//lanpr->line_color[0] = 0.22;
		//lanpr->line_color[1] = 0.29;
		//lanpr->line_color[2] = 0.53;
		//lanpr->line_color[3] = 1;

		//lanpr->background_color[0] = 0.59;
		//lanpr->background_color[1] = 0.90;
		//lanpr->background_color[2] = 0.51;
		//lanpr->background_color[3] = 1;

		lanpr->reloaded = 1;

		OneTime.InitComplete = 1;
	}

	/* SNAKE */

	DRW_texture_ensure_fullscreen_2D_multisample(&txl->depth, GPU_DEPTH_COMPONENT32F, 8, 0);
	DRW_texture_ensure_fullscreen_2D_multisample(&txl->color, GPU_RGBA32F, 8, 0);
	DRW_texture_ensure_fullscreen_2D_multisample(&txl->normal, GPU_RGBA32F, 8, 0);
	DRW_texture_ensure_fullscreen_2D_multisample(&txl->edge_intermediate, GPU_RGBA32F, 8, 0);

	DRW_texture_ensure_fullscreen_2D_multisample(&txl->ms_resolve_depth, GPU_DEPTH_COMPONENT32F, 8, 0);
	DRW_texture_ensure_fullscreen_2D_multisample(&txl->ms_resolve_color, GPU_RGBA32F, 8, 0);


	GPU_framebuffer_ensure_config(&fbl->passes, {
		GPU_ATTACHMENT_TEXTURE(txl->depth),
		GPU_ATTACHMENT_TEXTURE(txl->color),
		GPU_ATTACHMENT_TEXTURE(txl->normal),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});

	GPU_framebuffer_ensure_config(&fbl->edge_intermediate, {
		GPU_ATTACHMENT_TEXTURE(txl->depth),
		GPU_ATTACHMENT_TEXTURE(txl->edge_intermediate),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});

	GPU_framebuffer_ensure_config(&fbl->edge_thinning, {
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_TEXTURE(txl->color),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});


	if (!OneTime.multichannel_shader) {
		OneTime.multichannel_shader =
			GPU_shader_create(
				datatoc_gpu_shader_3D_normal_smooth_color_vert_glsl,
				datatoc_lanpr_snake_multichannel_fragment, NULL, NULL, NULL);

	}
	if (!OneTime.edge_detect_shader) {
		OneTime.edge_detect_shader =
			GPU_shader_create(
				datatoc_common_fullscreen_vert_glsl,
				datatoc_lanpr_snake_edge_fragment, NULL, NULL, NULL);

	}
	if (!OneTime.edge_thinning_shader) {
		OneTime.edge_thinning_shader =
			GPU_shader_create(
				datatoc_common_fullscreen_vert_glsl,
				datatoc_lanpr_image_peel_fragment, NULL, NULL, NULL);

	}
	if (!OneTime.snake_connection_shader) {
		OneTime.snake_connection_shader =
			GPU_shader_create(
				datatoc_lanpr_line_connection_vertex,
				datatoc_lanpr_line_connection_fragment,
				datatoc_lanpr_line_connection_geometry,
				NULL, NULL);
	}

	/* DPIX */
	lanpr_init_atlas_inputs(ved);

	/* SOFTWARE */
	if (!OneTime.software_shader) {
		OneTime.software_shader =
			GPU_shader_create(
				datatoc_lanpr_software_scale_compensate_vertex,
				datatoc_lanpr_atlas_preview_fragment,
				datatoc_lanpr_software_line_width_geometry,
				NULL, NULL);
	}

	GPU_framebuffer_ensure_config(&fbl->software_ms, {
		GPU_ATTACHMENT_TEXTURE(txl->ms_resolve_depth),
		GPU_ATTACHMENT_TEXTURE(txl->ms_resolve_color),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});

}
static void lanpr_engine_free(void){
	void *ved = OneTime.ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;

	DRW_pass_free(psl->color_pass);
	DRW_pass_free(psl->edge_intermediate);

	GPU_framebuffer_free(fbl->passes);
	GPU_framebuffer_free(fbl->edge_intermediate);
	GPU_framebuffer_free(fbl->edge_thinning);
	GPU_framebuffer_free(fbl->software_ms);

	DRW_texture_free(txl->depth);
	DRW_texture_free(txl->color);
	DRW_texture_free(txl->normal);
	DRW_texture_free(txl->edge_intermediate);
	DRW_texture_free(txl->ms_resolve_depth);
	DRW_texture_free(txl->ms_resolve_color);

	BLI_mempool_destroy(stl->g_data->mp_line_strip);
	BLI_mempool_destroy(stl->g_data->mp_line_strip_point);
	BLI_mempool_destroy(stl->g_data->mp_sample);
	BLI_mempool_destroy(stl->g_data->mp_batch_list);

	MEM_freeN(stl->g_data->line_result_8bit);
	MEM_freeN(stl->g_data->line_result);
	MEM_freeN(stl->g_data);

	lanpr_destroy_atlas(vedata);

	stl->g_data = 0;
}

static void lanpr_cache_init(void *vedata){

	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;

	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	if (!stl->g_data) {
		/* Alloc transient pointers */
		stl->g_data = MEM_callocN(sizeof(*stl->g_data), __func__);
		stl->g_data->mp_sample = BLI_mempool_create(sizeof(LANPR_TextureSample), 0, 512, BLI_MEMPOOL_NOP);
		stl->g_data->mp_line_strip = BLI_mempool_create(sizeof(LANPR_LineStrip), 0, 512, BLI_MEMPOOL_NOP);
		stl->g_data->mp_line_strip_point = BLI_mempool_create(sizeof(LANPR_LineStripPoint), 0, 1024, BLI_MEMPOOL_NOP);
		stl->g_data->mp_batch_list = BLI_mempool_create(sizeof(LANPR_BatchItem), 0, 128, BLI_MEMPOOL_NOP);
	}

	LANPR_PrivateData *pd = stl->g_data;

	const DRWContextState *draw_ctx = DRW_context_state_get();
	Scene *scene = DEG_get_evaluated_scene(draw_ctx->depsgraph);
	SceneLANPR *lanpr = &scene->lanpr;
	View3D *v3d = draw_ctx->v3d;
	Object *camera;
	if (v3d) {
		RegionView3D *rv3d = draw_ctx->rv3d;
		camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	}
	else {
		camera = scene->camera;
	}

	psl->color_pass = DRW_pass_create("Color Pass", DRW_STATE_WRITE_COLOR | DRW_STATE_DEPTH_LESS_EQUAL | DRW_STATE_WRITE_DEPTH);
	stl->g_data->multipass_shgrp = DRW_shgroup_create(OneTime.multichannel_shader, psl->color_pass);


	if (lanpr->master_mode == LANPR_MASTER_MODE_SNAKE) {
		struct GPUBatch *quad = DRW_cache_fullscreen_quad_get();

		psl->edge_intermediate = DRW_pass_create("Edge Detection", DRW_STATE_WRITE_COLOR);
		stl->g_data->edge_detect_shgrp = DRW_shgroup_create(OneTime.edge_detect_shader, psl->edge_intermediate);
		DRW_shgroup_uniform_texture_ref(stl->g_data->edge_detect_shgrp, "TexSample0", &txl->depth);
		DRW_shgroup_uniform_texture_ref(stl->g_data->edge_detect_shgrp, "TexSample1", &txl->color);
		DRW_shgroup_uniform_texture_ref(stl->g_data->edge_detect_shgrp, "TexSample2", &txl->normal);

		DRW_shgroup_uniform_float(stl->g_data->edge_detect_shgrp, "zNear", &stl->g_data->znear, 1);
		DRW_shgroup_uniform_float(stl->g_data->edge_detect_shgrp, "zfFar", &stl->g_data->zfar, 1);

		DRW_shgroup_uniform_float(stl->g_data->edge_detect_shgrp, "uValue0", &stl->g_data->normal_clamp, 1);// normal clamp
		DRW_shgroup_uniform_float(stl->g_data->edge_detect_shgrp, "uValue1", &stl->g_data->normal_strength, 1);// normal strength
		DRW_shgroup_uniform_float(stl->g_data->edge_detect_shgrp, "uValue2", &stl->g_data->depth_clamp, 1);// depth clamp
		DRW_shgroup_uniform_float(stl->g_data->edge_detect_shgrp, "uValue3", &stl->g_data->depth_strength, 1);// depth strength
		DRW_shgroup_call_add(stl->g_data->edge_detect_shgrp, quad, NULL);

		psl->edge_thinning = DRW_pass_create("Edge Thinning Stage 1", DRW_STATE_WRITE_COLOR);
		stl->g_data->edge_thinning_shgrp = DRW_shgroup_create(OneTime.edge_thinning_shader, psl->edge_thinning);
		DRW_shgroup_uniform_texture_ref(stl->g_data->edge_thinning_shgrp, "TexSample0", &dtxl->color);
		DRW_shgroup_uniform_int(stl->g_data->edge_thinning_shgrp, "Stage", &stl->g_data->stage, 1);
		DRW_shgroup_call_add(stl->g_data->edge_thinning_shgrp, quad, NULL);

	}
	elif(lanpr->master_mode == LANPR_MASTER_MODE_DPIX && lanpr->active_layer)
	{
		LANPR_LineLayer *ll = lanpr->active_layer;
		psl->dpix_transform_pass = DRW_pass_create("DPIX Transform Stage", DRW_STATE_WRITE_COLOR);
		stl->g_data->dpix_transform_shgrp = DRW_shgroup_create(OneTime.dpix_transform_shader, psl->dpix_transform_pass);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "vert0_tex", &txl->dpix_in_pl);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "vert1_tex", &txl->dpix_in_pr);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "face_normal0_tex", &txl->dpix_in_nl);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "face_normal1_tex", &txl->dpix_in_nr);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "edge_mask_tex", &txl->dpix_in_edge_mask);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "sample_step", &stl->g_data->dpix_sample_step, 1);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "is_perspective", &stl->g_data->dpix_is_perspective, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_transform_shgrp, "viewport", stl->g_data->dpix_viewport, 1);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "buffer_width", &stl->g_data->dpix_buffer_width, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_transform_shgrp, "crease_threshold", &lanpr->crease_threshold, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_transform_shgrp, "crease_fade_threshold", &lanpr->crease_fade_threshold, 1);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "enable_crease", &ll->enable_crease, 1);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "enable_material", &ll->enable_material_seperate, 1);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "enable_edge_mark", &ll->enable_edge_mark, 1);
		DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "enable_intersection", &ll->enable_intersection, 1);

		psl->dpix_preview_pass = DRW_pass_create("DPIX Preview", DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL);
		stl->g_data->dpix_preview_shgrp = DRW_shgroup_create(OneTime.dpix_preview_shader, psl->dpix_preview_pass);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_preview_shgrp, "vert0_tex", &txl->dpix_out_pl);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_preview_shgrp, "vert1_tex", &txl->dpix_out_pr);
		DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_preview_shgrp, "edge_mask_tex", &txl->dpix_in_edge_mask);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "viewport", stl->g_data->dpix_viewport, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "color", ll->color, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "crease_color", ll->crease_color, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "material_color", ll->material_color, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "edge_mark_color", ll->edge_mark_color, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "intersection_color", ll->intersection_color, 1);
		DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "background_color", lanpr->background_color, 1);
		//DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "line_color", ll->line_color, 1); //we have color
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "depth_offset", &stl->g_data->dpix_depth_offset, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "depth_width_influence", &lanpr->depth_width_influence, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "depth_width_curve", &lanpr->depth_width_curve, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "depth_alpha_influence", &lanpr->depth_alpha_influence, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "depth_alpha_curve", &lanpr->depth_alpha_curve, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "line_thickness", &ll->thickness, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "line_thickness_crease", &ll->thickness_crease, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "line_thickness_material", &ll->thickness_material, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "line_thickness_edge_mark", &ll->thickness_edge_mark, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "line_thickness_intersection", &ll->thickness_intersection, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "zNear", &stl->g_data->dpix_znear, 1);
		DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "zFar", &stl->g_data->dpix_zfar, 1);

		pd->begin_index = 0;
		int fsize = sizeof(float) * 4 * TNS_DPIX_TEXTURE_SIZE * TNS_DPIX_TEXTURE_SIZE;

		if (lanpr->reloaded) {
			pd->atlas_pl = MEM_callocN(fsize, "atlas_point_l");
			pd->atlas_pr = MEM_callocN(fsize, "atlas_point_r");
			pd->atlas_nl = MEM_callocN(fsize, "atlas_normal_l");
			pd->atlas_nr = MEM_callocN(fsize, "atlas_normal_l");
			pd->atlas_edge_mask = MEM_callocN(fsize, "atlas_edge_mask"); // should always be float

			pd->dpix_batch_list.first = pd->dpix_batch_list.last = 0;
			BLI_mempool_clear(pd->mp_batch_list);
		}
	} elif(lanpr->master_mode == LANPR_MASTER_MODE_SOFTWARE)
	{
		;
		/*LANPR_LineLayer *ll;
		for (ll = lanpr->line_layers.first; ll; ll = ll->next) {
			ll->shgrp = DRW_shgroup_create(OneTime.software_shader, psl->software_pass);
			DRW_shgroup_uniform_vec4(ll->shgrp, "color", ll->color, 1);
			DRW_shgroup_uniform_vec4(ll->shgrp, "crease_color", ll->crease_color, 1);
			DRW_shgroup_uniform_vec4(ll->shgrp, "material_color", ll->material_color, 1);
			DRW_shgroup_uniform_vec4(ll->shgrp, "edge_mark_color", ll->edge_mark_color, 1);
			DRW_shgroup_uniform_vec4(ll->shgrp, "intersection_color", ll->intersection_color, 1);
			DRW_shgroup_uniform_float(ll->shgrp, "thickness", &ll->thickness, 1);
			DRW_shgroup_uniform_float(ll->shgrp, "thickness_crease", &ll->thickness_crease, 1);
			DRW_shgroup_uniform_float(ll->shgrp, "thickness_material", &ll->thickness_material, 1);
			DRW_shgroup_uniform_float(ll->shgrp, "thickness_edge_mark", &ll->thickness_edge_mark, 1);
			DRW_shgroup_uniform_float(ll->shgrp, "thickness_intersection", &ll->thickness_intersection, 1);
			DRW_shgroup_uniform_vec4(ll->shgrp, "preview_viewport", stl->g_data->dpix_viewport, 1);
			DRW_shgroup_uniform_vec4(ll->shgrp, "output_viewport", stl->g_data->output_viewport, 1);
			if (ll->batch) DRW_shgroup_call_add(ll->shgrp, ll->batch, NULL);
		}*/
	}



}

static void lanpr_cache_populate(void *vedata, Object *ob){

	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData *pd = stl->g_data;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;

	if (!DRW_object_is_renderable(ob)) return;
	if (ob == draw_ctx->object_edit) return;
	if (ob->type != OB_MESH) return;

	struct GPUBatch *geom = DRW_cache_object_surface_get(ob);
	if (geom) {
		DRW_shgroup_call_object_add_no_cull(stl->g_data->multipass_shgrp, geom, ob);
	}

	if (lanpr->master_mode == LANPR_MASTER_MODE_DPIX && lanpr->active_layer) {
		int idx = pd->begin_index;
		if (lanpr->reloaded) {
			pd->begin_index = lanpr_feed_atlas_data_obj(vedata, pd->atlas_pl, pd->atlas_pr, pd->atlas_nl, pd->atlas_nr, pd->atlas_edge_mask, ob, idx);
			lanpr_feed_atlas_trigger_preview_obj(vedata, ob, idx);
		}

	}
}

static void lanpr_cache_finish(void *vedata){
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData *pd = stl->g_data;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;
	float mat[4][4];
	unit_m4(mat);

	LANPR_BatchItem *bi;
	for (bi = pd->dpix_batch_list.first; bi; bi = (void *)bi->Item.next) {
		DRW_shgroup_call_add(pd->dpix_transform_shgrp, bi->dpix_transform_batch, bi->ob->obmat);
		DRW_shgroup_call_add(pd->dpix_preview_shgrp, bi->dpix_preview_batch, 0);
	}

	if (lanpr->reloaded && lanpr->master_mode == LANPR_MASTER_MODE_DPIX && lanpr->active_layer) {
		if (lanpr->render_buffer) {
			lanpr_feed_atlas_data_intersection_cache(vedata, pd->atlas_pl, pd->atlas_pr, pd->atlas_nl, pd->atlas_nr, pd->atlas_edge_mask, pd->begin_index);
			lanpr_create_atlas_intersection_preview(vedata, pd->begin_index);
		}

		GPU_texture_update(txl->dpix_in_pl, GPU_DATA_FLOAT, pd->atlas_pl);
		GPU_texture_update(txl->dpix_in_pr, GPU_DATA_FLOAT, pd->atlas_pr);
		GPU_texture_update(txl->dpix_in_nl, GPU_DATA_FLOAT, pd->atlas_nl);
		GPU_texture_update(txl->dpix_in_nr, GPU_DATA_FLOAT, pd->atlas_nr);
		GPU_texture_update(txl->dpix_in_edge_mask, GPU_DATA_FLOAT, pd->atlas_edge_mask);

		MEM_freeN(pd->atlas_pl);
		MEM_freeN(pd->atlas_pr);
		MEM_freeN(pd->atlas_nl);
		MEM_freeN(pd->atlas_nr);
		MEM_freeN(pd->atlas_edge_mask);
		pd->atlas_pl = 0;
		lanpr->reloaded = 0;

	}

	if (lanpr->render_buffer && lanpr->render_buffer->DPIXIntersectionBatch) {
		DRW_shgroup_call_add(pd->dpix_transform_shgrp, lanpr->render_buffer->DPIXIntersectionTransformBatch, 0);
		DRW_shgroup_call_add(pd->dpix_preview_shgrp, lanpr->render_buffer->DPIXIntersectionBatch, 0);
	}
}

static void lanpr_draw_scene_exec(void *vedata, GPUFrameBuffer *dfb) {
	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_FramebufferList *fbl = ((LANPR_Data *)vedata)->fbl;
	LANPR_PrivateData *pd = stl->g_data;

	float clear_col[4] = { 1.0f, 0.0f, 0.0f, 1.0f };
	float clear_depth = 1.0f;
	uint clear_stencil = 0xFF;

	GPU_framebuffer_bind(fbl->passes);
	GPUFrameBufferBits clear_bits = GPU_DEPTH_BIT | GPU_COLOR_BIT;
	GPU_framebuffer_clear(fbl->passes, clear_bits, clear_col, clear_depth, clear_stencil);

	const DRWContextState *draw_ctx = DRW_context_state_get();
	Scene *scene = DEG_get_evaluated_scene(draw_ctx->depsgraph);
	SceneLANPR *lanpr = &scene->lanpr;
	View3D *v3d = draw_ctx->v3d;
	Object *camera;
	if (v3d) {
		RegionView3D *rv3d = draw_ctx->rv3d;
		camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	}
	else {
		camera = scene->camera;
	}

	if (lanpr->master_mode == LANPR_MASTER_MODE_DPIX) {
		DRW_draw_pass(psl->color_pass);
		lanpr_dpix_draw_scene(txl, fbl, psl, stl->g_data, lanpr, dfb);
	}
	elif(lanpr->master_mode == LANPR_MASTER_MODE_SNAKE)
	{
		DRW_draw_pass(psl->color_pass);
		lanpr_snake_draw_scene(txl, fbl, psl, stl->g_data, lanpr, dfb);
	}
	elif(lanpr->master_mode == LANPR_MASTER_MODE_SOFTWARE)
	{
		GPU_framebuffer_bind(fbl->software_ms);
		GPU_framebuffer_clear(fbl->software_ms, clear_bits, lanpr->background_color, clear_depth, clear_stencil);

		if (lanpr->render_buffer && lanpr->render_buffer->ChainDrawBatch) {
			LANPR_LineLayer* ll;
			for (ll = lanpr->line_layers.first; ll; ll = ll->next) {
				LANPR_RenderBuffer* rb;
				psl->software_pass = DRW_pass_create("Software Render Preview", DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL);
				rb = lanpr->render_buffer;
				rb->ChainShgrp = DRW_shgroup_create(OneTime.software_shader, psl->software_pass);
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "color", ll->color, 1);
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "crease_color", ll->crease_color, 1);
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "material_color", ll->material_color, 1);
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "edge_mark_color", ll->edge_mark_color, 1);
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "intersection_color", ll->intersection_color, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "thickness", &ll->thickness, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "thickness_crease", &ll->thickness_crease, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "thickness_material", &ll->thickness_material, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "thickness_edge_mark", &ll->thickness_edge_mark, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "thickness_intersection", &ll->thickness_intersection, 1);

				DRW_shgroup_uniform_int(rb->ChainShgrp, "occlusion_level_begin", &ll->qi_begin, 1);
				DRW_shgroup_uniform_int(rb->ChainShgrp, "occlusion_level_end", &ll->qi_end, 1);

				int texw = GPU_texture_width(txl->ms_resolve_color), texh = GPU_texture_height(txl->ms_resolve_color);;
				pd->output_viewport[2] = scene->r.xsch;
				pd->output_viewport[3] = scene->r.ysch;
				pd->dpix_viewport[2] = texw;
				pd->dpix_viewport[3] = texh;
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "preview_viewport", stl->g_data->dpix_viewport, 1);
				DRW_shgroup_uniform_vec4(rb->ChainShgrp, "output_viewport", stl->g_data->output_viewport, 1);

				float *tld = &lanpr->taper_left_distance, *tls = &lanpr->taper_left_strength,
					  *trd = &lanpr->taper_right_distance, *trs = &lanpr->taper_right_strength;

				DRW_shgroup_uniform_float(rb->ChainShgrp, "TaperLDist", tld, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "TaperLStrength", tls, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "TaperRDist", lanpr->use_same_taper ? tld : trd, 1);
				DRW_shgroup_uniform_float(rb->ChainShgrp, "TaperRStrength", lanpr->use_same_taper ? tls : trs, 1);

				//need to add component enable/disable option.
				DRW_shgroup_call_add(rb->ChainShgrp, lanpr->render_buffer->ChainDrawBatch, NULL);
				// debug purpose
				//DRW_draw_pass(psl->color_pass);
				//DRW_draw_pass(psl->color_pass);
				DRW_draw_pass(psl->software_pass);
			}
		}

		GPU_framebuffer_blit(fbl->software_ms, 0, dfb, 0, GPU_COLOR_BIT);

	}
}

static void lanpr_draw_scene(void *vedata) {
	DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();
	lanpr_draw_scene_exec(vedata, dfbl->default_fb);
}

void LANPR_render_cache(
	void *vedata, struct Object *ob,
	struct RenderEngine *engine, struct Depsgraph *UNUSED(depsgraph)){

	lanpr_cache_populate(vedata, ob);

}

static void lanpr_render_to_image(LANPR_Data *vedata, RenderEngine *engine, struct RenderLayer *render_layer, const rcti *rect){
	LANPR_StorageList *stl = vedata->stl;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	Scene *scene = DEG_get_evaluated_scene(draw_ctx->depsgraph);

	/* refered to eevee's code */

	/* Init default FB and render targets:
	 * In render mode the default framebuffer is not generated
	 * because there is no viewport. So we need to manually create it or
	 * not use it. For code clarity we just allocate it make use of it. */
	DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	DRW_texture_ensure_fullscreen_2D(&dtxl->depth, GPU_DEPTH_COMPONENT32F, 0);
	DRW_texture_ensure_fullscreen_2D(&dtxl->color, GPU_RGBA32F, 0);

	GPU_framebuffer_ensure_config(&dfbl->default_fb, {
		GPU_ATTACHMENT_TEXTURE(dtxl->depth),
		GPU_ATTACHMENT_TEXTURE(dtxl->color),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});
	scene->lanpr.reloaded = 1;

	lanpr_engine_init(vedata);
	lanpr_cache_init(vedata);
	DRW_render_object_iter(vedata, engine, draw_ctx->depsgraph, LANPR_render_cache);
	lanpr_cache_finish(vedata);

	float clear_col[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	float clear_depth = 1.0f;
	uint clear_stencil = 0xFF;
	GPUFrameBufferBits clear_bits = GPU_DEPTH_BIT | GPU_COLOR_BIT;

	GPU_framebuffer_bind(dfbl->default_fb);
	GPU_framebuffer_clear(dfbl->default_fb, clear_bits, clear_col, clear_depth, clear_stencil);

	lanpr_draw_scene_exec(vedata, dfbl->default_fb);

	// read it back so we can again display and save it.
	const char *viewname = RE_GetActiveRenderView(engine->re);
	RenderPass *rp = RE_pass_find_by_name(render_layer, RE_PASSNAME_COMBINED, viewname);
	GPU_framebuffer_bind(dfbl->default_fb);
	GPU_framebuffer_read_color(dfbl->default_fb,
	                           rect->xmin, rect->ymin,
	                           BLI_rcti_size_x(rect), BLI_rcti_size_y(rect),
	                           4, 0, rp->rect);

}

static void lanpr_view_update(void *vedata){
	//LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	//if (stl->g_data) {
	//	stl->g_data->view_updated = true;
	//}

	//our update flag is in SceneLANPR.
	const DRWContextState *draw_ctx = DRW_context_state_get();
	SceneLANPR *lanpr = &DEG_get_evaluated_scene(draw_ctx->depsgraph)->lanpr;
	lanpr->reloaded = 1; // very bad solution, this will slow down animation.
}

//static void lanpr_id_update(void *vedata, ID *id){
//	const DRWContextState *draw_ctx = DRW_context_state_get();
//    SceneLANPR *lanpr = &DEG_get_evaluated_scene(draw_ctx->depsgraph)->lanpr;
//
//	/* look at eevee_engine.c */
//	switch (GS(id->name)) {
//		case ID_OB:
//		    //seems doesn't need this one currently...
//			//eevee_id_object_update(vedata, (Object *)id);
//			lanpr->reloaded = 1;
//			break;
//		case ID_ME:
//		    lanpr->reloaded=1;
//			break;
//		default:
//			/* pass */
//			break;
//	}
//}


static const DrawEngineDataSize lanpr_data_size = DRW_VIEWPORT_DATA_SIZE(LANPR_Data);

DrawEngineType draw_engine_lanpr_type = {
	NULL, NULL,
	N_("LANPR"),
	&lanpr_data_size, // why should we have the "&" ?
	&lanpr_engine_init,
	&lanpr_engine_free,
	&lanpr_cache_init,
	&lanpr_cache_populate,
	&lanpr_cache_finish,
	NULL,//draw background
	&lanpr_draw_scene,//draw scene, looks like that not much difference except a camera overlay image.
	&lanpr_view_update,
	NULL, //&lanpr_id_update, wait till I figure out how to do this.
	&lanpr_render_to_image,
};

RenderEngineType DRW_engine_viewport_lanpr_type = {
	NULL, NULL,
	LANPR_ENGINE, N_("LANPR"), RE_INTERNAL,
	NULL,// update
	&DRW_render_to_image,// render to img
	NULL,// bake
	NULL,// doesn't seem to be what I thought it was... &lanpr_view_update,// view update
	NULL,// render to view
	NULL,// update in script
	NULL,// update in render pass
	&draw_engine_lanpr_type,
	{NULL, NULL, NULL}
};

