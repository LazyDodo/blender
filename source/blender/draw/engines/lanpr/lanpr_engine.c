#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
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



LANPROneTimeInit OneTime;


static void lanpr_engine_init(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	//LANPR_ViewLayerData *sldata = EEVEE_view_layer_data_ensure();
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	//if(!stl->g_data) stl->g_data = MEM_callocN(sizeof(*stl->g_data), __func__);

	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	SceneLANPR* lanpr = &draw_ctx->scene->lanpr;

	lanpr->reloaded = 1;

	if (!OneTime.InitComplete) {
		lanpr->depth_clamp = 0.01;
		lanpr->depth_strength = 800;
		lanpr->normal_clamp = 2;
		lanpr->normal_strength = 10;
		lanpr->line_thickness = 2;
		lanpr->taper_left_distance = 20;
		lanpr->taper_left_strength = 0.9;
		lanpr->taper_right_distance = 20;
		lanpr->taper_right_strength = 0.9;

		lanpr->line_color[0] = 0.22;
		lanpr->line_color[1] = 0.29;
		lanpr->line_color[2] = 0.53;
		lanpr->line_color[3] = 1;

		lanpr->background_color[0] = 0.59;
		lanpr->background_color[1] = 0.90;
		lanpr->background_color[2] = 0.51;
		lanpr->background_color[3] = 1;

		lanpr->reloaded = 1;

		OneTime.InitComplete=1;
	}


	/* Main Buffer */
	DRW_texture_ensure_fullscreen_2D(&txl->depth, GPU_DEPTH_COMPONENT32F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
	DRW_texture_ensure_fullscreen_2D(&txl->color, GPU_RGBA16F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
	DRW_texture_ensure_fullscreen_2D(&txl->normal, GPU_RGBA16F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
    DRW_texture_ensure_fullscreen_2D(&txl->edge_intermediate, GPU_RGBA16F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
	
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
		GPU_ATTACHMENT_TEXTURE(txl->depth),
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
			datatoc_lanpr_snake_multichannel_fragment,NULL,NULL,NULL);

    }
	if (!OneTime.edge_detect_shader) {
	OneTime.edge_detect_shader = 
		GPU_shader_create(
			datatoc_common_fullscreen_vert_glsl,
			datatoc_lanpr_snake_edge_fragment,NULL,NULL,NULL);

    }
	if (!OneTime.edge_thinning_shader) {
	OneTime.edge_thinning_shader = 
		GPU_shader_create(
			datatoc_common_fullscreen_vert_glsl,
			datatoc_lanpr_image_peel_fragment,NULL,NULL,NULL);

    }
	if (!OneTime.snake_connection_shader) {
	OneTime.snake_connection_shader = 
		GPU_shader_create(
			datatoc_lanpr_line_connection_vertex,
			datatoc_lanpr_line_connection_fragment,
			datatoc_lanpr_line_connection_geometry,
			NULL,NULL);
    }

	lanpr_init_atlas_inputs(ved);

}
static void lanpr_engine_free(void){
	void* ved= OneTime.ved;
    LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	//LANPR_ViewLayerData *sldata = EEVEE_view_layer_data_ensure();
    LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;

	DRW_pass_free(psl->color_pass);
	DRW_pass_free(psl->edge_intermediate);

	GPU_framebuffer_free(fbl->passes);
	GPU_framebuffer_free(fbl->edge_intermediate);
	GPU_framebuffer_free(fbl->edge_thinning);

	DRW_texture_free(txl->depth);
	DRW_texture_free(txl->color);
	DRW_texture_free(txl->normal);
	DRW_texture_free(txl->edge_intermediate);

	BLI_mempool_destroy(stl->g_data->mp_line_strip);
	BLI_mempool_destroy(stl->g_data->mp_line_strip_point);
	BLI_mempool_destroy(stl->g_data->mp_sample);

	MEM_freeN(stl->g_data->line_result_8bit);
	MEM_freeN(stl->g_data->line_result);
	MEM_freeN(stl->g_data);
	stl->g_data=0;
}

static void lanpr_cache_init(void *vedata){

	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;
	
	if (!stl->g_data) {
		/* Alloc transient pointers */
		stl->g_data = MEM_callocN(sizeof(*stl->g_data), __func__);
		stl->g_data->mp_sample = BLI_mempool_create(sizeof(LANPR_TextureSample), 0, 512, BLI_MEMPOOL_NOP);
		stl->g_data->mp_line_strip = BLI_mempool_create(sizeof(LANPR_LineStrip), 0, 512, BLI_MEMPOOL_NOP);
		stl->g_data->mp_line_strip_point = BLI_mempool_create(sizeof(LANPR_LineStripPoint), 0, 1024, BLI_MEMPOOL_NOP);
	}

	LANPR_PrivateData* pd = stl->g_data;

	psl->color_pass = DRW_pass_create("Color Pass", DRW_STATE_WRITE_COLOR | DRW_STATE_DEPTH_LESS_EQUAL | DRW_STATE_WRITE_DEPTH);
	stl->g_data->multipass_shgrp = DRW_shgroup_create(OneTime.multichannel_shader, psl->color_pass);


	struct Gwn_Batch *quad = DRW_cache_fullscreen_quad_get();


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
	DRW_shgroup_uniform_texture_ref(stl->g_data->edge_thinning_shgrp, "TexSample0", &txl->edge_intermediate);
	DRW_shgroup_uniform_int(stl->g_data->edge_thinning_shgrp, "Stage", &stl->g_data->stage, 1);
	DRW_shgroup_call_add(stl->g_data->edge_thinning_shgrp, quad, NULL);

	psl->edge_thinning_2 = DRW_pass_create("Edge Thinning Stage 2", DRW_STATE_WRITE_COLOR);
	stl->g_data->edge_thinning_shgrp_2 = DRW_shgroup_create(OneTime.edge_thinning_shader, psl->edge_thinning_2);
	DRW_shgroup_uniform_texture_ref(stl->g_data->edge_thinning_shgrp_2, "TexSample0", &txl->color);
	DRW_shgroup_uniform_int(stl->g_data->edge_thinning_shgrp_2, "Stage", &stl->g_data->stage, 1);
	DRW_shgroup_call_add(stl->g_data->edge_thinning_shgrp_2, quad, NULL);

	psl->dpix_transform_pass = DRW_pass_create("DPIX Transform Stage", DRW_STATE_WRITE_COLOR);
	stl->g_data->dpix_transform_shgrp = DRW_shgroup_create(OneTime.dpix_transform_shader, psl->dpix_transform_pass);
	DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "vert0_tex", &txl->dpix_in_pl);
	DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "vert1_tex", &txl->dpix_in_pr);
	DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "face_normal0_tex", &txl->dpix_in_nl);
	DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "face_normal1_tex", &txl->dpix_in_nr);
    DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "sample_step", &stl->g_data->dpix_sample_step, 1);
	DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "is_perspective", &stl->g_data->dpix_is_perspective, 1);
	DRW_shgroup_uniform_vec4(stl->g_data->dpix_transform_shgrp, "viewport", stl->g_data->dpix_viewport, 4);
    DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "buffer_width", &stl->g_data->dpix_buffer_width, 1);

	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;

	psl->dpix_preview_pass = DRW_pass_create("DPIX Preview", DRW_STATE_WRITE_COLOR|DRW_STATE_WRITE_DEPTH|DRW_STATE_DEPTH_LESS_EQUAL);
	stl->g_data->dpix_preview_shgrp = DRW_shgroup_create(OneTime.dpix_preview_shader, psl->dpix_preview_pass);
	DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_preview_shgrp, "vert0_tex", &txl->dpix_out_pl);
	DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_preview_shgrp, "vert1_tex", &txl->dpix_out_pr);
	DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "viewport", stl->g_data->dpix_viewport, 4);
	DRW_shgroup_uniform_vec4(stl->g_data->dpix_preview_shgrp, "color", lanpr->line_color, 4);
    DRW_shgroup_uniform_float(stl->g_data->dpix_preview_shgrp, "depth_offset", &stl->g_data->dpix_depth_offset, 1);
	
	pd->begin_index = 0;
	int tsize = sizeof(float) * 4 * TNS_DPIX_TEXTURE_SIZE*TNS_DPIX_TEXTURE_SIZE;
	if (!pd->atlas_pl) {
		pd->atlas_pl = MEM_callocN(tsize, "atlas_point_l");
		pd->atlas_pr = MEM_callocN(tsize, "atlas_point_r");
		pd->atlas_nl = MEM_callocN(tsize, "atlas_normal_l");
		pd->atlas_nr = MEM_callocN(tsize, "atlas_normal_l");
	}
	if(lanpr->reloaded){
	    memset(pd->atlas_pl, 0, tsize);
	    memset(pd->atlas_pr, 0, tsize);
	    memset(pd->atlas_nl, 0, tsize);
	    memset(pd->atlas_nr, 0, tsize);
	}
}

static void lanpr_cache_populate(void *vedata, Object *ob){
    
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData* pd = stl->g_data;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;
	
	if (!DRW_object_is_renderable(ob)) {
		return;
	}

	if (ob == draw_ctx->object_edit) {
		return;
	}

	struct Gwn_Batch *geom = DRW_cache_object_surface_get(ob);
	if (geom) {
        DRW_shgroup_call_object_add(stl->g_data->multipass_shgrp, geom, ob);
	}

	if(lanpr->reloaded){
		int idx = pd->begin_index;
		pd->begin_index = lanpr_feed_atlas_data_obj(vedata,pd->atlas_pl,pd->atlas_pr,pd->atlas_nl,pd->atlas_nr,ob,idx);
		lanpr_feed_atlas_trigger_preview_obj(vedata,ob,idx);
	}
}

static void lanpr_cache_finish(void *vedata){
    LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData* pd = stl->g_data;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;

	if(lanpr->reloaded){
		GPU_texture_update(txl->dpix_in_pl,pd->atlas_pl);
		GPU_texture_update(txl->dpix_in_pr,pd->atlas_pr);
		GPU_texture_update(txl->dpix_in_nl,pd->atlas_nl);
		GPU_texture_update(txl->dpix_in_nr,pd->atlas_nr);

		MEM_freeN(pd->atlas_pl);
		MEM_freeN(pd->atlas_pr);
		MEM_freeN(pd->atlas_nl);
		MEM_freeN(pd->atlas_nr);
		pd->atlas_pl = 0;
		lanpr->reloaded = 0;
	}
}

static void lanpr_draw_scene(void *vedata)
{
	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_FramebufferList *fbl = ((LANPR_Data *)vedata)->fbl;

	float clear_col[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	float clear_depth = 1.0f;
	uint clear_stencil = 0xFF;

	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
	DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();

    GPU_framebuffer_bind(fbl->passes);
    GPUFrameBufferBits clear_bits = GPU_DEPTH_BIT | GPU_COLOR_BIT;
    GPU_framebuffer_clear(fbl->passes, clear_bits, clear_col, clear_depth, clear_stencil);
   
    DRW_draw_pass(psl->color_pass);

	const DRWContextState *draw_ctx = DRW_context_state_get();
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;

	if(lanpr->master_mode == LANPR_MASTER_MODE_DPIX){
		lanpr_dpix_draw_scene(txl,fbl,psl,stl->g_data,lanpr);
	}else{
		lanpr_snake_draw_scene(txl,fbl,psl,stl->g_data,lanpr);
	}
}

static const DrawEngineDataSize lanpr_data_size = DRW_VIEWPORT_DATA_SIZE(LANPR_Data);

DrawEngineType draw_engine_lanpr_type = {
	NULL, NULL,
	N_("LANPR"),
	&lanpr_data_size,
	&lanpr_engine_init,
	&lanpr_engine_free,
	&lanpr_cache_init,
	&lanpr_cache_populate,
	&lanpr_cache_finish,
	NULL,//draw background
	lanpr_draw_scene,//draw scene, looks like that not much difference except a camera overlay image.
	NULL,
	NULL,
	NULL,
};

RenderEngineType DRW_engine_viewport_lanpr_type = {
	NULL, NULL,
	LANPR_ENGINE, N_("LANPR"), RE_INTERNAL,
	NULL,// update
	NULL,// render to img
	NULL,// bake
	NULL,// view update
	NULL,// render to view
	NULL,// update in script
	NULL,// update in render pass
	&draw_engine_lanpr_type,
	{NULL, NULL, NULL}
};