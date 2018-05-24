#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "NUL_TNS.h"
#include "DRW_render.h"
#include "BKE_object.h"
#include "DNA_camera_types.h"
#include "GPU_immediate.h"
#include "GPU_immediate_util.h"

#include <math.h>

extern char datatoc_common_fullscreen_vert_glsl[];
extern char datatoc_gpu_shader_3D_normal_smooth_color_vert_glsl[];
extern char datatoc_lanpr_snake_multichannel_fragment[];
extern char datatoc_lanpr_snake_edge_fragment[];
extern char datatoc_lanpr_image_peel_fragment[];

//==============================================================[ ATLAS / DPIX ]

// will be updated here very soon(-ish).....


//=====================================================================[ SNAKE ]


//==============================================[ ENGINE ]

typedef struct LANPROneTimeInit{
    GPUShader* multichannel_shader;
	GPUShader* edge_detect_shader;
	GPUShader* edge_thinning_shader;
	void* ved;
} LANPROneTimeInit;

LANPROneTimeInit OneTime;

static void lanpr_engine_init(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	//LANPR_ViewLayerData *sldata = EEVEE_view_layer_data_ensure();
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;

	/* Main Buffer */
	DRW_texture_ensure_fullscreen_2D(&txl->depth, GPU_DEPTH_COMPONENT32F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
	DRW_texture_ensure_fullscreen_2D(&txl->color, GPU_RGBA16F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
	DRW_texture_ensure_fullscreen_2D(&txl->normal, GPU_RGBA16F, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
    DRW_texture_ensure_fullscreen_2D(&txl->edge_intermediate, GPU_RGBA8, DRW_TEX_FILTER | DRW_TEX_MIPMAP);
	
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
}

static void lanpr_cache_init(void *vedata){

	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_TextureList *txl = ((LANPR_Data *)vedata)->txl;


	if (!stl->g_data) {
		/* Alloc transient pointers */
		stl->g_data = MEM_mallocN(sizeof(*stl->g_data), __func__);
	}


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
}

static void lanpr_cache_populate(void *vedata, Object *ob){
    
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	
	if (!DRW_object_is_renderable(ob)) {
		return;
	}

	const DRWContextState *draw_ctx = DRW_context_state_get();
	if (ob == draw_ctx->object_edit) {
		return;
	}

	struct Gwn_Batch *geom = DRW_cache_object_surface_get(ob);
	if (geom) {
        DRW_shgroup_call_object_add(stl->g_data->multipass_shgrp, geom, ob);
	}
}

static void lanpr_cache_finish(void *vedata){

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
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;

	stl->g_data->znear = camera? ((Camera*)camera->data)->clipsta:0.1;
    stl->g_data->zfar = camera? ((Camera*)camera->data)->clipend:100;
	stl->g_data->normal_clamp =    draw_ctx->scene->lanpr.normal_clamp;
	stl->g_data->normal_strength = draw_ctx->scene->lanpr.normal_strength;
	stl->g_data->depth_clamp =     draw_ctx->scene->lanpr.depth_clamp;
	stl->g_data->depth_strength =  draw_ctx->scene->lanpr.depth_strength;

	GPU_framebuffer_bind(fbl->edge_intermediate);
	//GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
	
	DRW_draw_pass(psl->edge_intermediate);
	
    stl->g_data->stage = 0;
    GPU_framebuffer_bind(fbl->edge_thinning);
	GPU_framebuffer_clear(fbl->edge_thinning, clear_bits, clear_col, clear_depth, clear_stencil);
    DRW_draw_pass(psl->edge_thinning);

	stl->g_data->stage = 1;
	GPU_framebuffer_bind(fbl->edge_intermediate);
	//GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
    DRW_draw_pass(psl->edge_thinning_2);

	GPU_framebuffer_bind(dfbl->default_fb);

	DRW_transform_to_display(txl->edge_intermediate);
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