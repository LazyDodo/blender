#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
#include "NUL_TNS.h"
#include "DRW_render.h"
#include "BKE_object.h"
#include "DNA_camera_types.h"
#include "GPU_immediate.h"
#include "GPU_immediate_util.h"
#include "GPU_framebuffer.h"

#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"


#include <math.h>

extern char datatoc_common_fullscreen_vert_glsl[];
extern char datatoc_gpu_shader_3D_normal_smooth_color_vert_glsl[];
extern char datatoc_lanpr_snake_multichannel_fragment[];
extern char datatoc_lanpr_snake_edge_fragment[];
extern char datatoc_lanpr_image_peel_fragment[];
extern char datatoc_lanpr_line_connection_vertex[];
extern char datatoc_lanpr_line_connection_fragment[];
extern char datatoc_lanpr_line_connection_geometry[];

//==============================================================[ ATLAS / DPIX ]

// will be updated here very soon(-ish).....


//=====================================================================[ SNAKE ]


//==============================================[ ENGINE ]

typedef struct LANPROneTimeInit{
    GPUShader* multichannel_shader;
	GPUShader* edge_detect_shader;
	GPUShader* edge_thinning_shader;
	GPUShader* snake_connection_shader;
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

	//if(!stl->g_data) stl->g_data = MEM_callocN(sizeof(*stl->g_data), __func__);

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
	if (!OneTime.snake_connection_shader) {
	OneTime.snake_connection_shader = 
		GPU_shader_create(
			datatoc_lanpr_line_connection_vertex,
			datatoc_lanpr_line_connection_fragment,
			datatoc_lanpr_line_connection_geometry,
			NULL,NULL);
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

int _TNS_ColOffsets[] = { -1,0,1,1,1,0,-1,-1 };
int _TNS_RowOffsets[] = { -1,-1,-1,0,1,1,1,0 };

int _TNS_Deviates[8][8] = {
	{ 0,1,2,3,4,3,2,1 },
	{ 1,0,1,2,3,4,3,2 },
	{ 2,1,0,1,2,3,4,3 },
	{ 3,2,1,0,1,2,3,4 },
	{ 4,3,2,1,0,1,2,3 },
	{ 3,4,3,2,1,0,1,2 },
	{ 2,3,4,3,2,1,0,1 },
	{ 1,2,3,4,3,2,1,0 }
};

#define TNS_CLAMP_TEXTURE_W(t,Col)\
	{if (Col >= t->width) Col = t->width - 1; if (Col < 0) Col = 0;}

#define TNS_CLAMP_TEXTURE_H(t,Row)\
	{if (Row >= t->height) Row = t->height - 1;if (Row < 0) Row = 0;}

#define TNS_CLAMP_TEXTURE_CONTINUE(t,Col,Row)\
	{if (Col >= t->width) continue; if (Col < 0) continue;\
     if (Row >= t->height) continue; if (Row < 0)continue; }


static LANPR_TextureSample* lanpr_any_uncovered_samples(LANPR_PrivateData* pd){
	return BLI_pophead(&pd->pending_samples);
}

int lanpr_direction_deviate(int From, int To) {
	return _TNS_Deviates[From - 1][To - 1];
}

int lanpr_detect_direction(LANPR_PrivateData* pd, int Col, int Row, int LastDirection) {
	int Deviate[9] = {100};
	int MinDeviate = 0;
	int i;
	LANPR_TextureSample* ts;

	for (i = 0; i < 8; i++) {
		TNS_CLAMP_TEXTURE_CONTINUE(pd, (_TNS_ColOffsets[i] + Col), (_TNS_RowOffsets[i] + Row));
		if (ts = pd->sample_table[(_TNS_ColOffsets[i] + Col) + (_TNS_RowOffsets[i] + Row) * pd->width]) {
			if (!LastDirection) return i + 1;
			Deviate[i+1] = lanpr_direction_deviate(i, LastDirection);
			if (!MinDeviate || Deviate[MinDeviate] > Deviate[i + 1]) MinDeviate = i + 1;
		}
	}

	return MinDeviate;
}

LANPR_LineStrip* lanpr_create_line_strip(LANPR_PrivateData* pd) {
	LANPR_LineStrip* ls = BLI_mempool_calloc(pd->mp_line_strip);
	return ls;
}
LANPR_LineStripPoint* lanpr_append_point(LANPR_PrivateData* pd, LANPR_LineStrip* ls, real X, real Y, real Z) {
	LANPR_LineStripPoint* lsp = BLI_mempool_calloc(pd->mp_line_strip_point);

	lsp->P[0] = X;
	lsp->P[1] = Y;
	lsp->P[2] = Z;

	BLI_addtail(&ls->points,lsp);

	ls->point_count++;

	return lsp;
}
LANPR_LineStripPoint* lanpr_push_point(LANPR_PrivateData* pd, LANPR_LineStrip* ls, real X, real Y, real Z) {
	LANPR_LineStripPoint* lsp = BLI_mempool_calloc(pd->mp_line_strip_point);

	lsp->P[0] = X;
	lsp->P[1] = Y;
	lsp->P[2] = Z;

	BLI_addhead(&ls->points, lsp);

	ls->point_count++;

	return lsp;
}

//void lanpr_remove_point(LANPR_LineStrip* ls,tnsLineStripPoint* lsp) {
//	lstRemoveItem(&ls->Points, lsp);
//	FreeMem(lsp);
//}

void lanpr_destroy_line_strip(LANPR_PrivateData* pd, LANPR_LineStrip* ls) {
	LANPR_LineStripPoint* lsp;
	while (lsp = BLI_pophead(&ls->points)) {
		BLI_mempool_free(pd->mp_line_strip_point, lsp);
	}
	BLI_mempool_free(pd->mp_line_strip, ls);
}

void lanpr_remove_sample(LANPR_PrivateData* pd, int Row, int Col) {
	LANPR_TextureSample* ts;
	ts = pd->sample_table[Row*pd->width + Col];
	pd->sample_table[Row*pd->width + Col] = 0;
	
	BLI_remlink(&pd->pending_samples, ts);
	ts->Item.prev=0;ts->Item.next=0;
	BLI_addtail(&pd->erased_samples, ts);
}

int lanpr_grow_snake_r(LANPR_PrivateData* pd, LANPR_LineStrip* ls, LANPR_LineStripPoint* ThisP, int Direction) {
	LANPR_LineStripPoint* NewP = ThisP,*p2;
	int Length = 5;
	int l = 0;
	int Deviate,Dir=Direction,NewDir;
	int AddPoint;
	int TX = NewP->P[0], TY = NewP->P[1];

	while (NewDir = lanpr_detect_direction(pd, TX, TY, Dir)) {
		AddPoint = 0;
		Deviate = lanpr_direction_deviate(NewDir, Dir);
		Dir = NewDir;

		l++;
		TX += _TNS_ColOffsets[NewDir-1];
		TY += _TNS_RowOffsets[NewDir-1];

		if (Deviate < 2) {
			lanpr_remove_sample(pd, TY, TX);
		}elif(Deviate < 3) {
			lanpr_remove_sample(pd, TY, TX);
			AddPoint = 1;
		}else {
			lanpr_remove_sample(pd, TY, TX);
			return;
		}

		if (AddPoint || l == Length) {
			p2 = lanpr_append_point(pd, ls, TX, TY, 0);
			NewP = p2;
			l = 0;
		}
	}
	if(TX!=ThisP->P[0] || TY!=ThisP->P[1])
		lanpr_append_point(pd, ls, TX, TY, 0);
}

int lanpr_grow_snake_l(LANPR_PrivateData* pd, LANPR_LineStrip* ls, LANPR_LineStripPoint* ThisP, int Direction) {
	LANPR_LineStripPoint* NewP = ThisP, *p2;
	int Length = 5;
	int l = 0;
	int Deviate, Dir = Direction, NewDir;
	int AddPoint;
	int TX = NewP->P[0], TY = NewP->P[1];

	while (NewDir = lanpr_detect_direction(pd, TX, TY, Dir)) {
		AddPoint = 0;
		Deviate = lanpr_direction_deviate(NewDir, Dir);
		Dir = NewDir;

		l++;
		TX += _TNS_ColOffsets[NewDir - 1];
		TY += _TNS_RowOffsets[NewDir - 1];

		if (Deviate < 2) {
			lanpr_remove_sample(pd, TY, TX);
		}elif(Deviate < 4) {
			lanpr_remove_sample(pd, TY, TX);
			AddPoint = 1;
		}
		else {
			lanpr_remove_sample(pd, TY, TX);
			return;
		}

		if (AddPoint || l == Length) {
			p2 = lanpr_push_point(pd, ls, TX, TY, 0);
			NewP = p2;
			l = 0;
		}
	}
	if(TX!=ThisP->P[0] || TY!=ThisP->P[1])
		lanpr_push_point(pd, ls, TX, TY, 0);
}

int lanpr_reverse_direction(int From) {
	From -= 4;
	if (From <= 0)From += 8;
	return From;
}

#define tMatDist2v(p1,p2)\
    sqrt(((p1)[0]-(p2)[0])*((p1)[0]-(p2)[0]) + ((p1)[1]-(p2)[1])*((p1)[1]-(p2)[1]))

#define tnsLinearItp(L,R,T)\
((L)*(1.0f - (T)) + (R)*(T))

void lanpr_texture_to_ndc(int x,int y, int w,int h, float* x_ndc, float* y_ndc){
    *x_ndc = tnsLinearItp(-1,1,(float)x/(float)w);
	*y_ndc = tnsLinearItp(-1,1, (float)y/(float)h);
}

Gwn_Batch *lanpr_get_snake_batch(LANPR_PrivateData* pd, LANPR_LineStrip* ls){
	int Count = ls->point_count;
	int i;
	LANPR_LineStripPoint* lsp,*plsp;
	u32bit *Index_adjacent;
	float* Verts;
	float* Lengths;
	float TotalLength=0;

	Index_adjacent = MEM_callocN(sizeof(unsigned int) * (Count - 1) * 4, "Index_adjacent buffer pre alloc");
	Verts = MEM_callocN(sizeof(float) * Count * 2, "Verts buffer pre alloc");
	Lengths = MEM_callocN(sizeof(float)* Count, "Length buffer pre alloc");

	Gwn_IndexBufBuilder elb;
	GWN_indexbuf_init_ex(&elb, GWN_PRIM_LINES_ADJ, (Count - 1) * 4, Count, true);

	for (i = 0; i < Count-1; i++) {
		Index_adjacent[i * 4 + 0] = i - 1;
		Index_adjacent[i * 4 + 1] = i;
		Index_adjacent[i * 4 + 2] = i + 1;
		Index_adjacent[i * 4 + 3] = i + 2;
		GWN_indexbuf_add_line_adj_verts(&elb, (i-1<0?0:i-1), i, i+1, (i+2>=Count-1?Count-1:i+2));
	}
	Index_adjacent[0] = 0;
	Index_adjacent[(Count - 1) * 4 - 1] = Count-1;

	i = 0;
	float xf,yf;
	for (lsp = (LANPR_LineStripPoint *)(ls->points.first); lsp; lsp = (LANPR_LineStripPoint *)(lsp->Item.next)) {
		lanpr_texture_to_ndc(lsp->P[0],lsp->P[1],pd->width, pd->height, &xf,&yf);
		Verts[i * 2 + 0] = xf;
		Verts[i * 2 + 1] = yf;
		if (plsp = (LANPR_LineStripPoint *)(lsp->Item.prev)) {
			TotalLength += tMatDist2v(plsp->P, lsp->P);
			Lengths[i] = TotalLength;
		}
		i++;
	}
	ls->total_length = TotalLength;

	static Gwn_VertFormat format = { 0 };
	static struct { uint pos, uvs; } attr_id;
	if (format.attrib_ct == 0) {
		attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
		attr_id.uvs = GWN_vertformat_attr_add(&format, "uvs", GWN_COMP_F32, 1, GWN_FETCH_FLOAT);
	}

	Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);
	GWN_vertbuf_data_alloc(vbo, Count);

	for (int i = 0; i < Count; ++i) {
		GWN_vertbuf_attr_set(vbo, attr_id.pos, i, &Verts[i*2]);
		GWN_vertbuf_attr_set(vbo, attr_id.uvs, i, &Lengths[i]);
	}

	MEM_freeN(Index_adjacent);
	MEM_freeN(Verts);
	MEM_freeN(Lengths);
	
	return GWN_batch_create_ex(GWN_PRIM_LINES_ADJ, vbo, GWN_indexbuf_build(&elb), GWN_USAGE_STREAM);
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
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;
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
	clear_bits = GPU_COLOR_BIT;
	GPU_framebuffer_clear(fbl->edge_thinning, clear_bits, clear_col, clear_depth, clear_stencil);
    DRW_draw_pass(psl->edge_thinning);

	stl->g_data->stage = 1;
	GPU_framebuffer_bind(fbl->edge_intermediate);
	//GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
    DRW_draw_pass(psl->edge_thinning_2);

	stl->g_data->stage = 0;
    GPU_framebuffer_bind(fbl->edge_thinning);
	GPU_framebuffer_clear(fbl->edge_thinning, clear_bits, clear_col, clear_depth, clear_stencil);
    DRW_draw_pass(psl->edge_thinning);

	stl->g_data->stage = 1;
	GPU_framebuffer_bind(fbl->edge_intermediate);
	//GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
    DRW_draw_pass(psl->edge_thinning_2);

	int texw = GPU_texture_width(txl->edge_intermediate) ,texh = GPU_texture_height(txl->edge_intermediate);;
    int size = texw*texh;
	int recreate=0;
	if(size != stl->g_data->width*stl->g_data->height) recreate =1;

	if(recreate){
		if(stl->g_data->line_result) MEM_freeN(stl->g_data->line_result);
		stl->g_data->line_result = MEM_callocN(sizeof(float) * size,"Texture readback buffer");

		if(stl->g_data->line_result_8bit) MEM_freeN(stl->g_data->line_result_8bit);
        stl->g_data->line_result_8bit = MEM_callocN(sizeof(unsigned char) * size,"Texture readback buffer 8bit");

		if(stl->g_data->sample_table) MEM_freeN(stl->g_data->sample_table);
		stl->g_data->sample_table = MEM_callocN(sizeof(void*) * size,"Texture readback buffer 8bit");

		stl->g_data->width = texw;
		stl->g_data->height = texh;
	}

	GPU_framebuffer_read_color(fbl->edge_intermediate,0,0,texw, texh,1,0, stl->g_data->line_result);

	float sample;
	int h, w;
	for (h = 0; h < texh; h++) {
		for (w = 0; w < texw; w++) {
			int index = h*texw + w;
			if ((sample = stl->g_data->line_result[index]) > 0.9) {
				stl->g_data->line_result_8bit[index] = 255;
				LANPR_TextureSample* ts = BLI_mempool_calloc(stl->g_data->mp_sample);
				BLI_addtail(&stl->g_data->pending_samples, ts);
				stl->g_data->sample_table[index] = ts;
				ts->X = w;
				ts->Y = h;
			}else{
				stl->g_data->sample_table[index] = 0;
			}
		}
	}

	LANPR_TextureSample *ts;
    LANPR_LineStrip* ls;
	LANPR_LineStripPoint* lsp;
	while(ts = lanpr_any_uncovered_samples(stl->g_data)){
		int Direction=0;
		LANPR_LineStripPoint tlsp = { 0 };

		tlsp.P[0] = ts->X;
		tlsp.P[1] = ts->Y;

		if (Direction = lanpr_detect_direction(stl->g_data, ts->X,ts->Y, Direction)) {
			BLI_addtail(&stl->g_data->line_strips, (ls = lanpr_create_line_strip(stl->g_data)));
			lsp = lanpr_append_point(stl->g_data, ls, ts->X, ts->Y, 0);
			lanpr_remove_sample(stl->g_data, ts->Y, ts->X);

			lanpr_grow_snake_r(stl->g_data, ls, lsp, Direction);

			lanpr_grow_snake_l(stl->g_data, ls, lsp, lanpr_reverse_direction(Direction));
		}

		//count++;
	}

	//GPU_framebuffer_bind()
	GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, lanpr->background_color, clear_depth, clear_stencil);


	for (ls = (LANPR_LineStrip *)(stl->g_data->line_strips.first); ls; ls = (LANPR_LineStrip *)(ls->Item.next)) {
		if (ls->point_count < 2) continue;

		float* tld = &lanpr->taper_left_distance, *tls = &lanpr->taper_left_strength,
			*trd = &lanpr->taper_right_distance, *trs = &lanpr->taper_right_strength;

		Gwn_Batch* snake_batch = lanpr_get_snake_batch(stl->g_data,ls);

		psl->snake_pass = DRW_pass_create("Snake Visualization Pass", DRW_STATE_WRITE_COLOR);
		stl->g_data->snake_shgrp = DRW_shgroup_create(OneTime.snake_connection_shader, psl->snake_pass);
		DRW_shgroup_uniform_float(stl->g_data->snake_shgrp, "LineWidth", &lanpr->line_thickness, 1);
		DRW_shgroup_uniform_float(stl->g_data->snake_shgrp, "TotalLength", &ls->total_length, 1);
		DRW_shgroup_uniform_float(stl->g_data->snake_shgrp, "TaperLDist", tld, 1);
		DRW_shgroup_uniform_float(stl->g_data->snake_shgrp, "TaperLStrength", tls, 1);
		DRW_shgroup_uniform_float(stl->g_data->snake_shgrp, "TaperRDist", lanpr->use_same_taper?tld:trd, 1);
		DRW_shgroup_uniform_float(stl->g_data->snake_shgrp, "TaperRStrength", lanpr->use_same_taper?tls:trs, 1);
		//lanpr->line_color[0] = lanpr->line_color[1] =lanpr->line_color[2] = lanpr->line_color[3] =1;
		DRW_shgroup_uniform_vec4(stl->g_data->snake_shgrp, "LineColor", lanpr->line_color, 1);

		DRW_shgroup_call_add(stl->g_data->snake_shgrp, snake_batch, NULL);
		DRW_draw_pass(psl->snake_pass);
	}
	

	BLI_mempool_clear(stl->g_data->mp_sample);
	BLI_mempool_clear(stl->g_data->mp_line_strip);
	BLI_mempool_clear(stl->g_data->mp_line_strip_point);

	stl->g_data->pending_samples.first = stl->g_data->pending_samples.last = 0;
	stl->g_data->erased_samples.first = stl->g_data->erased_samples.last = 0;
	stl->g_data->line_strips.first = stl->g_data->line_strips.last = 0;
	//stl->g_data->line_strip_point.first = stl->g_data->line_strip_point.last = 0;


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