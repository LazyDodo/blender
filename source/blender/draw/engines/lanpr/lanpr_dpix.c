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
#include "DNA_meshdata_types.h"
#include "BKE_customdata.h"
#include "DEG_depsgraph_query.h"
#include "GPU_draw.h"

#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"
#include "bmesh.h"

#include <math.h>

extern LANPR_SharedResource lanpr_share;
extern char datatoc_lanpr_dpix_project_passthrough_vert_glsl[];
extern char datatoc_lanpr_dpix_project_clip_frag_glsl[];
extern char datatoc_lanpr_dpix_preview_geom_glsl[];
extern char datatoc_lanpr_dpix_preview_frag_glsl[];

void lanpr_init_atlas_inputs(void *ved){
	lanpr_share.ved_viewport = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d && rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;

	if (lanpr->reloaded || !txl->dpix_in_pl) {
		DRW_texture_ensure_2D(&txl->dpix_in_pl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_pr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_nl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_nr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_edge_mask, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA8, 0);
		DRW_texture_ensure_2D(&txl->dpix_out_pl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_out_pr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_out_length, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
	}

	GPU_framebuffer_ensure_config(&fbl->dpix_transform, {
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_TEXTURE(txl->dpix_out_pl),
		GPU_ATTACHMENT_TEXTURE(txl->dpix_out_pr),
		GPU_ATTACHMENT_TEXTURE(txl->dpix_out_length),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});

	GPU_framebuffer_ensure_config(&fbl->dpix_preview, {
		GPU_ATTACHMENT_TEXTURE(txl->depth),
		GPU_ATTACHMENT_TEXTURE(txl->color),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});

	if (!lanpr_share.dpix_transform_shader) {
		lanpr_share.dpix_transform_shader =
			DRW_shader_create(
				datatoc_lanpr_dpix_project_passthrough_vert_glsl,
				NULL,
				datatoc_lanpr_dpix_project_clip_frag_glsl,
				NULL);
	}
	if (!lanpr_share.dpix_preview_shader) {
		lanpr_share.dpix_preview_shader =
			DRW_shader_create(
				datatoc_lanpr_dpix_project_passthrough_vert_glsl,
				datatoc_lanpr_dpix_preview_geom_glsl,
				datatoc_lanpr_dpix_preview_frag_glsl,
				NULL);
	}
}
void lanpr_destroy_atlas(void *ved){
	//no need to free things, no custom data.
}

int lanpr_feed_atlas_data_obj(void *vedata,
                              float *AtlasPointsL, float *AtlasPointsR,
                              float *AtlasFaceNormalL, float *AtlasFaceNormalR,
                              float *AtlasEdgeMask,
                              Object *ob, int begin_index) {
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;

	if (!DRW_object_is_renderable(ob)) return begin_index;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	if (ob == draw_ctx->object_edit) return begin_index;
	if (ob->type != OB_MESH) return begin_index;

	Mesh *me = ob->data;
	BMesh *bm;
	struct BMFace *f1, *f2;
	struct BMVert *v1, *v2;
	struct BMEdge *e;
	struct BMLoop *l1, *l2;
	FreestyleEdge *fe;
	int CanFindFreestyle = 0;
	int vert_count = me->totvert, edge_count = me->totedge, face_count = me->totface;
	int i, idx;

	const BMAllocTemplate allocsize = BMALLOC_TEMPLATE_FROM_ME(me);
	bm = BM_mesh_create(&allocsize,
	                    &((struct BMeshCreateParams) {.use_toolflags = true, }));
	BM_mesh_bm_from_me(bm, me, &((struct BMeshFromMeshParams) {.calc_face_normal = true, }));
	BM_mesh_elem_table_ensure(bm, BM_VERT | BM_EDGE | BM_FACE);

	if (CustomData_has_layer(&bm->edata, CD_FREESTYLE_EDGE)) {
		CanFindFreestyle = 1;
	}

	for (i = 0; i < edge_count; i++) {
		f1 = 0;
		f2 = 0;
		e = BM_edge_at_index(bm, i);
		v1 = e->v1;
		v2 = e->v2;
		l1 = e->l;
		l2 = e->l ? e->l->radial_next : 0;
		if (l1) f1 = l1->f;
		if (l2) f2 = l2->f;

		idx = (begin_index + i) * 4;

		AtlasPointsL[idx + 0] = v1->co[0];
		AtlasPointsL[idx + 1] = v1->co[1];
		AtlasPointsL[idx + 2] = v1->co[2];
		AtlasPointsL[idx + 3] = 1;

		AtlasPointsR[idx + 0] = v2->co[0];
		AtlasPointsR[idx + 1] = v2->co[1];
		AtlasPointsR[idx + 2] = v2->co[2];
		AtlasPointsR[idx + 3] = 1;

		if (CanFindFreestyle) {
			fe = CustomData_bmesh_get(&bm->edata, e->head.data, CD_FREESTYLE_EDGE);
			if (fe->flag & FREESTYLE_EDGE_MARK) AtlasEdgeMask[idx + 1] = 1; // channel G
		}

		if (f1) {
			AtlasFaceNormalL[idx + 0] = f1->no[0];
			AtlasFaceNormalL[idx + 1] = f1->no[1];
			AtlasFaceNormalL[idx + 2] = f1->no[2];
			AtlasFaceNormalL[idx + 3] = 1;
		}
		else {
			AtlasFaceNormalL[idx + 0] = 0;
			AtlasFaceNormalL[idx + 1] = 0;
			AtlasFaceNormalL[idx + 2] = 0;
			AtlasFaceNormalL[idx + 3] = 0;
		}

		if (f2 && f2 != f1) { // this is for edge condition
			AtlasFaceNormalR[idx + 0] = f2->no[0];
			AtlasFaceNormalR[idx + 1] = f2->no[1];
			AtlasFaceNormalR[idx + 2] = f2->no[2];
			AtlasFaceNormalR[idx + 3] = 1;

			if (f2->mat_nr != f1->mat_nr) AtlasEdgeMask[idx] = 1; // channel R

		}
		else {
			AtlasFaceNormalR[idx + 0] = 0;
			AtlasFaceNormalR[idx + 1] = 0;
			AtlasFaceNormalR[idx + 2] = 0;
			AtlasFaceNormalR[idx + 3] = 0;
		}

	}

	BM_mesh_free(bm);

	return begin_index + edge_count;
}

int lanpr_feed_atlas_data_intersection_cache(void *vedata,
                                             float *AtlasPointsL, float *AtlasPointsR,
                                             float *AtlasFaceNormalL, float *AtlasFaceNormalR,
                                             float *AtlasEdgeMask,
                                             int begin_index) {
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData *pd = stl->g_data;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;
	LANPR_RenderBuffer *rb = lanpr->render_buffer;
	nListItemPointer *lip;
	LANPR_RenderLine *rl;
	int i, idx;

	i = 0;

	if (!rb) return 0;

	for (lip = rb->IntersectionLines.first; lip; lip = lip->pNext) {
		rl = lip->p;

		idx = (begin_index + i) * 4;
		AtlasEdgeMask[idx + 2] = 1; // channel B

		AtlasPointsL[idx + 0] = rl->L->GLocation[0];
		AtlasPointsL[idx + 1] = rl->L->GLocation[1];
		AtlasPointsL[idx + 2] = rl->L->GLocation[2];
		AtlasPointsL[idx + 3] = 1;

		AtlasPointsR[idx + 0] = rl->R->GLocation[0];
		AtlasPointsR[idx + 1] = rl->R->GLocation[1];
		AtlasPointsR[idx + 2] = rl->R->GLocation[2];
		AtlasPointsR[idx + 3] = 1;

		AtlasFaceNormalL[idx + 0] = 0;
		AtlasFaceNormalL[idx + 1] = 0;
		AtlasFaceNormalL[idx + 2] = 1;
		AtlasFaceNormalL[idx + 3] = 0;

		AtlasFaceNormalR[idx + 0] = 0;
		AtlasFaceNormalR[idx + 1] = 0;
		AtlasFaceNormalR[idx + 2] = 1;
		AtlasFaceNormalR[idx + 3] = 0;

		i++;
	}

	return begin_index + i;
}

void lanpr_dpix_index_to_coord(int index, float *x, float *y){
	(*x) = tnsLinearItp(-1, 1, (float)(index % TNS_DPIX_TEXTURE_SIZE + 0.5) / (float)TNS_DPIX_TEXTURE_SIZE);
	(*y) = tnsLinearItp(-1, 1, (float)(index / TNS_DPIX_TEXTURE_SIZE + 0.5) / (float)TNS_DPIX_TEXTURE_SIZE);
}

void lanpr_dpix_index_to_coord_absolute(int index, float *x, float *y){
	(*x) = (float)(index % TNS_DPIX_TEXTURE_SIZE) + 0.5;
	(*y) = (float)(index / TNS_DPIX_TEXTURE_SIZE) + 0.5;
}

int lanpr_feed_atlas_trigger_preview_obj(void *vedata, Object *ob, int begin_index) {
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData *pd = stl->g_data;
	Mesh *me = ob->data;
	if (ob->type != OB_MESH) return begin_index;
	int edge_count = me->totedge;
	int i;
	float co[2];

	static GPUVertFormat format = { 0 };
	static struct { uint pos, uvs; } attr_id;
	if (format.attr_len == 0) {
		attr_id.pos = GPU_vertformat_attr_add(&format, "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	}

	static GPUVertFormat format2 = { 0 };
	static struct { uint pos, uvs; } attr_id2;
	if (format2.attr_len == 0) {
		attr_id2.pos = GPU_vertformat_attr_add(&format2, "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	}

	GPUVertBuf *vbo = GPU_vertbuf_create_with_format(&format);
	GPUVertBuf *vbo2 = GPU_vertbuf_create_with_format(&format2);
	GPU_vertbuf_data_alloc(vbo, edge_count);
	GPU_vertbuf_data_alloc(vbo2, edge_count);

	for (i = 0; i < edge_count; i++) {
		lanpr_dpix_index_to_coord(i + begin_index, &co[0], &co[1]);
		GPU_vertbuf_attr_set(vbo, attr_id.pos, i, co);
		lanpr_dpix_index_to_coord_absolute(i + begin_index, &co[0], &co[1]);
		GPU_vertbuf_attr_set(vbo2, attr_id2.pos, i, co);
	}

	GPUBatch *gb = GPU_batch_create_ex(GPU_PRIM_POINTS, vbo, 0, GPU_USAGE_STATIC | GPU_BATCH_OWNS_VBO);
	GPUBatch *gb2 = GPU_batch_create_ex(GPU_PRIM_POINTS, vbo2, 0, GPU_USAGE_STATIC | GPU_BATCH_OWNS_VBO);

	LANPR_BatchItem *bi = BLI_mempool_alloc(pd->mp_batch_list);
	BLI_addtail(&pd->dpix_batch_list, bi);
	bi->dpix_transform_batch = gb;
	bi->dpix_preview_batch = gb2;
	bi->ob = ob;

	return begin_index + edge_count;
}

void lanpr_create_atlas_intersection_preview(void *vedata, int begin_index) {
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData *pd = stl->g_data;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	SceneLANPR *lanpr = &draw_ctx->scene->lanpr;
	LANPR_RenderBuffer *rb = lanpr->render_buffer;
	float co[2];
	int i;

	if (!rb) return;

	if (rb->DPIXIntersectionBatch) GPU_batch_discard(rb->DPIXIntersectionBatch);
	rb->DPIXIntersectionBatch = 0;

	if (!rb->IntersectionCount) return;

	static GPUVertFormat format = { 0 };
	static struct { uint pos, uvs; } attr_id;
	if (format.attr_len == 0) {
		attr_id.pos = GPU_vertformat_attr_add(&format, "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	}
	static GPUVertFormat format2 = { 0 };
	static struct { uint pos, uvs; } attr_id2;
	if (format2.attr_len == 0) {
		attr_id2.pos = GPU_vertformat_attr_add(&format2, "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	}

	GPUVertBuf *vbo = GPU_vertbuf_create_with_format(&format);
	GPU_vertbuf_data_alloc(vbo, rb->IntersectionCount);

	GPUVertBuf *vbo2 = GPU_vertbuf_create_with_format(&format2);
	GPU_vertbuf_data_alloc(vbo2, rb->IntersectionCount);

	for (i = 0; i < rb->IntersectionCount; i++) {
		lanpr_dpix_index_to_coord(i + begin_index, &co[0], &co[1]);
		GPU_vertbuf_attr_set(vbo, attr_id.pos, i, co);
		lanpr_dpix_index_to_coord_absolute(i + begin_index, &co[0], &co[1]);
		GPU_vertbuf_attr_set(vbo2, attr_id2.pos, i, co);
	}
	rb->DPIXIntersectionTransformBatch = GPU_batch_create_ex(GPU_PRIM_POINTS, vbo, 0, GPU_USAGE_STATIC | GPU_BATCH_OWNS_VBO);
	rb->DPIXIntersectionBatch = GPU_batch_create_ex(GPU_PRIM_POINTS, vbo2, 0, GPU_USAGE_STATIC | GPU_BATCH_OWNS_VBO);
}


void lanpr_dpix_draw_scene(LANPR_TextureList *txl, LANPR_FramebufferList *fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr, GPUFrameBuffer *DefaultFB, int is_render) {
	float clear_col[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	float clear_depth = 1.0f;
	uint clear_stencil = 0xFF;

	if (!lanpr->active_layer) return; /* return early in case we don't have line layers. DPIX only use the first layer. */

	int texw = GPU_texture_width(txl->edge_intermediate), texh = GPU_texture_height(txl->edge_intermediate);;
	int tsize = texw * texh;

	const DRWContextState *draw_ctx = DRW_context_state_get();
	Scene *scene = DEG_get_evaluated_scene(draw_ctx->depsgraph);
	View3D *v3d = draw_ctx->v3d;
	Object *camera = 0;
	if (v3d) {
		RegionView3D *rv3d = draw_ctx->rv3d;
		camera = (rv3d && rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	}
	if (!camera) {
		camera = scene->camera;
	}
	if (is_render && !camera) return;

	pd->dpix_viewport[2] = texw;
	pd->dpix_viewport[3] = texh;
	pd->dpix_is_perspective = 1;
	pd->dpix_sample_step = 1;
	pd->dpix_buffer_width = TNS_DPIX_TEXTURE_SIZE;
	pd->dpix_depth_offset = 0.0001;
	pd->dpix_znear = camera ? ((Camera *)camera->data)->clipsta : v3d->near;
	pd->dpix_zfar = camera ? ((Camera *)camera->data)->clipend : v3d->far;

	GPU_point_size(1);
	//GPU_line_width(2);
	GPU_framebuffer_bind(fbl->dpix_transform);
	DRW_draw_pass(psl->dpix_transform_pass);

	GPU_framebuffer_bind(fbl->dpix_preview);
	GPUFrameBufferBits clear_bits = GPU_COLOR_BIT;
	GPU_framebuffer_clear(fbl->dpix_preview, clear_bits, lanpr->background_color, clear_depth, clear_stencil);
	DRW_draw_pass(psl->dpix_preview_pass);

	GPU_framebuffer_bind(DefaultFB);
	GPU_framebuffer_clear(DefaultFB, clear_bits, lanpr->background_color, clear_depth, clear_stencil);
	DRW_multisamples_resolve(txl->depth, txl->color, 1);
}
