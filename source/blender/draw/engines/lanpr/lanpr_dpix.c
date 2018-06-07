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

extern struct LANPROneTimeInit OneTime;
extern char datatoc_lanpr_atlas_project_passthrough_vertex[];
extern char datatoc_lanpr_atlas_project_clip_fragment[];
extern char datatoc_lanpr_atlas_preview_geometry[];
extern char datatoc_lanpr_atlas_preview_fragment[];

void lanpr_init_atlas_inputs(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	//LANPR_ViewLayerData *sldata = EEVEE_view_layer_data_ensure();
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	//txl->dpix_in_pl = 

	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
	SceneLANPR* lanpr=&draw_ctx->scene->lanpr;

	if(lanpr->reloaded || !txl->dpix_in_pl){
		DRW_texture_ensure_2D(&txl->dpix_in_pl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_pr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_nl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_in_nr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);

		DRW_texture_ensure_2D(&txl->dpix_out_pl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_out_pr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
		DRW_texture_ensure_2D(&txl->dpix_out_length, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGBA32F, 0);
	}


	/* Main Buffer */

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

    if (!OneTime.dpix_transform_shader) {
	OneTime.dpix_transform_shader = 
		GPU_shader_create(
			datatoc_lanpr_atlas_project_passthrough_vertex,
			datatoc_lanpr_atlas_project_clip_fragment,
			NULL,NULL,NULL);
    }
	if (!OneTime.dpix_preview_shader) {
	OneTime.dpix_preview_shader = 
		GPU_shader_create(
		    datatoc_lanpr_atlas_project_passthrough_vertex,
			datatoc_lanpr_atlas_preview_fragment,
			datatoc_lanpr_atlas_preview_geometry,
			//NULL,
			NULL,NULL);
    }
}
void lanpr_destroy_atlas(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PassList *psl = ((LANPR_Data *)vedata)->psl;

	DRW_pass_free(psl->dpix_transform_pass);
	DRW_pass_free(psl->dpix_preview_pass);

	GPU_framebuffer_free(fbl->dpix_transform);
	GPU_framebuffer_free(fbl->dpix_preview);

	DRW_texture_free(txl->dpix_in_pl);
	DRW_texture_free(txl->dpix_in_pr);
	DRW_texture_free(txl->dpix_in_nl);
	DRW_texture_free(txl->dpix_in_nr);
	DRW_texture_free(txl->dpix_out_pl);
	DRW_texture_free(txl->dpix_out_pr);
}

int lanpr_feed_atlas_data_obj(void* vedata,
	float* AtlasPointsL, float* AtlasPointsR,
	float* AtlasFaceNormalL, float* AtlasFaceNormalR,
	Object* ob, int BeginIndex) {
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;

	if (!DRW_object_is_renderable(ob)) return BeginIndex;
	const DRWContextState *draw_ctx = DRW_context_state_get();
	if (ob == draw_ctx->object_edit) return BeginIndex;
	if(ob->type != OB_MESH) return BeginIndex;

	Mesh* me = ob->data;
	BMesh* bm;
	struct BMFace *f1,*f2;
	struct BMVert *v1,*v2;
	struct BMEdge *e;
	struct BMLoop *l1,*l2;
    int vert_count = me->totvert, edge_count = me->totedge, face_count = me->totface;
	int i,idx;

	const BMAllocTemplate allocsize = BMALLOC_TEMPLATE_FROM_ME(me);
	bm = BM_mesh_create(&allocsize,
			            &((struct BMeshCreateParams){.use_toolflags = true,}));
	BM_mesh_bm_from_me(bm, me, &((struct BMeshFromMeshParams){.calc_face_normal = true,}));
	BM_mesh_elem_table_ensure(bm,BM_VERT|BM_EDGE|BM_FACE);
	
	for(i=0; i<edge_count; i++){
		f1=0;
		f2=0;
		e = BM_edge_at_index(bm,i);
		v1 = e->v1;
		v2 = e->v2;
		l1 = e->l;
		l2 = e->l?e->l->radial_next:0;
		if(l1) f1 = l1->f;
		if(l2) f2 = l2->f;

		idx = (BeginIndex+i)*4;

		AtlasPointsL[idx + 0] = v1->co[0];
		AtlasPointsL[idx + 1] = v1->co[1];
		AtlasPointsL[idx + 2] = v1->co[2];
		AtlasPointsL[idx + 3] = 1;

		AtlasPointsR[idx + 0] = v2->co[0];
		AtlasPointsR[idx + 1] = v2->co[1];
		AtlasPointsR[idx + 2] = v2->co[2];
		AtlasPointsR[idx + 3] = 1;

		if(f1){
			AtlasFaceNormalL[idx + 0] = f1->no[0];
			AtlasFaceNormalL[idx + 1] = f1->no[1];
			AtlasFaceNormalL[idx + 2] = f1->no[2];
			AtlasFaceNormalL[idx + 3] = 1;
		}else{
			AtlasFaceNormalL[idx + 0] = 0;
			AtlasFaceNormalL[idx + 1] = 0;
			AtlasFaceNormalL[idx + 2] = 0;
			AtlasFaceNormalL[idx + 3] = 0;
		}

		if(f2 && f2!=f1){ // this is for edge condition
			AtlasFaceNormalR[idx + 0] = f2->no[0];
			AtlasFaceNormalR[idx + 1] = f2->no[1];
			AtlasFaceNormalR[idx + 2] = f2->no[2];
			AtlasFaceNormalR[idx + 3] = 1;
		}else{
			AtlasFaceNormalR[idx + 0] = 0;
			AtlasFaceNormalR[idx + 1] = 0;
			AtlasFaceNormalR[idx + 2] = 0;
			AtlasFaceNormalR[idx + 3] = 0;
		}

	}

	BM_mesh_free(bm);
	
	return BeginIndex + edge_count;
}

void lanpr_dpix_index_to_coord(int index, float* x,float* y){
    (*x) = tnsLinearItp(-1,1,(float)(index % TNS_DPIX_TEXTURE_SIZE+0.5)/(float)TNS_DPIX_TEXTURE_SIZE);
	(*y) = tnsLinearItp(-1,1,(float)(index / TNS_DPIX_TEXTURE_SIZE+0.5)/(float)TNS_DPIX_TEXTURE_SIZE);
}
void lanpr_dpix_index_to_coord_absolute(int index, float* x,float* y){
	(*x) = (float)(index % TNS_DPIX_TEXTURE_SIZE)+0.5;
    (*y) = (float)(index / TNS_DPIX_TEXTURE_SIZE)+0.5;
}

int lanpr_feed_atlas_trigger_preview_obj(void* vedata, Object* ob, int BeginIndex) {
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	LANPR_PrivateData* pd = stl->g_data;
	Mesh* me = ob->data;
	if (ob->type != OB_MESH) return BeginIndex;
	int edge_count = me->totedge;
	int i;
	float co[2];

	static Gwn_VertFormat format = { 0 };
	static struct { uint pos, uvs; } attr_id;
	if (format.attrib_ct == 0) {
		attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
	}

	static Gwn_VertFormat format2 = { 0 };
	static struct { uint pos, uvs; } attr_id2;
	if (format2.attrib_ct == 0) {
		attr_id2.pos = GWN_vertformat_attr_add(&format2, "pos", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
	}

	Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);
	Gwn_VertBuf *vbo2 = GWN_vertbuf_create_with_format(&format2);
	GWN_vertbuf_data_alloc(vbo, edge_count);
	GWN_vertbuf_data_alloc(vbo2, edge_count);

	for(i=0;i<edge_count;i++){
        lanpr_dpix_index_to_coord(i+BeginIndex,&co[0],&co[1]);
		GWN_vertbuf_attr_set(vbo, attr_id.pos, i, co);
		lanpr_dpix_index_to_coord_absolute(i+BeginIndex,&co[0],&co[1]);
		GWN_vertbuf_attr_set(vbo2, attr_id2.pos, i, co);
	}
	
	Gwn_Batch* gb = GWN_batch_create_ex(GWN_PRIM_POINTS, vbo, 0, GWN_USAGE_STATIC|GWN_BATCH_OWNS_VBO);
    Gwn_Batch* gb2 = GWN_batch_create_ex(GWN_PRIM_POINTS, vbo2, 0, GWN_USAGE_STATIC|GWN_BATCH_OWNS_VBO);

	LANPR_BatchItem *bi = BLI_mempool_alloc(pd->mp_batch_list);
	BLI_addtail(&pd->dpix_batch_list,bi);
	bi->dpix_transform_batch = gb;
	bi->dpix_preview_batch = gb2;
	bi->ob = ob;

	return BeginIndex + edge_count;
}


void lanpr_dpix_draw_scene(LANPR_TextureList* txl, LANPR_FramebufferList * fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr){
    	float clear_col[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	    float clear_depth = 1.0f;
	    uint clear_stencil = 0xFF;

        DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
	    DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();

        int texw = GPU_texture_width(txl->edge_intermediate) ,texh = GPU_texture_height(txl->edge_intermediate);;
	    int tsize = texw*texh;
        
        pd->dpix_viewport[2] = texw;
		pd->dpix_viewport[3] = texh;
		pd->dpix_is_perspective = 1;
		pd->dpix_sample_step = 1;
		pd->dpix_buffer_width = TNS_DPIX_TEXTURE_SIZE;
		pd->dpix_depth_offset=0.0001;

        glPointSize(1);
		glLineWidth(2);
		GPU_framebuffer_bind(fbl->dpix_transform);
		//GPU_disable_program_point_size();
		DRW_draw_pass(psl->dpix_transform_pass);

		//GPU_framebuffer_bind(fbl->edge_intermediate);
		//DRW_draw_pass(psl->color_pass);// use depth

		glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH, GL_NICEST);
	    glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		GPU_framebuffer_bind(fbl->dpix_preview);
		GPUFrameBufferBits clear_bits = GPU_COLOR_BIT;
		GPU_framebuffer_clear(fbl->dpix_preview, clear_bits, lanpr->background_color, clear_depth, clear_stencil);
		DRW_draw_pass(psl->dpix_preview_pass);

		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);

		GPU_framebuffer_bind(dfbl->default_fb);
		//DRW_transform_to_display(txl->dpix_out_pl);
		DRW_transform_to_display(txl->color);
}
