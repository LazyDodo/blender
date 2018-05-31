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

typedef struct LANPROneTimeInit{
    
	/* Snake */

	GPUShader* multichannel_shader;
	GPUShader* edge_detect_shader;
	GPUShader* edge_thinning_shader;
	GPUShader* snake_connection_shader;

	/* DPIX */

	GPUShader* dpix_transform_shader;

	void* ved;
} LANPROneTimeInit;

//==============================================================[ ATLAS / DPIX ]


void lanpr_init_atlas_inputs(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;
	//LANPR_ViewLayerData *sldata = EEVEE_view_layer_data_ensure();
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

	txl->dpix_in_pl = 

	const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;


	/* Main Buffer */
	DRW_texture_ensure_2D(&txl->dpix_in_pl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);
	DRW_texture_ensure_2D(&txl->dpix_in_pr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);
	DRW_texture_ensure_2D(&txl->dpix_in_nl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);
	DRW_texture_ensure_2D(&txl->dpix_in_nr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);

	DRW_texture_ensure_2D(&txl->dpix_out_pl, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);
	DRW_texture_ensure_2D(&txl->dpix_out_pr, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);
	DRW_texture_ensure_2D(&txl->dpix_out_length, TNS_DPIX_TEXTURE_SIZE, TNS_DPIX_TEXTURE_SIZE, GPU_RGB32F, DRW_TEX_FILTER);

	GPU_framebuffer_ensure_config(&fbl->dpix_transform, {
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_TEXTURE(txl->dpix_out_pl),
		GPU_ATTACHMENT_TEXTURE(txl->dpix_out_pr),
		GPU_ATTACHMENT_TEXTURE(txl->dpix_out_length),
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE,
		GPU_ATTACHMENT_LEAVE
	});

	if (!OneTime.dpix_transform_shader) {
	OneTime.dpix_transform_shader = 
		GPU_shader_create(
			datatoc_gpu_shader_3D_normal_smooth_color_vert_glsl,
			datatoc_lanpr_snake_multichannel_fragment,NULL,NULL,NULL);
    }
}
void lanpr_destroy_atlas(void *ved){
	OneTime.ved = ved;
	LANPR_Data *vedata = (LANPR_Data *)ved;
	LANPR_TextureList *txl = vedata->txl;
	LANPR_FramebufferList *fbl = vedata->fbl;
	LANPR_StorageList *stl = ((LANPR_Data *)vedata)->stl;

	DRW_pass_free(psl->dpix_transform_pass);

	GPU_framebuffer_free(fbl->dpix_transform);

	DRW_texture_free(txl->dpix_in_pl);
	DRW_texture_free(txl->dpix_in_pr);
	DRW_texture_free(txl->dpix_in_nl);
	DRW_texture_free(txl->dpix_in_nr);
	DRW_texture_free(txl->dpix_out_pl);
	DRW_texture_free(txl->dpix_out_pr);
}

static Gwn_VertBuf *lanpr_mesh_get_tri_pos_and_normals_raw(
        MeshRenderData *rdata, const bool use_hide,
        Gwn_VertBuf **r_vbo)
{
	BLI_assert(rdata->types & (MR_DATATYPE_VERT | MR_DATATYPE_LOOPTRI | MR_DATATYPE_LOOP | MR_DATATYPE_POLY));

	if (*r_vbo == NULL) {
		static Gwn_VertFormat format = { 0 };
		static struct { uint pos, nor; } attr_id;
		if (format.attrib_ct == 0) {
			attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
			attr_id.nor = GWN_vertformat_attr_add(&format, "nor", GWN_COMP_I10, 3, GWN_FETCH_INT_TO_FLOAT_UNIT);
		}

		const int tri_len = mesh_render_data_looptri_len_get(rdata);

		Gwn_VertBuf *vbo = *r_vbo = GWN_vertbuf_create_with_format(&format);

		const int vbo_len_capacity = tri_len * 3;
		int vbo_len_used = 0;
		GWN_vertbuf_data_alloc(vbo, vbo_len_capacity);

		Gwn_VertBufRaw pos_step, nor_step;
		GWN_vertbuf_attr_get_raw_data(vbo, attr_id.pos, &pos_step);
		GWN_vertbuf_attr_get_raw_data(vbo, attr_id.nor, &nor_step);

		float (*lnors)[3] = rdata->loop_normals;

		if (rdata->edit_bmesh) {
			Gwn_PackedNormal *pnors_pack, *vnors_pack;

			if (lnors == NULL) {
				mesh_render_data_ensure_poly_normals_pack(rdata);
				mesh_render_data_ensure_vert_normals_pack(rdata);

				pnors_pack = rdata->poly_normals_pack;
				vnors_pack = rdata->vert_normals_pack;
			}

			for (int i = 0; i < tri_len; i++) {
				const BMLoop **bm_looptri = (const BMLoop **)rdata->edit_bmesh->looptris[i];
				const BMFace *bm_face = bm_looptri[0]->f;

				/* use_hide always for edit-mode */
				if (BM_elem_flag_test(bm_face, BM_ELEM_HIDDEN)) {
					continue;
				}

				if (lnors) {
					for (uint t = 0; t < 3; t++) {
						const float *nor = lnors[BM_elem_index_get(bm_looptri[t])];
						*((Gwn_PackedNormal *)GWN_vertbuf_raw_step(&nor_step)) = GWN_normal_convert_i10_v3(nor);
					}
				}
				else if (BM_elem_flag_test(bm_face, BM_ELEM_SMOOTH)) {
					for (uint t = 0; t < 3; t++) {
						*((Gwn_PackedNormal *)GWN_vertbuf_raw_step(&nor_step)) = vnors_pack[BM_elem_index_get(bm_looptri[t]->v)];
					}
				}
				else {
					const Gwn_PackedNormal *snor_pack = &pnors_pack[BM_elem_index_get(bm_face)];
					for (uint t = 0; t < 3; t++) {
						*((Gwn_PackedNormal *)GWN_vertbuf_raw_step(&nor_step)) = *snor_pack;
					}
				}

				/* TODO(sybren): deduplicate this and all the other places it's pasted to in this file. */
				if (rdata->edit_data && rdata->edit_data->vertexCos) {
					for (uint t = 0; t < 3; t++) {
						int vidx = BM_elem_index_get(bm_looptri[t]->v);
						const float *pos = rdata->edit_data->vertexCos[vidx];
						copy_v3_v3(GWN_vertbuf_raw_step(&pos_step), pos);
					}
				}
				else {
					for (uint t = 0; t < 3; t++) {
						copy_v3_v3(GWN_vertbuf_raw_step(&pos_step), bm_looptri[t]->v->co);
					}
				}
			}
		}
		else {
			if (lnors == NULL) {
				/* Use normals from vertex. */
				mesh_render_data_ensure_poly_normals_pack(rdata);
			}

			for (int i = 0; i < tri_len; i++) {
				const MLoopTri *mlt = &rdata->mlooptri[i];
				const MPoly *mp = &rdata->mpoly[mlt->poly];

				if (use_hide && (mp->flag & ME_HIDE)) {
					continue;
				}

				const uint vtri[3] = {
					rdata->mloop[mlt->tri[0]].v,
					rdata->mloop[mlt->tri[1]].v,
					rdata->mloop[mlt->tri[2]].v,
				};

				if (lnors) {
					for (uint t = 0; t < 3; t++) {
						const float *nor = lnors[mlt->tri[t]];
						*((Gwn_PackedNormal *)GWN_vertbuf_raw_step(&nor_step)) = GWN_normal_convert_i10_v3(nor);
					}
				}
				else if (mp->flag & ME_SMOOTH) {
					for (uint t = 0; t < 3; t++) {
						const MVert *mv = &rdata->mvert[vtri[t]];
						*((Gwn_PackedNormal *)GWN_vertbuf_raw_step(&nor_step)) = GWN_normal_convert_i10_s3(mv->no);
					}
				}
				else {
					const Gwn_PackedNormal *pnors_pack = &rdata->poly_normals_pack[mlt->poly];
					for (uint t = 0; t < 3; t++) {
						*((Gwn_PackedNormal *)GWN_vertbuf_raw_step(&nor_step)) = *pnors_pack;
					}
				}

				for (uint t = 0; t < 3; t++) {
					const MVert *mv = &rdata->mvert[vtri[t]];
					copy_v3_v3(GWN_vertbuf_raw_step(&pos_step), mv->co);
				}
			}
		}

		vbo_len_used = GWN_vertbuf_raw_used(&pos_step);
		BLI_assert(vbo_len_used == GWN_vertbuf_raw_used(&nor_step));

		if (vbo_len_capacity != vbo_len_used) {
			GWN_vertbuf_data_resize(vbo, vbo_len_used);
		}
	}
	return *r_vbo;
}

void DRW_cache_mesh_surface_get(Object *ob)
{
	BLI_assert(ob->type == OB_MESH);

	Mesh *me = ob->data;	
		
	MeshBatchCache *cache = mesh_batch_cache_get(me);

	if (cache->triangles_with_normals == NULL) {
		const int datatype = MR_DATATYPE_VERT | MR_DATATYPE_LOOPTRI | MR_DATATYPE_LOOP | MR_DATATYPE_POLY;
		MeshRenderData *rdata = mesh_render_data_create(me, datatype);

		cache->triangles_with_normals = GWN_batch_create(
		        GWN_PRIM_TRIS, mesh_batch_cache_get_tri_pos_and_normals(rdata, cache), NULL);

		mesh_render_data_free(rdata);
	}

}

int lanpr_feed_atlas_data_obj(
	float* AtlasPointsL, float* AtlasPointsR,
	float* AtlasFaceNormalL, float* AtlasFaceNormalR,
	tns3DObject* o, int Begin) {
	
	tnsVert* v;
	tnsEdge* e;
	tnsFace* f;
	tns3DObject* io;
	tnsMeshObject* mo;
	int NextBegin = Begin;
	int ThisCount=0;

	for (io = o->ChildObjects.pFirst; io; io = io->Item.pNext) {
		NextBegin = tns_FeedAtlasDataRecursive(
			AtlasPointsL, AtlasPointsR, AtlasFaceNormalL, AtlasFaceNormalR, io, NextBegin);
	}

	if (o->Type != TNS_OBJECT_MESH) return Begin;

	tObjRefreshGeometeryIndex(o);
	mo = o;

	for (e = mo->E.pFirst; e; e = e->Item.pNext) {
		int offset = ThisCount + NextBegin;
		AtlasPointsL[offset + 0] = e->VL->P[0];
		AtlasPointsL[offset + 1] = e->VL->P[1];
		AtlasPointsL[offset + 2] = e->VL->P[2];
		AtlasPointsR[offset + 0] = e->VR->P[0];
		AtlasPointsR[offset + 1] = e->VR->P[1];
		AtlasPointsR[offset + 2] = e->VR->P[2];
		AtlasFaceNormalL[offset + 0] = e->FL ? e->FL->FaceNormal[0] : 0;
		AtlasFaceNormalL[offset + 1] = e->FL ? e->FL->FaceNormal[1] : 0;
		AtlasFaceNormalL[offset + 2] = e->FL ? e->FL->FaceNormal[2] : 0;
		AtlasFaceNormalR[offset + 0] = e->FR ? e->FR->FaceNormal[0] : 0;
		AtlasFaceNormalR[offset + 1] = e->FR ? e->FR->FaceNormal[1] : 0;
		AtlasFaceNormalR[offset + 2] = e->FR ? e->FR->FaceNormal[2] : 0;
		ThisCount += 3;
	}

	return ThisCount + NextBegin;

}
void lanpr_feed_atlas_data(tnsScene* s) {
	tns3DObject* o;
	int NextBegin = 0;
	float* AtlasPointsL,*AtlasPointsR,
	       *AtlasFaceNormalL,*AtlasFaceNormalR;


	if (!s || !s->AtlasPointsL) return;

	int width = s->AtlasPointsL->Width;
	int count = width*width*3;

	AtlasPointsL = CreateNewBuffer(float, count);
	AtlasPointsR = CreateNewBuffer(float, count);
	AtlasFaceNormalL = CreateNewBuffer(float, count);
	AtlasFaceNormalR = CreateNewBuffer(float, count);

	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		NextBegin = tns_FeedAtlasDataRecursive(
			AtlasPointsL, AtlasPointsR, AtlasFaceNormalL, AtlasFaceNormalR, o, NextBegin);
	}

	//NextBegin = AtlasCount;

	tnsFeed2DRectTextureData(s->AtlasPointsL, AtlasPointsL);
	tnsFeed2DRectTextureData(s->AtlasPointsR, AtlasPointsR);
	tnsFeed2DRectTextureData(s->AtlasFaceNormalL, AtlasFaceNormalL);
	tnsFeed2DRectTextureData(s->AtlasFaceNormalR, AtlasFaceNormalR);

	FreeMem(AtlasPointsL);
	FreeMem(AtlasPointsR);
	FreeMem(AtlasFaceNormalL);
	FreeMem(AtlasFaceNormalR);
}

void lanpr_get_atlas_buffer_position(int index, float* x, float* y) {
	int px, py;
	px = index % TNS_ATLAS_DEFAULT_INPUT_WIDTH;
	py = index / TNS_ATLAS_DEFAULT_INPUT_WIDTH;

	*x = (float)px / TNS_ATLAS_DEFAULT_INPUT_WIDTH * 2 - 1;
	*y = (float)py / TNS_ATLAS_DEFAULT_INPUT_WIDTH * 2 - 1;
}

int lanpr_make_trigger_batch_recursive(int Begin, tns3DObject* o) {
	tnsMeshObject* mo = o,*io;
	int i;
	int Next = Begin;

	for (io = o->ChildObjects.pFirst; io; io = io->Base.Item.pNext) {
		Next = tns_MakeAtlasTriggerBatchRecursive(Next, io);
	}

	if (o->Type != TNS_OBJECT_MESH) return Begin;

	float* data = CreateNewBuffer(float, mo->numE*2);

	for (i = 0; i < mo->numE; i++) {
		tns_GetAtlasBufferPosition(Begin+i, &data[i*2], &data[i*2 + 1]);
	}

	mo->AtlasTriggerBatch = tnsCreateBatch(mo->numE, 2, data, 0);
	//mo->AtlasTriggerBatch->BeginElementOffset = Begin;

	Next = Begin + mo->numE;

	

	return Next;
}

void lanpr_calculate_view_direction(tnsScene* s, tnsCamera* c, tnsVector3d Vec) {
	tnsVector3d Direction = { 0,0,-1 };
	tnsVector3d Trans;
	tnsMatrix44d inv;
	tMatInverse44d(inv, c->Base.GlobalTransform);
	tMatApplyRotation43d(Trans, inv, Direction);
	//tMatVectorCopy3d(Trans, Vec);
	//tMatVectorMultiSelf3d(Trans, -1);
	tMatVectorCopy3d(Trans, Vec);

}

void lanpr_trigger_this_object(tns3DObject* o) {
	tnsMeshObject *mo = o;
	tnsShader* cs = T->AtlasTransformShader;

	if (o->Type != TNS_OBJECT_MESH || !mo->AtlasTriggerBatch) return;

	glBindBuffer(GL_ARRAY_BUFFER, mo->AtlasTriggerBatch->VBO);
	glEnableVertexAttribArray(cs->vertexIndex);
	glVertexAttribPointer(cs->vertexIndex, mo->AtlasTriggerBatch->Dimention, GL_FLOAT, 0, 0, 0);

	glPointSize(1);
	glDrawArrays(GL_POINTS, 0,
		//mo->AtlasTriggerBatch->BeginElementOffset*mo->AtlasTriggerBatch->Dimention,
		mo->AtlasTriggerBatch->NumVert*mo->AtlasTriggerBatch->Dimention);
}

void lanpr_trigger_atlas_transform_object(tns3DObject* o) {

	tnsMeshObject* mo = o, *io;
	int i;
	
	for (io = o->ChildObjects.pFirst; io; io = io->Base.Item.pNext) {
		tnsPushMatrix();
		tnsApplyObjectMatrix(io);
		tns_TriggerThisObject(io);
		if (o->ChildObjects.pFirst) {
			tns_TriggerAtlasTransformRecursive(io);
		}
		tnsPopMatrix();
	}
}
void lanpr_trigger_atlas_transform(tnsScene* s, n3DViewUiExtra* e) {
	tns3DObject* o;
	tnsShader* sd = T->AtlasTransformShader;
	//if (!s) return;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		tnsPushMatrix();
		tnsApplyObjectMatrix(o);
		tns_TriggerThisObject(o);
		if (o->ChildObjects.pFirst) {
			tns_TriggerAtlasTransformRecursive(o);
		}
		tnsPopMatrix();
	}

}

extern const GLuint TNS_ATTACHMENT_ARRAY_0_1_2[];

void lanpr_make_atlas_trigger_batch(tnsScene* s) {
	tns3DObject* o;
	int Begin = 0;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		Begin = tns_MakeAtlasTriggerBatchRecursive(Begin, o);
	}
}

void lanpr_atlas_draw_transform(n3DViewUiExtra* e) {
	tns3DObject* o; tnsScene* s = e->CurrentScene;
	tnsShader* sd = T->AtlasTransformShader;
	if (!s) return;
	int W = e->AtlasPointsOutL->Width, H = e->AtlasPointsOutL->Height;

	tnsUseShader(sd);
	tnsEnableShaderv(sd);

	tnsDrawToOffscreen(e->AtlasFBO, 3, TNS_ATTACHMENT_ARRAY_0_1_2);
	tnsViewportWithScissor(0, 0,W, H);
	//tnsOrtho(0, W, H, 0, -100, 100);

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);

	tnsActiveTexture(GL_TEXTURE0);
	tnsBindTexture0(s->AtlasPointsL);
	tnsActiveTexture(GL_TEXTURE1);
	tnsBindTexture1(s->AtlasPointsR);
	tnsActiveTexture(GL_TEXTURE2);
	tnsBindTexture2(s->AtlasFaceNormalL);
	tnsActiveTexture(GL_TEXTURE3);
	tnsBindTexture3(s->AtlasFaceNormalR);

	tnsVector3d Viewdir;
	real* camloc = e->ViewingCamera->Base.GLocation;

	tnsApplyCameraView(W, H, e->ViewingCamera);

	tns_CalculateViewDir(s, e->ViewingCamera, Viewdir);

	glUniform4f(sd->uniform0Index, 0, 0, e->OffScr->pColorTextures[0]->Width, e->OffScr->pColorTextures[0]->Height);
	glUniform3f(sd->uniform4Index, (GLfloat)camloc[0], (GLfloat)camloc[1], (GLfloat)camloc[2]);
	glUniform3f(sd->uniform5Index, (GLfloat)Viewdir[0], (GLfloat)Viewdir[1], (GLfloat)Viewdir[2]);
	glUniform1f(sd->uniform6Index, 1);
	glUniform1f(sd->uniform7Index, TNS_ATLAS_DEFAULT_INPUT_WIDTH);

	
	//tnsVertex2d(10, 10);
	//tnsVertex2d(10, 90);
	//tnsVertex2d(90, 90);
	//tnsVertex2d(90, 10);
	//tnsPackAs(GL_POINTS);
	//tnsFlush();

	tns_TriggerAtlasTransform(s, e);
}

void lanpr_atlas_draw_edge_preview(n3DViewUiExtra* e) {
	tnsEnableShaderv(T->AtlasPreviewShader);
	glDisableVertexAttribArray(T->AtlasPreviewShader->colorIndex);
	glVertexAttrib4f(T->AtlasPreviewShader->colorIndex, 1, 1, 1, 1);

	tnsActiveTexture(GL_TEXTURE0);
	tnsBindTexture0(e->AtlasPointsOutL);

	tnsActiveTexture(GL_TEXTURE1);
	tnsBindTexture1(e->AtlasPointsOutR);

	//tnsUseUiShader();
	//tnsUseShader(T->AtlasPreviewShader);

	//tnsColor4d(1, 1, 1, 1);
	//tnsVertex2d(0.10, 0.10);
	//tnsVertex2d(0.90, 0.10);
	//tnsVertex2d(0.90, 0.90);
	//tnsVertex2d(0.10, 0.90);
	//tnsPackAs(GL_LINE_LOOP);
	//tnsFlush();

	glUniform2f(T->AtlasPreviewShader->uniform1Index, e->OffScr->pColorTextures[0]->Width, e->OffScr->pColorTextures[0]->Height);
	glUniform1f(T->AtlasPreviewShader->uniform0Index, TNS_ATLAS_DEFAULT_INPUT_WIDTH);

	glLineWidth(2);
	glEnable(GL_LINE_SMOOTH);

	tns_TriggerAtlasTransform(e->CurrentScene, e);

	glLineWidth(1);
	glDisable(GL_LINE_SMOOTH);

}


//=====================================================================[ SNAKE ]


//==============================================[ ENGINE ]


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

	psl->dpix_transform_pass = DRW_pass_create("DPIX Transform Stage", DRW_STATE_WRITE_COLOR);
	stl->g_data->dpix_transform_shgrp = DRW_shgroup_create(OneTime.dpix_transform_shader, psl->dpix_transform_pass);
	//DRW_shgroup_uniform_texture_ref(stl->g_data->dpix_transform_shgrp, "TexSample0", &txl->color);
	//DRW_shgroup_uniform_int(stl->g_data->dpix_transform_shgrp, "Stage", &stl->g_data->stage, 1);
	//DRW_shgroup_call_add(stl->g_data->dpix_transform_shgrp, quad, NULL);
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

void lanpr_count_drawing_elements(LANPR_PrivateData* pd, int* vert_count, int* index_adjacent_count){
    int v_count=0;
	int e_count=0;
	LANPR_LineStrip* ls;
	for (ls = (LANPR_LineStrip *)(pd->line_strips.first); ls; ls = (LANPR_LineStrip *)(ls->Item.next)) {
		v_count += (ls->point_count);
		e_count += ((ls->point_count - 1) * 4);
	}
	*vert_count = v_count;
	*index_adjacent_count = e_count;
}

Gwn_Batch *lanpr_get_snake_batch(LANPR_PrivateData* pd){
	LANPR_LineStrip* ls;
	LANPR_LineStripPoint* lsp, *plsp;
	int i;
	u32bit *Index_adjacent;
	float* Verts;
	float* Lengths;
	float TotalLength=0;
	int v_count,e_count;

	lanpr_count_drawing_elements(pd,&v_count,&e_count);

	//Index_adjacent = MEM_callocN(sizeof(unsigned int) * e_count, "Index_adjacent buffer pre alloc");
	Verts = MEM_callocN(sizeof(float) * v_count * 2, "Verts buffer pre alloc");
	Lengths = MEM_callocN(sizeof(float)* v_count, "Length buffer pre alloc");

	Gwn_IndexBufBuilder elb;
	GWN_indexbuf_init_ex(&elb, GWN_PRIM_LINES_ADJ, e_count, v_count, true);

	int vert_offset=0;

	for (ls = (LANPR_LineStrip *)(pd->line_strips.first); ls; ls = (LANPR_LineStrip *)(ls->Item.next)) {
		for (i = 0; i < ls->point_count-1; i++) {
			int v1 = i+vert_offset-1;
			int v2 = i+vert_offset;
			int v3 = i+vert_offset+1;
			int v4 = i+vert_offset+2;
			if(v1<0) v1=0;
			if(v4>= v_count) v4 = v_count -1;
			GWN_indexbuf_add_line_adj_verts(&elb, v1,v2,v3,v4);
		}

		i = 0;
		float xf,yf;
		for (lsp = (LANPR_LineStripPoint *)(ls->points.first); lsp; lsp = (LANPR_LineStripPoint *)(lsp->Item.next)) {
			lanpr_texture_to_ndc(lsp->P[0],lsp->P[1],pd->width, pd->height, &xf,&yf);
			Verts[vert_offset*2 + i * 2 + 0] = xf;
			Verts[vert_offset*2 + i * 2 + 1] = yf;
			if (plsp = (LANPR_LineStripPoint *)(lsp->Item.prev)) {
				TotalLength += tMatDist2v(plsp->P, lsp->P);
				Lengths[vert_offset + i] = TotalLength;
			}
			i++;
		}
		ls->total_length = TotalLength;

		vert_offset+=(ls->point_count);
	}

	static Gwn_VertFormat format = { 0 };
	static struct { uint pos, uvs; } attr_id;
	if (format.attrib_ct == 0) {
		attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
		attr_id.uvs = GWN_vertformat_attr_add(&format, "uvs", GWN_COMP_F32, 1, GWN_FETCH_FLOAT);
	}

	Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);
	GWN_vertbuf_data_alloc(vbo, v_count);

	for (int i = 0; i < v_count; ++i) {
		GWN_vertbuf_attr_set(vbo, attr_id.pos, i, &Verts[i*2]);
		GWN_vertbuf_attr_set(vbo, attr_id.uvs, i, &Lengths[i]);
	}

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
	stl->g_data->normal_clamp =    lanpr->normal_clamp;
	stl->g_data->normal_strength = lanpr->normal_strength;
	stl->g_data->depth_clamp =     lanpr->depth_clamp;
	stl->g_data->depth_strength =  lanpr->depth_strength;
	//GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
	GPU_framebuffer_bind(fbl->edge_intermediate);
	DRW_draw_pass(psl->edge_intermediate);

	if((!lanpr->enable_vector_trace) && (!lanpr->display_thinning_result)){
		GPU_framebuffer_bind(dfbl->default_fb);
		DRW_transform_to_display(txl->edge_intermediate);
		return;
	}

	if(lanpr->display_thinning_result || lanpr->enable_vector_trace){
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

		if(!lanpr->enable_vector_trace){
			GPU_framebuffer_bind(dfbl->default_fb);
			DRW_transform_to_display(txl->edge_intermediate);
			return;
		}
	}
    

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

	float* tld = &lanpr->taper_left_distance, *tls = &lanpr->taper_left_strength,
		 *trd = &lanpr->taper_right_distance, *trs = &lanpr->taper_right_strength;

	Gwn_Batch* snake_batch = lanpr_get_snake_batch(stl->g_data);

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