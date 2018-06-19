#include "NUL4.h"
#include "NUL_Util.h"
#include "NUL_TNS.h"
#include "tinycthread.h"
#include "png.h"

#include <math.h>

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/

void lanpr_make_render_buffers_recursive(tns3DObject* o,real* MVMat,real* MVPMat,tnsRenderBuffer* rb, real HeightMultiply){
	tns3DObject* oc;
	tnsMeshObject* mo;
	tnsVert* v;
	tnsFace* f;
	tnsEdge* e;
	tnsRenderTriangle* rt;
	tnsMatrix44d NewMVP;
	tnsMatrix44d NewMV;
	tnsFrameBuffer* fb = rb->FrameBuffer;
	tnsRenderElementLinkNode* reln;
	tnsCamera* c = rb->Scene->ActiveCamera;
	tnsMaterial* m;

	//if (o->RenderTriangles) FreeMem(o->RenderTriangles);
	//if (o->RenderVertices) FreeMem(o->RenderVertices);

	tMatMultiply44d(NewMVP, MVPMat, o->SelfTransform);
	tMatMultiply44d(NewMV, MVMat, o->SelfTransform);

	if (o->Type == TNS_OBJECT_MESH) {
		mo = o;
		o->RenderVertices = CreateNewBuffer(tnsRenderVert, mo->numV);
		o->RenderTriangles = calloc(mo->TriangleCount, rb->TriangleSize);//CreateNewBuffer(tnsRenderTriangle, mo->TriangleCount);
		//o->RenderLines = CreateNewBuffer(tnsRenderLine, mo->TriangulatedEdgeCount);

		reln = lstAppendPointerStaticSized(&rb->VertexBufferPointers, &rb->RenderDataPool, o->RenderVertices,
			sizeof(tnsRenderElementLinkNode));
		reln->ElementCount = mo->numV;
		reln->ObjectRef = mo;

		reln = lstAppendPointerStaticSized(&rb->TriangleBufferPointers, &rb->RenderDataPool, o->RenderTriangles,
			sizeof(tnsRenderElementLinkNode));
		reln->ElementCount = mo->TriangleCount;
		reln->ObjectRef = mo;

		for (v = mo->V.pFirst; v; v = v->Item.pNext) {
			tRdrTransformRenderVert(v, o->RenderVertices, NewMV, NewMVP,c);//,HeightMultiply,c->ZMin,c->ZMax);
			tObjRecalculateVertNormal(v);
		}

		for (e = mo->E.pFirst; e; e = e->Item.pNext) {
			tnsRenderLine* rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
			rl->L = e->VL->RV;
			rl->R = e->VR->RV;
			e->RenderLine = rl;
			tnsRenderLineSegment* rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
			lstAppendItem(&rl->Segments, rls);
			lstAppendItem(&rb->AllRenderLines, rl);
		}


		rt = o->RenderTriangles;
		for (f = mo->F.pFirst; f; f = f->Item.pNext) {
			tObjRecalculateFaceAverageNormal(f);
			tObjSimpleTriangulateRender(f, rt, rb->TriangleSize, o->RenderVertices, &rb->AllRenderLines, &rt, &rb->RenderDataPool);
			// already done in the func above. tRdrCalculateRenderTriangleNormal(rt);
			tMatApplyNormalTransform43d(f->GNormal, MVMat, f->FaceNormal);
			tMatNormalizeSelf3d(f->GNormal);
			m = tnsGetIndexedMaterial(rb->Scene, f->MaterialID);
			//if(m) m->PreviewVCount += (f->TriangleCount*3);
		}
	}

	for (oc = o->ChildObjects.pFirst; oc; oc = oc->Item.pNext) {
		tRdrMakeRenderGeometryBuffersRecursive(oc, NewMV, NewMVP, rb, HeightMultiply);
	}
}
void lanpr_make_render_buffers(tnsScene* s,tnsCamera* c, tnsRenderBuffer* rb, int HeightMultiply) {
	tns3DObject* o;
	tnsMatrixStackItem* tmsi = tKnlGetCurrentMatStackItem();
	tnsShader* current_shader = 0;
	tnsMatrix44d proj,view,result, inv;
	tnsFrameBuffer* fb = rb->FrameBuffer;

	if (!c) return;

	memset(rb->Scene->MaterialPointers, 0, sizeof(void*) * 2048);

	real asp = ((real)fb->W / (real)fb->H);

	if (c->CameraType == TNS_CAMERA_PERSPECTIVE) {
		tMatMakePerspectiveMatrix44d(proj, c->FOV, asp, c->ZMin, c->ZMax);
	}elif(c->CameraType == TNS_CAMERA_ORTHO) {
		real w = c->OrthScale;
		tMatMakeOrthographicMatrix44d(proj, -w, w, -w / asp, w / asp, c->ZMin, c->ZMax);
	}

	tMatLoadIdentity44d(view);

	tObjApplySelfTransformMatrix(c, 0);
	tObjApplyGlobalTransformMatrixReverted(c);
	tMatInverse44d(inv, c->Base.GlobalTransform);
	tMatMultiply44d(result, proj, inv);
	memcpy(proj, result, sizeof(tnsMatrix44d));
	memcpy(fb->ViewProjection, proj, sizeof(tnsMatrix44d));

	tMatInverse44d(fb->VPInverse, fb->ViewProjection);

	void* a;
	while (a = lstPopPointer(&rb->TriangleBufferPointers)) FreeMem(a);
	while (a = lstPopPointer(&rb->VertexBufferPointers)) FreeMem(a);

	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		//tObjApplyGlobalTransformMatrixRecursive(o);
		lanpr_make_render_buffers_recursive(o, view, proj, rb, (real)HeightMultiply);
	}
}

static int lanpr_compute_feature_lines_exec(bContext *C, wmOperator *op, const wmEvent *event){
	tnsRenderBuffer* rb = a->This->EndInstance;
	thrd_t SeperateThread;

	if (!rb || !rb->Scene || !rb->Scene->ActiveCamera) {
		return OPERATOR_CANCELLED;
	}

	rb->TriangleSize = tRdrGetRenderTriangleSize(rb);

	thrd_create(&SeperateThread, THREAD_ComputeFeatureLines, rb);

	return OPERATOR_FINISHED;
}

static void lanpr_compute_feature_lines_cancel(bContext *C, wmOperator *op){
    return;
}


void SCENE_OT_lanpr_calculate_feature_lines(struct wmOperatorType* ot){

	/* identifiers */
	ot->name = "Calculate Feature Lines";
	ot->description = "LANPR calculates feature line in current scene";
	ot->idname = "SCENE_OT_lanpr_calculate";

	/* api callbacks */
	//ot->invoke = screen_render_invoke; /* why we need both invoke and exec? */
	//ot->modal = screen_render_modal;
	ot->cancel = lanpr_compute_feature_lines_cancel;
	ot->exec = lanpr_compute_feature_lines;
}