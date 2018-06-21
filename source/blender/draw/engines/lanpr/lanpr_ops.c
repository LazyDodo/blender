#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
#include "lanpr_all.h"
#include "lanpr_util.h"
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

#include "WM_types.h"
#include "WM_api.h"

#include <math.h>

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/

struct Object;

/* ====================================== base structures =========================================== */

#define TNS_BOUND_AREA_CROSSES(b1,b2)\
((b1)[0]<(b2)[1] && (b1)[1]>(b2)[0] && (b1)[3]<(b2)[2] && (b1)[2]>(b2)[3])

void lanpr_MakeInitialBoundingAreas(LANPR_RenderBuffer* rb) {
	int SpW = 20;
	int SpH = rb->H / (rb->W / SpW);
	int Row, Col;
	LANPR_BoundingArea* ba;
	real W = (real)rb->W;
	real H = (real)rb->H;
	real SpanW = (real)1 / SpW * 2.0;
	real SpanH = (real)1 / SpH * 2.0;

	rb->TileCountX = SpW;
	rb->TileCountY = SpH;
	rb->WidthPerTile = SpanW;
	rb->HeightPerTile = SpanH;

	rb->BoundingAreaCount = SpW * SpH;
	rb->InitialBoundingAreas = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_BoundingArea) * rb->BoundingAreaCount);

	for (Row = 0; Row < SpH; Row++) {
		for (Col = 0; Col < SpW; Col++) {
			ba = &rb->InitialBoundingAreas[Row * 20 + Col];

			ba->L = SpanW * Col - 1.0;
			ba->R = (Col == SpW - 1) ? 1.0 : (SpanW * (Col + 1) - 1.0);
			ba->U = 1.0 - SpanH * Row;
			ba->B = (Row == SpH - 1) ? -1.0 : (1.0 - SpanH * (Row + 1));

			ba->CX = (ba->L + ba->R) / 2;
			ba->CY = (ba->U + ba->B) / 2;

			if (Row) {
				lstAppendPointerStatic(&ba->UP, &rb->RenderDataPool, &rb->InitialBoundingAreas[(Row - 1) * 20 + Col]);
			}
			if (Col) {
				lstAppendPointerStatic(&ba->LP, &rb->RenderDataPool, &rb->InitialBoundingAreas[Row * 20 + Col - 1]);
			}
			if (Row != SpH -1) {
				lstAppendPointerStatic(&ba->BP, &rb->RenderDataPool, &rb->InitialBoundingAreas[(Row + 1) * 20 + Col]);
			}
			if (Col != SpW - 1) {
				lstAppendPointerStatic(&ba->RP, &rb->RenderDataPool, &rb->InitialBoundingAreas[Row * 20 + Col + 1]);
			}
		}
	}
}
void lanpr_ConnectNewBoundingAreas(LANPR_RenderBuffer* rb, LANPR_BoundingArea* Root) {
	LANPR_BoundingArea* ba = Root->Child, *tba;
	nListItemPointer* lip,*lip2,*lip3,*NextLip;
	nStaticMemoryPool* mph= &rb->RenderDataPool;

	lstAppendPointerStaticPool(mph, &ba[1].RP, &ba[0]);
	lstAppendPointerStaticPool(mph, &ba[0].LP, &ba[1]);
	lstAppendPointerStaticPool(mph, &ba[1].BP, &ba[2]);
	lstAppendPointerStaticPool(mph, &ba[2].UP, &ba[1]);
	lstAppendPointerStaticPool(mph, &ba[2].RP, &ba[3]);
	lstAppendPointerStaticPool(mph, &ba[3].LP, &ba[2]);
	lstAppendPointerStaticPool(mph, &ba[3].UP, &ba[0]);
	lstAppendPointerStaticPool(mph, &ba[0].BP, &ba[3]);

	for (lip = Root->LP.pFirst; lip; lip = lip->pNext) {
		tba = lip->p;
		if (ba[1].U > tba->B && ba[1].B < tba->U) { lstAppendPointerStaticPool(mph, &ba[1].LP, tba); lstAppendPointerStaticPool(mph, &tba->RP, &ba[1]); }
		if (ba[2].U > tba->B && ba[2].B < tba->U) { lstAppendPointerStaticPool(mph, &ba[2].LP, tba); lstAppendPointerStaticPool(mph, &tba->RP, &ba[2]); }
	}
	for (lip = Root->RP.pFirst; lip; lip = lip->pNext) {
		tba = lip->p;
		if (ba[0].U > tba->B && ba[0].B < tba->U) { lstAppendPointerStaticPool(mph, &ba[0].RP, tba); lstAppendPointerStaticPool(mph, &tba->LP, &ba[0]); }
		if (ba[3].U > tba->B && ba[3].B < tba->U) { lstAppendPointerStaticPool(mph, &ba[3].RP, tba); lstAppendPointerStaticPool(mph, &tba->LP, &ba[3]); }
	}
	for (lip = Root->UP.pFirst; lip; lip = lip->pNext) {
		tba = lip->p;
		if (ba[0].R > tba->L && ba[0].L < tba->R) { lstAppendPointerStaticPool(mph, &ba[0].UP, tba); lstAppendPointerStaticPool(mph, &tba->BP, &ba[0]); }
		if (ba[1].R > tba->L && ba[1].L < tba->R) { lstAppendPointerStaticPool(mph, &ba[1].UP, tba); lstAppendPointerStaticPool(mph, &tba->BP, &ba[1]); }
	}
	for (lip = Root->BP.pFirst; lip; lip = lip->pNext) {
		tba = lip->p;
		if (ba[2].R > tba->L && ba[2].L < tba->R) { lstAppendPointerStaticPool(mph, &ba[2].BP, tba); lstAppendPointerStaticPool(mph, &tba->UP, &ba[2]); }
		if (ba[3].R > tba->L && ba[3].L < tba->R) { lstAppendPointerStaticPool(mph, &ba[3].BP, tba); lstAppendPointerStaticPool(mph, &tba->UP, &ba[3]); }
	}
	for (lip = Root->LP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((LANPR_BoundingArea*)lip->p)->RP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((LANPR_BoundingArea*)lip->p)->RP, lip2);
				if (ba[1].U > tba->B && ba[1].B < tba->U)lstAppendPointerStaticPool(mph, &tba->RP, &ba[1]);
				if (ba[2].U > tba->B && ba[2].B < tba->U)lstAppendPointerStaticPool(mph, &tba->RP, &ba[2]);
			}
		}
	}
	for (lip = Root->RP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((LANPR_BoundingArea*)lip->p)->LP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((LANPR_BoundingArea*)lip->p)->LP, lip2);
				if (ba[0].U > tba->B && ba[0].B < tba->U)lstAppendPointerStaticPool(mph, &tba->LP, &ba[0]);
				if (ba[3].U > tba->B && ba[3].B < tba->U)lstAppendPointerStaticPool(mph, &tba->LP, &ba[3]);
			}
		}
	}
	for (lip = Root->UP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((LANPR_BoundingArea*)lip->p)->BP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((LANPR_BoundingArea*)lip->p)->BP, lip2);
				if (ba[0].R > tba->L && ba[0].L < tba->R)lstAppendPointerStaticPool(mph, &tba->UP, &ba[0]);
				if (ba[1].R > tba->L && ba[1].L < tba->R)lstAppendPointerStaticPool(mph, &tba->UP, &ba[1]);
			}
		}
	}
	for (lip = Root->BP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((LANPR_BoundingArea*)lip->p)->UP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((LANPR_BoundingArea*)lip->p)->UP, lip2);
				if (ba[2].R > tba->L && ba[2].L < tba->R)lstAppendPointerStaticPool(mph, &tba->BP, &ba[2]);
				if (ba[3].R > tba->L && ba[3].L < tba->R)lstAppendPointerStaticPool(mph, &tba->BP, &ba[3]);
			}
		}
	}
	while (lstPopPointerNoFree(&Root->LP));
	while (lstPopPointerNoFree(&Root->RP));
	while (lstPopPointerNoFree(&Root->UP));
	while (lstPopPointerNoFree(&Root->BP));
}
void lanpr_AssociateTriangleWithBoundingArea(LANPR_RenderBuffer* rb, LANPR_BoundingArea* RootBoundingArea, LANPR_RenderTriangle* rt, real* LRUB, int Recursive);
int lanpr_TriangleCalculateIntersectionsInTile(LANPR_RenderBuffer* rb, LANPR_RenderTriangle* rt, LANPR_BoundingArea* ba);

void lanpr_SplitBoundingArea(LANPR_RenderBuffer* rb, LANPR_BoundingArea* Root) {
	LANPR_BoundingArea* ba = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_BoundingArea) * 4);
	LANPR_RenderTriangle* rt;

	ba[0].L = Root->CX;
	ba[0].R = Root->R;
	ba[0].U = Root->U;
	ba[0].B = Root->CY;
	ba[0].CX = (ba[0].L + ba[0].R) / 2;
	ba[0].CY = (ba[0].U + ba[0].B) / 2;

	ba[1].L = Root->L;
	ba[1].R = Root->CX;
	ba[1].U = Root->U;
	ba[1].B = Root->CY;
	ba[1].CX = (ba[1].L + ba[1].R) / 2;
	ba[1].CY = (ba[1].U + ba[1].B) / 2;

	ba[2].L = Root->L;
	ba[2].R = Root->CX;
	ba[2].U = Root->CY;
	ba[2].B = Root->B;
	ba[2].CX = (ba[2].L + ba[2].R) / 2;
	ba[2].CY = (ba[2].U + ba[2].B) / 2;

	ba[3].L = Root->CX;
	ba[3].R = Root->R;
	ba[3].U = Root->CY;
	ba[3].B = Root->B;
	ba[3].CX = (ba[3].L + ba[3].R) / 2;
	ba[3].CY = (ba[3].U + ba[3].B) / 2;

	Root->Child = ba;

	lanpr_ConnectNewBoundingAreas(rb, Root);

	while (rt = lstPopPointerNoFree(&Root->AssociatedTriangles)) {
		LANPR_BoundingArea* ba = Root->Child;
		real B[4];
		B[0] = TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
		B[1] = TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
		B[2] = TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
		B[3] = TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[0].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[0], rt, B, 0);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[1].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[1], rt, B, 0);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[2].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[2], rt, B, 0);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[3].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[3], rt, B, 0);
	}

	rb->BoundingAreaCount += 3;
}
int lanpr_LineCrossesBoundingArea(LANPR_RenderBuffer* fb, tnsVector2d L, tnsVector2d R, LANPR_BoundingArea* ba) {
	real vx, vy;
	tnsVector4d Converted;
	real c1, c;

	if ((Converted[0] = (real)ba->L) > TNS_MAX2(L[0], R[0])) return 0;
	if ((Converted[1] = (real)ba->R) < TNS_MIN2(L[0], R[0])) return 0;
	if ((Converted[2] = (real)ba->B) > TNS_MAX2(L[1], R[1])) return 0;
	if ((Converted[3] = (real)ba->U) < TNS_MIN2(L[1], R[1])) return 0;

	vx = L[0] - R[0];
	vy = L[1] - R[1];

	c1 = vx * (Converted[2] - L[1]) - vy * (Converted[0] - L[0]);
	c = c1;

	c1 = vx * (Converted[2] - L[1]) - vy * (Converted[1] - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	c1 = vx * (Converted[3] - L[1]) - vy * (Converted[0] - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	c1 = vx * (Converted[3] - L[1]) - vy * (Converted[1] - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	return 0;
}
int lanpr_TriangleCoversBoundingArea(LANPR_RenderBuffer* fb, LANPR_RenderTriangle* rt, LANPR_BoundingArea* ba) {

	real vx, vy;
	tnsVector2d p1, p2, p3, p4;
	real
		*FBC1 = rt->V[0]->FrameBufferCoord,
		*FBC2 = rt->V[1]->FrameBufferCoord,
		*FBC3 = rt->V[2]->FrameBufferCoord;

	p3[0] = p1[0] = (real)ba->L;
	p2[1] = p1[1] = (real)ba->B;
	p2[0] = p4[0] = (real)ba->R;
	p3[1] = p4[1] = (real)ba->U;

	if (FBC1[0] >= p1[0] && FBC1[0] <= p2[0] && FBC1[1] >= p1[1] && FBC1[1] <= p3[1]) return 1;
	if (FBC2[0] >= p1[0] && FBC2[0] <= p2[0] && FBC2[1] >= p1[1] && FBC2[1] <= p3[1]) return 1;
	if (FBC3[0] >= p1[0] && FBC3[0] <= p2[0] && FBC3[1] >= p1[1] && FBC3[1] <= p3[1]) return 1;

	if (lanpr_PointInsideTrianglef(p1, FBC1, FBC2, FBC3) ||
		lanpr_PointInsideTrianglef(p2, FBC1, FBC2, FBC3) ||
		lanpr_PointInsideTrianglef(p3, FBC1, FBC2, FBC3) ||
		lanpr_PointInsideTrianglef(p4, FBC1, FBC2, FBC3)) return 1;

	if  (lanpr_LineCrossesBoundingArea(fb, FBC1, FBC2, ba)) return 1;
	elif(lanpr_LineCrossesBoundingArea(fb, FBC2, FBC3, ba)) return 1;
	elif(lanpr_LineCrossesBoundingArea(fb, FBC3, FBC1, ba)) return 1;

	return 0;
}
void lanpr_AssociateTriangleWithBoundingArea(LANPR_RenderBuffer* rb, LANPR_BoundingArea* RootBoundingArea, LANPR_RenderTriangle* rt, real* LRUB, int Recursive) {
	if (!lanpr_TriangleCoversBoundingArea(rb, rt, RootBoundingArea)) return;
	if (!RootBoundingArea->Child) {
		lstAppendPointerStaticPool(&rb->RenderDataPool, &RootBoundingArea->AssociatedTriangles, rt);
		RootBoundingArea->TriangleCount++;
		if (RootBoundingArea->TriangleCount > 200 && Recursive) {
			lanpr_SplitBoundingArea(rb, RootBoundingArea);
		}
		if(Recursive) lanpr_TriangleCalculateIntersectionsInTile(rb, rt, RootBoundingArea);
	}else {
		LANPR_BoundingArea* ba = RootBoundingArea->Child;
		real* B1 = LRUB;
		real B[4];
		if (!LRUB) {
			B[0] = TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
			B[1] = TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
			B[2] = TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
			B[3] = TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
			B1 = B;
		}
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[0].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[0], rt, B1, Recursive);
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[1].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[1], rt, B1, Recursive);
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[2].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[2], rt, B1, Recursive);
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[3].L)) lanpr_AssociateTriangleWithBoundingArea(rb, &ba[3], rt, B1, Recursive);
	}
}
int lanpr_GetTriangleBoundingTile(LANPR_RenderBuffer* rb, LANPR_RenderTriangle* rt,int* RowBegin,int* RowEnd,int* ColBegin,int* ColEnd) {
	real SpW = rb->WidthPerTile, SpH = rb->HeightPerTile;
	real B[4];

	if (!rt->F || !rt->V[0] || !rt->V[1] || !rt->V[2]) return 0;

	B[0] = TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
	B[1] = TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
	B[2] = TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
	B[3] = TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);

	if (B[0] > 1 || B[1] < -1 || B[2]>1 || B[3] < -1) return 0;

	(*ColBegin) = (int)((B[0] + 1.0) / SpW);
	(*ColEnd) = (int)((B[1] + 1.0) / SpW);
	(*RowEnd) = rb->TileCountY - (int)((B[2] + 1.0) / SpH) - 1;
	(*RowBegin) = rb->TileCountY - (int)((B[3] + 1.0) / SpH) - 1;

	if ((*ColEnd) >= rb->TileCountX) (*ColEnd) = rb->TileCountX-1;
	if ((*RowEnd) >= rb->TileCountY) (*RowEnd) = rb->TileCountY-1;
	if ((*ColBegin) < 0) (*ColBegin) = 0;
	if ((*RowBegin) < 0) (*RowBegin) = 0;

	return 1;
}
void lanpr_AddTriangles(LANPR_RenderBuffer* rb) {
	LANPR_RenderElementLinkNode* reln;
	LANPR_RenderTriangle* rt;
	LANPR_RenderTile* tile;
	tnsMatrix44d VP;
	Camera* c = ((Camera*)rb->Scene->ActiveCamera);
	int i, lim;
	int x1, x2, y1, y2;
	int r, co;
	//tnsMatrix44d proj, view, result, inv;
	//tMatMakePerspectiveMatrix44d(proj, c->FOV, (real)fb->W / (real)fb->H, c->ZMin, c->ZMax);
	//tMatLoadIdentity44d(view);
	//tObjApplySelfTransformMatrix(c, 0);
	//tObjApplyGlobalTransformMatrixReverted(c);
	//tMatInverse44d(inv, c->Base.GlobalTransform);
	//tMatMultiply44d(result, proj, inv);
	//memcpy(proj, result, sizeof(tnsMatrix44d));

	//tnsglobal_TriangleIntersectionCount = 0;

	//tnsset_RenderOverallProgress(rb, NUL_MH2);
	rb->CalculationStatus = TNS_CALCULATION_INTERSECTION;
	//nulThreadNotifyUsers("tns.render_buffer_list.calculation_status");

	for (reln = rb->TriangleBufferPointers.pFirst; reln; reln = reln->Item.pNext) {
		rt = reln->Pointer;
		lim = reln->ElementCount;
		for (i = 0; i < lim; i++) {
			if (rt->CullStatus) {
				rt++; continue;
			}
			if (lanpr_GetTriangleBoundingTile(rb, rt, &y1, &y2, &x1, &x2)) {
				for (co = x1; co <= x2; co++) {
					for (r = y1; r <= y2; r++) {
						lanpr_AssociateTriangleWithBoundingArea(rb, &rb->InitialBoundingAreas[r * 20 + co], rt, 0, 1);
					}
				}
			}else{
				;//throw away.
			}
			rt = (void*)(((BYTE*)rt) + rb->TriangleSize);
			//if (tnsglobal_TriangleIntersectionCount >= 2000) {
				//tnsset_PlusRenderIntersectionCount(rb, tnsglobal_TriangleIntersectionCount);
				//tnsglobal_TriangleIntersectionCount = 0;
			//}
		}
	}
	//tnsset_PlusRenderIntersectionCount(rb, tnsglobal_TriangleIntersectionCount);
}
LANPR_BoundingArea* lanpr_GetNextBoundingArea(LANPR_BoundingArea* This, LANPR_RenderLine* rl, real x, real y, real k, int PositiveX, int PositiveY, real* NextX, real* NextY) {
	real rx, ry, ux, uy, lx, ly, bx, by;
	real r1, r2, r;
	LANPR_BoundingArea* ba; nListItemPointer* lip;
	if (PositiveX) {
		rx = This->R;
		ry = y + k * (rx - x);
		if (PositiveY) {
			uy = This->U;
			ux = x + (uy - y) / k;
			r1 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], rx);
			r2 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], ux);
			if (TNS_MIN2(r1, r2) > 1) return 0;
			if (r1 <= r2) {
				for (lip = This->RP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->U >= ry && ba->B < ry) { *NextX = rx; *NextY = ry; return ba; }
				}
			}else {
				for (lip = This->UP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->R >= ux && ba->L < ux) { *NextX = ux; *NextY = uy; return ba; }
				}
			}
		}else {
			by = This->B;
			bx = x + (by - y) / k;
			r1 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], rx);
			r2 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], bx);
			if (TNS_MIN2(r1, r2) > 1) return 0;
			if (r1 <= r2) {
				for (lip = This->RP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->U >= ry && ba->B < ry) { *NextX = rx; *NextY = ry; return ba; }
				}
			}else {
				for (lip = This->BP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->R >= bx && ba->L < bx) { *NextX = bx; *NextY = by; return ba; }
				}
			}
		}
	}else {
		lx = This->L;
		ly = y + k * (lx - x);
		if (PositiveY) {
			uy = This->U;
			ux = x + (uy - y) / k;
			r1 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], lx);
			r2 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], ux);
			if (TNS_MIN2(r1, r2) > 1) return 0;
			if (r1 <= r2) {
				for (lip = This->LP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->U >= ly && ba->B < ly) { *NextX = lx; *NextY = ly; return ba; }
				}
			}else {
				for (lip = This->UP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->R >= ux && ba->L < ux) { *NextX = ux; *NextY = uy; return ba; }
				}
			}
		}else {
			by = This->B;
			bx = x + (by - y) / k;
			r1 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], lx);
			r2 = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], bx);
			if (TNS_MIN2(r1, r2) > 1) return 0;
			if (r1 <= r2) {
				for (lip = This->LP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->U >= ly && ba->B < ly) { *NextX = lx; *NextY = ly; return ba; }
				}
			}else {
				for (lip = This->BP.pFirst; lip; lip = lip->pNext) {
					ba = lip->p;
					if (ba->R >= bx && ba->L < bx) { *NextX = bx; *NextY = by; return ba; }
				}
			}
		}
	}
	return 0;
}

LANPR_BoundingArea* lanpr_GetBoundingArea(LANPR_RenderBuffer* rb, real x, real y) {
	LANPR_BoundingArea* iba;
	real SpW = rb->WidthPerTile, SpH = rb->HeightPerTile;
	int c = (int)((x + 1.0) / SpW);
	int r = rb->TileCountY - (int)((y + 1.0) / SpH) - 1;
	if (r < 0) r = 0;
	if (c < 0) c = 0;
	if (r >= rb->TileCountY) r = rb->TileCountY-1;
	if (c >= rb->TileCountX) c = rb->TileCountX-1;

	iba = &rb->InitialBoundingAreas[r * 20 + c];
	while (iba->Child) {
		if (x > iba->CX) {
			if (y > iba->CY) iba = &iba->Child[0];
			else iba = &iba->Child[3];
		}else {
			if (y > iba->CY) iba = &iba->Child[1];
			else iba = &iba->Child[2];
		}
	}
	return iba;
}
LANPR_BoundingArea* lanpr_GetFirstPossibleBoundingArea(LANPR_RenderBuffer* rb, LANPR_RenderLine* rl) {
	LANPR_BoundingArea* iba;
	real p[2] = { rl->L->FrameBufferCoord[0], rl->L->FrameBufferCoord[1] };
	tnsVector2d LU = { -1,1 }, RU = { 1,1 }, LB = { -1,-1 }, RB = { 1, -1 };
	real r=1,sr=1;

	if (p[0] > -1 && p[0]<1 && p[1] > -1 && p[1]<1) {
		return lanpr_GetBoundingArea(rb, p[0], p[1]);
	}else {
		if (lanpr_LineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, LU, RU, &sr) && sr<r && sr>0)r = sr;
		if (lanpr_LineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, LB, RB, &sr) && sr<r && sr>0)r = sr;
		if (lanpr_LineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, LB, LU, &sr) && sr<r && sr>0)r = sr;
		if (lanpr_LineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, RB, RU, &sr) && sr<r && sr>0)r = sr;
		lanpr_LinearInterpolate2dv(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, r, p);

		return lanpr_GetBoundingArea(rb, p[0], p[1]);
	}

	real SpW = rb->WidthPerTile, SpH = rb->HeightPerTile;

	return iba;
}


/* ======================================= geometry ============================================ */

void lanpr_CutLineIntegrated(LANPR_RenderBuffer* rb, LANPR_RenderLine* rl, real Begin, real End) {
	LANPR_RenderLineSegment* rls = rl->Segments.pFirst, *irls;
	LANPR_RenderLineSegment *BeginSegment = 0, *EndSegment = 0;
	LANPR_RenderLineSegment *ns = 0, *ns2 = 0;
	LANPR_RenderLineSegment *BeforeBegin, *AfterEnd;
	LANPR_RenderLineSegment *Next;

	if (TNS_DOUBLE_CLOSE_ENOUGH(Begin, End)) return;

	if (Begin != Begin)
		Begin = 0;
	if (End != End)
		End = 0;

	if (Begin > End) {
		real t = Begin;
		Begin = End;
		End = t;
	}

	for (rls = rl->Segments.pFirst; rls; rls = rls->Item.pNext) {
		if (TNS_DOUBLE_CLOSE_ENOUGH(rls->at, Begin)) {
			BeginSegment = rls;
			ns = BeginSegment;
			break;
		}
		if (!rls->Item.pNext) {
			break;
		}
		irls = rls->Item.pNext;
		if (irls->at > Begin && Begin > rls->at) {
			BeginSegment = irls;
			ns = memStaticAquireThread(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
			break;
		}
	}
	for (rls = BeginSegment; rls; rls = rls->Item.pNext) {
		if (TNS_DOUBLE_CLOSE_ENOUGH(rls->at, End)) {
			EndSegment = rls;
			ns2 = EndSegment;
			break;
		}
		//irls = rls->Item.pNext;
		if (rls->at > End) {
			EndSegment = rls;
			ns2 = memStaticAquireThread(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
			break;
		}
	}

	if (!ns) ns = memStaticAquireThread(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
	if (!ns2) ns2 = memStaticAquireThread(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));

	if (BeginSegment) {
		if (BeginSegment != ns) {
			ns->OccludeLevel = BeginSegment->Item.pPrev ? (irls = BeginSegment->Item.pPrev)->OccludeLevel : 0;
			lstInsertItemBefore(&rl->Segments, (void*)ns, (void*)BeginSegment);
		}
	}
	else {
		ns->OccludeLevel = (irls = rl->Segments.pLast)->OccludeLevel;
		lstAppendItem(&rl->Segments, ns);
	}
	if (EndSegment) {
		if (EndSegment != ns2) {
			ns2->OccludeLevel = EndSegment->Item.pPrev ? (irls = EndSegment->Item.pPrev)->OccludeLevel : 0;
			lstInsertItemBefore(&rl->Segments, (void*)ns2, (void*)EndSegment);
		}
	}
	else {
		ns2->OccludeLevel = (irls = rl->Segments.pLast)->OccludeLevel;
		lstAppendItem(&rl->Segments, ns2);
	}

	ns->at = Begin;
	ns2->at = End;

	for (rls = ns; rls && rls != ns2; rls = rls->Item.pNext) {
		rls->OccludeLevel++;
	}
}


int lanpr_MakeNextOcclusionTaskInfo(LANPR_RenderBuffer* rb, LANPR_RenderTaskInfo* rti) {
	nListItemPointer* p;
	int i;
	int res = 0;

	//EnterCriticalSection(&rb->csManagement);

	if (rb->ContourManaged) {
		p = rb->ContourManaged;
		rti->Contour = (void*)p;
		rti->ContourPointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->ContourManaged = p;
		rti->ContourPointers.pLast = p ? p->pPrev : rb->Contours.pLast;
		res = 1;
	}
	else {
		lstEmptyDirect(&rti->ContourPointers);
		rti->Contour = 0;
	}

	if (rb->IntersectionManaged) {
		p = rb->IntersectionManaged;
		rti->Intersection = (void*)p;
		rti->IntersectionPointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->IntersectionManaged = p;
		rti->IntersectionPointers.pLast = p ? p->pPrev : rb->IntersectionLines.pLast;
		res = 1;
	}
	else {
		lstEmptyDirect(&rti->IntersectionPointers);
		rti->Intersection = 0;
	}

	if (rb->CreaseManaged) {
		p = rb->CreaseManaged;
		rti->Crease = (void*)p;
		rti->CreasePointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->CreaseManaged = p;
		rti->CreasePointers.pLast = p ? p->pPrev : rb->CreaseLines.pLast;
		res = 1;
	}
	else {
		lstEmptyDirect(&rti->CreasePointers);
		rti->Crease = 0;
	}

	if (rb->MaterialManaged) {
		p = rb->MaterialManaged;
		rti->Material = (void*)p;
		rti->MaterialPointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->MaterialManaged = p;
		rti->MaterialPointers.pLast = p ? p->pPrev : rb->MaterialLines.pLast;
		res = 1;
	}
	else {
		lstEmptyDirect(&rti->MaterialPointers);
		rti->Material = 0;
	}

	//LeaveCriticalSection(&rb->csManagement);

	return res;
}
void lanpr_CalculateSingleLineOcclusion(LANPR_RenderBuffer* rb, LANPR_RenderLine* rl, int ThreadID) {
	real x = rl->L->FrameBufferCoord[0], y = rl->L->FrameBufferCoord[1];
	LANPR_BoundingArea* ba = lanpr_GetFirstPossibleBoundingArea(rb, rl);
	LANPR_BoundingArea* nba = ba;
	LANPR_RenderTriangleThread* rt;
	nListItemPointer* lip;
	LANPR_RenderVert* rv;
	Camera* c = rb->Scene->ActiveCamera;
	real l, r;
	real k = (rl->R->FrameBufferCoord[1] - rl->L->FrameBufferCoord[1]) / (rl->R->FrameBufferCoord[0] - rl->L->FrameBufferCoord[0] + 1e-30);
	int PositiveX = (rl->R->FrameBufferCoord[0] - rl->L->FrameBufferCoord[0])>0 ? 1 : 0;
	int PositiveY = (rl->R->FrameBufferCoord[1] - rl->L->FrameBufferCoord[1])>0 ? 1 : 0;

	while (nba) {

		for (lip = nba->AssociatedTriangles.pFirst; lip; lip = lip->pNext) {
			rt = lip->p;
			//if (rt->Testing[ThreadID] == rl || rl->L->IntersectWith == rt || rl->R->IntersectWith == rt) continue;
			rt->Testing[ThreadID] = rl;
			if (lanpr_TriangleLineImageSpaceIntersectTestOnlyV2(rt, rl, c, rb->ViewProjection, &l, &r)) {
				lanpr_CutLineIntegrated(rb, rl, l, r);
			}
		}
		nba = lanpr_GetNextBoundingArea(nba, rl, x, y, k, PositiveX, PositiveY, &x, &y);
	}
}
void THREAD_CalculateLineOcclusion(LANPR_RenderTaskInfo* rti) {
	LANPR_RenderBuffer* rb = rti->RenderBuffer;
	int ThreadId = rti->ThreadID;
	LANPR_RenderLine* rl;
	nListItemPointer* lip;
	int count = 0;

	while (lanpr_MakeNextOcclusionTaskInfo(rb, rti)) {

		for (lip = (void*)rti->Contour; lip&& lip->pPrev != rti->ContourPointers.pLast; lip = lip->pNext) {
			lanpr_CalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);

			count++;
		}
		//tnsset_PlusRenderContourProcessedCount(rb, count);
		count = 0;

		for (lip = (void*)rti->Crease; lip && lip->pPrev != rti->CreasePointers.pLast; lip = lip->pNext) {
			lanpr_CalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);
			count++;
		}
		//tnsset_PlusRenderCreaseProcessedCount(rb, count);
		count = 0;

		for (lip = (void*)rti->Intersection; lip&& lip->pPrev != rti->IntersectionPointers.pLast; lip = lip->pNext) {
			lanpr_CalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);
			count++;
		}
		//tnsset_PlusRenderIntersectionProcessedCount(rb, count);
		count = 0;

		for (lip = (void*)rti->Material; lip&& lip->pPrev != rti->MaterialPointers.pLast; lip = lip->pNext) {
			lanpr_CalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);
			count++;
		}
		//tnsset_PlusRenderMaterialProcessedCount(rb, count);
		count = 0;

	}
	//thrd_exit(0);
}

int lanpr_GetNormal(tnsVector3d v1, tnsVector3d v2, tnsVector3d v3, tnsVector3d n, tnsVector3d Pos) {
	tnsVector3d vec1, vec2;

	tMatVectorMinus3d(vec1, v2, v1);
	tMatVectorMinus3d(vec2, v3, v1);
	tMatVectorCross3d(n, vec1, vec2);
	tMatNormalizeSelf3d(n);
	if (Pos && (tMatDot3d(n, Pos, 1) < 0)) {
		tMatVectorMultiSelf3d(n, -1.0f);
		return 1;
	}
	return 0;
}

int lanpr_BoundBoxCrosses(tnsVector4d xxyy1, tnsVector4d xxyy2) {
	real XMax1 = TNS_MAX2(xxyy1[0], xxyy1[1]);
	real XMin1 = TNS_MIN2(xxyy1[0], xxyy1[1]);
	real YMax1 = TNS_MAX2(xxyy1[2], xxyy1[3]);
	real YMin1 = TNS_MIN2(xxyy1[2], xxyy1[3]);
	real XMax2 = TNS_MAX2(xxyy2[0], xxyy2[1]);
	real XMin2 = TNS_MIN2(xxyy2[0], xxyy2[1]);
	real YMax2 = TNS_MAX2(xxyy2[2], xxyy2[3]);
	real YMin2 = TNS_MIN2(xxyy2[2], xxyy2[3]);

	if (XMax1<XMin2 || XMin1>XMax2) return 0;
	if (YMax1<YMin2 || YMin1>YMax2) return 0;

	return 1;
}
int lanpr_LineCrossesTile(LANPR_RenderBuffer* fb, tnsVector2d L, tnsVector2d R, LANPR_RenderTile* Tile) {
	real vx, vy;
	tnsVector4d Converted;
	real c1, c;

	if ((Converted[0] = (real)Tile->SubX - (real)fb->W / 2) > TNS_MAX2(L[0], R[0])) return 0;
	if ((Converted[1] = (real)Tile->SubXLim - (real)fb->W / 2) < TNS_MIN2(L[0], R[0])) return 0;
	if ((Converted[2] = (real)Tile->SubY - (real)fb->H / 2) > TNS_MAX2(L[1], R[1])) return 0;
	if ((Converted[3] = (real)Tile->SubYLim - (real)fb->H / 2) < TNS_MIN2(L[1], R[1])) return 0;

	//LineBoundBox[0] = L[0];
	//LineBoundBox[1] = R[0];
	//LineBoundBox[2] = L[1];
	//LineBoundBox[3] = R[1];

	//if (!lanpr_BoundBoxCrosses(Converted, LineBoundBox))return 0;

	//tMatVectorMinus2d(vec, L, R);
	vx = L[0] - R[0];
	vy = L[1] - R[1];

	c1 = vx * (Converted[2] - L[1]) - vy * (Converted[0] - L[0]);
	c = c1;

	c1 = vx * (Converted[2] - L[1]) - vy * (Converted[1] - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	c1 = vx * (Converted[3] - L[1]) - vy * (Converted[0] - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	c1 = vx * (Converted[3] - L[1]) - vy * (Converted[1] - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	//c1 = vec[0] * (Converted[2] - L[1]) - vec[1] * (Converted[0] - L[0]);
	//if (c1*c < 0)return 1;

	//test[0] = Converted[1]; test[1] = Converted[2];
	//tMatVectorMinus2d(test, test, L);
	//tMatVectorCrossOnly3d(CrossResult, vec, test);
	//if (tMatDot3d(CrossResult, CrossResultSave, 0) < 0) 
	//	return 1;
	//else tMatVectorCopy3d(CrossResult, CrossResultSave);

	//test[0] = Converted[0]; test[1] = Converted[3];
	//tMatVectorMinus2d(test, test, L);
	//tMatVectorCrossOnly3d(CrossResult, vec, test);
	//if (tMatDot3d(CrossResult, CrossResultSave, 0) < 0)
	//	return 1;
	//else tMatVectorCopy3d(CrossResult, CrossResultSave);
	//
	//test[0] = Converted[1]; test[1] = Converted[3];
	//tMatVectorMinus2d(test, test, L);
	//tMatVectorCrossOnly3d(CrossResult, vec, test);
	//if (tMatDot3d(CrossResult, CrossResultSave, 0) < 0) 
	//	return 1;

	return 0;
}
int lanpr_PointInsideTrianglef(tnsVector2d v, tnsVector2d v0, tnsVector2d v1, tnsVector2d v2) {
	double cl, c;

	cl = (v0[0] - v[0]) * (v1[1] - v[1]) - (v0[1] - v[1]) * (v1[0] - v[0]);
	c = cl;

	cl = (v1[0] - v[0]) * (v2[1] - v[1]) - (v1[1] - v[1]) * (v2[0] - v[0]);
	if (c*cl <= 0)return 0;
	else c = cl;

	cl = (v2[0] - v[0]) * (v0[1] - v[1]) - (v2[1] - v[1]) * (v0[0] - v[0]);
	if (c*cl <= 0)return 0;
	else c = cl;

	cl = (v0[0] - v[0]) * (v1[1] - v[1]) - (v0[1] - v[1]) * (v1[0] - v[0]);
	if (c*cl <= 0)return 0;

	return 1;
}
int lanpr_PointOnLinef(tnsVector2d v, tnsVector2d v0, tnsVector2d v1) {
	real c1, c2;

	c1 = tMatGetLinearRatio(v0[0], v1[0], v[0]);
	c2 = tMatGetLinearRatio(v0[1], v1[1], v[1]);

	if (TNS_DOUBLE_CLOSE_ENOUGH(c1, c2) && c1 >= 0 && c1 <= 1) return 1;

	return 0;
}
int lanpr_PointTriangleRelation(tnsVector2d v, tnsVector2d v0, tnsVector2d v1, tnsVector2d v2) {
	double cl, c;
	real r;
	if (lanpr_PointOnLinef(v, v0, v1) || lanpr_PointOnLinef(v, v1, v2) || lanpr_PointOnLinef(v, v2, v0)) return 1;

	cl = (v0[0] - v[0]) * (v1[1] - v[1]) - (v0[1] - v[1]) * (v1[0] - v[0]);
	c = cl;

	cl = (v1[0] - v[0]) * (v2[1] - v[1]) - (v1[1] - v[1]) * (v2[0] - v[0]);
	if ((r = c*cl) < 0) return 0;
	elif(r == 0) return 1;
	else c = cl;

	cl = (v2[0] - v[0]) * (v0[1] - v[1]) - (v2[1] - v[1]) * (v0[0] - v[0]);
	if ((r = c*cl) < 0) return 0;
	elif(r == 0) return 1;
	else c = cl;

	cl = (v0[0] - v[0]) * (v1[1] - v[1]) - (v0[1] - v[1]) * (v1[0] - v[0]);
	if ((r = c*cl) < 0) return 0;
	elif(r == 0) return 1;

	return 2;
}
int lanpr_PointInsideTriangle3d(tnsVector3d v, tnsVector3d v0, tnsVector3d v1, tnsVector3d v2) {
	real cl, c;
	tnsVector3d L, R;
	tnsVector3d N1, N2;

	tMatVectorMinus3d(L, v1, v0);
	tMatVectorMinus3d(R, v, v1);
	tMatVectorCross3d(N1, L, R);

	tMatVectorMinus3d(L, v2, v1);
	tMatVectorMinus3d(R, v, v2);
	tMatVectorCross3d(N2, L, R);

	if (tMatDot3d(N1, N2, 0) < 0) return 0;

	tMatVectorMinus3d(L, v0, v2);
	tMatVectorMinus3d(R, v, v0);
	tMatVectorCross3d(N1, L, R);

	if (tMatDot3d(N1, N2, 0) < 0) return 0;

	tMatVectorMinus3d(L, v1, v0);
	tMatVectorMinus3d(R, v, v1);
	tMatVectorCross3d(N2, L, R);

	if (tMatDot3d(N1, N2, 0) < 0) return 0;

	return 1;
}
int lanpr_PointInsideTriangle3de(tnsVector3d v, tnsVector3d v0, tnsVector3d v1, tnsVector3d v2) {
	tnsVector3d L, R;
	tnsVector3d N1, N2;
	real d;

	tMatVectorMinus3d(L, v1, v0);
	tMatVectorMinus3d(R, v, v1);
	//tMatNormalizeSelf3d(L);
	//tMatNormalizeSelf3d(R);
	tMatVectorCross3d(N1, L, R);

	tMatVectorMinus3d(L, v2, v1);
	tMatVectorMinus3d(R, v, v2);
	//tMatNormalizeSelf3d(L);
	//tMatNormalizeSelf3d(R);
	tMatVectorCross3d(N2, L, R);

	if ((d = tMatDot3d(N1, N2, 0)) < 0)return 0;
	//if (d<DBL_EPSILON) return -1;

	tMatVectorMinus3d(L, v0, v2);
	tMatVectorMinus3d(R, v, v0);
	//tMatNormalizeSelf3d(L);
	//tMatNormalizeSelf3d(R);
	tMatVectorCross3d(N1, L, R);

	if ((d = tMatDot3d(N1, N2, 0)) < 0) return 0;
	//if (d<DBL_EPSILON) return -1;

	tMatVectorMinus3d(L, v1, v0);
	tMatVectorMinus3d(R, v, v1);
	//tMatNormalizeSelf3d(L);
	//tMatNormalizeSelf3d(R);
	tMatVectorCross3d(N2, L, R);

	if ((d = tMatDot3d(N1, N2, 0)) < 0) return 0;
	//if (d<DBL_EPSILON) return -1;

	return 1;
}

LANPR_RenderElementLinkNode* lanpr_NewCullTriangleSpace64(LANPR_RenderBuffer* rb) {
	LANPR_RenderElementLinkNode* reln;

	LANPR_RenderTriangle* RenderTriangles = calloc(64, rb->TriangleSize);//CreateNewBuffer(LANPR_RenderTriangle, 64);

	reln = lstAppendPointerStaticSized(&rb->TriangleBufferPointers, &rb->RenderDataPool, RenderTriangles,
		sizeof(LANPR_RenderElementLinkNode));
	reln->ElementCount = 64;
	reln->Additional = 1;

	return reln;
}
LANPR_RenderElementLinkNode* lanpr_NewCullPointSpace64(LANPR_RenderBuffer* rb) {
	LANPR_RenderElementLinkNode* reln;

	LANPR_RenderVert* RenderVertices = CreateNewBuffer(LANPR_RenderVert, 64);

	reln = lstAppendPointerStaticSized(&rb->VertexBufferPointers, &rb->RenderDataPool, RenderVertices,
		sizeof(LANPR_RenderElementLinkNode));
	reln->ElementCount = 64;
	reln->Additional = 1;

	return reln;
}
void lanpr_CalculateRenderTriangleNormal(LANPR_RenderTriangle* rt);
void lanpr_PostTriangle(LANPR_RenderTriangle* rt, LANPR_RenderTriangle* orig) {
	if (rt->V[0])tMatVectorAccum3d(rt->GC, rt->V[0]->FrameBufferCoord);
	if (rt->V[1])tMatVectorAccum3d(rt->GC, rt->V[1]->FrameBufferCoord);
	if (rt->V[2])tMatVectorAccum3d(rt->GC, rt->V[2]->FrameBufferCoord);
	tMatVectorMultiSelf3d(rt->GC, 1.0f / 3.0f);

	tMatVectorCopy3d(orig->GN, rt->GN);
}
void lanpr_CullTriangles(LANPR_RenderBuffer* rb) {
	LANPR_RenderLine* rl;
	LANPR_RenderTriangle* rt, *rt1;
	LANPR_RenderVert* rv;
	LANPR_RenderElementLinkNode* reln, *veln, *teln;
	LANPR_RenderLineSegment* rls;
	tnsVector4d p1, p2;
	real* MVInverse = rb->VPInverse;
	int i;
	real a;
	int VCount = 0, TCount = 0;
	Object* o;

	veln = lanpr_NewCullPointSpace64(rb);
	teln = lanpr_NewCullTriangleSpace64(rb);
	rv = &((LANPR_RenderVert*)veln->Pointer)[VCount];
	rt1 = (void*)(((BYTE*)teln->Pointer) + rb->TriangleSize*TCount);

	for (reln = rb->TriangleBufferPointers.pFirst; reln; reln = reln->Item.pNext) {
		i = 0;
		if (reln->Additional) continue;
		o = reln->ObjectRef;
		for (i; i < reln->ElementCount; i++) {
			int In1 = 0, In2 = 0, In3 = 0;
			rt = (void*)(((BYTE*)reln->Pointer) + rb->TriangleSize*i);
			if (rt->V[0]->FrameBufferCoord[3] < 0) In1 = 1;
			if (rt->V[1]->FrameBufferCoord[3] < 0) In2 = 1;
			if (rt->V[2]->FrameBufferCoord[3] < 0) In3 = 1;

			rt->RL[0]->ObjectRef = o;
			rt->RL[1]->ObjectRef = o;
			rt->RL[2]->ObjectRef = o;

			if (VCount > 60) {
				veln->ElementCount = VCount;
				veln = lanpr_NewCullPointSpace64(rb);
				VCount = 0;
			}

			if (TCount > 60) {
				teln->ElementCount = TCount;
				teln = lanpr_NewCullTriangleSpace64(rb);
				TCount = 0;
			}

			if ((!rt->RL[0]->Item.pNext && !rt->RL[0]->Item.pPrev) ||
				(!rt->RL[1]->Item.pNext && !rt->RL[1]->Item.pPrev) ||
				(!rt->RL[2]->Item.pNext && !rt->RL[2]->Item.pPrev)) {
				printf("'");
			}

			rv = &((LANPR_RenderVert*)veln->Pointer)[VCount];
			rt1 = &((LANPR_RenderTriangle*)teln->Pointer)[TCount];


			switch (In1 + In2 + In3) {
			case 0:
				continue;
			case 3:
				rt->CullStatus = TNS_CULL_DISCARD;
				continue;
			case 2:
				rt->CullStatus = TNS_CULL_USED;
				if (!In1) {
					a = rt->V[0]->FrameBufferCoord[2] / (rt->V[0]->FrameBufferCoord[2] - rt->V[2]->FrameBufferCoord[2]);
					rv[0].FrameBufferCoord[0] = (1 - a)*rt->V[0]->FrameBufferCoord[0] + a* rt->V[2]->FrameBufferCoord[0];
					rv[0].FrameBufferCoord[1] = (1 - a)*rt->V[0]->FrameBufferCoord[1] + a* rt->V[2]->FrameBufferCoord[1];
					rv[0].FrameBufferCoord[2] = 0;
					rv[0].FrameBufferCoord[3] = (1 - a)*rt->V[0]->FrameBufferCoord[3] + a* rt->V[2]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[0].GLocation, MVInverse, rv[0].FrameBufferCoord);

					a = rt->V[0]->FrameBufferCoord[2] / (rt->V[0]->FrameBufferCoord[2] - rt->V[1]->FrameBufferCoord[2]);
					rv[1].FrameBufferCoord[0] = (1 - a)*rt->V[0]->FrameBufferCoord[0] + a* rt->V[1]->FrameBufferCoord[0];
					rv[1].FrameBufferCoord[1] = (1 - a)*rt->V[0]->FrameBufferCoord[1] + a* rt->V[1]->FrameBufferCoord[1];
					rv[1].FrameBufferCoord[2] = 0;
					rv[1].FrameBufferCoord[3] = (1 - a)*rt->V[0]->FrameBufferCoord[3] + a* rt->V[1]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[1].GLocation, MVInverse, rv[1].FrameBufferCoord);

					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1->RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = rt->V[0];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1->RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[0];
					rl->R = &rv[0];
					rl->TL = rt->RL[2]->TL == rt ? rt1 : rt->RL[2]->TL;
					rl->TR = rt->RL[2]->TR == rt ? rt1 : rt->RL[2]->TR;
					rt1->RL[2] = rl;

					rt1->V[0] = rt->V[0];
					rt1->V[1] = &rv[1];
					rt1->V[2] = &rv[0];

					lanpr_PostTriangle(rt1, rt);

					VCount += 2;
					TCount += 1;
					continue;
				}elif(!In3) {
					a = rt->V[2]->FrameBufferCoord[2] / (rt->V[2]->FrameBufferCoord[2] - rt->V[0]->FrameBufferCoord[2]);
					rv[0].FrameBufferCoord[0] = (1 - a)*rt->V[2]->FrameBufferCoord[0] + a* rt->V[0]->FrameBufferCoord[0];
					rv[0].FrameBufferCoord[1] = (1 - a)*rt->V[2]->FrameBufferCoord[1] + a* rt->V[0]->FrameBufferCoord[1];
					rv[0].FrameBufferCoord[2] = 0;
					rv[0].FrameBufferCoord[3] = (1 - a)*rt->V[2]->FrameBufferCoord[3] + a* rt->V[0]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[0].GLocation, MVInverse, rv[0].FrameBufferCoord);

					a = rt->V[2]->FrameBufferCoord[2] / (rt->V[2]->FrameBufferCoord[2] - rt->V[1]->FrameBufferCoord[2]);
					rv[1].FrameBufferCoord[0] = (1 - a)*rt->V[2]->FrameBufferCoord[0] + a* rt->V[1]->FrameBufferCoord[0];
					rv[1].FrameBufferCoord[1] = (1 - a)*rt->V[2]->FrameBufferCoord[1] + a* rt->V[1]->FrameBufferCoord[1];
					rv[1].FrameBufferCoord[2] = 0;
					rv[1].FrameBufferCoord[3] = (1 - a)*rt->V[2]->FrameBufferCoord[3] + a* rt->V[1]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[1].GLocation, MVInverse, rv[1].FrameBufferCoord);

					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = &rv[1];
					rl->TL = rt1;
					rt1->RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = rt->V[2];
					rl->TL = rt->RL[1]->TL == rt ? rt1 : rt->RL[1]->TL;
					rl->TR = rt->RL[1]->TR == rt ? rt1 : rt->RL[1]->TR;
					rt1->RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[2];
					rl->R = &rv[0];
					rl->TL = rt->RL[2]->TL == rt ? rt1 : rt->RL[2]->TL;
					rl->TR = rt->RL[2]->TR == rt ? rt1 : rt->RL[2]->TR;
					rt1->RL[2] = rl;

					rt1->V[0] = &rv[1];
					rt1->V[1] = rt->V[2];
					rt1->V[2] = &rv[0];

					lanpr_PostTriangle(rt1, rt);

					VCount += 2;
					TCount += 1;
					continue;
				}elif(!In2) {
					a = rt->V[1]->FrameBufferCoord[2] / (rt->V[1]->FrameBufferCoord[2] - rt->V[0]->FrameBufferCoord[2]);
					rv[0].FrameBufferCoord[0] = (1 - a)*rt->V[1]->FrameBufferCoord[0] + a* rt->V[0]->FrameBufferCoord[0];
					rv[0].FrameBufferCoord[1] = (1 - a)*rt->V[1]->FrameBufferCoord[1] + a* rt->V[0]->FrameBufferCoord[1];
					rv[0].FrameBufferCoord[2] = 0;
					rv[0].FrameBufferCoord[3] = (1 - a)*rt->V[1]->FrameBufferCoord[3] + a* rt->V[0]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[0].GLocation, MVInverse, rv[0].FrameBufferCoord);

					a = rt->V[1]->FrameBufferCoord[2] / (rt->V[1]->FrameBufferCoord[2] - rt->V[2]->FrameBufferCoord[2]);
					rv[1].FrameBufferCoord[0] = (1 - a)*rt->V[1]->FrameBufferCoord[0] + a* rt->V[2]->FrameBufferCoord[0];
					rv[1].FrameBufferCoord[1] = (1 - a)*rt->V[1]->FrameBufferCoord[1] + a* rt->V[2]->FrameBufferCoord[1];
					rv[1].FrameBufferCoord[2] = 0;
					rv[1].FrameBufferCoord[3] = (1 - a)*rt->V[1]->FrameBufferCoord[3] + a* rt->V[2]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[1].GLocation, MVInverse, rv[1].FrameBufferCoord);

					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1->RL[2] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[1];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1->RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[1];
					rl->R = &rv[1];
					rl->TL = rt->RL[1]->TL == rt ? rt1 : rt->RL[1]->TL;
					rl->TR = rt->RL[1]->TR == rt ? rt1 : rt->RL[1]->TR;
					rt1->RL[1] = rl;

					rt1->V[0] = rt->V[1];
					rt1->V[1] = &rv[1];
					rt1->V[2] = &rv[0];

					lanpr_PostTriangle(rt1, rt);

					VCount += 2;
					TCount += 1;
					continue;
				}
				break;
			case 1:
				rt->CullStatus = TNS_CULL_USED;
				if (In1) {
					a = rt->V[0]->FrameBufferCoord[2] / (rt->V[0]->FrameBufferCoord[2] - rt->V[2]->FrameBufferCoord[2]);
					rv[0].FrameBufferCoord[0] = (1 - a)*rt->V[0]->FrameBufferCoord[0] + a* rt->V[2]->FrameBufferCoord[0];
					rv[0].FrameBufferCoord[1] = (1 - a)*rt->V[0]->FrameBufferCoord[1] + a* rt->V[2]->FrameBufferCoord[1];
					rv[0].FrameBufferCoord[2] = 0;
					rv[0].FrameBufferCoord[3] = (1 - a)*rt->V[0]->FrameBufferCoord[3] + a* rt->V[2]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[0].GLocation, MVInverse, rv[0].FrameBufferCoord);

					a = rt->V[0]->FrameBufferCoord[2] / (rt->V[0]->FrameBufferCoord[2] - rt->V[1]->FrameBufferCoord[2]);
					rv[1].FrameBufferCoord[0] = (1 - a)*rt->V[0]->FrameBufferCoord[0] + a* rt->V[1]->FrameBufferCoord[0];
					rv[1].FrameBufferCoord[1] = (1 - a)*rt->V[0]->FrameBufferCoord[1] + a* rt->V[1]->FrameBufferCoord[1];
					rv[1].FrameBufferCoord[2] = 0;
					rv[1].FrameBufferCoord[3] = (1 - a)*rt->V[0]->FrameBufferCoord[3] + a* rt->V[1]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[1].GLocation, MVInverse, rv[1].FrameBufferCoord);

					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1[0].RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[1];
					rl->TL = &rt1[0];
					rl->TR = &rt1[1];
					rt1[0].RL[2] = rl;
					rt1[1].RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[1];
					rl->R = &rv[1];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1[0].RL[0] = rl;

					rt1[0].V[0] = rt->V[1];
					rt1[0].V[1] = &rv[1];
					rt1[0].V[2] = &rv[0];

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[2];
					rl->R = &rv[0];
					rl->TL = rt->RL[2]->TL == rt ? rt1 : rt->RL[2]->TL;
					rl->TR = rt->RL[2]->TR == rt ? rt1 : rt->RL[2]->TR;
					rt1[1].RL[2] = rl;
					rt1[1].RL[1] = rt->RL[1];

					rt1[1].V[0] = &rv[0];
					rt1[1].V[1] = rt->V[1];
					rt1[1].V[2] = rt->V[2];

					lanpr_PostTriangle(&rt1[0], rt);
					lanpr_PostTriangle(&rt1[1], rt);

					VCount += 2;
					TCount += 2;
					continue;
				}elif(In2) {
					a = rt->V[1]->FrameBufferCoord[2] / (rt->V[1]->FrameBufferCoord[2] - rt->V[0]->FrameBufferCoord[2]);
					rv[0].FrameBufferCoord[0] = (1 - a)*rt->V[1]->FrameBufferCoord[0] + a* rt->V[0]->FrameBufferCoord[0];
					rv[0].FrameBufferCoord[1] = (1 - a)*rt->V[1]->FrameBufferCoord[1] + a* rt->V[0]->FrameBufferCoord[1];
					rv[0].FrameBufferCoord[2] = 0;
					rv[0].FrameBufferCoord[3] = (1 - a)*rt->V[1]->FrameBufferCoord[3] + a* rt->V[0]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[0].GLocation, MVInverse, rv[0].FrameBufferCoord);

					a = rt->V[1]->FrameBufferCoord[2] / (rt->V[1]->FrameBufferCoord[2] - rt->V[2]->FrameBufferCoord[2]);
					rv[1].FrameBufferCoord[0] = (1 - a)*rt->V[1]->FrameBufferCoord[0] + a* rt->V[2]->FrameBufferCoord[0];
					rv[1].FrameBufferCoord[1] = (1 - a)*rt->V[1]->FrameBufferCoord[1] + a* rt->V[2]->FrameBufferCoord[1];
					rv[1].FrameBufferCoord[2] = 0;
					rv[1].FrameBufferCoord[3] = (1 - a)*rt->V[1]->FrameBufferCoord[3] + a* rt->V[2]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[1].GLocation, MVInverse, rv[1].FrameBufferCoord);

					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[1]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1[0].RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[2];
					rl->TL = &rt1[0];
					rl->TR = &rt1[1];
					rt1[0].RL[2] = rl;
					rt1[1].RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[2];
					rl->R = &rv[1];
					rl->TL = rt->RL[1]->TL == rt ? rt1 : rt->RL[1]->TL;
					rl->TR = rt->RL[1]->TR == rt ? rt1 : rt->RL[1]->TR;
					rt1[0].RL[0] = rl;

					rt1[0].V[0] = rt->V[2];
					rt1[0].V[1] = &rv[1];
					rt1[0].V[2] = &rv[0];

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[0];
					rl->R = &rv[0];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1[1].RL[2] = rl;
					rt1[1].RL[1] = rt->RL[2];

					rt1[1].V[0] = &rv[0];
					rt1[1].V[1] = rt->V[2];
					rt1[1].V[2] = rt->V[0];

					lanpr_PostTriangle(&rt1[0], rt);
					lanpr_PostTriangle(&rt1[1], rt);

					VCount += 2;
					TCount += 2;
					continue;
				}elif(In3) {
					a = rt->V[2]->FrameBufferCoord[2] / (rt->V[2]->FrameBufferCoord[2] - rt->V[0]->FrameBufferCoord[2]);
					rv[0].FrameBufferCoord[0] = (1 - a)*rt->V[2]->FrameBufferCoord[0] + a* rt->V[0]->FrameBufferCoord[0];
					rv[0].FrameBufferCoord[1] = (1 - a)*rt->V[2]->FrameBufferCoord[1] + a* rt->V[0]->FrameBufferCoord[1];
					rv[0].FrameBufferCoord[2] = 0;
					rv[0].FrameBufferCoord[3] = (1 - a)*rt->V[2]->FrameBufferCoord[3] + a* rt->V[0]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[0].GLocation, MVInverse, rv[0].FrameBufferCoord);

					a = rt->V[2]->FrameBufferCoord[2] / (rt->V[2]->FrameBufferCoord[2] - rt->V[1]->FrameBufferCoord[2]);
					rv[1].FrameBufferCoord[0] = (1 - a)*rt->V[2]->FrameBufferCoord[0] + a* rt->V[1]->FrameBufferCoord[0];
					rv[1].FrameBufferCoord[1] = (1 - a)*rt->V[2]->FrameBufferCoord[1] + a* rt->V[1]->FrameBufferCoord[1];
					rv[1].FrameBufferCoord[2] = 0;
					rv[1].FrameBufferCoord[3] = (1 - a)*rt->V[2]->FrameBufferCoord[3] + a* rt->V[1]->FrameBufferCoord[3];
					tMatApplyTransform44dTrue(rv[1].GLocation, MVInverse, rv[1].FrameBufferCoord);

					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, (void*)rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1[0].RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[0];
					rl->TL = &rt1[0];
					rl->TR = &rt1[1];
					rt1[0].RL[2] = rl;
					rt1[1].RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[2];
					rl->R = &rv[1];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1[0].RL[0] = rl;

					rt1[0].V[0] = rt->V[0];
					rt1[0].V[1] = &rv[1];
					rt1[0].V[2] = &rv[0];

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = rt->V[1];
					rl->R = &rv[0];
					rl->TL = rt->RL[1]->TL == rt ? rt1 : rt->RL[1]->TL;
					rl->TR = rt->RL[1]->TR == rt ? rt1 : rt->RL[1]->TR;
					rt1[1].RL[2] = rl;
					rt1[1].RL[1] = rt->RL[1];

					rt1[1].V[0] = &rv[0];
					rt1[1].V[1] = rt->V[1];
					rt1[1].V[2] = rt->V[2];

					lanpr_PostTriangle(&rt1[0], rt);
					lanpr_PostTriangle(&rt1[1], rt);

					VCount += 2;
					TCount += 2;
					continue;
				}
				break;
			}
		}
		teln->ElementCount = TCount;
		veln->ElementCount = VCount;
	}
}
void lanpr_PerspectiveDivision(LANPR_RenderBuffer* rb) {
	LANPR_RenderVert* rv;
	LANPR_RenderElementLinkNode* reln;
	Camera* cam = rb->Scene->ActiveCamera;
	int i;

	if (cam->CameraType == TNS_CAMERA_ORTHO) return;

	for (reln = rb->VertexBufferPointers.pFirst; reln; reln = reln->Item.pNext) {
		rv = reln->Pointer;
		for (i = 0; i < reln->ElementCount; i++) {
			//if (rv->FrameBufferCoord[2] < -DBL_EPSILON) continue;
			tMatVectorMultiSelf3d(rv[i].FrameBufferCoord, 1 / rv[i].FrameBufferCoord[3]);
			rv[i].FrameBufferCoord[2] = cam->ZMin*cam->ZMax / (cam->ZMax - fabs(rv[i].FrameBufferCoord[2]) * (cam->ZMax - cam->ZMin));
		}
	}
}

void lanpr_TransformRenderVert(BMVert* V, LANPR_RenderVert* RVBuf, real* MVMat, real* MVPMat, Camera* Camera) {//real HeightMultiply, real ZMin, real ZMax) {
	LANPR_RenderVert* rv = &RVBuf[V->I];
	rv->V = V;
	V->RV = rv;
	tMatApplyTransform43d(rv->GLocation, MVMat, V->P);
	tMatApplyTransform44d(rv->FrameBufferCoord, MVPMat, V->P);

	//if(rv->FrameBufferCoord[2]>0)tMatVectorMultiSelf3d(rv->FrameBufferCoord, (1 / rv->FrameBufferCoord[3]));
	//else tMatVectorMultiSelf3d(rv->FrameBufferCoord, -rv->FrameBufferCoord[3]);
	//   rv->FrameBufferCoord[2] = Camera->ZMin* Camera->ZMax / (Camera->ZMax - fabs(rv->FrameBufferCoord[2]) * (Camera->ZMax - Camera->ZMin));
}
void lanpr_CalculateRenderTriangleNormal(LANPR_RenderTriangle* rt) {
	tnsVector3d L, R;
	tMatVectorMinus3d(L, rt->V[1]->GLocation, rt->V[0]->GLocation);
	tMatVectorMinus3d(R, rt->V[2]->GLocation, rt->V[0]->GLocation);
	tMatVectorCross3d(rt->GN, L, R);
	tMatNormalizeSelf3d(rt->GN);
}
void lanpr_MakeRenderGeometryBuffersRecursive(Object* o, real* MVMat, real* MVPMat, LANPR_RenderBuffer* rb, real HeightMultiply) {
	Object* oc;
	Mesh* mo;
	BMVert* v;
	BMFace* f;
	BMEdge* e;
	LANPR_RenderTriangle* rt;
	tnsMatrix44d NewMVP;
	tnsMatrix44d NewMV;
	LANPR_RenderBuffer* fb = rb->FrameBuffer;
	LANPR_RenderElementLinkNode* reln;
	Camera* c = rb->Scene->ActiveCamera;
	Material* m;

	//if (o->RenderTriangles) FreeMem(o->RenderTriangles);
	//if (o->RenderVertices) FreeMem(o->RenderVertices);

	tMatMultiply44d(NewMVP, MVPMat, o->SelfTransform);
	tMatMultiply44d(NewMV, MVMat, o->SelfTransform);

	if (o->Type == TNS_OBJECT_MESH) {
		mo = o;
		o->RenderVertices = CreateNewBuffer(LANPR_RenderVert, mo->numV);
		o->RenderTriangles = calloc(mo->TriangleCount, rb->TriangleSize);//CreateNewBuffer(LANPR_RenderTriangle, mo->TriangleCount);
																		 //o->RenderLines = CreateNewBuffer(LANPR_RenderLine, mo->TriangulatedEdgeCount);

		reln = lstAppendPointerStaticSized(&rb->VertexBufferPointers, &rb->RenderDataPool, o->RenderVertices,
			sizeof(LANPR_RenderElementLinkNode));
		reln->ElementCount = mo->numV;
		reln->ObjectRef = mo;

		reln = lstAppendPointerStaticSized(&rb->TriangleBufferPointers, &rb->RenderDataPool, o->RenderTriangles,
			sizeof(LANPR_RenderElementLinkNode));
		reln->ElementCount = mo->TriangleCount;
		reln->ObjectRef = mo;

		for (v = mo->V.pFirst; v; v = v->Item.pNext) {
			lanpr_TransformRenderVert(v, o->RenderVertices, NewMV, NewMVP, c);//,HeightMultiply,c->ZMin,c->ZMax);
			tObjRecalculateVertNormal(v);
		}

		for (e = mo->E.pFirst; e; e = e->Item.pNext) {
			LANPR_RenderLine* rl = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLine));
			rl->L = e->VL->RV;
			rl->R = e->VR->RV;
			e->RenderLine = rl;
			LANPR_RenderLineSegment* rls = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineSegment));
			lstAppendItem(&rl->Segments, rls);
			lstAppendItem(&rb->AllRenderLines, rl);
		}


		rt = o->RenderTriangles;
		for (f = mo->F.pFirst; f; f = f->Item.pNext) {
			tObjRecalculateFaceAverageNormal(f);
			tObjSimpleTriangulateRender(f, rt, rb->TriangleSize, o->RenderVertices, &rb->AllRenderLines, &rt, &rb->RenderDataPool);
			// already done in the func above. lanpr_CalculateRenderTriangleNormal(rt);
			tMatApplyNormalTransform43d(f->GNormal, MVMat, f->FaceNormal);
			tMatNormalizeSelf3d(f->GNormal);
			m = tnsGetIndexedMaterial(rb->Scene, f->MaterialID);
			//if(m) m->PreviewVCount += (f->TriangleCount*3);
		}
	}

	for (oc = o->ChildObjects.pFirst; oc; oc = oc->Item.pNext) {
		lanpr_MakeRenderGeometryBuffersRecursive(oc, NewMV, NewMVP, rb, HeightMultiply);
	}
}
void tnsMakeRenderGeometryBuffers(Scene* s, Camera* c, LANPR_RenderBuffer* rb, int HeightMultiply) {
	Object* o;
	tnsMatrix44d proj, view, result, inv;
	LANPR_RenderBuffer* fb = rb->FrameBuffer;

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
		tObjApplyGlobalTransformMatrixRecursive(o);
		lanpr_MakeRenderGeometryBuffersRecursive(o, view, proj, rb, (real)HeightMultiply);
	}
}

void tnsZeroGeomtryBuffers(tnsScene* s) {
	Object* o;
	Object* oc;
	if (!s) return;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		lanpr_ZeroGeomtryBuffersRecursive(o);
	}
}



void lanpr_generate_geom_buffer(struct Object *ob){
    
}

static int lanpr_compute_feature_lines_exec(struct bContext *C, struct wmOperator *op){
	Scene *scene = CTX_data_scene(C);
	SceneLANPR* lanpr = &scene->lanpr;
    LANPR_RenderBuffer* rb;
	
	/* need threading, later.... */

    rb = lanpr_CreateRenderBuffer(lanpr);

	lanpr_MakeInitialBoundingAreas(rb);


	return OPERATOR_FINISHED;
}

static void lanpr_compute_feature_lines_cancel(struct bContext *C, struct wmOperator *op){

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
	ot->exec = lanpr_compute_feature_lines_exec;
}



/* ====================================== render control ======================================= */

void lanpr_DestroyRenderData(LANPR_RenderBuffer* rb) {
	LANPR_RenderElementLinkNode* reln;

	rb->ContourCount = 0;
	rb->ContourManaged = 0;
	rb->IntersectionCount = 0;
	rb->IntersectionManaged = 0;
	rb->MaterialLineCount = 0;
	rb->MaterialManaged = 0;
	rb->CreaseCount = 0;
	rb->CreaseManaged = 0;
	rb->CalculationStatus = TNS_CALCULATION_IDLE;

	lstEmptyDirect(&rb->Contours);
	lstEmptyDirect(&rb->IntersectionLines);
	lstEmptyDirect(&rb->CreaseLines);
	lstEmptyDirect(&rb->MaterialLines);
	lstEmptyDirect(&rb->AllRenderLines);

	//tnsZeroGeomtryBuffers(rb->Scene);

	while (reln = lstPopItem(&rb->VertexBufferPointers)) {
		FreeMem(reln->Pointer);
	}

	while (reln = lstPopItem(&rb->TriangleBufferPointers)) {
		FreeMem(reln->Pointer);
	}

	memStaticDestroy(&rb->RenderDataPool);
}

LANPR_RenderBuffer* lanpr_CreateRenderBuffer(SceneLANPR* lanpr) {
	if (lanpr->render_buffer) {
		lanpr_DestroyRenderData(lanpr->render_buffer);
		mem_free(lanpr->render_buffer);
	}

	LANPR_RenderBuffer* rb = MEM_callocN(sizeof(LANPR_RenderBuffer), "creating LANPR render buffer");

	lanpr->render_buffer = rb;

	return rb;
}


/* internal */

LANPR_LineStyle* lanpr_new_line_layer(SceneLANPR* lanpr){
    LANPR_LineStyle* ls = MEM_callocN(sizeof(LANPR_LineStyle),"Line Style");
	BLI_addtail(&lanpr->line_style_layers,ls);
	lanpr->active_layer = ls;
	return ls;
}

static int lanpr_add_line_layer_exec(struct bContext *C, struct wmOperator *op){
	Scene *scene = CTX_data_scene(C);
	SceneLANPR* lanpr = &scene->lanpr;

    lanpr_new_line_layer(lanpr);

	return OPERATOR_FINISHED;
}


void SCENE_OT_lanpr_add_line_layer(struct wmOperatorType* ot){
	
	ot->name = "Add line layer";
	ot->description = "Add a new line layer";
	ot->idname = "SCENE_OT_lanpr_add_line_layer";

	ot->exec = lanpr_add_line_layer_exec;

}

static int lanpr_delete_line_layer_exec(struct bContext *C, struct wmOperator *op){

	return OPERATOR_FINISHED;
}


void SCENE_OT_lanpr_delete_line_layer(struct wmOperatorType* ot){
	
	ot->name = "Delete line layer";
	ot->description = "Delete selected line layer";
	ot->idname = "SCENE_OT_lanpr_delete_line_layer";

	ot->exec = lanpr_delete_line_layer_exec;
	
}