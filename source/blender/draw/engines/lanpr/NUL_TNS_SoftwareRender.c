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

extern tnsMain* T;
extern NUL MAIN;

static int tnsglobal_TriangleIntersectionCount;


void tnsset_RenderOverallProgress(tnsRenderBuffer* rb, real value);
void tnsset_PlusRenderContourCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderCreaseCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderMaterialCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderIntersectionCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderContourProcessedCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderCreaseProcessedCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderMaterialProcessedCount(tnsRenderBuffer* rb, int value);
void tnsset_PlusRenderIntersectionProcessedCount(tnsRenderBuffer* rb, int value);



void tRdrMakeFakeFrameBuffer(tnsRenderBuffer* rb, int W, int H, int Samples) {
	tnsFrameBuffer* fb = memAquireOnly(sizeof(tnsFrameBuffer));
	fb->SubPixelSample = Samples;
	fb->W = W;
	fb->H = H;
	rb->FrameBuffer = fb;
	fb->OutputAALevel = TNS_OUTPUT_AA_16;
}
void tRdrMakeFrameBuffer(tnsRenderBuffer* rb,int W,int H,int Samples) {
	tnsFrameBuffer* fb = memAquireOnly(sizeof(tnsFrameBuffer));
	//fb->Pixels = CreateNewBuffer(tnsRenderSubPixel, W*H*Samples*Samples);

	//if (!fb->Pixels) MessageBox(0, "�ڴ�ը��", "�ڴ�ը��", 0);

	fb->SubPixelSample = Samples;
	fb->W = W;
	fb->H = H;
	rb->FrameBuffer = fb;
}
void tRdrMakeRenderTiles(tnsRenderBuffer* rb,int TileW,int TileH) {
	tnsRenderTile* rt;
	tnsRenderTile* m;
	tnsFrameBuffer* fb = rb->FrameBuffer;
	int fbw=fb->W, fbh=fb->H;
	char sample = fb->SubPixelSample;
	int Rows = fbh / TileH;
	int Colums = fbw / TileW;
	int r, c;

	//if (!(rb->State & TNS_RENDERBUFFER_RASTERIZER_COMPLETE))return;

	if (Rows*TileH < fbh) 
		Rows += 1;
	if (Colums*TileW < fbw) 
		Colums += 1;

	rb->FrameBuffer->TileCountX = Colums;
	rb->FrameBuffer->TileCountY = Rows;
	rb->FrameBuffer->TileSizeW = TileW;
	rb->FrameBuffer->TileSizeH = TileH;

	rt = CreateNewBuffer(tnsRenderTile, Rows*Colums);
	rb->FrameBuffer->Tiles = rt;
	for (r = 0; r < Rows; r++) {
		for (c = 0; c < Colums; c++) {
			m = &TNS_TILE(rt,r,c,Colums);
			m->Row = r;
			m->Column = c;
			m->FirstPixel = TNS_FRAMEBUFFER_PIXEL((rb->FrameBuffer), r, c);
			m->SubX = c*TileW;
			m->SubY = r*TileH; 
			if (c != Colums - 1)m->SubXLim = m->SubX+TileW; else m->SubXLim = fbw;
			if (r != Rows - 1)m->SubYLim = m->SubY+TileH; else m->SubYLim = fbh;
			m->FX = (real)m->SubX / (real)fbw;
			m->FXLim = (real)m->SubXLim / (real)fbw;
			m->FY = (real)m->SubY / (real)fbh;
			m->FYLim = (real)m->SubYLim / (real)fbh;
			if (c == Colums - 1) m->FXLim = 1.0f;
			if (r == Rows - 1) m->FYLim = 1.0f;
		}
	}
}
void tRdrClearRenderTiles(tnsRenderBuffer* rb, int TileW, int TileH) {
	tnsRenderTile* rt;
	tnsRenderTile* m;
	tnsFrameBuffer* fb = rb->FrameBuffer;
	FreeMem(fb->Tiles);
}


void tRdrCutLineIntegrated(tnsRenderBuffer* rb, tnsRenderLine* rl, real Begin, real End) {
	tnsRenderLineSegment* rls = rl->Segments.pFirst, *irls;
	tnsRenderLineSegment *BeginSegment = 0, *EndSegment = 0;
	tnsRenderLineSegment *ns = 0, *ns2 = 0;
	tnsRenderLineSegment *BeforeBegin, *AfterEnd;
	tnsRenderLineSegment *Next;

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
			ns = memStaticAquireThread(&rb->RenderDataPool,sizeof(tnsRenderLineSegment));
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
			ns2 = memStaticAquireThread(&rb->RenderDataPool,sizeof(tnsRenderLineSegment));
			break;
		}
	}

	if (!ns) ns = memStaticAquireThread(&rb->RenderDataPool,sizeof(tnsRenderLineSegment));
	if (!ns2) ns2 = memStaticAquireThread(&rb->RenderDataPool,sizeof(tnsRenderLineSegment));

	if (BeginSegment) {
		if (BeginSegment != ns) {
			ns->OccludeLevel = BeginSegment->Item.pPrev ? (irls = BeginSegment->Item.pPrev)->OccludeLevel : 0;
			lstInsertItemBefore(&rl->Segments, ns, BeginSegment);
		}
	}
	else {
		ns->OccludeLevel = (irls = rl->Segments.pLast)->OccludeLevel;
		lstAppendItem(&rl->Segments, ns);
	}
	if (EndSegment) {
		if (EndSegment != ns2) {
			ns2->OccludeLevel = EndSegment->Item.pPrev ? (irls = EndSegment->Item.pPrev)->OccludeLevel : 0;
			lstInsertItemBefore(&rl->Segments, ns2, EndSegment);
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

	//for (rls = rl->Segments.pFirst; rls; rls = rls->Item.pNext) {
	//	if ((irls = rls->Item.pNext) && irls->at < rls->at) {
	//		tnsRenderLineSegment* deb;
	//		printf("WRONG!\n");
	//		for (deb = rl->Segments.pFirst; deb; deb = deb->Item.pNext) {
	//			printf("[%d-%d] ", deb->at, deb->OccludeLevel);
	//		}
	//		printf("\n");
	//		break;
	//	}
	//}
}


#define TNS_BOUND_AREA_CROSSES(b1,b2)\
((b1)[0]<(b2)[1] && (b1)[1]>(b2)[0] && (b1)[3]<(b2)[2] && (b1)[2]>(b2)[3])
void tRdrMakeInitialBoundingAreas(tnsRenderBuffer* rb) {
	int SpW = 20;
	int SpH = rb->FrameBuffer->H / (rb->FrameBuffer->W / SpW);
	int Row, Col;
	tnsBoundingArea* ba;
	real W = (real)rb->FrameBuffer->W;
	real H = (real)rb->FrameBuffer->H;
	real SpanW = (real)1 / SpW * 2.0;
	real SpanH = (real)1 / SpH * 2.0;

	rb->FrameBuffer->TileCountX = SpW;
	rb->FrameBuffer->TileCountY = SpH;
	rb->FrameBuffer->WidthPerTile = SpanW;
	rb->FrameBuffer->HeightPerTile = SpanH;

	rb->BoundingAreaCount = SpW * SpH;
	rb->InitialBoundingAreas = memStaticAquire(&rb->RenderDataPool, sizeof(tnsBoundingArea) * rb->BoundingAreaCount);

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
void tRdrConnectNewBoundingAreas(tnsRenderBuffer* rb, tnsBoundingArea* Root) {
	tnsBoundingArea* ba = Root->Child, *tba;
	nListItemPointer* lip,*lip2,*lip3,*NextLip;
	nMemoryPool* mph= &rb->RenderDataPool;

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
		for (lip2 = ((tnsBoundingArea*)lip->p)->RP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((tnsBoundingArea*)lip->p)->RP, lip2);
				if (ba[1].U > tba->B && ba[1].B < tba->U)lstAppendPointerStaticPool(mph, &tba->RP, &ba[1]);
				if (ba[2].U > tba->B && ba[2].B < tba->U)lstAppendPointerStaticPool(mph, &tba->RP, &ba[2]);
			}
		}
	}
	for (lip = Root->RP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((tnsBoundingArea*)lip->p)->LP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((tnsBoundingArea*)lip->p)->LP, lip2);
				if (ba[0].U > tba->B && ba[0].B < tba->U)lstAppendPointerStaticPool(mph, &tba->LP, &ba[0]);
				if (ba[3].U > tba->B && ba[3].B < tba->U)lstAppendPointerStaticPool(mph, &tba->LP, &ba[3]);
			}
		}
	}
	for (lip = Root->UP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((tnsBoundingArea*)lip->p)->BP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((tnsBoundingArea*)lip->p)->BP, lip2);
				if (ba[0].R > tba->L && ba[0].L < tba->R)lstAppendPointerStaticPool(mph, &tba->UP, &ba[0]);
				if (ba[1].R > tba->L && ba[1].L < tba->R)lstAppendPointerStaticPool(mph, &tba->UP, &ba[1]);
			}
		}
	}
	for (lip = Root->BP.pFirst; lip; lip = lip->pNext) {
		for (lip2 = ((tnsBoundingArea*)lip->p)->UP.pFirst; lip2; lip2 = NextLip) {
			NextLip = lip2->pNext;
			tba = lip2->p;
			if (tba == Root) {
				lstRemovePointerItemNoFree(&((tnsBoundingArea*)lip->p)->UP, lip2);
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
void tRdrAssociateTriangleWithBoundingArea(tnsRenderBuffer* rb, tnsBoundingArea* RootBoundingArea, tnsRenderTriangle* rt, real* LRUB, int Recursive);
int tRdrTriangleCalculateIntersectionsInTile(tnsRenderBuffer* rb, tnsRenderTriangle* rt, tnsBoundingArea* ba);

void tRdrSplitBoundingArea(tnsRenderBuffer* rb, tnsBoundingArea* Root) {
	tnsBoundingArea* ba = memStaticAquire(&rb->RenderDataPool, sizeof(tnsBoundingArea) * 4);
	tnsRenderTriangle* rt;

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

	tRdrConnectNewBoundingAreas(rb, Root);

	while (rt = lstPopPointerNoFree(&Root->AssociatedTriangles)) {
		tnsBoundingArea* ba = Root->Child;
		real B[4];
		B[0] = TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
		B[1] = TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
		B[2] = TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
		B[3] = TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[0].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[0], rt, B, 0);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[1].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[1], rt, B, 0);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[2].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[2], rt, B, 0);
		if (TNS_BOUND_AREA_CROSSES(B, &ba[3].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[3], rt, B, 0);
	}

	rb->BoundingAreaCount += 3;
}
int tRdrLineCrossesBoundingArea(tnsFrameBuffer* fb, tnsVector2d L, tnsVector2d R, tnsBoundingArea* ba) {
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
int tRdrTriangleCoversBoundingArea(tnsFrameBuffer* fb, tnsRenderTriangle* rt, tnsBoundingArea* ba) {

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

	if (tRdrPointInsideTrianglef(p1, FBC1, FBC2, FBC3) ||
		tRdrPointInsideTrianglef(p2, FBC1, FBC2, FBC3) ||
		tRdrPointInsideTrianglef(p3, FBC1, FBC2, FBC3) ||
		tRdrPointInsideTrianglef(p4, FBC1, FBC2, FBC3)) return 1;

	if (tRdrLineCrossesBoundingArea(fb, FBC1, FBC2, ba)) return 1;
	elif(tRdrLineCrossesBoundingArea(fb, FBC2, FBC3, ba)) return 1;
	elif(tRdrLineCrossesBoundingArea(fb, FBC3, FBC1, ba)) return 1;

	return 0;
}
void tRdrAssociateTriangleWithBoundingArea(tnsRenderBuffer* rb, tnsBoundingArea* RootBoundingArea, tnsRenderTriangle* rt, real* LRUB, int Recursive) {
	if (!tRdrTriangleCoversBoundingArea(rb->FrameBuffer, rt, RootBoundingArea)) return;
	if (!RootBoundingArea->Child) {
		lstAppendPointerStaticPool(&rb->RenderDataPool, &RootBoundingArea->AssociatedTriangles, rt);
		RootBoundingArea->TriangleCount++;
		if (RootBoundingArea->TriangleCount > 200 && Recursive) {
			tRdrSplitBoundingArea(rb, RootBoundingArea);
		}
		if(Recursive) tRdrTriangleCalculateIntersectionsInTile(rb, rt, RootBoundingArea);
	}else {
		tnsBoundingArea* ba = RootBoundingArea->Child;
		real* B1 = LRUB;
		real B[4];
		if (!LRUB) {
			B[0] = TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
			B[1] = TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
			B[2] = TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
			B[3] = TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
			B1 = B;
		}
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[0].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[0], rt, B1, Recursive);
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[1].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[1], rt, B1, Recursive);
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[2].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[2], rt, B1, Recursive);
		if (TNS_BOUND_AREA_CROSSES(B1, &ba[3].L)) tRdrAssociateTriangleWithBoundingArea(rb, &ba[3], rt, B1, Recursive);
	}
}
int tRdrGetTriangleBoundingTile(tnsRenderBuffer* rb, tnsRenderTriangle* rt,int* RowBegin,int* RowEnd,int* ColBegin,int* ColEnd) {
	real SpW = rb->FrameBuffer->WidthPerTile, SpH = rb->FrameBuffer->HeightPerTile;
	real B[4];

	if (!rt->F || !rt->V[0] || !rt->V[1] || !rt->V[2]) return 0;

	B[0] = TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
	B[1] = TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]);
	B[2] = TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);
	B[3] = TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]);

	if (B[0] > 1 || B[1] < -1 || B[2]>1 || B[3] < -1) return 0;

	(*ColBegin) = (int)((B[0] + 1.0) / SpW);
	(*ColEnd) = (int)((B[1] + 1.0) / SpW);
	(*RowEnd) = rb->FrameBuffer->TileCountY - (int)((B[2] + 1.0) / SpH) - 1;
	(*RowBegin) = rb->FrameBuffer->TileCountY - (int)((B[3] + 1.0) / SpH) - 1;

	if ((*ColEnd) >= rb->FrameBuffer->TileCountX) (*ColEnd) = rb->FrameBuffer->TileCountX-1;
	if ((*RowEnd) >= rb->FrameBuffer->TileCountY) (*RowEnd) = rb->FrameBuffer->TileCountY-1;
	if ((*ColBegin) < 0) (*ColBegin) = 0;
	if ((*RowBegin) < 0) (*RowBegin) = 0;

	return 1;
}
void tRdrAddTriangles(tnsRenderBuffer* rb) {
	tnsRenderElementLinkNode* reln;
	tnsRenderTriangle* rt;
	tnsFrameBuffer* fb = rb->FrameBuffer;
	tnsRenderTile* tile;
	tnsMatrix44d VP;
	tnsCamera* c = ((tnsCamera*)rb->Scene->ActiveCamera);
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

	tnsglobal_TriangleIntersectionCount = 0;

	tnsset_RenderOverallProgress(rb, NUL_MH2);
	rb->CalculationStatus = TNS_CALCULATION_INTERSECTION;
	nulThreadNotifyUsers("tns.render_buffer_list.calculation_status");

	for (reln = rb->TriangleBufferPointers.pFirst; reln; reln = reln->Item.pNext) {
		rt = reln->Pointer;
		lim = reln->ElementCount;
		for (i = 0; i < lim; i++) {
			if (rt->CullStatus) {
				rt++; continue;
			}
			if (tRdrGetTriangleBoundingTile(rb, rt, &y1, &y2, &x1, &x2)) {
				for (co = x1; co <= x2; co++) {
					for (r = y1; r <= y2; r++) {
						tRdrAssociateTriangleWithBoundingArea(rb, &rb->InitialBoundingAreas[r * 20 + co], rt, 0, 1);
					}
				}
			}else{
				;//throw away.
			}
			rt = ((BYTE*)rt) + rb->TriangleSize;
			if (tnsglobal_TriangleIntersectionCount >= 2000) {
				tnsset_PlusRenderIntersectionCount(rb, tnsglobal_TriangleIntersectionCount);
				tnsglobal_TriangleIntersectionCount = 0;
			}
		}
	}
	tnsset_PlusRenderIntersectionCount(rb, tnsglobal_TriangleIntersectionCount);
}
tnsBoundingArea* tRdrGetNextBoundingArea(tnsBoundingArea* This, tnsRenderLine* rl, real x, real y, real k, int PositiveX, int PositiveY, real* NextX, real* NextY) {
	real rx, ry, ux, uy, lx, ly, bx, by;
	real r1, r2, r;
	tnsBoundingArea* ba; nListItemPointer* lip;
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

tnsBoundingArea* tRdrGetBoundingArea(tnsRenderBuffer* rb, real x, real y) {
	tnsBoundingArea* iba;
	real SpW = rb->FrameBuffer->WidthPerTile, SpH = rb->FrameBuffer->HeightPerTile;
	int c = (int)((x + 1.0) / SpW);
	int r = rb->FrameBuffer->TileCountY - (int)((y + 1.0) / SpH) - 1;
	if (r < 0) r = 0;
	if (c < 0) c = 0;
	if (r >= rb->FrameBuffer->TileCountY) r = rb->FrameBuffer->TileCountY-1;
	if (c >= rb->FrameBuffer->TileCountX) c = rb->FrameBuffer->TileCountX-1;

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
tnsBoundingArea* tRdrGetFirstPossibleBoundingArea(tnsRenderBuffer* rb, tnsRenderLine* rl) {
	tnsBoundingArea* iba;
	real p[2] = { rl->L->FrameBufferCoord[0], rl->L->FrameBufferCoord[1] };
	tnsVector2d LU = { -1,1 }, RU = { 1,1 }, LB = { -1,-1 }, RB = { 1, -1 };
	real r=1,sr=1;

	if (p[0] > -1 && p[0]<1 && p[1] > -1 && p[1]<1) {
		return tRdrGetBoundingArea(rb, p[0], p[1]);
	}else {
		if (tRdrLineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, LU, RU, &sr) && sr<r && sr>0)r = sr;
		if (tRdrLineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, LB, RB, &sr) && sr<r && sr>0)r = sr;
		if (tRdrLineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, LB, LU, &sr) && sr<r && sr>0)r = sr;
		if (tRdrLineIntersectTest2d(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, RB, RU, &sr) && sr<r && sr>0)r = sr;
		tnsLinearInterpolate2dv(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, r, p);

		return tRdrGetBoundingArea(rb, p[0], p[1]);
	}

	real SpW = rb->FrameBuffer->WidthPerTile, SpH = rb->FrameBuffer->HeightPerTile;

	return iba;
}
int tRdrTriangleLineImageSpaceIntersectTestOnly(tnsRenderTriangle* rt, tnsRenderLine* rl, double* From, double* To);
int tRdrTriangleLineImageSpaceIntersectTestOnlyV2(tnsRenderTriangle* rt, tnsRenderLine* rl, tnsCamera* cam, tnsMatrix44d vp, double* From, double* To);

int tRdrMakeNextOcclusionTaskInfo(tnsRenderBuffer* rb, tnsRenderTaskInfo* rti) {
	nListItemPointer* p;
	int i;
	int res = 0;

	EnterCriticalSection(&rb->csManagement);

	if (rb->ContourManaged) {
		p = rb->ContourManaged;
		rti->Contour = p;
		rti->ContourPointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->ContourManaged = p;
		rti->ContourPointers.pLast = p ? p->pPrev : rb->Contours.pLast;
		res = 1;
	}else {
		lstEmptyDirect(&rti->ContourPointers);
		rti->Contour = 0;
	}

	if (rb->IntersectionManaged) {
		p = rb->IntersectionManaged;
		rti->Intersection = p;
		rti->IntersectionPointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->IntersectionManaged = p;
		rti->IntersectionPointers.pLast = p ? p->pPrev : rb->IntersectionLines.pLast;
		res = 1;
	}else {
		lstEmptyDirect(&rti->IntersectionPointers);
		rti->Intersection = 0;
	}

	if (rb->CreaseManaged) {
		p = rb->CreaseManaged;
		rti->Crease = p;
		rti->CreasePointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->CreaseManaged = p;
		rti->CreasePointers.pLast = p ? p->pPrev : rb->CreaseLines.pLast;
		res = 1;
	}else {
		lstEmptyDirect(&rti->CreasePointers);
		rti->Crease = 0;
	}

	if (rb->MaterialManaged) {
		p = rb->MaterialManaged;
		rti->Material = p;
		rti->MaterialPointers.pFirst = p;
		for (i = 0; i < TNS_THREAD_LINE_COUNT && p; i++) {
			p = p->pNext;
		}
		rb->MaterialManaged = p;
		rti->MaterialPointers.pLast = p ? p->pPrev : rb->MaterialLines.pLast;
		res = 1;
	}else {
		lstEmptyDirect(&rti->MaterialPointers);
		rti->Material = 0;
	}

	LeaveCriticalSection(&rb->csManagement);

	return res;
}
void tRdrCalculateSingleLineOcclusion(tnsRenderBuffer* rb, tnsRenderLine* rl, int ThreadID) {
	real x = rl->L->FrameBufferCoord[0], y = rl->L->FrameBufferCoord[1];
	tnsBoundingArea* ba = tRdrGetFirstPossibleBoundingArea(rb, rl);
	tnsBoundingArea* nba = ba;
	tnsRenderTriangleThread* rt;
	nListItemPointer* lip;
	tnsRenderVert* rv;
	tnsCamera* c = rb->Scene->ActiveCamera;
	real l, r;
	real k = (rl->R->FrameBufferCoord[1] - rl->L->FrameBufferCoord[1]) / (rl->R->FrameBufferCoord[0] - rl->L->FrameBufferCoord[0] + 1e-30);
	int PositiveX = (rl->R->FrameBufferCoord[0] - rl->L->FrameBufferCoord[0])>0 ? 1 : 0;
	int PositiveY = (rl->R->FrameBufferCoord[1] - rl->L->FrameBufferCoord[1])>0 ? 1 : 0;

	while (nba) {

		for (lip = nba->AssociatedTriangles.pFirst; lip; lip = lip->pNext) {
			rt = lip->p;
			if (rt->Testing[ThreadID] == rl || rl->L->IntersectWith == rt || rl->R->IntersectWith == rt) continue;
			rt->Testing[ThreadID] = rl;
			if (tRdrTriangleLineImageSpaceIntersectTestOnlyV2(rt, rl,c,rb->FrameBuffer->ViewProjection, &l, &r)) {
			    tRdrCutLineIntegrated(rb, rl, l, r);
			}
		}
		nba = tRdrGetNextBoundingArea(nba, rl, x, y, k, PositiveX, PositiveY, &x, &y);
	}
}
void THREAD_CalculateLineOcclusion(tnsRenderTaskInfo* rti) {
	tnsRenderBuffer* rb = rti->RenderBuffer;
	int ThreadId = rti->ThreadID;
	tnsRenderLine* rl;
	nListItemPointer* lip;
	int count=0;

	while (tRdrMakeNextOcclusionTaskInfo(rb, rti)) {

		for (lip = rti->Contour; lip&& lip->pPrev != rti->ContourPointers.pLast; lip = lip->pNext) {
			tRdrCalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);

			count++;
		}
		tnsset_PlusRenderContourProcessedCount(rb, count);
		count = 0;

		for (lip = rti->Crease; lip && lip->pPrev != rti->CreasePointers.pLast; lip = lip->pNext) {
			tRdrCalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);
			count++;
		}
		tnsset_PlusRenderCreaseProcessedCount(rb, count);
		count = 0;

		for (lip = rti->Intersection; lip&& lip->pPrev != rti->IntersectionPointers.pLast; lip = lip->pNext) {
			tRdrCalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);
			count++;
		}
		tnsset_PlusRenderIntersectionProcessedCount(rb, count);
		count = 0;

		for (lip = rti->Material; lip&& lip->pPrev != rti->MaterialPointers.pLast; lip = lip->pNext) {
			tRdrCalculateSingleLineOcclusion(rb, lip->p, rti->ThreadID);
			count++;
		}
		tnsset_PlusRenderMaterialProcessedCount(rb, count);
		count = 0;

	}
	thrd_exit(0);
}

int tRdrGetNormal(tnsVector3d v1, tnsVector3d v2, tnsVector3d v3, tnsVector3d n, tnsVector3d Pos) {
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

tnsRenderTile* tRdrFrameBufferCoordToTile(tnsFrameBuffer* fb, tnsVector2d CoordXY) {
	int r, c;
	int Row, Column;
	tnsVector2d Converted = { 0 };
	tnsRenderTile* rt = 0;

	tMatVectorCopy2d(CoordXY, Converted);

	Converted[0] += 1;
	Converted[1] += 1;
	Converted[0] /= 2;
	Converted[1] /= 2;

	if (fb->WidthPerTile < 0.00001) {
		fb->WidthPerTile = ((real)fb->TileSizeW / (real)fb->W);
		fb->HeightPerTile = ((real)fb->TileSizeH / (real)fb->H);
	}

	r = (int)(((real)Converted[1] / (real)fb->H) / (real)fb->HeightPerTile);
	c = (int)(((real)Converted[0] / (real)fb->W) / (real)fb->WidthPerTile);

	if (r < 0 || r >= fb->TileCountY) return 0;
	if (c < 0 || c >= fb->TileCountX) return 0;

	rt = &TNS_TILE(fb->Tiles, r, c, fb->TileCountX);

	Row = r;
	Column = c;

	if (TNS_IN_TILE_X(rt,(Converted[0]/fb->W))) {
		if (TNS_IN_TILE_Y(rt, (Converted[1]/fb->H))) {
		}else {
			if (r < fb->TileCountY - 1) Row = r + 1;
		}
	}else {
		if (TNS_IN_TILE_Y(rt, (Converted[1] / fb->H))) {
			if (c < fb->TileCountX - 1) Column = c + 1;
		}
		else {
			if (r < fb->TileCountY - 1) Row = r + 1;
			if (c < fb->TileCountX - 1) Column = c + 1;
		}
	}

	rt = &TNS_TILE(fb->Tiles, r, c, fb->TileCountX);
}
int tRdrPointExtendTouchesTileEdge(tnsFrameBuffer* fb, tnsRenderTile* FromTile, tnsVector2d FromPoint,tnsVector2d ToPoint, int* L, int* R, int* T, int* B,tnsVector2d CrossPoint) {
	tnsVector2d Vector;
	real tilt;

	tMatVectorMinus2d(Vector, ToPoint, FromPoint);

	tilt = Vector[1] / Vector[0];

	if (Vector[0]>0) {
		real YC = (FromTile->FXLim - FromPoint[0]) * tilt + FromPoint[1];
		if ((Vector[1]>0 && YC <= ToPoint[1])|| (Vector[1]<0 && YC >= ToPoint[1])) {
			if (TNS_IN_TILE_Y(FromTile, YC)) {
				*R = 1;
			}elif(Vector[1]>0) {
				real XC = (FromTile->FYLim - FromPoint[1]) / fabsf(tilt) + FromPoint[0];
				if (TNS_IN_TILE_X(FromTile, XC)) {
					*T = 1;
				}
			}else {
				real XC = (FromPoint[1] - FromTile->FY) / fabsf(tilt) + FromPoint[0];
				if (TNS_IN_TILE_X(FromTile, XC)) {
					*B = 1;
				}
			}
		}
		else return 0;
	}else {
		real YC = ( FromPoint[0] - FromTile->FX) * tilt + FromPoint[1];
		if ((Vector[1]>0 && YC <= ToPoint[1]) || (Vector[1]<0 && YC >= ToPoint[1])) {
			if (TNS_IN_TILE_Y(FromTile, YC)) {
				*L = 1;
			}elif(Vector[1] > 0) {
				real XC = FromPoint[0] - (FromTile->FYLim - FromPoint[1]) / fabsf(tilt);
				if (TNS_IN_TILE_X(FromTile, XC)) {
					*T = 1;
				}
			}else {
				real XC = FromPoint[0] - (FromPoint[1] - FromTile->FY) / fabsf(tilt);
				if (TNS_IN_TILE_X(FromTile, XC)) {
					*B = 1;
				}
			}
		}else return 0;
	}
	return 1;
}
int tRdrBoundBoxCrosses(tnsVector4d xxyy1, tnsVector4d xxyy2) {
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
int tRdrLineCrossesTile(tnsFrameBuffer* fb, tnsVector2d L, tnsVector2d R, tnsRenderTile* Tile) {
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

	//if (!tRdrBoundBoxCrosses(Converted, LineBoundBox))return 0;

	//tMatVectorMinus2d(vec, L, R);
	vx = L[0] - R[0];
	vy = L[1] - R[1];

	c1 = vx * (Converted[2]-L[1]) - vy * (Converted[0] -L[0]);
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
int tRdrPointInsideTrianglef(tnsVector2d v, tnsVector2d v0,tnsVector2d v1, tnsVector2d v2) {
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
int tRdrPointOnLinef(tnsVector2d v, tnsVector2d v0, tnsVector2d v1) {
	real c1, c2;

	c1 = tMatGetLinearRatio(v0[0], v1[0], v[0]);
	c2 = tMatGetLinearRatio(v0[1], v1[1], v[1]);

	if (TNS_DOUBLE_CLOSE_ENOUGH(c1, c2) && c1>=0 &&c1<=1) return 1;

	return 0;
}
int tRdrPointTriangleRelation(tnsVector2d v, tnsVector2d v0, tnsVector2d v1, tnsVector2d v2) {
	double cl, c;
	real r;
	if (tRdrPointOnLinef(v, v0, v1) || tRdrPointOnLinef(v, v1, v2) || tRdrPointOnLinef(v, v2, v0)) return 1;

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
int tRdrPointInsideTriangle3d(tnsVector3d v, tnsVector3d v0, tnsVector3d v1, tnsVector3d v2) {
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
int tRdrPointInsideTriangle3de(tnsVector3d v, tnsVector3d v0, tnsVector3d v1, tnsVector3d v2) {
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

	if ((d=tMatDot3d(N1, N2, 0)) < 0) return 0;
	//if (d<DBL_EPSILON) return -1;

	tMatVectorMinus3d(L, v1, v0);
	tMatVectorMinus3d(R, v, v1);
	//tMatNormalizeSelf3d(L);
	//tMatNormalizeSelf3d(R);
	tMatVectorCross3d(N2, L, R);

	if ((d=tMatDot3d(N1, N2, 0)) < 0) return 0;
	//if (d<DBL_EPSILON) return -1;

	return 1;
}

int tRdrTileInsideTriangle(tnsFrameBuffer* fb, tnsRenderTriangle* rt, tnsRenderTile* Tile) {
	tnsVector2d Converted = { 0 };

	Converted[0] = (real)Tile->SubX - (real)fb->W / 2;
	Converted[1] = (real)Tile->SubY - (real)fb->H / 2;
	if (!tRdrPointInsideTrianglef(Converted, rt->V[0]->FrameBufferCoord, rt->V[1]->FrameBufferCoord, rt->V[2]->FrameBufferCoord))return 0;

	//Converted[0] = (real)Tile->SubX - (real)fb->W / 2;
	Converted[1] = (real)Tile->SubYLim - (real)fb->H / 2;
	if (!tRdrPointInsideTrianglef(Converted, rt->V[0]->FrameBufferCoord, rt->V[1]->FrameBufferCoord, rt->V[2]->FrameBufferCoord))return 0;

	Converted[0] = (real)Tile->SubXLim - (real)fb->W / 2;
	//Converted[1] = (real)Tile->SubYLim - (real)fb->H / 2;
	if (!tRdrPointInsideTrianglef(Converted, rt->V[0]->FrameBufferCoord, rt->V[1]->FrameBufferCoord, rt->V[2]->FrameBufferCoord))return 0;

	//Converted[0] = (real)Tile->SubXLim - (real)fb->W / 2;
	Converted[1] = (real)Tile->SubY - (real)fb->H / 2;
	if (!tRdrPointInsideTrianglef(Converted, rt->V[0]->FrameBufferCoord, rt->V[1]->FrameBufferCoord, rt->V[2]->FrameBufferCoord))return 0;

	return 1;
}
int tRdrTriangleCoversTile(tnsFrameBuffer* fb, tnsRenderTriangle* rt, tnsRenderTile* Tile) {

	real vx, vy;
	tnsVector2d p1,p2,p3,p4;
	real
		*FBC1 = rt->V[0]->FrameBufferCoord,
		*FBC2 = rt->V[1]->FrameBufferCoord,
		*FBC3 = rt->V[2]->FrameBufferCoord;

	p3[0] = p1[0] = (real)Tile->SubX - (real)fb->W / 2;
	p2[1] = p1[1] = (real)Tile->SubY - (real)fb->H / 2;
	p2[0] = p4[0] = (real)Tile->SubXLim - (real)fb->W / 2;
	p3[1] = p4[1] = (real)Tile->SubYLim - (real)fb->H / 2;

	if (FBC1[0] >= p1[0] && FBC1[0]<=p2[0] && FBC1[1]>=p1[1] && FBC1[1] <= p3[1]) return 1;
	if (FBC2[0] >= p1[0] && FBC2[0]<=p2[0] && FBC2[1]>=p1[1] && FBC2[1] <= p3[1]) return 1;
	if (FBC3[0] >= p1[0] && FBC3[0]<=p2[0] && FBC3[1]>=p1[1] && FBC3[1] <= p3[1]) return 1;

	if (tRdrPointInsideTrianglef(p1, FBC1, FBC2, FBC3) ||
		tRdrPointInsideTrianglef(p2, FBC1, FBC2, FBC3) ||
		tRdrPointInsideTrianglef(p3, FBC1, FBC2, FBC3) ||
		tRdrPointInsideTrianglef(p4, FBC1, FBC2, FBC3)) return 1;

	if (tRdrLineCrossesTile(fb, FBC1, FBC2, Tile)) return 1;
	elif (tRdrLineCrossesTile(fb, FBC2, FBC3, Tile)) return 1;
	elif (tRdrLineCrossesTile(fb, FBC3, FBC1, Tile)) return 1;


	return 0;
}
int tRdrGetTriangleBoundTile(tnsFrameBuffer* fb, tnsRenderTriangle* triangle, int* XMin, int* XMax, int* YMin, int* YMax) {
	tnsRenderTile* rt;
	int x1, x2, x3, y1, y2, y3;
	int Row, Col;
	int t;

	rt = tRdrFrameBufferCoordToTile(fb, triangle->V[0]->FrameBufferCoord);
	if (!rt) return 0;
	x1 = rt->Column;
	y1 = rt->Row;

	rt = tRdrFrameBufferCoordToTile(fb, triangle->V[1]->FrameBufferCoord);
	if (!rt) return 0;
	x2 = rt->Column;
	y2 = rt->Row;

	rt = tRdrFrameBufferCoordToTile(fb, triangle->V[2]->FrameBufferCoord);
	if (!rt) return 0;
	x3 = rt->Column;
	y3 = rt->Row;

	*XMin = TNS_MIN3(x1, x2, x3);
	*XMax = TNS_MAX3(x1, x2, x3);
	*YMin = TNS_MIN3(y1, y2, y3);
	*YMax = TNS_MAX3(y1, y2, y3);

	return 1;
}

tnsRenderElementLinkNode* tRdrNewCullTriangleSpace64(tnsRenderBuffer* rb) {
	tnsRenderElementLinkNode* reln;

	tnsRenderTriangle* RenderTriangles = calloc(64, rb->TriangleSize);//CreateNewBuffer(tnsRenderTriangle, 64);

	reln = lstAppendPointerStaticSized(&rb->TriangleBufferPointers, &rb->RenderDataPool, RenderTriangles,
		sizeof(tnsRenderElementLinkNode));
	reln->ElementCount = 64;
	reln->Additional = 1;

	return reln;
}
tnsRenderElementLinkNode* tRdrNewCullPointSpace64(tnsRenderBuffer* rb) {
	tnsRenderElementLinkNode* reln;

	tnsRenderVert* RenderVertices = CreateNewBuffer(tnsRenderVert, 64);

	reln = lstAppendPointerStaticSized(&rb->VertexBufferPointers, &rb->RenderDataPool, RenderVertices,
		sizeof(tnsRenderElementLinkNode));
	reln->ElementCount = 64;
	reln->Additional = 1;

	return reln;
}
void tRdrCalculateRenderTriangleNormal(tnsRenderTriangle* rt);
void tRdrPostTriangle(tnsRenderTriangle* rt, tnsRenderTriangle* orig) {
	if (rt->V[0])tMatVectorAccum3d(rt->GC, rt->V[0]->FrameBufferCoord);
	if (rt->V[1])tMatVectorAccum3d(rt->GC, rt->V[1]->FrameBufferCoord);
	if (rt->V[2])tMatVectorAccum3d(rt->GC, rt->V[2]->FrameBufferCoord);
	tMatVectorMultiSelf3d(rt->GC, 1.0f / 3.0f);

	tMatVectorCopy3d(orig->GN, rt->GN);
}
void tRdrCullTriangles(tnsRenderBuffer* rb) {
	tnsRenderLine* rl;
	tnsRenderTriangle* rt, *rt1;
	tnsRenderVert* rv;
	tnsRenderElementLinkNode* reln,*veln,*teln;
	tnsRenderLineSegment* rls;
	tnsVector4d p1, p2;
	real* MVInverse = rb->FrameBuffer->VPInverse;
	int i;
	real a;
	int VCount = 0, TCount = 0;
	tns3DObject* o;

	veln = tRdrNewCullPointSpace64(rb);
	teln = tRdrNewCullTriangleSpace64(rb);
	rv = &((tnsRenderVert*)veln->Pointer)[VCount];
	rt1 = ((BYTE*)teln->Pointer) + rb->TriangleSize*TCount;

	for (reln = rb->TriangleBufferPointers.pFirst; reln; reln = reln->Item.pNext) {
		i = 0;
		if (reln->Additional) continue;
		o = reln->ObjectRef;
		for (i; i < reln->ElementCount; i++) {
			int In1 = 0, In2 = 0, In3 = 0;
			rt = ((BYTE*)reln->Pointer) + rb->TriangleSize*i;
			if (rt->V[0]->FrameBufferCoord[3] < 0) In1 = 1;
			if (rt->V[1]->FrameBufferCoord[3] < 0) In2 = 1;
			if (rt->V[2]->FrameBufferCoord[3] < 0) In3 = 1;

			rt->RL[0]->ObjectRef = o;
			rt->RL[1]->ObjectRef = o;
			rt->RL[2]->ObjectRef = o;

			if (VCount > 60) {
				veln->ElementCount = VCount;
				veln = tRdrNewCullPointSpace64(rb);
				VCount = 0;
			}

			if (TCount > 60) {
				teln->ElementCount = TCount;
				teln = tRdrNewCullTriangleSpace64(rb);
				TCount = 0;
			}

			if ((!rt->RL[0]->Item.pNext && !rt->RL[0]->Item.pPrev)||
				(!rt->RL[1]->Item.pNext && !rt->RL[1]->Item.pPrev)||
				(!rt->RL[2]->Item.pNext && !rt->RL[2]->Item.pPrev)) {
				printf("'");
			}

			rv = &((tnsRenderVert*)veln->Pointer)[VCount];
			rt1 = &((tnsRenderTriangle*)teln->Pointer)[TCount];

			
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

					lstRemoveItem(&rb->AllRenderLines, rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[2]);
					
					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1->RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = rt->V[0];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1->RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					tRdrPostTriangle(o,rt1,rt);

					VCount += 2;
					TCount += 1;
					continue;
				}elif (!In3) {
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

					lstRemoveItem(&rb->AllRenderLines, rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = &rv[1];
					rl->TL = rt1;
					rt1->RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = rt->V[2];
					rl->TL = rt->RL[1]->TL == rt ? rt1 : rt->RL[1]->TL;
					rl->TR = rt->RL[1]->TR == rt ? rt1 : rt->RL[1]->TR;
					rt1->RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					tRdrPostTriangle(rt1,rt);

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

					lstRemoveItem(&rb->AllRenderLines, rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1->RL[2] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[1];
					rl->TL = rt->RL[0]->TL == rt ? rt1 : rt->RL[0]->TL;
					rl->TR = rt->RL[0]->TR == rt ? rt1 : rt->RL[0]->TR;
					rt1->RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					tRdrPostTriangle(rt1, rt);

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

					lstRemoveItem(&rb->AllRenderLines, rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1[0].RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[1];
					rl->TL = &rt1[0];
					rl->TR = &rt1[1];
					rt1[0].RL[2] = rl;
					rt1[1].RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					tRdrPostTriangle(&rt1[0],rt);
					tRdrPostTriangle(&rt1[1],rt);

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

					lstRemoveItem(&rb->AllRenderLines, rt->RL[0]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[1]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1[0].RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[2];
					rl->TL = &rt1[0];
					rl->TR = &rt1[1];
					rt1[0].RL[2] = rl;
					rt1[1].RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					tRdrPostTriangle(&rt1[0], rt);
					tRdrPostTriangle(&rt1[1], rt);

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

					lstRemoveItem(&rb->AllRenderLines, rt->RL[1]);
					lstRemoveItem(&rb->AllRenderLines, rt->RL[2]);

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[1];
					rl->R = &rv[0];
					rl->TL = rt1;
					rt1[0].RL[1] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
					lstAppendItem(&rl->Segments, rls);
					lstAppendItem(&rb->AllRenderLines, rl);
					rl->L = &rv[0];
					rl->R = rt->V[0];
					rl->TL = &rt1[0];
					rl->TR = &rt1[1];
					rt1[0].RL[2] = rl;
					rt1[1].RL[0] = rl;

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					rl = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
					rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
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

					tRdrPostTriangle(&rt1[0], rt);
					tRdrPostTriangle(&rt1[1], rt);

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
void tRdrPerspectiveDivision(tnsRenderBuffer* rb) {
	tnsRenderVert* rv;
	tnsRenderElementLinkNode* reln;
	tnsCamera* cam = rb->Scene->ActiveCamera;
	int i;

	if (cam->CameraType == TNS_CAMERA_ORTHO) return;

	for (reln = rb->VertexBufferPointers.pFirst; reln; reln = reln->Item.pNext) {
		rv = reln->Pointer;
		for (i = 0; i < reln->ElementCount; i++) {
			//if (rv->FrameBufferCoord[2] < -DBL_EPSILON) continue;
			tMatVectorMultiSelf3d(rv[i].FrameBufferCoord, 1/rv[i].FrameBufferCoord[3]);
			rv[i].FrameBufferCoord[2] = cam->ZMin*cam->ZMax / (cam->ZMax - fabs(rv[i].FrameBufferCoord[2]) * (cam->ZMax - cam->ZMin));
		}
	}
}

void tRdrTransformRenderVert(tnsVert* V, tnsRenderVert* RVBuf, real* MVMat, real* MVPMat, tnsCamera* Camera){//real HeightMultiply, real ZMin, real ZMax) {
	tnsRenderVert* rv = &RVBuf[V->I];
	rv->V = V;
	V->RV = rv;
	tMatApplyTransform43d(rv->GLocation, MVMat, V->P);
	tMatApplyTransform44d(rv->FrameBufferCoord, MVPMat, V->P);

	//if(rv->FrameBufferCoord[2]>0)tMatVectorMultiSelf3d(rv->FrameBufferCoord, (1 / rv->FrameBufferCoord[3]));
	//else tMatVectorMultiSelf3d(rv->FrameBufferCoord, -rv->FrameBufferCoord[3]);
 //   rv->FrameBufferCoord[2] = Camera->ZMin* Camera->ZMax / (Camera->ZMax - fabs(rv->FrameBufferCoord[2]) * (Camera->ZMax - Camera->ZMin));
}
void tRdrCalculateRenderTriangleNormal(tnsRenderTriangle* rt) {
	tnsVector3d L, R;
	tMatVectorMinus3d(L, rt->V[1]->GLocation , rt->V[0]->GLocation);
	tMatVectorMinus3d(R, rt->V[2]->GLocation , rt->V[0]->GLocation);
	tMatVectorCross3d(rt->GN, L, R);
	tMatNormalizeSelf3d(rt->GN);
}
void tRdrMakeRenderGeometryBuffersRecursive(tns3DObject* o,real* MVMat,real* MVPMat,tnsRenderBuffer* rb, real HeightMultiply){
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
void tnsMakeRenderGeometryBuffers(tnsScene* s,tnsCamera* c, tnsRenderBuffer* rb, int HeightMultiply) {
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
		tObjApplyGlobalTransformMatrixRecursive(o);
		tRdrMakeRenderGeometryBuffersRecursive(o, view, proj, rb, (real)HeightMultiply);
	}
}
void tRdrZeroGeomtryBuffersRecursive(tns3DObject* o) {
	tns3DObject* oc;
	o->RenderTriangles = 0;
	o->RenderVertices = 0;
	for (oc = o->ChildObjects.pFirst; oc; oc = oc->Item.pNext) {
		tRdrZeroGeomtryBuffersRecursive(oc);
	}
}
void tnsZeroGeomtryBuffers(tnsScene* s) {
	tns3DObject* o;
	tns3DObject* oc;
	if (!s) return;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		tRdrZeroGeomtryBuffersRecursive(o);
	}
}

void tnsCreateRasterizerbuffer(tnsRenderBuffer* rb, int W, int H, int Samples, int TileW, int TileH) {
	tRdrMakeFrameBuffer(rb, W, H, Samples);
	tRdrMakeRenderTiles(rb, TileW, TileH);	
}
void tnsTransformGeomtryToRenderBuffer(tnsRenderBuffer* rb) {
	//if (rb->State&TNS_RENDERBUFFER_GEOMETRY_COMPLETE) return;

	if (!rb->Scene || !rb->Scene->ActiveCamera || !rb->FrameBuffer) return;

	tnsset_RenderOverallProgress(rb, 10);
	rb->CalculationStatus = TNS_CALCULATION_GEOMETRY;
	nulThreadNotifyUsers("tns.render_buffer_list.calculation_status");

	tnsMakeRenderGeometryBuffers(rb->Scene, rb->Scene->ActiveCamera, rb, rb->FrameBuffer->H);
}
tnsRenderBuffer* tnsCreateRenderBuffer(char* name) {
	tnsRenderBuffer* rb = memAquireHyper(sizeof(tnsRenderBuffer));

	rb->Scene = T->World.ActiveScene;

	strSafeSet(&rb->Name, name);

	lstAppendItem(&T->RenderBuffers, rb);

	return rb;
}

#define INTERSECT_SORT_MIN_TO_MAX_3(ia,ib,ic,lst)\
{\
lst[0] = TNS_MIN3_INDEX(ia,ib,ic);\
lst[1] = ( ((ia<=ib && ib<=ic)||(ic<=ib && ib<=ia)) ? 1 : (((ic<=ia && ia<=ib)||(ib<ia && ia<=ic))? 0 : 2));\
lst[2] = TNS_MAX3_INDEX(ia,ib,ic);\
}

//ia ib ic are ordered
#define INTERSECT_JUST_GREATER(is,order,num,index)\
{\
index = (num<is[order[0]]?order[0]:(num<is[order[1]]?order[1]:(num<is[order[2]]?order[2]:order[2])));\
}

//ia ib ic are ordered
#define INTERSECT_JUST_SMALLER(is,order,num,index)\
{\
index = (num>is[order[2]]?order[2]:(num>is[order[1]]?order[1]:(num>is[order[0]]?order[0]:order[0])));\
}

void tRdrGetInterpolatePoint2d(tnsVector2d L, tnsVector2d R, real Ratio, tnsVector2d Result) {
	Result[0] = tnsLinearItp(L[0], R[0], Ratio);
	Result[1] = tnsLinearItp(L[1], R[1], Ratio);
}
void tRdrGetInterpolatePoint3d(tnsVector2d L, tnsVector2d R, real Ratio, tnsVector2d Result) {
	Result[0] = tnsLinearItp(L[0], R[0], Ratio);
	Result[1] = tnsLinearItp(L[1], R[1], Ratio);
	Result[2] = tnsLinearItp(L[2], R[2], Ratio);
}

int tRdtGetZIntersectPoint(tnsVector3d TL, tnsVector3d TR, tnsVector3d LL, tnsVector3d LR, tnsVector3d IntersectResult) {
	//real lzl = 1 / LL[2], lzr = 1 / LR[2], tzl = 1 / TL[2], tzr = 1 / TR[2];
	real lzl = LL[2], lzr = LR[2], tzl = TL[2], tzr = TR[2];
	real l = tzl - lzl, r = tzr - lzr;
	real ratio;
	int rev = l < 0 ? -1 : 1;//-1:occlude left 1:occlude right

	if (l*r >= 0) {
		if (l == 0) IntersectResult[2] = r > 0 ? -1 : 1;
		else if (r == 0) IntersectResult[2] = l > 0 ? -1 : 1;
		else 
		IntersectResult[2] = rev;
		return 0;
	}
	l = fabsf(l);
	r = fabsf(r);
	ratio = l / (l + r);

	IntersectResult[2] = tnsLinearInterpolate(lzl, lzr, ratio);
	tRdrGetInterpolatePoint2d(LL, LR, ratio, IntersectResult);

	return rev;
}
tnsRenderVert* tRdrFindSharedVertex(tnsRenderLine* rl, tnsRenderTriangle* rt) {
	if (rl->L == rt->V[0] || rl->L == rt->V[1] || rl->L == rt->V[2]) return rl->L;
	elif(rl->R == rt->V[0] || rl->R == rt->V[1] || rl->R == rt->V[2]) return rl->R;
	else return 0;
}
void tRdrFindEdgeNoVertex(tnsRenderTriangle* rt, tnsRenderVert* rv,tnsVector3d L, tnsVector3d R) {
	if (rt->V[0] == rv) {
		tMatVectorCopy3d(rt->V[1]->FrameBufferCoord, L);
		tMatVectorCopy3d(rt->V[2]->FrameBufferCoord, R);
	}elif(rt->V[1] == rv) {
		tMatVectorCopy3d(rt->V[2]->FrameBufferCoord, L);
		tMatVectorCopy3d(rt->V[0]->FrameBufferCoord, R);
	}elif(rt->V[2] == rv) {
		tMatVectorCopy3d(rt->V[0]->FrameBufferCoord, L);
		tMatVectorCopy3d(rt->V[1]->FrameBufferCoord, R);
	}
}
tnsRenderLine* tRdrAnotherEdge(tnsRenderTriangle* rt, tnsRenderVert* rv) {
	if (rt->V[0] == rv) {
		return rt->RL[1];
	}elif(rt->V[1] == rv) {
		return rt->RL[2];
	}elif(rt->V[2] == rv) {
		return rt->RL[0];
	}
}

int tRdrShareEdge(tnsRenderTriangle* rt, tnsRenderVert* rv, tnsRenderLine* rl) {
	tnsRenderVert* another = rv == rl->L ? rl->R : rl->L;

	if (rt->V[0] == rv) {
		if (another == rt->V[1] || another == rt->V[2]) return 1;
		return 0;
	}elif(rt->V[1] == rv) {
		if (another == rt->V[0] || another == rt->V[2]) return 1;
		return 0;
	}elif(rt->V[2] == rv) {
		if (another == rt->V[0] || another == rt->V[1]) return 1;
		return 0;
	}
}
int tRdrShareEdgeDirect(tnsRenderTriangle* rt, tnsRenderLine* rl) {
	if (rt->RL[0] == rl|| rt->RL[1]==rl || rt->RL[2] == rl) 
		return 1;
	return 0;
}
int tRdrTriangleLineImageSpaceIntersectTestOnly(tnsRenderTriangle* rt, tnsRenderLine* rl, double* From, double* To) {
	double dl, dr;
	double ratio;
	double is[3];
	int order[3];
	int LCross, RCross;
	int a, b, c;
	int ret;
	int noCross = 0;
	tnsVector3d TL,TR,LL,LR;
	tnsVector3d IntersectResult;
	tnsRenderVert* Share;
	int StL=0, StR=0;
	int OccludeSide;

	double
		*LFBC = rl->L->FrameBufferCoord,
		*RFBC = rl->R->FrameBufferCoord,
		*FBC0 = rt->V[0]->FrameBufferCoord,
		*FBC1 = rt->V[1]->FrameBufferCoord,
		*FBC2 = rt->V[2]->FrameBufferCoord;

	//bound box.
	if (TNS_MIN3(FBC0[2], FBC1[2], FBC2[2]) > TNS_MAX2(LFBC[2], RFBC[2])) return 0;
	if (TNS_MAX3(FBC0[0], FBC1[0], FBC2[0]) < TNS_MIN2(LFBC[0], RFBC[0])) return 0;
	if (TNS_MIN3(FBC0[0], FBC1[0], FBC2[0]) > TNS_MAX2(LFBC[0], RFBC[0])) return 0;
	if (TNS_MAX3(FBC0[1], FBC1[1], FBC2[1]) < TNS_MIN2(LFBC[1], RFBC[1])) return 0;
	if (TNS_MIN3(FBC0[1], FBC1[1], FBC2[1]) > TNS_MAX2(LFBC[1], RFBC[1])) return 0;
	
	if (Share = tRdrFindSharedVertex(rl, rt)) {
		tnsVector3d CL, CR;
		double r;
		int status;
		double r2;

		//if (rl->IgnoreConnectedFace/* && tRdrShareEdge(rt, Share, rl)*/) 
			//return 0;

		tRdrFindEdgeNoVertex(rt, Share, CL, CR);
		status = tRdrLineIntersectTest2d(LFBC, RFBC, CL, CR, &r);

		//LL[2] = 1 / tnsLinearItp(1 / LFBC[2], 1 / RFBC[2], r);
		LL[0] = tnsLinearItp(LFBC[0], RFBC[0], r);
		LL[1] = tnsLinearItp(LFBC[1], RFBC[1], r);
		LL[2] = tnsLinearItp(LFBC[2], RFBC[2], r);

		r2 = tRdrGetLinearRatio(CL, CR, LL);
		//LR[2] = 1 / tnsLinearItp(1 / CL[2], 1 / CR[2], r2);
		LR[0] = tnsLinearItp(CL[0], CR[0], r2);
		LR[1] = tnsLinearItp(CL[1], CR[1], r2);
		LR[2] = tnsLinearItp(CL[2], CR[2], r2);


		if (LL[2] <= (LR[2]+0.000000001)) return 0;

		StL = tRdrPointInsideTrianglef(LFBC, FBC0, FBC1, FBC2);
		StR = tRdrPointInsideTrianglef(RFBC, FBC0, FBC1, FBC2);

		if ((StL && Share == rl->R) ||
			(StR && Share == rl->L)) {
			*From = 0;
			*To = 1;
			return 1;
		}

		if (!status) return 0;

		if (rl->L == Share) {
			*From = 0;
			*To = r;
		}else {
			*From = r;
			*To = 1;
		}
		return 1;

	}

	//printf("L(%d %d %d) R(%d %d %d)\n     1(%d %d %d)\n    2(%d %d %d)\n    3(%d %d %d)\n",
	//	rl->L->V->P[0], rl->L->V->P[1], rl->L->V->P[2],
	//	rl->R->V->P[0], rl->R->V->P[1], rl->R->V->P[2],
	//	rt->V[0]->V->P[0], rt->V[0]->V->P[1], rt->V[0]->V->P[2],
	//	rt->V[1]->V->P[0], rt->V[1]->V->P[1], rt->V[1]->V->P[2],
	//	rt->V[2]->V->P[0], rt->V[2]->V->P[1], rt->V[2]->V->P[2]);


	a = tRdrLineIntersectTest2d(LFBC, RFBC, FBC0, FBC1, &is[0]);
	b = tRdrLineIntersectTest2d(LFBC, RFBC, FBC1, FBC2, &is[1]);
	c = tRdrLineIntersectTest2d(LFBC, RFBC, FBC2, FBC0, &is[2]);

	if (!a && !b &&!c) {
		if (!(StL = tRdrPointTriangleRelation(LFBC, FBC0, FBC1, FBC2)) &&
			!(StR = tRdrPointTriangleRelation(RFBC, FBC0, FBC1, FBC2))) {
			return 0;//not occluding
		}
	}

	StL = tRdrPointTriangleRelation(LFBC, FBC0, FBC1, FBC2);
	StR = tRdrPointTriangleRelation(RFBC, FBC0, FBC1, FBC2);

	INTERSECT_SORT_MIN_TO_MAX_3(is[0], is[1], is[2], order);

	if (StL) {
		INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
		INTERSECT_JUST_GREATER(is, order, -DBL_TRIANGLE_LIM, RCross);
		//if (is[LCross]>=0 || is[RCross] >= 1) return 0;
	}elif(StR) {
		INTERSECT_JUST_SMALLER(is, order, 1.0f+ DBL_TRIANGLE_LIM, LCross);
		INTERSECT_JUST_GREATER(is, order, 1.0f- DBL_TRIANGLE_LIM, RCross);
		//if (is[LCross] <= 0 || is[RCross] <= 1) return 0;
	}else {
		if (a) {
			if (b) {
				LCross = is[0] < is[1] ? 0 : 1;
				RCross = is[0] < is[1] ? 1 : 0;
			}else {
				LCross = is[0] < is[2] ? 0 : 2;
				RCross = is[0] < is[2] ? 2 : 0;
			}
		}elif(c) {
			LCross = is[1] < is[2] ? 1 : 2;
			RCross = is[1] < is[2] ? 2 : 1;
		}else {
			return 0;
		}
		//if (rl->IgnoreConnectedFace/* && tRdrShareEdge(rt, Share, rl)*/)
		//	return 0;
		if (TNS_MAX3(FBC0[2], FBC1[2], FBC2[2]) <( TNS_MIN2(LFBC[2], RFBC[2])- 0.000001)) {
			*From = is[LCross];
			*To = is[RCross];
			TNS_CLAMP((*From), 0, 1);
			TNS_CLAMP((*To), 0, 1);
			return 1;
		}
	}

	LL[2] = tRdrGetLineZ(LFBC, RFBC, is[LCross]);
	LR[2] = tRdrGetLineZ(LFBC, RFBC, is[RCross]);
	tRdrGetInterpolatePoint2d(LFBC, RFBC, is[LCross], LL);
	tRdrGetInterpolatePoint2d(LFBC, RFBC, is[RCross], LR);

	TL[2] = tRdrGetLineZPoint(rt->V[LCross]->FrameBufferCoord, rt->V[(LCross > 1 ? 0 : (LCross + 1))]->FrameBufferCoord, LL);
	TR[2] = tRdrGetLineZPoint(rt->V[RCross]->FrameBufferCoord, rt->V[(RCross > 1 ? 0 : (RCross + 1))]->FrameBufferCoord, LR);
	tMatVectorCopy2d(LL, TL);
	tMatVectorCopy2d(LR, TR);

	if (OccludeSide = tRdtGetZIntersectPoint(TL, TR, LL, LR, IntersectResult)) {
		real r = tRdrGetLinearRatio(LFBC, RFBC, IntersectResult);
		if (OccludeSide > 0) {
			if (r > 1 /*|| r < 0*/) return 0;
			*From = TNS_MAX2(r,0);
			*To = TNS_MIN2(is[RCross], 1);
		}else {
			if (r < 0 /*|| r > 1*/) return 0;
			*From = TNS_MAX2(is[LCross], 0);
			*To = TNS_MIN2(r, 1);
		}
		//*From = TNS_MAX2(TNS_MAX2(r, is[LCross]), 0);
		//*To = TNS_MIN2(r, TNS_MIN2(is[RCross], 1));
	}elif(IntersectResult[2]<0){
		if ((!StL) && (!StR) && (a + b + c < 2) || is[LCross]>is[RCross]) return 0;
		*From = is[LCross];
		*To = is[RCross];
	}else return 0;

	TNS_CLAMP((*From), 0, 1);
	TNS_CLAMP((*To), 0, 1);

	//if ((TNS_FLOAT_CLOSE_ENOUGH(*From, 0) && TNS_FLOAT_CLOSE_ENOUGH(*To, 1)) ||
	//	(TNS_FLOAT_CLOSE_ENOUGH(*To, 0) && TNS_FLOAT_CLOSE_ENOUGH(*From, 1)) ||
	//	TNS_FLOAT_CLOSE_ENOUGH(*From, *To)) return 0;

	return 1;
}
int tRdrTriangleLineImageSpaceIntersectTestOnlyV2(tnsRenderTriangle* rt, tnsRenderLine* rl, tnsCamera* cam, tnsMatrix44d vp, double* From, double* To) {
	double dl, dr;
	double ratio;
	double is[3] = {0};
	int order[3];
	int LCross=-1, RCross=-1;
	int a, b, c;
	int ret;
	int noCross = 0;
	tnsVector3d TL, TR, LL, LR;
	tnsVector3d IntersectResult;
	tnsRenderVert* Share;
	int StL = 0, StR = 0,In;
	int OccludeSide;

	tnsVector3d LV;
	tnsVector3d RV;
	real* CV = &cam->RenderViewDir;
	real DotL, DotR, DotLA, DotRA;
	real DotF;
	tnsRenderVert* Result, *rv;
	tnsVector3d GLocation,Trans;
	real Cut = -1;
	int NextCut, NextCut1;
	int status;


	double
		*LFBC = rl->L->FrameBufferCoord,
		*RFBC = rl->R->FrameBufferCoord,
		*FBC0 = rt->V[0]->FrameBufferCoord,
		*FBC1 = rt->V[1]->FrameBufferCoord,
		*FBC2 = rt->V[2]->FrameBufferCoord;

	//printf("%f %f %f   %f %f\n", FBC0[2], FBC1[2], FBC2[2], LFBC[2], RFBC[2]);

	//bound box.
	if (TNS_MIN3(FBC0[2], FBC1[2], FBC2[2]) > TNS_MAX2(LFBC[2], RFBC[2])) 
		return 0;
	if (TNS_MAX3(FBC0[0], FBC1[0], FBC2[0]) < TNS_MIN2(LFBC[0], RFBC[0])) return 0;
	if (TNS_MIN3(FBC0[0], FBC1[0], FBC2[0]) > TNS_MAX2(LFBC[0], RFBC[0])) return 0;
	if (TNS_MAX3(FBC0[1], FBC1[1], FBC2[1]) < TNS_MIN2(LFBC[1], RFBC[1])) return 0;
	if (TNS_MIN3(FBC0[1], FBC1[1], FBC2[1]) > TNS_MAX2(LFBC[1], RFBC[1])) return 0;

	if (tRdrShareEdgeDirect(rt, rl))
		return 0;

	a = tRdrLineIntersectTest2d(LFBC, RFBC, FBC0, FBC1, &is[0]);
	b = tRdrLineIntersectTest2d(LFBC, RFBC, FBC1, FBC2, &is[1]);
	c = tRdrLineIntersectTest2d(LFBC, RFBC, FBC2, FBC0, &is[2]);

	INTERSECT_SORT_MIN_TO_MAX_3(is[0], is[1], is[2], order);

	tMatVectorMinus3d(LV, rl->L->GLocation, rt->V[0]->GLocation);
	tMatVectorMinus3d(RV, rl->R->GLocation, rt->V[0]->GLocation);

	if(cam->CameraType == TNS_CAMERA_PERSPECTIVE) tMatVectorMinus3d(CV, cam->Base.GLocation, rt->V[0]->GLocation);

	DotL = tMatDot3d(LV, rt->GN, 0);
	DotR = tMatDot3d(RV, rt->GN, 0);
	DotF = tMatDot3d(CV, rt->GN, 0);

	if (!DotF) return 0;

	//if ((Share = tRdrFindSharedVertex(rl, rt)) ||
	//	(!rl->TL && ((Share = tRdrShareEdgeDirect(rt, rl->L->IntersectingLine)?rl->L:
	//	(tRdrShareEdgeDirect(rt, rl->R->IntersectingLine)? rl->R:0))))) {

	//	tnsVector3d CL, CR;
	//	double r;
	//	double r2;
	//	tnsRenderVert* another = Share == rl->L ? rl->R : rl->L;
	//	int in,index;
	//	real cut;

	//	if (tRdrShareEdgeDirect(rt, rl))
	//		return 0;

	//	in = tRdrPointInsideTrianglef(another->FrameBufferCoord, rt->V[0]->FrameBufferCoord, rt->V[1]->FrameBufferCoord, rt->V[2]->FrameBufferCoord);

	//	if (DotR*DotF + DotL*DotF >= 0) return 0;

	//	if (!rl->TL) {
	//		if (in) {
	//			*From = 0;
	//			*To = 1;
	//			return 1;
	//		}else {
	//			if (another->IntersectingLine == rt->RL[0]) cut = is[0];
	//			elif(another->IntersectingLine == rt->RL[1]) cut = is[1];
	//			elif(another->IntersectingLine == rt->RL[2]) cut = is[2];
	//			if (Share == rl->L) {
	//				*From = 0;
	//				INTERSECT_JUST_GREATER(is, order, 0, index);
	//				if (!TNS_ABC(index)) return 0;
	//				if (is[index] < DBL_EPSILON) return 0;
	//				*To = is[index];
	//				//TNS_CLAMP(*From, 0, 1);
	//				//TNS_CLAMP(*To, 0, 1);
	//				return 1;
	//			}else {
	//				INTERSECT_JUST_SMALLER(is, order, 1, index);
	//				if (!TNS_ABC(index)) return 0;
	//				if (is[index] < DBL_EPSILON) return 0;
	//				*From = is[index];
	//				*To = 1;
	//				//TNS_CLAMP(*From, 0, 1);
	//				//TNS_CLAMP(*To, 0, 1);
	//				return 1;
	//			}
	//		}
	//	}

	//	if (Share == rl->L) {
	//		if (in) {
	//			*From = 0;
	//			*To = 1;
	//			return 1;
	//		}else {
	//			cut = TNS_MAX3(is[0], is[1], is[2]);
	//			if (!TNS_MAX3_INDEX_ABC(is[0], is[1], is[2])) return 0;
	//			if (cut <= 1 && cut > 0.000001) {
	//				*From = 0;
	//				*To = cut;
	//				return 1;
	//			}
	//			return 0;
	//		}
	//	}else {
	//		if (in) {
	//			*From = 0;
	//			*To = 1;
	//			return 1;
	//		}
	//		else {
	//			cut = TNS_MIN3(is[0], is[1], is[2]);
	//			if (!TNS_MIN3_INDEX_ABC(is[0], is[1], is[2])) return 0;
	//			if (cut <= 0.99999 && cut > 0.000001) {
	//				*From = cut;
	//				*To = 1;
	//				return 1;
	//			}
	//			return 0;
	//		}
	//	}

	//}

	if (!a && !b && !c) {
		if (!(StL = tRdrPointTriangleRelation(LFBC, FBC0, FBC1, FBC2)) &&
			!(StR = tRdrPointTriangleRelation(RFBC, FBC0, FBC1, FBC2))) {
			return 0;//not occluding
		}
	}

	StL = tRdrPointTriangleRelation(LFBC, FBC0, FBC1, FBC2);
	StR = tRdrPointTriangleRelation(RFBC, FBC0, FBC1, FBC2);


	//for (rv = rt->IntersectingVerts.pFirst; rv; rv = rv->Item.pNext) {
	//	if (rv->IntersectWith == rt && rv->IntersectingLine == rl) {
	//		Cut = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], rv->FrameBufferCoord[0]);
	//		break;
	//	}
	//}


	DotLA = fabs(DotL); if (DotLA < DBL_EPSILON) { DotLA = 0; DotL = 0; }
	DotRA = fabs(DotR); if (DotRA < DBL_EPSILON) { DotRA = 0; DotR = 0; }
	if (DotL - DotR == 0) Cut = 100000;
	else if (DotL*DotR <= 0) {
		Cut = DotLA / fabs(DotL - DotR);
	}else {
		Cut = fabs(DotR + DotL) / fabs(DotL - DotR);
		Cut = DotRA > DotLA ? 1 - Cut : Cut;
	}

	if (cam->CameraType == TNS_CAMERA_PERSPECTIVE) {
		tnsLinearInterpolate3dv(rl->L->GLocation, rl->R->GLocation, Cut, GLocation);
		tMatApplyTransform44d(Trans, vp, GLocation);
		tMatVectorMultiSelf3d(Trans, (1 / Trans[3])/**HeightMultiply/2*/);
	}else {
		tnsLinearInterpolate3dv(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord, Cut, Trans);
		//tMatApplyTransform44d(Trans, vp, GLocation);
	}

	//Trans[2] = tMatDist3dv(GLocation, cam->Base.GLocation);
	//Trans[2] = cam->ZMin*cam->ZMax / (cam->ZMax - fabs(Trans[2]) * (cam->ZMax - cam->ZMin));


	Cut = tMatGetLinearRatio(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], Trans[0]);

	In = tRdrPointInsideTrianglef(Trans, rt->V[0]->FrameBufferCoord, rt->V[1]->FrameBufferCoord, rt->V[2]->FrameBufferCoord);


	if (StL == 2) {
		if (StR == 2) {
			INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
			INTERSECT_JUST_GREATER(is, order, 1 - DBL_TRIANGLE_LIM, RCross);
		}elif(StR == 1) {
			INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
			INTERSECT_JUST_GREATER(is, order, 1 - DBL_TRIANGLE_LIM, RCross);
		}elif(StR == 0) {
			INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
			INTERSECT_JUST_GREATER(is, order, 0, RCross);
		}
	}elif(StL == 1) {
		if (StR == 2) {
			INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
			INTERSECT_JUST_GREATER(is, order, 1 - DBL_TRIANGLE_LIM, RCross);
		}elif(StR == 1) {
			INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
			INTERSECT_JUST_GREATER(is, order, 1 - DBL_TRIANGLE_LIM, RCross);
		}elif(StR == 0) {
			INTERSECT_JUST_GREATER(is, order, DBL_TRIANGLE_LIM, RCross);
			if (TNS_ABC(RCross) && is[RCross] > (DBL_TRIANGLE_LIM)) {
				INTERSECT_JUST_SMALLER(is, order, DBL_TRIANGLE_LIM, LCross);
			}else {
				INTERSECT_JUST_SMALLER(is, order, -DBL_TRIANGLE_LIM, LCross);
				INTERSECT_JUST_GREATER(is, order, -DBL_TRIANGLE_LIM, RCross);
			}
		}
	}elif(StL == 0) {
		if (StR == 2) {
			INTERSECT_JUST_SMALLER(is, order, 1 - DBL_TRIANGLE_LIM, LCross);
			INTERSECT_JUST_GREATER(is, order, 1 - DBL_TRIANGLE_LIM, RCross);
		}elif(StR == 1) {
			INTERSECT_JUST_SMALLER(is, order, 1 - DBL_TRIANGLE_LIM, LCross);
			if (TNS_ABC(LCross) && is[LCross] < (1 - DBL_TRIANGLE_LIM)) {
				INTERSECT_JUST_GREATER(is, order, 1 - DBL_TRIANGLE_LIM, RCross);
			}else {
				INTERSECT_JUST_SMALLER(is, order, 1 + DBL_TRIANGLE_LIM, LCross);
				INTERSECT_JUST_GREATER(is, order, 1 + DBL_TRIANGLE_LIM, RCross);
			}
		}elif(StR == 0) {
			INTERSECT_JUST_GREATER(is, order, 0, LCross);
			if (TNS_ABC(LCross) && is[LCross] > 0) {
				INTERSECT_JUST_GREATER(is, order, is[LCross], RCross);
			}else {
				INTERSECT_JUST_GREATER(is, order, is[LCross], LCross);
				INTERSECT_JUST_GREATER(is, order, is[LCross], RCross);
			}
		}
	}


	real LF = DotL*DotF, RF = DotR*DotF;
	//int CrossCount = a + b + c;
	//if (CrossCount == 2) {
	//	INTERSECT_JUST_GREATER(is, order, 0, LCross);
	//	if (!TNS_ABC(LCross)) INTERSECT_JUST_GREATER(is, order, is[LCross], LCross);
	//	INTERSECT_JUST_GREATER(is, order, is[LCross], RCross);
	//}elif(CrossCount == 1 || StL+StR==1) {
	//	if (StL) {
	//		INTERSECT_JUST_GREATER(is, order, DBL_TRIANGLE_LIM, RCross);
	//		INTERSECT_JUST_SMALLER(is, order, is[RCross], LCross);
	//	}elif(StR) {
	//		INTERSECT_JUST_SMALLER(is, order, 1 - DBL_TRIANGLE_LIM, LCross);
	//		INTERSECT_JUST_GREATER(is, order, is[LCross], RCross);
	//	}
	//}elif(CrossCount == 0) {
	//	INTERSECT_JUST_SMALLER(is, order, 0, LCross);
	//	INTERSECT_JUST_GREATER(is, order, 1, RCross);
	//}

	if (LF <= 0 && RF <= 0 && (DotL || DotR)) {

		*From = TNS_MAX2(0,is[LCross]);
		*To = TNS_MIN2(1,is[RCross]);
		if (*From >= *To) return 0;
		return 1;
	}elif(LF >= 0 && RF <= 0 && (DotL || DotR)) {
		*From = TNS_MAX2(Cut, is[LCross]);
		*To = TNS_MIN2(1, is[RCross]);
		if (*From >= *To) return 0;
		return 1; 
	}elif(LF <= 0 && RF >= 0 && (DotL || DotR)) {
		*From = TNS_MAX2(0, is[LCross]);
		*To = TNS_MIN2(Cut, is[RCross]);
		if (*From >= *To) return 0;
		return 1;
	}else 
		return 0;
	return 1;
}

tnsRenderVert* tRdrTriangleShareEdge(tnsRenderTriangle* l, tnsRenderTriangle* r) {
	if (l->RL[0] == r->RL[0]) return r->RL[0];
	if (l->RL[0] == r->RL[1]) return r->RL[1];
	if (l->RL[0] == r->RL[2]) return r->RL[2];
	if (l->RL[1] == r->RL[0]) return r->RL[0];
	if (l->RL[1] == r->RL[1]) return r->RL[1];
	if (l->RL[1] == r->RL[2]) return r->RL[2];
	if (l->RL[2] == r->RL[0]) return r->RL[0];
	if (l->RL[2] == r->RL[1]) return r->RL[1];
	if (l->RL[2] == r->RL[2]) return r->RL[2];
	return 0;
}
tnsRenderVert* tRdrTriangleSharePoint(tnsRenderTriangle* l, tnsRenderTriangle* r) {
	if (l->V[0] == r->V[0]) return r->V[0];
	if (l->V[0] == r->V[1]) return r->V[1];
	if (l->V[0] == r->V[2]) return r->V[2];
	if (l->V[1] == r->V[0]) return r->V[0];
	if (l->V[1] == r->V[1]) return r->V[1];
	if (l->V[1] == r->V[2]) return r->V[2];
	if (l->V[2] == r->V[0]) return r->V[0];
	if (l->V[2] == r->V[1]) return r->V[1];
	if (l->V[2] == r->V[2]) return r->V[2];
	return 0;
}



// tRdrAssociateLineWithTile


tnsRenderVert* tRdrTriangleLineIntersectionTest(tnsRenderLine* rl, tnsRenderTriangle* rt, tnsRenderTriangle* Testing, tnsRenderVert* Last) {
	tnsVector3d LV;
	tnsVector3d RV;
	real DotL, DotR;
	tnsRenderVert* Result, *rv;
	tnsVector3d GLocation;
	tnsRenderVert* L = rl->L, *R = rl->R;
	int result;

	int i;

	for (rv = Testing->IntersectingVerts.pFirst; rv; rv = rv->Item.pNext) {
		if (rv->IntersectWith == rt && rv->IntersectingLine == rl) 
			return rv;
	}


	tMatVectorMinus3d(LV, L->GLocation, Testing->V[0]->GLocation);
	tMatVectorMinus3d(RV, R->GLocation, Testing->V[0]->GLocation);

	DotL = tMatDot3d(LV, Testing->GN, 0);
	DotR = tMatDot3d(RV, Testing->GN, 0);

	if (DotL*DotR > 0 || (!DotL&&!DotR))
		return 0;

	DotL = fabs(DotL);
	DotR = fabs(DotR);

	tnsLinearInterpolate3dv(L->GLocation, R->GLocation, DotL / (DotL + DotR), GLocation);

	if (Last && TNS_DOUBLE_CLOSE_ENOUGH(Last->GLocation[0], GLocation[0])
		&& TNS_DOUBLE_CLOSE_ENOUGH(Last->GLocation[1], GLocation[1])
		&& TNS_DOUBLE_CLOSE_ENOUGH(Last->GLocation[2], GLocation[2])) {

		Last->IntersectintLine2 = rl;
		return 0;
	}

	if (!(result = tRdrPointInsideTriangle3de(GLocation, Testing->V[0]->GLocation, Testing->V[1]->GLocation, Testing->V[2]->GLocation)))
		return 0;
	/*elif(result < 0) {
		return 0;
	}*/

	

	Result = memAquireOnly(sizeof(tnsRenderVert));

	if (DotL > 0 ||DotR<0) Result->Positive = 1; else Result->Positive = 0;

	//Result->IntersectingOnFace = Testing;
	Result->EdgeUsed = 1;
	//Result->IntersectL = L;
	Result->V = R; //Caution!
	//Result->IntersectWith = rt;
	tMatVectorCopy3d(GLocation, Result->GLocation);

	lstAppendItem(&Testing->IntersectingVerts, Result);

	return Result;
}
tnsRenderLine* tRdrTriangleGenerateIntersectionLineOnly(tnsRenderBuffer* rb, tnsRenderTriangle* rt, tnsRenderTriangle* Testing) {
	tnsRenderVert* L=0, *R=0;
	tnsRenderVert** Next = &L;
	tnsRenderLine* Result;
	tnsRenderVert* E0T=0;
	tnsRenderVert* E1T=0;
	tnsRenderVert* E2T=0;
	tnsRenderVert* TE0=0;
	tnsRenderVert* TE1=0;
	tnsRenderVert* TE2=0;
	tnsFrameBuffer* fb = rb->FrameBuffer;
	tnsVector3d* cl = rb->Scene->ActiveCamera->GLocation;
	real ZMax = ((tnsCamera*)rb->Scene->ActiveCamera)->ZMax;
	real ZMin = ((tnsCamera*)rb->Scene->ActiveCamera)->ZMin;
	tnsRenderVert* Share = tRdrTriangleSharePoint(Testing, rt);

	if (Share) {
		tnsRenderVert* NewShare;
		tnsRenderLine* rl = tRdrAnotherEdge(rt, Share);

		L = NewShare = memStaticAquire(&rb->RenderDataPool, (sizeof(tnsRenderVert)));

		NewShare->Positive = 1;
		NewShare->EdgeUsed = 1;
		//NewShare->IntersectL = L;
		NewShare->V = R; //Caution!
						 //Result->IntersectWith = rt;
		tMatVectorCopy3d(Share->GLocation, NewShare->GLocation);

		R = tRdrTriangleLineIntersectionTest(rl, rt, Testing, 0);

		if (!R) {
			rl = tRdrAnotherEdge(Testing, Share);
			R = tRdrTriangleLineIntersectionTest(rl, Testing, rt, 0);
			if (!R) return 0;
			lstAppendItem(&Testing->IntersectingVerts, NewShare);
		}else{
			lstAppendItem(&rt->IntersectingVerts, NewShare);
		}

	}else {
					  E0T = tRdrTriangleLineIntersectionTest(rt->RL[0], rt, Testing, 0); if (E0T && (!(*Next))) { (*Next) = E0T; (*Next)->IntersectingLine = rt->RL[0];  Next = &R; }
					  E1T = tRdrTriangleLineIntersectionTest(rt->RL[1], rt, Testing, L); if (E1T && (!(*Next))) { (*Next) = E1T; (*Next)->IntersectingLine = rt->RL[1];  Next = &R; }
		if (!(*Next)) E2T = tRdrTriangleLineIntersectionTest(rt->RL[2], rt, Testing, L); if (E2T && (!(*Next))) { (*Next) = E2T; (*Next)->IntersectingLine = rt->RL[2];  Next = &R; }

		if (!(*Next))TE0 = tRdrTriangleLineIntersectionTest(Testing->RL[0], Testing, rt,L); if (TE0 && (!(*Next))) { (*Next) = TE0; (*Next)->IntersectingLine = Testing->RL[0]; Next = &R; }
		if (!(*Next))TE1 = tRdrTriangleLineIntersectionTest(Testing->RL[1], Testing, rt, L); if (TE1 && (!(*Next))) { (*Next) = TE1; (*Next)->IntersectingLine = Testing->RL[1]; Next = &R; }
		if (!(*Next))TE2 = tRdrTriangleLineIntersectionTest(Testing->RL[2], Testing, rt, L); if (TE2 && (!(*Next))) { (*Next) = TE2; (*Next)->IntersectingLine = Testing->RL[2]; Next = &R; }

		if (!(*Next)) return 0;
	}
	tMatApplyTransform44d(L->FrameBufferCoord, fb->ViewProjection, L->GLocation);
	tMatApplyTransform44d(R->FrameBufferCoord, fb->ViewProjection, R->GLocation);
	tMatVectorMultiSelf3d(L->FrameBufferCoord, (1 / L->FrameBufferCoord[3])/**HeightMultiply/2*/);
	tMatVectorMultiSelf3d(R->FrameBufferCoord, (1 / R->FrameBufferCoord[3])/**HeightMultiply/2*/);

	//L->FrameBufferCoord[2] = tMatDist3dv(L->GLocation, cl);
	L->FrameBufferCoord[2] = ZMin*ZMax / (ZMax - fabs(L->FrameBufferCoord[2]) * (ZMax - ZMin));
	//if (L->FrameBufferCoord[3] < 0) {
	//	L->FrameBufferCoord[1] *= -1;
	//	L->FrameBufferCoord[0] *= -1;
	//	L->FrameBufferCoord[2] *= -1;
	//}

	//R->FrameBufferCoord[2] = tMatDist3dv(R->GLocation, cl);
	R->FrameBufferCoord[2] = ZMin*ZMax / (ZMax - fabs(R->FrameBufferCoord[2]) * (ZMax - ZMin));
	//if (R->FrameBufferCoord[3] < 0) {
	//	R->FrameBufferCoord[1] *= -1;
	//	R->FrameBufferCoord[0] *= -1;
	//	R->FrameBufferCoord[2] *= -1;
	//}



	/*if (L->FrameBufferCoord[3] < 0) {
		L->FrameBufferCoord[2] *= -1;
		L->FrameBufferCoord[1] *= -1;
		L->FrameBufferCoord[0] *= -1;
	}
	L->FrameBufferCoord[2] /= L->FrameBufferCoord[3];
	L->FrameBufferCoord[3] = L->FrameBufferCoord[2];
	L->FrameBufferCoord[2] = ZMin*ZMax / (ZMax - L->FrameBufferCoord[2] * (ZMax - ZMin));

	if (R->FrameBufferCoord[3] < 0) {
		R->FrameBufferCoord[2] *= -1;
		R->FrameBufferCoord[1] *= -1;
		R->FrameBufferCoord[0] *= -1;
	}
	R->FrameBufferCoord[2] /= R->FrameBufferCoord[3];
	R->FrameBufferCoord[3] = R->FrameBufferCoord[2];
	R->FrameBufferCoord[2] = ZMin*ZMax / (ZMax - R->FrameBufferCoord[2] * (ZMax - ZMin));
	*/

	L->IntersectWith = rt;
	R->IntersectWith = Testing;

	///*((1 / rl->L->FrameBufferCoord[3])*rb->FrameBuffer->H / 2)

	Result = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLine));
	Result->L = L;
	Result->R = R;
	tnsRenderLineSegment* rls = memStaticAquire(&rb->RenderDataPool, sizeof(tnsRenderLineSegment));
	lstAppendItem(&Result->Segments,rls);
	lstAppendItem(&rb->AllRenderLines, Result);
	lstAppendPointerStatic(&rb->IntersectionLines, &rb->RenderDataPool, Result);

	tnsglobal_TriangleIntersectionCount++;

	//rb->IntersectionCount++;

	return Result;
}
int tRdrTriangleCalculateIntersectionsInTile(tnsRenderBuffer* rb, tnsRenderTriangle* rt, tnsBoundingArea* ba) {
	tnsVector3d n, c = { 0 };
	tnsVector3d TL, TR;
	tnsRenderTriangle* TestingTriangle;
	tnsRenderLine* TestingLine;
	tnsRenderLine* Result = 0;
	tnsRenderVert* rv;
	nListItemPointer* lip, *NextLip;
	real l, r;
	int a = 0;

	real
		*FBC0 = rt->V[0]->FrameBufferCoord,
		*FBC1 = rt->V[1]->FrameBufferCoord,
		*FBC2 = rt->V[2]->FrameBufferCoord;

	for (lip = ba->AssociatedTriangles.pFirst; lip; lip = NextLip) {
		NextLip = lip->pNext;
		TestingTriangle = lip->p;
		if (TestingTriangle == rt || TestingTriangle->Testing == rt || tRdrTriangleShareEdge(rt,TestingTriangle))
			continue;	
		TestingTriangle->Testing = rt;
		real
			*RFBC0 = TestingTriangle->V[0]->FrameBufferCoord,
			*RFBC1 = TestingTriangle->V[1]->FrameBufferCoord,
			*RFBC2 = TestingTriangle->V[2]->FrameBufferCoord;

		if (TNS_MIN3(FBC0[2], FBC1[2], FBC2[2]) > TNS_MAX3(RFBC0[2], RFBC1[2], RFBC2[2])) continue;
		if (TNS_MAX3(FBC0[2], FBC1[2], FBC2[2]) < TNS_MIN3(RFBC0[2], RFBC1[2], RFBC2[2])) continue;
		if (TNS_MIN3(FBC0[0], FBC1[0], FBC2[0]) > TNS_MAX3(RFBC0[0], RFBC1[0], RFBC2[0])) continue;
		if (TNS_MAX3(FBC0[0], FBC1[0], FBC2[0]) < TNS_MIN3(RFBC0[0], RFBC1[0], RFBC2[0])) continue;
		if (TNS_MIN3(FBC0[1], FBC1[1], FBC2[1]) > TNS_MAX3(RFBC0[1], RFBC1[1], RFBC2[1])) continue;
		if (TNS_MAX3(FBC0[1], FBC1[1], FBC2[1]) < TNS_MIN3(RFBC0[1], RFBC1[1], RFBC2[1])) continue;


		Result = tRdrTriangleGenerateIntersectionLineOnly(rb, rt, TestingTriangle);
	}

}


int tRdrLineCrossesFrame(tnsVector2d L, tnsVector2d R) {
	real vx, vy;
	tnsVector4d Converted;
	real c1, c;

	if (-1 > TNS_MAX2(L[0], R[0])) return 0;
	if (1 < TNS_MIN2(L[0], R[0])) return 0;
	if (-1 > TNS_MAX2(L[1], R[1])) return 0;
	if (1 < TNS_MIN2(L[1], R[1])) return 0;

	vx = L[0] - R[0];
	vy = L[1] - R[1];

	c1 = vx * (-1 - L[1]) - vy * (-1 - L[0]);
	c = c1;

	c1 = vx * (-1 - L[1]) - vy * (1 - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	c1 = vx * (1 - L[1]) - vy * (-1 - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	c1 = vx * (1 - L[1]) - vy * (1 - L[0]);
	if (c1*c <= 0)return 1;
	else c = c1;

	return 0;
}

void tRdrComputeViewVector(tnsRenderBuffer* rb) {
	tnsVector3d Direction = { 0,0,-1 };
	tnsVector3d Trans;
	tnsMatrix44d inv;
	tMatInverse44d(inv, rb->Scene->ActiveCamera->GlobalTransform);
	tMatApplyRotation43d(Trans, inv, Direction);
	tMatVectorCopy3d(Trans, rb->ViewVector);
	tMatVectorMultiSelf3d(Trans, -1);
	tMatVectorCopy3d(Trans, ((tnsCamera*)rb->Scene->ActiveCamera)->RenderViewDir);
}

void tRdrComputeSceneContours(tnsRenderBuffer* rb) {
	real* ViewVector = &rb->ViewVector;
	tnsEdge* e;
	real Dot1 = 0, Dot2 = 0;
	real Result;
	tnsVector4d GNormal;
	int Add = 0;
	tnsCamera* c = rb->Scene->ActiveCamera;
	tnsRenderLine* rl;
	int ContourCount=0;
	int CreaseCount=0;
	int MaterialCount=0;

	rb->OverallProgress = 20;
	rb->CalculationStatus= TNS_CALCULATION_CONTOUR;
	nulThreadNotifyUsers("tns.render_buffer_list.calculation_status");

	if (c->CameraType == TNS_CAMERA_ORTHO) {
		tRdrComputeViewVector(rb);
	}

	for (rl = rb->AllRenderLines.pFirst; rl; rl = rl->Item.pNext) {
		//if(rl->Testing)
		//if (!tRdrLineCrossesFrame(rl->L->FrameBufferCoord, rl->R->FrameBufferCoord)) 
		//	continue;

		Add = 0; Dot1 = 0; Dot2 = 0;

		if (c->CameraType == TNS_CAMERA_PERSPECTIVE) {
			tMatVectorMinus3d(ViewVector, rl->L->GLocation, c->Base.GLocation);
		}

		if (rl->TL) Dot1 = tMatDot3d(ViewVector, rl->TL->GN, 0); else Add = 1;
		if (rl->TR) Dot2 = tMatDot3d(ViewVector, rl->TR->GN, 0); else Add = 1;

		if (!Add) {
			if ((Result = Dot1*Dot2) <= 0) Add = 1;
			elif(tMatDot3d(rl->TL->GN, rl->TR->GN, 0) < rb->CreaseCos) Add = 2;
			elif(rl->TL && rl->TR && rl->TL->F && rl->TR->F && rl->TL->F->MaterialID != rl->TR->F->MaterialID) Add = 3;
		}

		if (Add==1) {
			lstAppendPointerStatic(&rb->Contours, &rb->RenderDataPool, rl);
			ContourCount++;
		}elif(Add == 2) {
			lstAppendPointerStatic(&rb->CreaseLines, &rb->RenderDataPool, rl);
			CreaseCount++;
		}elif(Add == 3) {
			lstAppendPointerStatic(&rb->MaterialLines, &rb->RenderDataPool, rl);
			MaterialCount++;
		}
		if (ContourCount >= 100000) {
			tnsset_PlusRenderContourCount(rb, ContourCount);
			ContourCount = 0;
		}
		if (CreaseCount >= 100000) {
			tnsset_PlusRenderCreaseCount(rb, CreaseCount);
			CreaseCount = 0;
		}
		if (MaterialCount >= 100000) {
			tnsset_PlusRenderMaterialCount(rb, MaterialCount);
			MaterialCount = 0;
		}
	}
	tnsset_PlusRenderContourCount(rb, ContourCount);
	tnsset_PlusRenderCreaseCount(rb, CreaseCount);
	tnsset_PlusRenderMaterialCount(rb, MaterialCount);
}


void tRdrClearRenderState(tnsRenderBuffer* rb) {
	rb->ContourCount = 0;
	rb->ContourManaged = 0;
	rb->IntersectionCount = 0;
	rb->IntersectionManaged = 0;
	rb->MaterialLineCount = 0;
	rb->MaterialManaged = 0;
	rb->CreaseCount = 0;
	rb->CreaseManaged = 0;
	rb->CalculationStatus = TNS_CALCULATION_IDLE;
}

void tnsComputeFeatureLines(tnsRenderBuffer* rb) {
	int r, c;
	int ThreadCount = rb->ThreadCount; if (ThreadCount <= 0) ThreadCount = 1;
	int i;
	int res;
	tnsRenderTaskInfo* rti = CreateNewBuffer(tnsRenderTaskInfo,ThreadCount);

	rb->CreaseCos = cos(TNS_PI - rb->CreaseAngle);

	tRdrMakeInitialBoundingAreas(rb);
	tRdrComputeSceneContours(rb);
	tRdrAddTriangles(rb);

	rb->ContourManaged = rb->Contours.pFirst;
	rb->CreaseManaged = rb->CreaseLines.pFirst;
	rb->IntersectionManaged = rb->IntersectionLines.pFirst;
	rb->MaterialManaged = rb->MaterialLines.pFirst;

	tnsset_RenderOverallProgress(rb, 70);
	rb->CalculationStatus = TNS_CALCULATION_OCCLUTION;
	nulThreadNotifyUsers("tns.render_buffer_list.calculation_status");

	for (i = 0; i < ThreadCount; i++) {
		rti[i].ThreadID = i;
		rti[i].RenderBuffer = rb;
		thrd_create(&rti[i].ThreadHandle, THREAD_CalculateLineOcclusion, &rti[i]);
	}

	for (i = 0; i < ThreadCount; i++) {
		thrd_join(rti[i].ThreadHandle, &res);
	}

	FreeMem(rti);
}

void tRdrRebuildRenderDrawCommand(tnsRenderBuffer* rb, tnsRenderDrawCommand* rdc);

int tRdrDrawEdgePreview(tnsRenderBuffer* rb, tnsRenderDrawCommand* OverrideLayer, tnsGroup* OverrideGroup, real ThicknessScale, n2DViewUiExtra* e, tnsOffscreen* Off);
void tRdrSaveRenderBufferPreviewAsImage(tnsRenderBuffer* rb, char* name, tnsRenderDrawCommand* OverrideLayer, tnsGroup* OverrideGroup) {
	png_structp png;
	png_infop pngInfo;
	png_colorp pngColor;
	FILE* f;
	char* filename[1024] = {0};
	int W = rb->FrameBuffer->W, H = rb->FrameBuffer->H;
	int W2 = W / 2, H2 = H / 2;

	// file prepare =======================================

	sprintf(filename, "%s", name);

	if (strcmp(strgetLastSegmentSeperateBy(filename, '.'), "png") &&
		strcmp(strgetLastSegmentSeperateBy(filename, '.'), "PNG"))
		strcat(filename, ".png");

	f = fopen(filename, "wb");
	if (!f) return;

	png = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	pngInfo = png_create_info_struct(png);

	png_init_io(png, f);

	png_set_IHDR(png, pngInfo,
		W, H, 8,
		PNG_COLOR_TYPE_RGBA,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_BASE,
		PNG_FILTER_TYPE_BASE);

	png_color_8 pngBit;
	pngBit.alpha = 8;
	pngBit.red = 8;
	pngBit.green = 8;
	pngBit.blue = 8;
	png_set_sBIT(png, pngInfo, &pngBit);

	png_write_info(png, pngInfo);

	png_set_shift(png, &pngBit);
	png_set_packing(png);
	//png_set_swap_alpha(png);//argb->rgba//�������

	// render ============================================

	tnsOffscreen* off = tnsCreate2DOffscreenWithDepthSupersample(GL_RGBA, W, H, 0, 0);
	tnsOffscreen* noff = tnsCreate2DOffscreenBasic(GL_RGBA, W, H, 0);

	tnsDrawToOffscreen(off, 1, 0);
	tnsViewportWithScissor(0, 0, W, H);
	tnsResetViewMatrix();
	tnsResetModelMatrix();
	tnsResetProjectionMatrix();
	tnsOrtho(-W2, W2, -H2, H2, 100, -100);
	glClearColor(rb->BackgroundColor[0], rb->BackgroundColor[1], rb->BackgroundColor[2],
		rb->OutputTransparent ? 0 : rb->BackgroundColor[3]);
	glClear(GL_COLOR_BUFFER_BIT);
	tnsUseUiShader();
	tRdrDrawEdgePreview(rb, OverrideLayer, OverrideGroup, 1.0, 0, 0);

	tnsDrawToOffscreenOnlyBind(noff, 1, 0);
	tnsReadFromOffscreen(off);
	tnsPassColorBetweenOffscreens(off, noff, 0, 0, W + 0, H + 0, 0, 0, W + 0, H + 0, GL_NEAREST);

	tnsReadFromOffscreen(noff);

	u8bit* ImageData = CreateNewBuffer(u8bit, W * H * 4);

	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ROW_LENGTH, W);
	glPixelStorei(GL_PACK_IMAGE_HEIGHT, H);

	glReadPixels(0, 0, W, H, GL_RGBA, GL_UNSIGNED_BYTE, ImageData);
	
	u8bit** Rows = CreateNewBuffer(u8bit*, H);
	int i;
	for (i = 0; i < H; i++) {
		Rows[H-i-1] = &ImageData[i * W * 4];
	}	
	//png_bytep
	png_write_image(png, Rows);
	png_write_end(png, pngInfo);
	png_destroy_write_struct(png, pngInfo);
	fclose(f);

	tnsDelete2DOffscreen(off);
	tnsDelete2DOffscreen(noff);

	FreeMem(ImageData);
	FreeMem(Rows);
	//new offscreen buffer
	//tile render
	//save as png
}

int tRdrGetRenderTriangleSize(tnsRenderBuffer* rb) {
	return sizeof(tnsRenderTriangle) + (sizeof(tnsRenderLine*)*rb->ThreadCount);
}

void tRdrRebuildAllCommand(tnsRenderBuffer* rb);

int THREAD_ComputeFeatureLines(tnsRenderBuffer* rb) {
	tnsTransformGeomtryToRenderBuffer(rb);
	tRdrCullTriangles(rb);
	tRdrPerspectiveDivision(rb);
	tnsComputeFeatureLines(rb);

	tnsset_RenderOverallProgress(rb, 100);
	rb->CalculationStatus = TNS_CALCULATION_FINISHED;
	nulThreadNotifyUsers("tns.render_buffer_list.calculation_status");

	return 0;
}


int ACTINV_LoadExchange(nActuator* a, nEvent* e) {

	nulInvoke(a, "NUL_file_dialog", e, 0, 0, 0);

	
	nulNotifyUsers("tns.render_buffer_list.draw_commands");

	return NUL_RUNNING;
}
int ACTMOD_LoadExchange(nActuatorIntern* a, nEvent* e) {
	char buf[1024] = {0};
	char* format;
	
	if (a->ConfirmData) {
		if (a->ConfirmData->Mode == NUL_CONFIRM_CANCEL) {
			nulConfirmSameDataIfAny(a);
			return NUL_CANCELED;
		}
		if (a->ConfirmData->Mode == NUL_CONFIRM_OK) {
			nulGetConfirmString(a, buf); 
			if (buf[0]) {
				format = strgetLastSegmentSeperateBy(buf, '.');
				if (strcmp(format, "lasdexchange")) {
					nulEnableMessagePanel(a, 0, "Opps", "This is not a lasdexchange file", e->x - 180, e->y - 40, 250, e);
				}else {
					tnsCreateScene("Loaded Scene");
					tnsLoadExchange(buf);
					nulEnableMessagePanel(a, 0, "Oh Yes", "Loaded Selected Exchange File", e->x - 180, e->y - 40, 250, e);
				}
			}
		}
	}


	//nulRedrawCurrentPanel();
	nulNotifyUsers("tns");

	return NUL_FINISHED;
}
int ACTCHK_AreThereAnyRenderBuffers(nPropPack* This, nStringSplitor* ss) {
	if (T->RenderBuffers.pFirst) return 1;
	return 0;
}
int ACTINV_CreateNewRenderBuffer(nActuator* a, nEvent* e) {
	char Name[128] = "Buffer";
	tnsRenderBuffer* irb=T->RenderBuffers.pFirst;
	while (irb) {
		if (strIsTheSame(Name, irb->Name->Ptr)) {
			strMakeDifferentName(Name);
			irb = T->RenderBuffers.pFirst;
			continue;
		}
		irb = irb->Item.pNext;
	}

	tnsRenderBuffer* rb = tnsCreateRenderBuffer(Name);

	tRdrMakeFakeFrameBuffer( rb, 1920, 1080, 2);

	rb->CreaseAngle = rad(140);
	rb->ThreadCount = 1;
	rb->MaxOccludeLevel = 4;
	rb->FrameBuffer->SubPixelSample = 4;

	rb->BackgroundColor[0] = rb->BackgroundColor[1] = rb->BackgroundColor[2] = 0.2;
	rb->BackgroundColor[3] = 1;

	rb->ShowFast = rb->ShowLine = rb->ShowMaterial = 1;

	InitializeCriticalSection(&rb->csData);
	InitializeCriticalSection(&rb->csInfo);
	InitializeCriticalSection(&rb->csManagement);
	InitializeCriticalSection(&rb->RenderDataPool.csMem);

	nulNotifyUsers("tns.render_buffer_list");

	return NUL_FINISHED;
}
int ACTINV_DestroyBufferData(nActuator* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;
	tnsRenderElementLinkNode* reln;

	if (!rb) return NUL_FINISHED;

	tRdrClearRenderState(rb);

	lstEmptyDirect(&rb->Contours);
	lstEmptyDirect(&rb->IntersectionLines);
	lstEmptyDirect(&rb->CreaseLines);
	lstEmptyDirect(&rb->MaterialLines);
	lstEmptyDirect(&rb->AllRenderLines);

	tnsZeroGeomtryBuffers(rb->Scene);

	while (reln = lstPopItem(&rb->VertexBufferPointers)) {
		FreeMem(reln->Pointer);
	}

	while (reln = lstPopItem(&rb->TriangleBufferPointers)) {
		FreeMem(reln->Pointer);
	}

	memStaticDestroy(&rb->RenderDataPool);

	return NUL_FINISHED;
}
int ACTINV_DestroyRenderBuffer(nActuator* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;

	ACTINV_DestroyBufferData(a, e);

	memFree(rb->FrameBuffer);

	if (T->ActiveRenderBuffer == rb) {
		T->ActiveRenderBuffer = 0;
	}

	lstRemoveItem(&T->RenderBuffers, rb);
	strSafeDestroy(&rb->Name);

	memFree(rb);

	nulNotifyUsers("tns.render_buffer_list");

	return NUL_FINISHED;
}

int ACTCHK_IsBufferGeomrtryReady(nPropPack* This, nStringSplitor* ss) {
	return 1;
	//if (This && ((tnsRenderBuffer*)This->EndInstance)->State & TNS_RENDERBUFFER_GEOMETRY_COMPLETE) return 1;
	//return 0;
}
int ACTINV_ComputeFeatureLines(nActuator* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;
	thrd_t SeperateThread;

	if (!rb->Scene || !rb->Scene->ActiveCamera) {
		nulEnableMessagePanel(a, 0, "��ż", "��û��ѡ����Ҫִ�м���ĳ���", e->x-150, e->y-50, 200, e);
		return NUL_CANCELED;
	}

	rb->TriangleSize = tRdrGetRenderTriangleSize(rb);

	thrd_create(&SeperateThread, THREAD_ComputeFeatureLines, rb);

	

	if (!rb) return NUL_FINISHED;

	return NUL_FINISHED;
}
tnsRenderDrawCommand* tnsNewLineSetOnly(tnsRenderBuffer* rb) {
	tnsRenderDrawCommand* rdc, *irdc;
	char Name[128] = "LineLayer";

	if (!rb) return NUL_FINISHED;

	rdc = memAquire(sizeof(tnsRenderDrawCommand));

	rdc->Thickness = 1;
	rdc->Color[0] = rdc->Color[1] = rdc->Color[2] = rdc->Color[3] = 1;

	rdc->ParentRB = rb;

	rdc->Transparency = 1;

	rdc->DrawThisCommand = 1;
	//rdc->DepthTest = 1;
	rdc->DrawContour = rdc->DrawCrease = rdc->DrawMaterialLines = rdc->DrawIntersections = 1;

	for (irdc = rb->DrawCommands.pFirst; irdc;) {
		if (strIsTheSame(Name, irdc->Name->Ptr)) {
			strMakeDifferentName(Name);
			irdc = rb->DrawCommands.pFirst;
			continue;
		}
		irdc = irdc->Item.pNext;
	}
	strSafeSet(&rdc->Name, Name);

	return rdc;
}
int ACTINV_NewLineSet(nActuator* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;
	tnsRenderDrawCommand* rdc;

	rdc = tnsNewLineSetOnly(rb);

	if (strArgumentMatch(a->ExtraInstructionsP, "position", "tail")) {
		lstAppendItem(&rb->DrawCommands, rdc);
	}
	else {
		lstPushItem(&rb->DrawCommands, rdc);
	}

	tRdrRebuildRenderDrawCommand(rb, rdc);

	nulNotifyUsers("tns.render_buffer_list.draw_commands");

	return NUL_FINISHED;
}
int ACTINV_DestroyLineSet(nActuator* a, nEvent* e) {
	tnsRenderDrawCommand* rdc = a->This->EndInstance;

	if (!rdc) return NUL_FINISHED;

	lstRemoveItem(&rdc->ParentRB->DrawCommands, rdc);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDeleteBuffers(1, &rdc->VBO);

	memFree(rdc);

	nulNotifyUsers("tns.render_buffer_list.draw_commands");

	return NUL_FINISHED;
}
int ACTINV_MoveLineSet(nActuator* a, nEvent* e) {
	tnsRenderDrawCommand* rdc = a->This->EndInstance;

	if (!rdc) return NUL_FINISHED;

	if (strArgumentMatch(a->ExtraInstructionsP, "direction", "up")) {
		lstMoveUp(&rdc->ParentRB->DrawCommands, rdc);
	}elif(strArgumentMatch(a->ExtraInstructionsP, "direction", "down")) {
		lstMoveDown(&rdc->ParentRB->DrawCommands, rdc);
	}

	nulNotifyUsers("tns.render_buffer_list.draw_commands");
	nulNotifyUsers("tns.render_buffer_list");


	return NUL_FINISHED;
}
int ACTINV_RebuildAllCommands(nActuator* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;

	tRdrRebuildAllCommand(rb);
	return NUL_FINISHED;
}
int ACTINV_AutoCreateLineSet(nActuator* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;
	tnsRenderDrawCommand* rdc;

	rdc = tnsNewLineSetOnly(rb);
	rdc->Thickness = 2;

	lstAppendItem(&rb->DrawCommands, rdc);

	rdc = tnsNewLineSetOnly(rb);
	rdc->OccludeBegin = 1;
	rdc->OccludeEnd = 1;
	rdc->Color[0] = 0.314;
	rdc->Color[1] = 0.596;
	rdc->Color[2] = 1;

	lstAppendItem(&rb->DrawCommands, rdc);

	rdc = tnsNewLineSetOnly(rb);
	rdc->OccludeBegin = 2;
	rdc->OccludeEnd = 2;
	rdc->Color[0] = 0.135;
	rdc->Color[1] = 0.304;
	rdc->Color[2] = 0.508;	
	

	lstAppendItem(&rb->DrawCommands, rdc);

	tRdrRebuildAllCommand(rb);

	nulNotifyUsers("tns.render_buffer_list.draw_commands");

	return NUL_FINISHED;
}

static char Message[] = "Please fill in these fields before exporting image:";
static char MessageFolder[] = "    - Output folder";
static char MessagePrefix[] = "    - File name prefix";
static char MessageConnector[] = "    - File name connector";
static char MessageLayerName[] = "    - One or more layers have empty/illegal names.";
static char MessageSuccess[] = "Sucessfully Saved Image(s).";
static char MessageHalfSuccess[] = "Some image(s) failed to save.";
static char MessageFailed[] = "No saving action performed.";

int ACTINV_SaveRenderBufferPreview(nActuatorIntern* a, nEvent* e) {
	tnsRenderBuffer* rb = a->This->EndInstance;
	tnsRenderDrawCommand* rdc;
	char FullPath[1024] = "";

	if (!rb) return;

	tnsFrameBuffer *fb = rb->FrameBuffer;

	if (fb->OutputMode == TNS_OUTPUT_MODE_COMBINED) {
		if ((!fb->OutputFolder || !fb->OutputFolder->Ptr) || (!fb->ImagePrefix || !fb->ImagePrefix->Ptr)) {
			nPanelMessageList List = {0};
			nulAddPanelMessage(&List, Message);
			if ((!fb->OutputFolder || !fb->OutputFolder->Ptr)) nulAddPanelMessage(&List, MessageFolder);
			if ((!fb->ImagePrefix || !fb->ImagePrefix->Ptr)) nulAddPanelMessage(&List, MessagePrefix);
			nulAddPanelMessage(&List, MessageFailed);
			nulEnableMultiMessagePanel(a, 0, "Caution", &List, e->x, e->y, 500, e);
			return NUL_FINISHED;
		}
		strcat(FullPath, fb->OutputFolder->Ptr);
		strcat(FullPath, fb->ImagePrefix->Ptr);
		tRdrSaveRenderBufferPreviewAsImage(rb, FullPath, 0, 0);
	}elif(fb->OutputMode == TNS_OUTPUT_MODE_PER_LAYER) {
		nPanelMessageList List = { 0 };
		int unnamed = 0;
		if ((!fb->OutputFolder || !fb->OutputFolder->Ptr) || (!fb->ImagePrefix || !fb->ImagePrefix->Ptr)|| (!fb->ImageNameConnector || !fb->ImageNameConnector->Ptr)) {
			nulAddPanelMessage(&List, Message);
			if ((!fb->OutputFolder||!fb->OutputFolder->Ptr)) nulAddPanelMessage(&List, MessageFolder);
			if ((!fb->ImagePrefix|| !fb->ImagePrefix->Ptr)) nulAddPanelMessage(&List, MessagePrefix);
			if ((!fb->ImageNameConnector|| !fb->ImageNameConnector->Ptr)) nulAddPanelMessage(&List, MessageConnector);
			nulAddPanelMessage(&List, MessageFailed);
			nulEnableMultiMessagePanel(a, 0, "Caution", &List, e->x, e->y, 500, e);
			return NUL_FINISHED;
		}
		for (rdc = rb->DrawCommands.pFirst; rdc; rdc = rdc->Item.pNext) {
			FullPath[0] = 0;
			if ((!rdc->Name || !rdc->Name->Ptr) && !unnamed) {
				nulAddPanelMessage(&List, MessageHalfSuccess);
				nulAddPanelMessage(&List, MessageLayerName);
				unnamed = 1;
				continue;
			}
			strcat(FullPath, fb->OutputFolder->Ptr);
			strcat(FullPath, fb->ImagePrefix->Ptr);
			strcat(FullPath, fb->ImageNameConnector->Ptr);
			strcat(FullPath, rdc->Name->Ptr);
			tRdrSaveRenderBufferPreviewAsImage(rb, FullPath, rdc, 0);
		}
		if(unnamed)nulEnableMultiMessagePanel(a, 0, "Caution", &List, e->x, e->y, 500, e);
	}

	return NUL_FINISHED;
}
int ACTINV_SaveSingleLayer(nActuator* a, nEvent* e) {
	tnsRenderDrawCommand* rdc = a->This->EndInstance;
	char FullPath[1024] = "";
	int fail = 0;

	if (!rdc)return;

	tnsFrameBuffer* fb = rdc->ParentRB->FrameBuffer;

	if (!fb) return;

	nPanelMessageList List = { 0 };

	if ((!fb->OutputFolder || !fb->OutputFolder->Ptr) || (!fb->ImagePrefix || !fb->ImagePrefix->Ptr) || (!fb->ImageNameConnector || !fb->ImageNameConnector->Ptr)) {
		nulAddPanelMessage(&List, Message);
		if ((!fb->OutputFolder || !fb->OutputFolder->Ptr)) nulAddPanelMessage(&List, MessageFolder);
		if ((!fb->ImagePrefix || !fb->ImagePrefix->Ptr)) nulAddPanelMessage(&List, MessagePrefix);
		if ((!fb->ImageNameConnector || !fb->ImageNameConnector->Ptr)) nulAddPanelMessage(&List, MessageConnector);
		fail = 1;
	}
	if (!rdc->Name || !rdc->Name->Ptr) {
		nulAddPanelMessage(&List, MessageHalfSuccess);
		nulAddPanelMessage(&List, MessageLayerName);
		fail = 1;
	}
	if (fail) {
		nulAddPanelMessage(&List, MessageFailed);
		nulEnableMultiMessagePanel(a, 0, "Caution", &List, e->x, e->y, 500, e);
		return NUL_FINISHED;
	}


	FullPath[0] = 0;
	strcat(FullPath, fb->OutputFolder->Ptr);
	strcat(FullPath, fb->ImagePrefix->Ptr);
	strcat(FullPath, fb->ImageNameConnector->Ptr);
	strcat(FullPath, rdc->Name->Ptr);
	tRdrSaveRenderBufferPreviewAsImage(rdc->ParentRB, FullPath, rdc, 0);
	

	return NUL_FINISHED;
}



long tRdrCountLeveledEdgeSegmentCount(nListHandle* LineList, int OccludeLevel,tnsGroup* OverrideGroup,int Exclusive) {
	nListItemPointer* lip;
	tnsRenderLine* rl;
	tnsRenderLineSegment* rls;
	tns3DObject* o;
	int not = 0;
	long Count = 0;
	for (lip = LineList->pFirst; lip; lip = lip->pNext) {
		rl = lip->p;
		o = rl->ObjectRef;
		for (rls = rl->Segments.pFirst; rls; rls = rls->Item.pNext) {
			if (OverrideGroup) {
				if (tnsGroupHaveObject(OverrideGroup, rl->ObjectRef) && Exclusive) continue;
				if (!tnsGroupHaveObject(OverrideGroup, rl->ObjectRef) && !Exclusive) continue;
			}
			if (rls->OccludeLevel == OccludeLevel)Count++;
		}
	}
	return Count;
}
long tRdrCountIntersectionSegmentCount(tnsRenderBuffer* rb) {
	tnsRenderLine* rl;
	tnsRenderLineSegment* rls;
	long Count = 0;
	for (rl = rb->IntersectionLines.pFirst; rl; rl = rl->Item.pNext) {
		Count++;
	}
	return Count;
}
void* tRdrMakeLeveledEdgeVertexArray(tnsRenderBuffer* rb, nListHandle* LineList, float* VertexArray, int OccludeLevel, tnsGroup* OverrideGroup, int Exclusive) {
	nListItemPointer* lip;
	tnsRenderLine* rl;
	tnsRenderLineSegment* rls,*irls;
	tns3DObject* o;
	real W = rb->FrameBuffer->W/2, H = rb->FrameBuffer->H/2;
	long i = 0;
	float* V = VertexArray;
	for (lip = LineList->pFirst; lip; lip = lip->pNext) {
		rl = lip->p;
		o = rl->ObjectRef;
		if (OverrideGroup) {
			if (tnsGroupHaveObject(OverrideGroup, rl->ObjectRef) && Exclusive) continue;
			if (!tnsGroupHaveObject(OverrideGroup, rl->ObjectRef) && !Exclusive) continue;
		}

		if(o) o->LineRenderingDone = 1;
		for (rls = rl->Segments.pFirst; rls; rls = rls->Item.pNext) {
			if (rls->OccludeLevel == OccludeLevel) {
				*V = tnsLinearInterpolate(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], rls->at) * W;
				V++;
				*V = tnsLinearInterpolate(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1], rls->at) * H;
				V++;
				*V = tnsLinearInterpolate(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], rls->Item.pNext ? (irls = rls->Item.pNext)->at : 1) * W;
				V++;
				*V = tnsLinearInterpolate(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1], rls->Item.pNext ? (irls = rls->Item.pNext)->at : 1) * H;
				V++;
			}
		}
	}
	return V;
}
u32bit tRdrMakeBoundingAreaVBORecursive(float* V, u32bit Begin, tnsBoundingArea* ba,float HalfW,float HalfH) {
	u32bit Index = Begin;
	if (ba->Child) {
		Index = tRdrMakeBoundingAreaVBORecursive(V, Index, &ba->Child[0], HalfW, HalfH);
		Index = tRdrMakeBoundingAreaVBORecursive(V, Index, &ba->Child[1], HalfW, HalfH);
		Index = tRdrMakeBoundingAreaVBORecursive(V, Index, &ba->Child[2], HalfW, HalfH);
		Index = tRdrMakeBoundingAreaVBORecursive(V, Index, &ba->Child[3], HalfW, HalfH);
		return Index;
	}else {
		float* v = &V[Begin];
		v[0] = ba->L*HalfW; v[1] = ba->U*HalfH;
		v[2] = ba->L*HalfW; v[3] = ba->B*HalfH;

		v[4] = ba->L*HalfW; v[5] = ba->B*HalfH;
		v[6] = ba->R*HalfW; v[7] = ba->B*HalfH;

		v[8] = ba->R*HalfW; v[9] = ba->B*HalfH;
		v[10] = ba->R*HalfW; v[11] = ba->U*HalfH;

		v[12] = ba->R*HalfW; v[13] = ba->U*HalfH;
		v[14] = ba->L*HalfW; v[15] = ba->U*HalfH;
		return Index + 16;
	}
}
void tRdrMakeBoundingAreaVBOs(tnsRenderBuffer* rb) {
	float* V = CreateNewBuffer(float, rb->BoundingAreaCount * 16);
	int r, c;
	u32bit Index = 0;

	if (rb->BaVBO) glDeleteBuffers(1, &rb->BaVBO);

	for (r = 0; r < rb->FrameBuffer->TileCountY; r++) {
		for (c = 0; c < rb->FrameBuffer->TileCountX; c++) {
			Index = tRdrMakeBoundingAreaVBORecursive(V, Index, &rb->InitialBoundingAreas[r * 20 + c], rb->FrameBuffer->W/2.0,rb->FrameBuffer->H/2.0);
		}
	}

	glGenBuffers(1, &rb->BaVBO);
	glBindBuffer(GL_ARRAY_BUFFER, rb->BaVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*rb->BoundingAreaCount * 16, V, GL_DYNAMIC_DRAW);
	FreeMem(V);
}
int tRdrCountMaterialTriangles(nListHandle* ObjectList, tnsMaterial* m, tnsGroup* OverrideGroup, int IsExclude) {
	float* V, *N;
	tnsRenderElementLinkNode* reln;
	tnsMeshObject* o;
	tnsGroup* g;
	int Count = 0;

	if (!m) return 0;

	for (o = ObjectList->pFirst; o; o = o->Base.Item.pNext) {

		Count += tRdrCountMaterialTriangles(&o->Base.ChildObjects, m, OverrideGroup, IsExclude);

		if (o->Base.Type != TNS_OBJECT_MESH) continue;

		if (OverrideGroup) {
			int Have = tnsGroupHaveObject(OverrideGroup, o);
			if (Have && IsExclude) continue;
			if (!Have && !IsExclude) continue;
		}

		tnsFace* f;
		for (f = o->F.pFirst; f; f = f->Item.pNext) {
			if (f->MaterialID == m->ID) Count+=(f->TriangleCount);
		}
	}

	return Count;
}
int tRdrMakeMaterialPoints(tnsRenderBuffer* rb, float* V, float* N, int Offset, nListHandle* ObjectList, tnsMaterial* m, tnsGroup* OverrideGroup, int IsExclude) {
	tnsRenderElementLinkNode* reln;
	tnsMeshObject* o;
	int ofst = Offset;
	real W = rb->FrameBuffer->W / 2, H = rb->FrameBuffer->H / 2;

	for (o = ObjectList->pFirst; o; o = o->Base.Item.pNext) {
		if (o->Base.Type != TNS_OBJECT_MESH) continue;

		ofst = tRdrMakeMaterialPoints(rb, V, N, ofst, &o->Base.ChildObjects, m, OverrideGroup, IsExclude);
		
		if (OverrideGroup) {
			int Have = tnsGroupHaveObject(OverrideGroup, o);
			if (Have && IsExclude) continue;
			if (!Have && !IsExclude) continue;
		}


		o->Base.MaterialRenderingDone = 1;

		tnsFace* f;
		for (f = o->F.pFirst; f; f = f->Item.pNext) {
			if (f->MaterialID == m->ID) {
				int i = 0;
				tnsLoopItem* li = f->Loop.pFirst;
				if (!li) continue;
				tnsLoopItem* nli = li->Item.pNext;
				tnsLoopItem* nnli = nli->Item.pNext;
				for (i; i < f->TriangleCount; i++) {
					if (!li->Begin->RV) return 0;

					V[ofst + i * 9] = li->Begin->RV->FrameBufferCoord[0]*W;
					V[ofst + i * 9 + 1] = li->Begin->RV->FrameBufferCoord[1]*H;
					V[ofst + i * 9 + 2] = li->Begin->RV->FrameBufferCoord[2];

					V[ofst + i * 9 + 3] = nli->Begin->RV->FrameBufferCoord[0] * W;
					V[ofst + i * 9 + 4] = nli->Begin->RV->FrameBufferCoord[1] * H;
					V[ofst + i * 9 + 5] = nli->Begin->RV->FrameBufferCoord[2];

					V[ofst + i * 9 + 6] = nnli->Begin->RV->FrameBufferCoord[0] * W;
					V[ofst + i * 9 + 7] = nnli->Begin->RV->FrameBufferCoord[1] * H;
					V[ofst + i * 9 + 8] = nnli->Begin->RV->FrameBufferCoord[2];

					N[ofst + i * 9] =     f->GNormal[0];
					N[ofst + i * 9 + 1] = f->GNormal[1];
					N[ofst + i * 9 + 2] = f->GNormal[2];

					N[ofst + i * 9 + 3] = f->GNormal[0];
					N[ofst + i * 9 + 4] = f->GNormal[1];
					N[ofst + i * 9 + 5] = f->GNormal[2];

					N[ofst + i * 9 + 6] = f->GNormal[0];
					N[ofst + i * 9 + 7] = f->GNormal[1];
					N[ofst + i * 9 + 8] = f->GNormal[2];

					li = li->Item.pNext;
					if (!li) li = f->Loop.pFirst;
					nli = li->Item.pNext;
					if (!nli) nli = f->Loop.pFirst;
					nnli = nli->Item.pNext;
					if (!nnli) nnli = f->Loop.pFirst;

					li = li->Item.pNext;
					if (!li) li = f->Loop.pFirst;
					nli = li->Item.pNext;
					if (!nli) nli = f->Loop.pFirst;
					nnli = nli->Item.pNext;
					if (!nnli) nnli = f->Loop.pFirst;
				}
				ofst += (f->TriangleCount)*9;
			}
		}
	}
	return ofst;
}


void tRdrRebuildRenderDrawCommand(tnsRenderBuffer* rb, tnsRenderDrawCommand* rdc) {
	int Count=0;
	int level;
	float* V, *tv, *N;;

	if (!rb || !rb->Scene) return;

	if (rdc->VBO) {
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glDeleteBuffers(1, &rdc->VBO);
	}
	if (rdc->NBO) {
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glDeleteBuffers(1, &rdc->NBO);
	}

	if (rdc->Type == TNS_COMMAND_LINE) {
		glGenBuffers(1, &rdc->VBO);
		glBindBuffer(GL_ARRAY_BUFFER, rdc->VBO);

		for (level = rdc->OccludeBegin; level <= rdc->OccludeEnd; level++) {
			if (rdc->DrawContour) Count += tRdrCountLeveledEdgeSegmentCount(&rb->Contours, level, rdc->OverrideGroup, rdc->ExcludeGroup);
			if (rdc->DrawIntersections) Count += tRdrCountLeveledEdgeSegmentCount(&rb->IntersectionLines, level, rdc->OverrideGroup, rdc->ExcludeGroup);
			if (rdc->DrawCrease) Count += tRdrCountLeveledEdgeSegmentCount(&rb->CreaseLines, level, rdc->OverrideGroup, rdc->ExcludeGroup);
			if (rdc->DrawMaterialLines) Count += tRdrCountLeveledEdgeSegmentCount(&rb->MaterialLines, level, rdc->OverrideGroup, rdc->ExcludeGroup);
		}

		rdc->VertCount = Count * 2;

		tv = V = CreateNewBuffer(float, 4 * Count);

		for (level = rdc->OccludeBegin; level <= rdc->OccludeEnd; level++) {
			if (rdc->DrawContour)tv = tRdrMakeLeveledEdgeVertexArray(rb, &rb->Contours, tv, level, rdc->OverrideGroup, rdc->ExcludeGroup);
			if (rdc->DrawIntersections)tv = tRdrMakeLeveledEdgeVertexArray(rb, &rb->IntersectionLines, tv, level, rdc->OverrideGroup, rdc->ExcludeGroup);
			if (rdc->DrawCrease)tv = tRdrMakeLeveledEdgeVertexArray(rb, &rb->CreaseLines, tv, level, rdc->OverrideGroup, rdc->ExcludeGroup);
			if (rdc->DrawMaterialLines)tv = tRdrMakeLeveledEdgeVertexArray(rb, &rb->MaterialLines, tv, level, rdc->OverrideGroup, rdc->ExcludeGroup);
		}

		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4 * Count, V, GL_DYNAMIC_DRAW);

		FreeMem(V);
		return;
	}

	if (rdc->Type == TNS_COMMAND_MATERIAL || rdc->Type == TNS_COMMAND_EDGE) {
		if (!rdc->MaterialRef) {
			rdc->VertCount = 0;
			return;
		}

		Count = tRdrCountMaterialTriangles(&rb->Scene->Objects.pFirst, rdc->MaterialRef, rdc->OverrideGroup, rdc->ExcludeGroup);

		if (Count) {
			rdc->VertCount = Count;

			V = CreateNewBuffer(float, 9 * Count);
			N = CreateNewBuffer(float, 9 * Count);

			tRdrMakeMaterialPoints(rb, V, N, 0, &rb->Scene->Objects.pFirst, rdc->MaterialRef, rdc->OverrideGroup, rdc->ExcludeGroup);

			glGenBuffers(1, &rdc->VBO);
			glBindBuffer(GL_ARRAY_BUFFER, rdc->VBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 9 * Count, V, GL_DYNAMIC_DRAW);

			glGenBuffers(1, &rdc->NBO);
			glBindBuffer(GL_ARRAY_BUFFER, rdc->NBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 9 * Count, N, GL_DYNAMIC_DRAW);

			FreeMem(V);
			FreeMem(N);

			return;
		}
			
	}
	
}
void tRdrRebuildAllCommand(tnsRenderBuffer* rb) {
	tnsRenderDrawCommand* rdc;
	if (!rb) return;
	tnsCleanObjectFinishMarks(rb->Scene);
	for (rdc = rb->DrawCommands.pLast; rdc; rdc = rdc->Item.pPrev) {
		tRdrRebuildRenderDrawCommand(rb, rdc);
	}
	nulNotifyUsers("tns.render_buffer_list");
}


void DebugTileGridDraw(tnsFrameBuffer* fb);
void DebugTriangulateDraw(tnsRenderBuffer* rb);
void tRdrDrawBoundingAreas(tnsRenderBuffer* rb);

void tRdrDrawRenderBufferFrame(tnsRenderBuffer* rb, nBoxedTheme* bt) {
	tnsScene* s = rb->Scene;
	tnsMaterial* m;
	tnsBranchedCommand * bc;
	tnsShader* cs = T->uiShader;
	int W2 = rb->FrameBuffer->W / 2, H2 = rb->FrameBuffer->H / 2;

	tnsUseUiShader();
	tnsColor4dv(bt->Disabled->BorderColor->RGBA);
	tnsVertex2d(W2, H2);
	tnsVertex2d(-W2, H2);
	tnsVertex2d(-W2, -H2);
	tnsVertex2d(W2, -H2);
	tnsPackAs(GL_LINE_LOOP);
	tnsFlush();
}

void tnsUseTextureSobelColorMSShader();
int ae;
int tRdrDrawEdgePreview(tnsRenderBuffer* rb, tnsRenderDrawCommand* OverrideLayer, tnsGroup* OverrideGroup,
	real ThicknessScale, n2DViewUiExtra* e, tnsOffscreen* Off){

	long Count;
	int Level;
	tnsShader* cs = T->uiShader;
	tnsRenderDrawCommand* rdc;
	int ThicknessUntrue;

	if (rb->CalculationStatus != TNS_CALCULATION_FINISHED) return;

	tnsEnableShaderv(cs);

	//DebugTriangulateDraw(rb);
	//tRdrDrawBoundingAreas(rb);

	glClear(GL_DEPTH_BUFFER_BIT);

	glEnableVertexAttribArray(cs->vertexIndex);
	glDisableVertexAttribArray(cs->colorIndex);

	//glEnable(GL_DEPTH_TEST);

	if (rb->OverrideDisplay == TNS_OVERRIDE_DISPLAY_HIDE) return;

	for (rdc = OverrideLayer ? OverrideLayer : rb->DrawCommands.pLast; rdc; rdc = rdc->Item.pPrev) {
		int Draw = rb->OverrideDisplay == TNS_OVERRIDE_DISPLAY_SWAP ? !rdc->DrawThisCommand : rdc->DrawThisCommand;
		if (!Draw && rb->OverrideDisplay!=TNS_OVERRIDE_DISPLAY_SHOW) continue;
		if (!((rb->ShowLine && rdc->Type == TNS_COMMAND_LINE) ||
			  (rb->ShowFast && rdc->Type == TNS_COMMAND_EDGE) ||
			  (rb->ShowMaterial && rdc->Type == TNS_COMMAND_MATERIAL))) continue;

		if ((rdc->Type == TNS_COMMAND_MATERIAL || rdc->Type == TNS_COMMAND_EDGE) && (!rdc->MaterialRef || !rdc->VertCount)) continue;

		if (e && (rdc->TransparencyMode == TNS_TRANSPARENCY_DRAW_LAYERED || rdc->Type == TNS_COMMAND_EDGE)) {
			tnsDrawToAllExtraAttachments(e->OffScr); ae = glGetError();
			glClearColor(1, 1, 1, 0);
			glClear(GL_COLOR_BUFFER_BIT);
		}elif(Off && rdc->TransparencyMode == TNS_TRANSPARENCY_DRAW_LAYERED) {
			tnsDrawToAllExtraAttachments(Off);
			glClearColor(1, 1, 1, 0);
			glClear(GL_COLOR_BUFFER_BIT);
		}

		if (rdc->DepthTest || rdc->Type == TNS_COMMAND_EDGE) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);

		if (rdc->ClearDepthBuffer) glClear(GL_DEPTH_BUFFER_BIT);
		
		if (rdc->Type == TNS_COMMAND_LINE) {
			real LineW = rdc->Thickness * ThicknessScale;
			if (LineW >= 10) ThicknessUntrue = NUL_LINE_WIDTH_WARNING_TOO_WIDE;
			elif(LineW <= 0.1) ThicknessUntrue = NUL_LINE_WIDTH_WARNING_TOO_THIN;
			glBindBuffer(GL_ARRAY_BUFFER, rdc->VBO);
			glVertexAttribPointer(cs->vertexIndex, 2, GL_FLOAT, 0, 0, 0);
			glVertexAttrib4dv(cs->colorIndex, rdc->Color);
			glLineWidth(LineW);
			if (rdc->UseStipple) {
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(rdc->StippleSize, rdc->StipplePattern);
			}
			else glDisable(GL_LINE_STIPPLE);

			glDrawArrays(GL_LINES, 0, rdc->VertCount * 2);
		}

		if (rdc->Type == TNS_COMMAND_MATERIAL || rdc->Type == TNS_COMMAND_EDGE) {
			if (rdc->Type == TNS_COMMAND_EDGE){
				tnsEnableShaderv(T->ExtraBuffersShader);
				glBindBuffer(GL_ARRAY_BUFFER, rdc->VBO);
				glEnableVertexAttribArray(T->ExtraBuffersShader->vertexIndex);
				glVertexAttribPointer(T->ExtraBuffersShader->vertexIndex, 3, GL_FLOAT, 0, 0, 0);
				glBindBuffer(GL_ARRAY_BUFFER, rdc->NBO);
				glEnableVertexAttribArray(T->ExtraBuffersShader->normalIndex);
				glVertexAttribPointer(T->ExtraBuffersShader->normalIndex, 3, GL_FLOAT, 0, 0, 0);
				glDisableVertexAttribArray(T->ExtraBuffersShader->colorIndex);
				glVertexAttrib4dv(T->ExtraBuffersShader->colorIndex, rdc->OverrideColor ? rdc->Color : rdc->MaterialRef->Color);
			}
			else {
				tnsEnableShaderv(cs);
				glBindBuffer(GL_ARRAY_BUFFER, rdc->VBO);
				glVertexAttribPointer(cs->vertexIndex, 3, GL_FLOAT, 0, 0, 0);
				glDisableVertexAttribArray(cs->colorIndex);
				glVertexAttrib4dv(cs->colorIndex, rdc->OverrideColor ? rdc->Color : rdc->MaterialRef->Color);
			}
			glDrawArrays(GL_TRIANGLES, 0, rdc->VertCount * 3);
		}

		if (e && (rdc->TransparencyMode == TNS_TRANSPARENCY_DRAW_LAYERED || rdc->Type == TNS_COMMAND_EDGE)) {
			int W = e->OffScr->pColorTextures[0]->Width, H = e->OffScr->pColorTextures[0]->Height;
			tnsDrawToOffscreen(e->OffScr, 1, 0);
			tnsViewportWithScissor(0, 0, W, H);
			tnsResetProjectionMatrix();
			tnsOrtho(0,W,H,0,-100,100);
			glDisable(GL_DEPTH_TEST);
			
			if (rdc->Type == TNS_COMMAND_EDGE) {
				tnsColor4d(1, 1, 1, rdc->Transparency);
				tnsEnableShaderv(T->SobelColorShader);
				//tnsUseTextureSobelColorMSShader();
				glUniform1f(T->SobelColorShader->uniform0Index, (float)rdc->NormalEdgeClamp);
				glUniform1f(T->SobelColorShader->uniform1Index, (float)rdc->NormalEdgeStrength);
				glUniform1f(T->SobelColorShader->uniform2Index, (float)rdc->DepthEdgeClamp);
				glUniform1f(T->SobelColorShader->uniform3Index, (float)rdc->DepthEdgeStrength);
				//int a = glGetError();
				tnsUseTexture1(e->OffScr->pColorTextures[1]);
				tnsUseTexture2(e->OffScr->pColorTextures[2]);
				tnsDraw2DTextureDirectly(e->OffScr->pDepthTexture, 0, 0, W, H);
				tnsFlush();
				tnsUnbindTexture1(e->OffScr->pColorTextures[1]);
				tnsUnbindTexture2(e->OffScr->pColorTextures[2]);
				tnsUnbindTexture0(e->OffScr->pDepthTexture);
			}else {
				tnsUseTextureMultisampleAlphaShader();
				tnsColor4d(0, 0, 0, rdc->Transparency);
				tnsDraw2DTextureDirectly(e->OffScr->pColorTextures[1], 0, 0, W, H);
				tnsFlush();
				tnsUnbindTexture0(e->OffScr->pColorTextures[1]);
			}

			tnsEnableShaderv(cs);
			glDisableVertexAttribArray(cs->uvIndex);
			glDisableVertexAttribArray(cs->colorIndex);

			tnsResetViewMatrix();
			tnsResetProjectionMatrix();
			tnsOrtho(e->PanX - W*e->ZoomX / 2, e->PanX + W*e->ZoomX / 2, e->PanY - e->ZoomY*H / 2, e->PanY + e->ZoomY*H / 2, 100, -100);
		}elif(Off && rdc->TransparencyMode == TNS_TRANSPARENCY_DRAW_LAYERED) {
			int W = Off->pColorTextures[0]->Width, H = Off->pColorTextures[0]->Height;
			tnsDrawToOffscreen(Off, 1, 0);
			tnsViewportWithScissor(0, 0, W, H);
			tnsResetProjectionMatrix();
			tnsOrtho(0, W, H, 0, -100, 100);
			glDisable(GL_DEPTH_TEST);
			tnsUseTextureMultisampleAlphaShader();
			tnsColor4d(0, 0, 0, rdc->Transparency);
			tnsDraw2DTextureDirectly(e->OffScr->pColorTextures[1], 0, 0, W, H);
			tnsFlush();
			tnsEnableShaderv(cs);
			tnsUnbindTexture0(e->OffScr->pColorTextures[1]);
			glDisableVertexAttribArray(cs->uvIndex);
			glDisableVertexAttribArray(cs->colorIndex);

			tnsResetViewMatrix();
			tnsResetProjectionMatrix();
			tnsOrtho(-W / 2, W / 2, -H / 2, H / 2, 100, -100);
		}

		if (OverrideLayer) break;
	}
	glDisable(GL_LINE_STIPPLE);
	glDisable(GL_DEPTH_TEST);
	glLineWidth(1.0);

	return ThicknessUntrue;
}
void tRdrDrawBoundingAreas(tnsRenderBuffer* rb) {
	tnsShader* cs = T->uiShader;

	if (rb->InitialBoundingAreas && !rb->BaVBO) tRdrMakeBoundingAreaVBOs(rb);

	tnsEnableShaderv(cs);

	glEnableVertexAttribArray(cs->vertexIndex);
	glDisableVertexAttribArray(cs->colorIndex);

	glBindBuffer(GL_ARRAY_BUFFER, rb->BaVBO);
	glVertexAttribPointer(cs->vertexIndex, 2, GL_FLOAT, 0, 0, 0);
	glVertexAttrib4f(cs->colorIndex, 0, 1, 1, 0.1);
	glDrawArrays(GL_LINES, 0, rb->BoundingAreaCount*16);

}
void tRdrDrawThisBoundingArea(tnsBoundingArea* ba, real HalfW, real HalfH) {
	tnsVertex2d(ba->L*HalfW, ba->U*HalfH);
	tnsVertex2d(ba->R*HalfW, ba->U*HalfH);
	tnsVertex2d(ba->R*HalfW, ba->B*HalfH);
	tnsVertex2d(ba->L*HalfW, ba->B*HalfH);
	tnsPackAs(GL_TRIANGLE_FAN);
}
void tRdrDrawBoundingAreaTriangles(tnsBoundingArea* ba, real HalfW, real HalfH) {
	tnsRenderTriangle* rt;
	nListItemPointer* lip;
	for (lip = ba->AssociatedTriangles.pFirst; lip; lip = lip->pNext) {
		rt = lip->p;
		tnsVertex2d(rt->V[0]->FrameBufferCoord[0] * HalfW, rt->V[0]->FrameBufferCoord[1] * HalfH);
		tnsVertex2d(rt->V[1]->FrameBufferCoord[0] * HalfW, rt->V[1]->FrameBufferCoord[1] * HalfH);
		tnsVertex2d(rt->V[2]->FrameBufferCoord[0] * HalfW, rt->V[2]->FrameBufferCoord[1] * HalfH);
	}

	tnsPackAs(GL_TRIANGLES);
}
void tRdrDrawLinkedBoundingAreas(tnsRenderBuffer* rb,tnsBoundingArea* ba) {
	nListItemPointer* lip;
	tnsBoundingArea* nba;
	real HalfW = (real)rb->FrameBuffer->W / 2 , HalfH = (real)rb->FrameBuffer->H / 2;

	tnsUseUiShader();
	
	tnsColor4d(0.3, 0.7, 0.9, 0.7);
	tRdrDrawThisBoundingArea(ba, HalfW, HalfH);

	tnsColor4d(0.7, 0.9, 0.3, 0.5);
	for (lip = ba->LP.pFirst; lip; lip = lip->pNext) {
		tRdrDrawThisBoundingArea(lip->p, HalfW, HalfH);
	}
	for (lip = ba->RP.pFirst; lip; lip = lip->pNext) {
		tRdrDrawThisBoundingArea(lip->p, HalfW, HalfH);
	}
	for (lip = ba->UP.pFirst; lip; lip = lip->pNext) {
		tRdrDrawThisBoundingArea(lip->p, HalfW, HalfH);
	}
	for (lip = ba->BP.pFirst; lip; lip = lip->pNext) {
		tRdrDrawThisBoundingArea(lip->p, HalfW, HalfH);
	}

	tnsColor4d(0.9, 0.3, 0.3, 0.3);
	tRdrDrawBoundingAreaTriangles(ba, HalfW, HalfH);

	tnsFlush();
}

void tnsMakeRenderLinePreviewArray(nListHandle* lst,real* V) {
	tnsRenderLine* rl;
	u32bit i = 0;
	for (rl = lst->pFirst; rl; rl = rl->Item.pNext) {
		V[i] = rl->L->FrameBufferCoord[0];
		V[i + 1] = rl->L->FrameBufferCoord[1];
		V[i + 2] = rl->R->FrameBufferCoord[0];
		V[i + 3] = rl->R->FrameBufferCoord[1];
		i += 4;
	}
}
void tnsFrameBufferTilesPreviewDraw(nBoxedTheme* bt, tnsFrameBuffer* fb, n2DViewUiExtra* e) {
	tnsRenderTile* rt;
	int r, c;

	//tnsUseUiShader();

	for (r = 0; r < fb->TileCountY; r++) {
		for (c = 0; c < fb->TileCountX; c++) {
			rt = &TNS_TILE(fb->Tiles, r, c, fb->TileCountX);
			if (rt->AssociatedTriangles.pFirst) {
				tnsColor4d(1, 0.7, 0, 0.3);
				tnsVertex2d(rt->SubX - fb->W / 2, rt->SubY - fb->H / 2);
				tnsVertex2d(rt->SubX - fb->W / 2, rt->SubYLim - fb->H / 2);
				tnsVertex2d(rt->SubXLim - fb->W / 2, rt->SubYLim - fb->H / 2);
				tnsVertex2d(rt->SubXLim - fb->W / 2, rt->SubY - fb->H / 2);
				tnsPackAs(GL_TRIANGLE_FAN);
			}
			if (rt->AssociatedLines.pFirst) {
				tnsColor4d(0.2, 1, 1, 1);
				tnsVertex2d(rt->SubX - fb->W / 2, rt->SubY - fb->H / 2);
				tnsVertex2d(rt->SubX - fb->W / 2, rt->SubYLim - fb->H / 2);
				tnsVertex2d(rt->SubXLim - fb->W / 2, rt->SubYLim - fb->H / 2);
				tnsVertex2d(rt->SubXLim - fb->W / 2, rt->SubY - fb->H / 2);
				tnsPackAs(GL_LINE_LOOP);
			}
		}
	}
	tnsFlush();
}
void tnsRenderBufferPreviewDrawTransparentGrid(nBoxedTheme* bt, tnsRenderBuffer* rb, n2DViewUiExtra* e) {
	real V[8] = { 0 };
	real UV[8] = { 0 };
	nUiItem* ui = e->ParentUi;
	tnsFrameBuffer* fb = rb->FrameBuffer;

	real w = fb->W/2;
	real h = fb->H / 2;
	real gw = w / 10 / e->ZoomX;
	real gh = h / 10 / e->ZoomY;

	real px = e->PanX / 10 / e->ZoomX, py = - e->PanY / 10 / e->ZoomY;
	
	tnsUseTransparentGridShader();

	if (e->ImageDrawAlpha == 1)tnsColor4d(0.5, 0.5, 0.5, 0.5);
	else if (e->ImageDrawAlpha == 2)tnsColor4d(0.9, 0.9, 0.9, 0.8);
	else return;

	tnsMakeQuad2d(V,
		-w, h, w, h,
		w, -h, -w, -h);
	tnsMakeQuad2d(UV, -gw - px, -gh - py, gw - px,-gh - py,gw - px,gh - py,-gw - px,gh - py);
	tnsTexCoordArray2d(UV, 4);
	tnsVertexArray2d(V, 4);
	tnsPackAs(GL_TRIANGLE_FAN);
	tnsFlush();
}
void tnsRenderBufferPreviewDraw(nBoxedTheme* bt, tnsRenderBuffer* rb,n2DViewUiExtra* e) {
	tnsShader* cs = T->uiShader;
	tnsBranchedCommand* bc;
	tnsBoundingArea* ba;
	real* V;

	if (!rb) return;

	if (!e->ClearBackground) {
		tnsClearColorv(rb->BackgroundColor);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}else {
		;
	}
	tnsRenderBufferPreviewDrawTransparentGrid(bt, rb, e);

	tnsEnableShaderv(cs);

	if (rb->InitialBoundingAreas && (ba = tRdrGetBoundingArea(rb, e->ClickedX / rb->FrameBuffer->W * 2, e->ClickedY / rb->FrameBuffer->H * 2))) {
		tRdrDrawLinkedBoundingAreas(rb, ba);
	}

	//tRdrDrawMaterialPreview(rb);

	if (rb->CalculationStatus){
		e->LineWidthWarning = tRdrDrawEdgePreview(rb, 0, 0, e->AdaptiveLineWidth ? 1.0f / e->ZoomX : 1.0f, e, 0);
	}


	if (e->ImageDrawBorder) {
		tRdrDrawRenderBufferFrame(rb, bt);
	}

}


}

