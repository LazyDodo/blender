#include "NUL4.h"

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/
#include "NUL_Util.h"
#include "NUL_TNS.h"
#include <math.h>



extern tnsMain* T;
extern NUL MAIN;


tnsVert* tObjCreateVertex(tnsMeshObject* mo, real x, real y, real z) {
	tnsVert* v = memAquire(sizeof(tnsVert));
	v->P[0] = x;
	v->P[1] = y;
	v->P[2] = z;

	v->N[2] = 1.0f;

	//v->C[3] = 1.0f;

	lstAppendItem(&mo->V, v);
	mo->numV++;

	return v;
}

tnsVert* tObjEdgeNextVertex(tnsEdge* e, tnsVert* v);

tnsEdge* tObjEdgeExist(tnsMeshObject* mo, tnsVert* v0, tnsVert* v1) {
	tnsEdge* e;
	tnsEdgeItem* ei;

	//for (e = mo->E.pFirst; e; e = e->Item.pNext) {
	//	if (((e->VL == v0) && (e->VR == v1)) || ((e->VR == v0) && (e->VL == v1))) return e;
	//}

	for (ei = v0->EdgeItems.pFirst; ei; ei = ei->Item.pNext) {
		if (tObjEdgeNextVertex(ei->E, v0) == v1) return ei->E;
	}
	return 0;
}

tnsVert* tObjEdgeShareVertex(tnsEdge* e0, tnsEdge* e1) {
	if (e0->VL == e1->VL) return e0->VL;
	if (e0->VL == e1->VR) return e0->VL;
	if (e0->VR == e1->VL) return e0->VR;
	if (e0->VR == e1->VR) return e0->VR;
	return 0;
}

tnsVert* tObjEdgeNextVertex(tnsEdge* e, tnsVert* v) {
	if (v == e->VL) return e->VR;
	else if (v == e->VR) return e->VL;
	else return 0;
}

tnsEdge*  tObjCreateEdge(tnsMeshObject* mo, tnsVert* v0, tnsVert* v1) {
	tnsEdge* e;
	tnsEdgeItem* ei;
	if ((!v0) || (!v1)) return 0;
	if (!(e = tObjEdgeExist(mo, v0, v1))) {
		e = memAquire(sizeof(tnsEdge));
		e->VL = v0; 
		e->VR = v1;

		lstAppendItem(&mo->E, e);
		mo->numE++;

		lstAppendPointer(&v0->EdgeItems, e);
		lstAppendPointer(&v1->EdgeItems, e);
	}
	return e;
}

tnsFace* tObjCreateEmptyFace(tnsMeshObject* mo) {
	tnsFace* f;
	if (!mo) return;
	f = memAquire(sizeof(tnsFace));

	lstAppendItem(&mo->F, f);

	mo->numF++;

	return f;
}

int tObjAppendFaceEdge(tnsMeshObject* mo, tnsFace* f,tnsEdge* e) {
	tnsLoopItem* li;
	tnsVert* v;
	if (e->FL && e->FR) return 0;
	if (!f->Loop.pFirst) {
		lstAppendPointerSized(&f->Loop, e,sizeof(tnsLoopItem));
		if (!e->FL) e->FL = f;
		else e->FR = f;
		return 1;
	}
	else if ((v=tObjEdgeShareVertex(((tnsLoopItem*)f->Loop.pLast)->E, e))) {
		if (f->Loop.pFirst == f->Loop.pLast) {
			li = f->Loop.pFirst;
			li->Begin = tObjEdgeNextVertex(li->E, v);
			lstAppendPointer(&li->Begin->FaceItems, f);
		}
		li = lstAppendPointerSized(&f->Loop, e, sizeof(tnsLoopItem));
		li->Begin = v;
		lstAppendPointer(&li->Begin->FaceItems, f);
		if (!e->FL) e->FL = f;
		else e->FR = f;
		return 1;
	}
	return 0;
}

void tObjRecalculateVertNormal(tnsVert* v) {
	tnsFace* f; tnsLoopItem* li;
	tnsVert* nv;
	nListItemPointer* lip;
	tnsVector3d v1 = { 0 }, v2 = { 0 }, vn = {0};
	v->N[0] = v->N[1] = v->N[2] = 0;
	for (lip = v->FaceItems.pFirst; lip; lip = lip->pNext) {
		f = lip->p;
		f->FaceNormal[0] = f->FaceNormal[1] = f->FaceNormal[2] = 0;
		for (li = f->Loop.pFirst; li; li = li->Item.pNext) {
			if (tObjEdgeNextVertex(li->E, li->Begin) == v) {
				nv = (li->Item.pNext? 
					tObjEdgeNextVertex(((tnsLoopItem*)li->Item.pNext)->E,v):
					tObjEdgeNextVertex(((tnsLoopItem*)f->Loop.pFirst)->E,v));

				tMatVectorMinus3d(v1, li->Begin->P, v->P);
				tMatVectorMinus3d(v2, nv->P, v->P);
				tMatVectorCross3d(vn, v2, v1);
				tMatNormalizeSelf3d(vn);
				tMatVectorAccum3d(v->N, vn);
				tMatVectorAccum3d(f->FaceNormal,vn);

				break;
			}
		}
		tMatNormalizeSelf3d(f->FaceNormal);
	}
	tMatNormalizeSelf3d(v->N);
}
void tObjRecalculateFaceAverageNormal(tnsFace* f) {
	tnsLoopItem* li = f->Loop.pFirst;
	tnsVert* vb,*nv;
	nListItemPointer* lip;
	tnsVector3d v1 = { 0 }, v2 = { 0 }, vn = { 0 }, cg = { 0 };
	int i = 0;

	f->FaceNormal[0] = f->FaceNormal[1] = f->FaceNormal[2] = 0;
	f->Center[0] = f->Center[1] = f->Center[2] = 0;

	for (li; li; li = li->Item.pNext) {
		tMatVectorAccum3d(f->Center, li->Begin->P);
		i++;
	}
	tMatVectorMultiSelf3d(f->Center, 1.0f / (real)i);

	i -= 2;
	f->TriangleCount = i < 0 ? 0 : i;

	for (li = f->Loop.pFirst; li; li = li->Item.pNext) {

		nv = li->Item.pNext ? ((tnsLoopItem*)li->Item.pNext)->Begin : ((tnsLoopItem*)f->Loop.pFirst)->Begin;

		tMatVectorMinus3d(v1, li->Begin->P, f->Center);
		tMatVectorMinus3d(v2, nv->P, f->Center);
		tMatVectorCross3d(vn, v2, v1);
		tMatVectorAccum3d(f->FaceNormal, vn);
	}
	tMatNormalizeSelf3d(f->FaceNormal);
}

int tObjFinishFaceLoop(tnsMeshObject* mo, tnsFace* f) {
	if (!f->Loop.pFirst) {
		return 0;
	}
	if (tObjEdgeNextVertex(((tnsLoopItem*)f->Loop.pLast)->E, ((tnsLoopItem*)f->Loop.pLast)->Begin) == ((tnsLoopItem*)f->Loop.pFirst)->Begin) {
		tnsLoopItem* li = f->Loop.pFirst; if (!li->Item.pNext) {return 0;  }
		tObjRecalculateFaceAverageNormal(f);
		return 1;
	}
	return 0;
}
int tObjEmptyFaceLoop(tnsMeshObject* mo, tnsFace* f) {
	tnsLoopItem* li;
	while (li = lstPopItem(&f->Loop)) {
		memFree(li);
	}
	return 1;
}


tnsFace* tObjCreateQuadFace(tnsMeshObject* mo, tnsVert* v0, tnsVert* v1, tnsVert* v2, tnsVert* v3) {
	tnsEdge* e0,*e1,*e2,*e3;
	tnsFace* f;
	if (!mo) return 0;
	e0 = tObjCreateEdge(mo, v0, v1);
	e1 = tObjCreateEdge(mo, v1, v2);
	e2 = tObjCreateEdge(mo, v2, v3);
	e3 = tObjCreateEdge(mo, v3, v0);

	f = tObjCreateEmptyFace(mo);
	tObjAppendFaceEdge(mo, f, e0);
	tObjAppendFaceEdge(mo, f, e1);
	tObjAppendFaceEdge(mo, f, e2);
	tObjAppendFaceEdge(mo, f, e3);
	tObjFinishFaceLoop(mo,f);

	tObjRecalculateVertNormal(v0);
	tObjRecalculateVertNormal(v1);
	tObjRecalculateVertNormal(v2);
	tObjRecalculateVertNormal(v3);

	return f;
}

real tObjGetAngleAlreadyNormalized2v3d(tnsVector3d v1, tnsVector3d v2, tnsVector3d Positive) {
	tnsVector3d res;
	tMatVectorCross3d(res, v1, v2);
	real rad = TNS_PI - acos(tMatDot3d(v1, v2, 0));
	//rad = rad < 0 ? -rad : rad;

	if (Positive) {
		if (tMatDot3d(res, Positive, 0) < 0) {
			return rad;
		}else return 2 * TNS_PI - rad;
	}
	return rad;
}
real tObjGetLoopBeginAngle(tnsFace* f, tnsLoopItem* li, tnsVector3d Positive) {
	tnsVector3d v1, v2;
	tnsVector3d res;
	if (!li->Item.pPrev) {
		tMatVectorMinus3d(v1, li->Begin->P, ((tnsLoopItem*)f->Loop.pLast)->Begin->P);
		tMatVectorMinus3d(v2, ((tnsLoopItem*)li->Item.pNext)->Begin->P, li->Begin->P);
	}elif(!li->Item.pNext) {
		tMatVectorMinus3d(v1, li->Begin->P, ((tnsLoopItem*)li->Item.pPrev)->Begin->P);
		tMatVectorMinus3d(v2, ((tnsLoopItem*)f->Loop.pFirst)->Begin->P, li->Begin->P);
	}
	else {
		tMatVectorMinus3d(v1, li->Begin->P, ((tnsLoopItem*)li->Item.pPrev)->Begin->P);
		tMatVectorMinus3d(v2, ((tnsLoopItem*)li->Item.pNext)->Begin->P, li->Begin->P);
	}
	tMatNormalizeSelf3d(v1);
	tMatNormalizeSelf3d(v2);

	return tObjGetAngleAlreadyNormalized2v3d(v1, v2, Positive);
}
static tnsRenderLine *temp_LastTriangulateLine;
int tObjEarCutTrangulatePickWhom(tnsFace* f, tnsTriangulateNode* tn, int Index,int F,int FF,int L,int LL, int TotalVerts) {
	int i;
	tnsLoopItem
		*lif = tn[F].LoopItem,
		*lil = tn[L].LoopItem,
		*lic = tn[Index].LoopItem;
	int Fail = 0, Fail1 = 0, Fail2 = 0;

	for (i = 0; i < TotalVerts; i++) {
		if (tn[i].Picked) continue;
		tnsLoopItem* lit = tn[i].LoopItem;
		if (lit == lif || lit == lic || lit == lil) continue;
		if (tRdrPointInsideTriangle3d(lit->Begin->P, lic->Begin->P, lil->Begin->P, lif->Begin->P)) {
			Fail = 1; break;
		}
	}
	if (!Fail) return Index;

	lif = tn[FF].LoopItem;
	lil = tn[Index].LoopItem;
	lic = tn[F].LoopItem;
	for (i = 0; i < TotalVerts; i++) {
		if (tn[i].Picked) continue;
		tnsLoopItem* lit = tn[i].LoopItem;
		if (lit == lif || lit == lic || lit == lil) continue;
		if (tRdrPointInsideTriangle3d(lit->Begin->P, lic->Begin->P, lil->Begin->P, lif->Begin->P)) {
			Fail1 = 1; break;
		}
	}
	//if (!Fail1) return F;

	lif = tn[Index].LoopItem;
	lil = tn[LL].LoopItem;
	lic = tn[L].LoopItem;
	for (i = 0; i < TotalVerts; i++) {
		if (tn[i].Picked) continue;
		tnsLoopItem* lit = tn[i].LoopItem;
		if (lit == lif || lit == lic || lit == lil) continue;
		if (tRdrPointInsideTriangle3d(lit->Begin->P, lic->Begin->P, lil->Begin->P, lif->Begin->P)) {
			Fail2 = 1; break;
		}
	}
	if (Fail1 && Fail2) return Index;
	if (Fail1) return L;
	if (Fail2) return F;
	if (!Fail1 && !Fail2) return tn[L].Angle < tn[F].Angle ? L : F;

	//if (!Fail2) return L;

	return Index;
}
void tObjEarCutTrangulatePickNode(tnsFace* f, tnsTriangulateNode* tn, int Index,int TotalVerts,int Remaining, u32bit* Index3, tnsRenderLine** Line3, nListHandle* Lines, nMemoryPool* RenderMemory) {
	int FowardItem = Index-1, LatterItem = Index+1, FF, LL;
	tnsVector3d v1, v2;
	int Pick;

	for (FowardItem; FowardItem >= -1; FowardItem--) {
		if (FowardItem == -1) {FowardItem = TotalVerts; continue;}
		if (!tn[FowardItem].Picked) {
			break;
		}
		if (FowardItem == Index) break;
	}
	for (FF = FowardItem-1; FF >= -1; FF--) {
		if (FF == -1) { FF = TotalVerts; continue; }
		if (!tn[FF].Picked) {
			break;
		}
		if (FF == Index) break;
	}
	for (LatterItem; LatterItem <= TotalVerts; LatterItem++) {
		if (LatterItem == TotalVerts) {LatterItem = -1; continue;}
		if (!tn[LatterItem].Picked) {
			break;
		}
		if (LatterItem == Index) break;
	}
	for (LL = LatterItem+1; LL <= TotalVerts; LL++) {
		if (LL == TotalVerts) { LL = -1; continue; }
		if (!tn[LL].Picked) {
			break;
		}
		if (LL == Index) break;
	}

	if (FF == Index) 
		FF = LatterItem;
	if (LL == Index) 
		LL = FowardItem;

	//Should make a better solution
	Pick = tObjEarCutTrangulatePickWhom(f, tn, Index, FowardItem, FF, LatterItem, LL, TotalVerts);
	if (Pick != Index) return tObjEarCutTrangulatePickNode(f, tn, Pick, TotalVerts, Remaining, Index3, Line3, Lines, RenderMemory);

	Index3[0] = tn[Index].LoopItem->Begin->I;
	Index3[1] = tn[LatterItem].LoopItem->Begin->I;
	Index3[2] = tn[FowardItem].LoopItem->Begin->I;

	if (RenderMemory) {
		tnsRenderLine* rl;
		tnsRenderLineSegment* rls;
		if (Remaining > 3) {
			rl = memStaticAquire(RenderMemory, sizeof(tnsRenderLine));
			rls = memStaticAquire(RenderMemory, sizeof(tnsRenderLineSegment));
			temp_LastTriangulateLine = rl;
			lstAppendItem(&rl->Segments, rls);
			lstAppendItem(Lines, rl);
		}else {
			rl = temp_LastTriangulateLine;
		}
		if (Remaining > 3) {
			Line3[0] = tn[Index].FowardRL;
			Line3[1] = rl;
			Line3[2] = tn[Index].BackwardRL;
			tn[FowardItem].FowardRL = rl;
			tn[LatterItem].BackwardRL = rl;
			rl->R = tn[LatterItem].LoopItem->Begin->RV;
			rl->L = tn[FowardItem].LoopItem->Begin->RV;
		}else {
			Line3[0] = tn[Index].FowardRL;
			Line3[1] = tn[LatterItem].FowardRL;
			Line3[2] = tn[Index].BackwardRL;
		}
		//printf("L %d-%d %d-%d %d-%d\n", Line3[0]->L->V->I, Line3[0]->R->V->I, Line3[1]->L->V->I, Line3[1]->R->V->I, Line3[2]->L->V->I, Line3[2]->R->V->I);
	}

	tn[Index].Picked = 1;

	if (Remaining == 3) {
		tn[LatterItem].Picked = 1;
		tn[FowardItem].Picked = 1;
		return;
	}

	tMatVectorMinus3d(v1,  tn[FowardItem].LoopItem->Begin->P, tn[FF].LoopItem->Begin->P);
	tMatVectorMinus3d(v2, tn[LatterItem].LoopItem->Begin->P, tn[FowardItem].LoopItem->Begin->P);
	tMatNormalizeSelf3d(v1);
	tMatNormalizeSelf3d(v2);
	tn[FowardItem].Angle = tObjGetAngleAlreadyNormalized2v3d(v1, v2, f->FaceNormal);

	tMatVectorMinus3d(v1, tn[LatterItem].LoopItem->Begin->P, tn[FowardItem].LoopItem->Begin->P);
	tMatVectorMinus3d(v2, tn[LL].LoopItem->Begin->P, tn[LatterItem].LoopItem->Begin->P);
	tMatNormalizeSelf3d(v1);
	tMatNormalizeSelf3d(v2);
	tn[LatterItem].Angle = tObjGetAngleAlreadyNormalized2v3d(v1, v2, f->FaceNormal);
}
int tObjEarCutTrangulateRecursive(tnsFace* f, tnsTriangulateNode* tn, u32bit* Out_IndexList, tnsRenderLine** Out_LineList, nListHandle* Lines, int TotalVerts,int Remaining, nMemoryPool* RenderMemory) {
	int i = 0;
	int SmallestAngle = -100;
	u32bit Index3[3] = { 0 };
	tnsRenderLine* Line3[3] = { 0 };

	for (i; i < TotalVerts; i++) {
		if (tn[i].Picked) continue;
		//printf("tn[%d].Angle: %lf\n", tn[i].LoopItem->Begin->I, deg(tn[i].Angle));
		if (SmallestAngle < -10 || tn[SmallestAngle].Angle > tn[i].Angle) {
			SmallestAngle = i;
		}
	}
	//printf("\n");

	tObjEarCutTrangulatePickNode(f, tn, SmallestAngle, TotalVerts, Remaining, Index3, Line3, Lines, RenderMemory);

	if(Remaining>3)
		tObjEarCutTrangulateRecursive(f, tn, Out_IndexList, Out_LineList, Lines, TotalVerts, Remaining-1, RenderMemory);

	Out_IndexList[(Remaining - 3) * 3 ] = Index3[0];
	Out_IndexList[(Remaining - 3) * 3 + 1] = Index3[1];
	Out_IndexList[(Remaining - 3) * 3 + 2] = Index3[2];

	if (Out_LineList) {
		Out_LineList[(Remaining - 3) * 3] = Line3[0];
		Out_LineList[(Remaining - 3) * 3 + 1] = Line3[1];
		Out_LineList[(Remaining - 3) * 3 + 2] = Line3[2];
	}
}
int tObjEarCutTrangulate(tnsFace* f, tnsTriangulateNode* tn, u32bit* Out_IndexList, tnsRenderLine** Out_LineList, nListHandle* Lines, nMemoryPool* RenderMemory) {
	int i = 0;
	tnsLoopItem* li=f->Loop.pFirst;

	if (!tn || !Out_IndexList) return 0;

	for (i; i < f->TriangleCount+2; i++) {
		tn[i].Angle = tObjGetLoopBeginAngle(f, li, f->FaceNormal);
		tn[i].LoopItem = li;
		tn[i].Picked = 0;
		tn[i].FowardRL   = li->E->RenderLine;
		tn[i].BackwardRL = li->Item.pPrev ? ((tnsLoopItem*)li->Item.pPrev)->E->RenderLine : ((tnsLoopItem*)f->Loop.pLast)->E->RenderLine;
		li = li->Item.pNext;
	}

	tObjEarCutTrangulateRecursive(f, tn, Out_IndexList, Out_LineList, Lines, f->TriangleCount+2, f->TriangleCount+2, RenderMemory);
}
void tRdrCalculateRenderTriangleNormal(tnsRenderTriangle* rt);
void tObjAssirnRenderLineWithTriangle(tnsRenderTriangle* rt) {
	if (!rt->RL[0]->TL)
		rt->RL[0]->TL = rt;
	elif(!rt->RL[0]->TR)
		rt->RL[0]->TR = rt;

	if (!rt->RL[1]->TL)
		rt->RL[1]->TL = rt;
	elif(!rt->RL[1]->TR)
		rt->RL[1]->TR = rt;

	if (!rt->RL[2]->TL)
		rt->RL[2]->TL = rt;  
	elif(!rt->RL[2]->TR)
		rt->RL[2]->TR = rt;
}
void tObjSimpleTriangulateRender(tnsFace* F, tnsRenderTriangle* TBuf, int TriangleSize, tnsRenderVert* VBuf, nListHandle* LBuf, tnsRenderTriangle** Next, nMemoryPool* RenderData) {
	tnsLoopItem* fli = F->Loop.pFirst;	
	tnsLoopItem* lli = F->Loop.pLast;
	tnsVert* v0;
	tnsLoopItem* li;
	int i = 0;
	real* v1, *v2, *v3;
	u32bit* sb = T->SharedTEBuf;
	real d;
	tnsRenderTriangle* rt,*rt2;

	if (!fli) return;

	v0 = fli->Begin;

	if (F->TriangleCount > 2) {
		tObjEarCutTrangulate(F, T->SharedTN, T->SharedTEBuf, T->SharedRLBuf, LBuf, RenderData);

		for (i = 0; i < F->TriangleCount; i++) {

			rt = ((BYTE*)TBuf) + TriangleSize*i;

			rt->V[0] = &VBuf[T->SharedTEBuf[i * 3]];
			rt->V[1] = &VBuf[T->SharedTEBuf[i * 3 + 1]];
			rt->V[2] = &VBuf[T->SharedTEBuf[i * 3 + 2]];

			rt->RL[0] = T->SharedRLBuf[i * 3];
			rt->RL[1] = T->SharedRLBuf[i * 3 + 1];
			rt->RL[2] = T->SharedRLBuf[i * 3 + 2];

			rt->F = F;

			tRdrCalculateRenderTriangleNormal(rt);
			if (d = tMatDot3d(rt->GN, F->GNormal, 0) < 0)
				tMatVectorMultiSelf3d(rt->GN, -1);

			tMatVectorAccum3d(rt->GC, rt->V[0]->FrameBufferCoord);
			tMatVectorAccum3d(rt->GC, rt->V[1]->FrameBufferCoord);
			tMatVectorAccum3d(rt->GC, rt->V[2]->FrameBufferCoord);
			tMatVectorMultiSelf3d(rt->GC, 1.0f / 3.0f);

			//tObjGetNormal(TBuf[i].V[0]->V->P, TBuf[i].V[0]->V->P, TBuf[i].V[1]->V->P, TBuf[i].V[2]->V->P, TBuf[i].Normal);
		}
		for (i = 0; i < F->TriangleCount; i++) {
			tnsRenderTriangle* irt = ((BYTE*)TBuf) + TriangleSize*i;
			tObjAssirnRenderLineWithTriangle(irt);
		}
		(*Next) = ((BYTE*)TBuf) + TriangleSize*F->TriangleCount;
	}
	else {
		rt = ((BYTE*)TBuf) + TriangleSize*i;

		//if (F->TriangleCount == 1) {
		rt->V[0] = &VBuf[v0->I];
		rt->V[1] = &VBuf[((tnsLoopItem*)fli->Item.pNext)->Begin->I];
		rt->V[2] = &VBuf[((tnsLoopItem*)lli)->Begin->I];

		rt->RL[0] = fli->E->RenderLine;
		rt->RL[1] = ((tnsLoopItem*)fli->Item.pNext)->E->RenderLine;
		rt->RL[2] = lli->E->RenderLine;

		rt->F = F;

		tRdrCalculateRenderTriangleNormal(rt);
		if (d = tMatDot3d(rt->GN, F->GNormal, 0) < 0)
			tMatVectorMultiSelf3d(rt->GN, -1);


		tMatVectorAccum3d(rt->GC, rt->V[0]->FrameBufferCoord);
		tMatVectorAccum3d(rt->GC, rt->V[1]->FrameBufferCoord);
		tMatVectorAccum3d(rt->GC, rt->V[2]->FrameBufferCoord);
		tMatVectorMultiSelf3d(rt->GC, 1.0f / 3.0f);

		i++;

		rt2 = ((BYTE*)TBuf) + TriangleSize*i;

		if (F->TriangleCount == 2) {
			tnsRenderLine* NewRL = memStaticAquire(RenderData, sizeof(tnsRenderLine));
			lstAppendItem(LBuf, NewRL);
			tnsRenderLineSegment* rls = memStaticAquire(RenderData, sizeof(tnsRenderLineSegment));
			lstAppendItem(&NewRL->Segments, rls);

			rt2->V[0] = &VBuf[((tnsLoopItem*)fli->Item.pNext)->Begin->I];
			rt2->V[1] = &VBuf[((tnsLoopItem*)lli->Item.pPrev)->Begin->I];
			rt2->V[2] = &VBuf[((tnsLoopItem*)lli)->Begin->I];

			rt2->RL[0] = ((tnsLoopItem*)fli->Item.pNext)->E->RenderLine;
			rt2->RL[1] = ((tnsLoopItem*)lli->Item.pPrev)->E->RenderLine;
			rt2->RL[2] = NewRL;

			rt->RL[1] = NewRL;

			rt->F = F;
			rt2->F = F;

			NewRL->L = rt2->V[2];
			NewRL->R = rt2->V[0];

			tObjAssirnRenderLineWithTriangle(rt);
			tObjAssirnRenderLineWithTriangle(rt2);

			tRdrCalculateRenderTriangleNormal(rt2);
			if (d = tMatDot3d(rt2->GN, F->GNormal, 0) < 0)
				tMatVectorMultiSelf3d(rt2->GN, -1);


			tMatVectorAccum3d(rt2->GC, rt2->V[0]->FrameBufferCoord);
			tMatVectorAccum3d(rt2->GC, rt2->V[1]->FrameBufferCoord);
			tMatVectorAccum3d(rt2->GC, rt2->V[2]->FrameBufferCoord);
			tMatVectorMultiSelf3d(rt2->GC, 1.0f / 3.0f);

			i++;
		}else {
			tObjAssirnRenderLineWithTriangle(rt);
		}

		(*Next) = ((BYTE*)TBuf) + TriangleSize*i;
	}
}
void tObjSimpleTiangulateCommand(tnsFace* F, u32bit* EBuf, u32bit* I,int* CommandCount) {
	tnsLoopItem* fli = F->Loop.pFirst;
	tnsVert* v1;
	tnsLoopItem* li;
	if (!fli) return;
	v1= fli->Begin;
	if (!v1)return;

	if (F->TriangleCount > 2) {
		tObjEarCutTrangulate(F, T->SharedTN, &EBuf[(*I) * 3], 0, 0, 0);
		//tObjEarCutTrangulate(F, T->SharedTN, T->SharedTEBuf, 0, 0, 0);

		(*I) += F->TriangleCount;
		(*CommandCount) += F->TriangleCount;
	}
	else {
		for (li = fli->Item.pNext; li->Item.pNext; li = li->Item.pNext) {
			EBuf[(*I) * 3] = v1->I;
			EBuf[(*I) * 3 + 1] = li->Begin->I;
			EBuf[(*I) * 3 + 2] = ((tnsLoopItem*)li->Item.pNext)->Begin->I;
			(*I) += 1;
			(*CommandCount)++;
		}
	}
}

void tObjRefreshGeometeryIndex(tnsMeshObject* mo) {
	tnsVert* v;
	tnsEdge* e;
	tnsFace* f;
	long I = 0;
	float* VBuf;
	float* NBuf;
	GLuint* EBuf;
	int CommandCount = 0;

	if (mo->Batch) {
		tnsDeleteBatch(mo->Batch);
		mo->Batch = 0;
	}

	//if (mo->numV) {
		VBuf = CreateNewBuffer(float, mo->numV * 3);
		NBuf = CreateNewBuffer(float, mo->numV * 3);
	//}

	for (v = mo->V.pFirst; v; v = v->Item.pNext) {
		v->I = I;
		VBuf[I * 3] = v->P[0];
		VBuf[I * 3 + 1] = v->P[1];
		VBuf[I * 3 + 2] = v->P[2];
		NBuf[I * 3] = v->N[0];
		NBuf[I * 3 + 1] = v->N[1];
		NBuf[I * 3 + 2] = v->N[2];
		I++;
	}

	mo->Batch = tnsCreateBatch(mo->numV, 3, VBuf, NBuf);

	mo->TriangleCount = 0;
	mo->TriangulatedEdgeCount = 0;
	for (f = mo->F.pFirst; f; f = f->Item.pNext) {
		tObjRecalculateFaceAverageNormal(f);
		mo->TriangleCount += f->TriangleCount;
		mo->TriangulatedEdgeCount += f->TriangleCount - 1;
	}

	for (e = mo->E.pFirst; e; e = e->Item.pNext) {
		mo->TriangulatedEdgeCount++;
	}

	EBuf = CreateNewBuffer(GLuint, mo->TriangleCount * 3);

	I = 0;
	//int MaxI = mo->TriangleCount * 3;
	for (f = mo->F.pFirst; f; f = f->Item.pNext) {
		//if (I*3 >= MaxI) break;
		tObjSimpleTiangulateCommand(f, EBuf, &I,&CommandCount);
	}

	tnsCreateCommand(mo->Batch, CommandCount, 3, GL_TRIANGLES,TNS_INTERNAL_MODE_WIRE, EBuf);

	FreeMem(EBuf);
	FreeMem(VBuf);
	FreeMem(NBuf);
}

void tObjCreateMeshPlane(tns3DObject* o, real AtX, real AtY, real AtZ) {
	tnsMeshObject* mo = o;
	tnsVert* v0, *v1, *v2, *v3;

	v0 = tObjCreateVertex(o, AtX - 1, AtY - 1, AtZ);
	v1 = tObjCreateVertex(o, AtX - 1, AtY + 1, AtZ);
	v2 = tObjCreateVertex(o, AtX + 1, AtY + 1, AtZ);
	v3 = tObjCreateVertex(o, AtX + 1, AtY - 1, AtZ);

	tObjCreateQuadFace(o, v0, v1, v2, v3);
}

tns3DObject* tnsCreateMeshObjectPlane(tnsScene* Scene, char* Name, real AtX, real AtY, real AtZ) {
	tnsMeshObject* o;
	tnsWorld* w = &T->World;

	if (!Scene) return 0;

	o = memAquireHyper(sizeof(tnsMeshObject));
	tObjInitObjectBase(o, Scene, Name, TNS_OBJECT_MESH,
		AtX, AtY, AtZ, 0, 0, 0, 1.0f, TNS_ROTATION_XYZ_EULER, 1.0f);

	tObjCreateMeshPlane(o, AtX, AtY, AtZ);

	tObjRefreshGeometeryIndex(o);

	return o;
} 


tnsMeshObject* tnsLoadObjectFromFile(char* FileName) {
	FILE* file = fopen(FileName,"r");
	int NumVert, NumFace;
	tnsVert** V;
	int i,j;
	real Buf[3] = { 0 };
	real NBuf[3] = { 0 };
	int Index;
	int L, R, Beg;
	int NumEdges;
	tnsMeshObject* o;
	tnsEdge* e;
	tnsFace* f;
	tnsVert* v;

	o = memAquireHyper(sizeof(tnsMeshObject));
	tObjInitObjectBase(o, T->World.ActiveScene, FileName, TNS_OBJECT_MESH,
		0, 0, 0, 0, 0, 0, 1.0f, TNS_ROTATION_XYZ_EULER, 1.0f);

	fscanf(file, "%d", &NumVert);

	V = CreateNewBuffer(tnsVert*, NumVert);

	for (i = 0; i < NumVert; i++) {
		fscanf(file, "%d %lf %lf %lf %lf %lf %lf", &Index, &Buf[0], &Buf[1], &Buf[2], &NBuf[0], &NBuf[1], &NBuf[2]);
		V[i] = tObjCreateVertex(o, Buf[0], Buf[1], Buf[2]);
		V[i]->N[0] = NBuf[0];
		V[i]->N[1] = NBuf[1];
		V[i]->N[2] = NBuf[2];
	}

	fscanf(file, "%d", &NumFace);

	int Empty = 0;

	for (i = 0; i < NumFace; i++) {
		fscanf(file, "%d", &NumEdges);
		f = tObjCreateEmptyFace(o);
		Empty = 0;
		for (j = 0; j <= NumEdges; j++) {
			if(j != NumEdges) fscanf(file, "%d", &R);
			else L = Beg;
			if (j == 0) {
				Beg = L = R; continue;
			}else {
				e = tObjCreateEdge(o, V[L], V[R]);
				if (!tObjAppendFaceEdge(o, f, e)) { Empty = 1; }
				L = R;
				continue;
			}
		}
		if(!Empty)tObjFinishFaceLoop(o, f);
		else tObjEmptyFaceLoop(o, f);
	}

	//for (v = o->V.pFirst; v; v = v->Item.pNext) {
	//	tObjRecalculateVertNormal(v);
	//}

	tObjRefreshGeometeryIndex(o);

	FreeMem(V);

	return o;
}




int tExcReadKeyword(FILE* f) {
	char kwd[128] = { 0 };
	fscanf(f, "%s", kwd);

	while (kwd[0] == '#') {
		fscanf(f, "%[^\n]\n", kwd);
		fscanf(f, "%s", kwd);
	}

	if (!strcmp(kwd, "ObjectCount"))         return TNS_KEYWORD_OBJECT_COUNT;
	if (!strcmp(kwd, "ObjectName"))          return TNS_KEYWORD_OBJECT_NAME;
	if (!strcmp(kwd, "ObjectType"))          return TNS_KEYWORD_OBJECT_TYPE;
	if (!strcmp(kwd, "ObjectLocation"))      return TNS_KEYWORD_OBJECT_LOCATION;
	if (!strcmp(kwd, "ObjectRotation"))      return TNS_KEYWORD_OBJECT_ROTATION;
	if (!strcmp(kwd, "ObjectScale"))         return TNS_KEYWORD_OBJECT_SCALE;
	if (!strcmp(kwd, "ObjectChildrenCount")) return TNS_KEYWORD_OBJECT_CHILDREN_COUNT;
	if (!strcmp(kwd, "ObjectChildrenNames")) return TNS_KEYWORD_OBJECT_CHILDREN_NAMES;
	//if (!strcmp(kwd, "MeshCouont"))          return TNS_KEYWORD_MESH_COUNT;
	//if (!strcmp(kwd, "MeshName"))            return TNS_KEYWORD_MESH_NAME;
	if (!strcmp(kwd, "MeshVertexCount"))     return TNS_KEYWORD_MESH_VERTEX_COUNT;
	if (!strcmp(kwd, "MeshVertices"))        return TNS_KEYWORD_MESH_VERTICES;
	if (!strcmp(kwd, "MeshLoopCount"))       return TNS_KEYWORD_MESH_LOOP_COUNT;
	if (!strcmp(kwd, "MeshTopology"))        return TNS_KEYWORD_MESH_TOPOLOGY;
	if (!strcmp(kwd, "MeshUVCount"))         return TNS_KEYWORD_MESH_UV_COUNT;
	if (!strcmp(kwd, "MeshUVName"))          return TNS_KEYWORD_MESH_UV_NAME;
	if (!strcmp(kwd, "MaterialCount"))       return TNS_KEYWORD_MATERIAL_COUNT;
	if (!strcmp(kwd, "MaterialName"))        return TNS_KEYWORD_MATERIAL_NAME;
	if (!strcmp(kwd, "MaterialColor"))        return TNS_KEYWORD_MATERIAL_COLOR;
	if (!strcmp(kwd, "MaterialEnd"))        return TNS_KEYWORD_MATERIAL_END;

	if (!strcmp(kwd, "GroupCount"))        return TNS_KEYWORD_GROUP_COUNT;
	if (!strcmp(kwd, "GroupName"))        return TNS_KEYWORD_GROUP_NAME;
	if (!strcmp(kwd, "ObjectGroupCount"))        return TNS_KEYWORD_OBJECT_GROUP_COUNT;


	return -1;
}
tns3DObject* tExcGenerateTypedObject(FILE* file) {
	char Type[128] = { 0 };
	char kwd[32] = { 0 };

	int NumVert, NumFace;
	tnsVert** V;
	int i, j;
	real Buf[3] = { 0 };
	real NBuf[3] = { 0 };
	int Index;
	int L, R, Beg;
	int NumEdges;
	int NumUV;
	int UVItemCount;
	//tnsMeshObject* o;
	tnsEdge* e;
	tnsFace* f;
	tnsVert* v;
	int MatInfo;
	fscanf(file, "%s", Type);
	real a;
	int CamType = 0;

	tns3DObject* o;

	if (!strcmp(Type, "MESH")) {
		o = memAquireHyper(sizeof(tnsMeshObject));
		fscanf(file, "%s", kwd);//vcount
		fscanf(file, "%d", &NumVert);
		V = CreateNewBuffer(tnsVert*, NumVert);
		fscanf(file, "%s", kwd);//vlist
		for (i = 0; i < NumVert; i++) {
			fscanf(file, "%d %lf %lf %lf %lf %lf %lf", &Index, &Buf[0], &Buf[1], &Buf[2], &NBuf[0], &NBuf[1], &NBuf[2]);
			V[i] = tObjCreateVertex(o, Buf[0], Buf[1], Buf[2]);
			V[i]->N[0] = NBuf[0];
			V[i]->N[1] = NBuf[1];
			V[i]->N[2] = NBuf[2];
			V[i]->I = Index;
		}
		fscanf(file, "%s", kwd);//fcount
		fscanf(file, "%d", &NumFace);
		if (NumFace) {
			fscanf(file, "%s", kwd);//flist
			int Empty = 0;

			for (i = 0; i < NumFace; i++) {
				fscanf(file, "%d", &NumEdges);
				f = tObjCreateEmptyFace(o);
				Empty = 0;
				for (j = 0; j <= NumEdges; j++) {
					if (j != NumEdges) fscanf(file, "%d", &R);
					else L = Beg;
					if (j == 0) {
						Beg = L = R; continue;
					}
					else {
						e = tObjCreateEdge(o, V[L], V[R]);
						if (!tObjAppendFaceEdge(o, f, e)) { Empty = 1; }
						L = R;
						continue;
					}
				}
				if (!Empty)tObjFinishFaceLoop(o, f);
				else tObjEmptyFaceLoop(o, f);

				fscanf(file, "%d", &MatInfo);
				f->MaterialID = MatInfo;
			}
		}

		fscanf(file, "%s", kwd);//uvcount
		fscanf(file, "%d\n", &NumUV);
		for (i = 0; i < NumUV; i++) {
			fscanf(file, "%s", kwd);//uvname
			fscanf(file, "%[^\n]\n", kwd);//name
			fscanf(file, "%s", kwd);//tag
			fscanf(file, "%d", &UVItemCount);//count

			int j;
			for (j = 0; j < UVItemCount; j++) {
				fscanf(file, "%lf %lf\n", &a, &a);//count
			}
		}

		o->Type = TNS_OBJECT_MESH;
		tObjRefreshGeometeryIndex(o);
		FreeMem(V);

	}elif(!strcmp(Type, "EMPTY")) {
		o = memAquireHyper(sizeof(tns3DObject));
		o->Type = TNS_OBJECT_PLACEHOLDER;

	}elif(!strcmp(Type, "CAMERA")) {
		o = memAquireHyper(sizeof(tnsCamera));
		int active;

		fscanf(file, "%s", kwd);//fov
		fscanf(file, "%lf", &((tnsCamera*)o)->FOV);

		fscanf(file, "%s", kwd);//near
		fscanf(file, "%lf", &((tnsCamera*)o)->ZMin);

		fscanf(file, "%s", kwd);//far
		fscanf(file, "%lf", &((tnsCamera*)o)->ZMax);

		((tnsCamera*)o)->FocusDistance = 1;

		fscanf(file, "%s", kwd);//type
		if (strIsTheSame(kwd, "CameraType")){

			fscanf(file, "%s", kwd);//type
			if (strIsTheSame(kwd, "ORTHO")) {
				CamType = 1;
			}elif(strIsTheSame(kwd, "PERSP")) {
				CamType = 0;
			}

			((tnsCamera*)o)->CameraType = CamType;

			fscanf(file, "%s", kwd);//orth scale
			fscanf(file, "%lf", &((tnsCamera*)o)->OrthScale);
			fscanf(file, "%s", kwd);//active
		}
		
		fscanf(file, "%d", &active);
		if (active) {
			if(T->World.ActiveScene) T->World.ActiveScene->ActiveCamera=o;
		}

		o->Type = TNS_OBJECT_CAMERA;
	}else {
		o = memAquireHyper(sizeof(tns3DObject));
		o->Type = TNS_OBJECT_PLACEHOLDER;
	}

	return o;
}

int tExcLoadObject(FILE* f) {
	tnsVector3d Loc;
	tnsVector3d Rot;
	tnsVector3d Scale;
	tns3DObject* o;
	char Name[128] = { 0 };
	char NameC[128] = { 0 };
	int Type=-1;
	int ChildrenCount,i,GroupCount;
	int Finished = 0;
	nListHandle Children = {0};

	while (!Finished) {
		switch (tExcReadKeyword(f)) {
		case TNS_KEYWORD_OBJECT_NAME:
			fscanf(f, " %[^\n]\n", Name);
			break;
		case TNS_KEYWORD_OBJECT_LOCATION:
			fscanf(f, "%lf %lf %lf", &Loc[0], &Loc[1], &Loc[2]);
			break;
		case TNS_KEYWORD_OBJECT_ROTATION:
			fscanf(f, "%lf %lf %lf", &Rot[0], &Rot[1], &Rot[2]);
			break;
		case TNS_KEYWORD_OBJECT_SCALE:
			fscanf(f, "%lf %lf %lf", &Scale[0], &Scale[1], &Scale[2]);
			break;
		case TNS_KEYWORD_OBJECT_GROUP_COUNT:
			//wait for o
			fscanf(f, "%d", &GroupCount);
			break;
		case TNS_KEYWORD_OBJECT_CHILDREN_COUNT:
			fscanf(f, "%d", &ChildrenCount);
			fscanf(f, "%s\n", NameC);//tag
			for (i = 0; i < ChildrenCount; i++) {
				tnsChildObjectReadTemp* cort = memAquireOnly(sizeof(tnsChildObjectReadTemp));
				fscanf(f, "%[^\n]\n", &cort->Name);
				lstAppendItem(&Children, cort);
			}
			break;
		case TNS_KEYWORD_OBJECT_TYPE:
			o = tExcGenerateTypedObject(f);

			tObjInitObjectBase(o, T->World.ActiveScene, Name, o->Type,
				Loc[0], Loc[1], Loc[2], Rot[0], Rot[1], Rot[2], 1.0f, TNS_ROTATION_XYZ_EULER, Scale[0]);

			lstCopyHandle(&o->ChildReadTemp, &Children);

			//here
			for (i = 0; i < GroupCount; i++) {
				fscanf(f, "%[^\n]\n", NameC);
				tnsAddToGroup(o, tnsFindGroup(T->World.ActiveScene, NameC));
			}

			Finished = 1;
			break;
		}
	}
}

int tExcLoadMaterial(FILE* f) {
	tnsVector4d Color;
	char Name[128] = { 0 };
	int Type = -1;
	int Finished = 0;
	tnsMaterial* m;

	while (!Finished) {
		switch (tExcReadKeyword(f)) {
		case TNS_KEYWORD_MATERIAL_NAME:
			fscanf(f, " %[^\n]\n", Name);
			m = tnsCreateMaterial(T->World.ActiveScene, Name);

			break;
		case TNS_KEYWORD_MATERIAL_COLOR:
			fscanf(f, "%lf %lf %lf %lf", &m->Color[0], &m->Color[1], &m->Color[2], &m->Color[3]);
			tMatVectorCopy4d(m->Color, m->LinearColor);
			m->LinearColor[3] = 1;
			tnsLinearToLog(m->LinearColor, MAIN.Gamma);
			break;

		case TNS_KEYWORD_MATERIAL_END:
			Finished = 1;
			break;
		}
	}

}

int tExcPostReadRecursive(tns3DObject* po) {
	tnsScene* s = T->World.ActiveScene;
	tns3DObject* o, *NextObject, *co;
	tnsChildObjectReadTemp* cort;
	tns3DObject* From;

	for (o = po->ChildObjects.pFirst; o; o = NextObject) {
		NextObject = o->Item.pNext;

		while (cort = lstPopItem(&o->ChildReadTemp)) {

			co = tnsFindObject(s, cort->Name, &From);

			if (co) {
				if (co == NextObject && NextObject) NextObject = NextObject->Item.pNext;
				if (From) {
					lstRemoveItem(&From->ChildObjects, co);
				}
				else {
					lstRemoveItem(&s->Objects, co);
				}
				lstAppendItem(&o->ChildObjects, co);
				co->ParentObject = o;
			}
			memFree(cort);
		}

		tExcPostReadRecursive(o);
	}
}
int tExcPostRead() {
	tnsScene* s = T->World.ActiveScene;
	tns3DObject* o, *NextObject,*co;
	tnsChildObjectReadTemp* cort;
	tns3DObject* From;

	for (o = s->Objects.pFirst; o; o = NextObject) {
		NextObject = o->Item.pNext;

		while (cort = lstPopItem(&o->ChildReadTemp)) {

			co = tnsFindObject(s, cort->Name, &From);

			if (co) {
				if (co == NextObject && NextObject) NextObject = NextObject->Item.pNext;
				if(From){
					lstRemoveItem(&From->ChildObjects, co);
				}else {
					lstRemoveItem(&s->Objects, co);
				}
				lstAppendItem(&o->ChildObjects, co);
				co->ParentObject = o;
			}			
			memFree(cort);
		}
		tExcPostReadRecursive(o);
	}

	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		tObjApplyGlobalTransformMatrixRecursive(o);
	}

	tns_MakeAtlasTriggerBatch(s);
	tns_FeedAtlasData(s);
}

int tnsLoadExchange(char* FileName) {
	FILE* f = fopen(FileName, "r");

	int kw;

	int i,count;
	char Name[128] = { 0 };

	if (!f) return 0;

	while (kw = tExcReadKeyword(f)) {
		if (ftell(f) == EOF || kw < 0) break;

		switch (kw) {
		case TNS_KEYWORD_OBJECT_COUNT:
			fscanf(f, "%d", &count);
			for (i = 0; i < count; i++) {
				tExcLoadObject(f);
			}
			break;
		case TNS_KEYWORD_MATERIAL_COUNT:
			fscanf(f, "%d", &count);
			for (i = 0; i < count; i++) {
				tExcLoadMaterial(f);
			}
			break;
		case TNS_KEYWORD_GROUP_COUNT:
			fscanf(f, "%d\n", &count);
			for (i = 0; i < count; i++) {
				fscanf(f, "%[^\n]\n", Name);
				tnsCreateGroup(T->World.ActiveScene, Name);
			}
			break;
		}
	}	

	fclose(f);

	tExcPostRead();
	
	return 1;
}


void tnsDrawMeshObject(tnsMeshObject* mo, int OverrideDisplayMode) {
	int Mode = OverrideDisplayMode ? OverrideDisplayMode : mo->Base.DrawMode;
	tnsShader* cs = T->TEST_MatcapShader;
	tnsBranchedCommand* bc;

	if (!mo->Batch) return;
	

	//glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

	glBindBuffer(GL_ARRAY_BUFFER, mo->Batch->VBO);
	glEnableVertexAttribArray(cs->vertexIndex);
	glVertexAttribPointer(cs->vertexIndex, mo->Batch->Dimention, GL_FLOAT, 0, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, mo->Batch->NBO);
	glEnableVertexAttribArray(cs->normalIndex);
	glVertexAttribPointer(cs->normalIndex, mo->Batch->Dimention, GL_FLOAT, 0, 0, 0);

	for (bc = mo->Batch->Branches.pFirst; bc; bc = bc->Item.pNext) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bc->EBO);
		glDrawElements(bc->DrawAs, bc->ElementCount*bc->Dimention, GL_UNSIGNED_INT, 0);
	}

	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	
}