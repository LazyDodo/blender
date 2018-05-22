#include <gl/glew.h>
#include "NUL4.h"
#include <math.h>

extern tnsMain* T;
extern NUL MAIN;


//==============================================================[ ATLAS ]

void tns_InitAtlasInput(tnsScene* s, int width) {
	if (s->AtlasPointsL) return;

	s->AtlasPointsL = tnsCreate2DRectTextureAdvanced(GL_RGBA32F, GL_RGB, GL_FLOAT, width, width, 0);
	s->AtlasPointsR = tnsCreate2DRectTextureAdvanced(GL_RGBA32F, GL_RGB, GL_FLOAT, width, width, 0);
	s->AtlasFaceNormalL = tnsCreate2DRectTextureAdvanced(GL_RGBA32F, GL_RGB, GL_FLOAT, width, width, 0);
	s->AtlasFaceNormalR = tnsCreate2DRectTextureAdvanced(GL_RGBA32F, GL_RGB, GL_FLOAT, width, width, 0);
}
void tns_DestroyAtlasInputs(tnsScene* s) {
	tnsDeleteTexture(s->AtlasPointsL); s->AtlasPointsL = 0;
	tnsDeleteTexture(s->AtlasPointsR); s->AtlasPointsL = 0;
	tnsDeleteTexture(s->AtlasFaceNormalL); s->AtlasFaceNormalL = 0;
	tnsDeleteTexture(s->AtlasFaceNormalR); s->AtlasFaceNormalR = 0;
}
void tns_InitAtlasOutputFBO(n3DViewUiExtra* e) {
	tnsOffscreen* o;
	e->AtlasFBO = o = tnsCreateOffscreenHandle();
	tnsAttach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT0, e->AtlasPointsOutL);
	tnsAttach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT1, e->AtlasPointsOutR);
	tnsAttach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT2, e->AtlasEdgeLengthAccumDepth);
}
void tns_InitAtlasOutputs(n3DViewUiExtra* e, int width) {
	if (e->AtlasFBO) return;

	e->AtlasPointsOutL = tnsCreate2DRectTextureAdvanced(GL_RGBA32F, GL_RGBA, GL_FLOAT, width, width, 0);
	e->AtlasPointsOutR = tnsCreate2DRectTextureAdvanced(GL_RGBA32F, GL_RGBA, GL_FLOAT, width, width, 0);
	e->AtlasEdgeLengthAccumDepth = tnsCreate2DRectTextureAdvanced(GL_RGBA32F,GL_RGBA, GL_FLOAT, width, width, 0);

	tns_InitAtlasOutputFBO(e);
}
void tns_DestroyAtlasOutputs(n3DViewUiExtra* e) {
	tnsDelete2DOffscreen(e->AtlasFBO);
	e->AtlasPointsOutL = 0;
	e->AtlasPointsOutR = 0;
	e->AtlasEdgeLengthAccumDepth = 0;
	e->AtlasFBO = 0;
}

void tObjRefreshGeometeryIndex(tnsMeshObject* mo);

int tns_FeedAtlasDataRecursive(
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
void tns_FeedAtlasData(tnsScene* s) {
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

void tns_GetAtlasBufferPosition(int index, float* x, float* y) {
	int px, py;
	px = index % TNS_ATLAS_DEFAULT_INPUT_WIDTH;
	py = index / TNS_ATLAS_DEFAULT_INPUT_WIDTH;

	*x = (float)px / TNS_ATLAS_DEFAULT_INPUT_WIDTH * 2 - 1;
	*y = (float)py / TNS_ATLAS_DEFAULT_INPUT_WIDTH * 2 - 1;
}

int tns_MakeAtlasTriggerBatchRecursive(int Begin, tns3DObject* o) {
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

void tns_CalculateViewDir(tnsScene* s, tnsCamera* c, tnsVector3d Vec) {
	tnsVector3d Direction = { 0,0,-1 };
	tnsVector3d Trans;
	tnsMatrix44d inv;
	tMatInverse44d(inv, c->Base.GlobalTransform);
	tMatApplyRotation43d(Trans, inv, Direction);
	//tMatVectorCopy3d(Trans, Vec);
	//tMatVectorMultiSelf3d(Trans, -1);
	tMatVectorCopy3d(Trans, Vec);

}

void tns_TriggerThisObject(tns3DObject* o) {
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

void tns_TriggerAtlasTransformRecursive(tns3DObject* o) {

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
void tns_TriggerAtlasTransform(tnsScene* s, n3DViewUiExtra* e) {
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

void tns_MakeAtlasTriggerBatch(tnsScene* s) {
	tns3DObject* o;
	int Begin = 0;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		Begin = tns_MakeAtlasTriggerBatchRecursive(Begin, o);
	}
}

void tns_AtlasDrawTransform(n3DViewUiExtra* e) {
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

void tns_AtlasDrawEdgePreview(n3DViewUiExtra* e) {
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


tnsTextureSample* tns_AnyUncoveredSample(tnsTexture* t) {
	int row, col;
	u8bit* sample = 0;
	tnsTextureSample* ts = lstPopItem(&t->PendingSamples);

	return ts;
}

int tns_PointCloseToLine(real X, real Y, tnsVector2d From, tnsVector2d To, real Distance) {
	tnsVector3d P1 = { 0 }, P2 = { 0 }, P0 = { 0 };
	tnsVector3d Line12, Line01, Line02;
	tnsVector3d C1, C2;
	real Len1,Len2, Dot1, Dot2;

	P0[0] = X; P0[1] = Y;
	tMatVectorCopy2d(From, P1);
	tMatVectorCopy2d(To, P2);

	tMatVectorMinus3d(Line12, P2, P1);
	tMatVectorMinus3d(Line01, P1, P0);
	tMatVectorMinus3d(Line02, P2, P0);
	tMatNormalizeSelf3d(Line12);

	tMatVectorCross3d(C1, Line01, Line12);
	tMatVectorCross3d(C2, Line02, Line12);

	if ((Len1 = tMatLength3d(C1)) > Distance || (Len2 = tMatLength3d(C2)) > Distance) return 0;

	if (tMatDot2d(Line12, Line01, 0)>0 && tMatDist2dv(P0, P1)>Distance) return 0;
	if (tMatDot2d(Line12, Line02, 0)<0 && tMatDist2dv(P0, P2)>Distance) return 0;

	return 1;
}

void tns_RemoveSample(tnsTexture* t, int Row, int Col) {
	tnsTextureSample* ts;
	ts = t->SamplePtr[Row*t->Width + Col];
	t->SamplePtr[Row*t->Width + Col] = 0;
	lstRemoveItem(&t->PendingSamples, ts);
	lstAppendItem(&t->ErasedSamples, ts);
}


//Directions
//(0,0)
//    1 2 3
//    8 0 4
//    7 6 5

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


int tns_ReverseDirection(int From) {
	From -= 4;
	if (From <= 0)From += 8;
	return From;
}

int tns_DirectionDeviate(int From, int To) {
	return _TNS_Deviates[From - 1][To - 1];
}


int tns_DetectDirection(tnsTexture* t, int Col, int Row, int LastDirection) {
	int Deviate[9] = {100};
	int MinDeviate = 0;
	int i;
	tnsTextureSample* ts;

	for (i = 0; i < 8; i++) {
		TNS_CLAMP_TEXTURE_CONTINUE(t, (_TNS_ColOffsets[i] + Col), (_TNS_RowOffsets[i] + Row));
		if (ts = t->SamplePtr[(_TNS_ColOffsets[i] + Col) + (_TNS_RowOffsets[i] + Row) * t->Width]) {
			if (!LastDirection) return i + 1;
			Deviate[i+1] = tns_DirectionDeviate(i, LastDirection);
			if (!MinDeviate || Deviate[MinDeviate] > Deviate[i + 1]) MinDeviate = i + 1;
		}
	}
	//if (MinDeviate) {
	//	tns_RemoveSample(ts, _TNS_ColOffsets[Deviate[MinDeviate]], _TNS_RowOffsets[Deviate[MinDeviate]]);

	//}

	return MinDeviate;
}

int tns_GrowSnakeR(tnsTexture* t, tnsLineStrip* ls, tnsLineStripPoint* ThisP, int Direction) {
	tnsLineStripPoint* NewP = ThisP,*p2;
	int Length = 5;
	int l = 0;
	int Deviate,Dir=Direction,NewDir;
	int AddPoint;
	int TX = NewP->P[0], TY = NewP->P[1];

	while (NewDir = tns_DetectDirection(t, TX, TY, Dir)) {
		AddPoint = 0;
		Deviate = tns_DirectionDeviate(NewDir, Dir);
		Dir = NewDir;

		l++;
		TX += _TNS_ColOffsets[NewDir-1];
		TY += _TNS_RowOffsets[NewDir-1];

		if (Deviate < 2) {
			tns_RemoveSample(t, TY, TX);
		}elif(Deviate < 3) {
			tns_RemoveSample(t, TY, TX);
			AddPoint = 1;
		}else {
			tns_RemoveSample(t, TY, TX);
			return;
		}

		if (AddPoint || l == Length) {
			p2 = tnsAppendPoint(ls, TX, TY, 0);
			NewP = p2;
			l = 0;
		}
	}
}

int tns_GrowSnakeL(tnsTexture* t, tnsLineStrip* ls, tnsLineStripPoint* ThisP, int Direction) {
	tnsLineStripPoint* NewP = ThisP, *p2;
	int Length = 5;
	int l = 0;
	int Deviate, Dir = Direction, NewDir;
	int AddPoint;
	int TX = NewP->P[0], TY = NewP->P[1];

	while (NewDir = tns_DetectDirection(t, TX, TY, Dir)) {
		AddPoint = 0;
		Deviate = tns_DirectionDeviate(NewDir, Dir);
		Dir = NewDir;

		l++;
		TX += _TNS_ColOffsets[NewDir - 1];
		TY += _TNS_RowOffsets[NewDir - 1];

		if (Deviate < 2) {
			tns_RemoveSample(t, TY, TX);
		}elif(Deviate < 4) {
			tns_RemoveSample(t, TY, TX);
			AddPoint = 1;
		}
		else {
			tns_RemoveSample(t, TY, TX);
			return;
		}

		if (AddPoint || l == Length) {
			p2 = tnsPushPoint(ls, TX, TY, 0);
			NewP = p2;
			l = 0;
		}
	}
}

void tns_CreateSnakeOffscreen(n3DViewUiExtra* e,int W,int H) {
	e->SnakeReadFBO = tnsCreate2DOffscreenBasic(GL_RGBA, W, H, 0);
}
void tns_DestroySnakeOffscreen(n3DViewUiExtra* e, int W, int H) {
	if (!e->SnakeReadFBO) return;
	tnsDelete2DOffscreen(e->SnakeReadFBO);
	e->SnakeReadFBO = 0;
}

void tns_BuildSnakes(n3DViewUiExtra* e){
	int r, c;
	tnsTexture* t = e->SnakeReadFBO->pColorTextures[0];
	tnsLineStrip* ls;
	tnsLineStripPoint* lsp;
	tnsTextureSample* ts;
	int count = 0;
	//tnsFilterKernel Filter = {0};
	//tnsGetGaussianKernel(&Filter, 5, 1);

	tnsReadbackTexture(t);

	while (ls = lstPopItem(&t->LineStrips)) {
		tnsDestroyLineStrip(ls);
	}

	while (ts = tns_AnyUncoveredSample(t)) {
		int Direction=0;
		tnsLineStripPoint tlsp = { 0 };

		tlsp.P[0] = ts->X;
		tlsp.P[1] = ts->Y;

		if (Direction = tns_DetectDirection(t, ts->X,ts->Y, Direction)) {
			lstAppendItem(&t->LineStrips, ls = tnsCreateLineStrip());
			lsp = tnsAppendPoint(ls, ts->X, ts->Y, 0);
			tns_RemoveSample(t, ts->Y, ts->X);

			tns_GrowSnakeR(t, ls, lsp, Direction);

			tns_GrowSnakeL(t, ls, lsp, tns_ReverseDirection(Direction));
		}

		count++;
		
		//if (count > 1000) return;
	}
}


void tns_DrawThisSnake(tnsLineStrip* ls) {
	int Count = ls->PointCount;
	int i;
	tnsLineStripPoint* lsp,*plsp;
	u32bit *Index_adjacent;
	float* Verts;
	float* Lengths;
	GLuint VBuf;
	GLuint IBuf;
	GLuint LBuf;

	float TotalLength=0;

	Index_adjacent = CreateNewBuffer(u32bit, (Count - 1) * 4);
	Verts = CreateNewBuffer(float, Count * 2);
	Lengths = CreateNewBuffer(float, Count);

	for (i = 0; i < Count-1; i++) {
		Index_adjacent[i * 4 + 0] = i - 1;
		Index_adjacent[i * 4 + 1] = i;
		Index_adjacent[i * 4 + 2] = i + 1;
		Index_adjacent[i * 4 + 3] = i + 2;
	}
	Index_adjacent[0] = 0;
	Index_adjacent[(Count - 1) * 4 - 1] = Count-1;

	i = 0;
	for (lsp = ls->Points.pFirst; lsp; lsp = lsp->Item.pNext) {
		Verts[i * 2 + 0] = lsp->P[0];
		Verts[i * 2 + 1] = lsp->P[1];
		if (plsp = lsp->Item.pPrev) {
			TotalLength += tMatDist2dv(plsp->P, lsp->P);
			Lengths[i] = TotalLength;
		}
		i++;
	}

	glGenBuffers(1, &VBuf);
	glBindBuffer(GL_ARRAY_BUFFER, VBuf);
	glBufferData(GL_ARRAY_BUFFER, Count*2 * sizeof(GLfloat), Verts, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(T->LineConnectionShader->vertexIndex);
	glVertexAttribPointer(T->LineConnectionShader->vertexIndex, 2, GL_FLOAT, 0, 0, 0);

	glGenBuffers(1, &LBuf);
	glBindBuffer(GL_ARRAY_BUFFER, LBuf);
	glBufferData(GL_ARRAY_BUFFER, Count * sizeof(GLfloat), Lengths, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(T->LineConnectionShader->uvIndex);
	glVertexAttribPointer(T->LineConnectionShader->uvIndex, 1, GL_FLOAT, 0, 0, 0);

	glGenBuffers(1, &IBuf);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBuf);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, (Count - 1) * 4 * sizeof(GLuint), Index_adjacent, GL_DYNAMIC_DRAW);

	FreeMem(Index_adjacent); FreeMem(Verts); FreeMem(Lengths);

	glDisableVertexAttribArray(T->LineConnectionShader->colorIndex);
	glVertexAttrib4f(T->LineConnectionShader->colorIndex, 1, 1, 0, 1);

	glUniform1f(T->LineConnectionShader->uniform0Index, 1);
	glUniform1f(T->LineConnectionShader->uniform1Index, TotalLength);

	glUniform1f(T->LineConnectionShader->uniform2Index, 00);
	glUniform1f(T->LineConnectionShader->uniform3Index, 00);

	glUniform1f(T->LineConnectionShader->uniform4Index, 0.0);
	glUniform1f(T->LineConnectionShader->uniform5Index, 0.0);


	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElements(GL_LINES_ADJACENCY, (Count - 1) * 4, GL_UNSIGNED_INT, 0);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glDeleteBuffers(1, &VBuf);
	glDeleteBuffers(1, &IBuf);
	glDeleteBuffers(1, &LBuf);
}

void tns_DrawSnakes(n3DViewUiExtra* e) {
	tnsTexture* t = e->SnakeReadFBO->pColorTextures[0];
	tnsLineStrip* ls;
	tnsLineStripPoint* lsp;
	int Count=0;
	int VCount;

	tnsEnableShaderv(T->LineConnectionShader);

	for (ls = t->LineStrips.pFirst; ls; ls = ls->Item.pNext) {
		if (ls->PointCount < 2) continue;

		tns_DrawThisSnake(ls);
		
	}
}

