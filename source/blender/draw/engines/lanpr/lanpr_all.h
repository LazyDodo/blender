#pragma once

#include "lanpr_util.h"
#include "BLI_mempool.h"
#include "BLI_utildefines.h"
//#include "GPU_framebuffer.h"
#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"
#include "DNA_listBase.h"
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

#include "BLI_threads.h"

#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"
#include "bmesh.h"

#include "WM_types.h"
#include "WM_api.h"



#define LANPR_ENGINE "BLENDER_LANPR"

#define TNS_PI 3.1415926535897932384626433832795
#define deg(r) r / TNS_PI * 180.0
#define rad(d) d *TNS_PI / 180.0

#define tMatDist2v(p1, p2) \
	sqrt(((p1)[0] - (p2)[0]) * ((p1)[0] - (p2)[0]) + ((p1)[1] - (p2)[1]) * ((p1)[1] - (p2)[1]))

#define tnsLinearItp(L, R, T) \
	((L) * (1.0f - (T)) + (R)*(T))


typedef struct LANPROneTimeInit {

	/* Snake */
	GPUShader *multichannel_shader;
	GPUShader *edge_detect_shader;
	GPUShader *edge_thinning_shader;
	GPUShader *snake_connection_shader;

	/* DPIX */
	GPUShader *dpix_transform_shader;
	GPUShader *dpix_preview_shader;

	/* Software */
	GPUShader *software_shader;

	void *ved;


	/* for default value assignment */

	int InitComplete;

} LANPROneTimeInit;

#define TNS_DPIX_TEXTURE_SIZE 2048

typedef struct LANPR_TextureSample {
	Link Item;
	int X, Y;
	float Z;    // for future usage
} LANPR_TextureSample;

typedef struct LANPR_LineStripPoint {
	Link Item;
	float P[3];
} LANPR_LineStripPoint;

typedef struct LANPR_LineStrip {
	Link Item;
	ListBase points;
	int point_count;
	float total_length;
}LANPR_LineStrip;

typedef struct LANPR_PassList {
	/* Snake */
	struct DRWPass *depth_pass;
	struct DRWPass *color_pass;
	struct DRWPass *normal_pass;
	struct DRWPass *edge_intermediate;
	struct DRWPass *edge_thinning;
	struct DRWPass *snake_pass;

	/* DPIX */
	struct DRWPass *dpix_transform_pass;
	struct DRWPass *dpix_preview_pass;

	/* SOFTWARE */
	struct DRWPass *software_pass;

} LANPR_PassList;

typedef struct LANPR_FramebufferList {

	/* Snake */
	struct GPUFrameBuffer *passes;
	struct GPUFrameBuffer *edge_intermediate;
	struct GPUFrameBuffer *edge_thinning;

	/* DPIX */
	struct GPUFrameBuffer *dpix_transform;
	struct GPUFrameBuffer *dpix_preview;

	/* Software */
	struct GPUFrameBuffer *software_ms;

} LANPR_FramebufferList;

typedef struct LANPR_TextureList {

	struct GPUTexture *color;
	struct GPUTexture *normal;
	struct GPUTexture *depth;
	struct GPUTexture *edge_intermediate;

	struct GPUTexture *dpix_in_pl;/* point L */
	struct GPUTexture *dpix_in_pr;/* point R */
	struct GPUTexture *dpix_in_nl;/* normal L */
	struct GPUTexture *dpix_in_nr;/* normal R */
	struct GPUTexture *dpix_in_edge_mask;/* RGBA, R:Material, G: Freestyle Edge Mark, BA:Reserved for future usage */

	struct GPUTexture *dpix_out_pl;
	struct GPUTexture *dpix_out_pr;
	struct GPUTexture *dpix_out_length;

	/* multisample resolve */
	struct GPUTexture *ms_resolve_depth;
	struct GPUTexture *ms_resolve_color;

} LANPR_TextureList;

typedef struct LANPR_PrivateData {
	DRWShadingGroup *multipass_shgrp;
	DRWShadingGroup *edge_detect_shgrp;
	DRWShadingGroup *edge_thinning_shgrp;
	DRWShadingGroup *snake_shgrp;

	DRWShadingGroup *dpix_transform_shgrp;
	DRWShadingGroup *dpix_preview_shgrp;

	// moved into line layer.
	//DRWShadingGroup *software_shgrp;

	//snake

	float normal_clamp;
	float normal_strength;
	float depth_clamp;
	float depth_strength;

	float zfar;
	float znear;

	int stage;//thinning

	float         *line_result;
	unsigned char *line_result_8bit;
	int width, height;          // if not match recreate buffer.
	void         **sample_table;

	BLI_mempool *mp_sample;
	BLI_mempool *mp_line_strip;
	BLI_mempool *mp_line_strip_point;
	BLI_mempool *mp_batch_list;

	ListBase pending_samples;
	ListBase erased_samples;
	ListBase line_strips;

	// dpix data

	void *atlas_pl;
	void *atlas_pr;
	void *atlas_nl;
	void *atlas_nr;
	void *atlas_edge_mask;

	int begin_index;

	int dpix_sample_step;
	int dpix_is_perspective;
	float dpix_viewport[4];
	float output_viewport[4];
	int dpix_buffer_width;
	float dpix_depth_offset;

	float dpix_znear;
	float dpix_zfar;

	// drawing

	unsigned v_buf;
	unsigned i_buf;
	unsigned l_buf;

	GPUVertFormat snake_gwn_format;
	GPUBatch *snake_batch;

	ListBase dpix_batch_list;

} LANPR_PrivateData;

typedef struct LANPR_StorageList {
	LANPR_PrivateData *g_data;
} LANPR_StorageList;

typedef struct LANPR_BatchItem {
	Link Item;
	GPUBatch *dpix_transform_batch;
	GPUBatch *dpix_preview_batch;
	Object *ob;
} LANPR_BatchItem;

typedef struct LANPR_Data {
	void *engine_type;
	LANPR_FramebufferList *fbl;
	LANPR_TextureList *txl;
	LANPR_PassList *psl;
	LANPR_StorageList *stl;
} LANPR_Data;



/* Below ported from NUL_TNS.h */

#define TNS_THREAD_LINE_COUNT 10000

#define TNS_CALCULATION_IDLE         0
#define TNS_CALCULATION_GEOMETRY     1
#define TNS_CALCULATION_CONTOUR      2
#define TNS_CALCULATION_INTERSECTION 3
#define TNS_CALCULATION_OCCLUTION    4
#define TNS_CALCULATION_FINISHED     100

typedef struct LANPR_RenderTaskInfo {
	//thrd_t           ThreadHandle;

	struct LANPR_RenderBuffer *RenderBuffer;
	int ThreadID;

	struct nListItemPointer *Contour;
	nListHandle ContourPointers;

	struct nListItemPointer *Intersection;
	nListHandle IntersectionPointers;

	struct nListItemPointer *Crease;
	nListHandle CreasePointers;

	struct nListItemPointer *Material;
	nListHandle MaterialPointers;

} LANPR_RenderTaskInfo;

typedef struct LANPR_RenderBuffer {
	struct LANPR_RenderBuffer *prev, *next;

	//nSafeString*       Name;

	//tnsFrameBuffer*    FrameBuffer;
	//now we move frame buffer content here
	int W, H;
	int SubPixelSample;//1,2,3,4, Use Squared Value.
	int TileSizeW, TileSizeH;
	int TileCountX, TileCountY;
	real WidthPerTile, HeightPerTile;
	tnsMatrix44d ViewProjection;
	tnsMatrix44d VPInverse;

	nSafeString *OutputFolder;   //end with a slash;
	nSafeString *ImagePrefix;
	nSafeString *ImageNameConnector;

	int OutputMode;
	int OutputAALevel;

	struct LANPR_BoundingArea *InitialBoundingAreas;
	u32bit BoundingAreaCount;
	//u32bit             BaVBO;
	//u32bit           BaFillVBO;

	nListHandle VertexBufferPointers;
	nListHandle LineBufferPointers;
	nListHandle TriangleBufferPointers;
	nListHandle AllRenderLines;

	nListHandle IntersectingVertexBuffer;

	struct GPUBatch *DPIXIntersectionTransformBatch;
	struct GPUBatch *DPIXIntersectionBatch;

	/* use own-implemented one */
	nStaticMemoryPool RenderDataPool;

	Material           *MaterialPointers[2048];

	//render status

	real ViewVector[3];

	int TriangleSize;

	u32bit ContourCount;
	u32bit ContourProcessed;
	nListItemPointer *ContourManaged;
	nListHandle Contours;

	u32bit IntersectionCount;
	u32bit IntersectionProcessed;
	nListItemPointer *IntersectionManaged;
	nListHandle IntersectionLines;

	u32bit CreaseCount;
	u32bit CreaseProcessed;
	nListItemPointer *CreaseManaged;
	nListHandle CreaseLines;

	u32bit MaterialLineCount;
	u32bit MaterialProcessed;
	nListItemPointer *MaterialManaged;
	nListHandle MaterialLines;

	u32bit EdgeMarkCount;
	u32bit EdgeMarkProcessed;
	nListItemPointer *EdgeMarkManaged;
	nListHandle EdgeMarks;

	nListHandle Chains;
	GPUBatch*  ChainDrawBatch;
	DRWShadingGroup* ChainShgrp;

	SpinLock csInfo;
	SpinLock csData;
	SpinLock csManagement;

	//settings

	//int             OutputTransparent;
	//real            BackgroundColor[4];

	int MaxOccludeLevel;
	real CreaseAngle;
	real CreaseCos;
	int CreaseAllowOverride;
	int ThreadCount;

	real OverallProgress;
	int CalculationStatus;

	int DrawMaterialPreview;
	real MaterialTransparency;

	int ShowLine;
	int ShowFast;
	int ShowMaterial;
	int OverrideDisplay;

	//nListHandle DrawCommands; /* not used now */

	//tnsRenderBufferPreviewNode RenderPreview[32];

	struct Scene *Scene;
	struct Object *Camera;

	int calculate_intersections;
	int size;

	//tnsRenderTriangles are in mesh object.
}LANPR_RenderBuffer;


#define TNS_CULL_DISCARD 2
#define TNS_CULL_USED    1

typedef struct LANPR_RenderTriangle {
	nListItem Item;
	struct LANPR_RenderVert *V[3];
	struct LANPR_RenderLine *RL[3];
	real GN[3];
	real GC[3];
	//struct BMFace *F;
	short MaterialID;
	nListHandle IntersectingVerts;
	char CullStatus;
	struct LANPR_RenderTriangle *Testing;   //Should Be tRT** Testing[NumOfThreads]
}LANPR_RenderTriangle;

typedef struct LANPR_RenderTriangleThread {
	struct LANPR_RenderTriangle Base;
	struct LANPR_RenderLine *Testing[128];    //max thread support;
}LANPR_RenderTriangleThread;

typedef struct LANPR_RenderElementLinkNode {
	nListItem Item;
	void *Pointer;
	int ElementCount;
	void *ObjectRef;
	char Additional;
}LANPR_RenderElementLinkNode;

typedef struct LANPR_RenderLineSegment {
	nListItem Item;
	//real     Begin, End;  // 0->At[L] 1->At[R]
	real at;
	u8bit OccludeLevel;    //after
}LANPR_RenderLineSegment;

typedef struct LANPR_RenderVert {
	nListItem Item;
	real GLocation[4];
	real FrameBufferCoord[4];
	int FrameBufferCoordi[2];
	struct BMVert *V;              //Used As R When Intersecting
	struct LANPR_RenderLine *IntersectingLine;
	struct LANPR_RenderLine *IntersectintLine2;
	struct LANPR_RenderTriangle *IntersectWith;     //   Positive 1         Negative 0
	//tnsRenderTriangle* IntersectingOnFace;       //         <|               |>
	char Positive;                    //                 L---->|----->R	 L---->|----->R
	char EdgeUsed;                    //                      <|		       |>
}LANPR_RenderVert;

#define LANPR_EDGE_FLAG_EDGE_MARK    1
#define LANPR_EDGE_FLAG_CONTOUR      2
#define LANPR_EDGE_FLAG_CREASE       4
#define LANPR_EDGE_FLAG_MATERIAL     8
#define LANPR_EDGE_FLAG_INTERSECTION 16
#define LANPR_EDGE_FLAG_FLOATING     32 // floating edge, unimplemented yet
#define LANPR_EDGE_FLAG_CHAIN_PICKED 64

#define LANPR_EDGE_FLAG_ALL_TYPE     0x3f

typedef struct LANPR_RenderLine {
	nListItem Item;
	struct LANPR_RenderVert *L, *R;
	struct LANPR_RenderTriangle *TL, *TR;
	nListHandle Segments;
	//tnsEdge*       Edge;//should be edge material
	//tnsRenderTriangle* Testing;//Should Be tRT** Testing[NumOfThreads]
	char MinOcclude;
	char Flags; // also for line type determination on chainning
	struct Object *ObjectRef;
}LANPR_RenderLine;

typedef struct LANPR_RenderLineChain {
	nListItem   Item;
	nListHandle Chain;
	//int         SegmentCount;  // we count before draw cmd.
	float       Length;          // calculated before draw cmd.
}LANPR_RenderLineChain;

typedef struct LANPR_RenderLineChainItem {
	nListItem   Item;
	float       pos[3]; // need z value for fading
	char        LineType;      //      style of [1]       style of [2]
	char        OccludeLevel;  // [1]--------------->[2]---------------->[3]--....
}LANPR_RenderLineChainItem;

typedef struct LANPR_BoundingArea {
	real L, R, U, B;
	real CX, CY;

	struct LANPR_BoundingArea *Child;//1,2,3,4 quadrant

	nListHandle LP;
	nListHandle RP;
	nListHandle UP;
	nListHandle BP;

	int TriangleCount;
	nListHandle LinkedTriangles;
	nListHandle LinkedLines;
	nListHandle LinkedChains;//reserved for multithread chainning
}LANPR_BoundingArea;


#define TNS_COMMAND_LINE     0
#define TNS_COMMAND_MATERIAL 1
#define TNS_COMMAND_EDGE     2

#define TNS_TRANSPARENCY_DRAW_SIMPLE  0
#define TNS_TRANSPARENCY_DRAW_LAYERED 1

#define TNS_OVERRIDE_ONLY                     0
#define TNS_OVERRIDE_EXCLUDE                  1
//#define TNS_OVERRIDE_ALL_OTHERS_OUTSIDE_GROUP 2
//#define TNS_OVERRIDE_ALL_OTHERS_IN_GROUP      3
//#define TNS_OVERRIDE_ALL_OTHERS               4

extern RenderEngineType DRW_engine_viewport_lanpr_type;



#define tnsLinearItp(L, R, T) \
	((L) * (1.0f - (T)) + (R)*(T))


#define TNS_TILE(tile, r, c, CCount) \
	tile[r * CCount + c]

#define TNS_CLAMP(a, Min, Max) \
	a = a < Min ? Min : (a > Max ? Max : a)

#define TNS_MAX2_INDEX(a, b) \
	(a > b ? 0 : 1)

#define TNS_MIN2_INDEX(a, b) \
	(a < b ? 0 : 1)

#define TNS_MAX3_INDEX(a, b, c) \
	(a > b ? (a > c ? 0 : (b > c ? 1 : 2)) : (b > c ? 1 : 2))

#define TNS_MIN3_INDEX(a, b, c) \
	(a < b ? (a < c ? 0 : (b < c ? 1 : 2)) : (b < c ? 1 : 2))

#define TNS_MAX3_INDEX_ABC(x, y, z) \
	(x > y ? (x > z ? a : (y > z ? b : c)) : (y > z ? b : c))

#define TNS_MIN3_INDEX_ABC(x, y, z) \
	(x < y ? (x < z ? a : (y < z ? b : c)) : (y < z ? b : c))

#define TNS_ABC(index) \
	(index == 0 ? a : (index == 1 ? b : c))


#define TNS_DOUBLE_CLOSE_ENOUGH(a, b) \
	(((a) + DBL_EDGE_LIM) >= (b) && ((a) - DBL_EDGE_LIM) <= (b))

//#define TNS_DOUBLE_CLOSE_ENOUGH(a,b)\
// //(((a)+0.00000000001)>=(b) && ((a)-0.0000000001)<=(b))


#define TNS_FLOAT_CLOSE_ENOUGH_WIDER(a, b) \
	(((a) + 0.0000001) >= (b) && ((a) - 0.0000001) <= (b))

#define TNS_IN_TILE_X(RenderTile, Fx) \
	(RenderTile->FX <= Fx && RenderTile->FXLim >= Fx)

#define TNS_IN_TILE_Y(RenderTile, Fy) \
	(RenderTile->FY <= Fy && RenderTile->FYLim >= Fy)


#define TNS_IN_TILE(RenderTile, Fx, Fy) \
	(TNS_IN_TILE_X(RenderTile, Fx) && TNS_IN_TILE_Y(RenderTile, Fy))

__inline void tMatConvert44df(tnsMatrix44d from, tnsMatrix44f to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
	to[3] = from[3];
	to[4] = from[4];
	to[5] = from[5];
	to[6] = from[6];
	to[7] = from[7];
	to[8] = from[8];
	to[9] = from[9];
	to[10] = from[10];
	to[11] = from[11];
	to[12] = from[12];
	to[13] = from[13];
	to[14] = from[14];
	to[15] = from[15];
}

__inline int lanpr_TrangleLineBoundBoxTest(LANPR_RenderTriangle *rt, LANPR_RenderLine *rl) {
	if (MAX3(rt->V[0]->FrameBufferCoord[2], rt->V[1]->FrameBufferCoord[2], rt->V[2]->FrameBufferCoord[2]) > MIN2(rl->L->FrameBufferCoord[2], rl->R->FrameBufferCoord[2])) return 0;
	if (MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]) < MIN2(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0])) return 0;
	if (MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]) > MAX2(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0])) return 0;
	if (MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]) < MIN2(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1])) return 0;
	if (MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]) > MAX2(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1])) return 0;
	return 1;
}

__inline double tMatGetLinearRatio(real L, real R, real FromL);
__inline int lanpr_LineIntersectTest2d(tnsVector2d a1, tnsVector2d a2, tnsVector2d b1, tnsVector2d b2, double *aRatio) {
	double k1, k2;
	double x;
	double y;
	double Ratio;
	double xDiff = (a2[0] - a1[0]);// +DBL_EPSILON;
	double xDiff2 = (b2[0] - b1[0]);

	if (xDiff == 0) {
		if (xDiff2 == 0) {
			*aRatio = 0;
			return 0;
		}
		double r2 = tMatGetLinearRatio(b1[0], b2[0], a1[0]);
		y = tnsLinearItp(b1[1], b2[1], r2);
		*aRatio = Ratio = tMatGetLinearRatio(a1[1], a2[1], y);
	}
	else {
		if (xDiff2 == 0) {
			Ratio = tMatGetLinearRatio(a1[0], a2[0], b1[0]);
			//y = tnsLinearItp(a1[1], a2[1], r2);
			*aRatio = Ratio;
		}
		else {
			k1 = (a2[1] - a1[1]) / xDiff;
			k2 = (b2[1] - b1[1]) / xDiff2;

			if ((k1 == k2))
				return 0;

			x = (a1[1] - b1[1] - k1 * a1[0] + k2 * b1[0]) / (k2 - k1);

			Ratio = (x - a1[0]) / xDiff;

			*aRatio = Ratio;
		}
	}



	if (b1[0] == b2[0]) {
		y = tnsLinearItp(a1[1], a2[1], Ratio);
		if (y > MAX2(b1[1], b2[1]) || y < MIN2(b1[1], b2[1])) return 0;
	}
	else
	if (Ratio <= 0 || Ratio > 1 ||
	    (b1[0] > b2[0] && x > b1[0]) ||
	    (b1[0] < b2[0] && x < b1[0]) ||
	    (b2[0] > b1[0] && x > b2[0]) ||
	    (b2[0] < b1[0] && x < b2[0]))
		return 0;

	return 1;
}
__inline double lanpr_GetLineZ(tnsVector3d L, tnsVector3d R, real Ratio) {
	//double z = 1 / tnsLinearItp(1 / L[2], 1 / R[2], Ratio);
	double z = tnsLinearItp(L[2], R[2], Ratio);
	return z;
}
__inline double lanpr_GetLineZPoint(tnsVector3d L, tnsVector3d R, tnsVector3d FromL) {
	double r = (FromL[0] - L[0]) / (R[0] - L[0]);
	return tnsLinearItp(L[2], R[2], r);
	//return 1 / tnsLinearItp(1 / L[2], 1 / R[2], r);
}
__inline double lanpr_GetLinearRatio(tnsVector3d L, tnsVector3d R, tnsVector3d FromL) {
	double r = (FromL[0] - L[0]) / (R[0] - L[0]);
	return r;
}

__inline double tMatGetLinearRatio(real L, real R, real FromL) {
	double r = (FromL - L) / (R - L);
	return r;
}
__inline void tMatVectorMinus2d(tnsVector2d result, tnsVector2d l, tnsVector2d r) {
	result[0] = l[0] - r[0];
	result[1] = l[1] - r[1];
}

__inline void tMatVectorMinus3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[0] - r[0];
	result[1] = l[1] - r[1];
	result[2] = l[2] - r[2];
}
__inline void tMatVectorSubtract3d(tnsVector3d l, tnsVector3d r) {
	l[0] = l[0] - r[0];
	l[1] = l[1] - r[1];
	l[2] = l[2] - r[2];
}
__inline void tMatVectorPlus3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[0] + r[0];
	result[1] = l[1] + r[1];
	result[2] = l[2] + r[2];
}
__inline void tMatVectorAccum3d(tnsVector3d l, tnsVector3d r) {
	l[0] = l[0] + r[0];
	l[1] = l[1] + r[1];
	l[2] = l[2] + r[2];
}
__inline void tMatVectorAccum2d(tnsVector2d l, tnsVector2d r) {
	l[0] = l[0] + r[0];
	l[1] = l[1] + r[1];
}
__inline void tMatVectorNegate3d(tnsVector3d result, tnsVector3d l) {
	result[0] = -l[0];
	result[1] = -l[1];
	result[2] = -l[2];
}
__inline void tMatVectorNegateSelf3d(tnsVector3d l) {
	l[0] = -l[0];
	l[1] = -l[1];
	l[2] = -l[2];
}
__inline void tMatVectorCopy2d(tnsVector2d from, tnsVector2d to) {
	to[0] = from[0];
	to[1] = from[1];
}
__inline void tMatVectorCopy3d(tnsVector3d from, tnsVector3d to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
}
__inline void tMatVectorCopy4d(tnsVector4d from, tnsVector4d to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
	to[3] = from[3];
}
__inline void tMatVectorMultiSelf4d(tnsVector3d from, real num) {
	from[0] *= num;
	from[1] *= num;
	from[2] *= num;
	from[3] *= num;
}
__inline void tMatVectorMultiSelf3d(tnsVector3d from, real num) {
	from[0] *= num;
	from[1] *= num;
	from[2] *= num;
}
__inline void tMatVectorMultiSelf2d(tnsVector3d from, real num) {
	from[0] *= num;
	from[1] *= num;
}

__inline real tMatDirectionToRad(tnsVector2d Dir) {
	real arcc = acos(Dir[0]);
	real arcs = asin(Dir[1]);

	if (Dir[0] >= 0) {
		if (Dir[1] >= 0) return arcc;
		else return TNS_PI * 2 - arcc;
	}
	else {
		if (Dir[1] >= 0) return arcs + TNS_PI / 2;
		else return TNS_PI + arcs;
	}
}


__inline void tMatVectorConvert4fd(tnsVector4f from, tnsVector4d to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
	to[3] = from[3];
}

__inline void tMatVectorConvert3fd(tnsVector3f from, tnsVector3d to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
}


int lanpr_PointInsideTrianglef(tnsVector2d v, tnsVector2d v0, tnsVector2d v1, tnsVector2d v2);
real lanpr_LinearInterpolate(real L, real R, real T);
void lanpr_LinearInterpolate2dv(real *L, real *R, real T, real *Result);
void lanpr_LinearInterpolate3dv(real *L, real *R, real T, real *Result);



// functions

// dpix

void lanpr_init_atlas_inputs(void *ved);
void lanpr_destroy_atlas(void *ved);
int lanpr_feed_atlas_data_obj(void *vedata,
                              float *AtlasPointsL, float *AtlasPointsR,
                              float *AtlasFaceNormalL, float *AtlasFaceNormalR,
                              float *AtlasEdgeMask,
                              Object *ob, int BeginIndex);

int lanpr_feed_atlas_data_intersection_cache(void *vedata,
                                             float *AtlasPointsL, float *AtlasPointsR,
                                             float *AtlasFaceNormalL, float *AtlasFaceNormalR,
                                             float *AtlasEdgeMask,
                                             int BeginIndex);

int lanpr_feed_atlas_trigger_preview_obj(void *vedata, Object *ob, int BeginIndex);
void lanpr_create_atlas_intersection_preview(void *vedata, int BeginIndex);

void lanpr_dpix_draw_scene(LANPR_TextureList *txl, LANPR_FramebufferList *fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr, GPUFrameBuffer *DefaultFB);

void lanpr_snake_draw_scene(LANPR_TextureList *txl, LANPR_FramebufferList *fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr, GPUFrameBuffer *DefaultFB);
