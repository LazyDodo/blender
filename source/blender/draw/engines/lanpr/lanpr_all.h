#pragma once

#include "NUL_Util.h"
#include "BLI_mempool.h"
#include "GPU_framebuffer.h"
#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"



#define LANPR_ENGINE "BLENDER_LANPR"

#define TNS_PI 3.1415926535897932384626433832795
#define deg(r) r/TNS_PI*180.0
#define rad(d) d*TNS_PI/180.0

#define tMatDist2v(p1,p2)\
    sqrt(((p1)[0]-(p2)[0])*((p1)[0]-(p2)[0]) + ((p1)[1]-(p2)[1])*((p1)[1]-(p2)[1]))

#define tnsLinearItp(L,R,T)\
((L)*(1.0f - (T)) + (R)*(T))


typedef struct LANPROneTimeInit{
    
	/* Snake */

	GPUShader* multichannel_shader;
	GPUShader* edge_detect_shader;
	GPUShader* edge_thinning_shader;
	GPUShader* snake_connection_shader;

	/* DPIX */

	GPUShader* dpix_transform_shader;
	GPUShader* dpix_preview_shader;

	void* ved;


    /* for default value assignment */

    int        InitComplete;
    
} LANPROneTimeInit;

#define TNS_DPIX_TEXTURE_SIZE 2048

typedef struct LANPR_TextureSample {
	Link      Item;
	int       X,Y;
	float     Z;// for future usage
} LANPR_TextureSample;

typedef struct LANPR_LineStripPoint {
	Link     Item;
	float P[3];
} LANPR_LineStripPoint;

typedef struct LANPR_LineStrip{
	Link     Item;
	ListBase points;
	int      point_count;
	float    total_length;
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

} LANPR_PassList;

typedef struct LANPR_FramebufferList {

	/* Snake */
	struct GPUFrameBuffer *passes;
	struct GPUFrameBuffer *edge_intermediate;
	struct GPUFrameBuffer *edge_thinning;

    /* DPIX */
	struct GPUFrameBuffer *dpix_transform;
	struct GPUFrameBuffer *dpix_preview;

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

} LANPR_TextureList;

typedef struct LANPR_PrivateData {
	DRWShadingGroup *multipass_shgrp;
	DRWShadingGroup *edge_detect_shgrp;
	DRWShadingGroup *edge_thinning_shgrp;
    DRWShadingGroup *snake_shgrp;
	
	DRWShadingGroup *dpix_transform_shgrp;
	DRWShadingGroup *dpix_preview_shgrp;

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
	int            width,height;// if not match recreate buffer.
	void         **sample_table;

	BLI_mempool*  mp_sample;
	BLI_mempool*  mp_line_strip;
	BLI_mempool*  mp_line_strip_point;
	BLI_mempool*  mp_batch_list;
	
	ListBase      pending_samples;
	ListBase      erased_samples;
    ListBase      line_strips;

	// dpix data

	void*         atlas_pl;
	void*         atlas_pr;
	void*         atlas_nl;
	void*         atlas_nr;
	void*         atlas_edge_mask;

	int           begin_index;

	int           dpix_sample_step;
	int           dpix_is_perspective;
	float         dpix_viewport[4];
    int           dpix_buffer_width;
	float         dpix_depth_offset;

	float         dpix_znear;
	float         dpix_zfar;

	// drawing

	unsigned        v_buf;
	unsigned        i_buf;
	unsigned        l_buf;
    
	Gwn_VertFormat   snake_gwn_format;
	Gwn_Batch*       snake_batch;

	ListBase         dpix_batch_list;

} LANPR_PrivateData;

typedef struct LANPR_StorageList {
	LANPR_PrivateData *g_data;
} LANPR_StorageList;

typedef struct LANPR_BatchItem {
	Link             Item;
	Gwn_Batch*       dpix_transform_batch;
	Gwn_Batch*       dpix_preview_batch;
	Object*          ob;
} LANPR_BatchItem;

typedef struct LANPR_Data {
	void *engine_type;
	LANPR_FramebufferList *fbl;
	LANPR_TextureList *txl;
	LANPR_PassList *psl;
	LANPR_StorageList *stl;
} LANPR_Data;



typedef struct tnsRenderTaskInfo {
	//thrd_t           ThreadHandle;

	//tnsRenderBuffer* RenderBuffer;
	int              ThreadID;

	//nListItemPointer* Contour;
	//nListHandle       ContourPointers;

	//nListItemPointer* Intersection;
	//nListHandle       IntersectionPointers;

	//nListItemPointer* Crease;
	//nListHandle       CreasePointers;

	//nListItemPointer* Material;
	//nListHandle       MaterialPointers;

} tnsRenderTaskInfo;

/* Below ported from NUL_TNS.h */

typedef struct LANPR_RenderBuffer {
	struct LANPR_RenderBuffer *prev, *next;

	//nSafeString*       Name;

	//tnsFrameBuffer*    FrameBuffer;

	//tnsBoundingArea*   InitialBoundingAreas;
	//u32bit             BoundingAreaCount;
	//u32bit             BaVBO;
	//u32bit           BaFillVBO;

	ListBase           VertexBufferPointers;
	ListBase           TriangleBufferPointers;
	
	ListBase           AllRenderLines;

	//nListHandle     IntersectingVertexBuffer;

    /* BLI_xx equ? */
	//nStaticMemoryPool  RenderDataPool;

	//render status

	real            ViewVector[3];

	int             TriangleSize;

	u32bit          ContourCount;
	u32bit          ContourProcessed;
	//nListItemPointer* ContourManaged;
	ListBase        Contours;

	u32bit          IntersectionCount;
	u32bit          IntersectionProcessed;
	//nListItemPointer* IntersectionManaged;
	ListBase        IntersectionLines;

	u32bit          CreaseCount;
	u32bit          CreaseProcessed;
	//nListItemPointer* CreaseManaged;
	ListBase        CreaseLines;

	u32bit          MaterialLineCount;
	u32bit          MaterialProcessed;
	//nListItemPointer* MaterialManaged;
	ListBase        MaterialLines;

	//CRITICAL_SECTION csInfo;
	//CRITICAL_SECTION csData;
	//CRITICAL_SECTION csManagement;

	//settings

	//int             OutputTransparent;
	//real            BackgroundColor[4];

	int             MaxOccludeLevel;
	real            CreaseAngle;
	real            CreaseCos;
	int             CreaseAllowOverride;
	int             ThreadCount;
	
	real            OverallProgress;
	int             CalculationStatus;

	int             DrawMaterialPreview;
	real            MaterialTransparency;

	int             ShowLine;
	int             ShowFast;
	int             ShowMaterial;
	int             OverrideDisplay;

	ListBase        DrawCommands;

	//tnsRenderBufferPreviewNode RenderPreview[32];

	Scene*          Scene;
	//tnsCamera* Camera;

	//tnsRenderTriangles are in mesh object.
}LANPR_RenderBuffer;


#define TNS_CULL_DISCARD 2
#define TNS_CULL_USED    1

typedef struct LANPR_RenderTriangle {
	Link              Item;
	struct LANPR_RenderVert* V[3];
	struct LANPR_RenderLine* RL[3];
	real              GN[3];
	real              GC[3];
	struct BMFace*           F;
	ListBase          IntersectingVerts;
	char              CullStatus;
	struct LANPR_RenderLine* Testing;	//Should Be tRT** Testing[NumOfThreads]
}LANPR_RenderTriangle;

typedef struct LANPR_RenderTriangleThread {
	struct LANPR_RenderTriangle Base;
	struct LANPR_RenderLine*    Testing[128]; //max thread support;
}LANPR_RenderTriangleThread;

typedef struct LANPR_RenderElementLinkNode {
	Link      Item;
	void*     Pointer;
	int       ElementCount;
	void*     ObjectRef;
	char      Additional;
}LANPR_RenderElementLinkNode;

typedef struct LANPR_tnsRenderLineSegment {
	Link      Item;
	//real     Begin, End;  // 0->At[L] 1->At[R]
	real      at;
	u8bit     OccludeLevel;//after
	int       PreviewIndex;
}LANPR_tnsRenderLineSegment;

struct LANPR_RenderVert{
	Link   Item;
	real GLocation[4];
	real FrameBufferCoord[4];
	int FrameBufferCoordi[2];
	struct BMVert*    V;           //Used As R When Intersecting
	struct LANPR_RenderLine*     IntersectingLine;
	struct LANPR_RenderVert*     IntersectintLine2;
	struct LANPR_RenderTriangle* IntersectWith;     //   Positive 1         Negative 0
	//tnsRenderTriangle* IntersectingOnFace;       //         <|               |>
	char        Positive;             //                 L---->|----->R	 L---->|----->R
	char        EdgeUsed;             //                      <|		       |>
}LANPR_RenderVert;

typedef struct LANPR_RenderLine {
	Link            Item;
	struct LANPR_RenderVert *L, *R;
	struct LANPR_RenderTriangle *TL, *TR;
	ListBase        Segments;
	//tnsEdge*       Edge;//should be edge material
	//tnsRenderTriangle* Testing;//Should Be tRT** Testing[NumOfThreads]
	char            MinOcclude;
	struct Object*         ObjectRef;
	//char            IgnoreConnectedFace;
	//char            CullStatus;
}LANPR_RenderLine;

typedef struct LANPR_BoundingArea {
	real        L, R, U, B;
	real        CX, CY;
	
	struct LANPR_BoundingArea* Child;//1,2,3,4 quadrant

	ListBase    LP;
	ListBase    RP;
	ListBase    UP;
	ListBase    BP;

	int         TriangleCount;
	ListBase    AssociatedTriangles;
}LANPR_BoundingArea;

typedef struct LANPR_RenderSubPixel {
	real                  Depth;
	struct LANPR_RenderTriangle* BelongTo;
	real                  Weight[3];  //belongto->vp 1 2 3
}LANPR_RenderSubPixel;

typedef struct LANPR_RenderTile {
	int                Row, Column;
	int                SubX, SubY, SubXLim, SubYLim;//lower Left Corner As 0
	real               FX, FY, FXLim, FYLim;  //ratio;
	//LANPR_RenderSubPixel* FirstPixel;            //lower Left Corner As 0
	ListBase           AssociatedTriangles;   //lstptrs
	ListBase           AssociatedLines;       //lstptrs
	char               Rendered;
}LANPR_RenderTile;


extern RenderEngineType DRW_engine_viewport_lanpr_type;



#define LANPR_MASTER_MODE_DPIX         0
#define LANPR_MASTER_MODE_SNAKE        1

#define LANPR_POST_PROCESSING_DISABLED 0
#define LANPR_POST_PROCESSING_ENABLED  1

#define LANPR_USE_DIFFERENT_TAPER      0
#define LANPR_USE_SAME_TAPER           1

#define LANPR_DISABLE_TIP_EXTEND       0
#define LANPR_ENABLE_TIP_EXTEND        1



#define TNS_TILE(tile,r,c,CCount)\
tile[r*CCount+c]

#define TNS_CLAMP(a,Min,Max)\
a=a<Min?Min:(a>Max?Max:a)

#define TNS_MAX2(a,b)\
(a>b?a:b)

#define TNS_MIN2(a,b)\
(a<b?a:b)

#define TNS_MAX3(a,b,c)\
(a>TNS_MAX2(b,c)?a:TNS_MAX2(b,c))

#define TNS_MIN3(a,b,c)\
(a<TNS_MIN2(b,c)?a:TNS_MIN2(b,c))

#define TNS_MAX2_INDEX(a,b)\
(a>b?0:1)

#define TNS_MIN2_INDEX(a,b)\
(a<b?0:1)

#define TNS_MAX3_INDEX(a,b,c)\
(a>b?(a>c?0:(b>c?1:2)):(b>c?1:2))

#define TNS_MIN3_INDEX(a,b,c)\
(a<b?(a<c?0:(b<c?1:2)):(b<c?1:2))

#define TNS_MAX3_INDEX_ABC(x,y,z)\
(x>y?(x>z?a:(y>z?b:c)):(y>z?b:c))

#define TNS_MIN3_INDEX_ABC(x,y,z)\
(x<y?(x<z?a:(y<z?b:c)):(y<z?b:c))

#define TNS_ABC(index)\
(index==0?a:(index==1?b:c))


#define TNS_DOUBLE_CLOSE_ENOUGH(a,b)\
(((a)+DBL_EDGE_LIM)>=(b) && ((a)-DBL_EDGE_LIM)<=(b))

//#define TNS_DOUBLE_CLOSE_ENOUGH(a,b)\
//(((a)+0.00000000001)>=(b) && ((a)-0.0000000001)<=(b))


#define TNS_FLOAT_CLOSE_ENOUGH_WIDER(a,b)\
(((a)+0.0000001)>=(b) && ((a)-0.0000001)<=(b))


#define TNS_FRAMEBUFFER_PIXEL(FrameBuffer,Row,Column)\
&((FrameBuffer)->Pixels[Row*FrameBuffer->TileSizeW*FrameBuffer->W*FrameBuffer->SubPixelSample + Column*FrameBuffer->H*FrameBuffer->TileSizeH*FrameBuffer->SubPixelSample])

#define TNS_IN_TILE_X(RenderTile,Fx)\
(RenderTile->FX<=Fx && RenderTile->FXLim>=Fx)

#define TNS_IN_TILE_Y(RenderTile,Fy)\
(RenderTile->FY<=Fy && RenderTile->FYLim>=Fy)

#define TNS_IN_TILE(RenderTile,Fx,Fy)\
(TNS_IN_TILE_X(RenderTile,Fx) && TNS_IN_TILE_Y(RenderTile,Fy))



// functions 

// dpix

void lanpr_init_atlas_inputs(void *ved);
void lanpr_destroy_atlas(void *ved);
int lanpr_feed_atlas_data_obj(void* vedata,
	float* AtlasPointsL, float* AtlasPointsR,
	float* AtlasFaceNormalL, float* AtlasFaceNormalR,
	float* AtlasEdgeMask,
	Object* ob, int BeginIndex);
int lanpr_feed_atlas_trigger_preview_obj(void* vedata, Object* ob, int BeginIndex);

void lanpr_dpix_draw_scene(LANPR_TextureList* txl, LANPR_FramebufferList * fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr);


//snake

void lanpr_snake_draw_scene(LANPR_TextureList* txl, LANPR_FramebufferList * fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr);
