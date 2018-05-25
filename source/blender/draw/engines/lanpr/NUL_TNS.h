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
	struct DRWPass *depth_pass;
	struct DRWPass *color_pass;
	struct DRWPass *normal_pass;
	struct DRWPass *edge_intermediate;
	struct DRWPass *edge_thinning;
	struct DRWPass *edge_thinning_2;
	struct DRWPass *snake_pass;
} LANPR_PassList;

typedef struct LANPR_FramebufferList {
	struct GPUFrameBuffer *passes;
	struct GPUFrameBuffer *edge_intermediate;
	struct GPUFrameBuffer *edge_thinning;
	struct GPUFrameBuffer *on_screen;
	//and something...
} LANPR_FramebufferList;

typedef struct LANPR_TextureList {
	struct GPUTexture *color;
    struct GPUTexture *normal;
	struct GPUTexture *depth;
	struct GPUTexture *edge_intermediate;
} LANPR_TextureList;

typedef struct LANPR_PrivateData {
	DRWShadingGroup *multipass_shgrp;
	DRWShadingGroup *edge_detect_shgrp;
	DRWShadingGroup *edge_thinning_shgrp;
	DRWShadingGroup *edge_thinning_shgrp_2;
    DRWShadingGroup *snake_shgrp;
	
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
	
	ListBase      pending_samples;
	ListBase      erased_samples;
    ListBase      line_strips;

	// drawing

	unsigned        v_buf;
	unsigned        i_buf;
	unsigned        l_buf;
    
	Gwn_VertFormat   snake_gwn_format;
	Gwn_Batch*       snake_batch;

} LANPR_PrivateData;

typedef struct LANPR_StorageList {
	LANPR_PrivateData *g_data;
} LANPR_StorageList;

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


extern RenderEngineType DRW_engine_viewport_lanpr_type;




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


void tnsClearAll();
void tnsClearColorv(real* rgba);
void tnsClearColor(real r, real g, real b, real a);
void tnsSwitchToCurrentWindowContext(void* wnd);