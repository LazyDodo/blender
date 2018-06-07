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
	struct DRWPass *edge_thinning_2;
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
	
	struct GPUTexture *dpix_out_pl;
	struct GPUTexture *dpix_out_pr;
	struct GPUTexture *dpix_out_length;

} LANPR_TextureList;

typedef struct LANPR_PrivateData {
	DRWShadingGroup *multipass_shgrp;
	DRWShadingGroup *edge_detect_shgrp;
	DRWShadingGroup *edge_thinning_shgrp;
	DRWShadingGroup *edge_thinning_shgrp_2;
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
	
	ListBase      pending_samples;
	ListBase      erased_samples;
    ListBase      line_strips;

	// dpix data

	void*         atlas_pl;
	void*         atlas_pr;
	void*         atlas_nl;
	void*         atlas_nr;	

	int           begin_index;

	int           dpix_sample_step;
	int           dpix_is_perspective;
	float         dpix_viewport[4];
    int           dpix_buffer_width;
	float         dpix_depth_offset;

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



// functions 

// dpix

void lanpr_init_atlas_inputs(void *ved);
void lanpr_destroy_atlas(void *ved);
int lanpr_feed_atlas_data_obj(void* vedata,
	float* AtlasPointsL, float* AtlasPointsR,
	float* AtlasFaceNormalL, float* AtlasFaceNormalR,
	Object* ob, int BeginIndex);
void lanpr_feed_atlas_trigger_preview_obj(void* vedata, Object* ob, int BeginIndex);

void lanpr_dpix_draw_scene(LANPR_TextureList* txl, LANPR_FramebufferList * fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr);


//snake

void lanpr_snake_draw_scene(LANPR_TextureList* txl, LANPR_FramebufferList * fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr);
