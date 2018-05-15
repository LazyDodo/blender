/*
 * Copyright 2017, Blender Foundation.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor(s): Antonio Vazquez
 *
 */

/** \file blender/draw/engines/gpencil/gpencil_engine.h
 *  \ingroup draw
 */

#ifndef __GPENCIL_ENGINE_H__
#define __GPENCIL_ENGINE_H__

#include "GPU_batch.h"

struct tGPspoint;
struct bGPDstroke;
struct ModifierData;
struct GPENCIL_Data;
struct GPENCIL_StorageList;
struct Object;
struct MaterialGPencilStyle;
struct RenderEngine;
struct RenderLayer;

 /* TODO: these could be system parameter in userprefs screen */
#define GPENCIL_MAX_GP_OBJ 256 

#define GPENCIL_CACHE_BLOCK_SIZE 8 
#define GPENCIL_MAX_SHGROUPS 65536
#define GPENCIL_MIN_BATCH_SLOTS_CHUNK 16

#define GPENCIL_COLOR_SOLID   0
#define GPENCIL_COLOR_TEXTURE 1
#define GPENCIL_COLOR_PATTERN 2

#define GP_SIMPLIFY(ts) ((ts->gpencil_simplify & GP_TOOL_FLAG_SIMPLIFY))
#define GP_SIMPLIFY_ONPLAY(playing) (((playing == true) && (ts->gpencil_simplify & GP_TOOL_FLAG_SIMPLIFY_ON_PLAY)) || ((ts->gpencil_simplify & GP_TOOL_FLAG_SIMPLIFY_ON_PLAY) == 0))
#define GP_SIMPLIFY_FILL(ts, playing) ((GP_SIMPLIFY_ONPLAY(playing) && (GP_SIMPLIFY(ts)) && (ts->gpencil_simplify & GP_TOOL_FLAG_SIMPLIFY_VIEW_FILL))) 
#define GP_SIMPLIFY_MODIF(ts, playing) ((GP_SIMPLIFY_ONPLAY(playing) && (GP_SIMPLIFY(ts)) && (ts->gpencil_simplify & GP_TOOL_FLAG_SIMPLIFY_VIEW_MODIF))) 

#define GP_IS_CAMERAVIEW ((rv3d != NULL) && (rv3d->persp == RV3D_CAMOB && v3d->camera))

 /* *********** OBJECTS CACHE *********** */

/* used to save gpencil objects */
typedef struct tGPencilObjectCache {
	struct Object *ob;
	int init_grp, end_grp;
	int idx;  /*original index, can change after sort */

	float zdepth;  /* z-depth value to sort gp object */
	bool temp_ob;  /* flag to tag temporary objects that must be removed after drawing loop */
} tGPencilObjectCache;

  /* *********** LISTS *********** */

typedef struct GPENCIL_shgroup {
	int s_clamp;
	int stroke_style;
	int color_type;
	int mode;
	int t_mix;
	int t_flip;
	int t_clamp;
	int fill_style;
	int keep_size;
	float obj_scale;
	struct DRWShadingGroup *shgrps_fill;
	struct DRWShadingGroup *shgrps_stroke;
} GPENCIL_shgroup;

typedef struct GPENCIL_Storage {
	int shgroup_id; /* total elements */
	float unit_matrix[4][4];
	int stroke_style;
	int color_type;
	int mode;
	int xray;
	int keep_size;
	float obj_scale;
	float pixfactor;
	int playing;
	bool is_render;
	int uselines;
	float gridsize[2];
	float gridcolor[3];
	const float *pixsize;
	float render_pixsize;
	int tonemapping;

	/* simplify settings*/
	bool simplify_fill;
	bool simplify_modif;
	bool simplify_vfx;

	/* Render Matrices and data */
	float persmat[4][4], persinv[4][4];
	float viewmat[4][4], viewinv[4][4];
	float winmat[4][4], wininv[4][4];
	float view_vecs[2][4]; /* vec4[2] */
} GPENCIL_Storage;

typedef struct GPENCIL_StorageList {
	struct GPENCIL_Storage *storage;
	struct g_data *g_data;
	struct GPENCIL_shgroup *shgroups;
} GPENCIL_StorageList;

typedef struct GPENCIL_PassList {
	struct DRWPass *stroke_pass;
	struct DRWPass *edit_pass;
	struct DRWPass *drawing_pass;
	struct DRWPass *mix_pass;
	struct DRWPass *mix_pass_noblend;
	struct DRWPass *painting_pass;
	struct DRWPass *paper_pass;
} GPENCIL_PassList;

typedef struct GPENCIL_FramebufferList {
	struct GPUFrameBuffer *main;
	struct GPUFrameBuffer *temp_fb_a;
	struct GPUFrameBuffer *painting_fb;
} GPENCIL_FramebufferList;

typedef struct GPENCIL_TextureList {
	struct GPUTexture *texture;
} GPENCIL_TextureList;

typedef struct GPENCIL_Data {
	void *engine_type; /* Required */
	struct GPENCIL_FramebufferList *fbl;
	struct GPENCIL_TextureList *txl;
	struct GPENCIL_PassList *psl;
	struct GPENCIL_StorageList *stl;

	/* render textures */
	struct GPUTexture *render_depth_tx;
	struct GPUTexture *render_color_tx;
} GPENCIL_Data;

/* *********** STATIC *********** */
typedef struct g_data {
	struct DRWShadingGroup *shgrps_edit_point;
	struct DRWShadingGroup *shgrps_edit_line;
	struct DRWShadingGroup *shgrps_drawing_stroke;
	struct DRWShadingGroup *shgrps_drawing_fill;

	/* for buffer only one batch is nedeed because the drawing is only of one stroke */
	Gwn_Batch *batch_buffer_stroke;
	Gwn_Batch *batch_buffer_fill;

	int gp_cache_used; /* total objects in cache */
	int gp_cache_size; /* size of the cache */
	struct tGPencilObjectCache *gp_object_cache;

	int session_flag;

} g_data; /* Transient data */

/* flags for fast drawing support */
typedef enum eGPsession_Flag {
	GP_DRW_PAINT_HOLD     = (1 << 0),
	GP_DRW_PAINT_IDLE     = (1 << 1),
	GP_DRW_PAINT_FILLING  = (1 << 2),
	GP_DRW_PAINT_READY    = (1 << 3),
	GP_DRW_PAINT_PAINTING = (1 << 4),
} eGPsession_Flag;

typedef struct GPENCIL_e_data {
	/* general drawing shaders */
	struct GPUShader *gpencil_fill_sh;
	struct GPUShader *gpencil_stroke_sh;
	struct GPUShader *gpencil_point_sh;
	struct GPUShader *gpencil_edit_point_sh;
	struct GPUShader *gpencil_line_sh;
	struct GPUShader *gpencil_drawing_fill_sh;
	struct GPUShader *gpencil_fullscreen_sh;
	struct GPUShader *gpencil_simple_fullscreen_sh;
	struct GPUShader *gpencil_painting_sh;
	struct GPUShader *gpencil_paper_sh;

	/* textures */
	struct GPUTexture *painting_depth_tx;
	struct GPUTexture *painting_color_tx;

	struct GPUTexture *gpencil_blank_texture;
	
	/* runtime pointers texture */
	struct GPUTexture *input_depth_tx;
	struct GPUTexture *input_color_tx;

	/* working textures */
	struct GPUTexture *temp_color_tx_a;
	struct GPUTexture *temp_depth_tx_a;

} GPENCIL_e_data; /* Engine data */

/* Gwn_Batch Cache */
typedef struct GpencilBatchCache {
	/* For normal strokes, a variable number of batch can be needed depending of number of strokes.
	   It could use the stroke number as total size, but when activate the onion skining, the number
	   can change, so the size is changed dinamically.
	 */
	Gwn_Batch **batch_stroke;
	Gwn_Batch **batch_fill;
	Gwn_Batch **batch_edit;
	Gwn_Batch **batch_edlin;

	/* settings to determine if cache is invalid */
	bool is_dirty;
	bool is_editmode;
	int cache_frame;

	/* keep information about the size of the cache */
	int cache_size;  /* total batch slots available */
	int cache_idx;   /* current slot index */
} GpencilBatchCache;

/* general drawing functions */
struct DRWShadingGroup *DRW_gpencil_shgroup_stroke_create(struct GPENCIL_e_data *e_data, struct GPENCIL_Data *vedata, struct DRWPass *pass, struct GPUShader *shader, 
														  struct Object *ob, struct bGPdata *gpd, struct MaterialGPencilStyle *gp_style, int id, bool onion);
void DRW_gpencil_populate_datablock(struct GPENCIL_e_data *e_data, void *vedata, struct Scene *scene, struct Object *ob, struct ToolSettings *ts, struct bGPdata *gpd);
void DRW_gpencil_populate_buffer_strokes(struct GPENCIL_e_data *e_data, void *vedata, struct ToolSettings *ts, struct Object *ob);
void DRW_gpencil_populate_multiedit(struct GPENCIL_e_data *e_data, void *vedata, struct Scene *scene, struct Object *ob, struct bGPdata *gpd);
void DRW_gpencil_triangulate_stroke_fill(struct bGPDstroke *gps);

/* create geometry functions */
struct Gwn_Batch *DRW_gpencil_get_point_geom(struct bGPDstroke *gps, short thickness, const float ink[4]);
struct Gwn_Batch *DRW_gpencil_get_stroke_geom(struct bGPDframe *gpf, struct bGPDstroke *gps, short thickness, const float ink[4]);
struct Gwn_Batch *DRW_gpencil_get_fill_geom(struct Object *ob, struct bGPDstroke *gps, const float color[4]);
struct Gwn_Batch *DRW_gpencil_get_edit_geom(struct bGPDstroke *gps, float alpha, short dflag);
struct Gwn_Batch *DRW_gpencil_get_edlin_geom(struct bGPDstroke *gps, float alpha, short dflag);
struct Gwn_Batch *DRW_gpencil_get_buffer_stroke_geom(struct bGPdata *gpd, float matrix[4][4], short thickness);
struct Gwn_Batch *DRW_gpencil_get_buffer_fill_geom(struct bGPdata *gpd);
struct Gwn_Batch *DRW_gpencil_get_buffer_point_geom(struct bGPdata *gpd, float matrix[4][4], short thickness);

/* object cache functions */
struct tGPencilObjectCache *gpencil_object_cache_add(struct tGPencilObjectCache *cache_array, struct Object *ob,
							  bool is_temp, int *gp_cache_size, int *gp_cache_used);

/* geometry batch cache functions */
void gpencil_batch_cache_check_free_slots(struct Object *ob);
struct GpencilBatchCache *gpencil_batch_cache_get(struct Object *ob, int cfra);

/* modifier functions */
void gpencil_instance_modifiers(struct GPENCIL_StorageList *stl, struct Object *ob);

/* main functions */
void GPENCIL_engine_init(void *vedata);
void GPENCIL_cache_init(void *vedata);
void GPENCIL_cache_populate(void *vedata, struct Object *ob);
void GPENCIL_cache_finish(void *vedata);
void GPENCIL_draw_scene(void *vedata);

/* render */
void GPENCIL_render_init(struct GPENCIL_Data *ved, struct RenderEngine *engine, struct Depsgraph *depsgraph);
void GPENCIL_render_to_image(void *vedata, struct RenderEngine *engine, struct RenderLayer *render_layer, const rcti *rect);

#endif /* __GPENCIL_ENGINE_H__ */
