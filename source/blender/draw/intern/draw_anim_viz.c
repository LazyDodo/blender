/*
 * ***** BEGIN GPL LICENSE BLOCK *****
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
 * The Original Code is Copyright (C) 2009/2018 by the Blender Foundation.
 *
 * Contributor(s): Joshua Leung
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/draw/intern/draw_anim_viz.c
 *  \ingroup draw
 */


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "BLI_sys_types.h"

#include "DNA_anim_types.h"
#include "DNA_armature_types.h"
#include "DNA_scene_types.h"
#include "DNA_screen_types.h"
#include "DNA_view3d_types.h"
#include "DNA_object_types.h"

#include "BLI_math.h"
#include "BLI_dlrbTree.h"

#include "BKE_animsys.h"
#include "BKE_action.h"

#include "ED_keyframes_draw.h"

#include "UI_resources.h"

#include "DEG_depsgraph_query.h"

#include "DRW_engine.h"
#include "DRW_render.h"

#include "GPU_shader.h"
#include "GPU_immediate.h"
#include "GPU_matrix.h"

#include "draw_common.h"
#include "draw_manager_text.h"

#include "draw_mode_engines.h"

extern struct GPUUniformBuffer *globals_ubo; /* draw_common.c */

/* ********************************* Lists ************************************** */
/* All lists are per viewport specific datas.
 * They are all free when viewport changes engines
 * or is free itself.
 */

/* XXX: How to show frame numbers, etc.?  Currently only doing the dots and lines */
typedef struct MPATH_PassList {
	struct DRWPass *lines;
	struct DRWPass *points;
} MPATH_PassList;

typedef struct MPATH_StorageList {
	struct MPATH_PrivateData *g_data;
} MPATH_StorageList;

typedef struct MPATH_Data {
	void *engine_type;
	DRWViewportEmptyList *fbl;
	DRWViewportEmptyList *txl;
	MPATH_PassList *psl;
	MPATH_StorageList *stl;
} MPATH_Data;

struct {
	GPUShader *mpath_line_sh;
	GPUShader *mpath_points_sh;
} e_data = {0};

/* *************************** Path Cache *********************************** */

/* Just convert the CPU cache to GPU cache. */
static Gwn_VertBuf *mpath_vbo_get(bMotionPath *mpath)
{
	if (!mpath->points_vbo) {
		Gwn_VertFormat format = {0};
		/* Match structure of bMotionPathVert. */
		uint pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
		GWN_vertformat_attr_add(&format, "flag", GWN_COMP_I32, 1, GWN_FETCH_INT);
		mpath->points_vbo = GWN_vertbuf_create_with_format(&format);
		GWN_vertbuf_data_alloc(mpath->points_vbo, mpath->length);

		/* meh... a useless memcpy. */
		Gwn_VertBufRaw raw_data;
		GWN_vertbuf_attr_get_raw_data(mpath->points_vbo, pos, &raw_data);
		memcpy(GWN_vertbuf_raw_step(&raw_data), mpath->points, sizeof(bMotionPathVert) * mpath->length);
	}
	return mpath->points_vbo;
}

static Gwn_Batch *mpath_batch_line_get(bMotionPath *mpath)
{
	if (!mpath->batch_line) {
		mpath->batch_line = GWN_batch_create(GWN_PRIM_LINE_STRIP, mpath_vbo_get(mpath), NULL);
	}
	return mpath->batch_line;
}

static Gwn_Batch *mpath_batch_points_get(bMotionPath *mpath)
{
	if (!mpath->batch_points) {
		mpath->batch_points = GWN_batch_create(GWN_PRIM_POINTS, mpath_vbo_get(mpath), NULL);
	}
	return mpath->batch_points;
}

/* *************************** Motion Path Drawing ****************************** */

/* TODO:
 * - options to draw paths with lines
 * - include support for editing the path verts */

/* Set up drawing environment for drawing motion paths */
void draw_motion_paths_init(View3D *v3d, ARegion *ar) 
{
	RegionView3D *rv3d = ar->regiondata;
	
	if (v3d->zbuf) glDisable(GL_DEPTH_TEST);
	
	gpuPushMatrix();
	gpuLoadMatrix(rv3d->viewmat);
}

/* set color
 * - more intense for active/selected bones, less intense for unselected bones
 * - black for before current frame, green for current frame, blue for after current frame
 * - intensity decreases as distance from current frame increases
 *
 * If the user select custom color, the color is replaced for the color selected in UI panel
 * - 75% Darker color is used for previous frames
 * - 50% Darker color for current frame
 * - User selected color for next frames
 */
static void set_motion_path_color(Scene *scene, bMotionPath *mpath, int i, short sel, int sfra, int efra,
	float prev_color[3], float frame_color[3], float next_color[3], unsigned color)
{
	int frame = sfra + i;
	int blend_base = (abs(frame - CFRA) == 1) ? TH_CFRAME : TH_BACK; /* "bleed" cframe color to ease color blending */
	unsigned char ubcolor[3];

#define SET_INTENSITY(A, B, C, min, max) (((1.0f - ((C - B) / (C - A))) * (max - min)) + min)
	float intensity;  /* how faint */

	if (frame < CFRA) {
		if (mpath->flag & MOTIONPATH_FLAG_CUSTOM) {
			/* Custom color: previous frames color is darker than current frame */
			rgb_float_to_uchar(ubcolor, prev_color);
		}
		else {
			/* black - before cfra */
			if (sel) {
				/* intensity = 0.5f; */
				intensity = SET_INTENSITY(sfra, i, CFRA, 0.25f, 0.75f);
			}
			else {
				/* intensity = 0.8f; */
				intensity = SET_INTENSITY(sfra, i, CFRA, 0.68f, 0.92f);
			}

			UI_GetThemeColorBlend3ubv(TH_WIRE, blend_base, intensity, ubcolor);
		}
	}
	else if (frame > CFRA) {
		if (mpath->flag & MOTIONPATH_FLAG_CUSTOM) {
			/* Custom color: next frames color is equal to user selected color */
			rgb_float_to_uchar(ubcolor, next_color);
		}
		else {
			/* blue - after cfra */
			if (sel) {
				/* intensity = 0.5f; */
				intensity = SET_INTENSITY(CFRA, i, efra, 0.25f, 0.75f);
			}
			else {
				/* intensity = 0.8f; */
				intensity = SET_INTENSITY(CFRA, i, efra, 0.68f, 0.92f);
			}

			UI_GetThemeColorBlend3ubv(TH_BONE_POSE, blend_base, intensity, ubcolor);
		}
	}
	else {
		if (mpath->flag & MOTIONPATH_FLAG_CUSTOM) {
			/* Custom color: current frame color is slightly darker than user selected color */
			rgb_float_to_uchar(ubcolor, frame_color);
		}
		else {
			/* green - on cfra */
			if (sel) {
				intensity = 0.5f;
			}
			else {
				intensity = 0.99f;
			}
			UI_GetThemeColorBlendShade3ubv(TH_CFRAME, TH_BACK, intensity, 10, ubcolor);
		}
	}

	immAttrib3ubv(color, ubcolor);

#undef SET_INTENSITY
}

/* Draw the given motion path for an Object or a Bone 
 *  - assumes that the viewport has already been initialized properly
 *    i.e. draw_motion_paths_init() has been called
 */
void draw_motion_path_instance(Scene *scene, 
                               Object *ob, bPoseChannel *pchan,
                               bAnimVizSettings *avs, bMotionPath *mpath)
{
	//RegionView3D *rv3d = ar->regiondata;
	bMotionPathVert *mpv, *mpv_start;
	int i, stepsize = avs->path_step;
	int sfra, efra, sind, len;
	float prev_color[3];
	float frame_color[3];
	float next_color[3];

	/* Custom color - Previous frames: color is darker than current frame */
	mul_v3_v3fl(prev_color, mpath->color, 0.25f);

	/* Custom color - Current frame: color is slightly darker than user selected color */
	mul_v3_v3fl(frame_color, mpath->color, 0.5f);

	/* Custom color - Next frames: color is equal to user selection */
	copy_v3_v3(next_color, mpath->color);

	/* Save old line width */
	GLfloat old_width;
	glGetFloatv(GL_LINE_WIDTH, &old_width);
	
	/* get frame ranges */
	if (avs->path_type == MOTIONPATH_TYPE_ACFRA) {
		/* With "Around Current", we only choose frames from around 
		 * the current frame to draw.
		 */
		sfra = CFRA - avs->path_bc;
		efra = CFRA + avs->path_ac;
	}
	else {
		/* Use the current display range */
		sfra = avs->path_sf;
		efra = avs->path_ef;
	}
	
	/* no matter what, we can only show what is in the cache and no more 
	 * - abort if whole range is past ends of path
	 * - otherwise clamp endpoints to extents of path
	 */
	if (sfra < mpath->start_frame) {
		/* start clamp */
		sfra = mpath->start_frame;
	}
	if (efra > mpath->end_frame) {
		/* end clamp */
		efra = mpath->end_frame;
	}
	
	if ((sfra > mpath->end_frame) || (efra < mpath->start_frame)) {
		/* whole path is out of bounds */
		return;
	}
	
	len = efra - sfra;
	
	if ((len <= 0) || (mpath->points == NULL)) {
		return;
	}
	
	/* get pointers to parts of path */
	sind = sfra - mpath->start_frame;
	mpv_start = (mpath->points + sind);

	/* draw curve-line of path */
	/* Draw lines only if line drawing option is enabled */
	if (mpath->flag & MOTIONPATH_FLAG_LINES) {
		/* set line thickness */
		glLineWidth(mpath->line_thickness);

		Gwn_VertFormat *format = immVertexFormat();
		unsigned int pos = GWN_vertformat_attr_add(format, "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
		unsigned int color = GWN_vertformat_attr_add(format, "color", GWN_COMP_U8, 3, GWN_FETCH_INT_TO_FLOAT_UNIT);

		immBindBuiltinProgram(GPU_SHADER_3D_SMOOTH_COLOR);

		immBegin(GWN_PRIM_LINE_STRIP, len);
		for (i = 0, mpv = mpv_start; i < len; i++, mpv++) {
			short sel = (pchan) ? (pchan->bone->flag & BONE_SELECTED) : (ob->flag & SELECT);

			/* Set color */
			set_motion_path_color(scene, mpath, i, sel, sfra, efra, prev_color, frame_color, next_color, color);

			/* draw a vertex with this color */
			immVertex3fv(pos, mpv->co);
		}
		immEnd();

		immUnbindProgram();

		/* back to old line thickness */
		glLineWidth(old_width);
	}

	unsigned int pos = GWN_vertformat_attr_add(immVertexFormat(), "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);

	immBindBuiltinProgram(GPU_SHADER_3D_UNIFORM_COLOR);

	/* Point must be bigger than line thickness */
	glPointSize(mpath->line_thickness + 1.0);
	
	/* draw little black point at each frame */
	immUniformColor3ub(0, 0, 0);

	immBegin(GWN_PRIM_POINTS, len);

	for (i = 0, mpv = mpv_start; i < len; i++, mpv++) {
		immVertex3fv(pos, mpv->co);
	}

	immEnd();

	/* Draw little white dots at each framestep value or replace with custom color */
	if (mpath->flag & MOTIONPATH_FLAG_CUSTOM) {
		immUniformColor3fv(mpath->color);
	}
	else {
		immUniformThemeColor(TH_TEXT_HI);
	}

	immBegin(GWN_PRIM_POINTS, (len + stepsize - 1) / stepsize);

	for (i = 0, mpv = mpv_start; i < len; i += stepsize, mpv += stepsize) {
		immVertex3fv(pos, mpv->co);
	}

	immEnd();
	
	/* Draw big green dot where the current frame is 
	 * NOTE: this is only done when keyframes are shown, since this adds similar types of clutter
	 */
	if ((avs->path_viewflag & MOTIONPATH_VIEW_KFRAS) &&
	    (sfra < CFRA) && (CFRA <= efra)) 
	{
		glPointSize(mpath->line_thickness + 5.0);
		immUniformThemeColor(TH_CFRAME);

		immBegin(GWN_PRIM_POINTS, 1);

		mpv = mpv_start + (CFRA - sfra);
		immVertex3fv(pos, mpv->co);

		immEnd();
	}

	immUnbindProgram();
	
	/* XXX, this isn't up to date but probably should be kept so. */
	invert_m4_m4(ob->imat, ob->obmat);
	
#if 0 /* FIXME!!! How can we draw text in the new 3d Viewport */

	/* Draw frame numbers at each framestep value */
	if (avs->path_viewflag & MOTIONPATH_VIEW_FNUMS) {
		unsigned char col[4];
		UI_GetThemeColor3ubv(TH_TEXT_HI, col);
		col[3] = 255;
		
		for (i = 0, mpv = mpv_start; i < len; i += stepsize, mpv += stepsize) {
			int frame = sfra + i;
			char numstr[32];
			size_t numstr_len;
			float co[3];
			
			/* only draw framenum if several consecutive highlighted points don't occur on same point */
			if (i == 0) {
				numstr_len = sprintf(numstr, " %d", frame);
				mul_v3_m4v3(co, ob->imat, mpv->co);
				view3d_cached_text_draw_add(co, numstr, numstr_len,
				                            0, V3D_CACHE_TEXT_WORLDSPACE | V3D_CACHE_TEXT_ASCII, col);
			}
			else if ((i >= stepsize) && (i < len - stepsize)) {
				bMotionPathVert *mpvP = (mpv - stepsize);
				bMotionPathVert *mpvN = (mpv + stepsize);
				
				if ((equals_v3v3(mpv->co, mpvP->co) == 0) || (equals_v3v3(mpv->co, mpvN->co) == 0)) {
					numstr_len = sprintf(numstr, " %d", frame);
					mul_v3_m4v3(co, ob->imat, mpv->co);
					view3d_cached_text_draw_add(co, numstr, numstr_len,
					                            0, V3D_CACHE_TEXT_WORLDSPACE | V3D_CACHE_TEXT_ASCII, col);
				}
			}
		}
	}
	
	/* Keyframes - dots and numbers */
	if (avs->path_viewflag & MOTIONPATH_VIEW_KFRAS) {
		unsigned char col[4];
		
		AnimData *adt = BKE_animdata_from_id(&ob->id);
		DLRBT_Tree keys;
		
		/* build list of all keyframes in active action for object or pchan */
		BLI_dlrbTree_init(&keys);
		
		if (adt) {
			/* it is assumed that keyframes for bones are all grouped in a single group
			 * unless an option is set to always use the whole action
			 */
			if ((pchan) && (avs->path_viewflag & MOTIONPATH_VIEW_KFACT) == 0) {
				bActionGroup *agrp = BKE_action_group_find_name(adt->action, pchan->name);
				
				if (agrp) {
					agroup_to_keylist(adt, agrp, &keys, NULL);
					BLI_dlrbTree_linkedlist_sync(&keys);
				}
			}
			else {
				action_to_keylist(adt, adt->action, &keys, NULL);
				BLI_dlrbTree_linkedlist_sync(&keys);
			}
		}
		
		/* Draw slightly-larger yellow dots at each keyframe */
		UI_GetThemeColor3ubv(TH_VERTEX_SELECT, col);
		col[3] = 255;
		
		/* point must be bigger than line */
		glPointSize(mpath->line_thickness + 3.0);

		pos = GWN_vertformat_attr_add(immVertexFormat(), "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);

		immBindBuiltinProgram(GPU_SHADER_3D_UNIFORM_COLOR);
		immUniformColor3ubv(col);
		
		immBeginAtMost(GWN_PRIM_POINTS, len);

		for (i = 0, mpv = mpv_start; i < len; i++, mpv++) {
			int    frame = sfra + i; 
			float mframe = (float)(frame);
			
			if (BLI_dlrbTree_search_exact(&keys, compare_ak_cfraPtr, &mframe)) {
				immVertex3fv(pos, mpv->co);
			}
		}

		immEnd();

		immUnbindProgram();
		
		/* Draw frame numbers of keyframes  */
		if (avs->path_viewflag & MOTIONPATH_VIEW_KFNOS) {
			float co[3];
			for (i = 0, mpv = mpv_start; i < len; i++, mpv++) {
				float mframe = (float)(sfra + i);
				
				if (BLI_dlrbTree_search_exact(&keys, compare_ak_cfraPtr, &mframe)) {
					char numstr[32];
					size_t numstr_len;
					
					numstr_len = sprintf(numstr, " %d", (sfra + i));
					mul_v3_m4v3(co, ob->imat, mpv->co);
					view3d_cached_text_draw_add(co, numstr, numstr_len,
					                            0, V3D_CACHE_TEXT_WORLDSPACE | V3D_CACHE_TEXT_ASCII, col);
				}
			}
		}
		
		BLI_dlrbTree_free(&keys);
	}
#endif
}

/* Clean up drawing environment after drawing motion paths */
void draw_motion_paths_cleanup(View3D *v3d)
{
	if (v3d->zbuf) glEnable(GL_DEPTH_TEST);
	gpuPopMatrix();
}

/* ******************* Populate Cache from MotionPath Sources ******************* */

/* Motion Paths from Object */


/* *************************** Draw Engine Entrypoints ************************** */

static void MPATH_engine_init(void *UNUSED(vedata))
{
}

static void MPATH_engine_free(void)
{
}

/* Here init all passes and shading groups
 * Assume that all Passes are NULL */
static void MPATH_cache_init(void *vedata)
{
	MPATH_PassList *psl = ((MPATH_Data *)vedata)->psl;

	{
		DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS;
		psl->lines = DRW_pass_create("Motionpath Line Pass", state);
	}

	{
		DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS | DRW_STATE_POINT;
		psl->points = DRW_pass_create("Motionpath Point Pass", state);
	}
}

static void MPATH_cache_motion_path(MPATH_PassList *psl, Scene *scene,
                                    Object *ob, bPoseChannel *pchan,
                                    bAnimVizSettings *avs, bMotionPath *mpath)
{
	const DRWContextState *draw_ctx = DRW_context_state_get();
	struct DRWTextStore *dt = DRW_text_cache_ensure();
	int txt_flag = DRW_TEXT_CACHE_GLOBALSPACE | DRW_TEXT_CACHE_ASCII;
	int stepsize = avs->path_step;
	int sfra, efra, sind, len;
	int cfra = (int)DEG_get_ctime(draw_ctx->depsgraph);
	bool sel = (pchan) ? (pchan->bone->flag & BONE_SELECTED) : (ob->flag & SELECT);
	bool show_keyframes = (avs->path_viewflag & MOTIONPATH_VIEW_KFRAS) != 0;
	bMotionPathVert *mpv, *mpv_start;

	/* get frame ranges */
	if (avs->path_type == MOTIONPATH_TYPE_ACFRA) {
		/* With "Around Current", we only choose frames from around
		 * the current frame to draw.
		 */
		sfra = cfra - avs->path_bc;
		efra = cfra + avs->path_ac + 1;
	}
	else {
		/* Use the current display range */
		sfra = avs->path_sf;
		efra = avs->path_ef;
	}

	/* no matter what, we can only show what is in the cache and no more
	 * - abort if whole range is past ends of path
	 * - otherwise clamp endpoints to extents of path
	 */
	if (sfra < mpath->start_frame) {
		/* start clamp */
		sfra = mpath->start_frame;
	}
	if (efra > mpath->end_frame) {
		/* end clamp */
		efra = mpath->end_frame;
	}

	if ((sfra > mpath->end_frame) || (efra < mpath->start_frame)) {
		/* whole path is out of bounds */
		return;
	}

	len = efra - sfra;

	if ((len <= 0) || (mpath->points == NULL)) {
		return;
	}

	sind = sfra - mpath->start_frame;
	mpv_start = (mpath->points + sind);

	bool use_custom_col = (mpath->flag & MOTIONPATH_FLAG_CUSTOM) != 0;

	/* draw curve-line of path */
	/* Draw lines only if line drawing option is enabled */
	if (mpath->flag & MOTIONPATH_FLAG_LINES) {
		DRWShadingGroup *shgrp = DRW_shgroup_create(mpath_line_shader_get(), psl->lines);
		DRW_shgroup_uniform_int_copy(shgrp, "frameCurrent", cfra);
		DRW_shgroup_uniform_int_copy(shgrp, "frameStart", sfra);
		DRW_shgroup_uniform_int_copy(shgrp, "frameEnd", efra);
		DRW_shgroup_uniform_int_copy(shgrp, "cacheStart", mpath->start_frame);
		DRW_shgroup_uniform_int_copy(shgrp, "lineThickness", mpath->line_thickness);
		DRW_shgroup_uniform_bool_copy(shgrp, "selected", sel);
		DRW_shgroup_uniform_bool_copy(shgrp, "useCustomColor", use_custom_col);
		DRW_shgroup_uniform_vec2(shgrp, "viewportSize", DRW_viewport_size_get(), 1);
		DRW_shgroup_uniform_block(shgrp, "globalsBlock", globals_ubo);
		if (use_custom_col) {
			DRW_shgroup_uniform_vec3(shgrp, "customColor", mpath->color, 1);
		}
		/* Only draw the required range. */
		DRW_shgroup_call_range_add(shgrp, mpath_batch_line_get(mpath), NULL, sind, len);
	}

	/* Draw points. */
	DRWShadingGroup *shgrp = DRW_shgroup_create(mpath_points_shader_get(), psl->points);
	DRW_shgroup_uniform_int_copy(shgrp, "frameCurrent", cfra);
	DRW_shgroup_uniform_int_copy(shgrp, "cacheStart", mpath->start_frame);
	DRW_shgroup_uniform_int_copy(shgrp, "pointSize", mpath->line_thickness);
	DRW_shgroup_uniform_int_copy(shgrp, "stepSize", stepsize);
	DRW_shgroup_uniform_bool_copy(shgrp, "selected", sel);
	DRW_shgroup_uniform_bool_copy(shgrp, "showKeyFrames", show_keyframes);
	DRW_shgroup_uniform_bool_copy(shgrp, "useCustomColor", use_custom_col);
	DRW_shgroup_uniform_block(shgrp, "globalsBlock", globals_ubo);
	if (use_custom_col) {
		DRW_shgroup_uniform_vec3(shgrp, "customColor", mpath->color, 1);
	}
	/* Only draw the required range. */
	DRW_shgroup_call_range_add(shgrp, mpath_batch_points_get(mpath), NULL, sind, len);

	/* Draw frame numbers at each framestep value */
	bool show_kf_no = (avs->path_viewflag & MOTIONPATH_VIEW_KFNOS) != 0;
	if ((avs->path_viewflag & (MOTIONPATH_VIEW_FNUMS)) || (show_kf_no && show_keyframes)) {
		int i;
		unsigned char col[4], col_kf[4];
		UI_GetThemeColor3ubv(TH_TEXT_HI, col);
		UI_GetThemeColor3ubv(TH_VERTEX_SELECT, col_kf);
		col[3] = col_kf[3] = 255;

		for (i = 0, mpv = mpv_start; i < len; i += stepsize, mpv += stepsize) {
			int frame = sfra + i;
			char numstr[32];
			size_t numstr_len;
			float co[3];
			bool is_keyframe = (mpv->flag & MOTIONPATH_VERT_KEY) != 0;

			if ((show_keyframes && show_kf_no && is_keyframe) ||
			    ((avs->path_viewflag & MOTIONPATH_VIEW_FNUMS) && (i == 0)))
			{
				numstr_len = sprintf(numstr, " %d", frame);
				mul_v3_m4v3(co, ob->imat, mpv->co);
				DRW_text_cache_add(dt, co, numstr, numstr_len, 0, txt_flag, (is_keyframe) ? col_kf : col);
			}
			else if (avs->path_viewflag & MOTIONPATH_VIEW_FNUMS) {
				bMotionPathVert *mpvP = (mpv - stepsize);
				bMotionPathVert *mpvN = (mpv + stepsize);
				/* only draw framenum if several consecutive highlighted points don't occur on same point */
				if ((equals_v3v3(mpv->co, mpvP->co) == 0) || (equals_v3v3(mpv->co, mpvN->co) == 0)) {
					numstr_len = sprintf(numstr, " %d", frame);
					mul_v3_m4v3(co, ob->imat, mpv->co);
					DRW_text_cache_add(dt, co, numstr, numstr_len, 0, txt_flag, col);
				}
			}
		}
	}
}

/* Add geometry to shading groups. Execute for each objects */
static void MPATH_cache_populate(void *vedata, Object *ob)
{
	MPATH_PassList *psl = ((MPATH_Data *)vedata)->psl;
	const DRWContextState *draw_ctx = DRW_context_state_get();

	if ((draw_ctx->v3d->overlay.flag & V3D_OVERLAY_MOTION_PATHS) == 0) {
		return;
	}

	if (ob->type == OB_ARMATURE) {
		if (DRW_pose_mode_armature(ob, draw_ctx->obact)) {
			for (bPoseChannel *pchan = ob->pose->chanbase.first; pchan; pchan = pchan->next) {
				if (pchan->mpath) {
					MPATH_cache_motion_path(psl, scene_eval, ob, pchan, &ob->pose->avs, pchan->mpath);
				}
			}
		}
	}
	else {
		if (ob->mpath) {
			MPATH_cache_motion_path(psl, scene_eval, ob, NULL, &ob->avs, ob->mpath);
		}
	}
}

/* Draw time! Control rendering pipeline from here */
static void MPATH_draw_scene(void *vedata)
{
	MPATH_PassList *psl = ((MPATH_Data *)vedata)->psl;
	DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();
	DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
	
	MULTISAMPLE_SYNC_ENABLE(dfbl, dtxl)

	DRW_draw_pass(psl->lines);
	DRW_draw_pass(psl->points);

	MULTISAMPLE_SYNC_DISABLE(dfbl, dtxl)
}

/* *************************** Draw Engine Defines ****************************** */

static const DrawEngineDataSize MPATH_data_size = DRW_VIEWPORT_DATA_SIZE(MPATH_Data);

DrawEngineType draw_engine_motion_path_type = {
	NULL, NULL,
	N_("MotionPath"),
	&MPATH_data_size,
	&MPATH_engine_init,
	&MPATH_engine_free,
	&MPATH_cache_init,
	&MPATH_cache_populate,
	NULL,
	NULL,
	&MPATH_draw_scene,
	NULL,
	NULL,
};

