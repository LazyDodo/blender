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
 * The Original Code is Copyright (C) 2017, Blender Foundation
 * This is a new part of Blender
 *
 * Contributor(s): Antonio Vazquez
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 * Operators for creating new Grease Pencil primitives (boxes, circles, ...)
 */

 /** \file blender/editors/gpencil/gpencil_primitive.c
  *  \ingroup edgpencil
  */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_utildefines.h"
#include "BLI_math.h"

#include "BLT_translation.h"

#include "DNA_brush_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"
#include "DNA_screen_types.h"
#include "DNA_space_types.h"
#include "DNA_view3d_types.h"

#include "BKE_brush.h"
#include "BKE_colortools.h"
#include "BKE_context.h"
#include "BKE_deform.h"
#include "BKE_global.h"
#include "BKE_gpencil.h"
#include "BKE_main.h"
#include "BKE_material.h"
#include "BKE_paint.h"
#include "BKE_report.h"

#include "UI_interface.h"

#include "WM_api.h"
#include "WM_types.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "ED_gpencil.h"
#include "ED_object.h"
#include "ED_screen.h"
#include "ED_view3d.h"
#include "ED_space_api.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "gpencil_intern.h"

#define MIN_EDGES 2
#define MAX_EDGES 128

#define IDLE 0
#define IN_PROGRESS 1
#define IN_CURVE_EDIT_BEZIER 2

#define CURVE_CP_NONE 0
#define CURVE_CP_T 1
#define CURVE_CP_1 2
#define CURVE_CP_2 3
#define CURVE_CP_B 4

  /* ************************************************ */
  /* Core/Shared Utilities */

/* clear the session buffers (call this before AND after a paint operation) */
static void gp_session_validatebuffer(tGPDprimitive *p)
{
	bGPdata *gpd = p->gpd;

	/* clear memory of buffer (or allocate it if starting a new session) */
	if (gpd->runtime.sbuffer) {
		memset(gpd->runtime.sbuffer, 0, sizeof(tGPspoint) * GP_STROKE_BUFFER_MAX);
	}
	else {
		gpd->runtime.sbuffer = MEM_callocN(sizeof(tGPspoint) * GP_STROKE_BUFFER_MAX, "gp_session_strokebuffer");
	}

	/* reset indices */
	gpd->runtime.sbuffer_size = 0;

	/* reset flags */
	gpd->runtime.sbuffer_sflag = 0;
	gpd->runtime.sbuffer_sflag |= GP_STROKE_3DSPACE;
	if (p->cyclic) {
		gpd->runtime.sbuffer_sflag |= GP_STROKE_CYCLIC;
	}
}

static void gp_init_colors(tGPDprimitive *p)
{
	bGPdata *gpd = p->gpd;
	Brush *brush = p->brush;

	Material *ma = NULL;
	MaterialGPencilStyle *gp_style = NULL;

	/* use brush material */
	ma = BKE_gpencil_get_material_from_brush(brush);

	/* if no brush defaults, get material and color info
	 * NOTE: Ensures that everything we need will exist...
	 */
	if ((ma == NULL) || (ma->gp_style == NULL)) {
		BKE_gpencil_material_ensure(p->bmain, p->ob);

		/* assign always the first material to the brush */
		p->mat = give_current_material(p->ob, 1);
		brush->gpencil_settings->material = p->mat;
	}
	else {
		p->mat = ma;
	}

	/* check if the material is already on object material slots and add it if missing */
	if (BKE_gpencil_get_material_index(p->ob, p->mat) == 0) {
		BKE_object_material_slot_add(p->bmain, p->ob);
		assign_material(p->bmain, p->ob, ma, p->ob->totcol, BKE_MAT_ASSIGN_USERPREF);
	}

	/* assign color information to temp data */
	gp_style = p->mat->gp_style;
	if (gp_style) {

		/* set colors */
		copy_v4_v4(gpd->runtime.scolor, gp_style->stroke_rgba);
		copy_v4_v4(gpd->runtime.sfill, gp_style->fill_rgba);
		/* add some alpha to make easy the filling without hide strokes */
		if (gpd->runtime.sfill[3] > 0.8f) {
			gpd->runtime.sfill[3] = 0.8f;
		}

		gpd->runtime.mode = (short)gp_style->mode;
		gpd->runtime.bstroke_style = gp_style->stroke_style;
		gpd->runtime.bfill_style = gp_style->fill_style;
	}
}

  /* Poll callback for primitive operators */
static bool gpencil_primitive_add_poll(bContext *C)
{
	/* only 3D view */
	ScrArea *sa = CTX_wm_area(C);
	if (sa && sa->spacetype != SPACE_VIEW3D) {
		return 0;
	}

	/* need data to create primitive */
	bGPdata *gpd = CTX_data_gpencil_data(C);
	if (gpd == NULL) {
		return 0;
	}

	/* only in edit and paint modes
	 * - paint as it's the "drawing/creation mode"
	 * - edit as this is more of an atomic editing operation
	 *   (similar to copy/paste), and also for consistency
	 */
	if ((gpd->flag & (GP_DATA_STROKE_PAINTMODE | GP_DATA_STROKE_EDITMODE)) == 0) {
		CTX_wm_operator_poll_msg_set(C, "Primitives can only be added in Draw or Edit modes");
		return 0;
	}

	/* don't allow operator to function if the active layer is locked/hidden
	 * (BUT, if there isn't an active layer, we are free to add new layer when the time comes)
	 */
	bGPDlayer *gpl = BKE_gpencil_layer_getactive(gpd);
	if ((gpl) && (gpl->flag & (GP_LAYER_LOCKED | GP_LAYER_HIDE))) {
		CTX_wm_operator_poll_msg_set(C, "Primitives cannot be added as active layer is locked or hidden");
		return 0;
	}

	return 1;
}


/* ****************** Primitive Interactive *********************** */

/* Helper: Create internal strokes primitives data */
static void gp_primitive_set_initdata(bContext *C, tGPDprimitive *tgpi)
{
	ToolSettings *ts = CTX_data_tool_settings(C);
	Depsgraph *depsgraph = CTX_data_depsgraph(C);
	int cfra_eval = (int)DEG_get_ctime(depsgraph);

	bGPDlayer *gpl = CTX_data_active_gpencil_layer(C);

	/* if brush doesn't exist, create a new one */
	Paint *paint = &ts->gp_paint->paint;
	/* if not exist, create a new one */
	if (paint->brush == NULL) {
		/* create new brushes */
		BKE_brush_gpencil_presets(C);
	}
	tgpi->brush = paint->brush;

	/* if layer doesn't exist, create a new one */
	if (gpl == NULL) {
		gpl = BKE_gpencil_layer_addnew(tgpi->gpd, DATA_("Primitives"), true);
	}
	tgpi->gpl = gpl;

	/* create a new temporary frame */
	tgpi->gpf = MEM_callocN(sizeof(bGPDframe), "Temp bGPDframe");
	tgpi->gpf->framenum = tgpi->cframe = cfra_eval;

	/* create new temp stroke */
	bGPDstroke *gps = MEM_callocN(sizeof(bGPDstroke), "Temp bGPDstroke");
	gps->thickness = 2.0f;
	gps->inittime = 0.0f;

	/* enable recalculation flag by default */
	gps->flag |= GP_STROKE_RECALC_CACHES;
	gps->flag &= ~GP_STROKE_SELECT;
	/* the polygon must be closed, so enabled cyclic */
	if (tgpi->type != GP_STROKE_LINE && tgpi->type != GP_STROKE_ARC) {
		gps->flag |= GP_STROKE_CYCLIC;
	}
	else {
		gps->flag &= ~GP_STROKE_CYCLIC;
	}
	gps->flag |= GP_STROKE_3DSPACE;

	gps->mat_nr = BKE_gpencil_get_material_index(tgpi->ob, tgpi->mat) - 1;

	/* allocate memory for storage points, but keep empty */
	gps->totpoints = 0;
	gps->points = MEM_callocN(sizeof(bGPDspoint), "gp_stroke_points");
	/* initialize triangle memory to dummy data */
	gps->tot_triangles = 0;
	gps->triangles = NULL;
	gps->flag |= GP_STROKE_RECALC_CACHES;

	/* add to strokes */
	BLI_addtail(&tgpi->gpf->strokes, gps);
}

/* ----------------------- */
/* Drawing Callbacks */

/* Drawing callback for modal operator in 3d mode */
static void gpencil_primitive_draw_3d(const bContext *C, ARegion *UNUSED(ar), void *arg)
{
	tGPDprimitive *tgpi = (tGPDprimitive *)arg;
	ED_gp_draw_primitives(C, tgpi, REGION_DRAW_POST_VIEW);
}

/* ----------------------- */

/* Helper: Draw status message while the user is running the operator */
static void gpencil_primitive_status_indicators(bContext *C, tGPDprimitive *tgpi)
{
	Scene *scene = tgpi->scene;
	char status_str[UI_MAX_DRAW_STR];
	char msg_str[UI_MAX_DRAW_STR];

	if (tgpi->type == GP_STROKE_BOX) {
		BLI_strncpy(msg_str, IFACE_("Rectangle: ESC/RMB to cancel, LMB set origin, Enter/LMB to confirm, Shift to square, Alt to center"), UI_MAX_DRAW_STR);
	}
	else if (tgpi->type == GP_STROKE_LINE) {
		BLI_strncpy(msg_str, IFACE_("Line: ESC/RMB to cancel, LMB set origin, Enter/LMB to confirm, Alt to center"), UI_MAX_DRAW_STR);
	}
	else if (tgpi->type == GP_STROKE_ARC) {
		BLI_strncpy(msg_str, IFACE_("Arc: ESC/RMB to cancel, Enter/LMB to confirm, WHEEL/+- to adjust edge number, Shift to square, Alt to center, F to flip, C to Close"), UI_MAX_DRAW_STR);
	}
	else if (tgpi->type == GP_STROKE_BEZIER) {
		BLI_strncpy(msg_str, IFACE_("Bezier: ESC/RMB to cancel, Enter/LMB to confirm, WHEEL/+- to adjust edge number, Shift to square, Alt to center, C to Close"), UI_MAX_DRAW_STR);
	}
	else {
		BLI_strncpy(msg_str, IFACE_("Circle: ESC/RMB to cancel, Enter/LMB to confirm, WHEEL/+- to adjust edge number, Shift to square, Alt to center"), UI_MAX_DRAW_STR);
	}

	if (tgpi->type == GP_STROKE_CIRCLE || tgpi->type == GP_STROKE_ARC) {
		if (hasNumInput(&tgpi->num)) {
			char str_offs[NUM_STR_REP_LEN];

			outputNumInput(&tgpi->num, str_offs, &scene->unit);
			BLI_snprintf(status_str, sizeof(status_str), "%s: %s", msg_str, str_offs);
		}
		else {
			if (tgpi->flag == IN_PROGRESS) {
				BLI_snprintf(
					status_str, sizeof(status_str), "%s: %d (%d, %d) (%d, %d)", msg_str, (int)tgpi->tot_edges,
					tgpi->top[0], tgpi->top[1], tgpi->bottom[0], tgpi->bottom[1]);
			}
			else {
				BLI_snprintf(
					status_str, sizeof(status_str), "%s: %d (%d, %d)", msg_str, (int)tgpi->tot_edges,
					tgpi->bottom[0], tgpi->bottom[1]);
			}
		}
	}
	else {
		if (tgpi->flag == IN_PROGRESS) {
			BLI_snprintf(
				status_str, sizeof(status_str), "%s: (%d, %d) (%d, %d)", msg_str,
				tgpi->top[0], tgpi->top[1], tgpi->bottom[0], tgpi->bottom[1]);
		}
		else {
			BLI_snprintf(
				status_str, sizeof(status_str), "%s: (%d, %d)", msg_str,
				tgpi->bottom[0], tgpi->bottom[1]);
		}
	}
	ED_workspace_status_text(C, status_str);
}

/* ----------------------- */

/* create a rectangle */
static void gp_primitive_rectangle(tGPDprimitive *tgpi, tGPspoint *points2D)
{
	BLI_assert(tgpi->tot_edges == 4);

	int i = tgpi->tot_stored_edges;

	points2D[i].x = (float)tgpi->top[0];
	points2D[i].y = (float)tgpi->top[1];

	points2D[i + 1].x = (float)tgpi->bottom[0];
	points2D[i + 1].y = (float)tgpi->top[1];

	points2D[i + 2].x = (float)tgpi->bottom[0];
	points2D[i + 2].y = (float)tgpi->bottom[1];

	points2D[i + 3].x = (float)tgpi->top[0];
	points2D[i + 3].y = (float)tgpi->bottom[1];
}

/* create a line */
static void gp_primitive_line(tGPDprimitive *tgpi, tGPspoint *points2D)
{
	BLI_assert(tgpi->tot_edges == 2);

	int i = tgpi->tot_stored_edges;

	points2D[i].x = (float)tgpi->top[0];
	points2D[i].y = (float)tgpi->top[1];

	points2D[i + 1].x = (float)tgpi->bottom[0];
	points2D[i + 1].y = (float)tgpi->bottom[1];
}

/* create an arc */
static void gp_primitive_arc(tGPDprimitive *tgpi, tGPspoint *points2D)
{
	const int totpoints = (tgpi->tot_edges + tgpi->tot_stored_edges);
	const float step = M_PI_2 / (float)(tgpi->tot_edges - 1);
	float length[2];
	int start[2];
	int end[2];
	float a = 0.0f;

	start[0] = tgpi->top[0];
	start[1] = tgpi->top[1];
	end[0] = tgpi->bottom[0];
	end[1] = tgpi->bottom[1];

	if (tgpi->flip) {
		SWAP(int, end[0], start[0]);
		SWAP(int, end[1], start[1]);
	}

	length[0] = end[0] - start[0];
	length[1] = end[1] - start[1];

	for (int i = tgpi->tot_stored_edges; i < totpoints; i++) {
		tGPspoint *p2d = &points2D[i];
		p2d->x = (start[0] + sinf(a) * length[0]);
		p2d->y = (end[1] - cosf(a) * length[1]);
		a += step;
	}
}

/* create a bezier */
static void gp_primitive_bezier(tGPDprimitive *tgpi, tGPspoint *points2D)
{
	const int totpoints = (tgpi->tot_edges + tgpi->tot_stored_edges);
	const float step = 1.0f / (float)(tgpi->tot_edges - 1);
	float cp1[2];
	float cp2[2];
	float cp3[2];
	float cp4[2];
	float a = 0.0f;

	copy_v2fl_v2i(cp1, tgpi->top);
	copy_v2fl_v2i(cp2, tgpi->bezcp1);
	copy_v2fl_v2i(cp3, tgpi->bezcp2);
	copy_v2fl_v2i(cp4, tgpi->bottom);

	for (int i = tgpi->tot_stored_edges; i < totpoints; i++) {
		tGPspoint *p2d = &points2D[i];
		interp_v2_v2v2v2v2_cubic(&p2d->x, cp1, cp2, cp3, cp4, a);
		a += step;
	}
}

/* create a circle */
static void gp_primitive_circle(tGPDprimitive *tgpi, tGPspoint *points2D)
{
	const int totpoints = (tgpi->tot_edges + tgpi->tot_stored_edges);
	const float step = (2.0f * M_PI) / (float)(tgpi->tot_edges);
	float center[2];
	float radius[2];
	float a = 0.0f;

	/* TODO: Use math-lib functions for these? */
	center[0] = tgpi->top[0] + ((tgpi->bottom[0] - tgpi->top[0]) / 2.0f);
	center[1] = tgpi->top[1] + ((tgpi->bottom[1] - tgpi->top[1]) / 2.0f);
	radius[0] = fabsf(((tgpi->bottom[0] - tgpi->top[0]) / 2.0f));
	radius[1] = fabsf(((tgpi->bottom[1] - tgpi->top[1]) / 2.0f));

	for (int i = tgpi->tot_stored_edges; i < totpoints; i++) {
		tGPspoint *p2d = &points2D[i];
		p2d->x = (center[0] + cosf(a) * radius[0]);
		p2d->y = (center[1] + sinf(a) * radius[1]);
		a += step;
	}
}

/* Helper: Update shape of the stroke */
static void gp_primitive_update_strokes(bContext *C, tGPDprimitive *tgpi)
{
	ToolSettings *ts = tgpi->scene->toolsettings;
	bGPdata *gpd = tgpi->gpd;
	bGPDstroke *gps = tgpi->gpf->strokes.first;
	GP_Sculpt_Settings *gset = &ts->gp_sculpt;

	/* realloc points to new size */
	/* TODO: only do this if the size has changed? */
	gps->points = MEM_reallocN(gps->points, sizeof(bGPDspoint) * (tgpi->tot_edges + tgpi->tot_stored_edges));
	if (gps->dvert != NULL) {
		gps->dvert = MEM_reallocN(gps->dvert, sizeof(MDeformVert) * (tgpi->tot_edges + tgpi->tot_stored_edges));
	}
	gps->totpoints = (tgpi->tot_edges + tgpi->tot_stored_edges);

	/* compute screen-space coordinates for points */
	tgpi->points = MEM_reallocN(tgpi->points, sizeof(tGPspoint) * (tgpi->tot_edges + tgpi->tot_stored_edges));
	tGPspoint *points2D = tgpi->points;
	switch (tgpi->type) {
		case GP_STROKE_BOX:
			tgpi->cyclic = true;
			gp_primitive_rectangle(tgpi, points2D);
			break;
		case GP_STROKE_LINE:
			tgpi->cyclic = false;
			gp_primitive_line(tgpi, points2D);
			break;
		case GP_STROKE_CIRCLE:
			tgpi->cyclic = true;
			gp_primitive_circle(tgpi, points2D);
			break;
		case GP_STROKE_ARC:
			gp_primitive_arc(tgpi, points2D);
			break;
		case GP_STROKE_BEZIER:
			gp_primitive_bezier(tgpi, points2D);
		default:
			break;
	}

	if (ELEM(tgpi->type, GP_STROKE_ARC, GP_STROKE_BEZIER)) {
		if (tgpi->cyclic)
			gps->flag |= GP_STROKE_CYCLIC;
		else
			gps->flag &= ~GP_STROKE_CYCLIC;
	}	

	/* convert screen-coordinates to 3D coordinates */
	gp_session_validatebuffer(tgpi);
	gp_init_colors(tgpi);
	if (gset->flag & GP_SCULPT_SETT_FLAG_PRIMITIVE_CURVE) {
		curvemapping_initialize(ts->gp_sculpt.cur_primitive);
	}
	
	for (int i = 0; i < gps->totpoints; i++) {
		bGPDspoint *pt = &gps->points[i];
		tGPspoint *p2d = &points2D[i];

		/* Copy points to buffer */
		tGPspoint *tpt = ((tGPspoint *)(gpd->runtime.sbuffer) + gpd->runtime.sbuffer_size);
		tpt->x = p2d->x;
		tpt->y = p2d->y;

		/* calc pressure */
		float pressure = 1.0;
		if (ELEM(tgpi->type, GP_STROKE_ARC, GP_STROKE_BEZIER)) {
			if (gset->flag & GP_SCULPT_SETT_FLAG_PRIMITIVE_CURVE) {
				/* normalize value to evaluate curve */
				float value = (float)i / (gps->totpoints - 1);
				float curvef = curvemapping_evaluateF(gset->cur_primitive, 0, value);
				pressure = 1.0f * curvef;
				CLAMP_MIN(pressure, 0.1f);
			}
		}

		tpt->pressure = pressure;
		tpt->strength = tgpi->brush->gpencil_settings->draw_strength;
		
		tpt->time = p2d->time;
		tpt->uv_fac = p2d->uv_fac;
		tpt->uv_rot = p2d->uv_rot;

		gpd->runtime.sbuffer_size++;

		/* convert screen-coordinates to 3D coordinates */
		gp_stroke_convertcoords_tpoint(tgpi->scene, tgpi->ar, tgpi->ob, tgpi->gpl, p2d, NULL, &pt->x);

		pt->pressure = pressure;
		pt->strength = tgpi->brush->gpencil_settings->draw_strength;
		pt->time = 0.0f;
		pt->flag = 0;

		if (gps->dvert != NULL) {
			MDeformVert *dvert = &gps->dvert[i];
			dvert->totweight = 0;
			dvert->dw = NULL;
		}
	}

	/* if axis locked, reproject to plane locked */
	if (tgpi->lock_axis > GP_LOCKAXIS_VIEW) {
		bGPDspoint *tpt = gps->points;
		float origin[3];
		ED_gp_get_drawing_reference(tgpi->scene, tgpi->ob, tgpi->gpl,
			ts->gpencil_v3d_align, origin);

		for (int i = 0; i < gps->totpoints; i++, tpt++) {
			ED_gp_project_point_to_plane(tgpi->ob, tgpi->rv3d, origin,
				ts->gp_sculpt.lock_axis - 1,
				tpt);
		}
	}

	/* if parented change position relative to parent object */
	for (int i = 0; i < gps->totpoints; i++) {
		bGPDspoint *pt = &gps->points[i];
		gp_apply_parent_point(tgpi->depsgraph, tgpi->ob, tgpi->gpd, tgpi->gpl, pt);
	}

	/* force fill recalc */
	gps->flag |= GP_STROKE_RECALC_CACHES;

	DEG_id_tag_update(&gpd->id, ID_RECALC_COPY_ON_WRITE);
	DEG_id_tag_update(&gpd->id, ID_RECALC_TRANSFORM | ID_RECALC_GEOMETRY);
	WM_event_add_notifier(C, NC_GPENCIL | NA_EDITED, NULL);
}

/* add new segment to curve */
static void gpencil_primitive_add_segment(tGPDprimitive *tgpi)
{
	tgpi->tot_stored_edges += tgpi->tot_edges;
}

/* Update screen and stroke */
static void gpencil_primitive_update(bContext *C, wmOperator *op, tGPDprimitive *tgpi)
{
	/* update indicator in header */
	gpencil_primitive_status_indicators(C, tgpi);
	/* apply... */
	tgpi->type = RNA_enum_get(op->ptr, "type");
	tgpi->tot_edges = RNA_int_get(op->ptr, "edges");
	/* update points position */
	gp_primitive_update_strokes(C, tgpi);
}

/* ----------------------- */

static void gpencil_primitive_interaction_begin(tGPDprimitive *tgpi, const wmEvent *event)
{
	copy_v2_v2_int(tgpi->origin, event->mval);
	copy_v2_v2_int(tgpi->top, event->mval);
	copy_v2_v2_int(tgpi->bottom, event->mval);
	copy_v2_v2_int(tgpi->bezcp1, event->mval);
	copy_v2_v2_int(tgpi->bezcp2, event->mval);
}

/* Exit and free memory */
static void gpencil_primitive_exit(bContext *C, wmOperator *op)
{
	tGPDprimitive *tgpi = op->customdata;
	bGPdata *gpd = tgpi->gpd;

	/* don't assume that operator data exists at all */
	if (tgpi) {
		/* remove drawing handler */
		if (tgpi->draw_handle_3d) {
			ED_region_draw_cb_exit(tgpi->ar->type, tgpi->draw_handle_3d);
		}

		/* clear status message area */
		ED_workspace_status_text(C, NULL);

		MEM_SAFE_FREE(tgpi->points);
		/* finally, free memory used by temp data */
		BKE_gpencil_free_strokes(tgpi->gpf);
		MEM_SAFE_FREE(tgpi->gpf);
		MEM_freeN(tgpi);
	}

	/* free stroke buffer */
	if ((gpd != NULL) && (gpd->runtime.sbuffer)) {
		MEM_SAFE_FREE(gpd->runtime.sbuffer);
		gpd->runtime.sbuffer = NULL;

		/* clear flags */
		gpd->runtime.sbuffer_size = 0;
		gpd->runtime.sbuffer_sflag = 0;
	}

	DEG_id_tag_update(&gpd->id, ID_RECALC_TRANSFORM | ID_RECALC_GEOMETRY);
	WM_event_add_notifier(C, NC_GPENCIL | NA_EDITED, NULL);

	/* clear pointer */
	op->customdata = NULL;
}

/* Init new temporary primitive data */
static void gpencil_primitive_init(bContext *C, wmOperator *op)
{
	ToolSettings *ts = CTX_data_tool_settings(C);
	bGPdata *gpd = CTX_data_gpencil_data(C);
	Main *bmain = CTX_data_main(C);
	Scene *scene = CTX_data_scene(C);
	Depsgraph *depsgraph = CTX_data_depsgraph(C);
	int cfra_eval = (int)DEG_get_ctime(depsgraph);

	/* create temporary operator data */
	tGPDprimitive *tgpi = MEM_callocN(sizeof(tGPDprimitive), "GPencil Primitive Data");
	op->customdata = tgpi;

	tgpi->points = MEM_callocN(sizeof(tGPspoint), "gp primitive points2D");

	/* set current scene and window info */
	tgpi->bmain = CTX_data_main(C);
	tgpi->scene = scene;
	tgpi->ob = CTX_data_active_object(C);
	tgpi->sa = CTX_wm_area(C);
	tgpi->ar = CTX_wm_region(C);
	tgpi->rv3d = tgpi->ar->regiondata;
	tgpi->v3d = tgpi->sa->spacedata.first;
	tgpi->depsgraph = CTX_data_depsgraph(C);
	tgpi->win = CTX_wm_window(C);

	/* set current frame number */
	tgpi->cframe = cfra_eval;

	/* set GP datablock */
	tgpi->gpd = gpd;

	/* getcolor info */
	tgpi->mat = BKE_gpencil_material_ensure(bmain, tgpi->ob);

	/* set parameters */
	tgpi->type = RNA_enum_get(op->ptr, "type");

	if(ELEM(tgpi->type, GP_STROKE_ARC, GP_STROKE_BEZIER))
		tgpi->curve = true;
	else
		tgpi->curve = false;

	/* set default edge count */
	if (tgpi->type == GP_STROKE_CIRCLE) {
		RNA_int_set(op->ptr, "edges", 64);
	}
	else if (tgpi->curve) {
		RNA_int_set(op->ptr, "edges", 32);
	}
	else if (tgpi->type == GP_STROKE_BOX) {
		RNA_int_set(op->ptr, "edges", 4);
	}
	else { /* LINE */
		RNA_int_set(op->ptr, "edges", 2);
	}

	tgpi->tot_stored_edges = 0;
	tgpi->tot_edges = RNA_int_get(op->ptr, "edges");
	tgpi->flag = IDLE;

	tgpi->lock_axis = ts->gp_sculpt.lock_axis;

	/* set temp layer, frame and stroke */
	gp_primitive_set_initdata(C, tgpi);
}

/* ----------------------- */

/* Invoke handler: Initialize the operator */
static int gpencil_primitive_invoke(bContext *C, wmOperator *op, const wmEvent *event)
{
	wmWindow *win = CTX_wm_window(C);
	bGPdata *gpd = CTX_data_gpencil_data(C);
	tGPDprimitive *tgpi = NULL;

	/* initialize operator runtime data */
	gpencil_primitive_init(C, op);
	tgpi = op->customdata;

	const bool is_modal = RNA_boolean_get(op->ptr, "wait_for_input");
	if (!is_modal) {
		tgpi->flag = IN_PROGRESS;
		gpencil_primitive_interaction_begin(tgpi, event);
	}

	/* if in tools region, wait till we get to the main (3d-space)
	 * region before allowing drawing to take place.
	 */
	op->flag |= OP_IS_MODAL_CURSOR_REGION;

	/* Enable custom drawing handlers */
	tgpi->draw_handle_3d = ED_region_draw_cb_activate(tgpi->ar->type, gpencil_primitive_draw_3d, tgpi, REGION_DRAW_POST_VIEW);

	/* set cursor to indicate modal */
	WM_cursor_modal_set(win, BC_CROSSCURSOR);

	/* update sindicator in header */
	gpencil_primitive_status_indicators(C, tgpi);
	DEG_id_tag_update(&gpd->id, ID_RECALC_TRANSFORM | ID_RECALC_GEOMETRY);
	WM_event_add_notifier(C, NC_GPENCIL | NA_EDITED, NULL);

	/* add a modal handler for this operator */
	WM_event_add_modal_handler(C, op);

	return OPERATOR_RUNNING_MODAL;
}

/* Helper to complete a primitive */
static void gpencil_primitive_interaction_end(bContext *C, wmOperator *op, wmWindow *win, tGPDprimitive *tgpi)
{
	bGPDframe *gpf;
	bGPDstroke *gps;

	ToolSettings *ts = tgpi->scene->toolsettings;

	const int def_nr = tgpi->ob->actdef - 1;
	const bool have_weight = (bool)BLI_findlink(&tgpi->ob->defbase, def_nr);

	/* return to normal cursor and header status */
	ED_workspace_status_text(C, NULL);
	WM_cursor_modal_restore(win);

	/* insert keyframes as required... */
	gpf = BKE_gpencil_layer_getframe(tgpi->gpl, tgpi->cframe, GP_GETFRAME_ADD_NEW);

	/* prepare stroke to get transferred */
	gps = tgpi->gpf->strokes.first;
	if (gps) {
		gps->thickness = tgpi->brush->size;
		gps->flag |= GP_STROKE_RECALC_CACHES;
		gps->tot_triangles = 0;
	}

	/* transfer stroke from temporary buffer to the actual frame */
	BLI_movelisttolist(&gpf->strokes, &tgpi->gpf->strokes);
	BLI_assert(BLI_listbase_is_empty(&tgpi->gpf->strokes));

	/* add weights if required */
	if ((ts->gpencil_flags & GP_TOOL_FLAG_CREATE_WEIGHTS) && (have_weight)) {
		BKE_gpencil_dvert_ensure(gps);
		for (int i = 0; i < gps->totpoints; i++) {
			MDeformVert *ve = &gps->dvert[i];
			MDeformWeight *dw = defvert_verify_index(ve, def_nr);
			if (dw) {
				dw->weight = ts->vgroup_weight;
			}

		}
	}

	DEG_id_tag_update(&tgpi->gpd->id, ID_RECALC_COPY_ON_WRITE);
	DEG_id_tag_update(&tgpi->gpd->id, ID_RECALC_TRANSFORM | ID_RECALC_GEOMETRY);

	/* clean up temp data */
	gpencil_primitive_exit(C, op);
}

/* Helper to set bezier cp */
static void gpencil_primitive_set_midpoint(tGPDprimitive *tgpi)
{
	float midpoint[2], start[2], end[2];
	copy_v2fl_v2i(start, tgpi->top);
	copy_v2fl_v2i(end, tgpi->bottom);
	mid_v2_v2v2(midpoint, start, end);
	round_v2i_v2fl(tgpi->bezcp1, midpoint);
	copy_v2_v2_int(tgpi->bezcp2, tgpi->bezcp1);
}

/* Helper to square a primitive */
static void gpencil_primitive_to_square(tGPDprimitive *tgpi, const int x, const int y)
{
	int w = abs(x);
	int h = abs(y);
	if ((x > 0 && y > 0) || (x < 0 && y < 0)) {
		if (w > h)
			tgpi->bottom[1] = tgpi->origin[1] + x;
		else
			tgpi->bottom[0] = tgpi->origin[0] + y;
	}
	else {
		if (w > h)
			tgpi->bottom[1] = tgpi->origin[1] - x;
		else
			tgpi->bottom[0] = tgpi->origin[0] - y;
	}
}

#define MOVE_NONE 0
#define MOVE_ENDS 1
#define MOVE_CP 2

/* bezier event handling */
static void gpencil_primitive_bezier_event_handling(bContext *C, wmOperator *op, wmWindow *win, const wmEvent *event, tGPDprimitive *tgpi) {
	/* calculate nearest point then set cursor */
	int move = MOVE_NONE;
	float a = len_v2v2_int(event->mval, tgpi->top);
	float b = len_v2v2_int(event->mval, tgpi->bottom);

	if (tgpi->flag == IN_CURVE_EDIT_BEZIER) {
		if ((a < 10 && tgpi->tot_stored_edges == 0) || b < 10) {
			move = MOVE_ENDS;
			WM_cursor_modal_set(win, BC_RING_CURSOR);
		}
		else {
			move = MOVE_CP;
			WM_cursor_modal_set(win, BC_HANDCURSOR);
		}
	}

	switch (event->type) {
	case MOUSEMOVE:
		if ((event->val == KM_PRESS) && tgpi->sel_cp != CURVE_CP_NONE) {
			if (tgpi->sel_cp == CURVE_CP_T && tgpi->tot_stored_edges == 0) {
				copy_v2_v2_int(tgpi->top, event->mval);
			}
			else if (tgpi->sel_cp == CURVE_CP_B) {
				copy_v2_v2_int(tgpi->bottom, event->mval);
			}
			else if (tgpi->sel_cp == CURVE_CP_1) {
				int dx = (event->mval[0] - tgpi->mvalo[0]);
				int dy = (event->mval[1] - tgpi->mvalo[1]);
				tgpi->bezcp1[0] += dx;
				tgpi->bezcp1[1] += dy;
				if (event->ctrl)
					copy_v2_v2_int(tgpi->bezcp2, tgpi->bezcp1);
			}
			else if (tgpi->sel_cp == CURVE_CP_2) {
				int dx = (event->mval[0] - tgpi->mvalo[0]);
				int dy = (event->mval[1] - tgpi->mvalo[1]);
				tgpi->bezcp2[0] += dx;
				tgpi->bezcp2[1] += dy;
				if (event->ctrl)
					copy_v2_v2_int(tgpi->bezcp1, tgpi->bezcp2);
			}
			copy_v2_v2_int(tgpi->mvalo, event->mval);

			/* update screen */
			gpencil_primitive_update(C, op, tgpi);
		}
		else if ((event->val == KM_PRESS)) {
			gpencil_primitive_set_midpoint(tgpi);
			gpencil_primitive_update(C, op, tgpi);
		}
		break;
	case LEFTMOUSE:
		if ((event->val == KM_PRESS)) {
			copy_v2_v2_int(tgpi->mvalo, event->mval);
			/* find nearest cp based on stroke end points */
			if (move == MOVE_ENDS)
				tgpi->sel_cp = (a < b) ? CURVE_CP_T : CURVE_CP_B;
			else if (move == MOVE_CP)
				tgpi->sel_cp = (a < b) ? CURVE_CP_1 : CURVE_CP_2;
			else
				tgpi->sel_cp = CURVE_CP_NONE;
			break;
		}
		else if ((event->val == KM_RELEASE) && (tgpi->flag == IN_PROGRESS)) {
			/* set control points and enter edit mode */
			tgpi->flag = IN_CURVE_EDIT_BEZIER;
			gpencil_primitive_set_midpoint(tgpi);
			copy_v2_v2_int(tgpi->mvalo, event->mval);
			gpencil_primitive_update(C, op, tgpi);
		}
		else {
			tgpi->sel_cp = CURVE_CP_NONE;
		}
		break;
	case AKEY:
		if (tgpi->flag == IN_CURVE_EDIT_BEZIER) {
			tgpi->flag = IN_PROGRESS;
			gpencil_primitive_add_segment(tgpi);
			copy_v2_v2_int(tgpi->top, tgpi->bottom);
			copy_v2_v2_int(tgpi->origin, tgpi->top);
			gpencil_primitive_set_midpoint(tgpi);
			copy_v2_v2_int(tgpi->mvalo, event->mval);
		}
		break;
	}
}

/* Modal handler: Events handling during interactive part */
static int gpencil_primitive_modal(bContext *C, wmOperator *op, const wmEvent *event)
{
	tGPDprimitive *tgpi = op->customdata;
	wmWindow *win = CTX_wm_window(C);
	const bool has_numinput = hasNumInput(&tgpi->num);

	/* bezier event handling */
	if (tgpi->type == GP_STROKE_BEZIER) {
		gpencil_primitive_bezier_event_handling(C, op, win, event, tgpi);
	}

	switch (event->type) {
		case LEFTMOUSE:
			if ((event->val == KM_PRESS) && (tgpi->flag == IDLE)) {
				/* start drawing primitive */
				/* TODO: Ignore if not in main region yet */
				tgpi->flag = IN_PROGRESS;
				gpencil_primitive_interaction_begin(tgpi, event);
			}
			else if ((event->val == KM_RELEASE) && (tgpi->flag == IN_PROGRESS) && (tgpi->type != GP_STROKE_BEZIER)) {
				/* stop drawing primitive */
				tgpi->flag = IDLE;
				gpencil_primitive_interaction_end(C, op, win, tgpi);
				/* done! */
				return OPERATOR_FINISHED;
			}
			else {
				if (G.debug & G_DEBUG) {
					printf("GP Add Primitive Modal: LEFTMOUSE %d, Status = %d\n", event->val, tgpi->flag);
				}
			}
			break;
		case SPACEKEY:  /* confirm */
		case RETKEY:
		{
			tgpi->flag = IDLE;
			gpencil_primitive_interaction_end(C, op, win, tgpi);
			/* done! */
			return OPERATOR_FINISHED;
		}
		case RIGHTMOUSE:
		if (tgpi->flag == IN_CURVE_EDIT_BEZIER) {
			tgpi->flag = IDLE;
			gpencil_primitive_update(C, op, tgpi);
			gpencil_primitive_interaction_end(C, op, win, tgpi);
			/* done! */
			return OPERATOR_FINISHED;
		}
		case ESCKEY:
		{
			/* return to normal cursor and header status */
			ED_workspace_status_text(C, NULL);
			WM_cursor_modal_restore(win);

			/* clean up temp data */
			gpencil_primitive_exit(C, op);

			/* canceled! */
			return OPERATOR_CANCELLED;
		}
		case CKEY:
		{
			if ((event->val == KM_RELEASE) && tgpi->type == GP_STROKE_ARC) {
				tgpi->cyclic ^= 1;

				/* update screen */
				gpencil_primitive_update(C, op, tgpi);
			}
			break;
		}
		case FKEY:
		{
			if ((event->val == KM_RELEASE) && tgpi->type == GP_STROKE_ARC) {
				tgpi->flip ^= 1;

				/* update screen */
				gpencil_primitive_update(C, op, tgpi);
			}
			break;
		}
		case PADPLUSKEY:
		case WHEELUPMOUSE:
		{
			if ((event->val != KM_RELEASE) && (tgpi->type == GP_STROKE_CIRCLE || tgpi->type == GP_STROKE_ARC)) {
				tgpi->tot_edges = tgpi->tot_edges + 1;
				CLAMP(tgpi->tot_edges, MIN_EDGES, MAX_EDGES);
				RNA_int_set(op->ptr, "edges", tgpi->tot_edges);

				/* update screen */
				gpencil_primitive_update(C, op, tgpi);
			}
			break;
		}
		case PADMINUS:
		case WHEELDOWNMOUSE:
		{
			if ((event->val != KM_RELEASE) && (tgpi->type == GP_STROKE_CIRCLE || tgpi->type == GP_STROKE_ARC)) {
				tgpi->tot_edges = tgpi->tot_edges - 1;
				CLAMP(tgpi->tot_edges, MIN_EDGES, MAX_EDGES);
				RNA_int_set(op->ptr, "edges", tgpi->tot_edges);

				/* update screen */
				gpencil_primitive_update(C, op, tgpi);
			}
			break;
		}
		case MOUSEMOVE: /* calculate new position */
		{
			if (tgpi->flag == IN_CURVE_EDIT_BEZIER) {
				break;
			}

			/* only handle mousemove if not doing numinput */
			if (has_numinput == false) {
				/* update position of mouse */
				tgpi->bottom[0] = event->mval[0];
				tgpi->bottom[1] = event->mval[1];
				tgpi->top[0] = tgpi->origin[0];
				tgpi->top[1] = tgpi->origin[1];
				if (tgpi->flag == IDLE) {
					tgpi->origin[0] = event->mval[0];
					tgpi->origin[1] = event->mval[1];
				}
				/* Keep square if shift key */
				if (event->shift) {
					int x = tgpi->bottom[0] - tgpi->origin[0];
					int y = tgpi->bottom[1] - tgpi->origin[1];
					if (tgpi->type == GP_STROKE_LINE || tgpi->curve) {
						float angle = fabsf(atan2f((float)y, (float)x));
						if (angle < 0.4f || angle > (M_PI - 0.4f)) {
							tgpi->bottom[1] = tgpi->origin[1];
						}
						else if (angle > (M_PI_2 - 0.4f) && angle < (M_PI_2 + 0.4f)) {
							tgpi->bottom[0] = tgpi->origin[0];
						}
						else {
							gpencil_primitive_to_square(tgpi, x, y);
						}
					}
					else {
						gpencil_primitive_to_square(tgpi, x, y);
					}
				}
				/* Center primitive if alt key */
				if (event->alt) {
					tgpi->top[0] = tgpi->origin[0] - (tgpi->bottom[0] - tgpi->origin[0]);
					tgpi->top[1] = tgpi->origin[1] - (tgpi->bottom[1] - tgpi->origin[1]);
				}
				/* update screen */
				gpencil_primitive_update(C, op, tgpi);
			}
			break;
		}
		default:
		{
			if (tgpi->flag != IN_CURVE_EDIT_BEZIER && (event->val == KM_PRESS) && handleNumInput(C, &tgpi->num, event)) {
				float value;

				/* Grab data from numeric input, and store this new value (the user see an int) */
				value = tgpi->tot_edges;
				applyNumInput(&tgpi->num, &value);
				tgpi->tot_edges = value;

				CLAMP(tgpi->tot_edges, MIN_EDGES, MAX_EDGES);
				RNA_int_set(op->ptr, "edges", tgpi->tot_edges);

				/* update screen */
				gpencil_primitive_update(C, op, tgpi);

				break;
			}
			else {
				/* unhandled event - allow to pass through */
				return OPERATOR_RUNNING_MODAL | OPERATOR_PASS_THROUGH;
			}
		}
	}

	/* still running... */
	return OPERATOR_RUNNING_MODAL;
}

/* Cancel handler */
static void gpencil_primitive_cancel(bContext *C, wmOperator *op)
{
	/* this is just a wrapper around exit() */
	gpencil_primitive_exit(C, op);
}

void GPENCIL_OT_primitive(wmOperatorType *ot)
{
	static EnumPropertyItem primitive_type[] = {
		{GP_STROKE_BOX, "BOX", 0, "Box", ""},
		{GP_STROKE_LINE, "LINE", 0, "Line", ""},
		{GP_STROKE_CIRCLE, "CIRCLE", 0, "Circle", ""},
		{GP_STROKE_ARC, "ARC", 0, "Arc", ""},
		{GP_STROKE_BEZIER, "BEZIER", 0, "Bezier", ""},
		{0, NULL, 0, NULL, NULL}
	};

	/* identifiers */
	ot->name = "Grease Pencil Shapes";
	ot->idname = "GPENCIL_OT_primitive";
	ot->description = "Create predefined grease pencil stroke shapes";

	/* callbacks */
	ot->invoke = gpencil_primitive_invoke;
	ot->modal = gpencil_primitive_modal;
	ot->cancel = gpencil_primitive_cancel;
	ot->poll = gpencil_primitive_add_poll;

	/* flags */
	ot->flag = OPTYPE_UNDO | OPTYPE_BLOCKING;

	/* properties */
	PropertyRNA *prop;

	RNA_def_int(ot->srna, "edges", 4, MIN_EDGES, MAX_EDGES, "Edges", "Number of polygon edges", MIN_EDGES, MAX_EDGES);
	RNA_def_enum(ot->srna, "type", primitive_type, GP_STROKE_BOX, "Type", "Type of shape");

	prop = RNA_def_boolean(ot->srna, "wait_for_input", true, "Wait for Input", "");
	RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);
}

/* *************************************************************** */
