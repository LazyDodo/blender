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
 * The Original Code is Copyright (C) 2014 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Germano Cavalcante
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file snap3d_gizmo.c
 *  \ingroup edgizmolib
 *
 * \name Snap Gizmo
 *
 * 3D Gizmo
 *
 * \brief Simple snap gizmo which exposes the location of the snap point through the `ret_location` parameter.
 */

#include "BLI_math.h"

#include "DNA_scene_types.h"

#include "BKE_context.h"

#include "GPU_immediate.h"

#include "ED_view3d.h"
#include "ED_gizmo_library.h"
#include "ED_screen.h"
#include "ED_transform_snap_object_context.h"

#include "UI_resources.h" /* icons */

#include "RNA_access.h"
#include "RNA_define.h"

#include "WM_types.h"
#include "WM_api.h"

/* own includes */
#include "../gizmo_geometry.h"
#include "../gizmo_library_intern.h"


typedef struct SnapGizmo3D {
	wmGizmo gizmo;

	/* We could have other snap contexts, for now only support 3D view. */
	struct SnapObjectContext *snap_context_v3d;

	void *last_operator;
} SnapGizmo3D;


/* -------------------------------------------------------------------- */

static void gizmo_snap_setup(wmGizmo *gz)
{
	/* Default properties */
	WM_gizmo_target_property_def_rna(gz, "snap_elements", gz->ptr, "ret_snap_elements", -1);
	WM_gizmo_target_property_def_rna(gz, "location", gz->ptr, "ret_location", -1);
}

static void gizmo_snap_draw(const bContext *UNUSED(C), wmGizmo *gz)
{
	if ((gz->state & WM_GIZMO_STATE_HIGHLIGHT) == 0) {
		return;
	}

	float location[3];
	wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "location");
	WM_gizmo_target_property_float_get_array(gz, gz_prop, location);
	//RNA_float_get_array(gizmo_snap->gizmo.ptr, "location", location);

	uint pos = GPU_vertformat_attr_add(immVertexFormat(), "pos", GPU_COMP_F32, 3, GPU_FETCH_FLOAT);
	immBindBuiltinProgram(GPU_SHADER_3D_POINT_FIXED_SIZE_UNIFORM_COLOR);

	immUniformColor3f(1.0, 1.0, 1.0);

	immBegin(GPU_PRIM_POINTS, 1);
	immVertex3fv(pos, location);
	immEnd();

	immUnbindProgram();
}

static int gizmo_snap_test_select(
        bContext *C, wmGizmo *gz, const int mval[2])
{
	SnapGizmo3D *gizmo_snap = (SnapGizmo3D *)gz;
	wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "snap_elements");

	int snap_elements = RNA_property_enum_get(&gz_prop->ptr, gz_prop->prop);
	snap_elements &= ~(SCE_SNAP_MODE_INCREMENT | SCE_SNAP_MODE_VOLUME);
	if (!snap_elements) {
		return -1;
	}

	void *last_operator = CTX_wm_manager(C)->operators.last;
	if (last_operator != gizmo_snap->last_operator) {
		/* Something has changed since the last time.
		 * Has the mesh been changed?
		 * In the doubt we will clear the snap context. */
		if (gizmo_snap->snap_context_v3d) {
			ED_transform_snap_object_context_destroy(gizmo_snap->snap_context_v3d);
			gizmo_snap->snap_context_v3d = NULL;
		}
		gizmo_snap->last_operator = last_operator;

		/* return early to ensure that all objects were updated in time. */
		return -1;
	}

	ARegion *ar = CTX_wm_region(C);
	View3D *v3d = CTX_wm_view3d(C);
	if (gizmo_snap->snap_context_v3d == NULL) {
		gizmo_snap->snap_context_v3d = ED_transform_snap_object_context_create_view3d(
		        NULL, CTX_data_scene(C), CTX_data_depsgraph(C), 0, ar, v3d);
	}

	gz_prop = WM_gizmo_target_property_find(gz, "location");
	const float mval_fl[2] = {UNPACK2(mval)};
	float dist_px = 12.0f * U.pixelsize;
	float co[3];
	if (!ED_transform_snap_object_project_view3d(
	        gizmo_snap->snap_context_v3d,
	        snap_elements,
	        &(const struct SnapObjectParams){
	            .snap_select = SNAP_ALL,
	            .use_object_edit_cage = true,
	            .use_occlusion_test = true,
	        },
	        mval_fl, &dist_px,
	        co, NULL))
	{
		RegionView3D *rv3d = ar->regiondata;
		ED_view3d_win_to_3d(v3d, ar, rv3d->ofs, mval_fl, co);
	}

	WM_gizmo_target_property_float_set_array(C, gz, gz_prop, co);
	ED_region_tag_redraw(CTX_wm_region(C));

	return 0;
}

static int gizmo_snap_modal(
        bContext *UNUSED(C), wmGizmo *UNUSED(gz), const wmEvent *UNUSED(event),
        eWM_GizmoFlagTweak UNUSED(tweak_flag))
{
	return OPERATOR_RUNNING_MODAL;
}

static int gizmo_snap_invoke(
        bContext *UNUSED(C), wmGizmo *UNUSED(gz), const wmEvent *UNUSED(event))
{
	return OPERATOR_RUNNING_MODAL;
}

static void gizmo_snap_free(wmGizmo *gz)
{
	SnapGizmo3D *gizmo_snap = (SnapGizmo3D *)gz;
	if (gizmo_snap->snap_context_v3d) {
		ED_transform_snap_object_context_destroy(gizmo_snap->snap_context_v3d);
		gizmo_snap->snap_context_v3d = NULL;
	}
}

static void GIZMO_GT_snap_3d(wmGizmoType *gzt)
{
	/* identifiers */
	gzt->idname = "GIZMO_GT_snap_3d";

	/* api callbacks */
	gzt->setup = gizmo_snap_setup;
	gzt->draw = gizmo_snap_draw;
	gzt->test_select = gizmo_snap_test_select;
	gzt->modal = gizmo_snap_modal;
	gzt->invoke = gizmo_snap_invoke;
	gzt->free = gizmo_snap_free;

	gzt->struct_size = sizeof(SnapGizmo3D);

	/* Copy of EnumPropertyItem located in "rna_scene.c" */
	static EnumPropertyItem rna_enum_snap_element_items[] = {
		{SCE_SNAP_MODE_INCREMENT, "INCREMENT", ICON_SNAP_INCREMENT, "Increment", "Snap to increments of grid"},
		{SCE_SNAP_MODE_VERTEX, "VERTEX", ICON_SNAP_VERTEX, "Vertex", "Snap to vertices"},
		{SCE_SNAP_MODE_EDGE, "EDGE", ICON_SNAP_EDGE, "Edge", "Snap to edges"},
		{SCE_SNAP_MODE_FACE, "FACE", ICON_SNAP_FACE, "Face", "Snap to faces"},
		{SCE_SNAP_MODE_VOLUME, "VOLUME", ICON_SNAP_VOLUME, "Volume", "Snap to volume"},
		{0, NULL, 0, NULL, NULL}
	};

	RNA_def_enum_flag(
	        gzt->srna, "ret_snap_elements", rna_enum_snap_element_items,
	        SCE_SNAP_MODE_VERTEX | SCE_SNAP_MODE_EDGE | SCE_SNAP_MODE_FACE,
	        "Snap Elements", "");

	RNA_def_float_vector(gzt->srna, "ret_location", 3, NULL, FLT_MIN, FLT_MAX, "Location", "Snap Point Location", FLT_MIN, FLT_MAX);

	WM_gizmotype_target_property_def(gzt, "snap_elements", PROP_ENUM, 1);
	WM_gizmotype_target_property_def(gzt, "location", PROP_FLOAT, 3);
}

void ED_gizmotypes_snap_3d(void)
{
	WM_gizmotype_append(GIZMO_GT_snap_3d);
}

/** \} */
