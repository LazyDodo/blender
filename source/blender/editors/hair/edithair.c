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
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/hair/edithair.c
 *  \ingroup edhair
 */

#include "DNA_hair_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_array_utils.h"
#include "BLI_blenlib.h"

#include "BKE_context.h"
#include "BKE_global.h"
#include "BKE_hair.h"
#include "BKE_library.h"
#include "BKE_main.h"
#include "BKE_report.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_hair.h"
#include "ED_object.h"
#include "ED_screen.h"
#include "ED_types.h"
#include "ED_util.h"
#include "ED_view3d.h"

#include "hair_intern.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

/********************** Load/Make/Free ********************/

void ED_hair_edithair_make(Object *obedit)
{
	HairSystem *hsys = obedit->data;

	ED_hair_edithair_free(obedit);

	hsys->edithair = MEM_callocN(sizeof(EditHair), "edithair");

	BKE_hair_curve_data_copy(&hsys->edithair->curve_data, &hsys->curve_data, 0);
}

void ED_hair_edithair_load(Object *obedit)
{
	HairSystem *hsys = obedit->data;
	
	BKE_hair_curve_data_free(&hsys->curve_data);
	BKE_hair_curve_data_copy(&hsys->curve_data, &hsys->edithair->curve_data, 0);
}

void ED_hair_edithair_free(Object *ob)
{
	HairSystem *hsys = ob->data;
	
	if (hsys->edithair) {
		BKE_hair_curve_data_free(&hsys->edithair->curve_data);
		
		MEM_freeN(hsys->edithair);
		hsys->edithair = NULL;
	}
}

bool ED_hair_poll_object(bContext *C)
{
	Object *ob = ED_object_context(C);
	return ob && ob->type == OB_HAIR;
}


bool ED_hair_poll_editmode(bContext *C)
{
	Object *ob = ED_object_context(C);
	return ob && ob->type == OB_HAIR && (ob->mode & OB_MODE_EDIT);
}

bool ED_hair_poll_view3d(bContext *C)
{
	if (!ED_hair_poll_editmode) {
		return false;
	}

	ScrArea *sa = CTX_wm_area(C);
	ARegion *ar = CTX_wm_region(C);
	return (sa && sa->spacetype == SPACE_VIEW3D) &&
	       (ar && ar->regiontype == RGN_TYPE_WINDOW);
}

void ED_hair_init_view3d(bContext *C, ViewContext *vc)
{
	ED_view3d_viewcontext_init(C, vc);

	if (V3D_IS_ZBUF(vc->v3d)) {
		if (vc->v3d->flag & V3D_INVALID_BACKBUF) {
			/* needed or else the draw matrix can be incorrect */
			view3d_operator_needs_opengl(C);

			ED_view3d_backbuf_validate(vc);
			/* we may need to force an update here by setting the rv3d as dirty
			 * for now it seems ok, but take care!:
			 * rv3d->depths->dirty = 1; */
			ED_view3d_depth_update(vc->ar);
		}
	}
}
