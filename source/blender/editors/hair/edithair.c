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

	hsys->edithair->pattern = BKE_hair_pattern_copy(hsys->pattern, 0);
	BKE_hair_curve_data_copy(&hsys->edithair->curve_data, &hsys->curve_data, 0);
}

void ED_hair_edithair_load(Object *obedit)
{
	HairSystem *hsys = obedit->data;
	
	BKE_hair_pattern_free(hsys->pattern);
	BKE_hair_curve_data_free(&hsys->curve_data);
	hsys->pattern = BKE_hair_pattern_copy(hsys->edithair->pattern, 0);
	BKE_hair_curve_data_copy(&hsys->curve_data, &hsys->edithair->curve_data, 0);
}

void ED_hair_edithair_free(Object *ob)
{
	HairSystem *hsys = ob->data;
	
	if (hsys->edithair) {
		BKE_hair_pattern_free(hsys->edithair->pattern);
		BKE_hair_curve_data_free(&hsys->edithair->curve_data);
		
		MEM_freeN(hsys->edithair);
		hsys->edithair = NULL;
	}
}

int ED_hair_object_poll(bContext *C)
{
	Object *ob = ED_object_context(C);
	return ob && ob->type == OB_HAIR;
}
