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

/** \file blender/editors/hair/hair_ops.c
 *  \ingroup edhair
 */


#include <stdlib.h>
#include <math.h>

#include "DNA_hair_types.h"
#include "DNA_scene_types.h"

#include "RNA_access.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_hair.h"
#include "ED_object.h"
#include "ED_screen.h"
#include "ED_transform.h"

#include "hair_intern.h"

/************************* registration ****************************/

void ED_operatortypes_hair(void)
{
	WM_operatortype_append(HAIR_OT_select_all);

	WM_operatortype_append(HAIR_OT_add_test_hair);
}

void ED_operatormacros_hair(void)
{
	wmOperatorType *ot;
	wmOperatorTypeMacro *otmacro;

	UNUSED_VARS(ot, otmacro);
}

void ED_keymap_hair(wmKeyConfig *keyconf)
{
	wmKeyMap *keymap;
	wmKeyMapItem *kmi;

	UNUSED_VARS(kmi);
	
	keymap = WM_keymap_ensure(keyconf, "Hair", 0, 0);
	keymap->poll = ED_operator_edithair;
	
	ED_keymap_proportional_cycle(keyconf, keymap);
	ED_keymap_proportional_editmode(keyconf, keymap, true);
}
