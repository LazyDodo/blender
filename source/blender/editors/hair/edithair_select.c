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

/** \file blender/editors/hair/edithair_test.c
 *  \ingroup edhair
 */

#include "MEM_guardedalloc.h"

#include "BLI_math.h"

#include "DNA_hair_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_context.h"
#include "BKE_hair.h"
#include "BKE_hair_iterators.h"
#include "BKE_mesh_sample.h"

#include "DEG_depsgraph.h"

#include "ED_hair.h"
#include "ED_screen.h"
#include "ED_select_utils.h"
#include "ED_view3d.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "WM_api.h"
#include "WM_types.h"

#include "UI_resources.h"

#include "BLT_translation.h"

#include "hair_intern.h"  /* own include */

/************************ de select all operator ************************/

static bool hair_has_selected_follicles(EditHair *edit)
{
	const HairFollicle *follicle;
	HairIterator iter;
	BKE_HAIR_ITER_FOLLICLES(follicle, &iter, edit->pattern) {
		if (follicle->flag & HAIR_FOLLICLE_SELECT) {
			return true;
		}
	}
	return false;
}

static void hair_follicle_select_action_apply(HairFollicle *follicle, int action)
{
	switch (action) {
		case SEL_SELECT:
			if (!(follicle->flag & HAIR_FOLLICLE_SELECT)) {
				follicle->flag |= HAIR_FOLLICLE_SELECT;
			}
			break;
		case SEL_DESELECT:
			if (follicle->flag & HAIR_FOLLICLE_SELECT) {
				follicle->flag &= ~HAIR_FOLLICLE_SELECT;
			}
			break;
		case SEL_INVERT:
			if (!(follicle->flag & HAIR_FOLLICLE_SELECT)) {
				follicle->flag |= HAIR_FOLLICLE_SELECT;
			}
			else {
				follicle->flag &= ~HAIR_FOLLICLE_SELECT;
			}
			break;
	}
}

static int hair_select_all_exec(bContext *C, wmOperator *op)
{
	const ToolSettings *settings = CTX_data_tool_settings(C);
	Object *obedit = CTX_data_edit_object(C);;
	HairSystem *hsys = obedit->data;
	EditHair *edit = hsys->edithair;
	int action = RNA_enum_get(op->ptr, "action");

	if (action == SEL_TOGGLE) {
		switch (settings->hair_edit_settings.select_mode) {
			case HAIR_SELECT_FOLLICLES: {
				action = hair_has_selected_follicles(edit) ? SEL_DESELECT : SEL_SELECT;
				break;
			}
			case HAIR_SELECT_VERTICES: {
				BLI_assert(false);
				break;
			}
			case HAIR_SELECT_TIPS: {
				BLI_assert(false);
				break;
			}
		}
	}

	switch (settings->hair_edit_settings.select_mode) {
		case HAIR_SELECT_FOLLICLES: {
			HairFollicle *follicle;
			HairIterator iter;
			BKE_HAIR_ITER_FOLLICLES(follicle, &iter, edit->pattern) {
				hair_follicle_select_action_apply(follicle, action);
			}
			break;
		}
		case HAIR_SELECT_VERTICES: {
			BLI_assert(false);
			break;
		}
		case HAIR_SELECT_TIPS: {
			BLI_assert(false);
			break;
		}
	}

	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_SELECT);
	DEG_id_tag_update(obedit->data, DEG_TAG_SELECT_UPDATE);
	WM_event_add_notifier(C, NC_GEOM|ND_SELECT|NA_SELECTED, obedit);

	return OPERATOR_FINISHED;
}

void HAIR_OT_select_all(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "(De)select All";
	ot->idname = "HAIR_OT_select_all";
	ot->description = "(De)select all hair points";

	/* api callbacks */
	ot->exec = hair_select_all_exec;
	ot->poll = ED_operator_edithair;

	/* flags */
	ot->flag = OPTYPE_REGISTER|OPTYPE_UNDO;

	WM_operator_properties_select_all(ot);
}
