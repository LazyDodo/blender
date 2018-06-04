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

/** \file blender/editors/groom/editgroom_select.c
 *  \ingroup edgroom
 */

#include "DNA_groom_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_math.h"

#include "BKE_context.h"
#include "BKE_groom.h"
#include "BKE_report.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_screen.h"
#include "ED_types.h"
#include "ED_view3d.h"
#include "ED_groom.h"

#include "groom_intern.h"

#include "RNA_access.h"
#include "RNA_define.h"

bool ED_groom_select_check_regions(const EditGroom *edit)
{
	for (GroomRegion* region = edit->regions.first; region; region = region->next)
	{
		if (region->flag & GM_REGION_SELECT)
		{
			return true;
		}
	}
	
	return false;
}

bool ED_groom_select_check_curves(const EditGroom *edit)
{
	for (GroomRegion* region = edit->regions.first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		GroomSection *section = bundle->sections;
		for (int i = 0; i < bundle->totsections; ++i, ++section)
		{
			if (section->flag & GM_SECTION_SELECT) {
				return true;
			}
		}
	}
	
	return false;
}

bool ED_groom_select_check_sections(const EditGroom *edit)
{
	for (GroomRegion* region = edit->regions.first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		GroomSectionVertex *vertex = bundle->verts;
		for (int i = 0; i < bundle->totverts; ++i, ++vertex)
		{
			if (vertex->flag & GM_VERTEX_SELECT)
			{
				return true;
			}
		}
	}
	
	return false;
}

void ED_groom_select_regions(EditGroom *edit, EditGroomSelectCb select_cb, void *userdata)
{
	for (GroomRegion* region = edit->regions.first; region; region = region->next)
	{
		const bool select = select_cb(userdata, region->flag & GM_REGION_SELECT);
		if (select)
		{
			region->flag |= GM_REGION_SELECT;
		}
		else
		{
			region->flag &= ~GM_REGION_SELECT;
		}
	}
}

void ED_groom_select_curves(EditGroom *edit, EditGroomSelectCb select_cb, void *userdata)
{
	for (GroomRegion* region = edit->regions.first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		GroomSection *section = bundle->sections;
		for (int i = 0; i < bundle->totsections; ++i, ++section)
		{
			const bool select = select_cb(userdata, section->flag & GM_SECTION_SELECT);
			if (select)
			{
				section->flag |= GM_SECTION_SELECT;
			}
			else
			{
				section->flag &= ~GM_SECTION_SELECT;
			}
		}
	}
}

void ED_groom_select_sections(EditGroom *edit, EditGroomSelectCb select_cb, void *userdata)
{
	for (GroomRegion* region = edit->regions.first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		GroomSectionVertex *vertex = bundle->verts;
		for (int i = 0; i < bundle->totverts; ++i, ++vertex)
		{
			const bool select = select_cb(userdata, vertex->flag & GM_VERTEX_SELECT);
			if (select)
			{
				vertex->flag |= GM_VERTEX_SELECT;
			}
			else
			{
				vertex->flag &= ~GM_VERTEX_SELECT;
			}
		}
	}
}

static bool groom_select_all_cb(void *UNUSED(userdata), bool UNUSED(is_selected))
{
	return true;
}

static bool groom_deselect_all_cb(void *UNUSED(userdata), bool UNUSED(is_selected))
{
	return false;
}

static bool groom_select_swap_cb(void *UNUSED(userdata), bool is_selected)
{
	return !is_selected;
}

static bool groom_has_selected(EditGroom *edit, GroomEditMode mode)
{
	switch (mode)
	{
		case GM_EDIT_MODE_REGIONS:
			return ED_groom_select_check_regions(edit);
		case GM_EDIT_MODE_CURVES:
			return ED_groom_select_check_curves(edit);
		case GM_EDIT_MODE_SECTIONS:
			return ED_groom_select_check_sections(edit);
	}
	return false;
}

static int de_select_all_exec(bContext *C, wmOperator *op)
{
	GroomEditMode mode = CTX_data_tool_settings(C)->groom_edit_settings.mode;
	Object *obedit = CTX_data_edit_object(C);
	Groom *groom = obedit->data;
	int action = RNA_enum_get(op->ptr, "action");

	EditGroomSelectCb cb;
	switch (action) {
		case SEL_SELECT:
			cb = groom_select_all_cb;
			break;
		case SEL_DESELECT:
			cb = groom_deselect_all_cb;
			break;
		case SEL_INVERT:
			cb = groom_select_swap_cb;
			break;
		case SEL_TOGGLE:
		{
			if (groom_has_selected(groom->editgroom, mode)) {
				cb = groom_deselect_all_cb;
			}
			else
			{
				cb = groom_select_all_cb;
			}
		}
	}

	switch (mode)
	{
		case GM_EDIT_MODE_REGIONS:
			ED_groom_select_regions(groom->editgroom, cb, NULL);
			break;
		case GM_EDIT_MODE_CURVES:
			ED_groom_select_curves(groom->editgroom, cb, NULL);
			break;
		case GM_EDIT_MODE_SECTIONS:
			ED_groom_select_sections(groom->editgroom, cb, NULL);
			break;
	}

	WM_event_add_notifier(C, NC_GEOM | ND_SELECT, obedit->data);

	return OPERATOR_FINISHED;
}

void GROOM_OT_select_all(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "(De)select All";
	ot->idname = "GROOM_OT_select_all";
	ot->description = "(De)select all control points";

	/* api callbacks */
	ot->exec = de_select_all_exec;
	ot->poll = ED_operator_editgroom;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	/* properties */
	WM_operator_properties_select_all(ot);
}

/****************************** Mouse Selection *************************/

static void select_pick_findnearest_cb(
        void *userdata,
        GroomRegion *region,
        GroomSection *section,
        GroomSectionVertex *vertex,
        const float screen_co[2])
{
	struct
	{
		GroomRegion *region;
		GroomSection *section;
		GroomSectionVertex *vertex;
		float dist;
		bool select;
		float mval_fl[2];
	} *data = userdata;

	float dist_test = len_manhattan_v2v2(data->mval_fl, screen_co);
	
	/* bias towards unselected items */
	if (data->select &&
	    ((vertex && vertex->flag & GM_VERTEX_SELECT) ||
	     (section && section->flag & GM_SECTION_SELECT) ||
	     (region && region->flag & GM_REGION_SELECT)))
	{
		dist_test += 5.0f;
	}

	if (dist_test < data->dist) {
		data->dist = dist_test;
		data->region = region;
		data->section = section;
		data->vertex = vertex;
	}
}

static void groom_set_region_select_flags(Groom *groom, int flag)
{
	for (GroomRegion* region = groom->editgroom->regions.first; region; region = region->next)
	{
		region->flag = (region->flag & ~GM_REGION_SELECT) | (flag & GM_REGION_SELECT);
	}
}

static void groom_set_curve_select_flags(Groom *groom, int flag)
{
	for (GroomRegion* region = groom->editgroom->regions.first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		GroomSection *section = bundle->sections;
		for (int i = 0; i < bundle->totsections; ++i, ++section)
		{
			section->flag = (section->flag & ~GM_SECTION_SELECT) | (flag & GM_SECTION_SELECT);
		}
	}
}

static void groom_set_section_select_flags(Groom *groom, int flag)
{
	for (GroomRegion* region = groom->editgroom->regions.first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		GroomSectionVertex *vertex = bundle->verts;
		for (int i = 0; i < bundle->totverts; ++i, ++vertex)
		{
			vertex->flag = (vertex->flag & ~GM_VERTEX_SELECT) | (flag & GM_VERTEX_SELECT);
		}
	}
}

bool ED_groom_select_pick(bContext *C, const int mval[2], bool extend, bool deselect, bool toggle)
{
	ViewContext vc;
	ED_view3d_viewcontext_init(C, &vc);
	Groom *groom = vc.obedit->data;

	struct
	{
		GroomRegion *region;
		GroomSection *section;
		GroomSectionVertex *vertex;
		float dist;
		bool select;
		float mval_fl[2];
	} data = {NULL};

	data.dist = ED_view3d_select_dist_px();
	data.select = true;
	data.mval_fl[0] = mval[0];
	data.mval_fl[1] = mval[1];

	ED_view3d_init_mats_rv3d(vc.obedit, vc.rv3d);
	groom_foreachScreenVert(&vc, select_pick_findnearest_cb, &data, V3D_PROJ_TEST_CLIP_DEFAULT);

	bool found = false;
	if (data.vertex)
	{
		if (extend)
		{
			data.vertex->flag |= GM_VERTEX_SELECT;
		}
		else if (deselect)
		{
			data.vertex->flag &= ~GM_VERTEX_SELECT;
		}
		else if (toggle)
		{
			data.vertex->flag ^= GM_VERTEX_SELECT;
		}
		else
		{
			/* deselect all other verts */
			groom_set_section_select_flags(groom, 0);
			data.vertex->flag |= GM_VERTEX_SELECT;
		}
		
		if (data.vertex->flag & GM_VERTEX_SELECT)
		{
			/* set active section */
			groom_set_region_select_flags(groom, 0);
			groom_set_curve_select_flags(groom, 0);
			data.section->flag |= GM_SECTION_SELECT;
			data.region->flag |= GM_REGION_SELECT;
		}
		
		found = true;
	}
	else if (data.section)
	{
		if (extend)
		{
			data.section->flag |= GM_SECTION_SELECT;
		}
		else if (deselect)
		{
			data.section->flag &= ~GM_SECTION_SELECT;
		}
		else if (toggle)
		{
			data.section->flag ^= GM_SECTION_SELECT;
		}
		else
		{
			/* deselect all other sections */
			groom_set_curve_select_flags(groom, 0);
			data.section->flag |= GM_SECTION_SELECT;
		}
		
		if (data.section->flag & GM_SECTION_SELECT)
		{
			/* set active region */
			groom_set_region_select_flags(groom, 0);
			data.region->flag |= GM_REGION_SELECT;
		}
		
		found = true;
	}
	else if (data.region)
	{
		if (extend)
		{
			data.region->flag |= GM_REGION_SELECT;
		}
		else if (deselect)
		{
			data.region->flag &= ~GM_REGION_SELECT;
		}
		else if (toggle)
		{
			data.region->flag ^= GM_REGION_SELECT;
		}
		else
		{
			/* deselect all other regions */
			groom_set_region_select_flags(groom, 0);
			data.region->flag |= GM_REGION_SELECT;
		}
		
		found = true;
	}

	if (found)
	{
		WM_event_add_notifier(C, NC_GEOM | ND_SELECT, vc.obedit->data);
		return true;
	}

	return false;
}
