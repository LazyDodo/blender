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

/** \file blender/editors/groom/editgroom_region.c
 *  \ingroup edgroom
 */

#include "DNA_groom_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_math.h"

#include "BLT_translation.h"

#include "BKE_context.h"
#include "BKE_groom.h"
#include "BKE_library.h"
#include "BKE_object_facemap.h"

#include "DEG_depsgraph.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_object.h"
#include "ED_screen.h"
#include "ED_util.h"
#include "ED_view3d.h"
#include "ED_groom.h"

#include "UI_resources.h"
#include "UI_interface.h"

#include "groom_intern.h"

/* GROOM_OT_region_add */

static void region_add_set_bundle_curve(GroomRegion *region, const float loc[3], const float rot[3][3], float length)
{
	GroomBundle *bundle = &region->bundle;
	
	bundle->totsections = 2;
	bundle->sections = MEM_callocN(sizeof(GroomSection) * bundle->totsections, "groom bundle sections");
	
	madd_v3_v3v3fl(bundle->sections[0].center, loc, rot[2], 0.0f);
	madd_v3_v3v3fl(bundle->sections[1].center, loc, rot[2], length);
}

static int region_add_poll(bContext *C)
{
	if (!ED_groom_object_poll(C))
	{
		return false;
	}
	
	/* We want a scalp object to make this useful */
	Object *ob = ED_object_context(C);
	Groom *groom = ob->data;
	return groom->scalp_object != NULL;
}

static int region_add_exec(bContext *C, wmOperator *op)
{
	const Depsgraph *depsgraph = CTX_data_depsgraph(C);
	Object *ob = ED_object_context(C);
	Groom *groom = ob->data;
	char scalp_facemap_name[MAX_VGROUP_NAME];
	RNA_string_get(op->ptr, "scalp_facemap", scalp_facemap_name);
	if (scalp_facemap_name[0] == '\0' ||
	    !BKE_object_facemap_find_name(groom->scalp_object, scalp_facemap_name))
	{
		return OPERATOR_CANCELLED;
	}

	WM_operator_view3d_unit_defaults(C, op);

	float loc[3], rot[3];
	unsigned int layer;
	if (!ED_object_add_generic_get_opts(C, op, 'Z', loc, rot, NULL, &layer, NULL))
		return OPERATOR_CANCELLED;

	float mat[4][4];
	ED_object_new_primitive_matrix(C, ob, loc, rot, mat);

	GroomRegion *region = MEM_callocN(sizeof(GroomRegion), "groom region");
	ListBase *regions = (groom->editgroom ? &groom->editgroom->regions : &groom->regions);
	BLI_addtail(regions, region);

	float scalp_loc[3];
	float scalp_rot[3][3];
	zero_v3(scalp_loc);
	unit_m3(scalp_rot);
	
	if (BKE_groom_set_region_scalp_facemap(groom, region, scalp_facemap_name))
	{
		const struct Mesh *scalp = BKE_groom_get_scalp(depsgraph, groom);
		BLI_assert(scalp != NULL);
		
		if (BKE_groom_region_bind(depsgraph, groom, region, true))
		{
			BKE_groom_calc_region_transform_on_scalp(region, scalp, scalp_loc, scalp_rot);
		}
	}
	
	region_add_set_bundle_curve(region, scalp_loc, scalp_rot, 1.0f);
	BKE_groom_region_reset_shape(depsgraph, groom, region);
	
	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);
	DEG_id_tag_update(&ob->id, OB_RECALC_DATA);

	return OPERATOR_FINISHED;
}

static void region_add_draw(bContext *C, wmOperator *op)
{
	uiLayout *layout = op->layout;
	Object *ob = ED_object_context(C);
	Groom *groom = ob->data;
	PointerRNA scalp_ob_ptr;
	RNA_id_pointer_create(&groom->scalp_object->id, &scalp_ob_ptr);

	if (groom->scalp_object)
	{
		uiItemPointerR(layout, op->ptr, "scalp_facemap", &scalp_ob_ptr, "face_maps", NULL, ICON_NONE);
	}
	else
	{
		uiItemR(layout, op->ptr, "scalp_facemap", 0, NULL, ICON_NONE);
	}
}

void GROOM_OT_region_add(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Add Region";
	ot->description = "Add a new region to the groom object";
	ot->idname = "GROOM_OT_region_add";

	/* api callbacks */
	ot->exec = region_add_exec;
	ot->poll = region_add_poll;
	ot->invoke = WM_operator_props_popup_confirm;
	ot->ui = region_add_draw;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	ED_object_add_generic_props(ot, false);
	RNA_def_string(ot->srna, "scalp_facemap", NULL, MAX_VGROUP_NAME, "Scalp Facemap Name", "Facemap to which to bind the new region");
}

/* GROOM_OT_region_bind */

static int region_bind_poll(bContext *C)
{
	if (!ED_operator_scene_editable(C))
	{
		return 0;
	}
	
	Object *ob = ED_object_context(C);
	Groom *groom = ob->data;
	if (groom->editgroom)
	{
		return 0;
	}
	
	return 1;
}

static int region_bind_exec(bContext *C, wmOperator *op)
{
	const Depsgraph *depsgraph = CTX_data_depsgraph(C);
	Object *ob = ED_object_context(C);
	Groom *groom = ob->data;
	const bool force_rebind = RNA_boolean_get(op->ptr, "force_rebind");

	GroomRegion *region = CTX_data_pointer_get_type(C, "groom_region", &RNA_GroomRegion).data;
	if (!region)
	{
		region = BLI_findlink(&groom->regions, groom->active_region);
		if (!region)
		{
			return OPERATOR_CANCELLED;
		}
	}

	BKE_groom_region_bind(depsgraph, groom, region, force_rebind);

	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);
	DEG_id_tag_update(&ob->id, OB_RECALC_DATA);

	return OPERATOR_FINISHED;
}

void GROOM_OT_region_bind(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Bind Region";
	ot->description = "Bind a groom bundle to its scalp region";
	ot->idname = "GROOM_OT_region_bind";

	/* api callbacks */
	ot->exec = region_bind_exec;
	ot->poll = region_bind_poll;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	RNA_def_boolean(ot->srna, "force_rebind", true, "Force Rebind",
	                "Force rebinding of the groom region even if a binding already exists");
}
