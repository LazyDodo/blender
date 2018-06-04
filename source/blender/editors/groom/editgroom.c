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

/** \file blender/editors/groom/editgroom.c
 *  \ingroup edgroom
 */

#include "DNA_groom_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_array_utils.h"
#include "BLI_blenlib.h"

#include "BKE_context.h"
#include "BKE_global.h"
#include "BKE_groom.h"
#include "BKE_library.h"
#include "BKE_main.h"
#include "BKE_report.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_groom.h"
#include "ED_object.h"
#include "ED_screen.h"
#include "ED_types.h"
#include "ED_util.h"

#include "groom_intern.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

/********************** Load/Make/Free ********************/

static void groom_regions_free(ListBase *regions)
{
	for (GroomRegion *region = regions->first; region; region = region->next)
	{
		if (region->scalp_samples)
		{
			MEM_freeN(region->scalp_samples);
		}
		
		GroomBundle *bundle = &region->bundle;
		BKE_groom_bundle_curve_cache_clear(bundle);
		
		if (bundle->verts)
		{
			MEM_freeN(bundle->verts);
		}
		if (bundle->sections)
		{
			MEM_freeN(bundle->sections);
		}
		if (bundle->guides)
		{
			MEM_freeN(bundle->guides);
		}
		if (bundle->guide_shape_weights)
		{
			MEM_freeN(bundle->guide_shape_weights);
		}
	}
	BLI_freelistN(regions);
}

static void groom_regions_copy(ListBase *regions_dst, ListBase *regions_src)
{
	BLI_duplicatelist(regions_dst, regions_src);
	for (GroomRegion *region = regions_dst->first; region; region = region->next)
	{
		if (region->scalp_samples)
		{
			region->scalp_samples = MEM_dupallocN(region->scalp_samples);
		}
		
		GroomBundle *bundle = &region->bundle;
		if (bundle->curvecache)
		{
			bundle->curvecache = MEM_dupallocN(bundle->curvecache);
		}
		if (bundle->sections)
		{
			bundle->sections = MEM_dupallocN(bundle->sections);
		}
		if (bundle->verts)
		{
			bundle->verts = MEM_dupallocN(bundle->verts);
		}
		if (bundle->guides)
		{
			bundle->guides = MEM_dupallocN(bundle->guides);
		}
		if (bundle->guide_shape_weights)
		{
			bundle->guide_shape_weights = MEM_dupallocN(bundle->guide_shape_weights);
		}
	}
}

void ED_groom_editgroom_make(Object *obedit)
{
	Groom *groom = obedit->data;

	ED_groom_editgroom_free(obedit);

	groom->editgroom = MEM_callocN(sizeof(EditGroom), "editgroom");
	groom_regions_copy(&groom->editgroom->regions, &groom->regions);
}

void ED_groom_editgroom_load(Object *obedit)
{
	Groom *groom = obedit->data;
	
	groom_regions_free(&groom->regions);
	groom_regions_copy(&groom->regions, &groom->editgroom->regions);
}

void ED_groom_editgroom_free(Object *ob)
{
	Groom *groom = ob->data;
	
	if (groom->editgroom) {
		groom_regions_free(&groom->editgroom->regions);
		
		MEM_freeN(groom->editgroom);
		groom->editgroom = NULL;
	}
}
