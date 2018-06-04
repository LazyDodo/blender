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

/** \file blender/makesrna/intern/rna_groom.c
 *  \ingroup RNA
 */

#include <stdlib.h>

#include "DNA_groom_types.h"
#include "DNA_scene_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_string_utils.h"
#include "BLI_utildefines.h"

#include "BLT_translation.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "rna_internal.h"

#include "WM_types.h"
#include "DNA_object_types.h"


#ifdef RNA_RUNTIME

#include "BLI_listbase.h"
#include "BLI_math.h"

#include "WM_api.h"

#include "BKE_groom.h"
#include "BKE_material.h"
#include "BKE_object_facemap.h"

#include "DEG_depsgraph.h"

static void UNUSED_FUNCTION(rna_Groom_update)(Main *UNUSED(bmain), Scene *UNUSED(scene), PointerRNA *UNUSED(ptr))
{
	WM_main_add_notifier(NC_GROOM | NA_EDITED, NULL);
}

static void rna_Groom_update_data(Main *UNUSED(bmain), Scene *UNUSED(scene), PointerRNA *ptr)
{
	DEG_id_tag_update(ptr->id.data, 0);
	WM_main_add_notifier(NC_GROOM | ND_DATA, ptr->id.data);
}

static int rna_GroomBundle_is_bound_get(PointerRNA *ptr)
{
	GroomBundle *bundle = (GroomBundle *)ptr->data;
	return (bundle->scalp_region != NULL);
}

static void rna_GroomBundle_scalp_facemap_name_set(PointerRNA *ptr, const char *value)
{
	Groom *groom = (Groom *)ptr->id.data;
	GroomBundle *bundle = (GroomBundle *)ptr->data;
	
	if (groom->scalp_object)
	{
		bFaceMap *fm = BKE_object_facemap_find_name(groom->scalp_object, value);
		if (fm) {
			/* no need for BLI_strncpy_utf8, since this matches an existing facemap */
			BLI_strncpy(bundle->scalp_facemap_name, value, sizeof(bundle->scalp_facemap_name));
			/* Bind to the region right away */
			BKE_groom_bundle_bind(groom, bundle, true);
			return;
		}
	}
	
	bundle->scalp_facemap_name[0] = '\0';
	/* Unbind from region */
	BKE_groom_bundle_unbind(bundle);
}

static PointerRNA rna_Groom_active_region_get(PointerRNA *ptr)
{
	Groom *groom = (Groom *)ptr->id.data;
	PointerRNA r_ptr;
	RNA_pointer_create(&groom->id, &RNA_GroomRegion, BLI_findlink(&groom->regions, groom->active_region), &r_ptr);
	return r_ptr;
}

static int rna_Groom_active_region_index_get(PointerRNA *ptr)
{
	Groom *groom = (Groom *)ptr->id.data;
	return groom->active_region;
}

static void rna_Groom_active_region_index_set(PointerRNA *ptr, int value)
{
	Groom *groom = (Groom *)ptr->id.data;
	groom->active_region = value;
}

static void rna_Groom_active_region_index_range(
        PointerRNA *ptr,
        int *min,
        int *max,
        int *UNUSED(softmin),
        int *UNUSED(softmax))
{
	Groom *groom = (Groom *)ptr->id.data;
	*min = 0;
	*max = max_ii(0, BLI_listbase_count(&groom->regions) - 1);
}

static const EnumPropertyItem *rna_Groom_material_slot_itemf(
        bContext *C,
        PointerRNA *UNUSED(ptr),
        PropertyRNA *UNUSED(prop),
        bool *r_free)
{
	Object *ob = CTX_data_pointer_get(C, "object").data;
	Material *ma;
	EnumPropertyItem *item = NULL;
	EnumPropertyItem tmp = {0, "", 0, "", ""};
	int totitem = 0;
	int i;

	if (ob && ob->totcol > 0) {
		for (i = 1; i <= ob->totcol; i++) {
			ma = give_current_material(ob, i);
			tmp.value = i;
			tmp.icon = ICON_MATERIAL_DATA;
			if (ma) {
				tmp.name = ma->id.name + 2;
				tmp.identifier = tmp.name;
			}
			else {
				tmp.name = "Default Material";
				tmp.identifier = tmp.name;
			}
			RNA_enum_item_add(&item, &totitem, &tmp);
		}
	}
	else {
		tmp.value = 1;
		tmp.icon = ICON_MATERIAL_DATA;
		tmp.name = "Default Material";
		tmp.identifier = tmp.name;
		RNA_enum_item_add(&item, &totitem, &tmp);
	}

	RNA_enum_item_end(&item, &totitem);
	*r_free = true;

	return item;
}

#else

static void rna_def_groom_bundle(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;
	
	srna = RNA_def_struct(brna, "GroomBundle", NULL);
	RNA_def_struct_sdna(srna, "GroomBundle");
	RNA_def_struct_ui_text(srna, "Groom Bundle", "Bundle of hair originating from a scalp region");
	
	prop = RNA_def_property(srna, "is_bound", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_funcs(prop, "rna_GroomBundle_is_bound_get", NULL);
	RNA_def_property_clear_flag(prop, PROP_EDITABLE);
	RNA_def_property_ui_text(prop, "Bound", "Bundle was successfully bound to a scalp region");
	RNA_def_property_update(prop, NC_GROOM | ND_DRAW, NULL);
	
	prop = RNA_def_property(srna, "scalp_facemap", PROP_STRING, PROP_NONE);
	RNA_def_property_string_sdna(prop, NULL, "scalp_facemap_name");
	RNA_def_property_flag(prop, PROP_NEVER_UNLINK);
	RNA_def_property_ui_text(prop, "Scalp Vertex Group", "Face map name of the scalp region");
	RNA_def_property_string_funcs(prop, NULL, NULL, "rna_GroomBundle_scalp_facemap_name_set");
	RNA_def_property_update(prop, NC_GROOM | ND_DRAW, "rna_Groom_update_data");
	
	prop = RNA_def_property(srna, "guides_count", PROP_INT, PROP_NONE);
	RNA_def_property_range(prop, 1, INT_MAX);
	RNA_def_property_ui_range(prop, 1, 1000, 1, -1);
	RNA_def_property_ui_text(prop, "Guides Count", "Number of hair guides to generate in this region");
	RNA_def_property_update(prop, NC_GROOM | ND_DRAW, "rna_Groom_update_data");
}

static void rna_def_groom_region(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;
	
	srna = RNA_def_struct(brna, "GroomRegion", NULL);
	RNA_def_struct_sdna(srna, "GroomRegion");
	RNA_def_struct_ui_text(srna, "Groom Region", "Region on the scalp");
	
	prop = RNA_def_property(srna, "bundle", PROP_POINTER, PROP_NONE);
	RNA_def_property_pointer_sdna(prop, NULL, "bundle");
	RNA_def_property_struct_type(prop, "GroomBundle");
	RNA_def_property_clear_flag(prop, PROP_EDITABLE);
	RNA_def_property_ui_text(prop, "Bundle", "Bundle geometry");
}

/* groom.regions */
static void rna_def_groom_regions(BlenderRNA *brna, PropertyRNA *cprop)
{
	StructRNA *srna;
	PropertyRNA *prop;
	
	RNA_def_property_srna(cprop, "GroomRegions");
	srna = RNA_def_struct(brna, "GroomRegions", NULL);
	RNA_def_struct_sdna(srna, "Groom");
	RNA_def_struct_ui_text(srna, "Groom Regions", "Collection of groom regions");
	
	prop = RNA_def_property(srna, "active", PROP_POINTER, PROP_NONE);
	RNA_def_property_struct_type(prop, "GroomRegion");
	RNA_def_property_pointer_funcs(prop, "rna_Groom_active_region_get", NULL, NULL, NULL);
	RNA_def_property_ui_text(prop, "Active Groom Region", "Active groom region being displayed");
	RNA_def_property_update(prop, NC_GROOM | ND_DRAW, NULL);
	
	prop = RNA_def_property(srna, "active_index", PROP_INT, PROP_UNSIGNED);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_int_funcs(prop, "rna_Groom_active_region_index_get",
	                           "rna_Groom_active_region_index_set",
	                           "rna_Groom_active_region_index_range");
	RNA_def_property_ui_text(prop, "Active Groom Region Index", "Index of active groom region");
	RNA_def_property_update(prop, NC_GROOM | ND_DRAW, NULL);
}

static void rna_def_groom(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;
	
	static const EnumPropertyItem material_slot_items[] = {
		{0, "DUMMY", 0, "Dummy", ""},
		{0, NULL, 0, NULL, NULL}
	};
	
	srna = RNA_def_struct(brna, "Groom", "ID");
	RNA_def_struct_sdna(srna, "Groom");
	RNA_def_struct_ui_text(srna, "Groom", "Guide curve geometry for hair");
	RNA_def_struct_ui_icon(srna, ICON_NONE);
	
	/* Animation Data */
	rna_def_animdata_common(srna);
	
	prop = RNA_def_property(srna, "regions", PROP_COLLECTION, PROP_NONE);
	RNA_def_property_collection_sdna(prop, NULL, "regions", NULL);
	RNA_def_property_struct_type(prop, "GroomRegion");
	RNA_def_property_ui_text(prop, "Regions", "Regions of hair");
	rna_def_groom_regions(brna, prop);
	
	prop = RNA_def_property(srna, "curve_resolution", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "curve_res");
	RNA_def_property_range(prop, 1, 1024);
	RNA_def_property_ui_range(prop, 1, 64, 1, -1);
	RNA_def_property_ui_text(prop, "Curve Resolution", "Curve subdivisions per segment");
	RNA_def_property_update(prop, 0, "rna_Groom_update_data");
	
	prop = RNA_def_property(srna, "hair_system", PROP_POINTER, PROP_NONE);
	RNA_def_property_ui_text(prop, "Hair", "Hair data");
	RNA_def_property_clear_flag(prop, PROP_EDITABLE);
	
	prop = RNA_def_property(srna, "hair_draw_settings", PROP_POINTER, PROP_NONE);
	RNA_def_property_ui_text(prop, "Hair Draw Settings", "Hair draw settings");
	RNA_def_property_clear_flag(prop, PROP_EDITABLE);
	
	prop = RNA_def_property(srna, "scalp_object", PROP_POINTER, PROP_NONE);
	RNA_def_property_pointer_sdna(prop, NULL, "scalp_object");
	RNA_def_property_ui_text(prop, "Scalp Object", "Surface for attaching hairs");
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_update(prop, 0, "rna_Groom_update_data");
	
	prop = RNA_def_property(srna, "material_index", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "material_index");
	RNA_def_property_range(prop, 1, 32767);
	RNA_def_property_ui_text(prop, "Material Index", "Index of material slot used for rendering hair fibers");
	RNA_def_property_update(prop, 0, "rna_Groom_update_data");

	prop = RNA_def_property(srna, "material_slot", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "material_index");
	RNA_def_property_enum_items(prop, material_slot_items);
	RNA_def_property_enum_funcs(prop, NULL, NULL, "rna_Groom_material_slot_itemf");
	RNA_def_property_ui_text(prop, "Material Slot", "Material slot used for rendering particles");
	RNA_def_property_update(prop, 0, "rna_Groom_update_data");
	
	UNUSED_VARS(prop);
}

void RNA_def_groom(BlenderRNA *brna)
{
	rna_def_groom(brna);
	rna_def_groom_region(brna);
	rna_def_groom_bundle(brna);
}

#endif
