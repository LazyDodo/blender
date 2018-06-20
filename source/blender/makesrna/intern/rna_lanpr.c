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
 * Contributor(s): Blender Foundation (2008).
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/makesrna/intern/rna_lanpr.c
 *  \ingroup RNA
 */

#include <stdio.h>
#include <stdlib.h>

#include "BLI_utildefines.h"
#include "BLI_string_utils.h"

#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "rna_internal.h"

#include "DNA_lanpr_types.h"
#include "DNA_material_types.h"
#include "DNA_texture_types.h"

#include "WM_types.h"
#include "WM_api.h"

void RNA_def_lanpr(BlenderRNA *brna){

    StructRNA* srna;
	PropertyRNA* prop;

        /* line style layer */

	static const EnumPropertyItem lanpr_line_component_modes[] = {
	    {0, "NORMAL", 0, "Normal", "Normal, display all selected lines"},
        {1, "OBJECT", 0, "Object", "Display lines for selected object"},
		{2, "MATERIAL", 0, "Material", "Display lines that touches specifi material"},
        {3, "COLLECTION", 0, "Collection", "Display lines in specific collections"},
	    {0, NULL, 0, NULL, NULL}
    };

    srna = RNA_def_struct(brna, "LANPR_LineStyleComponent", NULL);
	RNA_def_struct_sdna(srna, "LANPR_LineStyleComponent");
	RNA_def_struct_ui_text(srna, "Line Style Component", "LANPR_LineStyleComponent");

//	prop = RNA_def_property(srna, "component_mode", PROP_ENUM, PROP_NONE);
//	RNA_def_property_enum_items(prop, lanpr_line_component_modes);
//	RNA_def_property_enum_default(prop, 0);
//	RNA_def_property_ui_text(prop, "Mode", "Limit the range of displayed lines");
//
	srna = RNA_def_struct(brna, "LANPR_LineStyle", NULL);
	RNA_def_struct_sdna(srna, "LANPR_LineStyle");
	RNA_def_struct_ui_text(srna, "Line Style", "LANPR_LineStyle layer");
//
//	prop = RNA_def_property(srna, "line_thickness", PROP_FLOAT, PROP_FACTOR);
//	RNA_def_property_float_default(prop, 1.0f);
//	RNA_def_property_ui_text(prop, "Thickness", "Master Thickness");
//	RNA_def_property_ui_range(prop, 0.0f, 30.0f, 0.01, 2);
//
	prop = RNA_def_property(srna, "comp", PROP_COLLECTION, PROP_NONE);
	RNA_def_property_collection_sdna(prop, NULL, "components", NULL);
	RNA_def_property_struct_type(prop, "LANPR_LineStyleComponent");
	RNA_def_property_ui_text(prop, "Components", "Line Layer Components");

}