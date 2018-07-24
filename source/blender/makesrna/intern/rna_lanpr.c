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

#ifdef RNA_RUNTIME



#else

void RNA_def_lanpr(BlenderRNA *brna){

    StructRNA* srna;
	PropertyRNA* prop;

        /* line style layer */

	static const EnumPropertyItem lanpr_line_component_modes[] = {
	    {0, "ALL", 0, "All", "Select All lines, lines are already selected are not affected"},
        {1, "OBJECT", 0, "Object", "Display lines for selected object"},
		{2, "MATERIAL", 0, "Material", "Display lines that touches specific material"},
        {3, "COLLECTION", 0, "Collection", "Display lines in specific collections"},
	    {0, NULL, 0, NULL, NULL}
    };

    srna = RNA_def_struct(brna, "LANPR_LineLayerComponent", NULL);
	RNA_def_struct_sdna(srna, "LANPR_LineLayerComponent");
	RNA_def_struct_ui_text(srna, "Line Layer Component", "LANPR_LineLayerComponent");

	prop = RNA_def_property(srna, "component_mode", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, lanpr_line_component_modes);
	RNA_def_property_enum_default(prop, 0);
	RNA_def_property_ui_text(prop, "Mode", "Limit the range of displayed lines");

    prop = RNA_def_property(srna, "object_select", PROP_POINTER, PROP_NONE);
    RNA_def_property_struct_type(prop, "Object");
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_ui_text(prop, "Object", "Display lines for selected object");

    prop = RNA_def_property(srna, "material_select", PROP_POINTER, PROP_NONE);
    RNA_def_property_struct_type(prop, "Material");
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_ui_text(prop, "Material", "Display lines that touches specific material");

    prop = RNA_def_property(srna, "collection_select", PROP_POINTER, PROP_NONE);
    RNA_def_property_struct_type(prop, "Collection");
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_ui_text(prop, "Collection", "Display lines in specific collections");
    


	srna = RNA_def_struct(brna, "LANPR_LineLayer", NULL);
	RNA_def_struct_sdna(srna, "LANPR_LineLayer");
	RNA_def_struct_ui_text(srna, "Line Layer", "LANPR_LineLayer");

	// removed for mow
    //prop = RNA_def_property(srna, "use_differnt_style", PROP_BOOLEAN, PROP_NONE);
	//RNA_def_property_ui_text(prop, "Different Style", "Use different line styles for differnt line types");

    prop = RNA_def_property(srna, "use_qi_range", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_ui_text(prop, "QI Range", "Use QI Range (occlusion levels) to select lines");

	prop = RNA_def_property(srna, "enable_contour", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_ui_text(prop, "Enable Contour", "Draw contour lines");

    prop = RNA_def_property(srna, "enable_crease", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_ui_text(prop, "Enable Crease", "Draw crease lines");

	prop = RNA_def_property(srna, "enable_edge_mark", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_ui_text(prop, "Enable Edge Mark", "Draw edge marks");

	prop = RNA_def_property(srna, "enable_material_seperate", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_ui_text(prop, "Enable Material Lines", "Draw material seperators");

	prop = RNA_def_property(srna, "enable_intersection", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_ui_text(prop, "Enable intersection Lines", "Draw intersection lines");

	prop = RNA_def_property(srna, "qi_begin", PROP_INT, PROP_NONE);
	RNA_def_property_int_default(prop, 0);
	RNA_def_property_ui_text(prop, "QI Begin", "QI Begin");
	RNA_def_property_range(prop, 0, 128);

	prop = RNA_def_property(srna, "qi_end", PROP_INT, PROP_NONE);
	RNA_def_property_int_default(prop, 0);
	RNA_def_property_ui_text(prop, "QI End", "QI End");
	RNA_def_property_range(prop, 0, 128);


	prop = RNA_def_property(srna, "thickness", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Thickness", "Master Thickness");
	RNA_def_property_ui_range(prop, 0.0f, 30.0f, 0.01, 2);

    prop = RNA_def_property(srna, "thickness_crease", PROP_FLOAT, PROP_FACTOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Thickness", "Crease Thickness");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.01, 2);

    prop = RNA_def_property(srna, "thickness_material", PROP_FLOAT, PROP_FACTOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Thickness", "Material Thickness");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.01, 2);

    prop = RNA_def_property(srna, "thickness_edge_mark", PROP_FLOAT, PROP_FACTOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Thickness", "Edge Mark Thickness");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.01, 2);

    prop = RNA_def_property(srna, "thickness_intersection", PROP_FLOAT, PROP_FACTOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Thickness", "Edge Mark Thickness");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.01, 2);

    prop = RNA_def_property(srna, "color", PROP_FLOAT, PROP_COLOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_array(prop, 3);
	RNA_def_property_ui_text(prop, "Color", "Master Color");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1, 2);

    prop = RNA_def_property(srna, "crease_color", PROP_FLOAT, PROP_COLOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_array(prop, 3);
	RNA_def_property_ui_text(prop, "Crease Color", "Drawing crease lines using this color");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1, 2);

	prop = RNA_def_property(srna, "material_color", PROP_FLOAT, PROP_COLOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_array(prop, 3);
	RNA_def_property_ui_text(prop, "Material Line Color", "Drawing material seperate lines using this color");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1, 2);

	prop = RNA_def_property(srna, "edge_mark_color", PROP_FLOAT, PROP_COLOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_array(prop, 3);
	RNA_def_property_ui_text(prop, "Edge Mark Color", "Drawing edge marks using this color");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1, 2);

    prop = RNA_def_property(srna, "intersection_color", PROP_FLOAT, PROP_COLOR);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_array(prop, 3);
	RNA_def_property_ui_text(prop, "Edge Mark Color", "Drawing edge marks using this color");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1, 2);

	prop = RNA_def_property(srna, "components", PROP_COLLECTION, PROP_NONE);
	RNA_def_property_collection_sdna(prop, NULL, "components", NULL);
	RNA_def_property_struct_type(prop, "LANPR_LineLayerComponent");
	RNA_def_property_ui_text(prop, "Components", "Line Layer Components");

}

#endif