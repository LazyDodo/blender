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
 * Contributor(s): Jörg Müller.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/makesrna/intern/rna_volume.c
 *  \ingroup RNA
 */

#include <stdlib.h>

#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "rna_internal.h"

#include "DNA_volume_types.h"

#include "BLI_math_base.h"

#ifdef RNA_RUNTIME

#include "BKE_volume.h"

static void rna_VolumeGrid_update(Main *UNUSED(bmain), Scene *UNUSED(scene), PointerRNA *UNUSED(ptr))
{
}

static void rna_Volume_update(Main *bmain, Scene *UNUSED(scene), PointerRNA *ptr)
{
	Volume *volume = ptr->data;
	BKE_volume_reload(bmain, volume);
}

static void rna_VolumeGrids_active_grid_index_range(
        PointerRNA *ptr, int *min, int *max, int *UNUSED(softmin), int *UNUSED(softmax))
{
	Volume *volume = (Volume *)ptr->data;

	*min = 0;
	*max = max_ii(0, BLI_listbase_count(&volume->grids) - 1);
}

static int rna_VolumeGrids_active_grid_index_get(PointerRNA *ptr)
{
	Volume *volume = (Volume *)ptr->data;
	return 0; // TODO
}

static void rna_VolumeGrids_active_grid_index_set(PointerRNA *ptr, int value)
{
	Volume *volume = (Volume *)ptr->data;
	// TODO
}

#else

static void rna_def_volume_grid(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;

	srna = RNA_def_struct(brna, "VolumeGrid", NULL);
	RNA_def_struct_ui_text(srna, "Volume Grid", "3D volume grid");
	RNA_def_struct_ui_icon(srna, ICON_VOLUME);

	prop = RNA_def_property(srna, "name", PROP_STRING, PROP_NONE);
	RNA_def_property_ui_text(prop, "Name", "Volume grid name");
	RNA_def_property_update(prop, 0, "rna_VolumeGrid_update");
}

static void rna_def_volume_grids(BlenderRNA *brna,  PropertyRNA *cprop)
{
	StructRNA *srna;
	PropertyRNA *prop;

	RNA_def_property_srna(cprop, "VolumeGrids");
	srna = RNA_def_struct(brna, "VolumeGrids", NULL);
	RNA_def_struct_sdna(srna, "Volume");
	RNA_def_struct_ui_text(srna, "Volume Grids", "3D volume grids");

	prop = RNA_def_property(srna, "active_index", PROP_INT, PROP_UNSIGNED);
	RNA_def_property_int_funcs(prop, "rna_VolumeGrids_active_grid_index_get",
	                           "rna_VolumeGrids_active_grid_index_set",
	                           "rna_VolumeGrids_active_grid_index_range");
	RNA_def_property_ui_text(prop, "Active Grid Index", "Index of active volume grid");
}

static void rna_def_volume(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;

	srna = RNA_def_struct(brna, "Volume", "ID");
	RNA_def_struct_ui_text(srna, "Volume", "Volume data-block for 3D volume grids");
	RNA_def_struct_ui_icon(srna, ICON_VOLUME);

	prop = RNA_def_property(srna, "filepath", PROP_STRING, PROP_FILEPATH);
	RNA_def_property_ui_text(prop, "File Path", "Volume sample file used by this Volume data-block");
	RNA_def_property_update(prop, 0, "rna_Volume_update");

	prop = RNA_def_property(srna, "packed_file", PROP_POINTER, PROP_NONE);
	RNA_def_property_pointer_sdna(prop, NULL, "packedfile");
	RNA_def_property_ui_text(prop, "Packed File", "");

	prop = RNA_def_property(srna, "grids", PROP_COLLECTION, PROP_NONE);
	RNA_def_property_struct_type(prop, "VolumeGrid");
	RNA_def_property_ui_text(prop, "Grids", "3D volume grids");
	rna_def_volume_grids(brna, prop);


	/* common */
	rna_def_animdata_common(srna);
}


void RNA_def_volume(BlenderRNA *brna)
{
	rna_def_volume_grid(brna);
	rna_def_volume(brna);
}

#endif

