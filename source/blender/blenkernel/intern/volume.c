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
 * Contributor(s): Brecht Van Lommel.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/volume.c
 *  \ingroup bke
 */

#include "MEM_guardedalloc.h"

#include "DNA_object_types.h"
#include "DNA_sound_types.h"
#include "DNA_volume_types.h"

#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_string.h"
#include "BLI_utildefines.h"

#include "BKE_animsys.h"
#include "BKE_global.h"
#include "BKE_library.h"
#include "BKE_library_query.h"
#include "BKE_library_remap.h"
#include "BKE_main.h"
#include "BKE_packedFile.h"
#include "BKE_volume.h"

#ifdef WITH_OPENVDB
#include "openvdb_capi.h"
#endif

/* OpenVDB integration */

static void volume_create_openvdb(Volume *volume)
{
	if (!volume->filepath[0]) {
		return;
	}

#ifdef WITH_OPENVDB
	OpenVDBReader *reader = OpenVDBReader_create();
	OpenVDBReader_open(reader, volume->filepath);

	size_t num_grids = OpenVDBReader_num_grids(reader);

	for (int i = 0; i < num_grids; i++) {
		VolumeGrid *grid = MEM_callocN(sizeof(VolumeGrid), "VolumeGrid");
		const char *name = OpenVDBReader_grid_name(reader, i);
		BLI_strncpy(grid->name, name, sizeof(grid->name));
		BLI_addtail(&volume->grids, grid);
	}

	OpenVDBReader_free(reader);
#endif
}

static void volume_free_openvdb(Volume *volume)
{
#ifdef WITH_OPENVDB
	BLI_freelistN(&volume->grids);
#endif
}

/* Volume datablock */

void BKE_volume_init(Volume *volume)
{
	BLI_assert(MEMCMP_STRUCT_OFS_IS_ZERO(volume, id));

	volume->filepath[0] = '\0';
	volume->packedfile = NULL;
	BLI_listbase_clear(&volume->grids);
	volume->flag = 0;
}

void *BKE_volume_add(Main *bmain, const char *name)
{
	Volume *volume;

	volume = BKE_libblock_alloc(bmain, ID_VO, name, 0);

	BKE_volume_init(volume);

	return volume;
}

void BKE_volume_copy_data(Main *bmain, Volume *volume_dst, const Volume *UNUSED(volume_src), const int UNUSED(flag))
{
	if (volume_dst->packedfile) {
		volume_dst->packedfile = dupPackedFile(volume_dst->packedfile);
	}

	BLI_listbase_clear(&volume_dst->grids);
	BKE_volume_reload(bmain, volume_dst);
}

Volume *BKE_volume_copy(Main *bmain, const Volume *volume)
{
	Volume *volume_copy;
	BKE_id_copy_ex(bmain, &volume->id, (ID **)&volume_copy, 0, false);
	return volume_copy;
}

void BKE_volume_make_local(Main *bmain, Volume *volume, const bool lib_local)
{
	BKE_id_make_local_generic(bmain, &volume->id, true, lib_local);
}

void BKE_volume_free(Volume *volume)
{
	BKE_animdata_free((ID *)volume, false);
	volume_free_openvdb(volume);
}

void BKE_volume_reload(Main *UNUSED(bmain), Volume *volume)
{
	volume_free_openvdb(volume);
	volume_create_openvdb(volume);
}
