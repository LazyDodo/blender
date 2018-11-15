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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Copyright (C) 2018 by Martin Felke.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */


#include "MEM_guardedalloc.h"
#include "BKE_fracture.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_pointcache.h"

#include "BLI_listbase.h"
#include "BLI_string.h"

#include "DNA_fracture_types.h"
#include "DNA_object_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_scene_types.h"
#include "DNA_rigidbody_types.h"

#include "DEG_depsgraph_query.h"

#include "limits.h"
#include "PIL_time.h"

void BKE_fracture_dynamic_do(FractureModifierData *fmd, Object* ob, Scene* scene, Depsgraph* depsgraph, Main* bmain)
{
	int count = 0;
	int totpoint = 0;
	FractureID *fid = fmd->shared->fracture_ids.first;

	while(fid) {
		FracPointCloud points;
		if (!fid->mi->fractured) {
			points = BKE_fracture_points_get(depsgraph, fmd, ob, fid->mi);
			totpoint += points.totpoints;
			MEM_freeN(points.points);
			points.totpoints = 0;
		}
		fid = fid->next;
	}

	fid = (FractureID*)fmd->shared->fracture_ids.first;
	while(fid){
		if (!fid->mi->fractured) {
			BKE_fracture_do(fmd, fid->mi, ob, depsgraph, bmain, scene, false);
		}
		fid->mi->fractured = true;
		count++;

		BLI_remlink(&fmd->shared->fracture_ids, fid);
		fid = (FractureID*)fmd->shared->fracture_ids.first;
	}

	fmd->shared->fracture_ids.first = NULL;
	fmd->shared->fracture_ids.last = NULL;
}
