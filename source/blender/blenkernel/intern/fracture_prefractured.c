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

#include "BLI_listbase.h"
#include "BLI_math.h"

#include "BKE_fracture.h"
#include "BKE_pointcache.h"
#include "BKE_rigidbody.h"
#include "BKE_global.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_main.h"

#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_fracture_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

#include "DEG_depsgraph_query.h" //must be after DNA_object_types.h, else it fails to compile

#include "PIL_time.h"

#if 0
Mesh *BKE_fracture_prefractured_apply(FractureModifierData *fmd, Object *ob, Mesh *derivedData, Depsgraph* depsgraph)
{
    bool do_refresh = (fmd->auto_execute) || (fmd->dm_group && fmd->use_constraint_group && fmd->shared->refresh_constraints);
    Scene *scene = fmd->scene;

    Mesh *final_dm = derivedData;
  //  Mesh *group_dm = BKE_fracture_group_dm(fmd, derivedData, ob, do_refresh || fmd->refresh);

    if (do_refresh) {
        fmd->shared->refresh = true;
    }

    if (fmd->shared->refresh)
    {
        BKE_fracture_refresh(fmd, ob, derivedData, depsgraph);
    }

    /* TODO_5, get rid of fmd->dm and perhaps of fmd->visible_mesh (BMESH!) too, the latter should be runtime data for creating islands ONLY */
    /* we should ideally only have one cached derivedmesh */
    if (fmd->shared->dm && fmd->shared->frac_mesh && (fmd->shared->dm->totpoly > 0)) {
        final_dm = BKE_fracture_prefractured_do(fmd, ob, fmd->shared->dm, derivedData, NULL, 0, scene, depsgraph);
    }
    else {
        final_dm = BKE_fracture_prefractured_do(fmd, ob, derivedData, derivedData, NULL, 0, scene, depsgraph);
    }

    return final_dm;
}
#endif

MeshIsland *BKE_fracture_mesh_island_create(Mesh* me, Main* bmain, Scene *scene, Object *ob)
{
	MeshIsland *mi = MEM_callocN(sizeof(MeshIsland), "mesh_island");
	mi->mesh = me;

	INIT_MINMAX(mi->min, mi->max);
	BKE_mesh_minmax(mi->mesh, mi->min, mi->max);
	BKE_fracture_mesh_center_centroid_area(mi->mesh, mi->centroid);
	unit_qt(mi->rot);

	mi->rigidbody = BKE_rigidbody_create_shard(bmain, scene, ob, NULL, mi);
	mi->rigidbody->type = RBO_TYPE_ACTIVE;
	mi->rigidbody->mesh_island_index = 0; // set when adding !!!!
	BKE_rigidbody_calc_shard_mass(ob, mi);

	return mi;
}


Mesh* BKE_fracture_apply(FractureModifierData *fmd, Object *ob, Mesh *me_orig, Depsgraph* depsgraph)
{
	Scene *scene = DEG_get_evaluated_scene(depsgraph);

	Main* bmain = G.main;
	Mesh* me_assembled = NULL;
	Mesh *me_final = NULL;
	Mesh *me = me_orig; //BKE_fracture_mesh_copy(me_orig, ob);
	//Object* ob = DEG_get_evaluated_object(depsgraph, obj);

	if (fmd->shared->refresh)
	{
		MeshIsland *mi = NULL;
		Mesh *me_tmp = NULL;

		// HACK
		fmd->shared->scene = scene;

		/*free old stuff here */
		BKE_fracture_constraints_free(fmd, scene);
		BKE_fracture_dynamic_free(fmd, scene);

		/* if refresh and if having packdata, assemble an inputmesh here (override basemesh) */
		if (fmd->shared->pack_storage.first)
		{
			//treat this as basemesh ?
			me_tmp = BKE_fracture_assemble_mesh_from_islands(fmd, &fmd->shared->pack_storage, ob);
		}
		else {
			me_tmp = me;
		}

		mi = BKE_fracture_mesh_island_create(me_tmp, bmain, scene, ob);
		mi->id = 0;
		mi->rigidbody->mesh_island_index = 0;
		BLI_addtail(&fmd->shared->mesh_islands, mi);

		/*if refresh, perform fracture */
		BKE_fracture_do(fmd, mi, ob, depsgraph, bmain, scene);

		/* refresh means full reset for dynamic too, so rebuild here one state*/
		BKE_fracture_dynamic_new_entries_add(fmd, scene, false);

		fmd->shared->refresh_constraints = true;
		fmd->shared->refresh_autohide = true;
	}
	else if (fmd->shared->refresh_dynamic)
	{
		/*if dynamic event, push state to fracture sequence*/
		BKE_fracture_dynamic_do(fmd, ob, scene, depsgraph, bmain);
	}

	/* assemble mesh from transformed meshislands */
	if (fmd->shared->mesh_islands.first)
	{
		me_assembled = BKE_fracture_assemble_mesh_from_islands(fmd, &fmd->shared->mesh_islands, ob);
		//DEG_id_tag_update(&ob->id, DEG_TAG_COPY_ON_WRITE);
	}
	else {
		me_assembled = me;
	}

	fmd->shared->mesh_cached = me_assembled;

	/*if refresh constraints, build constraints */
	if (fmd->shared->refresh_constraints)
	{
		if (!fmd->shared->refresh)
			// do not free twice
			BKE_fracture_constraints_free(fmd, scene);

		BKE_fracture_external_constraints_setup(fmd, scene, ob);

		//if (!fmd->shared->refresh)
			/* update scene here in case only the constraints updated, dont update twice*/
		//	DEG_id_tag_update(&ob->id, DEG_TAG_COPY_ON_WRITE);
	}

	fmd->shared->refresh = false;

	/*if refresh autohide, initialize its structures */
	if (fmd->shared->refresh_autohide)
	{
		fmd->distortion_cached = false;
		BKE_fracture_autohide_refresh(fmd, ob, me_assembled);

		if (fmd->autohide_dist > 0 && !fmd->distortion_cached) {
			BKE_fracture_automerge_refresh(fmd, me_assembled);
		}

#if 0
		DEG_id_tag_update(DEG_get_original_id(&scene->id), DEG_TAG_COPY_ON_WRITE);
		DEG_id_tag_update(DEG_get_original_id(&ob->id), DEG_TAG_COPY_ON_WRITE);
		DEG_id_tag_update(DEG_get_original_id(&((Mesh*)(ob->data))->id), DEG_TAG_COPY_ON_WRITE);
#endif
	}

	/* if autohide / automerge etc perform postprocess */
	if (fmd->autohide_dist > 0 || fmd->automerge_dist > 0 || fmd->use_centroids || fmd->use_vertices)
	{
		//printf("Autohide \n");
		me_final = BKE_fracture_autohide_do(fmd, me_assembled, ob, scene);
		BKE_fracture_mesh_free(me_assembled);
		me_assembled = NULL;
	}
	else {
		me_final = me_assembled;
		if (!fmd->fix_normals) {
			BKE_mesh_calc_normals(me_final);
		}
	}

	/* return output mesh */
	return me_final;
}
