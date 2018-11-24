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

MeshIsland *BKE_fracture_mesh_island_create(Mesh* me, Main* bmain, Scene *scene, Object *ob, int frame)
{
	int endframe = scene->rigidbody_world->shared->pointcache->endframe;
	MeshIsland *mi = MEM_callocN(sizeof(MeshIsland), "mesh_island");
	mi->mesh = me;

	INIT_MINMAX(mi->min, mi->max);
	BKE_mesh_minmax(mi->mesh, mi->min, mi->max);
	BKE_fracture_mesh_center_centroid_area(mi->mesh, mi->centroid);
	unit_qt(mi->rot);

	mi->startframe = frame;
	mi->endframe = endframe;

	if (endframe >= frame) {
		frame = endframe - frame + 1;
		mi->locs = MEM_callocN(sizeof (float) * 3 *frame, "mi->locs");
		mi->rots = MEM_callocN(sizeof (float) * 4 *frame, "mi->rots");
		mi->vels = MEM_callocN(sizeof (float) * 3 *frame, "mi->vels");
		mi->aves = MEM_callocN(sizeof (float) * 3 *frame, "mi->aves");
	}

	mi->rigidbody = BKE_rigidbody_create_shard(bmain, scene, ob, NULL, mi);
	mi->rigidbody->type = RBO_TYPE_ACTIVE;
	mi->rigidbody->flag |= (RBO_FLAG_NEEDS_VALIDATE | RBO_FLAG_NEEDS_RESHAPE);
	BKE_rigidbody_calc_shard_mass(ob, mi);

	return mi;
}

static bool handle_initial_shards(FractureModifierData* fmd, Object* ob, Depsgraph *depsgraph, Main* bmain, Scene* scene, int frame)
{
	MeshIsland *mi = NULL;
	MeshIsland** mi_tmp = NULL;
	int i = 0;

	int count = 0;
	for (mi = fmd->shared->mesh_islands.first; mi; mi = mi->next)
	{
		if (!BKE_fracture_meshisland_check_frame(fmd, mi, frame) && mi->id > 0)
		{
			count++;
		}
	}

	if (count == 0) {
		return false;
	}

	mi_tmp = MEM_callocN(sizeof(MeshIsland*) * count, "mi_tmp");
	for (mi = fmd->shared->mesh_islands.first; mi; mi = mi->next)
	{
		if (!BKE_fracture_meshisland_check_frame(fmd, mi, frame) && mi->id > 0)
		{
			mi_tmp[i] = mi;
			i++;
		}
	}

	/*decouple from listbase because it will continue growing ... */
	for (mi = mi_tmp[0], i = 0; i < count; i++)
	{
		BKE_fracture_do(fmd, mi_tmp[i], ob, depsgraph, bmain, scene, true);
		mi_tmp[i]->endframe = frame;
	}

	MEM_freeN(mi_tmp);

	return true;
}

static void do_initial_prefracture(FractureModifierData* fmd, Object* ob, Depsgraph *depsgraph, Main* bmain,
                            Scene* scene, int frame, Mesh* me)
{
	MeshIsland *mi = NULL;
	Mesh *me_tmp = NULL;

	BKE_fracture_meshislands_free(fmd, scene);
	me_tmp = BKE_fracture_mesh_copy(me, ob);

	mi = BKE_fracture_mesh_island_create(me_tmp, bmain, scene, ob, frame);
	mi->id = 0;
	BLI_addtail(&fmd->shared->mesh_islands, mi);

	if (fmd->shared->last_islands) {
		MEM_freeN(fmd->shared->last_islands);
		fmd->shared->last_islands = NULL;
		fmd->shared->last_expected_islands = 0;
	}

	fmd->shared->last_islands = MEM_callocN(sizeof(MeshIsland*), "island initial");
	fmd->shared->last_islands[0] = mi;
	fmd->shared->last_expected_islands = 1;

	BKE_fracture_do(fmd, mi, ob, depsgraph, bmain, scene, true);

	//if ((fmd->point_source & MOD_FRACTURE_CUSTOM) == 0)
	//	mi->endframe = frame;
}


Mesh* BKE_fracture_apply(FractureModifierData *fmd, Object *ob, Mesh *me_orig, Depsgraph* depsgraph)
{
	Scene *scene = DEG_get_input_scene(depsgraph);
	float ctime = BKE_scene_frame_get(scene);
	int frame = (int)ctime;
	RigidBodyWorld* rbw = scene->rigidbody_world;

	Main* bmain = G.main;
	Mesh* me_assembled = NULL;
	Mesh *me_final = NULL;
	Mesh *me = me_orig;

	if (fmd->auto_execute && !BKE_rigidbody_check_sim_running(rbw, ctime)) {
		if (ob->rigidbody_object)
			fmd->shared->refresh = true;
	}

	if (fmd->shared->refresh /*|| fmd->shared->reset_shards*/)
	{
		/*reset_shards called from readfile.c; refresh from operator */

		/*free old stuff here */
		BKE_fracture_constraints_free(fmd, scene);

		/*keep shards at packing and at dynamic refresh */
		if (fmd->dm_group)
		{
			if (!handle_initial_shards(fmd, ob, depsgraph, bmain, scene, frame))
			{
				do_initial_prefracture(fmd, ob, depsgraph, bmain, scene, frame, me);
			}
		}
		else {
			if (BLI_listbase_is_empty(&fmd->shared->mesh_islands) || !fmd->use_dynamic) {
				/*rebuild shards after loading and prefracture refresh*/
				do_initial_prefracture(fmd, ob, depsgraph, bmain, scene, frame, me);
			}
		}

		fmd->shared->refresh_constraints = true;
		fmd->shared->refresh_autohide = true;
		//fmd->shared->reset_shards = false;
	}
	else if (fmd->shared->refresh_dynamic) {
		//handle_initial_shards(fmd, ob, depsgraph, bmain, scene, frame);

		/*if dynamic event, push state to fracture sequence*/
		BKE_fracture_dynamic_do(fmd, ob, scene, depsgraph, bmain);
		fmd->shared->refresh_dynamic = false;
		fmd->shared->refresh_constraints = true;
		fmd->shared->refresh_autohide = true;
	}

	//if (scene->rigidbody_world)
	//	BKE_rigidbody_modifier_update(scene, ob, scene->rigidbody_world, false, depsgraph);

	if (fmd->use_animated_mesh && fmd->anim_mesh_ob)
	{
		/*update bound island positions to follow bind object, after physics ran */
		BKE_fracture_animated_loc_rot(fmd, ob, false, depsgraph);
	}

	//BKE_rigidbody_modifier_sync((ModifierData*)fmd, ob, scene, ctime);

	/* assemble mesh from transformed meshislands */
	if (fmd->shared->mesh_islands.first) {
		me_assembled = BKE_fracture_assemble_mesh_from_islands(fmd, ob, ctime);
	}
	else {
		me_assembled = BKE_fracture_mesh_copy(me, ob);
	}

	fmd->shared->mesh_cached = me_assembled;

	/*if refresh constraints, build constraints */
	if (fmd->shared->refresh_constraints) {
		if (!fmd->shared->refresh)
			// do not free twice
			BKE_fracture_constraints_free(fmd, scene);

		BKE_fracture_external_constraints_setup(fmd, scene, ob);
	}

	fmd->shared->refresh = false;

	/*if refresh autohide, initialize its structures */
	if (fmd->shared->refresh_autohide) {
		fmd->distortion_cached = false;
		BKE_fracture_autohide_refresh(fmd, ob, me_assembled);

		if (fmd->autohide_dist > 0 && !fmd->distortion_cached) {
			BKE_fracture_automerge_refresh(fmd, me_assembled);
		}
	}

	if (me_assembled->totvert == 0) {
		/* just return the original mesh in case our mesh is empty */
		BKE_fracture_mesh_free(me_assembled);
		me_final = me_orig;
	}
	else
	{
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
	}

	fmd->last_frame = frame;

	/* return output mesh */
	return me_final;
}
