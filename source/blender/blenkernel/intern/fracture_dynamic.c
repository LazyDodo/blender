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
	int frame = (int)BKE_scene_frame_get(scene);

	{
		int count = 0;
		int totpoint = 0;
		FractureID *fid = fmd->shared->fracture_ids.first;

		while(fid) {
			FracPointCloud points;
			if (!fid->fractured) {
				points = BKE_fracture_points_get(depsgraph, fmd, ob, fid->mi);
				totpoint += points.totpoints;
				MEM_freeN(points.points);
				points.totpoints = 0;
			}
			fid = fid->next;
		}

		if (totpoint == 0)
		{
			fmd->update_dynamic = false;
		}
		else
		{
			if (fmd->update_dynamic)
			{
				if (!fmd->is_dynamic_external)
				{
					BKE_fracture_constraints_free(fmd, scene);
				}

				printf("ADD NEW 2: %s \n", ob->id.name);
				fmd->update_dynamic = false;
			}

			fid = (FractureID*)fmd->shared->fracture_ids.first;
			while(fid){
				if (!fid->fractured) {
					BKE_fracture_do(fmd, fid->mi, ob, depsgraph, bmain, scene);
				}
				fid->fractured = true;
				count++;

				fid = fid->next;
			}

#if 0
			if (count > 0)
			{
				printf("REFRESH: %s \n", ob->id.name);
				scene->rigidbody_world->flag |= RBW_FLAG_OBJECT_CHANGED;
				scene->rigidbody_world->flag |= RBW_FLAG_REFRESH_MODIFIERS;
			}
#endif
		}
	}
}

#if 0
static void fracture_dynamic_initialize(FractureModifierData *fmd, Object *ob, Mesh *dm, char (**names)[66], Scene* scene)
{
	/*TODO_1 refresh, move to BKE and just call from operator for prefractured case*/

	printf("ADD NEW 1: %s \n", ob->id.name);
	if ((fmd->last_frame == INT_MAX)) //TODO event if we restart at cache startframe
	{
		if (fmd->shared->reset_shards)
		{
			BKE_fracture_simulation_free(fmd, true, true, scene);
			BKE_fracture_free(fmd, true, true, scene);
			fmd->last_frame = 1;
		}
		else
		{
			BKE_fracture_simulation_free(fmd, false, true, scene);
			fmd->last_frame = 1;
		}

		//try to exec handlers only ONCE
		if (fmd->shared->frac_mesh == NULL) {
			// fire a callback which can then load external data right NOW
			//BLI_callback_exec(G.main, &ob->id, BLI_CB_EVT_FRACTURE_REFRESH);
			if (fmd->shared->frac_mesh != NULL) {
				bool tmp = fmd->shards_to_islands;

				fmd->shards_to_islands = false;
				BKE_fracture_assemble_mesh_from_islands(fmd, true, false);
				fmd->shards_to_islands = tmp;

				BKE_fracture_update_visual_mesh(fmd, ob, true);

				//store names here... gahhhh that is so clumsy
				if (names) {
					MeshIsland *mi;
					int i = 0, count = 0;
					count = BLI_listbase_count(&fmd->shared->meshIslands);

					(*names) = MEM_callocN(sizeof(char*) * 66 * count, "names");
					for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
						//names[i] = MEM_callocN(sizeof(char*) * 66, "names");
						BLI_snprintf((*names)[i], sizeof(mi->name), "%s", mi->name);
						i++;
					}
					//mi_count = i;
				}

				fracture_dynamic_new_entries_add(fmd, dm, ob, scene);
			}
		}
	}
	/* here we just create the fracmesh, in dynamic case we add the first sequence entry as well */
	if (fmd->shared->frac_mesh == NULL) {
		 fracture_dynamic_new_entries_add(fmd, dm, ob, scene);
	}
}
#endif

#if 0
bool BKE_fracture_dynamic_lookup_mesh_state(FractureModifierData *fmd, int frame)
{
	bool forward = false;
	bool backward = false;

	backward = ((fmd->last_frame > frame) && fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->prev);
	forward = ((fmd->last_frame < frame) && (fmd->shared->current_mi_entry) && (fmd->shared->current_mi_entry->next != NULL) &&
			   (fmd->shared->current_mi_entry->is_new == false));

	if (backward)
	{
		while (fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->prev &&
			   frame <= fmd->shared->current_mi_entry->prev->frame)
		{
			printf("Jumping backward because %d is smaller than %d\n", frame, fmd->shared->current_mi_entry->prev->frame);
			get_prev_entries(fmd);
		}
	}
	else if (forward)
	{
		while ((fmd->shared->current_mi_entry) && (fmd->shared->current_mi_entry->next != NULL) &&
			   (fmd->shared->current_mi_entry->is_new == false) &&
			   frame > fmd->shared->current_mi_entry->frame)
		{
			printf("Jumping forward because %d is greater than %d\n", frame, fmd->shared->current_mi_entry->frame);
			get_next_entries(fmd);
		}
	}

	return forward || backward;
}

static MeshIslandSequence* fracture_dynamic_meshisland_sequence_add(FractureModifierData* fmd,
                                                                    float frame, bool is_new)
{
	MeshIslandSequence *msq = MEM_mallocN(sizeof(MeshIslandSequence), "meshisland_sequence_add");
	msq->frame = (int)frame;

	if (!is_new) {
		msq->meshIslands.first = NULL;
		msq->meshIslands.last = NULL;
		msq->meshIslands = fmd->shared->mesh_islands;
		fmd->shared->refresh = false;
		msq->is_new = false;
	}
	else {
		msq->meshIslands.first = NULL;
		msq->meshIslands.last = NULL;
		msq->is_new = true;
	}

	BLI_addtail(&fmd->shared->meshIsland_sequence, msq);

	return msq;
}

void BKE_fracture_dynamic_new_entries_add(FractureModifierData* fmd, Scene *scene, bool is_new)
{
	int frame = (int)BKE_scene_frame_get(scene);
	int end = 250; //TODO might be problematic with longer sequences, take proper end value ?

	if (scene->rigidbody_world)
	{
		end = scene->rigidbody_world->shared->pointcache->endframe;
	}

	if (fmd->shared->current_mi_entry) {
		fmd->shared->current_mi_entry->frame = frame;
	}

	fmd->shared->current_mi_entry = fracture_dynamic_meshisland_sequence_add(fmd, end, is_new);
	fmd->shared->mesh_islands = fmd->shared->current_mi_entry->meshIslands;
}

static void get_next_entries(FractureModifierData *fmd)
{
	if (fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->next)
	{ // && fmd->current_mi_entry->next->is_new == false) {
		fmd->shared->current_mi_entry = fmd->shared->current_mi_entry->next;
		fmd->shared->mesh_islands = fmd->shared->current_mi_entry->meshIslands;
	}
}

static void get_prev_entries(FractureModifierData *fmd)
{
	if (fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->prev)
	{
		fmd->shared->current_mi_entry = fmd->shared->current_mi_entry->prev;
		fmd->shared->mesh_islands = fmd->shared->current_mi_entry->meshIslands;
	}
}
#endif
