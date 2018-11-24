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

#include <stdio.h>

#include "MEM_guardedalloc.h"

#include "BKE_fracture.h"
#include "BKE_mesh.h"
#include "BKE_deform.h"
#include "BKE_material.h"
#include "BKE_main.h"
#include "BKE_modifier.h"
#include "BKE_rigidbody.h"
#include "BKE_pointcache.h"
#include "BKE_object.h"

#include "BLI_listbase.h"
#include "BLI_ghash.h"
#include "BLI_math.h"
#include "BLI_math_vector.h"
#include "BLI_string.h"

#include "DNA_object_types.h"
#include "DNA_modifier_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_fracture_types.h"
#include "DNA_material_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

static Shard* fracture_object_to_island(FractureModifierData *fmd, Object *own, Object* target, Scene* scene);

static int fracture_collect_defgrp(Object* o, Object* ob, int defstart, GHash** def_index_map)
{
	bDeformGroup *vgroup, *ngroup;
	int k = 0;

	/* create vertexgroups on new object, if they dont exist already there*/
	for (vgroup = o->defbase.first; vgroup; vgroup = vgroup->next) {
		int index = defgroup_name_index(ob, vgroup->name);
		int key = defstart + k;

		if (index == -1) {
			// old group index + defstart to make it somehow linearized
			ngroup = MEM_callocN(sizeof(bDeformGroup), "collect deformGroup");
			memcpy(ngroup, vgroup, sizeof(bDeformGroup));
			BLI_addtail(&ob->defbase, ngroup);
			index = BLI_listbase_count(&ob->defbase)-1;
		}

		if (!BLI_ghash_haskey(*def_index_map, POINTER_FROM_INT(key)))
			BLI_ghash_insert(*def_index_map, POINTER_FROM_INT(key), POINTER_FROM_INT(index));

		k++;
	}

	if (ob->defbase.first && ob->actdef == 0)
		ob->actdef = 1;

	return k;
}

short BKE_fracture_collect_materials(Main* bmain, Object* o, Object* ob, int matstart, GHash** mat_index_map)
{
	short *totcolp = NULL;
	Material ***matarar = NULL;
	int j;

	/* append materials to target object, if not existing yet */
	totcolp = give_totcolp(o);
	matarar = give_matarar(o);

	for (j = 0; j < (*totcolp); j++)
	{
		void *key;
		int index = BKE_object_material_slot_find_index(ob, (*matarar)[j]);
		if (index == 0)
		{
			index = ob->totcol+1;
			assign_material(bmain, ob, (*matarar)[j], index, BKE_MAT_ASSIGN_USERPREF);
		}

		key = POINTER_FROM_INT(matstart+j);
		if (!BLI_ghash_haskey(*mat_index_map, key))
			BLI_ghash_insert(*mat_index_map, key, POINTER_FROM_INT(index));
	}

	return (*totcolp);
}

Shard* BKE_fracture_mesh_island_add(Main* bmain, FractureModifierData *fmd, Object* own, Object *target, Scene *scene)
{
	Shard *mi;
	short totcol = 0, totdef = 0;
	float loc[3], quat[4];
	int endframe = scene->rigidbody_world->shared->pointcache->endframe;
	//better: leave open... or update when changing endframe / startframe TODO FIX

	if (own->type != OB_MESH || !own->data)
		return NULL;

	if (target->type != OB_MESH || !target->data)
		return NULL;

	//lets see whether we need to add loc here too XXX TODO
	mat4_to_loc_quat(loc, quat, target->obmat);

	mi = fracture_object_to_island(fmd, own, target, scene);
	mi->endframe = endframe;

	copy_v3_v3(mi->loc, loc);
	copy_qt_qt(mi->rot, quat);
	copy_v3_v3(mi->loc, loc);

	//mi->rigidbody = BKE_rigidbody_create_shard(bmain, scene, own, target, mi);
	BLI_strncpy(mi->name, target->id.name + 2, MAX_ID_NAME - 2);

	//handle materials
	if (!fmd->shared->material_index_map)
	{
		fmd->shared->material_index_map = BLI_ghash_int_new("mat_index_map");
		fmd->shared->matstart = 1;
	}

	totcol = BKE_fracture_collect_materials(bmain, target, own, fmd->shared->matstart, &fmd->shared->material_index_map);
	if (totcol < 0)
		totcol = 0;
	fmd->shared->matstart += totcol;
	mi->totcol = totcol;

	/*XXXXX TODO deal with material deletion, and reorder (in material code) */

	//handle vertexgroups, too
	if (!fmd->shared->defgrp_index_map)
	{
		fmd->shared->defgrp_index_map = BLI_ghash_int_new("defgrp_index_map");
		fmd->shared->defstart = 0;
	}

	totdef = fracture_collect_defgrp(target, own, fmd->shared->defstart, &fmd->shared->defgrp_index_map);
	if (totdef < 0)
		totdef = 0;
	fmd->shared->defstart += totdef;
	mi->totdef = totdef;

	//XXX TODO handle UVs, shapekeys and more ?
//	fracture_collect_uv_tex(target, own);

	//add shard to pack storage
//	pack_storage_add(fmd, mi); //TODO FIX LOOKUP

	return mi;
}

void BKE_fracture_mesh_island_free(Shard *mi, Scene* scene)
{
	if (mi->mesh) {
		BKE_fracture_mesh_free(mi->mesh);
		mi->mesh = NULL;
	}

	if (mi->rigidbody) {
		if (scene) {
			BKE_rigidbody_remove_shard(scene, mi);
		}
		MEM_freeN(mi->rigidbody->shared);
		MEM_freeN(mi->rigidbody);
		mi->rigidbody = NULL;
	}

	if (mi->participating_constraints != NULL) {
		int i = 0;
		for (i = 0; i < mi->participating_constraint_count; i++)
		{
			RigidBodyShardCon *con = mi->participating_constraints[i];
			if (con) {
				con->mi1 = NULL;
				con->mi2 = NULL;
			}
		}

		MEM_freeN(mi->participating_constraints);
		mi->participating_constraints = NULL;
		mi->participating_constraint_count = 0;
	}

	if (mi->rots) {
		MEM_freeN(mi->rots);
		mi->rots = NULL;
	}

	if (mi->locs) {
		MEM_freeN(mi->locs);
		mi->locs = NULL;
	}

	if (mi->vels) {
		MEM_freeN(mi->vels);
		mi->vels = NULL;
	}

	if (mi->aves) {
		MEM_freeN(mi->aves);
		mi->aves = NULL;
	}

	if (mi->neighbors) {
		MEM_freeN(mi->neighbors);
		mi->neighbors = NULL;
		mi->neighbor_count = 0;
	}

	MEM_freeN(mi);
	mi = NULL;
}

void BKE_fracture_mesh_island_remove(FractureModifierData *fmd, Shard *mi, Scene* scene)
{
	if (BLI_listbase_is_single(&fmd->shared->shards))
	{
		BKE_fracture_mesh_island_remove_all(fmd, scene);
		return;
	}

	if (mi)
	{
		int i;

		BLI_remlink(&fmd->shared->shards, mi);
//		pack_storage_remove(fmd, mi, scene);

		for (i = 0; i < mi->participating_constraint_count; i++)
		{
			RigidBodyShardCon *con = mi->participating_constraints[i];
			BLI_remlink(&fmd->shared->constraints, con);
			BKE_rigidbody_remove_shard_con(scene->rigidbody_world, con);
		}

		BKE_fracture_mesh_island_free(mi, scene);
	}
}

void BKE_fracture_mesh_island_remove_all(FractureModifierData *fmd, Scene* scene)
{
	Shard *mi;

	//free pack storage
//	pack_storage_remove_all(fmd, scene);

	//free all constraints first
	BKE_fracture_constraints_free(fmd, scene);

	//free all meshislands
	while (fmd->shared->shards.first) {
		mi = fmd->shared->shards.first;
		BLI_remlink(&fmd->shared->shards, mi);
		BKE_fracture_mesh_island_free(mi, scene);
	}

	fmd->shared->shards.first = NULL;
	fmd->shared->shards.last = NULL;

	if (fmd->shared->material_index_map)
	{
		BLI_ghash_free(fmd->shared->material_index_map, NULL, NULL);
		fmd->shared->material_index_map = NULL;
		fmd->shared->matstart = 1;
	}

	if (fmd->shared->defgrp_index_map)
	{
		BLI_ghash_free(fmd->shared->defgrp_index_map, NULL, NULL);
		fmd->shared->defgrp_index_map = NULL;
		fmd->shared->defstart = 0;
	}
}

static Shard* fracture_object_to_island(FractureModifierData* fmd, Object *own, Object* target, Scene* scene)
{
	Shard *mi;
	Mesh *me = target->runtime.mesh_eval;
	MVert *mv;
	int v;
	SpaceTransform trans;
	float mat[4][4];
	int frame = (int)BKE_scene_frame_get(scene);

	if (me)
	{	//fallback if no derivedFinal available
		me = BKE_fracture_mesh_copy(me, own);
	}
	else {
		me = BKE_fracture_mesh_copy(target->data, own);
	}

	unit_m4(mat);
	BLI_space_transform_from_matrices(&trans, target->obmat, mat);

	// create temp shard -> that necessary at all ?
	mi = BKE_fracture_mesh_island_create(me, scene, own, frame);
	BLI_addtail(&fmd->shared->shards, mi);

	for (v = 0, mv = mi->mesh->mvert; v < mi->mesh->totvert; v++, mv++)
	{
		//shrink the shard ? (and take centroid diff into account here, too)
		BLI_space_transform_apply(&trans, mv->co);
	}

	BLI_space_transform_apply(&trans, mi->loc);

	return mi;
}
