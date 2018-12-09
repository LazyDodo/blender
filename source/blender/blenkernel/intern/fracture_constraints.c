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

#include "BKE_collection.h"
#include "BKE_fracture.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_pointcache.h"
#include "BKE_rigidbody.h"

#include "BLI_kdtree.h"
#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_math_matrix.h"
#include "BLI_ghash.h"
#include "BLI_string.h"

#include "DNA_fracture_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_collection_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

#include "DEG_depsgraph_query.h"

#include "PIL_time.h"

static void search_tree_based(FractureModifierData *rmd, Shard *mi, Shard **meshIslands, KDTree **combined_tree,
							  float co[3], Object *ob, Scene *scene);

static void remove_participants(RigidBodyShardCon* con, Shard *mi);


static int get_object_index(Scene *scene, Object *ob) {
	RigidBodyWorld *rbw = scene->rigidbody_world;
	Object *obb, *obbj;
	int i = 0;
	if (rbw && rbw->group) {

		obb = DEG_get_original_object(ob);
		FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rbw->group, obj)
		{
			obbj = DEG_get_original_object(obj);
			if (obb == obbj) {
				return i;
			}
			i++;
		}
		FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
	}

	return -1;
}

bool BKE_fracture_check_valid_shard(FractureModifierData *fmd, Shard *mi, Scene *scene)
{
	//why for gods sake is the correct frame only in the original scene ?
	Scene* sc = (Scene*)DEG_get_original_id(&scene->id);
	int frame = BKE_scene_frame_get(sc);

	return !BKE_fracture_meshisland_check_frame(fmd, mi, frame);
}

static int count_valid_shards(FractureModifierData *fmd, Scene *scene) {
	int count = 0;
	Shard *mi;

	for (mi = fmd->shared->shards.first; mi; mi = mi->next) {
		if (BKE_fracture_check_valid_shard(fmd, mi, scene)) {
			count++;
		}
	}

	return count;
}

static int prepareConstraintSearch(FractureModifierData *rmd, Shard ***mesh_islands, KDTree **combined_tree, Object *obj,
								   MVert** mverts, Scene *scene)
{
	Shard *mi;
	int i = 0, ret = 0;
	int islands = 0;

	if (rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY))
	{
		FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rmd->pack_group, ob)
		{
			if (obj != ob)
			{
				FractureModifierData *fmdi = (FractureModifierData *)modifiers_findByType(ob, eModifierType_Fracture);
				if (fmdi) {
					islands += count_valid_shards(fmdi, scene);
				}
			}
		}
		FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
	}
	else {

		islands = count_valid_shards(rmd, scene);
	}

	*mesh_islands = MEM_reallocN(*mesh_islands, islands * sizeof(Shard *));

	if (rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY))
	{

		FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rmd->pack_group, ob)
		{
			if (obj != ob)
			{
				FractureModifierData *fmdi = (FractureModifierData *)modifiers_findByType(ob, eModifierType_Fracture);
				if (fmdi) {
					for (mi = fmdi->shared->shards.first; mi; mi = mi->next) {
						mi->object_index = get_object_index(scene, ob);
						if (BKE_fracture_check_valid_shard(fmdi, mi, scene))
						{
							(*mesh_islands)[i] = mi;
							i++;
						}
					}
				}

				//j++;
			}
		}
		FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
	}
	else {

		for (mi = rmd->shared->shards.first; mi; mi = mi->next) {
			mi->object_index = get_object_index(scene, obj);
			if (BKE_fracture_check_valid_shard(rmd, mi, scene))
			{
				(*mesh_islands)[i] = mi;
				i++;
			}
		}
	}

	if (rmd->constraint_target == MOD_FRACTURE_CENTROID)
	{
		*combined_tree = BLI_kdtree_new(islands);
		for (i = 0; i < islands; i++) {
			float obj_centr[3];
			if (!(rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY)))
			{
				mul_v3_m4v3(obj_centr, obj->obmat, (*mesh_islands)[i]->loc);
			}
			else {
				//TODO, take care of recursive collections here !!!!
				CollectionObject *go = BLI_findlink(&scene->rigidbody_world->group->gobject,
				                                    (*mesh_islands)[i]->object_index);

				mul_v3_m4v3(obj_centr, go->ob->obmat, (*mesh_islands)[i]->loc);
				//copy_v3_v3(obj_centr, (*mesh_islands)[i]->loc);
			}

			BLI_kdtree_insert(*combined_tree, i, obj_centr);
		}

		BLI_kdtree_balance(*combined_tree);
		ret = islands;
	}
	else if (rmd->constraint_target == MOD_FRACTURE_VERTEX)
	{
		int totvert = 0;
		MVert *mvert = NULL;
		MVert *mv;

		if (rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY))
		{
			mvert = MEM_mallocN(sizeof(MVert), "mvert");

			FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rmd->pack_group, ob)
			{
				if (obj != ob)
				{
					float imat[4][4];
					FractureModifierData *fmdi = (FractureModifierData *)modifiers_findByType(ob, eModifierType_Fracture);
					if (fmdi && fmdi->shared->mesh_cached) {
						int v = fmdi->shared->mesh_cached->totvert;
						int x = 0;

						//invert_m4_m4(imat, go->ob->obmat);
						copy_m4_m4(imat, ob->obmat);

						mvert = MEM_reallocN(mvert, sizeof(MVert) * (totvert + v));
						memcpy(mvert + totvert, fmdi->shared->mesh_cached->mvert, v * sizeof(MVert));

						for (x = totvert; x < totvert + v; x++)
						{
							mul_v3_m4v3(mvert[x].co, imat, mvert[x].co);
						}

						totvert += v;
					}
				}
			}
			FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
		}
		else if (rmd && rmd->shared->mesh_cached) {
			totvert = rmd->shared->mesh_cached->totvert;
			mvert = rmd->shared->mesh_cached->mvert;
		}

		if (totvert > 0)
		{
			*combined_tree = BLI_kdtree_new(totvert);
			for (i = 0, mv = mvert; i < totvert; i++, mv++) {
				float co[3];
				if (rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY))
				{
					copy_v3_v3(co, mv->co);
				}
				else {
					mul_v3_m4v3(co, obj->obmat, mv->co);
				}

				BLI_kdtree_insert(*combined_tree, i, co);
			}

			BLI_kdtree_balance(*combined_tree);
			ret = totvert;
			*mverts = mvert;
		}
	}

	return ret;
}

static void create_constraints(FractureModifierData *rmd, Object *ob, Scene *scene)
{
	KDTree *coord_tree = NULL;
	Shard **mesh_islands = MEM_mallocN(sizeof(Shard *), "mesh_islands");
	int count, i = 0;
	Shard *mi;
	MVert *mvert = NULL;

	float max_mass = 0.0f;
	for (mi = rmd->shared->shards.first; mi; mi = mi->next)
	{
		if (mi->rigidbody->mass > max_mass)
			max_mass = mi->rigidbody->mass;
	}

	count = prepareConstraintSearch(rmd, &mesh_islands, &coord_tree, ob, &mvert, scene);

	for (i = 0; i < count; i++) {
		if (rmd->constraint_target == MOD_FRACTURE_CENTROID) {
			float co[3];
			if (rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY))
			{
				CollectionObject *go = BLI_findlink(&scene->rigidbody_world->group->gobject, mesh_islands[i]->object_index);
				mul_v3_m4v3(co, go->ob->obmat, mesh_islands[i]->loc);
			}
			else {
				copy_v3_v3(co, mesh_islands[i]->loc);
			}

			search_tree_based(rmd, mesh_islands[i], mesh_islands, &coord_tree, co, ob, scene);
		}
		else if (rmd->constraint_target == MOD_FRACTURE_VERTEX) {
			Shard *mii = NULL;
			mii = BLI_ghash_lookup(rmd->shared->vertex_island_map, POINTER_FROM_INT(i));
			search_tree_based(rmd, mii, mesh_islands, &coord_tree, mvert[i].co, ob, scene);
		}
	}

	if (coord_tree != NULL) {
		BLI_kdtree_free(coord_tree);
		coord_tree = NULL;
	}

	MEM_freeN(mesh_islands);

	if (rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY) && mvert)
	{	//was copied from modifiers... so remove now
		MEM_freeN(mvert);
	}
}

static void search_tree_based(FractureModifierData *rmd, Shard *mi, Shard **meshIslands,
							  KDTree **combined_tree, float co[3], Object *ob, Scene *scene)
{
	int r = 0, limit = 0, i = 0;
	KDTreeNearest *n3 = NULL;
	float dist, obj_centr[3];

	limit = rmd->constraint_limit;
	dist = rmd->contact_dist;

	if (!(rmd->pack_group && (rmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY)))
	{
		mul_v3_m4v3(obj_centr, ob->obmat, co);
	}
	else {
		copy_v3_v3(obj_centr, co);
	}

	r = BLI_kdtree_range_search(*combined_tree, obj_centr, &n3, dist);

	/* use centroid dist based approach here, together with limit */
	for (i = 0; i < r; i++) {
		Shard *mi2 = NULL;

		if (rmd->constraint_target == MOD_FRACTURE_CENTROID) {
			mi2 = meshIslands[(n3 + i)->index];
		}
		else if(rmd->constraint_target == MOD_FRACTURE_VERTEX) {
			int index = (n3 + i)->index;
			mi2 = BLI_ghash_lookup(rmd->shared->vertex_island_map, POINTER_FROM_INT(index));
		}
		if ((mi != mi2) && (mi2 != NULL)) {
			float thresh = rmd->breaking_threshold;
			int con_type = rmd->constraint_type;
			int startframe = 1;

			if ((i >= limit) && (limit > 0)) {
				break;
			}

			if (scene->rigidbody_world) {
				startframe = scene->rigidbody_world->shared->pointcache->startframe;
			}

			if (rmd->dynamic_new_constraints == MOD_FRACTURE_MIXED_DYNAMIC_CONSTRAINTS) {
				//only build between old and new (old = mi startframe == pointcache startframe)
				if (mi->startframe > startframe && mi2->startframe > startframe) {
					continue;
				}
			}

			if (rmd->dynamic_new_constraints == MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS)
			{	//only build between old
				if (mi->startframe > startframe || mi2->startframe > startframe) {
					continue;
				}
			}

			if (rmd->contact_size > 0.0f) {
				//limit to shards above the size, atleast one shard must be bigger
				if (BKE_rigidbody_check_island_size(rmd, mi, rmd->contact_size) &&
				    BKE_rigidbody_check_island_size(rmd, mi2, rmd->contact_size))
				{
					continue;
				}
			}

			/* does some more checks */
			BKE_fracture_meshislands_connect(scene, rmd, mi, mi2, con_type, thresh);
		}
	}

	if (n3 != NULL) {
		MEM_freeN(n3);
		n3 = NULL;
	}
}

void BKE_fracture_meshislands_connect(Scene* sc, FractureModifierData *fmd, Shard *mi1, Shard *mi2,
											 short con_type, float thresh)
{
	int con_found = false;
	RigidBodyShardCon *con;
	bool ok = mi1 && mi1->rigidbody;
	bool con_group = fmd->flag & MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY;
	ok = ok && mi2 && mi2->rigidbody;
	ok = ok && (fmd->flag & MOD_FRACTURE_USE_CONSTRAINTS);
	ok = ok && ((!(fmd->pack_group && con_group) && (mi1->object_index == mi2->object_index))||
		 (fmd->pack_group && con_group && (mi1->object_index != mi2->object_index)));

	if (ok) {
		/* search local constraint list instead of global one !!! saves lots of time */
		int i;
		for (i = 0; i < mi1->participating_constraint_count; i++) {
			con = mi1->participating_constraints[i];
			if (con && ((con->mi1 == mi2) || (con->mi2 == mi2))) {
				con_found = true;
				break;
			}
		}

		if (!con_found) {
			for (i = 0; i < mi2->participating_constraint_count; i++) {
				con = mi2->participating_constraints[i];
				if (con && ((con->mi1 == mi1) || (con->mi2 == mi1))) {
					con_found = true;
					break;
				}
			}
		}
	}

	if (!con_found && ok) {
		BKE_fracture_constraint_create(sc, fmd, mi1, mi2, con_type, thresh);
	}
}

void BKE_fracture_mesh_constraint_remove(FractureModifierData *fmd, RigidBodyShardCon* con, Scene *scene)
{
	remove_participants(con, con->mi1);
	remove_participants(con, con->mi2);

	BLI_remlink(&fmd->shared->constraints, con);
	BKE_rigidbody_remove_shard_con(scene->rigidbody_world, con);
	MEM_freeN(con);
}

void BKE_fracture_mesh_constraint_remove_all(FractureModifierData *fmd, Scene *scene)
{
	BKE_fracture_constraints_free(fmd, scene->rigidbody_world);
}

static void remove_participants(RigidBodyShardCon* con, Shard *mi)
{
	RigidBodyShardCon **cons;
	/* Probably wrong, would need to shrink array size... listbase would have been better here */
	/* info not necessary so omit */
	int count = 0;

	if (!mi)
		return;

	count = mi->participating_constraint_count;

	if (count > 1)
	{
		int i, j = 0;
		mi->participating_constraint_count--;
		cons = MEM_callocN(sizeof(RigidBodyShardCon *) * (count-1) , "temp cons");
		for (i = 0; i < mi->participating_constraint_count; i++)
		{
			if (mi->participating_constraints[i] != con)
			{
				cons[j] = con;
				j++;
			}
		}

		MEM_freeN(mi->participating_constraints);
		mi->participating_constraints = cons;
	}
}

void BKE_fracture_constraint_create(Scene* scene, FractureModifierData* fmd, Shard *mi1, Shard *mi2, short con_type, float thresh)
{
	RigidBodyShardCon *rbsc;
	rbsc = BKE_rigidbody_create_shard_constraint(scene, con_type, false); // !fmd->use_dynamic);
	rbsc->mi1 = mi1;
	rbsc->mi2 = mi2;

	BLI_snprintf(rbsc->name, 64, "%s-%s", rbsc->mi1->name, rbsc->mi2->name);

	if (thresh == 0 || !(fmd->flag & MOD_FRACTURE_USE_BREAKING)) {
		rbsc->flag &= ~RBC_FLAG_USE_BREAKING;
	}

	if (!(fmd->flag & MOD_FRACTURE_USE_CONSTRAINT_COLLISION)) {
		rbsc->flag |= RBC_FLAG_DISABLE_COLLISIONS;
	}
	else {
		rbsc->flag &= ~RBC_FLAG_DISABLE_COLLISIONS;
	}

	rbsc->flag |= RBC_FLAG_NEEDS_VALIDATE;

	if ((mi1->cluster_index != -1) && (mi2->cluster_index != -1) &&
		(mi1->cluster_index == mi2->cluster_index))
	{
		if (fmd->cluster_count > 1) {
			rbsc->breaking_threshold = fmd->cluster_breaking_threshold;
		}
		else {
			rbsc->breaking_threshold = thresh;
		}
	}
	else
	{
		if ((mi1->cluster_index != -1) && (mi2->cluster_index != -1) &&
			(mi1->cluster_index != mi2->cluster_index))
		{
			/* set a different type of constraint between clusters */
			rbsc->type = fmd->cluster_constraint_type;
		}
		rbsc->breaking_threshold = thresh;
	}

	if (fmd->thresh_defgrp_name[0]) {
		/* modify maximum threshold by minimum weight */
		rbsc->breaking_threshold = thresh * MIN2(mi1->thresh_weight, mi2->thresh_weight);
	}

	BLI_addtail(&fmd->shared->constraints, rbsc);

	if ((mi1->object_index == mi2->object_index))
	{
		/* store constraints per meshisland too, to allow breaking percentage */
		if (mi1->participating_constraints == NULL) {
			mi1->participating_constraints = MEM_callocN(sizeof(RigidBodyShardCon *), "part_constraints_mi1");
			mi1->participating_constraint_count = 0;
		}
		mi1->participating_constraints = MEM_reallocN(mi1->participating_constraints,
													  sizeof(RigidBodyShardCon *) * (mi1->participating_constraint_count + 1));
		mi1->participating_constraints[mi1->participating_constraint_count] = rbsc;
		mi1->participating_constraint_count++;

		if (mi2->participating_constraints == NULL) {
			mi2->participating_constraints = MEM_callocN(sizeof(RigidBodyShardCon *), "part_constraints_mi2");
			mi2->participating_constraint_count = 0;
		}
		mi2->participating_constraints = MEM_reallocN(mi2->participating_constraints,
													  sizeof(RigidBodyShardCon *) * (mi2->participating_constraint_count + 1));
		mi2->participating_constraints[mi2->participating_constraint_count] = rbsc;
		mi2->participating_constraint_count++;
	}
}

static void do_cluster_count(FractureModifierData *fmd, Object *obj)
{
	int k = 0;
	KDTree *tree;
	Shard *mi, **seeds;
	int seed_count, group_count = 0;
	float mat[4][4];
	CollectionObject *go = NULL;

	int mi_count;
	invert_m4_m4(mat, obj->obmat);

	/* zero clusters or one mean no clusters, all shards keep free */
	if (fmd->cluster_count < 1 && !fmd->cluster_group) {
		return;
	}

	/*initialize cluster "colors" -> membership of meshislands to clusters, initally all shards are "free" */
	for (mi = fmd->shared->shards.first; mi; mi = mi->next ) {
		mi->cluster_index = -1;
	}

	mi_count = BLI_listbase_count(&fmd->shared->shards);
	seed_count = (fmd->cluster_count > mi_count ? mi_count : fmd->cluster_count);
	//seed_count = fmd->cluster_count;

	if (fmd->cluster_group)
	{
		 group_count = BLI_listbase_count(&fmd->cluster_group->gobject);
	}

	seeds = MEM_mallocN(sizeof(Shard *) * (seed_count), "seeds");
	tree = BLI_kdtree_new(seed_count + group_count);
	//shard_tree = BLI_kdtree_new(mi_count);

	/* pick n seed locations, randomly scattered over the object */
	for (k = 0; k < seed_count; k++) {
		int which_index = k * (int)(mi_count / seed_count);
		Shard *which = (Shard *)BLI_findlink(&fmd->shared->shards, which_index);
		which->cluster_index = k;
		//print_v3("INSERT", which->loc);
		BLI_kdtree_insert(tree, k, which->loc);
		seeds[k] = which;
	}

	/*add the group here */
	if (fmd->cluster_group) {
		for (k = seed_count, go = fmd->cluster_group->gobject.first; go; k++, go = go->next)
		{
			float loc[3];

			mul_v3_m4v3(loc, mat, go->ob->loc);

			//print_v3("INSERT", loc);
			BLI_kdtree_insert(tree, k, loc);
		}
	}

	BLI_kdtree_balance(tree);


	/* assign each shard to its closest center */
	for (mi = fmd->shared->shards.first; mi; mi = mi->next ) {
		KDTreeNearest n;
		int index;
		Object *obb = NULL;
		bool in_shape = true;

		index = BLI_kdtree_find_nearest(tree, mi->loc, &n);
		if (index >= seed_count) {
			//TODO, recursive collection find by index
			CollectionObject *cob = NULL;
			cob = BLI_findlink(&fmd->cluster_group->gobject, index - seed_count);
			if (cob)
				obb = cob->ob;
		}

		/* test whether within shape if seed is a mesh object*/
		if (obb && obb->type == OB_MESH)
		{
			/* use geometry of meshes */
			MVert* mvert = NULL;
			MLoop* mloop = NULL;
			MPoly* mpoly = NULL, *mp = NULL;
			Mesh *dm = obb->runtime.mesh_eval;
			int totpoly, p;
			float mi_diff[3], normalized_mi_diff[3], len_mi;

			sub_v3_v3v3(mi_diff, mi->loc, obb->loc);
			normalize_v3_v3(normalized_mi_diff, mi_diff);
			len_mi = len_v3(mi_diff);

			if (!dm) {
				dm = obb->data;
			}

			mloop = dm->mloop;
			mvert = dm->mvert;
			mpoly = dm->mpoly;
			totpoly = dm->totpoly;

			for (p = 0, mp = mpoly; p < totpoly; p++, mp++)
			{
				//third point is always ob->loc
				float centr[3] = {0, 0, 0}, co[3], co2[3], centr_diff[3], len, dot, normalized[3];
				int i;
				in_shape = true;

				for (i = 0; i < mp->totloop; i++)
				{
					MVert *v, *v2;
					int index, index2;
					float co_ob[3], co2_ob[3], normal[3];

					index = mloop[i+mp->loopstart].v;
					if (i < mp->totloop-1)
					{
						index2 = mloop[i+1+mp->loopstart].v;
					}
					else {
						index2 = mloop[mp->loopstart].v;
					}

					v = mvert + index;
					v2 = mvert + index2;
					add_v3_v3(centr, v->co);

					mul_v3_m4v3(co, obb->obmat, v->co);
					mul_v3_m4v3(co2, obb->obmat, v2->co);

					//check triangle to objcentr... (as part of "frustum" this poly forms)
					sub_v3_v3v3(co_ob, co, obb->loc);
					sub_v3_v3v3(co2_ob, co2, obb->loc);
					cross_v3_v3v3(normal, co_ob, co2_ob);

					if (dot_v3v3(mi_diff, normal) < 0)
					{
						in_shape = false;
						//break;
					}
				}

				if (mp->totloop > 0) {
					mul_v3_fl(centr, 1.0f/(float)mp->totloop);
				}
				mul_m4_v3(obb->obmat, centr);

				if (in_shape)
				{
					in_shape = false;
					sub_v3_v3v3(centr_diff, centr, obb->loc);
					normalize_v3_v3(normalized, centr_diff);
					dot = dot_v3v3(normalized, normalized_mi_diff);
					len = len_v3(centr_diff);

					//print_v3("CENTR", centr_diff);
					//print_v3("MI", mi_diff);
					//dot2 = dot_v3v3(centr_diff, mi_diff);

					if ((len > len_mi * dot) && (dot > 0))
					{
						//printf("LEN_MI, LEN, DOT, DOT2 (OK) : %f %f %f %f\n", len_mi, len, dot, dot2);
						in_shape = true;
					}
				}

				if (in_shape)
					break;
			}
		}

		if (in_shape)
		{
			mi->cluster_index = index < seed_count ? seeds[index]->cluster_index : index;
		}
	}

	BLI_kdtree_free(tree);
	MEM_freeN(seeds);
}

static void do_clusters(FractureModifierData *fmd, Object* obj)
{
	/*grow clusters from all meshIslands */
	do_cluster_count(fmd, obj);
}

/* used for external python api */
RigidBodyShardCon *BKE_fracture_mesh_constraint_create(Scene *scene, FractureModifierData *fmd,
													 Shard *mi1, Shard *mi2, short con_type)
{
	RigidBodyShardCon *rbsc;

	if (mi1 == NULL || mi2 == NULL)
		return NULL;

	rbsc = BKE_rigidbody_create_shard_constraint(scene, con_type, false);
	rbsc->mi1 = mi1;
	rbsc->mi2 = mi2;

	/* moved default meshconstraint pos calculation here to creation, so you can override it later on*/
	/* do this for all constraints */
	/* location for fixed constraints doesnt matter, so keep old setting */
	if (rbsc->type == RBC_TYPE_FIXED) {
		copy_v3_v3(rbsc->pos, rbsc->mi1->rigidbody->pos);
	}
	else {
		/* else set location to center */
		add_v3_v3v3(rbsc->pos, rbsc->mi1->rigidbody->pos, rbsc->mi2->rigidbody->pos);
		mul_v3_fl(rbsc->pos, 0.5f);
	}

	copy_qt_qt(rbsc->orn, rbsc->mi1->rigidbody->orn);

	BLI_addtail(&fmd->shared->constraints, rbsc);

	/* store constraints per meshisland too, to allow breaking percentage */
	if (mi1->participating_constraints == NULL) {
		mi1->participating_constraints = MEM_mallocN(sizeof(RigidBodyShardCon *), "part_constraints_mi1");
		mi1->participating_constraints[0] = rbsc;
		mi1->participating_constraint_count = 1;
	}
	else
	{
		mi1->participating_constraints = MEM_reallocN(mi1->participating_constraints, sizeof(RigidBodyShardCon *) * (mi1->participating_constraint_count + 1));
		mi1->participating_constraints[mi1->participating_constraint_count] = rbsc;
		mi1->participating_constraint_count++;
	}

	if (mi2->participating_constraints == NULL) {
		mi2->participating_constraints = MEM_mallocN(sizeof(RigidBodyShardCon *), "part_constraints_mi2");
		mi2->participating_constraints[0] = rbsc;
		mi2->participating_constraint_count = 1;
	}
	else
	{
		mi2->participating_constraints = MEM_reallocN(mi2->participating_constraints, sizeof(RigidBodyShardCon *) * (mi2->participating_constraint_count + 1));
		mi2->participating_constraints[mi2->participating_constraint_count] = rbsc;
		mi2->participating_constraint_count++;
	}

	return rbsc;
}

void BKE_fracture_constraints_refresh(FractureModifierData *fmd, Object *ob, Scene* scene)
{
	double start = PIL_check_seconds_timer();

	do_clusters(fmd, ob);
	printf("Clustering done, %g\n", PIL_check_seconds_timer() - start);

	start = PIL_check_seconds_timer();

	if (fmd->flag & MOD_FRACTURE_USE_CONSTRAINTS) {
		//if (fmd->dynamic_new_constraints != MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS) {
		create_constraints(fmd, ob, scene);
	}

	printf("Building constraints done, %g\n", PIL_check_seconds_timer() - start);
	printf("Constraints: %d\n", BLI_listbase_count(&fmd->shared->constraints));
}

void BKE_fracture_constraints_free(FractureModifierData *fmd, RigidBodyWorld *rbw)
{
	Shard *mi = NULL;
	RigidBodyShardCon *rbsc = NULL;

	for (mi = fmd->shared->shards.first; mi; mi = mi->next) {
		if (mi->participating_constraints != NULL && mi->participating_constraint_count > 0) {
			int i;
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
	}

	while (fmd->shared->constraints.first) {
		rbsc = fmd->shared->constraints.first;
		BLI_remlink(&fmd->shared->constraints, rbsc);
		if (rbw)
			BKE_rigidbody_remove_shard_con(rbw, rbsc);
		MEM_freeN(rbsc);
		rbsc = NULL;
	}

	fmd->shared->constraints.first = NULL;
	fmd->shared->constraints.last = NULL;
}

