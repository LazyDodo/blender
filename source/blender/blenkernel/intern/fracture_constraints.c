
#include "MEM_guardedalloc.h"

#include "BKE_fracture.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
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
#include "DNA_group_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

#include "PIL_time.h"

static void search_tree_based(FractureModifierData *rmd, MeshIsland *mi, MeshIsland **meshIslands, KDTree **combined_tree,
							  float co[3], Object *ob, Scene *scene);

static void remove_participants(RigidBodyShardCon* con, MeshIsland *mi);


static int prepareConstraintSearch(FractureModifierData *rmd, MeshIsland ***mesh_islands, KDTree **combined_tree, Object *ob,
								   MVert** mverts)
{
	MeshIsland *mi;
	int i = 0, ret = 0;
	int islands = 0;

	if (rmd->dm_group && rmd->use_constraint_group)
	{
		CollectionObject *go;

		for (go = rmd->dm_group->gobject.first; go; go = go->next)
		{
			if ( go->ob != ob)
			{
				FractureModifierData *fmdi = (FractureModifierData *)modifiers_findByType(go->ob, eModifierType_Fracture);
				if (fmdi) {
					islands += BLI_listbase_count(&fmdi->shared->meshIslands);
				}
			}
		}
	}
	else {

		islands = BLI_listbase_count(&rmd->shared->meshIslands);
	}

	*mesh_islands = MEM_reallocN(*mesh_islands, islands * sizeof(MeshIsland *));

	if (rmd->dm_group && rmd->use_constraint_group)
	{
		CollectionObject *go;
		int j = 0;

		for (go = rmd->dm_group->gobject.first; go; go = go->next)
		{
			if ( go->ob != ob)
			{
				FractureModifierData *fmdi = (FractureModifierData *)modifiers_findByType(go->ob, eModifierType_Fracture);
				if (fmdi) {
					for (mi = fmdi->shared->meshIslands.first; mi; mi = mi->next) {
						mi->object_index = j;
						(*mesh_islands)[i] = mi;
						i++;
					}
				}

				j++;
			}
		}
	}
	else {

		for (mi = rmd->shared->meshIslands.first; mi; mi = mi->next) {
			mi->object_index = -1;
			(*mesh_islands)[i] = mi;
			i++;
		}
	}

	if (rmd->constraint_target == MOD_FRACTURE_CENTROID)
	{
		*combined_tree = BLI_kdtree_new(islands);
		for (i = 0; i < islands; i++) {
			float obj_centr[3];
			mul_v3_m4v3(obj_centr, ob->obmat, (*mesh_islands)[i]->centroid);
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

		if (rmd->dm_group && rmd->use_constraint_group)
		{
			CollectionObject *go;
			mvert = MEM_mallocN(sizeof(MVert), "mvert");

			for (go = rmd->dm_group->gobject.first; go; go = go->next)
			{
				if ( go->ob != ob)
				{
					float imat[4][4];
					FractureModifierData *fmdi = (FractureModifierData *)modifiers_findByType(go->ob, eModifierType_Fracture);
					if (fmdi && fmdi->shared->visible_mesh_cached) {
						int v = fmdi->shared->visible_mesh_cached->totvert;
						int x = 0;

						//invert_m4_m4(imat, go->ob->obmat);
						copy_m4_m4(imat, go->ob->obmat);

						mvert = MEM_reallocN(mvert, sizeof(MVert) * (totvert + v));
						memcpy(mvert + totvert, fmdi->shared->visible_mesh_cached->mvert, v * sizeof(MVert));

						for (x = totvert; x < totvert + v; x++)
						{
							mul_v3_m4v3(mvert[x].co, imat, mvert[x].co);
						}

						totvert += v;
					}
				}
			}
		}
		else if (rmd && rmd->shared->visible_mesh_cached) {
			totvert = rmd->shared->visible_mesh_cached->totvert;
			mvert = rmd->shared->visible_mesh_cached->mvert;
		}

		if (totvert > 0)
		{

			*combined_tree = BLI_kdtree_new(totvert);
			for (i = 0, mv = mvert; i < totvert; i++, mv++) {
				float co[3];
				if (rmd->dm_group && rmd->use_constraint_group)
				{
					copy_v3_v3(co, mv->co);
				}
				else {
					mul_v3_m4v3(co, ob->obmat, mv->co);
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
	MeshIsland **mesh_islands = MEM_mallocN(sizeof(MeshIsland *), "mesh_islands");
	int count, i = 0;
	MeshIsland *mi;
	MVert *mvert = NULL;

	float max_mass = 0.0f;
	for (mi = rmd->shared->meshIslands.first; mi; mi = mi->next)
	{
		if (mi->rigidbody->mass > max_mass)
			max_mass = mi->rigidbody->mass;
	}

	count = prepareConstraintSearch(rmd, &mesh_islands, &coord_tree, ob, &mvert);

	for (i = 0; i < count; i++) {
		if (rmd->constraint_target == MOD_FRACTURE_CENTROID) {
			search_tree_based(rmd, mesh_islands[i], mesh_islands, &coord_tree, NULL, ob, scene);
		}
		else if (rmd->constraint_target == MOD_FRACTURE_VERTEX) {
			MeshIsland *mii = NULL;
			mii = BLI_ghash_lookup(rmd->shared->vertex_island_map, SET_INT_IN_POINTER(i));
			search_tree_based(rmd, mii, mesh_islands, &coord_tree, mvert[i].co, ob, scene);
		}
	}

	if (coord_tree != NULL) {
		BLI_kdtree_free(coord_tree);
		coord_tree = NULL;
	}

	MEM_freeN(mesh_islands);

	if (rmd->dm_group && rmd->use_constraint_group && mvert)
	{	//was copied from modifiers... so remove now
		MEM_freeN(mvert);
	}
}

static void search_tree_based(FractureModifierData *rmd, MeshIsland *mi, MeshIsland **meshIslands,
							  KDTree **combined_tree, float co[3], Object *ob, Scene *scene)
{
	int r = 0, limit = 0, i = 0;
	KDTreeNearest *n3 = NULL;
	float dist, obj_centr[3];

	limit = rmd->constraint_limit;
	dist = rmd->contact_dist;
	//factor = rmd->mass_threshold_factor;

	if ((rmd->fracture_mode == MOD_FRACTURE_DYNAMIC) &&
		(rmd->dynamic_new_constraints != MOD_FRACTURE_ALL_DYNAMIC_CONSTRAINTS))
	{
		Shard* s = BKE_fracture_shard_find(&rmd->shared->frac_mesh->shard_map, mi->id);
		if (s->parent_id > -1) {
			return;
		}
	}

	if (rmd->constraint_target == MOD_FRACTURE_CENTROID) {
		mul_v3_m4v3(obj_centr, ob->obmat, mi->centroid);
	}
	else if (rmd->constraint_target == MOD_FRACTURE_VERTEX){
		if (!(rmd->dm_group && rmd->use_constraint_group))
		{
			mul_v3_m4v3(obj_centr, ob->obmat, co);
		}
		else {
			copy_v3_v3(obj_centr, co);
		}
	}

	r = BLI_kdtree_range_search(*combined_tree, obj_centr, &n3, dist);

	/* use centroid dist based approach here, together with limit */
	for (i = 0; i < r; i++) {
		MeshIsland *mi2 = NULL;

		if (rmd->constraint_target == MOD_FRACTURE_CENTROID) {
			mi2 = meshIslands[(n3 + i)->index];
		}
		else if(rmd->constraint_target == MOD_FRACTURE_VERTEX) {
			int index = (n3 + i)->index;
			mi2 = BLI_ghash_lookup(rmd->shared->vertex_island_map, SET_INT_IN_POINTER(index));
		}
		if ((mi != mi2) && (mi2 != NULL)) {
			float thresh = rmd->breaking_threshold;
			int con_type = rmd->use_compounds ? RBC_TYPE_COMPOUND : rmd->constraint_type;

			if ((i >= limit) && (limit > 0)) {
				break;
			}

			if ((rmd->fracture_mode == MOD_FRACTURE_DYNAMIC))
			{
				Shard* s1 = BKE_fracture_shard_find(&rmd->shared->frac_mesh->shard_map, mi->id);
				Shard* s2 = BKE_fracture_shard_find(&rmd->shared->frac_mesh->shard_map, mi2->id);

				if (rmd->dynamic_new_constraints == MOD_FRACTURE_MIXED_DYNAMIC_CONSTRAINTS) {
					//only build between old and new
					if (s1->parent_id > -1 && s2->parent_id > -1) {
						continue;
					}
				}
				else if (rmd->dynamic_new_constraints == MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS){
					// dont build at all
					if (s2->parent_id > -1) {
						continue;
					}
				}
			}

			BKE_fracture_constraint_create(scene, rmd, mi, mi2, con_type, thresh);
		}
	}

	if (n3 != NULL) {
		MEM_freeN(n3);
		n3 = NULL;
	}
}

void BKE_fracture_meshislands_connect(Scene* sc, FractureModifierData *fmd, MeshIsland *mi1, MeshIsland *mi2,
											 int con_type, float thresh)
{
	int con_found = false;
	RigidBodyShardCon *con;
	bool ok = mi1 && mi1->rigidbody;
	ok = ok && mi2 && mi2->rigidbody;
	ok = ok && fmd->use_constraints;
	ok = ok && ((!(fmd->dm_group && fmd->use_constraint_group) && (mi1->object_index == -1) && (mi2->object_index == -1))||
		 (fmd->dm_group && fmd->use_constraint_group && (mi1->object_index != mi2->object_index)));

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
		BKE_fracture_meshislands_connect(sc, fmd, mi1, mi2, con_type, thresh);
	}
}

void BKE_fracture_mesh_constraint_remove(FractureModifierData *fmd, RigidBodyShardCon* con, Scene *scene)
{
	remove_participants(con, con->mi1);
	remove_participants(con, con->mi2);

	BLI_remlink(&fmd->shared->meshConstraints, con);
	BKE_rigidbody_remove_shard_con(scene, con);
	MEM_freeN(con);
	if (fmd->constraint_count > 0)
	{
		fmd->constraint_count--;
	}
}

void BKE_fracture_mesh_constraint_remove_all(FractureModifierData *fmd, Scene *scene)
{
	BKE_fracture_constraints_free(fmd, scene);
	fmd->constraint_count = 0;
}

static void remove_participants(RigidBodyShardCon* con, MeshIsland *mi)
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

void BKE_fracture_constraint_create(Scene* scene, FractureModifierData* fmd, MeshIsland *mi1, MeshIsland *mi2, short con_type, float thresh)
{
	RigidBodyShardCon *rbsc;
	rbsc = BKE_rigidbody_create_shard_constraint(scene, con_type, fmd->fracture_mode != MOD_FRACTURE_DYNAMIC);
	rbsc->mi1 = mi1;
	rbsc->mi2 = mi2;

	BLI_snprintf(rbsc->name, 64, "%s-%s", rbsc->mi1->name, rbsc->mi2->name);

	if (thresh == 0 || fmd->use_breaking == false) {
		rbsc->flag &= ~RBC_FLAG_USE_BREAKING;
	}

	if (!fmd->use_constraint_collision) {
		rbsc->flag |= RBC_FLAG_DISABLE_COLLISIONS;
	}
	else {
		rbsc->flag &= ~RBC_FLAG_DISABLE_COLLISIONS;
	}

	if ((mi1->particle_index != -1) && (mi2->particle_index != -1) &&
		(mi1->particle_index == mi2->particle_index))
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
		if ((mi1->particle_index != -1) && (mi2->particle_index != -1) &&
			(mi1->particle_index != mi2->particle_index))
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

	BLI_addtail(&fmd->shared->meshConstraints, rbsc);

	if ((mi1->object_index == -1) && (mi2->object_index == -1))
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
	MeshIsland *mi, **seeds;
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
	for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next ) {
		mi->particle_index = -1;
	}

	mi_count = BLI_listbase_count(&fmd->shared->meshIslands);
	seed_count = (fmd->cluster_count > mi_count ? mi_count : fmd->cluster_count);
	//seed_count = fmd->cluster_count;

	if (fmd->cluster_group)
	{
		 group_count = BLI_listbase_count(&fmd->cluster_group->gobject);
	}

	seeds = MEM_mallocN(sizeof(MeshIsland *) * (seed_count), "seeds");
	tree = BLI_kdtree_new(seed_count + group_count);

	/* pick n seed locations, randomly scattered over the object */
	for (k = 0; k < seed_count; k++) {
		int which_index = k * (int)(mi_count / seed_count);
		MeshIsland *which = (MeshIsland *)BLI_findlink(&fmd->shared->meshIslands, which_index);
		which->particle_index = k;
		print_v3("INSERT", which->centroid);
		BLI_kdtree_insert(tree, k, which->centroid);
		seeds[k] = which;
	}

	/*add the group here */
	if (fmd->cluster_group) {
		for (k = seed_count, go = fmd->cluster_group->gobject.first; go; k++, go = go->next)
		{
			float loc[3];

			mul_v3_m4v3(loc, mat, go->ob->loc);

			print_v3("INSERT", loc);
			BLI_kdtree_insert(tree, k, loc);
		}
	}

	BLI_kdtree_balance(tree);

	/* assign each shard to its closest center */
	for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next ) {
		KDTreeNearest n;
		int index;

		index = BLI_kdtree_find_nearest(tree, mi->centroid, &n);
		mi->particle_index = index < seed_count ? seeds[index]->particle_index : index;
	}

	BLI_kdtree_free(tree);
	MEM_freeN(seeds);
}

static void do_clusters(FractureModifierData *fmd, Object* obj)
{
	/*grow clusters from all meshIslands */
	do_cluster_count(fmd, obj);
}

RigidBodyShardCon *BKE_fracture_mesh_constraint_create(Scene *scene, FractureModifierData *fmd,
													 MeshIsland *mi1, MeshIsland *mi2, short con_type)
{
	RigidBodyShardCon *rbsc;

	if (mi1 == NULL || mi2 == NULL)
		return NULL;

	rbsc = BKE_rigidbody_create_shard_constraint(scene, con_type, false);
	rbsc->mi1 = mi1;
	rbsc->mi2 = mi2;

	if (fmd->fracture_mode == MOD_FRACTURE_EXTERNAL)
	{
		/* disable breaking flag here by default, only enable later via python if necessary */
		rbsc->flag &= ~RBC_FLAG_USE_BREAKING;

		/* also delete all other "default" flags here, let them being overriden from python too */
		//rbsc->flag &= ~RBC_FLAG_ENABLED;
		rbsc->flag |= RBC_FLAG_DISABLE_COLLISIONS;

#if 0
		/* and dont allow to let constrained objects collide per default, as with regular constraints */
		rbsc->flag |= RBC_FLAG_DISABLE_COLLISIONS;
#endif
	}

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

	if (BLI_listbase_is_empty(&fmd->shared->meshConstraints))
	{
		fmd->constraint_count = 0;
	}

	BLI_addtail(&fmd->shared->meshConstraints, rbsc);
	fmd->constraint_count++;

#if 0
	if (index > -1)
	{
		rbsc->id = index;
	}
	else
	{
		rbsc->id = fmd->constraint_count-1;
	}
#endif

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

	if (fmd->use_constraints) {
		int count = 0;

		/* fire a callback which can then load external constraint data right NOW */
	   // BLI_callback_exec(G.main, &ob->id, BLI_CB_EVT_FRACTURE_CONSTRAINT_REFRESH);

		/*if we loaded constraints, dont create other ones now */
		count = BLI_listbase_count(&fmd->shared->meshConstraints);

		if (count == 0 || fmd->dynamic_new_constraints != MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS) {
			create_constraints(fmd, ob, scene); /* check for actually creating the constraints inside*/
		}
	}

	printf("Building constraints done, %g\n", PIL_check_seconds_timer() - start);
	printf("Constraints: %d\n", BLI_listbase_count(&fmd->shared->meshConstraints));
}

void BKE_fracture_constraints_free(FractureModifierData *fmd, Scene *scene)
{
	MeshIsland *mi = NULL;
	RigidBodyShardCon *rbsc = NULL;

	//hmm after loading the pointers might be out of sync...
	if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC) {
		if (fmd->shared->current_mi_entry) {
			fmd->shared->meshIslands = fmd->shared->current_mi_entry->meshIslands;
		}
		else {
			fmd->shared->meshIslands.first = NULL;
			fmd->shared->meshIslands.last = NULL;
		}
	}

	for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
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

	while (fmd->shared->meshConstraints.first) {
		rbsc = fmd->shared->meshConstraints.first;
		BLI_remlink(&fmd->shared->meshConstraints, rbsc);

		if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
		{
			BKE_rigidbody_remove_shard_con(scene, rbsc);
		}
		MEM_freeN(rbsc);
		rbsc = NULL;
	}

	fmd->shared->meshConstraints.first = NULL;
	fmd->shared->meshConstraints.last = NULL;
}

