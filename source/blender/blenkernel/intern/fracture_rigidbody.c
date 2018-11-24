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

#include "BLI_utildefines.h"
#include "BLI_kdtree.h"
#include "BLI_math.h"
#include "BLI_listbase.h"
#include "BLI_threads.h"

#include "DNA_defs.h"
#include "DNA_object_types.h"
#include "DNA_modifier_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_scene_types.h"
#include "DNA_ID.h"
#include "DNA_fracture_types.h"

#include "DEG_depsgraph_build.h"
#include "DEG_depsgraph_query.h"

#ifdef WITH_BULLET
#include "RBI_api.h"
#endif

#include "BKE_collection.h"
#include "BKE_fracture.h"
#include "BKE_rigidbody.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"
#include "BKE_object.h"
#include "BKE_modifier.h"
#include "BKE_global.h"
#include "BKE_pointcache.h"

/*====================================================================================================================*/

/* Fracture Modifier stuff */

#ifdef WITH_BULLET

/* Fracture Modifier related prototypes */

//static void resetDynamic(RigidBodyWorld *rbw, bool do_reset_always, bool override_bind, struct Depsgraph *depsgraph);
static void check_fracture(rbContactPoint *cp, Scene *scene);
static void test_deactivate_rigidbody(RigidBodyOb *rbo, Shard *mi);
static float box_volume(float size[3]);


void BKE_rigidbody_activate(RigidBodyOb* rbo, RigidBodyWorld *UNUSED(rbw), Shard *mi, Object *ob)
{
	RigidBodyShardCon *con;
	int i;

	if( !BKE_rigidbody_activate_by_size_check(ob, mi))
	{
		return;
	}

	if (rbo->flag & RBO_FLAG_KINEMATIC && rbo->type == RBO_TYPE_ACTIVE)
	{
		rbo->flag &= ~RBO_FLAG_KINEMATIC;
		rbo->flag |= RBO_FLAG_NEEDS_VALIDATE;

		//propagate trigger on impact / activation
		if (ob->rigidbody_object->flag & RBO_FLAG_PROPAGATE_TRIGGER) {
			rbo->flag |= RBO_FLAG_PROPAGATE_TRIGGER;
		}

		//RB_dworld_remove_body(rbw->physics_world, rbo->physics_object);
		RB_body_set_mass(rbo->shared->physics_object, RBO_GET_MASS(rbo));
		RB_body_set_kinematic_state(rbo->shared->physics_object, false);
		//RB_dworld_add_body(rbw->physics_world, rbo->physics_object, rbo->col_groups, mi, ob);
		RB_body_activate(rbo->shared->physics_object);
	}

	if (mi && ob->rigidbody_object->flag & RBO_FLAG_CONSTRAINT_DISSOLVE) {
		for (i = 0; i < mi->participating_constraint_count; i++) {
			bool different_cluster = false;
			bool dissolve_plastic = (ob->rigidbody_object->flag & RBO_FLAG_PLASTIC_DISSOLVE);
			con = mi->participating_constraints[i];

			different_cluster = ((con->mi1->cluster_index != con->mi2->cluster_index) ||
								((con->mi1->cluster_index == -1) && (con->mi2->cluster_index == -1)));

			if (con->physics_constraint && different_cluster) {
				if (dissolve_plastic) {
					con->flag |= RBC_FLAG_PLASTIC_ACTIVE;
				}

				if (con->breaking_threshold >= 0)
				{
					RB_constraint_set_enabled(con->physics_constraint, false);
				}
			}
		}
	}
}

bool BKE_rigidbody_modifier_active(FractureModifierData *rmd)
{
	return ((rmd != NULL) && (rmd->modifier.mode & (eModifierMode_Realtime | eModifierMode_Render)) //&&
			/*(rmd->shared->refresh == false || rmd->fracture_mode == MOD_FRACTURE_DYNAMIC)*/);
}

static void calc_dist_angle(RigidBodyShardCon *con, float *dist, float *angle, bool exact)
{
	float q1[4], q2[4], qdiff[4], axis[3];

	if (con == NULL || con->mi1 == NULL || con->mi2 == NULL ||
		con->mi1->rigidbody == NULL || con->mi2->rigidbody == NULL)
	{
		*dist = 0;
		*angle = 0;
		return;
	}

	sub_v3_v3v3(axis, con->mi1->rigidbody->pos, con->mi2->rigidbody->pos);
	*dist = len_v3(axis);

	copy_qt_qt(q1, con->mi1->rigidbody->orn);
	copy_qt_qt(q2, con->mi2->rigidbody->orn);

	if (exact)
	{
		float iquat1[4], iquat2[4];
		invert_qt_qt(iquat1, con->mi1->rot);
		invert_qt_qt(iquat2, con->mi2->rot);
		mul_qt_qtqt(q1, q1, iquat1);
		mul_qt_qtqt(q2, q2, iquat2);
		rotation_between_quats_to_quat(qdiff, q1, q2);
		normalize_qt(qdiff);
		*angle = 2.0f * saacos(qdiff[0]);
		if (!isfinite(*angle)) {
			*angle = 0.0f;
		}
	}
	else
	{
		//XXX TODO probably very wrong here
		invert_qt(q1);
		mul_qt_qtqt(qdiff, q1, q2);
		quat_to_axis_angle(axis, angle, qdiff);
	}
}

void BKE_rigidbody_start_dist_angle(RigidBodyShardCon *con, bool exact, bool both)
{
	/* store starting angle and distance per constraint*/
	float dist, angle;
	calc_dist_angle(con, &dist, &angle, exact);

	if (both)
	{
		con->start_dist = dist;
		con->start_angle = angle;
		//printf("Start Values(dist, angle) %f %f %f %f\n", con->start_dist, con->start_angle, dist, angle);
	}

	con->start_dist_deform = dist;
	con->start_angle_deform = angle;
}

float BKE_rigidbody_calc_max_con_mass(Object *ob)
{
	FractureModifierData *rmd;
	ModifierData *md;
	RigidBodyShardCon *con;
	float max_con_mass = 0, con_mass;

	for (md = ob->modifiers.first; md; md = md->next) {
		if (md->type == eModifierType_Fracture) {
			rmd = (FractureModifierData *)md;
			for (con = rmd->shared->constraints.first; con; con = con->next) {
				if ((con->mi1 != NULL && con->mi1->rigidbody != NULL) &&
					(con->mi2 != NULL && con->mi2->rigidbody != NULL)) {
					con_mass = con->mi1->rigidbody->mass + con->mi2->rigidbody->mass;
					if (con_mass > max_con_mass) {
						max_con_mass = con_mass;
					}
				}
			}

			return max_con_mass;
		}
	}

	return 0;
}

float BKE_rigidbody_calc_min_con_dist(Object *ob)
{
	FractureModifierData *rmd;
	ModifierData *md;
	RigidBodyShardCon *con;
	float min_con_dist = FLT_MAX, con_dist, con_vec[3];

	for (md = ob->modifiers.first; md; md = md->next) {
		if (md->type == eModifierType_Fracture) {
			rmd = (FractureModifierData *)md;
			for (con = rmd->shared->constraints.first; con; con = con->next) {
				if ((con->mi1 != NULL && con->mi1->rigidbody != NULL) &&
					(con->mi2 != NULL && con->mi2->rigidbody != NULL)) {
					sub_v3_v3v3(con_vec, con->mi1->loc, con->mi2->loc);
					con_dist = len_v3(con_vec);
					if (con_dist < min_con_dist) {
						min_con_dist = con_dist;
					}
				}
			}

			return min_con_dist;
		}
	}

	return FLT_MAX;
}


void BKE_rigidbody_calc_threshold(float max_con_mass, FractureModifierData *rmd, RigidBodyShardCon *con) {

	float max_thresh, thresh = 0.0f, con_mass;
	if ((max_con_mass == 0) && (rmd->flag & MOD_FRACTURE_USE_MASS_DEP_THRESHOLDS)) {
		return;
	}

	if ((con->mi1 == NULL) || (con->mi2 == NULL)) {
		return;
	}

	max_thresh = thresh = rmd->breaking_threshold;
	if ((con->mi1->rigidbody != NULL) && (con->mi2->rigidbody != NULL)) {

		con_mass = con->mi1->rigidbody->mass + con->mi2->rigidbody->mass;
		if (rmd->flag & MOD_FRACTURE_USE_MASS_DEP_THRESHOLDS)
		{
			thresh = (con_mass / max_con_mass) * max_thresh;
		}

		con->breaking_threshold = thresh;
	}
}

static float box_volume(float size[3])
{
	float volume;

	volume = size[0] * size[1] * size[2];
	if (volume == 0 && size[0] == 0) {
		volume = size[1] * size[2];
	}
	else if (volume == 0 && size[1] == 0) {
		volume = size[0] * size[2];
	}
	else if (volume == 0 && size[2] == 0) {
		volume = size[0] * size[1];
	}

	return volume;
}

/* helper function to calculate volume of rigidbody object */
float BKE_rigidbody_calc_volume_dm(Mesh *dm, RigidBodyOb *rbo, Object* ob)
{
	float loc[3]  = {0.0f, 0.0f, 0.0f};
	float size[3]  = {1.0f, 1.0f, 1.0f};
	float scale[3] = {1.0f, 1.0f, 1.0f};
	float radius = 1.0f;
	float height = 1.0f;

	float volume = 0.0f;

	/* if automatically determining dimensions, use the Object's boundbox
	 *	- assume that all quadrics are standing upright on local z-axis
	 *	- assume even distribution of mass around the Object's pivot
	 *	  (i.e. Object pivot is centralised in boundbox)
	 *	- boundbox gives full width
	 */
	/* XXX: all dimensions are auto-determined now... later can add stored settings for this*/
	BKE_fracture_mesh_boundbox_calc(dm, loc, size);
	mul_v3_fl(size, 2.0f);

	/* also take object scale into account */
	if (ob)
	{
		mat4_to_size(scale, ob->obmat);
		mul_v3_v3(size, scale);
	}

	if (ELEM(rbo->shape, RB_SHAPE_CAPSULE, RB_SHAPE_CYLINDER, RB_SHAPE_CONE)) {
		/* take radius as largest x/y dimension, and height as z-dimension */
		radius = MAX2(size[0], size[1]) * 0.5f;
		height = size[2];
	}
	else if (rbo->shape == RB_SHAPE_SPHERE) {
		/* take radius to the the largest dimension to try and encompass everything */
		radius = max_fff(size[0], size[1], size[2]) * 0.5f;
	}

	/* calculate volume as appropriate  */
	switch (rbo->shape) {

		case RB_SHAPE_SPHERE:
			volume = 4.0f / 3.0f * (float)M_PI * radius * radius * radius;
			break;

		/* for now, assume that capsule is close enough to a cylinder... */
		case RB_SHAPE_CAPSULE:
		case RB_SHAPE_CYLINDER:
			volume = (float)M_PI * radius * radius * height;
			break;

		case RB_SHAPE_CONE:
			volume = (float)M_PI / 3.0f * radius * radius * height;
			break;

		/* for now, all mesh shapes are just treated as boxes...
		 * NOTE: this may overestimate the volume, but other methods are overkill
		 */
		case RB_SHAPE_BOX:
			volume = box_volume(size);
			break;

		case RB_SHAPE_CONVEXH:
		case RB_SHAPE_TRIMESH:
		{
			if (!ob) {
				/* for quick shard mass (re)calculation approximate,
				 * else for many shards on complex geometries this becomes
				 * very slow */
				volume = box_volume(size);
			}
			else
			{
				//for initial mass calculation take exact values
				MVert *mvert = dm->mvert;
				int totvert = dm->totvert;
				MLoop *mloop = dm->mloop;
				MLoopTri *mlooptri = BKE_mesh_runtime_looptri_ensure(dm);
				int tottri = BKE_mesh_runtime_looptri_len(dm);

				BKE_mesh_calc_volume(mvert, totvert, mlooptri, tottri, mloop, &volume, NULL);

				if (volume == 0.0f)
					volume = box_volume(size);
			}
			break;
		}

#if 0 // XXX: not defined yet
		case RB_SHAPE_COMPOUND:
			volume = 0.0f;
			break;
#endif
	}

	/* return the volume calculated */
	return volume;
}

void BKE_rigidbody_calc_shard_mass(Object *ob, Shard *mi)
{
	Mesh *dm_ob = ob->runtime.mesh_eval, *dm_mi = mi->mesh;
	float vol_mi = 0, mass_mi = 0, vol_ob = 0, mass_ob = 0;

	if (dm_ob == NULL) {
		/* fallback method */
		if (ob->type == OB_MESH) {
			/* if we have a mesh, determine its volume */
			dm_ob = ob->data;
			vol_ob = BKE_rigidbody_calc_volume_dm(dm_ob, ob->rigidbody_object, ob);
		}
		else {
			/* else get object boundbox as last resort */
			float dim[3];
			BKE_object_dimensions_get(ob, dim);
			vol_ob = dim[0] * dim[1] * dim[2];
		}
	}
	else
	{
		vol_ob = BKE_rigidbody_calc_volume_dm(dm_ob, ob->rigidbody_object, ob);
	}

	mass_ob = ob->rigidbody_object->mass;

	if (vol_ob > 0) {
		dm_mi = mi->mesh;
		vol_mi = BKE_rigidbody_calc_volume_dm(dm_mi, mi->rigidbody, NULL);
		mass_mi = (vol_mi / vol_ob) * mass_ob;
		mi->rigidbody->mass = mass_mi;
	}

	if (mi->rigidbody->type == RBO_TYPE_ACTIVE) {
		if (mi->rigidbody->mass == 0)
			mi->rigidbody->mass = 0.001;  /* set a minimum mass for active objects */
	}

	/* only active bodies need mass update */
	if ((mi->rigidbody->shared->physics_object) && (mi->rigidbody->type == RBO_TYPE_ACTIVE)) {
		RB_body_set_mass(mi->rigidbody->shared->physics_object, RBO_GET_MASS(mi->rigidbody));
	}
}


/* --------------------- */

/* Create new physics sim collision shape for object and store it,
 * or remove the existing one first and replace...
 */
void BKE_rigidbody_validate_sim_shard_shape(Shard *mi, Object *ob, short rebuild)
{
	Mesh* me_phys = NULL;
	MVert *mv = NULL;
	RigidBodyOb *rbo = mi->rigidbody;
	rbCollisionShape *new_shape = NULL;
	float size[3] = {1.0f, 1.0f, 1.0f}, loc[3] = {0.0f, 0.0f, 0.0f};
	float radius = 1.0f;
	float height = 1.0f;
	float capsule_height;
	float hull_margin = 0.0f;
	bool can_embed = true;
	bool has_volume;
	float min[3], max[3], margin;
	int v = 0;

	/* sanity check */
	if (rbo == NULL)
		return;

	/* don't create a new shape if we already have one and don't want to rebuild it */
	if (rbo->shared->physics_shape && !rebuild)
		return;

	/* if automatically determining dimensions, use the Object's boundbox
	 *	- assume that all quadrics are standing upright on local z-axis
	 *	- assume even distribution of mass around the Object's pivot
	 *	  (i.e. Object pivot is centralized in boundbox)
	 *	- boundbox gives full width
	 */
	// XXX: all dimensions are auto-determined now... later can add stored settings for this
	/* get object dimensions without scaling */

	/*correct physics mesh location, must be around zero */
	me_phys = BKE_fracture_mesh_copy(mi->mesh, ob);
	for (v = 0, mv = me_phys->mvert; v < me_phys->totvert; mv++, v++)
	{
		sub_v3_v3(mv->co, mi->loc);
	}

	INIT_MINMAX(min, max);
	if (!BKE_mesh_minmax(mi->mesh, min, max)) {
		min[0] = min[1] = min[2] = -1.0f;
		max[0] = max[1] = max[2] = 1.0f;
	}

	mid_v3_v3v3(loc, min, max);
	size[0] = (max[0] - min[0]) / 2.0f;
	size[1] = (max[1] - min[1]) / 2.0f;
	size[2] = (max[2] - min[2]) / 2.0f;

	if (ELEM(rbo->shape, RB_SHAPE_CAPSULE, RB_SHAPE_CYLINDER, RB_SHAPE_CONE)) {
		/* take radius as largest x/y dimension, and height as z-dimension */
		radius = MAX2(size[0], size[1]);
		height = size[2];
	}
	else if (rbo->shape == RB_SHAPE_SPHERE) {

		/* take radius to the largest dimension to try and encompass everything */
		radius = (rbo->flag & RBO_FLAG_USE_MARGIN) ? min_fff(size[0], size[1], size[2]) :
				 max_fff(size[0], size[1], size[2]);
	}

	/* create new shape */
	switch (rbo->shape) {
		case RB_SHAPE_BOX:
			new_shape = RB_shape_new_box(size[0], size[1], size[2]);
			break;

		case RB_SHAPE_SPHERE:
			margin = (rbo->flag & RBO_FLAG_USE_MARGIN) ? rbo->margin : 0.0f;

			new_shape = RB_shape_new_sphere(margin);
			break;

		case RB_SHAPE_CAPSULE:
			capsule_height = (height - radius) * 2.0f;
			new_shape = RB_shape_new_capsule(radius, (capsule_height > 0.0f) ? capsule_height : 0.0f);
			break;
		case RB_SHAPE_CYLINDER:
			new_shape = RB_shape_new_cylinder(radius, height);
			break;

		case RB_SHAPE_CONE:
			new_shape = RB_shape_new_cone(radius, height * 2.0f);
			break;

		case RB_SHAPE_CONVEXH:
			/* try to embed collision margin */
			has_volume = (MIN3(size[0], size[1], size[2]) > 0.0f);

			if (!(rbo->flag & RBO_FLAG_USE_MARGIN) && has_volume)
				hull_margin = 0.04f;
			new_shape = BKE_rigidbody_get_shape_convexhull_from_mesh(me_phys, hull_margin, &can_embed);
			if (!(rbo->flag & RBO_FLAG_USE_MARGIN))
				rbo->margin = (can_embed && has_volume) ? 0.04f : 0.0f;      /* RB_TODO ideally we shouldn't directly change the margin here */
			break;
		case RB_SHAPE_TRIMESH:
		{
			new_shape = BKE_rigidbody_get_shape_trimesh_from_mesh(ob, me_phys);
			break;
		}
		case RB_SHAPE_COMPOUND:
		{
			new_shape = RB_shape_new_compound();
			break;
		}
	}
	/* assign new collision shape if creation was successful */
	if (new_shape) {
		margin = RBO_GET_MARGIN(rbo);
		if (rbo->shared->physics_shape)
			RB_shape_delete(rbo->shared->physics_shape);
		rbo->shared->physics_shape = new_shape;

		RB_shape_set_margin(rbo->shared->physics_shape, margin);
	}
	else { /* otherwise fall back to box shape */
		rbo->shape = RB_SHAPE_BOX;
		BKE_rigidbody_validate_sim_shard_shape(mi, ob, true);
	}

	BKE_fracture_mesh_free(me_phys);
}

/* --------------------- */

/* Create physics sim representation of shard given RigidBody settings
 * < rebuild: even if an instance already exists, replace it
 */
void BKE_rigidbody_validate_sim_shard(RigidBodyWorld *rbw, Shard *mi, Object *ob, FractureModifierData *fmd,
                                      short rebuild, int transfer_speeds, float isize[3])
{
	RigidBodyOb *rbo = (mi) ? mi->rigidbody : NULL;
	float loc[3];
	float rot[4];

	/* sanity checks:
	 *	- object doesn't have RigidBody info already: then why is it here?
	 */
	if (rbo == NULL)
		return;

	if (!rbw || !rbw->shared || !rbw->shared->physics_world) {
		return;
	}

	/* at validation, reset frame count as well */

	/* make sure collision shape exists */
	/* FIXME we shouldn't always have to rebuild collision shapes when rebuilding objects, but it's needed for constraints to update correctly */
	if (rbo->shared->physics_shape == NULL || rebuild)
		BKE_rigidbody_validate_sim_shard_shape(mi, ob, true);

	if (rbo->shared->physics_object) {
		if (rebuild == false /*|| mi->rigidbody->flag & RBO_FLAG_KINEMATIC_REBUILD*/)
			RB_dworld_remove_body(rbw->shared->physics_world, rbo->shared->physics_object);
	}

	if (!rbo->shared->physics_object || rebuild) {
		float size[3];

		/* remove rigid body if it already exists before creating a new one */
		if (rbo->shared->physics_object) {
			RB_body_delete(rbo->shared->physics_object);
		}

		copy_v3_v3(loc, rbo->pos);
		copy_qt_qt(rot, rbo->orn);
		copy_v3_v3(size, isize);

		mul_v3_v3(size, ob->size);
		rbo->shared->physics_object = RB_body_new(rbo->shared->physics_shape, loc, rot);

		RB_body_set_friction(rbo->shared->physics_object, rbo->friction);
		RB_body_set_restitution(rbo->shared->physics_object, rbo->restitution);

		RB_body_set_damping(rbo->shared->physics_object, rbo->lin_damping, rbo->ang_damping);
		RB_body_set_sleep_thresh(rbo->shared->physics_object, rbo->lin_sleep_thresh, rbo->ang_sleep_thresh);
		RB_body_set_activation_state(rbo->shared->physics_object, rbo->flag & RBO_FLAG_USE_DEACTIVATION);

		if (rbo->type == RBO_TYPE_PASSIVE || rbo->flag & RBO_FLAG_START_DEACTIVATED)
			RB_body_deactivate(rbo->shared->physics_object);


		RB_body_set_linear_factor(rbo->shared->physics_object,
								  (ob->protectflag & OB_LOCK_LOCX) == 0,
								  (ob->protectflag & OB_LOCK_LOCY) == 0,
								  (ob->protectflag & OB_LOCK_LOCZ) == 0);
		RB_body_set_angular_factor(rbo->shared->physics_object,
								   (ob->protectflag & OB_LOCK_ROTX) == 0,
								   (ob->protectflag & OB_LOCK_ROTY) == 0,
								   (ob->protectflag & OB_LOCK_ROTZ) == 0);

		RB_body_set_mass(rbo->shared->physics_object, RBO_GET_MASS(rbo));
		RB_body_set_kinematic_state(rbo->shared->physics_object, rbo->flag & RBO_FLAG_KINEMATIC || rbo->flag & RBO_FLAG_DISABLED);

		if (transfer_speeds)
		{
			RB_body_set_linear_velocity(rbo->shared->physics_object, rbo->lin_vel);
			RB_body_set_angular_velocity(rbo->shared->physics_object, rbo->ang_vel);
		}
	}
	else if (mi->rigidbody->flag & RBO_FLAG_KINEMATIC_REBUILD )
	{
		RB_body_deactivate(rbo->shared->physics_object);
		RB_body_set_mass(rbo->shared->physics_object, 0.0f);
		RB_body_set_kinematic_state(rbo->shared->physics_object, true);
	}

	if (rbw && rbw->shared->physics_world && rbo->shared->physics_object)
	{
		RB_dworld_add_body(rbw->shared->physics_world, rbo->shared->physics_object, rbo->col_groups, mi, ob);
	}

	rbo->flag &= ~RBO_FLAG_NEEDS_VALIDATE;
}

/* --------------------- */

Shard* BKE_rigidbody_closest_meshisland_to_point(FractureModifierData* fmd, Object *ob, Object *ob2, Scene* scene)
{
	Shard *mi, **mi_array = NULL;
	KDTree *tree;
	KDTreeNearest *n = NULL;
	int count = 0;
	int index = 0;
	float loc[3], min[3], max[3], vec[3] = {1, 1, 1};
	int i = 0, r = 0;

	count = BLI_listbase_count(&fmd->shared->shards);
	tree = BLI_kdtree_new(count);
	mi_array = MEM_mallocN(sizeof(Shard*) * count, "mi_array find_closest_meshisland");

	for (mi = fmd->shared->shards.first; mi; mi = mi->next) {
		mul_v3_m4v3(loc, ob->obmat, mi->loc);
		BLI_kdtree_insert(tree, i, loc);
		mi_array[i] = mi;
		i++;
	}

	BLI_kdtree_balance(tree);

	index = BLI_kdtree_find_nearest(tree, ob2->loc, n);

	//create "aabb"
	mul_v3_v3(vec, ob2->size);
	sub_v3_v3v3(min, ob2->loc, vec);
	add_v3_v3v3(max, ob2->loc, vec);

	if (index == -1) {
		MEM_freeN(mi_array);
		BLI_kdtree_free(tree);
		return NULL;
	}

	if (index >= count) {
		index = count-1;
	}

	mi = mi_array[index];

	//do a range search and clip against "aabb" of empty scale for additional inner constraint
	r = BLI_kdtree_range_search(tree, ob2->loc, &n, max_fff(UNPACK3(ob2->size)));
	for (i = 0; i < r; i++)
	{
		float co[3];
		copy_v3_v3(co, n[i].co);

		if ((co[0] > min[0] && co[0] < max[0]) &&
			(co[1] > min[1] && co[1] < max[1]) &&
			(co[2] > min[2] && co[2] < max[2]))
		{
			Shard* mi2;
			int ind = n[i].index;

			if (ind >= count) {
				ind = count-1;
			}

			mi2 = mi_array[ind];

			//connect ?
			BKE_fracture_constraint_create(scene, fmd, mi, mi2, fmd->constraint_type, fmd->breaking_threshold);
		}
	}

	if (n) {
		MEM_freeN(n);
	}

	MEM_freeN(mi_array);
	BLI_kdtree_free(tree);

	return mi;
}

static void rigidbody_set_springs_active(RigidBodyShardCon *rbc, bool active)
{
	if (rbc && rbc->physics_constraint && rbc->type == RBC_TYPE_6DOF_SPRING)
	{
		if (active) //XXX TEST purpose only
		{
			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->flag & RBC_FLAG_USE_SPRING_X);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->spring_stiffness_x);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->spring_damping_x);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->flag & RBC_FLAG_USE_SPRING_Y);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->spring_stiffness_y);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->spring_damping_y);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->flag & RBC_FLAG_USE_SPRING_Z);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->spring_stiffness_z);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->spring_damping_z);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->flag & RBC_FLAG_USE_SPRING_ANG_X);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->spring_stiffness_ang_x);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->spring_damping_ang_x);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->flag & RBC_FLAG_USE_SPRING_ANG_Y);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->spring_stiffness_ang_y);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->spring_damping_ang_y);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->flag & RBC_FLAG_USE_SPRING_ANG_Z);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->spring_stiffness_ang_z);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->spring_damping_ang_z);
		}
		else
		{
			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->flag & RBC_FLAG_USE_SPRING_X);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, 0);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->spring_damping_x);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->flag & RBC_FLAG_USE_SPRING_Y);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, 0);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->spring_damping_y);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->flag & RBC_FLAG_USE_SPRING_Z);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, 0);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->spring_damping_z);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->flag & RBC_FLAG_USE_SPRING_ANG_X);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, 0);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->spring_damping_ang_x);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->flag & RBC_FLAG_USE_SPRING_ANG_Y);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, 0);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->spring_damping_ang_y);

			RB_constraint_set_spring_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->flag & RBC_FLAG_USE_SPRING_ANG_Z);
			RB_constraint_set_stiffness_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, 0);
			RB_constraint_set_damping_6dof_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->spring_damping_ang_z);
		}
	}
}

static void rigidbody_create_shard_physics_constraint(FractureModifierData* fmd, Object* ob, RigidBodyShardCon *rbc, RigidBodyWorld *rbw)
{
	float loc[3];
	float rot[4];
	float lin_lower;
	float lin_upper;
	float ang_lower;
	float ang_upper;
	rbRigidBody *rb1;
	rbRigidBody *rb2;

	if (rbc && rbc->mi1 && rbc->mi2)
	{
		rb1 = rbc->mi1->rigidbody->shared->physics_object;
		rb2 = rbc->mi2->rigidbody->shared->physics_object;
	}
	else
	{
		return;
	}

	add_v3_v3v3(rbc->pos, rbc->mi1->rigidbody->pos, rbc->mi2->rigidbody->pos);
	mul_v3_fl(rbc->pos, 0.5f);

	copy_qt_qt(rbc->orn, rbc->mi1->rigidbody->orn);
	copy_v3_v3(loc, rbc->pos);
	copy_qt_qt(rot, rbc->orn);

	if (rb1 && rb2) {
		switch (rbc->type) {
			case RBC_TYPE_POINT:
				rbc->physics_constraint = RB_constraint_new_point(loc, rb1, rb2);
				break;
			case RBC_TYPE_FIXED:
				rbc->physics_constraint = RB_constraint_new_fixed(loc, rot, rb1, rb2);
				break;
			case RBC_TYPE_HINGE:
				rbc->physics_constraint = RB_constraint_new_hinge(loc, rot, rb1, rb2);
				if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_Z) {
					RB_constraint_set_limits_hinge(rbc->physics_constraint, rbc->limit_ang_z_lower, rbc->limit_ang_z_upper);
				}
				else
					RB_constraint_set_limits_hinge(rbc->physics_constraint, 0.0f, -1.0f);
				break;
			case RBC_TYPE_SLIDER:
				rbc->physics_constraint = RB_constraint_new_slider(loc, rot, rb1, rb2);
				if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_X)
					RB_constraint_set_limits_slider(rbc->physics_constraint, rbc->limit_lin_x_lower, rbc->limit_lin_x_upper);
				else
					RB_constraint_set_limits_slider(rbc->physics_constraint, 0.0f, -1.0f);
				break;
			case RBC_TYPE_PISTON:
				rbc->physics_constraint = RB_constraint_new_piston(loc, rot, rb1, rb2);
				if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_X) {
					lin_lower = rbc->limit_lin_x_lower;
					lin_upper = rbc->limit_lin_x_upper;
				}
				else {
					lin_lower = 0.0f;
					lin_upper = -1.0f;
				}
				if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_X) {
					ang_lower = rbc->limit_ang_x_lower;
					ang_upper = rbc->limit_ang_x_upper;
				}
				else {
					ang_lower = 0.0f;
					ang_upper = -1.0f;
				}
				RB_constraint_set_limits_piston(rbc->physics_constraint, lin_lower, lin_upper, ang_lower, ang_upper);
				break;
			case RBC_TYPE_6DOF_SPRING:
				rbc->physics_constraint = RB_constraint_new_6dof_spring(loc, rot, rb1, rb2);

				if ((rbc->plastic_angle < 0.0f) && (rbc->plastic_dist < 0.0f))
				{
					/* no plastic mode */
					rigidbody_set_springs_active(rbc, true);
				}
				else
				{
					/*plastic mode, activate depending on flag */
					/* mark immediate activation, so we dont activate again */

					if (rbc->flag & RBC_FLAG_USE_PLASTIC)
					{
						rbc->flag |= RBC_FLAG_PLASTIC_ACTIVE;
						rigidbody_set_springs_active(rbc, true);
					}
					else
					{
						rbc->flag &= ~RBC_FLAG_PLASTIC_ACTIVE;
						rigidbody_set_springs_active(rbc, false);
					}
				}

				RB_constraint_set_equilibrium_6dof_spring(rbc->physics_constraint);

			/* fall through */
			case RBC_TYPE_6DOF:
				if (rbc->type == RBC_TYPE_6DOF)     /* a litte awkward but avoids duplicate code for limits */
					rbc->physics_constraint = RB_constraint_new_6dof(loc, rot, rb1, rb2);

				if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_X)
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->limit_lin_x_lower, rbc->limit_lin_x_upper);
				else
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_LIN_X, 0.0f, -1.0f);

				if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_Y)
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->limit_lin_y_lower, rbc->limit_lin_y_upper);
				else
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_LIN_Y, 0.0f, -1.0f);

				if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_Z)
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->limit_lin_z_lower, rbc->limit_lin_z_upper);
				else
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_LIN_Z, 0.0f, -1.0f);

				if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_X)
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->limit_ang_x_lower, rbc->limit_ang_x_upper);
				else
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_ANG_X, 0.0f, -1.0f);

				if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_Y)
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->limit_ang_y_lower, rbc->limit_ang_y_upper);
				else
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_ANG_Y, 0.0f, -1.0f);

				if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_Z)
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->limit_ang_z_lower, rbc->limit_ang_z_upper);
				else
					RB_constraint_set_limits_6dof(rbc->physics_constraint, RB_LIMIT_ANG_Z, 0.0f, -1.0f);
				break;
			case RBC_TYPE_MOTOR:
				rbc->physics_constraint = RB_constraint_new_motor(loc, rot, rb1, rb2);

				RB_constraint_set_enable_motor(rbc->physics_constraint, rbc->flag & RBC_FLAG_USE_MOTOR_LIN, rbc->flag & RBC_FLAG_USE_MOTOR_ANG);
				RB_constraint_set_max_impulse_motor(rbc->physics_constraint, rbc->motor_lin_max_impulse, rbc->motor_ang_max_impulse);
				RB_constraint_set_target_velocity_motor(rbc->physics_constraint, rbc->motor_lin_target_velocity, rbc->motor_ang_target_velocity);
				break;
			case RBC_TYPE_COMPOUND:
				rbc->physics_constraint = RB_constraint_new_compound(rb1, rb2);
				break;
		}
	}
	else { /* can't create constraint without both rigid bodies */
		return;
	}

	RB_constraint_set_enabled(rbc->physics_constraint, rbc->flag & RBC_FLAG_ENABLED);

	if (rbc->flag & RBC_FLAG_USE_BREAKING)
		RB_constraint_set_breaking_threshold(rbc->physics_constraint, rbc->breaking_threshold);
	else
		RB_constraint_set_breaking_threshold(rbc->physics_constraint, FLT_MAX);

	if (rbc->flag & RBC_FLAG_OVERRIDE_SOLVER_ITERATIONS)
		RB_constraint_set_solver_iterations(rbc->physics_constraint, rbc->num_solver_iterations);
	else
		RB_constraint_set_solver_iterations(rbc->physics_constraint, -1);

	if (rbc->physics_constraint)
	{
		RB_constraint_set_id(rbc->physics_constraint, rbc->name);
	}
}

/* Create physics sim representation of constraint given rigid body constraint settings
 * < rebuild: even if an instance already exists, replace it
 */
void BKE_rigidbody_validate_sim_shard_constraint(RigidBodyWorld *rbw, FractureModifierData *fmd, Object* ob, RigidBodyShardCon *rbc, short rebuild)
{
	/* sanity checks:
	 *	- object should have a rigid body constraint
	 *  - rigid body constraint should have at least one constrained object
	 */
	if (rbc == NULL) {
		return;
	}

	if (ELEM(NULL, rbc->mi1, rbc->mi2)) {
		if (rbc->physics_constraint) {
			RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
			RB_constraint_delete(rbc->physics_constraint);
			rbc->physics_constraint = NULL;
		}
		return;
	}

	if (rbc->mi1 == rbc->mi2) {
		//owch, this happened... better check for sanity  here
		return;
	}

	if (ELEM(NULL, rbc->mi1->rigidbody, rbc->mi2->rigidbody)) {
		return;
	}

	if (rbc->physics_constraint) {
		if (rebuild == false)
		{
			if (!(rbc->flag & RBC_FLAG_USE_KINEMATIC_DEACTIVATION))
			{
				RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
			}
		}
	}
	if (rbc->physics_constraint == NULL || rebuild || (rbc->flag & RBC_FLAG_USE_KINEMATIC_DEACTIVATION) || (rbc->flag & RBC_FLAG_NEEDS_VALIDATE)) {

		/* remove constraint if it already exists before creating a new one */
		if (rbc->physics_constraint) {
			RB_constraint_delete(rbc->physics_constraint);
			rbc->physics_constraint = NULL;
		}

		rigidbody_create_shard_physics_constraint(fmd, ob, rbc, rbw);
	}

	if ((rbw && rbw->shared->physics_world && rbc->physics_constraint)) {
		RB_dworld_add_constraint(rbw->shared->physics_world, rbc->physics_constraint, rbc->flag & RBC_FLAG_DISABLE_COLLISIONS);
	}

	rbc->flag &= ~RBC_FLAG_USE_KINEMATIC_DEACTIVATION;
	rbc->flag &= ~RBC_FLAG_NEEDS_VALIDATE;
}

static bool colgroup_check(int group1, int group2)
{
	int i = 0;
	for (i = 0; i < 20; i++)
	{
		int v1, v2;
		v1 = (group1 & (1 << i));
		v2 = (group2 & (1 << i));

		//printf("%d, %d, %d\n", i, v1, v2);

		if ((v1 > 0) && (v1 == v2))
		{
			return true;
		}
	}

	return false;
}

static bool do_activate(Object* ob, Object *ob2, Shard *mi_compare, RigidBodyWorld *rbw, Shard *mi_trigger, bool activate)
{
	FractureModifierData *fmd;
	bool valid = true;
	bool antiValid = ob2->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER;
	bool wouldActivate = false;
	Shard *mi;

	fmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);
	valid = valid && (fmd != NULL);
	antiValid = antiValid && (fmd != NULL);

	valid = valid && (ob->rigidbody_object->flag & RBO_FLAG_IS_TRIGGERED);
	valid = valid && ((ob2->rigidbody_object->flag & RBO_FLAG_IS_TRIGGER) || ((ob2->rigidbody_object->flag & RBO_FLAG_PROPAGATE_TRIGGER) &&
			((mi_trigger) && (mi_trigger->rigidbody->flag & RBO_FLAG_PROPAGATE_TRIGGER))));

#if 0
	/*prefer dynamic trigger over trigger, and allow activation after queue is empty only (everything fractured) */
	if (mi_trigger && mi_trigger->rigidbody->flag & RBO_FLAG_DYNAMIC_TRIGGER) {
		//valid = valid && BLI_listbase_is_empty(&fmd->shared->fracture_ids);
		//allow "would activate" so a collision is registered, but dont actually activate here
		return false;
	}
#endif

	if (valid || antiValid)
	{
		for (mi = fmd->shared->shards.first; mi; mi = mi->next)
		{
			bool dissolve = ob->rigidbody_object->flag & RBO_FLAG_CONSTRAINT_DISSOLVE;

			bool same_cluster = ((mi->cluster_index != -1) &&
								(mi->cluster_index == mi_compare->cluster_index));

			bool different_cluster = !same_cluster && dissolve;
			bool activate_broken = fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN;

			RigidBodyOb* rbo = mi->rigidbody;
			if ((((rbo->flag & RBO_FLAG_KINEMATIC) || different_cluster) &&
				 ((mi_compare == mi) || (same_cluster && !dissolve && activate_broken))) && valid)
			{
				if (rbo->shared->physics_object && activate) {
					BKE_rigidbody_activate(rbo, rbw, mi, ob);
				}

				wouldActivate = true;
			}

			if ((mi_compare == mi) && antiValid && activate)
			{
				if (rbo->shared->physics_object) {
					test_deactivate_rigidbody(rbo, mi);
				}
			}
		}
	}
	else if (!fmd)
	{
		RigidBodyOb* rbo = ob->rigidbody_object;
		valid = ob2->rigidbody_object->flag & RBO_FLAG_IS_TRIGGER;
		antiValid = ob2->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER;

		if (rbo && valid && activate)
		{
			if (activate)
				BKE_rigidbody_activate(rbo, rbw, NULL, ob);

			wouldActivate = true;
		}

		if (rbo && antiValid && activate)
		{
			test_deactivate_rigidbody(rbo, NULL);
		}
	}

	return wouldActivate;
}

static int check_colgroup_ghost(Object* ob1, Object *ob2)
{
	int ret = 0;
	ret = colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups);
	return ret && (!(ob1->rigidbody_object->flag & RBO_FLAG_IS_GHOST) && !(ob2->rigidbody_object->flag & RBO_FLAG_IS_GHOST));
}

#if 0
static void fake_dynamic_collide(Object *ob1, Object *ob2, MeshIsland *mi1, MeshIsland *mi2, Scene* scene)
{
	if (((ob1->rigidbody_object->flag & RBO_FLAG_IS_GHOST) && (ob1->rigidbody_object->flag & RBO_FLAG_DYNAMIC_TRIGGER)) ||
		((ob1->rigidbody_object->flag & RBO_FLAG_IS_GHOST) && (ob2->rigidbody_object->flag & RBO_FLAG_DYNAMIC_TRIGGER)) ||
		((ob2->rigidbody_object->flag & RBO_FLAG_IS_GHOST) && (ob1->rigidbody_object->flag & RBO_FLAG_DYNAMIC_TRIGGER)) ||
		((ob2->rigidbody_object->flag & RBO_FLAG_IS_GHOST) && (ob2->rigidbody_object->flag & RBO_FLAG_DYNAMIC_TRIGGER)))
	{
		rbContactPoint point;
		point.contact_force = 0;

		if (mi1 || mi2) {
			if (mi1) {
				copy_v3_v3(point.contact_pos_world_onA, mi1->centroid);
				point.contact_islandA = mi1;
				point.contact_objectA = ob1;
			}
			else {
				copy_v3_v3(point.contact_pos_world_onA, ob1->loc);
				point.contact_islandA = NULL;
				point.contact_objectA = ob1;
			}

			if (mi2) {
				copy_v3_v3(point.contact_pos_world_onB, mi2->centroid);
				point.contact_islandB = mi2;
				point.contact_objectB = ob2;
			}
			else {
				copy_v3_v3(point.contact_pos_world_onB, ob2->loc);
				point.contact_islandB = NULL;
				point.contact_objectB = ob2;
			}

			check_fracture(&point, scene);
		}
	}
}
#endif

static bool check_constraint_island(FractureModifierData* fmd, Shard *mi1, Shard *mi2)
{
	if (mi1 && mi2 && fmd && (!(fmd->flag & MOD_FRACTURE_USE_CONSTRAINT_COLLISION) ||
	                          fmd->flag & MOD_FRACTURE_USE_SELF_COLLISION))
	{
		float dist_sq = len_squared_v3v3(mi1->loc, mi2->loc);
		bool self_collision = fmd->flag & MOD_FRACTURE_USE_SELF_COLLISION;
		bool constraint_collision = fmd->flag & MOD_FRACTURE_USE_CONSTRAINT_COLLISION;

		bool is_near = len_squared_v3v3(mi1->rigidbody->pos, mi2->rigidbody->pos) < dist_sq;
		bool same_island = mi1->constraint_index == mi2->constraint_index;
		bool same_near_self = same_island && self_collision && is_near;
		bool diff_clust_near_self = mi1->cluster_index != mi2->cluster_index && self_collision && is_near;
		bool regular_case = mi1 && mi2 && !constraint_collision && !same_island;

		if (self_collision)
		{
			if (mi1->rigidbody->shared->physics_shape)
			{
				RB_shape_set_margin(mi1->rigidbody->shared->physics_shape, is_near ? 0.0f : RBO_GET_MARGIN(mi1->rigidbody));
			}

			if (mi2->rigidbody->shared->physics_shape)
			{
				RB_shape_set_margin(mi2->rigidbody->shared->physics_shape, is_near ? 0.0f : RBO_GET_MARGIN(mi2->rigidbody));
			}
		}

		//collide if: same island and near, different cluster if clustered and same island and near

		return regular_case || same_near_self || diff_clust_near_self;

	}

	if (fmd)
	{
		return !(fmd->pack_group && (fmd->flag & MOD_FRACTURE_USE_CONSTRAINT_COLLISION));
	}

	return true;
}

/* this allows partial object activation, only some shards will be activated, called from bullet(!) */
int BKE_rigidbody_filter_callback(void* scene, void* island1, void* island2, void *blenderOb1, void* blenderOb2, bool activate)
{
	//pass scene here, in order to hopefully get the original one from DEG
	Shard* mi1, *mi2;
	Scene *sc = (Scene*)scene, *sc_orig;

	RigidBodyWorld *rbw;
	Object* ob1, *ob2;
	bool validOb = true, check_activate = false;

	// oh man... the pleasures of CoW.., mooo
	sc_orig = (Scene*)DEG_get_original_id(&sc->id);
	rbw = sc_orig->rigidbody_world;

	mi1 = (Shard*)island1;
	mi2 = (Shard*)island2;

	ob1 = (Object*)blenderOb1;
	ob2 = (Object*)blenderOb2;

	FractureModifierData *fmd1 = (FractureModifierData*)modifiers_findByType(ob1, eModifierType_Fracture);
	FractureModifierData *fmd2 = (FractureModifierData*)modifiers_findByType(ob2, eModifierType_Fracture);

	if (rbw == NULL)
	{
		/* just check for ghost flags here, do not activate anything */
		return check_colgroup_ghost(ob1, ob2);
	}

	if ((mi1 != NULL) && (mi2 != NULL)) {
		validOb = (ob1 != ob2 && colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups) &&
				  ((mi1->rigidbody->flag & RBO_FLAG_KINEMATIC) || (mi2->rigidbody->flag & RBO_FLAG_KINEMATIC)) &&
				  ((mi1->rigidbody->type == RBO_TYPE_ACTIVE) && (mi2->rigidbody->type == RBO_TYPE_ACTIVE)));
	}
	else if ((mi1 == NULL) && (mi2 != NULL)) {
		validOb = (colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups) &&
				  ((ob1->rigidbody_object->flag & RBO_FLAG_KINEMATIC) || (mi2->rigidbody->flag & RBO_FLAG_KINEMATIC)) &&
				  ((ob1->rigidbody_object->type == RBO_TYPE_ACTIVE) && (mi2->rigidbody->type == RBO_TYPE_ACTIVE)));
	}
	else if ((mi1 != NULL) && (mi2 == NULL)) {
		validOb = (colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups) &&
				  ((mi1->rigidbody->flag & RBO_FLAG_KINEMATIC) || (ob2->rigidbody_object->flag & RBO_FLAG_KINEMATIC)) &&
				  ((mi1->rigidbody->type == RBO_TYPE_ACTIVE) && (ob2->rigidbody_object->type == RBO_TYPE_ACTIVE)));
	}
	else
	{
		validOb = (colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups) &&
				  ((ob1->rigidbody_object->flag & RBO_FLAG_KINEMATIC) || (ob2->rigidbody_object->flag & RBO_FLAG_KINEMATIC)) &&
				  ((ob1->rigidbody_object->type == RBO_TYPE_ACTIVE) && (ob2->rigidbody_object->type == RBO_TYPE_ACTIVE)));
	}

	if (validOb || ((colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups) &&
				   ((ob1->rigidbody_object->flag & RBO_FLAG_CONSTRAINT_DISSOLVE) ||
					(ob2->rigidbody_object->flag & RBO_FLAG_CONSTRAINT_DISSOLVE) ||
					(ob1->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER) ||
					(ob2->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER)))))
	{
		//override for 2 regular rigidbodies to enable ghost trigger functionality; else bullet wont call this again here with "activate == true"


		if (ob1->rigidbody_object->flag & RBO_FLAG_IS_TRIGGERED)
		{
			bool override = activate || (ob2->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER);
			check_activate = do_activate(ob1, ob2, mi1, rbw, mi2, override);
		}

		if (ob2->rigidbody_object->flag & RBO_FLAG_IS_TRIGGERED)
		{
			bool override = activate || (ob1->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER);
			check_activate = do_activate(ob2, ob1, mi2, rbw, mi1, override);
		}
	}

	//if ghost is involved, and dynafrac trigger is enabled, try to call check_fracture manually here, without forces and with centroid as contact point
//	fake_dynamic_collide(ob1, ob2, mi1, mi2, rbw);
//	fake_dynamic_collide(ob2, ob1, mi2, mi1, rbw);

	validOb = (check_colgroup_ghost(ob1, ob2) && ((check_constraint_island(fmd1, mi1, mi2) &&
				  check_constraint_island(fmd2, mi2, mi1)) || (ob1 != ob2)));

	//return activate ? validOb : check_activate || validOb;
	return validOb;
}

static bool can_break(Object* collider, Object* ob, bool limit)
{
	//allow limit impact only on initial shard and 1st level shards ?
	if (collider && collider->rigidbody_object && (collider->rigidbody_object->flag & RBO_FLAG_DYNAMIC_TRIGGER))
	{
		if (limit && (collider == ob)) {
			return false;
		}

		//dont allow limit impact with ground
		if (collider->rigidbody_object->type == RBO_TYPE_PASSIVE) {
			return false;
		}

		return true;
	}

	return !limit;
}

bool BKE_rigidbody_activate_by_size_check(Object *ob, Shard *mi)
{
	FractureModifierData *fmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);

	if (!fmd) {
		return true;
	}

	if (ob->rigidbody_object->flag & (RBO_FLAG_IS_TRIGGERED | RBO_FLAG_KINEMATIC) &&
	        (fmd->flag & MOD_FRACTURE_USE_DYNAMIC))
	{
		//try to keep bigger shards in place
		return BKE_rigidbody_check_island_size(fmd, mi, fmd->dynamic_activation_size);
	}

	return true;
}

bool BKE_rigidbody_check_island_size(FractureModifierData *fmd, Shard *mi, float check_size)
{
	FractureQueueEntry *fid;
	float size = fmd->dynamic_min_size, diff[3], min[3], max[3];
	INIT_MINMAX(min, max);
	BKE_mesh_minmax(mi->mesh, min, max);

	sub_v3_v3v3(diff, max, min);

	if (check_size > -1.0f) {
		size = check_size;

		if ((diff[max_axis_v3(diff)] < size))
		{
			return true;
		}

		return false;
	}
	else {
		if (diff[max_axis_v3(diff)] < size)
		{
			return false;
		}

		for (fid = fmd->shared->dynamic_fracture_queue.first; fid; fid = fid->next)
		{
			if (fid->mi->id == mi->id)
			{
				return false;
			}
		}

		return true;
	}
}

static bool check_constraints(FractureModifierData *fmd, Shard *mi, RigidBodyWorld *rbw) {
	//count broken constraints
	RigidBodyShardCon *con;
	CollectionObject* go;
	int i = 0, broken = 0;
	float percentage;

	if (rbw && rbw->constraints) {
		//delete outer constraints here, to be sure
		for (go = rbw->constraints->gobject.first; go; go = go->next) {
			RigidBodyCon *rbc = go->ob->rigidbody_constraint;
			FractureModifierData *fmd1 = (FractureModifierData*)modifiers_findByType(rbc->ob1, eModifierType_Fracture);
			FractureModifierData *fmd2 = (FractureModifierData*)modifiers_findByType(rbc->ob2, eModifierType_Fracture);

			if ((rbc && rbc->physics_constraint && rbw && rbw->shared->physics_world) && ((fmd1 == fmd) || (fmd2 == fmd))) {
				RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
				RB_constraint_delete(rbc->physics_constraint);
				rbc->physics_constraint = NULL;
			}
		}
	}

	if (mi->participating_constraint_count == 0) {
		return true;
	}

	for (i = 0; i < mi->participating_constraint_count; i++) {
		con = mi->participating_constraints[i];
		if (con->physics_constraint && !RB_constraint_is_enabled(con->physics_constraint)) {
			broken++;
		}
	}

	percentage = (float)broken / (float)mi->participating_constraint_count;

	if ((percentage * 100) >= fmd->dynamic_percentage) {
		for (i = 0; i < mi->participating_constraint_count; i++) {
			con = mi->participating_constraints[i];
			if (con->physics_constraint) {
				RB_constraint_set_enabled(con->physics_constraint, false);
			}
		}

		return true;
	}

	return false;
}

static void check_fracture_meshisland(FractureModifierData *fmd, Shard *mi, Object* ob1, Object* ob2,
                                      RigidBodyWorld *rbw, float contact_pos[3], float force)
{
	bool canbreak = false;
	bool limit_impact = fmd->flag & MOD_FRACTURE_USE_LIMIT_IMPACT;

	if (mi->rigidbody->mass > 0) {
		force = force / mi->rigidbody->mass;
	}

	canbreak = (force > fmd->dynamic_force) || ((limit_impact && ob2) && can_break(ob2, ob1, limit_impact));

	if (canbreak && check_constraints(fmd, mi, rbw))
	{
		float size[3] =  {1.0f, 1.0f, 1.0f};

		if (ob1 == ob2 || (ob2 && ob2->rigidbody_object && ob1->rigidbody_object->type == RBO_TYPE_PASSIVE)) {
			//todo calculate shard...
			size[0] = size[1] = size[2] = -1.0f;
		}
		else if (ob2) {
			BKE_object_dimensions_get(ob2, size);
		}

		copy_v3_v3(mi->impact_loc, contact_pos);
		copy_v3_v3(mi->impact_size, size);

		/*only fracture on new entries, this is necessary because after loading a file
		 *the pointcache thinks it is empty and a fracture is attempted ! */
		if (BKE_rigidbody_check_island_size(fmd, mi, -2.0f) )
		{
			if (mi->fractured == false)
			{
				FractureQueueEntry* fid = MEM_mallocN(sizeof(FractureQueueEntry), "callback_fractureid");
				fid->mi = mi;
				BLI_addtail(&fmd->shared->dynamic_fracture_queue, fid);
				fmd->shared->flag |= MOD_FRACTURE_REFRESH_DYNAMIC;
				printf("FRACTURE : %d\n", mi->id);
			}
		}
	}
}


static void check_fracture(rbContactPoint* cp, Scene *scene)
{
	Object* ob1 = NULL, *ob2 = NULL;
	Shard *mi1 = NULL, *mi2 = NULL;
	FractureModifierData *fmd1, *fmd2;
	RigidBodyWorld *rbw = NULL;
	float force;
	int frame;

	if (!cp || !scene) {
		return;
	}

	frame = (int)BKE_scene_frame_get(scene);
	force = cp->contact_force;
	ob1 = (Object*)cp->contact_objectA;
	ob2 = (Object*)cp->contact_objectB;

	if (!ob1 || !ob2) {
		return;
	}

	fmd1 = (FractureModifierData*)modifiers_findByType(ob1, eModifierType_Fracture);
	if (fmd1 && (fmd1->flag & MOD_FRACTURE_USE_DYNAMIC))
	{
		mi1 = (Shard*)cp->contact_islandA;
		check_fracture_meshisland(fmd1, mi1, ob1, ob2, rbw, cp->contact_pos_world_onA, force);
	}

	fmd2 = (FractureModifierData*)modifiers_findByType(ob2, eModifierType_Fracture);
	if (fmd2 && (fmd2->flag & MOD_FRACTURE_USE_DYNAMIC))
	{
		mi2 = (Shard*)cp->contact_islandB;
		check_fracture_meshisland(fmd2, mi2, ob2, ob1, rbw, cp->contact_pos_world_onB, force);
	}

	//free contact point ?
	cp = NULL;
}

static ThreadMutex dynamic_lock = BLI_MUTEX_INITIALIZER;
void BKE_rigidbody_contact_callback(rbContactPoint* cp, void* sc)
{
	Scene* scene = (Scene*)DEG_get_original_id(sc);

	BLI_mutex_lock(&dynamic_lock);
	check_fracture(cp,scene);
	BLI_mutex_unlock(&dynamic_lock);
}

void BKE_rigidbody_id_callback(void* island, int* objectId, int* islandId)
{
	Shard *mi = (Shard*)island;
	//RigidBodyWorld *rbw = (RigidBodyWorld*)world;

	*objectId = -1;
	*islandId = -1;

	if (mi)
	{
		*objectId = mi->object_index;
		*islandId = mi->id;
	}
}

#if 0
static void rigidbody_passive_fake_hook(MeshIsland *mi, float co[3])
{
	//no reshape necessary as vertcount didnt change, but update rbo->pos / orn ? according to change of 1st vertex
	//fake hook system
	if (mi->rigidbody->type == RBO_TYPE_PASSIVE &&
		mi->rigidbody->shared->physics_object && !(mi->rigidbody->flag & RBO_FLAG_KINEMATIC))
	{
		float oldloc[3], loc[3], diff[3], pos[3];
		copy_v3_v3(oldloc, mi->mesh->mvert->co); //1st vertex

		//this location comes from the final DM, which might be changed by hook modifiers for example
		//XXX TODO maybe need a proper switch for this behavior, too
		copy_v3_v3(loc, co);
		sub_v3_v3v3(diff, oldloc, loc);
		//sub_v3_v3(diff, mi->centroid);

		//RB_body_get_position(mi->rigidbody->physics_object, pos);
		copy_v3_v3(pos, mi->rigidbody->pos);
		//print_v3("Pos:", pos);
		//print_v3("Diff", diff);

		sub_v3_v3(pos, diff);
		RB_body_set_kinematic_state(mi->rigidbody->shared->physics_object, true);

		//XXX TODO how to handle rotation properly ? and omit if kinematic, else it will interfere
		//copy_v3_v3(mi->rigidbody->pos, pos);
		RB_body_set_loc_rot(mi->rigidbody->shared->physics_object, pos, mi->rigidbody->orn);
		//BKE_rigidbody_update_cell(mi, ob, pos, mi->rigidbody->orn, fmd, -1);
	}
}
#endif

void BKE_rigidbody_shard_validate(RigidBodyWorld *rbw, Shard *mi, Object *ob, FractureModifierData* fmd,
                                  int rebuild, int transfer_speed, float size[3], float ctime)
{
	if (mi == NULL || mi->rigidbody == NULL) {
		return;
	}

	if (BKE_fracture_meshisland_check_frame(fmd, mi, (int)ctime)) {
		RigidBodyOb *rbo = mi->rigidbody;

		if (rbw && rbo && rbo->shared->physics_object && (mi->startframe < (int)ctime))
		{
			int i = 0;
			for (i = 0; i < mi->participating_constraint_count; i++)
			{
				/* dont forget removing constraints if any */
				RigidBodyShardCon *con = mi->participating_constraints[i];
				BKE_rigidbody_remove_shard_con(rbw, con);

				//do not re-validate
				con->flag &= ~(RBC_FLAG_NEEDS_VALIDATE);
			}

			if (rbw->shared->physics_world && rbo->shared->physics_object)
				RB_dworld_remove_body(rbw->shared->physics_world, rbo->shared->physics_object);

			if (rbo->shared->physics_object) {
				RB_body_delete(rbo->shared->physics_object);
				rbo->shared->physics_object = NULL;
			}

			if (rbo->shared->physics_shape) {
				RB_shape_delete(rbo->shared->physics_shape);
				rbo->shared->physics_shape = NULL;
			}

			//do not re-validate
			rbo->flag &= ~(RBO_FLAG_NEEDS_VALIDATE | RBO_FLAG_NEEDS_RESHAPE);
		}
	}

	if (rebuild /*|| (mi->rigidbody->flag & RBO_FLAG_KINEMATIC_REBUILD)*/) {
		/* World has been rebuilt so rebuild object */
		BKE_rigidbody_validate_sim_shard(rbw, mi, ob, fmd, true, transfer_speed, size);
	}
	else if (mi->rigidbody->flag & RBO_FLAG_NEEDS_VALIDATE) {
		BKE_rigidbody_validate_sim_shard(rbw, mi, ob, fmd, false, transfer_speed, size);
	}

	/* refresh shape... */
	if (mi->rigidbody->shared->physics_object && (mi->rigidbody->flag & RBO_FLAG_NEEDS_RESHAPE)) {
		/* mesh/shape data changed, so force shape refresh */
		BKE_rigidbody_validate_sim_shard_shape(mi, ob, true);
		/* now tell RB sim about it */
		// XXX: we assume that this can only get applied for active/passive shapes that will be included as rigidbodies
		RB_body_set_collision_shape(mi->rigidbody->shared->physics_object, mi->rigidbody->shared->physics_shape);
	}
	mi->rigidbody->flag &= ~(RBO_FLAG_NEEDS_VALIDATE | RBO_FLAG_NEEDS_RESHAPE);
}

static void activateCluster(Shard *mi, int cluster_index, RigidBodyWorld *rbw, Object *ob) {
	RigidBodyShardCon *con;
	int i = 0;
	for (i = 0; i < mi->participating_constraint_count; i++)
	{
		con = mi->participating_constraints[i];
		if (con->physics_constraint && con->mi1->cluster_index == cluster_index)
		{
			if (con->mi1->rigidbody->flag & RBO_FLAG_KINEMATIC) {
				BKE_rigidbody_activate(con->mi1->rigidbody, rbw, con->mi1, ob);
			}
		}

		if (con->physics_constraint && con->mi2->cluster_index == cluster_index)
		{
			if (con->mi2->rigidbody->flag & RBO_FLAG_KINEMATIC) {
				BKE_rigidbody_activate(con->mi2->rigidbody, rbw, con->mi2, ob);
			}
		}
	}
}

static void set_constraint_index(FractureModifierData *fmd, RigidBodyShardCon *con)
{
	if (con->physics_constraint)
	{
		if (!RB_constraint_is_enabled(con->physics_constraint))
		{
			fmd->constraint_island_count++;
			con->mi1->constraint_index = fmd->constraint_island_count;

			fmd->constraint_island_count++;
			con->mi2->constraint_index = fmd->constraint_island_count;
		}
		else
		{
			con->mi1->constraint_index = fmd->constraint_island_count;
			con->mi2->constraint_index = fmd->constraint_island_count;
		}
	}
}

static void handle_breaking_percentage(FractureModifierData* fmd, Object *ob, Shard *mi, RigidBodyWorld *rbw, int breaking_percentage)
{
	int broken_cons = 0, cons = 0, i = 0, cluster_cons = 0, broken_cluster_cons = 0;
	RigidBodyShardCon *con;

	cons = mi->participating_constraint_count;
	/* calc ratio of broken cons here, per MeshIsland and flag the rest to be broken too*/
	for (i = 0; i < cons; i++) {
		con = mi->participating_constraints[i];
		if (con) {
			if (fmd->cluster_breaking_percentage > 0)
			{
				/*only count as broken if between clusters!*/
				if (con->mi1->cluster_index != con->mi2->cluster_index)
				{
					cluster_cons++;

					if (con->physics_constraint)
					{
						if (!RB_constraint_is_enabled(con->physics_constraint)) {
							broken_cluster_cons++;
						}
					}
				}
			}

			if (con->physics_constraint && !RB_constraint_is_enabled(con->physics_constraint)) {
				broken_cons++;
			}
		}
	}

	if (cluster_cons > 0) {
		if ((float)broken_cluster_cons / (float)cluster_cons * 100 >= fmd->cluster_breaking_percentage) {
			for (i = 0; i < cons; i++) {
				con = mi->participating_constraints[i];
				if (con && con->mi1->cluster_index != con->mi2->cluster_index) {
					if (fmd->flag & MOD_FRACTURE_USE_BREAKING)
					{
						if (con->physics_constraint) {

							RB_constraint_set_enabled(con->physics_constraint, false);
							/*if (con->mi1->rigidbody->flag & RBO_FLAG_KINEMATIC ||
								con->mi2->rigidbody->flag & RBO_FLAG_KINEMATIC ) */
							if (fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN)
							{
								BKE_rigidbody_activate(con->mi1->rigidbody, rbw, con->mi1, ob);
								activateCluster(con->mi1, con->mi1->cluster_index, rbw, ob);

								BKE_rigidbody_activate(con->mi2->rigidbody, rbw, con->mi2, ob);
								activateCluster(con->mi2, con->mi2->cluster_index, rbw, ob);
							}
						}
					}
				}
			}
		}
	}

	if (cons > 0 && breaking_percentage > 0) {
		if ((float)broken_cons / (float)cons * 100 >= breaking_percentage) {
			/* break all cons if over percentage */
			for (i = 0; i < cons; i++) {
				con = mi->participating_constraints[i];
				if (con && (fmd->flag & MOD_FRACTURE_USE_BREAKING))
				{
					if (con->physics_constraint) {
						RB_constraint_set_enabled(con->physics_constraint, false);
						if (fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN)
						{
							BKE_rigidbody_activate(con->mi1->rigidbody, rbw, con->mi1, ob);
							BKE_rigidbody_activate(con->mi2->rigidbody, rbw, con->mi2, ob);
						}
					}
				}
			}
		}
	}
}

static void test_deactivate_rigidbody(RigidBodyOb *rbo, Shard* mi)
{
	if (rbo->shared->physics_object && ((rbo->flag & RBO_FLAG_KINEMATIC) == 0)) {
		float lin_vel[3], ang_vel[3];

		RB_body_get_linear_velocity(rbo->shared->physics_object, lin_vel);
		RB_body_get_angular_velocity(rbo->shared->physics_object, ang_vel);

		if (((len_squared_v3(lin_vel) < (rbo->lin_sleep_thresh * rbo->lin_sleep_thresh))) ||
		   ((len_squared_v3(ang_vel) < (rbo->ang_sleep_thresh * rbo->ang_sleep_thresh))))
		{
			int i = 0;

			rbo->flag |= RBO_FLAG_KINEMATIC;
			rbo->flag |= RBO_FLAG_KINEMATIC_REBUILD;
			rbo->flag |= RBO_FLAG_NEEDS_VALIDATE;

			if ((mi != NULL) && (rbo->flag & RBO_FLAG_PROPAGATE_TRIGGER))
			{
				for (i = 0; i < mi->participating_constraint_count; i++)
				{
					RigidBodyShardCon *con;
					con = mi->participating_constraints[i];

					if (con && con->physics_constraint && RB_constraint_is_enabled(con->physics_constraint))
					{
						RigidBodyOb *rb1 = con->mi1->rigidbody;
						RigidBodyOb *rb2 = con->mi2->rigidbody;

						RB_constraint_set_enabled(con->physics_constraint, false);

						rb1->flag |= RBO_FLAG_KINEMATIC;
						rb1->flag |= RBO_FLAG_KINEMATIC_REBUILD;
						rb1->flag |= RBO_FLAG_NEEDS_VALIDATE;

						rb2->flag |= RBO_FLAG_KINEMATIC;
						rb2->flag |= RBO_FLAG_KINEMATIC_REBUILD;
						rb2->flag |= RBO_FLAG_NEEDS_VALIDATE;
					}
				}
			}
		}
	}
}

static void deform_constraint(FractureModifierData *fmd, Object *ob, RigidBodyShardCon* rbsc, RigidBodyWorld *rbw)
{
	float thresh = RB_constraint_get_breaking_threshold(rbsc->physics_constraint);
	float weakening = 1.0f - fmd->deform_weakening;

	RB_dworld_remove_constraint(rbw->shared->physics_world, rbsc->physics_constraint);

	BKE_rigidbody_start_dist_angle(rbsc, true, false);
	BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, true);

	RB_constraint_set_breaking_threshold(rbsc->physics_constraint, thresh * weakening);

	RB_body_deactivate(rbsc->mi1->rigidbody->shared->physics_object);
	RB_body_deactivate(rbsc->mi2->rigidbody->shared->physics_object);
}

static void handle_deform_angle(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
								float anglediff, float weight, float deform_angle)
{
	if ((fmd->deform_angle > 0 || ((fmd->flag & MOD_FRACTURE_USE_DEFORM_ANGLE_WEIGHTED) && weight > 0)) &&
		(anglediff > deform_angle))
	{
		/* if we have cluster breaking angle, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_deform_angle > 0 && rbsc->mi1->cluster_index == rbsc->mi2->cluster_index) ||
			 fmd->cluster_deform_angle == 0)
		{
			deform_constraint(fmd, ob, rbsc, rbw);
		}
	}

	if ((fmd->cluster_deform_angle > 0) && (rbsc->mi1->cluster_index != rbsc->mi2->cluster_index)
		&& anglediff > fmd->cluster_deform_angle)
	{
		deform_constraint(fmd, ob, rbsc, rbw);
	}
}

static void handle_deform_dist(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
								float distdiff, float weight, float deform_dist)
{
	if ((fmd->deform_distance > 0 || ((fmd->flag & MOD_FRACTURE_USE_DEFORM_DISTANCE_WEIGHTED) && weight > 0)) &&
		(distdiff > deform_dist))
	{
		/* if we have cluster breaking angle, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_deform_distance > 0 && rbsc->mi1->cluster_index == rbsc->mi2->cluster_index) ||
			 fmd->cluster_deform_distance == 0)
		{
			deform_constraint(fmd, ob, rbsc, rbw);
		}
	}

	if ((fmd->cluster_deform_distance > 0) && (rbsc->mi1->cluster_index != rbsc->mi2->cluster_index)
		&& distdiff > fmd->cluster_deform_distance)
	{
		deform_constraint(fmd, ob, rbsc, rbw);
	}
}


static void handle_breaking_angle(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
								  float anglediff, float weight, float breaking_angle)
{
	if ((fmd->breaking_angle > 0 || ((fmd->flag & MOD_FRACTURE_USE_BREAKING_ANGLE_WEIGHTED) && weight > 0)) &&
		(anglediff > breaking_angle))
	{
		/* if we have cluster breaking angle, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_breaking_angle > 0 && rbsc->mi1->cluster_index == rbsc->mi2->cluster_index &&
			 rbsc->mi1->cluster_index != -1) || fmd->cluster_breaking_angle == 0)
		{
			if (fmd->flag & MOD_FRACTURE_USE_BREAKING)
			{
				//break constraint
				if (rbsc->physics_constraint) {
					RB_constraint_set_enabled(rbsc->physics_constraint, false);
					if (fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN)
					{
						BKE_rigidbody_activate(rbsc->mi1->rigidbody, rbw, rbsc->mi1, ob);
						BKE_rigidbody_activate(rbsc->mi2->rigidbody, rbw, rbsc->mi2, ob);
					}
				}
			}
		}
	}

	if ((fmd->cluster_breaking_angle > 0) && (rbsc->mi1->cluster_index != rbsc->mi2->cluster_index)
		&& anglediff > fmd->cluster_breaking_angle)
	{
		if (fmd->flag & MOD_FRACTURE_USE_BREAKING)
		{
			if (rbsc->physics_constraint) {
				RB_constraint_set_enabled(rbsc->physics_constraint, false);
				if (fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN)
				{
					BKE_rigidbody_activate(rbsc->mi1->rigidbody, rbw, rbsc->mi1, ob);
					BKE_rigidbody_activate(rbsc->mi2->rigidbody, rbw, rbsc->mi2, ob);
				}
			}
		}
	}
}

static void handle_breaking_distance(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
									 float distdiff, float weight, float breaking_distance)
{
	if ((fmd->breaking_distance > 0 || ((fmd->flag &MOD_FRACTURE_USE_BREAKING_DISTANCE_WEIGHTED) && weight > 0)) &&
		(distdiff > breaking_distance))
	{
		/* if we have cluster breaking distance, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_breaking_distance > 0 && rbsc->mi1->cluster_index == rbsc->mi2->cluster_index &&
			 rbsc->mi1->cluster_index != -1) || fmd->cluster_breaking_distance == 0)
		{
			if (fmd->flag & MOD_FRACTURE_USE_BREAKING)
			{
				if (rbsc->physics_constraint) {
					RB_constraint_set_enabled(rbsc->physics_constraint, false);
					if (fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN)
					{
						BKE_rigidbody_activate(rbsc->mi1->rigidbody, rbw, rbsc->mi1, ob);
						BKE_rigidbody_activate(rbsc->mi2->rigidbody, rbw, rbsc->mi2, ob);
					}
				}
			}
		}
	}

	if ((fmd->cluster_breaking_distance > 0) && (rbsc->mi1->cluster_index != rbsc->mi2->cluster_index)
		&& distdiff > fmd->cluster_breaking_distance)
	{
		if (fmd->flag & MOD_FRACTURE_USE_BREAKING)
		{
			if (rbsc->physics_constraint) {
				RB_constraint_set_enabled(rbsc->physics_constraint, false);
				if (fmd->flag & MOD_FRACTURE_USE_ACTIVATE_BROKEN)
				{
					BKE_rigidbody_activate(rbsc->mi1->rigidbody, rbw, rbsc->mi1, ob);
					BKE_rigidbody_activate(rbsc->mi2->rigidbody, rbw, rbsc->mi2, ob);
				}
			}
		}
	}
}

static void enable_plastic(RigidBodyShardCon *rbsc)
{
	if (!(rbsc->flag & RBC_FLAG_PLASTIC_ACTIVE) && rbsc->plastic_dist >= 0.0f && rbsc->plastic_angle >= 0.0f)
	{
		if (rbsc->physics_constraint)
		{
			/* activate only once */
			rbsc->flag |= RBC_FLAG_PLASTIC_ACTIVE;
			rigidbody_set_springs_active(rbsc, true);
			RB_constraint_set_equilibrium_6dof_spring(rbsc->physics_constraint);
			RB_constraint_set_enabled(rbsc->physics_constraint, true);
		}
	}
}

static void handle_plastic_breaking(RigidBodyShardCon *rbsc, RigidBodyWorld* rbw, short laststeps, float lastscale)
{
	float dist, angle, distdiff, anglediff;
	bool exceededAngle = false, exceededDist = false;

	/*match breaking threshold according to timescale and steps */
	if (rbsc->physics_constraint)
	{
		float step_ratio = (float)rbw->steps_per_second / (float)laststeps;
		float time_ratio = lastscale / rbw->time_scale;

		/*in case of generic without linear locks, ignore time ratio*/
		if ((rbsc->type == RBC_TYPE_6DOF || rbsc->type == RBC_TYPE_6DOF_SPRING) &&
			((rbsc->flag & (RBC_FLAG_USE_LIMIT_LIN_X | RBC_FLAG_USE_LIMIT_LIN_Y | RBC_FLAG_USE_LIMIT_LIN_Z)) == 0))
		{
			time_ratio = 1.0f;
		}

		RB_constraint_set_breaking_threshold(rbsc->physics_constraint, (rbsc->breaking_threshold / step_ratio) * time_ratio);
	}

	calc_dist_angle(rbsc, &dist, &angle, true);

	/* note, relative change in percentage is expected, -1 disables */
	distdiff = fabs(1.0f -(rbsc->start_dist/dist));

	/* TODO, ensure rigidbody orn is equal to quaternion of object !!! */
	// The construct "asin(sin(x))" is a triangle function to achieve a seamless rotation loop from input
	anglediff = asin(sin(fabs(rbsc->start_angle - angle) * 0.5f));

	//printf("Dist, Angle: %f %f %f %f %f %f\n", rbsc->start_dist, rbsc->start_angle, dist, angle, distdiff, anglediff);

	exceededAngle = ((rbsc->breaking_angle > 0.0f) && (anglediff > rbsc->breaking_angle));
	exceededDist = ((rbsc->breaking_dist > 0.0f) && (distdiff > (rbsc->breaking_dist + (anglediff / M_PI))));

	if (exceededDist || exceededAngle)
	{
		if (rbsc->type == RBC_TYPE_6DOF_SPRING)
		{
			enable_plastic(rbsc);
		}
		else if (rbsc->physics_constraint)
		{
			/* break regular connections */
			RB_constraint_set_enabled(rbsc->physics_constraint, false);
		}
	}

	exceededAngle = ((rbsc->plastic_angle > 0.0f) && (anglediff > rbsc->plastic_angle));
	exceededDist = ((rbsc->plastic_dist > 0.0f) && (distdiff > (rbsc->plastic_dist + (anglediff / M_PI))));

	/* break plastic connections */
	if ((exceededDist || exceededAngle))
	{
		if (rbsc->type == RBC_TYPE_6DOF_SPRING && rbsc->flag & RBC_FLAG_PLASTIC_ACTIVE)
		{
			if (rbsc->physics_constraint)
			{
				rigidbody_set_springs_active(rbsc, false);
				RB_constraint_set_enabled(rbsc->physics_constraint, rbsc->flag & RBC_FLAG_ENABLED);
			}
		}
	}
}

static void handle_regular_breaking(FractureModifierData *fmd, Object *ob, RigidBodyWorld *rbw, RigidBodyShardCon *rbsc,
                                    float max_con_mass)
{
	bool breaking_distance_weighted = fmd->flag & MOD_FRACTURE_USE_BREAKING_DISTANCE_WEIGHTED;
	bool breaking_angle_weighted = fmd->flag & MOD_FRACTURE_USE_BREAKING_ANGLE_WEIGHTED;
	bool deform_distance_weighted = fmd->flag & MOD_FRACTURE_USE_BREAKING_ANGLE_WEIGHTED;
	bool deform_angle_weighted = fmd->flag & MOD_FRACTURE_USE_DEFORM_ANGLE_WEIGHTED;

	float weight = MIN2(rbsc->mi1->thresh_weight, rbsc->mi2->thresh_weight);
	float breaking_angle = breaking_angle_weighted ? fmd->breaking_angle * weight : fmd->breaking_angle;
	float breaking_distance = breaking_distance_weighted ? fmd->breaking_distance * weight : fmd->breaking_distance;
	float deform_angle = deform_angle_weighted ? fmd->deform_angle * weight : fmd->deform_angle;
	float deform_distance = deform_distance_weighted ? fmd->deform_distance * weight : fmd->deform_distance;
	float dist, angle, distdiff, anglediff;

	if (fmd->flag & MOD_FRACTURE_USE_MASS_DEP_THRESHOLDS) {
		BKE_rigidbody_calc_threshold(max_con_mass, fmd, rbsc);
	}

	if ((((fmd->breaking_angle) > 0) || (breaking_angle_weighted && weight > 0) ||
		(fmd->breaking_distance > 0) || (breaking_distance_weighted && weight > 0) ||
		 (fmd->cluster_breaking_angle > 0 || (fmd->cluster_breaking_distance > 0))) /*&& !rebuild*/)
	{
		calc_dist_angle(rbsc, &dist, &angle, false);
		anglediff = fabs(angle - rbsc->start_angle);
		distdiff = fabs(dist - rbsc->start_dist);

		/* handle breaking */
		handle_breaking_angle(fmd, ob, rbsc, rbw, anglediff, weight, breaking_angle);
		handle_breaking_distance(fmd, ob, rbsc, rbw, distdiff, weight, breaking_distance);
	}

	if ((((fmd->deform_angle) > 0) || (deform_angle_weighted && weight > 0) ||
		(fmd->deform_distance > 0) || (deform_distance_weighted && weight > 0) ||
		 (fmd->cluster_deform_angle > 0 || (fmd->cluster_deform_distance > 0))) /*&& !rebuild*/)
	{
		if (rbsc->physics_constraint && RB_constraint_is_enabled(rbsc->physics_constraint))
		{
			calc_dist_angle(rbsc, &dist, &angle, false);
			anglediff = fabs(angle - rbsc->start_angle_deform);
			distdiff = fabs(dist - rbsc->start_dist_deform);

			/* handle deform */
			handle_deform_angle(fmd, ob, rbsc, rbw, anglediff, weight, deform_angle);
			handle_deform_dist(fmd, ob, rbsc, rbw, distdiff, weight, deform_distance);
		}
	}
}

static void handle_solver_iterations(RigidBodyWorld *rbw, FractureModifierData *fmd, RigidBodyShardCon *rbsc)
{
	int iterations;

	if (fmd->solver_iterations_override == 0) {
		iterations = rbw->num_solver_iterations;
	}
	else {
		if (rbsc && rbsc->mi1 && rbsc->mi2)
		{
			if ((rbsc->mi1->cluster_index != -1) && (rbsc->mi1->cluster_index == rbsc->mi2->cluster_index)) {
				iterations = fmd->cluster_solver_iterations_override;
			}
			else {
				iterations = fmd->solver_iterations_override;
			}
		}
		else {
			iterations = rbw->num_solver_iterations;
		}
	}

	/* dont automatically enable in External mode */
	if (iterations > 0) {
		rbsc->flag |= RBC_FLAG_OVERRIDE_SOLVER_ITERATIONS;
		rbsc->num_solver_iterations = iterations;
	}
}

bool BKE_rigidbody_modifier_update(Scene* scene, Object* ob, RigidBodyWorld *rbw, bool rebuild, Depsgraph *depsgraph)
{
	Shard *mi;
	RigidBodyShardCon *rbsc;
	short laststeps = rbw->steps_per_second;
	float lastscale = rbw->time_scale;
	int i = 0;
	FractureModifierData *fmd = NULL;
	//Scene *sc = (Scene*)DEG_get_original_id(&scene->id);
	float frame = BKE_scene_frame_get(scene); //maybe get original scene for correct frame

	fmd = (FractureModifierData*) modifiers_findByType(ob, eModifierType_Fracture);

	if (BKE_rigidbody_modifier_active(fmd)) {
		float max_con_mass = 0;
		bool is_empty = BLI_listbase_is_empty(&fmd->shared->shards);
		int count = 0, brokencount = 0, plastic = 0;
		float size[3] = {1.0f, 1.0f, 1.0f};
		float bbsize[3];
		float locbb[3];

		BKE_object_where_is_calc(depsgraph, scene, ob);
		fmd->constraint_island_count = 1;

		for (mi = fmd->shared->shards.first; mi; mi = mi->next) {
			if (mi->rigidbody == NULL) {
				continue;
			}
			else {  /* as usual, but for each shard now, and no constraints*/
				/* perform simulation data updates as tagged */
				/* refresh object... */
				int do_rebuild = rebuild;

				BKE_mesh_boundbox_calc(mi->mesh, locbb, bbsize);

				if ((rbw->flag & RBW_FLAG_REBUILD_CONSTRAINTS) && !(fmd->flag & MOD_FRACTURE_USE_DYNAMIC))
				{
					mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
					mi->rigidbody->flag &= ~RBO_FLAG_PROPAGATE_TRIGGER;
				}


				if (fmd->flag & MOD_FRACTURE_USE_BREAKING)
				{
					bool breaking_percentage_weighted = fmd->flag & MOD_FRACTURE_USE_BREAKING_PERCENTAGE_WEIGHTED;
					float weight = mi->thresh_weight;
					int breaking_percentage = breaking_percentage_weighted ? (fmd->breaking_percentage * weight) :
																				  fmd->breaking_percentage;

					if (fmd->breaking_percentage > 0 || (breaking_percentage_weighted && weight > 0) ||
						(fmd->cluster_breaking_percentage > 0))
					{
						handle_breaking_percentage(fmd, ob, mi, rbw, breaking_percentage);
					}
				}

				BKE_rigidbody_shard_validate(rbw, is_empty ? NULL : mi, ob, fmd, do_rebuild,
											 fmd->flag & MOD_FRACTURE_USE_DYNAMIC, bbsize, frame);

				mi->constraint_index = mi->id;

			}

			BKE_rigidbody_update_sim_ob(scene, rbw, ob, mi->rigidbody, mi->loc, mi, size, fmd, depsgraph);
		}

		if (fmd->flag & MOD_FRACTURE_USE_MASS_DEP_THRESHOLDS) {
			max_con_mass = BKE_rigidbody_calc_max_con_mass(ob);
		}

		for (rbsc = fmd->shared->constraints.first; rbsc; rbsc = rbsc->next) {

			//sanity check
			if (!rbsc || !rbsc->mi1 || !rbsc->mi2)
				continue;

			if (rebuild)
			{
				rbsc->start_angle = 0.0f;
				rbsc->start_dist = 0.0f;
			}

			/* needs probably taken into account BEFORE validation already, instead of after it
			 * (causing a delay by one frame with old value and thus different simulation behavior at start) */
			handle_solver_iterations(rbw, fmd, rbsc);

			if (rbsc->physics_constraint && !(RB_constraint_is_enabled(rbsc->physics_constraint)))
			{
				brokencount++;
			}

			if (rbsc->type == RBC_TYPE_6DOF_SPRING && (rbsc->flag & RBC_FLAG_PLASTIC_ACTIVE))
			{
				plastic++;
			}

			count++;


			if (rbsc->physics_constraint && rbw && (rbw->flag & RBW_FLAG_REBUILD_CONSTRAINTS) && !rebuild) {
				//printf("Rebuilding constraints\n");
#if 0
				if (!fmd->is_dynamic_external) {
					RB_constraint_set_enabled(rbsc->physics_constraint, rbsc->flag & RBC_FLAG_ENABLED);
					rbsc->flag |= RBC_FLAG_NEEDS_VALIDATE;
				}
#endif

				if (rbsc->type == RBC_TYPE_6DOF_SPRING)
				{
					if (rbsc->plastic_angle >= 0.0f || rbsc->plastic_dist >= 0.0f)
					{
						/*reset plastic constraints with immediate activation*/
						if (rbsc->flag & RBC_FLAG_USE_PLASTIC)
						{
							if (!(rbsc->flag & RBC_FLAG_PLASTIC_ACTIVE))
							{
								rbsc->flag |= RBC_FLAG_PLASTIC_ACTIVE;
								rigidbody_set_springs_active(rbsc, true);
								RB_constraint_set_enabled(rbsc->physics_constraint, true);
								if (rbsc->physics_constraint)
									RB_constraint_set_equilibrium_6dof_spring(rbsc->physics_constraint);
							}
						}
						else
						{
							rigidbody_set_springs_active(rbsc, false);
							RB_constraint_set_enabled(rbsc->physics_constraint, false);
							rbsc->flag &= ~RBC_FLAG_PLASTIC_ACTIVE;
						}
					}
				}
			}

			if (rebuild) {
				/* World has been rebuilt so rebuild constraint */
				BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, true);
				BKE_rigidbody_start_dist_angle(rbsc, true, true);
				//TODO ensure evaluation on transform change too
			}
			else if (rbsc->flag & RBC_FLAG_NEEDS_VALIDATE && !(fmd->flag & MOD_FRACTURE_USE_DYNAMIC)) {
				BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, false);
			}
			else if (fmd->flag & MOD_FRACTURE_USE_DYNAMIC) {
				if (rbsc->mi1 && rbsc->mi2) {
					if (!BKE_fracture_meshisland_check_frame(fmd, rbsc->mi1, frame) &&
					    !BKE_fracture_meshisland_check_frame(fmd, rbsc->mi2, frame))
					{
						rbsc->flag |= RBC_FLAG_NEEDS_VALIDATE;
						BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, false);
					}
				}
			}

			if (!rebuild)
			{
				handle_regular_breaking(fmd, ob, rbw, rbsc, max_con_mass);
			}

			if (!rebuild)
			{
				handle_plastic_breaking(rbsc, rbw, laststeps, lastscale);
			}

			if (rbsc->physics_constraint)
			{
				RigidBodyOb *rbo1 = rbsc->mi1->rigidbody;
				RigidBodyOb *rbo2 = rbsc->mi2->rigidbody;

				if ((rbo1->force_thresh > 0 || rbo2->force_thresh > 0))
				{
					if (RB_constraint_get_applied_impulse(rbsc->physics_constraint) >= rbo1->force_thresh + rbo2->force_thresh)
					{
						//TODO, should be the actual objects, not just "ob"... can differ in case of external constraints...
						RB_constraint_set_enabled(rbsc->physics_constraint, false);
						BKE_rigidbody_activate(rbo1, rbw, NULL, ob);
						BKE_rigidbody_activate(rbo2, rbw, NULL, ob);
					}
				}
			}

			set_constraint_index(fmd, rbsc);

			rbsc->flag &= ~RBC_FLAG_NEEDS_VALIDATE;
			lastscale = rbw->time_scale;
			laststeps = rbw->steps_per_second;

			i++;
		}

		printf("Constraints: Frame %d , Total %d,  Intact %d,  Broken %d, Plastic %d\n", (int)frame, count,
		       count-brokencount, brokencount, plastic);

		return true;
	}
	else
	{
		return false;
	}
}

bool BKE_rigidbody_modifier_sync(ModifierData *md, Object *ob, Scene *scene, float ctime)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	bool modFound = false;
	FractureModifierData *fmd = NULL;
	Shard *mi;
	RigidBodyOb *rbo;
	float size[3] = {1, 1, 1};
	float centr[3];

	if (md->type == eModifierType_Fracture) {

		fmd = (FractureModifierData *)md;

		if (BKE_rigidbody_modifier_active(fmd)) {
			modFound = true;

			invert_m4_m4(ob->imat, ob->obmat);

			for (mi = fmd->shared->shards.first; mi; mi = mi->next) {

				rbo = mi->rigidbody;
				if (!rbo || !ob->rigidbody_object) {
					continue;
				}

				/* use rigid body transform after cache start frame if objects is not being transformed */
				if (BKE_rigidbody_check_sim_running(rbw, ctime) && !(ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ)) {

					/* keep original transform when the simulation is muted */
					if (rbw->flag & RBW_FLAG_MUTED)
						break;
				}
				/* otherwise set rigid body transform to current obmat*/
				else //if (!(ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ))
				{
					mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);
					mat4_to_size(size, ob->obmat);
					copy_v3_v3(centr, mi->loc);
					mul_v3_v3(centr, size);
					mul_qt_v3(rbo->orn, centr);
					add_v3_v3(rbo->pos, centr);
					mul_qt_qtqt(rbo->orn, rbo->orn, mi->rot);

					zero_v3(rbo->lin_vel);
					zero_v3(rbo->ang_vel);
				}
			}

			if ((ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ) ||
				((ob->rigidbody_object) && (ob->rigidbody_object->flag & RBO_FLAG_KINEMATIC)))
			{
				if (ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ && rbw) {
					RigidBodyShardCon *con;

					BKE_rigidbody_cache_reset(scene);
					/* re-enable all constraints as well */
					for (con = fmd->shared->constraints.first; con; con = con->next) {
						if (con->physics_constraint)
							RB_constraint_set_enabled(con->physics_constraint, con->flag & RBC_FLAG_ENABLED);
					}
				}
			}

			return modFound;
		}
	}

	return modFound;
}

bool BKE_restoreKinematic(RigidBodyWorld *rbw, bool override_bind)
{
	CollectionObject *go;
	bool did_it = false;

	/*restore kinematic state of shards if object is kinematic*/
	for (go = rbw->group->gobject.first; go; go = go->next)	{
		bool kinematic = false, triggered = false;

		if ((go->ob) && (go->ob->rigidbody_object)) {
			kinematic = go->ob->rigidbody_object->flag & RBO_FLAG_KINEMATIC;
			triggered = go->ob->rigidbody_object->flag & RBO_FLAG_IS_TRIGGERED;

			FractureModifierData *fmd = (FractureModifierData*)modifiers_findByType(go->ob, eModifierType_Fracture);
			if (fmd && triggered)
			{
				Shard* mi;
				for (mi = fmd->shared->shards.first; mi; mi = mi->next)
				{
					if (mi->rigidbody)
					{
						mi->rigidbody->flag &= ~RBO_FLAG_KINEMATIC_REBUILD;
						if (kinematic)
						{
							bool use_animated_mesh = fmd->flag & MOD_FRACTURE_USE_ANIMATED_MESH;
							if (!use_animated_mesh || (use_animated_mesh &&
							   (override_bind || mi->rigidbody->flag & RBO_FLAG_KINEMATIC_BOUND)))
							{
								if (use_animated_mesh && override_bind)
								{
									mi->rigidbody->flag |= RBO_FLAG_KINEMATIC_BOUND;
								}

								mi->rigidbody->flag |= RBO_FLAG_KINEMATIC;
							}
						}
						else
						{
							if (fmd->flag & MOD_FRACTURE_USE_ANIMATED_MESH) {
								mi->rigidbody->flag |= RBO_FLAG_KINEMATIC_BOUND;
							}

							//might happen if being hit by a stop trigger, remove kinematic here in this case
							mi->rigidbody->flag &= ~RBO_FLAG_KINEMATIC;
						}
						mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
						did_it = true;
					}
				}
			}
			else if (!fmd && triggered)
			{	/* restore regular triggered objects back to kinematic at all, they very likely were kinematic before...
				 * user has to disable triggered if behavior is not desired */
				go->ob->rigidbody_object->flag &= ~RBO_FLAG_KINEMATIC_REBUILD;
				go->ob->rigidbody_object->flag |= RBO_FLAG_KINEMATIC;
				go->ob->rigidbody_object->flag |= RBO_FLAG_NEEDS_VALIDATE;
				did_it = true;
			}
		}
	}

	return did_it;
}

/* Add rigid body settings to the specified shard */
RigidBodyOb *BKE_rigidbody_create_shard(Object *ob, Object *target, Shard *mi)
{
	RigidBodyOb *rbo;
	float centr[3], size[3];

	/* sanity checks
	 *	- rigidbody world must exist
	 *	- shard must exist
	 *	- cannot add rigid body if it already exists
	 */
	if (mi == NULL || (mi->rigidbody != NULL))
		return NULL;

	if (ob->type != OB_MESH) {
		return NULL;
	}

	if ((ob->type == OB_MESH) && (((Mesh *)ob->data)->totvert == 0)) {
		return NULL;
	}

	if (!ob->rigidbody_object) {
		return NULL;
	}

	/* since we are always member of an object, dupe its settings,
	 * create new settings data, and link it up */
	if (target && target->rigidbody_object)
	{
		rbo = BKE_rigidbody_copy_object(target, 0);

		mat4_to_loc_quat(rbo->pos, rbo->orn, target->obmat);
		zero_v3(rbo->lin_vel);
		zero_v3(rbo->ang_vel);
	}
	else
	{
		/* regular FM case */
		rbo = BKE_rigidbody_copy_object(ob, 0);
		//rbo->type = mi->passive_weight > 0.01f ? RBO_TYPE_PASSIVE : RBO_TYPE_ACTIVE;

		/* set initial transform */
		mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);
		mat4_to_size(size, ob->obmat);

		//add initial "offset" (centroid), maybe subtract ob->obmat ?? (not sure)
		copy_v3_v3(centr, mi->loc);
		mul_v3_v3(centr, size);
		mul_qt_v3(rbo->orn, centr);
		add_v3_v3(rbo->pos, centr);
		zero_v3(rbo->lin_vel);
		zero_v3(rbo->ang_vel);
	}

	/* return this object */
	return rbo;
}

/* Add rigid body constraint to the specified object */
RigidBodyShardCon *BKE_rigidbody_create_shard_constraint(Scene *scene, short type, bool reset)
{
	RigidBodyShardCon *rbc;

	/* sanity checks
	 *	- rigidbody world must exist
	 *	- object must exist
	 *	- cannot add constraint if it already exists
	 */

	/* create new settings data, and link it up */
	rbc = MEM_callocN(sizeof(RigidBodyShardCon), "RigidBodyShardCon");

	/* set default settings */
	rbc->type = type;

	rbc->mi1 = NULL;
	rbc->mi2 = NULL;

	rbc->flag |= RBC_FLAG_ENABLED;
	rbc->flag &= ~RBC_FLAG_DISABLE_COLLISIONS;
	rbc->flag |= RBC_FLAG_USE_BREAKING;

	rbc->breaking_threshold = 10.0f; /* no good default here, just use 10 for now */
	rbc->num_solver_iterations = 10; /* 10 is Bullet default */

	rbc->limit_lin_x_lower = -1.0f;
	rbc->limit_lin_x_upper = 1.0f;
	rbc->limit_lin_y_lower = -1.0f;
	rbc->limit_lin_y_upper = 1.0f;
	rbc->limit_lin_z_lower = -1.0f;
	rbc->limit_lin_z_upper = 1.0f;
	rbc->limit_ang_x_lower = -M_PI_4;
	rbc->limit_ang_x_upper = M_PI_4;
	rbc->limit_ang_y_lower = -M_PI_4;
	rbc->limit_ang_y_upper = M_PI_4;
	rbc->limit_ang_z_lower = -M_PI_4;
	rbc->limit_ang_z_upper = M_PI_4;

	rbc->spring_damping_x = 0.5f;
	rbc->spring_damping_y = 0.5f;
	rbc->spring_damping_z = 0.5f;
	rbc->spring_damping_ang_x = 0.5f;
	rbc->spring_damping_ang_y = 0.5f;
	rbc->spring_damping_ang_z = 0.5f;
	rbc->spring_stiffness_x = 10.0f;
	rbc->spring_stiffness_y = 10.0f;
	rbc->spring_stiffness_z = 10.0f;
	rbc->spring_stiffness_ang_x = 10.0f;
	rbc->spring_stiffness_ang_y = 10.0f;
	rbc->spring_stiffness_ang_z = 10.0f;

	rbc->motor_lin_max_impulse = 1.0f;
	rbc->motor_lin_target_velocity = 1.0f;
	rbc->motor_ang_max_impulse = 1.0f;
	rbc->motor_ang_target_velocity = 1.0f;
	strcpy(rbc->name, "");
	zero_v3(rbc->pos);
	unit_qt(rbc->orn);

	/* flag cache as outdated */
	if (reset)
		BKE_rigidbody_cache_reset(scene);

	/* return this object */
	return rbc;
}

void BKE_rigidbody_remove_shard_con(RigidBodyWorld *rbw, RigidBodyShardCon *con)
{
	if (rbw && rbw->shared->physics_world && con && con->physics_constraint) {
		RB_dworld_remove_constraint(rbw->shared->physics_world, con->physics_constraint);
		RB_constraint_delete(con->physics_constraint);
		con->physics_constraint = NULL;
	}
}

void BKE_rigidbody_remove_shard(Scene *scene, Shard *mi)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	int i = 0;

	/* rbw can be NULL directly after linking / appending objects without their original scenes
	 * if an attempt to refracture is done then, this would crash here with null pointer access */
	if (mi->rigidbody != NULL && rbw != NULL) {

		RigidBodyShardCon *con;

		for (i = 0; i < mi->participating_constraint_count; i++) {
			con = mi->participating_constraints[i];
			BKE_rigidbody_remove_shard_con(rbw, con);
		}

		if (rbw->shared->physics_world && mi->rigidbody && mi->rigidbody->shared->physics_object)
			RB_dworld_remove_body(rbw->shared->physics_world, mi->rigidbody->shared->physics_object);

		if (mi->rigidbody->shared->physics_object) {
			RB_body_delete(mi->rigidbody->shared->physics_object);
			mi->rigidbody->shared->physics_object = NULL;
		}

		if (mi->rigidbody->shared->physics_shape) {
			RB_shape_delete(mi->rigidbody->shared->physics_shape);
			mi->rigidbody->shared->physics_shape = NULL;
		}
	}
}

bool BKE_rigidbody_remove_modifier(RigidBodyWorld* rbw, ModifierData *md, Object *ob)
{
	RigidBodyShardCon *con;
	Shard *mi;
	FractureModifierData *fmd;
	bool modFound = false;

	if (md->type == eModifierType_Fracture)
	{
		fmd = (FractureModifierData *)md;
		modFound = true;
		CollectionObject *go;

		for (con = fmd->shared->constraints.first; con; con = con->next) {
			if (rbw && rbw->shared->physics_world && con->physics_constraint) {
				RB_dworld_remove_constraint(rbw->shared->physics_world, con->physics_constraint);
				RB_constraint_delete(con->physics_constraint);
				con->physics_constraint = NULL;
			}
		}

		/*if we are part of a connected object, delete the parent's constraints here too*/
		for (go = rbw->group->gobject.first; go; go = go->next)
		{
			FractureModifierData *fmdi = (FractureModifierData*)modifiers_findByType(go->ob, eModifierType_Fracture);
			if (fmdi && ob != go->ob)
			{
				for (con = fmdi->shared->constraints.first; con; con = con->next) {
					if (rbw && rbw->shared->physics_world && con->physics_constraint) {
						RB_dworld_remove_constraint(rbw->shared->physics_world, con->physics_constraint);
						RB_constraint_delete(con->physics_constraint);
						con->physics_constraint = NULL;
					}
				}
			}
		}


		for (mi = fmd->shared->shards.first; mi; mi = mi->next) {
			if (mi->rigidbody != NULL) {
				if (rbw->shared->physics_world && mi->rigidbody && mi->rigidbody->shared->physics_object)
					RB_dworld_remove_body(rbw->shared->physics_world, mi->rigidbody->shared->physics_object);
				if (mi->rigidbody->shared->physics_object) {
					RB_body_delete(mi->rigidbody->shared->physics_object);
					mi->rigidbody->shared->physics_object = NULL;
				}

				if (mi->rigidbody->shared->physics_shape) {
					RB_shape_delete(mi->rigidbody->shared->physics_shape);
					mi->rigidbody->shared->physics_shape = NULL;
				}

				MEM_freeN(mi->rigidbody->shared);
				MEM_freeN(mi->rigidbody);
				mi->rigidbody = NULL;
			}
		}
	}

	return modFound;
}

#endif
