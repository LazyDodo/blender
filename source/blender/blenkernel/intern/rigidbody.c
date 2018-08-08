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
 * The Original Code is Copyright (C) 2013 Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Joshua Leung, Sergej Reich, Martin Felke
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file rigidbody.c
 *  \ingroup blenkernel
 *  \brief Blender-side interface and methods for dealing with Rigid Body simulations
 */

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_callbacks.h"
#include "BLI_math.h"
#include "BLI_kdtree.h"
#include "BLI_rand.h"
#include "BLI_sort.h"
#include "BLI_threads.h"
#include "BLI_utildefines.h"

#ifdef WITH_BULLET
#  include "RBI_api.h"
#endif

#include "DEG_depsgraph.h"

#include "DNA_fracture_types.h"
#include "DNA_ID.h"
#include "DNA_group_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_object_force_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

#include "BKE_collection.h"
#include "BKE_effect.h"
#include "BKE_fracture.h"
#include "BKE_global.h"
#include "BKE_layer.h"
#include "BKE_library.h"
#include "BKE_library_query.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"
#include "BKE_object.h"
#include "BKE_pointcache.h"
#include "BKE_rigidbody.h"
#include "BKE_modifier.h"
#include "BKE_scene.h"
#include "PIL_time.h"
#include "bmesh.h"


#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

/* ************************************** */
/* Memory Management */

/* Freeing Methods --------------------- */

#ifdef WITH_BULLET
//static void rigidbody_update_ob_array(RigidBodyWorld *rbw);

#else
static void RB_dworld_remove_constraint(void *UNUSED(world), void *UNUSED(con)) {}
static void RB_dworld_remove_body(void *UNUSED(world), void *UNUSED(body)) {}
static void RB_dworld_delete(void *UNUSED(world)) {}
static void RB_body_delete(void *UNUSED(body)) {}
static void RB_shape_delete(void *UNUSED(shape)) {}
static void RB_constraint_delete(void *UNUSED(con)) {}

#endif

/* Free rigidbody world */
void BKE_rigidbody_free_world(Scene *scene)
{
	bool is_orig = (scene->id.tag & LIB_TAG_COPIED_ON_WRITE) == 0;
	RigidBodyWorld *rbw = scene->rigidbody_world;
	scene->rigidbody_world = NULL;

	/* sanity check */
	if (!rbw)
		return;

	if (is_orig && rbw->shared->physics_world) {
		/* free physics references, we assume that all physics objects in will have been added to the world */
		if (rbw->constraints) {
			FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rbw->constraints, object)
			{
				if (object->rigidbody_constraint) {
					RigidBodyCon *rbc = object->rigidbody_constraint;
					if (rbc->physics_constraint) {
						RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
					}
				}
			}
			FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
		}

		if (rbw->group) {
			FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rbw->group, object)
			{
				BKE_rigidbody_free_object(object, rbw);
			}
			FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
		}
		/* free dynamics world */
		RB_dworld_delete(rbw->shared->physics_world);
	}

	if (rbw->objects)
		MEM_freeN(rbw->objects);

	if (rbw->cache_index_map) {
		MEM_freeN(rbw->cache_index_map);
		rbw->cache_index_map = NULL;
	}

	if (rbw->cache_offset_map) {
		MEM_freeN(rbw->cache_offset_map);
		rbw->cache_offset_map = NULL;
	}


	if (is_orig) {
		/* free cache */
		BKE_ptcache_free_list(&(rbw->shared->ptcaches));
		rbw->shared->pointcache = NULL;

		MEM_freeN(rbw->shared);
	}

	/* free effector weights */
	if (rbw->effector_weights)
		MEM_freeN(rbw->effector_weights);

	/* free rigidbody world itself */
	MEM_freeN(rbw);
}

/* Free RigidBody settings and sim instances */
void BKE_rigidbody_free_object(Object *ob, RigidBodyWorld *rbw)
{
	bool is_orig = (ob->id.tag & LIB_TAG_COPIED_ON_WRITE) == 0;
	RigidBodyOb *rbo = (ob) ? ob->rigidbody_object : NULL;

	/* sanity check */
	if (rbo == NULL)
		return;

	/* free physics references */
	if (is_orig) {
		if (rbo->shared->physics_object) {
			BLI_assert(rbw);
			if (rbw) {
				/* We can only remove the body from the world if the world is known.
				 * The world is generally only unknown if it's an evaluated copy of
				 * an object that's being freed, in which case this code isn't run anyway. */
				RB_dworld_remove_body(rbw->shared->physics_world, rbo->shared->physics_object);
			}

			RB_body_delete(rbo->shared->physics_object);
			rbo->shared->physics_object = NULL;
		}

		if (rbo->shared->physics_shape) {
			RB_shape_delete(rbo->shared->physics_shape);
			rbo->shared->physics_shape = NULL;
		}

		MEM_freeN(rbo->shared);
	}

	/* free data itself */
	MEM_freeN(rbo);
	ob->rigidbody_object = NULL;
}

/* Free RigidBody constraint and sim instance */
void BKE_rigidbody_free_constraint(Object *ob)
{
	RigidBodyCon *rbc = (ob) ? ob->rigidbody_constraint : NULL;

	/* sanity check */
	if (rbc == NULL)
		return;

	/* free physics reference */
	if (rbc->physics_constraint) {
		RB_constraint_delete(rbc->physics_constraint);
		rbc->physics_constraint = NULL;
	}

	/* free data itself */
	MEM_freeN(rbc);
	ob->rigidbody_constraint = NULL;
}

#ifdef WITH_BULLET

/* Copying Methods --------------------- */

/* These just copy the data, clearing out references to physics objects.
 * Anything that uses them MUST verify that the copied object will
 * be added to relevant groups later...
 */

RigidBodyOb *BKE_rigidbody_copy_object(const Object *ob, const int flag)
{
	RigidBodyOb *rboN = NULL;

	if (ob->rigidbody_object) {
		/* just duplicate the whole struct first (to catch all the settings) */
		rboN = MEM_dupallocN(ob->rigidbody_object);

		if ((flag & LIB_ID_CREATE_NO_MAIN) == 0) {
			/* This is a regular copy, and not a CoW copy for depsgraph evaluation */
			rboN->shared = MEM_callocN(sizeof(*rboN->shared), "RigidBodyOb_Shared");
		}

		/* tag object as needing to be verified */
		rboN->flag |= RBO_FLAG_NEEDS_VALIDATE;
	}

	/* return new copy of settings */
	return rboN;
}

RigidBodyCon *BKE_rigidbody_copy_constraint(const Object *ob, const int UNUSED(flag))
{
	RigidBodyCon *rbcN = NULL;

	if (ob->rigidbody_constraint) {
		/* just duplicate the whole struct first (to catch all the settings) */
		rbcN = MEM_dupallocN(ob->rigidbody_constraint);

		/* tag object as needing to be verified */
		rbcN->flag |= RBC_FLAG_NEEDS_VALIDATE;

		/* clear out all the fields which need to be revalidated later */
		rbcN->physics_constraint = NULL;
	}

	/* return new copy of settings */
	return rbcN;
}

/* ************************************** */
/* Setup Utilities - Validate Sim Instances */

/* get the appropriate evaluated mesh based on rigid body mesh source */
static Mesh *rigidbody_get_mesh(Object *ob)
{
	switch (ob->rigidbody_object->mesh_source) {
		case RBO_MESH_DEFORM:
			return ob->runtime.mesh_deform_eval;
		case RBO_MESH_FINAL:
			return ob->runtime.mesh_eval;
		case RBO_MESH_BASE:
			/* This mesh may be used for computing looptris, which should be done
			 * on the original; otherwise every time the CoW is recreated it will
			 * have to be recomputed. */
			BLI_assert(ob->rigidbody_object->mesh_source == RBO_MESH_BASE);
			return ob->runtime.mesh_orig;
	}

	/* Just return something sensible so that at least Blender won't crash. */
	BLI_assert(!"Unknown mesh source");
	return ob->runtime.mesh_eval;
}

rbCollisionShape *BKE_rigidbody_get_shape_convexhull_from_mesh(Mesh *dm, float margin, bool *can_embed)
{
	rbCollisionShape *shape = NULL;
	int totvert = dm->totvert;
	MVert *mvert = dm->mvert;

	if (dm && totvert) {
		shape = RB_shape_new_convex_hull((float *)mvert, sizeof(MVert), totvert, margin, can_embed);
	}
	else {
		printf("ERROR: no vertices to define Convex Hull collision shape with\n");
	}

	return shape;
}

/* create collision shape of mesh - triangulated mesh
 * returns NULL if creation fails.
 */
rbCollisionShape *BKE_rigidbody_get_shape_trimesh_from_mesh(Object *ob, MeshIsland* mi)
{
	rbCollisionShape *shape = NULL;
	Mesh *dm = NULL;
	MVert *mvert;
	const MLoopTri *looptri;
	int totvert;
	int tottri;
	const MLoop *mloop;

	if (mi && mi->physics_mesh) {
		dm = mi->physics_mesh;
	}
	else if (ob->type == OB_MESH) {
		dm = rigidbody_get_mesh(ob);
	}
	else {
		printf("ERROR: cannot make Triangular Mesh collision shape for non-Mesh object\n");
		return NULL;
	}

	/* ensure mesh validity, then grab data */
	if (dm == NULL)
		return NULL;

	mvert   = dm->mvert;
	totvert = dm->totvert;
	looptri = BKE_mesh_runtime_looptri_ensure(dm);
	tottri = BKE_mesh_runtime_looptri_len(dm);
	mloop = dm->mloop;

	/* sanity checking - potential case when no data will be present */
	if ((totvert == 0) || (tottri == 0)) {
		printf("WARNING: no geometry data converted for Mesh Collision Shape (ob = %s)\n", ob->id.name + 2);
	}
	else {
		rbMeshData *mdata;
		int i;

		/* init mesh data for collision shape */
		mdata = RB_trimesh_data_new(tottri, totvert);

		RB_trimesh_add_vertices(mdata, (float *)mvert, totvert, sizeof(MVert));

		/* loop over all faces, adding them as triangles to the collision shape
		 * (so for some faces, more than triangle will get added)
		 */
		if (mvert && looptri) {
			for (i = 0; i < tottri; i++) {
				/* add first triangle - verts 1,2,3 */
				const MLoopTri *lt = &looptri[i];
				int vtri[3];

				vtri[0] = mloop[lt->tri[0]].v;
				vtri[1] = mloop[lt->tri[1]].v;
				vtri[2] = mloop[lt->tri[2]].v;

				RB_trimesh_add_triangle_indices(mdata, i, UNPACK3(vtri));
			}
		}

		RB_trimesh_finish(mdata);

		/* construct collision shape
		 *
		 * These have been chosen to get better speed/accuracy tradeoffs with regards
		 * to limitations of each:
		 *    - BVH-Triangle Mesh: for passive objects only. Despite having greater
		 *                         speed/accuracy, they cannot be used for moving objects.
		 *    - GImpact Mesh:      for active objects. These are slower and less stable,
		 *                         but are more flexible for general usage.
		 */
		if (ob->rigidbody_object->type == RBO_TYPE_PASSIVE) {
			shape = RB_shape_new_trimesh(mdata);
		}
		else {
			shape = RB_shape_new_gimpact_mesh(mdata);
		}
	}

	return shape;
}

/* Create new physics sim collision shape for object and store it,
 * or remove the existing one first and replace...
 */
static void rigidbody_validate_sim_shape(Object *ob, bool rebuild)
{
	RigidBodyOb *rbo = ob->rigidbody_object;
	rbCollisionShape *new_shape = NULL;
	BoundBox *bb = NULL;
	float size[3] = {1.0f, 1.0f, 1.0f};
	float radius = 1.0f;
	float height = 1.0f;
	float capsule_height;
	float hull_margin = 0.0f;
	bool can_embed = true;
	bool has_volume;

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
	 */
	// XXX: all dimensions are auto-determined now... later can add stored settings for this
	/* get object dimensions without scaling */
	bb = BKE_object_boundbox_get(ob);
	if (bb) {
		size[0] = (bb->vec[4][0] - bb->vec[0][0]);
		size[1] = (bb->vec[2][1] - bb->vec[0][1]);
		size[2] = (bb->vec[1][2] - bb->vec[0][2]);
	}
	mul_v3_fl(size, 0.5f);

	if (ELEM(rbo->shape, RB_SHAPE_CAPSULE, RB_SHAPE_CYLINDER, RB_SHAPE_CONE)) {
		/* take radius as largest x/y dimension, and height as z-dimension */
		radius = MAX2(size[0], size[1]);
		height = size[2];
	}
	else if (rbo->shape == RB_SHAPE_SPHERE) {
		/* take radius to the largest dimension to try and encompass everything */
		radius = MAX3(size[0], size[1], size[2]);
	}

	/* create new shape */
	switch (rbo->shape) {
		case RB_SHAPE_BOX:
			new_shape = RB_shape_new_box(size[0], size[1], size[2]);
			break;

		case RB_SHAPE_SPHERE:
			new_shape = RB_shape_new_sphere(radius + RBO_GET_MARGIN(rbo));
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
			/* try to emged collision margin */
			has_volume = (MIN3(size[0], size[1], size[2]) > 0.0f);

			if (!(rbo->flag & RBO_FLAG_USE_MARGIN) && has_volume)
				hull_margin = 0.04f;
			if (ob->type == OB_MESH && ob->data) {
				new_shape = BKE_rigidbody_get_shape_convexhull_from_mesh((Mesh *)ob->data, hull_margin, &can_embed);
			}
			else {
				printf("ERROR: cannot make Convex Hull collision shape for non-Mesh object\n");
			}

			if (!(rbo->flag & RBO_FLAG_USE_MARGIN))
				rbo->margin = (can_embed && has_volume) ? 0.04f : 0.0f;  /* RB_TODO ideally we shouldn't directly change the margin here */
			break;
		case RB_SHAPE_TRIMESH:
			new_shape = BKE_rigidbody_get_shape_trimesh_from_mesh(ob, NULL);
			break;
	}
	/* use box shape if we can't fall back to old shape */
	if (new_shape == NULL && rbo->shared->physics_shape == NULL) {
		new_shape = RB_shape_new_box(size[0], size[1], size[2]);
	}
	/* assign new collision shape if creation was successful */
	if (new_shape) {
		if (rbo->shared->physics_shape) {
			RB_shape_delete(rbo->shared->physics_shape);
		}
		rbo->shared->physics_shape = new_shape;
		RB_shape_set_margin(rbo->shared->physics_shape, RBO_GET_MARGIN(rbo));
	}
}

/* --------------------- */

/* helper function to calculate volume of rigidbody object */
// TODO: allow a parameter to specify method used to calculate this?
void BKE_rigidbody_calc_volume(Object *ob, float *r_vol)
{
	RigidBodyOb *rbo = ob->rigidbody_object;

	float size[3]  = {1.0f, 1.0f, 1.0f};
	float radius = 1.0f;
	float height = 1.0f;

	float volume = 0.0f;

	/* if automatically determining dimensions, use the Object's boundbox
	 *	- assume that all quadrics are standing upright on local z-axis
	 *	- assume even distribution of mass around the Object's pivot
	 *	  (i.e. Object pivot is centralized in boundbox)
	 *	- boundbox gives full width
	 */
	// XXX: all dimensions are auto-determined now... later can add stored settings for this
	BKE_object_dimensions_get(ob, size);

	if (ELEM(rbo->shape, RB_SHAPE_CAPSULE, RB_SHAPE_CYLINDER, RB_SHAPE_CONE)) {
		/* take radius as largest x/y dimension, and height as z-dimension */
		radius = MAX2(size[0], size[1]) * 0.5f;
		height = size[2];
	}
	else if (rbo->shape == RB_SHAPE_SPHERE) {
		/* take radius to the largest dimension to try and encompass everything */
		radius = max_fff(size[0], size[1], size[2]) * 0.5f;
	}

	/* calculate volume as appropriate  */
	switch (rbo->shape) {
		case RB_SHAPE_BOX:
			volume = size[0] * size[1] * size[2];
			break;

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

		case RB_SHAPE_CONVEXH:
		case RB_SHAPE_TRIMESH:
		{
			if (ob->type == OB_MESH) {
				Mesh *mesh = rigidbody_get_mesh(ob);
				MVert *mvert;
				const MLoopTri *lt = NULL;
				int totvert, tottri = 0;
				const MLoop *mloop = NULL;

				/* ensure mesh validity, then grab data */
				if (mesh == NULL)
					return;

				mvert   = mesh->mvert;
				totvert = mesh->totvert;
				lt = BKE_mesh_runtime_looptri_ensure(mesh);
				tottri = mesh->runtime.looptris.len;
				mloop = mesh->mloop;

				if (totvert > 0 && tottri > 0) {
					BKE_mesh_calc_volume(mvert, totvert, lt, tottri, mloop, &volume, NULL);
				}
			}
			else {
				/* rough estimate from boundbox as fallback */
				/* XXX could implement other types of geometry here (curves, etc.) */
				volume = size[0] * size[1] * size[2];
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
	if (r_vol) *r_vol = volume;
}

void BKE_rigidbody_calc_center_of_mass(Object *ob, float r_center[3])
{
	RigidBodyOb *rbo = ob->rigidbody_object;

	float size[3]  = {1.0f, 1.0f, 1.0f};
	float height = 1.0f;

	zero_v3(r_center);

	/* if automatically determining dimensions, use the Object's boundbox
	 *	- assume that all quadrics are standing upright on local z-axis
	 *	- assume even distribution of mass around the Object's pivot
	 *	  (i.e. Object pivot is centralized in boundbox)
	 *	- boundbox gives full width
	 */
	// XXX: all dimensions are auto-determined now... later can add stored settings for this
	BKE_object_dimensions_get(ob, size);

	/* calculate volume as appropriate  */
	switch (rbo->shape) {
		case RB_SHAPE_BOX:
		case RB_SHAPE_SPHERE:
		case RB_SHAPE_CAPSULE:
		case RB_SHAPE_CYLINDER:
			break;

		case RB_SHAPE_CONE:
			/* take radius as largest x/y dimension, and height as z-dimension */
			height = size[2];
			/* cone is geometrically centered on the median,
			 * center of mass is 1/4 up from the base
			 */
			r_center[2] = -0.25f * height;
			break;

		case RB_SHAPE_CONVEXH:
		case RB_SHAPE_TRIMESH:
		{
			if (ob->type == OB_MESH) {
				Mesh *mesh = rigidbody_get_mesh(ob);
				MVert *mvert;
				const MLoopTri *looptri;
				int totvert, tottri;
				const MLoop *mloop;

				/* ensure mesh validity, then grab data */
				if (mesh == NULL)
					return;

				mvert   = mesh->mvert;
				totvert = mesh->totvert;
				looptri = BKE_mesh_runtime_looptri_ensure(mesh);
				tottri = mesh->runtime.looptris.len;
				mloop = mesh->mloop;

				if (totvert > 0 && tottri > 0) {
					BKE_mesh_calc_volume(mvert, totvert, looptri, tottri, mloop, NULL, r_center);
				}
			}
			break;
		}

#if 0 // XXX: not defined yet
		case RB_SHAPE_COMPOUND:
			volume = 0.0f;
			break;
#endif
	}
}

/* --------------------- */

/**
 * Create physics sim representation of object given RigidBody settings
 *
 * \param rebuild Even if an instance already exists, replace it
 */
static void rigidbody_validate_sim_object(RigidBodyWorld *rbw, Object *ob, bool rebuild, bool transfer_speeds)
{
	RigidBodyOb *rbo = (ob) ? ob->rigidbody_object : NULL;
	float loc[3];
	float rot[4];
	float size[3] = {1.0f, 1.0f, 1.0f};

	/* sanity checks:
	 *	- object doesn't have RigidBody info already: then why is it here?
	 */
	if (rbo == NULL)
		return;

	/* make sure collision shape exists */
	/* FIXME we shouldn't always have to rebuild collision shapes when rebuilding objects, but it's needed for constraints to update correctly */
	if (rbo->shared->physics_shape == NULL || rebuild)
		rigidbody_validate_sim_shape(ob, true);

	if (rbo->shared->physics_object) {
		RB_dworld_remove_body(rbw->shared->physics_world, rbo->shared->physics_object);
	}
	if (!rbo->shared->physics_object || rebuild) {
		/* remove rigid body if it already exists before creating a new one */
		if (rbo->shared->physics_object) {
			RB_body_delete(rbo->shared->physics_object);
		}

		mat4_to_loc_quat(loc, rot, ob->obmat);

		rbo->shared->physics_object = RB_body_new(rbo->shared->physics_shape, loc, rot, false, 0.0f, 0.0f, 0.0f, 0.0f, size);

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
	}

	if (rbw && rbw->shared->physics_world)
		RB_dworld_add_body(rbw->shared->physics_world, rbo->shared->physics_object, rbo->col_groups, NULL, ob, -1);

}

/* --------------------- */

static void rigidbody_constraint_init_spring(
		RigidBodyCon *rbc, void (*set_spring)(rbConstraint *, int, int),
		void (*set_stiffness)(rbConstraint *, int, float), void (*set_damping)(rbConstraint *, int, float))
{
	set_spring(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->flag & RBC_FLAG_USE_SPRING_X);
	set_stiffness(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->spring_stiffness_x);
	set_damping(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->spring_damping_x);

	set_spring(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->flag & RBC_FLAG_USE_SPRING_Y);
	set_stiffness(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->spring_stiffness_y);
	set_damping(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->spring_damping_y);

	set_spring(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->flag & RBC_FLAG_USE_SPRING_Z);
	set_stiffness(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->spring_stiffness_z);
	set_damping(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->spring_damping_z);

	set_spring(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->flag & RBC_FLAG_USE_SPRING_ANG_X);
	set_stiffness(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->spring_stiffness_ang_x);
	set_damping(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->spring_damping_ang_x);

	set_spring(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->flag & RBC_FLAG_USE_SPRING_ANG_Y);
	set_stiffness(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->spring_stiffness_ang_y);
	set_damping(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->spring_damping_ang_y);

	set_spring(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->flag & RBC_FLAG_USE_SPRING_ANG_Z);
	set_stiffness(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->spring_stiffness_ang_z);
	set_damping(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->spring_damping_ang_z);
}

static void rigidbody_constraint_set_limits(
		RigidBodyCon *rbc, void (*set_limits)(rbConstraint *, int, float, float))
{
	if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_X)
		set_limits(rbc->physics_constraint, RB_LIMIT_LIN_X, rbc->limit_lin_x_lower, rbc->limit_lin_x_upper);
	else
		set_limits(rbc->physics_constraint, RB_LIMIT_LIN_X, 0.0f, -1.0f);

	if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_Y)
		set_limits(rbc->physics_constraint, RB_LIMIT_LIN_Y, rbc->limit_lin_y_lower, rbc->limit_lin_y_upper);
	else
		set_limits(rbc->physics_constraint, RB_LIMIT_LIN_Y, 0.0f, -1.0f);

	if (rbc->flag & RBC_FLAG_USE_LIMIT_LIN_Z)
		set_limits(rbc->physics_constraint, RB_LIMIT_LIN_Z, rbc->limit_lin_z_lower, rbc->limit_lin_z_upper);
	else
		set_limits(rbc->physics_constraint, RB_LIMIT_LIN_Z, 0.0f, -1.0f);

	if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_X)
		set_limits(rbc->physics_constraint, RB_LIMIT_ANG_X, rbc->limit_ang_x_lower, rbc->limit_ang_x_upper);
	else
		set_limits(rbc->physics_constraint, RB_LIMIT_ANG_X, 0.0f, -1.0f);

	if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_Y)
		set_limits(rbc->physics_constraint, RB_LIMIT_ANG_Y, rbc->limit_ang_y_lower, rbc->limit_ang_y_upper);
	else
		set_limits(rbc->physics_constraint, RB_LIMIT_ANG_Y, 0.0f, -1.0f);

	if (rbc->flag & RBC_FLAG_USE_LIMIT_ANG_Z)
		set_limits(rbc->physics_constraint, RB_LIMIT_ANG_Z, rbc->limit_ang_z_lower, rbc->limit_ang_z_upper);
	else
		set_limits(rbc->physics_constraint, RB_LIMIT_ANG_Z, 0.0f, -1.0f);
}

/**
 * Create physics sim representation of constraint given rigid body constraint settings
 *
 * < rebuild: even if an instance already exists, replace it
 */
static void rigidbody_validate_sim_constraint(Scene* scene, Object *ob, bool rebuild)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	RigidBodyCon *rbc = (ob) ? ob->rigidbody_constraint : NULL;
	float loc[3];
	float rot[4];
	float lin_lower;
	float lin_upper;
	float ang_lower;
	float ang_upper;

	/* sanity checks:
	 *	- object should have a rigid body constraint
	 *  - rigid body constraint should have at least one constrained object
	 */
	if (rbc == NULL) {
		return;
	}

	if (ELEM(NULL, rbc->ob1, rbc->ob1->rigidbody_object, rbc->ob2, rbc->ob2->rigidbody_object)) {
		if (rbc->physics_constraint) {
			RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
			RB_constraint_delete(rbc->physics_constraint);
			rbc->physics_constraint = NULL;
		}
		return;
	}

	if (rbc->physics_constraint && rebuild == false) {
		RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
	}
	if (rbc->physics_constraint == NULL || rebuild) {
		rbRigidBody *rb1 = NULL, *rb2 = NULL;
		FractureModifierData *fmd1 = NULL, *fmd2 = NULL;
		MeshIsland *mi1, *mi2;

		/*add logic to find closest shards of 2 FM objects, or to attach one shard to one regular */
		/*for now, just take the first meshislands... later find the closest ones or the closest to center*/
		fmd1 = (FractureModifierData*)modifiers_findByType(rbc->ob1, eModifierType_Fracture);
		fmd2 = (FractureModifierData*)modifiers_findByType(rbc->ob2, eModifierType_Fracture);

		if (fmd1 && fmd2) {
			mi1 = BKE_rigidbody_closest_meshisland_to_point(fmd1, rbc->ob1, ob, scene, rbc);
			mi2 = BKE_rigidbody_closest_meshisland_to_point(fmd2, rbc->ob2, ob, scene, rbc);

			if (mi1 && mi2) {
				rb1 = mi1->rigidbody->shared->physics_object;
				rb2 = mi2->rigidbody->shared->physics_object;
			}
		}
		else if (fmd1) {
			mi1 = BKE_rigidbody_closest_meshisland_to_point(fmd1, rbc->ob1, ob, scene, rbc);
			if (mi1) {
				rb1 = mi1->rigidbody->shared->physics_object;
				rb2 = rbc->ob2->rigidbody_object->shared->physics_object;
			}
		}
		else if (fmd2) {
			mi2 = BKE_rigidbody_closest_meshisland_to_point(fmd2, rbc->ob2, ob, scene, rbc);
			if (mi2)
			{
				rb2 = mi2->rigidbody->shared->physics_object;
				rb1 = rbc->ob1->rigidbody_object->shared->physics_object;
			}
		}
		else {
			rb1 = rbc->ob1->rigidbody_object->shared->physics_object;
			rb2 = rbc->ob2->rigidbody_object->shared->physics_object;
		}

		/* remove constraint if it already exists before creating a new one */
		if (rbc->physics_constraint) {
			RB_constraint_delete(rbc->physics_constraint);
			rbc->physics_constraint = NULL;
		}

		mat4_to_loc_quat(loc, rot, ob->obmat);

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
					if (rbc->spring_type == RBC_SPRING_TYPE2) {
						rbc->physics_constraint = RB_constraint_new_6dof_spring2(loc, rot, rb1, rb2);

						rigidbody_constraint_init_spring(rbc, RB_constraint_set_spring_6dof_spring2, RB_constraint_set_stiffness_6dof_spring2, RB_constraint_set_damping_6dof_spring2);

						RB_constraint_set_equilibrium_6dof_spring2(rbc->physics_constraint);

						rigidbody_constraint_set_limits(rbc, RB_constraint_set_limits_6dof_spring2);
					}
					else {
						rbc->physics_constraint = RB_constraint_new_6dof_spring(loc, rot, rb1, rb2);

						rigidbody_constraint_init_spring(rbc, RB_constraint_set_spring_6dof_spring, RB_constraint_set_stiffness_6dof_spring, RB_constraint_set_damping_6dof_spring);

						RB_constraint_set_equilibrium_6dof_spring(rbc->physics_constraint);

						rigidbody_constraint_set_limits(rbc, RB_constraint_set_limits_6dof);
					}
					break;
				case RBC_TYPE_6DOF:
					rbc->physics_constraint = RB_constraint_new_6dof(loc, rot, rb1, rb2);

					rigidbody_constraint_set_limits(rbc, RB_constraint_set_limits_6dof);
					break;
				case RBC_TYPE_MOTOR:
					rbc->physics_constraint = RB_constraint_new_motor(loc, rot, rb1, rb2);

					RB_constraint_set_enable_motor(rbc->physics_constraint, rbc->flag & RBC_FLAG_USE_MOTOR_LIN, rbc->flag & RBC_FLAG_USE_MOTOR_ANG);
					RB_constraint_set_max_impulse_motor(rbc->physics_constraint, rbc->motor_lin_max_impulse, rbc->motor_ang_max_impulse);
					RB_constraint_set_target_velocity_motor(rbc->physics_constraint, rbc->motor_lin_target_velocity, rbc->motor_ang_target_velocity);
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
	}

	if (rbw && rbw->shared->physics_world && rbc->physics_constraint) {
		RB_dworld_add_constraint(rbw->shared->physics_world, rbc->physics_constraint, rbc->flag & RBC_FLAG_DISABLE_COLLISIONS);
	}

	if (rbc->physics_constraint)
	{
		/*
		char id[64];
		char *rest = id;
		char *ptr;
		char ptr2[64] = "0";

		strncpy(id, ob->id.name + 2, strlen(ob->id.name) - 2);

		ptr = strtok_r(id, ".", &rest);
		while (ptr != NULL)
		{
			strncpy(ptr2, ptr, strlen(ptr));
			ptr = strtok_r(NULL, ".", &rest);
		}*/

		RB_constraint_set_id(rbc->physics_constraint, ob->id.name + 2);
	}
}



/* --------------------- */

/* Create physics sim world given RigidBody world settings */
// NOTE: this does NOT update object references that the scene uses, in case those aren't ready yet!
void BKE_rigidbody_validate_sim_world(Scene *scene, RigidBodyWorld *rbw, bool rebuild)
{
	/* sanity checks */
	if (rbw == NULL)
		return;

	/* create new sim world */
	if (rebuild || rbw->shared->physics_world == NULL) {
		if (rbw->shared->physics_world)
			RB_dworld_delete(rbw->shared->physics_world);
		rbw->shared->physics_world = RB_dworld_new(scene->physics_settings.gravity, rbw, scene,
												   BKE_rigidbody_filter_callback, BKE_rigidbody_contact_callback,
												   BKE_rigidbody_id_callback, NULL);
	}

	RB_dworld_set_solver_iterations(rbw->shared->physics_world, rbw->num_solver_iterations);
	RB_dworld_set_split_impulse(rbw->shared->physics_world, rbw->flag & RBW_FLAG_USE_SPLIT_IMPULSE);
}

/* ************************************** */
/* Setup Utilities - Create Settings Blocks */

/* Set up RigidBody world */
RigidBodyWorld *BKE_rigidbody_create_world(Scene *scene)
{
	/* try to get whatever RigidBody world that might be representing this already */
	RigidBodyWorld *rbw;

	/* sanity checks
	 *	- there must be a valid scene to add world to
	 *	- there mustn't be a sim world using this group already
	 */
	if (scene == NULL)
		return NULL;

	/* create a new sim world */
	rbw = MEM_callocN(sizeof(RigidBodyWorld), "RigidBodyWorld");
	rbw->shared = MEM_callocN(sizeof(*rbw->shared), "RigidBodyWorld_Shared");

	/* set default settings */
	rbw->effector_weights = BKE_add_effector_weights(NULL);

	rbw->ltime = PSFRA;

	rbw->time_scale = 1.0f;

	rbw->steps_per_second = 60; /* Bullet default (60 Hz) */
	rbw->num_solver_iterations = 10; /* 10 is bullet default */

	rbw->shared->pointcache = BKE_ptcache_add(&(rbw->shared->ptcaches));
	rbw->shared->pointcache->step = 1;
	rbw->flag &=~ RBW_FLAG_OBJECT_CHANGED;
	rbw->flag &=~ RBW_FLAG_REFRESH_MODIFIERS;

	rbw->objects = MEM_mallocN(sizeof(Object *), "objects");
	rbw->cache_index_map = MEM_mallocN(sizeof(RigidBodyOb *), "cache_index_map");
	rbw->cache_offset_map = MEM_mallocN(sizeof(int), "cache_offset_map");

	/* return this sim world */
	return rbw;
}

/* Add rigid body settings to the specified shard */
RigidBodyOb *BKE_rigidbody_create_shard(Main* bmain, Scene *scene, Object *ob, Object *target, MeshIsland *mi)
{
	RigidBodyOb *rbo;
	RigidBodyWorld *rbw = BKE_rigidbody_get_world(scene);
	float centr[3], size[3], mat[4][4];

	/* sanity checks
	 *	- rigidbody world must exist
	 *	- shard must exist
	 *	- cannot add rigid body if it already exists
	 */
	if (mi == NULL || (mi->rigidbody != NULL))
		return NULL;

	if (ob->type != OB_MESH && ob->type != OB_FONT && ob->type != OB_CURVE && ob->type != OB_SURF) {
		return NULL;
	}

	if ((ob->type == OB_MESH) && (((Mesh *)ob->data)->totvert == 0)) {
		return NULL;
	}

	/* Add rigid body world and group if they don't exist for convenience */
	if (rbw == NULL) {
		rbw = BKE_rigidbody_create_world(scene);
		BKE_rigidbody_validate_sim_world(scene, rbw, false);
		scene->rigidbody_world = rbw;
	}
	if (rbw->group == NULL) {
		rbw->group = BKE_collection_add(bmain, NULL, "RigidBodyWorld");
	}

	/* make rigidbody object settings */
	if (ob->rigidbody_object == NULL) {
		ob->rigidbody_object = BKE_rigidbody_create_object(scene, ob, RBO_TYPE_ACTIVE, NULL);
	}
	else {
		ob->rigidbody_object->type = RBO_TYPE_ACTIVE;
		ob->rigidbody_object->flag |= RBO_FLAG_NEEDS_VALIDATE;
	}

	if (!BKE_collection_has_object(rbw->group, ob)) {
		BKE_collection_object_add(bmain, rbw->group, ob);
	}

		DEG_id_tag_update(&ob->id, OB_RECALC_OB);

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
		rbo->type = mi->ground_weight > 0.01f ? RBO_TYPE_PASSIVE : RBO_TYPE_ACTIVE;

		/* set initial transform */
		mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);
		mat4_to_size(size, ob->obmat);

		//add initial "offset" (centroid), maybe subtract ob->obmat ?? (not sure)
		copy_v3_v3(centr, mi->centroid);
		mul_v3_v3(centr, size);
		mul_qt_v3(rbo->orn, centr);
		add_v3_v3(rbo->pos, centr);
		zero_v3(rbo->lin_vel);
		zero_v3(rbo->ang_vel);
	}

	/* return this object */
	return rbo;
}

RigidBodyWorld *BKE_rigidbody_world_copy(RigidBodyWorld *rbw, const int flag)
{
	RigidBodyWorld *rbw_copy = MEM_dupallocN(rbw);

	if (rbw->effector_weights) {
		rbw_copy->effector_weights = MEM_dupallocN(rbw->effector_weights);
	}
	if ((flag & LIB_ID_CREATE_NO_USER_REFCOUNT) == 0) {
		id_us_plus((ID *)rbw_copy->group);
		id_us_plus((ID *)rbw_copy->constraints);
	}

	if ((flag & LIB_ID_CREATE_NO_MAIN) == 0) {
		/* This is a regular copy, and not a CoW copy for depsgraph evaluation */
		rbw_copy->shared = MEM_callocN(sizeof(*rbw_copy->shared), "RigidBodyWorld_Shared");
		BKE_ptcache_copy_list(&rbw_copy->shared->ptcaches, &rbw->shared->ptcaches, LIB_ID_COPY_CACHES);
		rbw_copy->shared->pointcache = rbw_copy->shared->ptcaches.first;
	}

	rbw_copy->objects = NULL;
	rbw_copy->numbodies = 0;
	BKE_rigidbody_update_ob_array(rbw_copy, false);

	return rbw_copy;
}

void BKE_rigidbody_world_groups_relink(RigidBodyWorld *rbw)
{
	ID_NEW_REMAP(rbw->group);
	ID_NEW_REMAP(rbw->constraints);
	ID_NEW_REMAP(rbw->effector_weights->group);
}

void BKE_rigidbody_world_id_loop(RigidBodyWorld *rbw, RigidbodyWorldIDFunc func, void *userdata)
{
	CollectionObject *go;

	func(rbw, (ID **)&rbw->group, userdata, IDWALK_CB_NOP);
	func(rbw, (ID **)&rbw->constraints, userdata, IDWALK_CB_NOP);
	func(rbw, (ID **)&rbw->effector_weights->group, userdata, IDWALK_CB_NOP);

	/* in regular blender one rigidbody is equivalent to one object, but in FM build you can have more
	 * rigidbodies in the world than objects in the object array... thats why loop over the group of objects
	 * of the rigidbody world... the rbw->objects array would contain group objects only anyway, also dont forget
	 * regular constraints if there are any */

	/*if (rbw->objects) {
		int i;
		int count = BLI_listbase_count(&rbw->group->gobject);
		for (i = 0; i < count; i++) {
			func(rbw, (ID **)&rbw->objects[i], userdata, IDWALK_CB_NOP);
		}
	}*/

	if (rbw->group) {
		for (go = rbw->group->gobject.first; go; go = go->next) {
			func(rbw, (ID **)&go->ob, userdata, IDWALK_CB_NOP);
		}
	}

	if (rbw->constraints) {
		for (go = rbw->constraints->gobject.first; go; go = go->next) {
			func(rbw, (ID **)&go->ob, userdata, IDWALK_CB_NOP);
		}
	}
}

/* Add rigid body settings to the specified object */
RigidBodyOb *BKE_rigidbody_create_object(Scene *scene, Object *ob, short type, MeshIsland *mi)
{
	RigidBodyOb *rbo;
	RigidBodyWorld *rbw = scene->rigidbody_world;
	FractureModifierData *fmd = NULL;

	/* sanity checks
	 *	- rigidbody world must exist
	 *	- object must exist
	 *	- cannot add rigid body if it already exists
	 */
	if (ob == NULL || (ob->rigidbody_object != NULL))
		return NULL;

	/* create new settings data, and link it up */
	rbo = MEM_callocN(sizeof(RigidBodyOb), "RigidBodyOb");
	rbo->shared = MEM_callocN(sizeof(*rbo->shared), "RigidBodyOb_Shared");

	if (mi != NULL && mi->rigidbody != NULL)
	{
		rbo->flag = mi->rigidbody->flag;
		rbo->shared->physics_object = NULL;
		rbo->shared->physics_shape = NULL;

		rbo->flag |= RBO_FLAG_NEEDS_VALIDATE;

		rbo->type = mi->rigidbody->type;
		rbo->mass = mi->rigidbody->mass;

		rbo->friction = mi->rigidbody->friction;
		rbo->restitution = mi->rigidbody->restitution;

		rbo->margin = mi->rigidbody->margin;

		rbo->lin_sleep_thresh = mi->rigidbody->lin_sleep_thresh;
		rbo->ang_sleep_thresh = mi->rigidbody->ang_sleep_thresh;
		rbo->force_thresh = mi->rigidbody->force_thresh;

		rbo->lin_damping = mi->rigidbody->lin_damping;
		rbo->ang_damping = mi->rigidbody->ang_damping;

		rbo->col_groups = mi->rigidbody->col_groups;

		rbo->shape = mi->rigidbody->shape;
		rbo->mesh_source = mi->rigidbody->mesh_source;
		rbo->meshisland_index = mi->rigidbody->meshisland_index;
		rbo->is_fractured = true;
		copy_v3_v3(rbo->pos, mi->rigidbody->pos);
		copy_qt_qt(rbo->orn, mi->rigidbody->orn);
		//mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);
	}
	else
	{
		/* set default settings */
		rbo->type = type;

		rbo->mass = 1.0f;

		rbo->friction = 0.5f; /* best when non-zero. 0.5 is Bullet default */
		rbo->restitution = 0.0f; /* best when zero. 0.0 is Bullet default */

		rbo->margin = 0.04f; /* 0.04 (in meters) is Bullet default */

		rbo->lin_sleep_thresh = 0.4f; /* 0.4 is half of Bullet default */
		rbo->ang_sleep_thresh = 0.5f; /* 0.5 is half of Bullet default */
		rbo->force_thresh = 0.0f; /*dont activate by force by default */

		rbo->lin_damping = 0.04f; /* 0.04 is game engine default */
		rbo->ang_damping = 0.1f; /* 0.1 is game engine default */

		rbo->col_groups = 1;

		/* use triangle meshes for passive objects
		 * use convex hulls for active objects since dynamic triangle meshes are very unstable
		 */
		if (type == RBO_TYPE_ACTIVE)
			rbo->shape = RB_SHAPE_CONVEXH;
		else
			rbo->shape = RB_SHAPE_TRIMESH;

		rbo->mesh_source = RBO_MESH_DEFORM;

		/* set initial transform */
		mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);

		rbo->meshisland_index = -1;
		rbo->is_fractured = false;
	}

	zero_v3(rbo->lin_vel);
	zero_v3(rbo->ang_vel);

	fmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);
	if (fmd && fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{	//keep cache here
		return rbo;
	}

	/* flag cache as outdated */
	BKE_rigidbody_cache_reset(rbw);

	/* return this object */
	return rbo;
}

/* Add rigid body constraint to the specified object */
RigidBodyCon *BKE_rigidbody_create_constraint(Scene *scene, Object *ob, short type, RigidBodyShardCon *con)
{
	RigidBodyCon *rbc;
	RigidBodyWorld *rbw = scene->rigidbody_world;

	/* sanity checks
	 *	- rigidbody world must exist
	 *	- object must exist
	 *	- cannot add constraint if it already exists
	 */
	if (ob == NULL || (ob->rigidbody_constraint != NULL))
		return NULL;

	/* create new settings data, and link it up */
	rbc = MEM_callocN(sizeof(RigidBodyCon), "RigidBodyCon");

	/* set default settings */
	rbc->type = type;

	rbc->ob1 = NULL;
	rbc->ob2 = NULL;

	rbc->flag |= RBC_FLAG_ENABLED;
	rbc->flag |= RBC_FLAG_DISABLE_COLLISIONS;

	if (con)
	{
		rbc->flag = con->flag;
		rbc->breaking_threshold = con->breaking_threshold;
		rbc->num_solver_iterations = con->num_solver_iterations;

		rbc->limit_lin_x_lower = con->limit_lin_x_lower;
		rbc->limit_lin_x_upper = con->limit_lin_x_upper;
		rbc->limit_lin_y_lower = con->limit_lin_y_lower;
		rbc->limit_lin_y_upper = con->limit_lin_y_upper;
		rbc->limit_lin_z_lower = con->limit_lin_z_lower;
		rbc->limit_lin_z_upper = con->limit_lin_z_upper;
		rbc->limit_ang_x_lower = con->limit_ang_x_lower;
		rbc->limit_ang_x_upper = con->limit_ang_x_upper;
		rbc->limit_ang_y_lower = con->limit_ang_y_lower;
		rbc->limit_ang_y_upper = con->limit_ang_y_upper;
		rbc->limit_ang_z_lower = con->limit_ang_z_lower;
		rbc->limit_ang_z_upper = con->limit_ang_z_upper;

		rbc->spring_damping_x = con->spring_damping_x;
		rbc->spring_damping_y = con->spring_damping_y;
		rbc->spring_damping_z = con->spring_damping_z;
		rbc->spring_damping_ang_x = con->spring_damping_ang_x;
		rbc->spring_damping_ang_y = con->spring_damping_ang_y;
		rbc->spring_damping_ang_z = con->spring_damping_ang_z;
		rbc->spring_stiffness_x = con->spring_stiffness_x;
		rbc->spring_stiffness_y = con->spring_stiffness_y;
		rbc->spring_stiffness_z = con->spring_stiffness_z;
		rbc->spring_stiffness_ang_x = con->spring_stiffness_ang_x;
		rbc->spring_stiffness_ang_y = con->spring_stiffness_ang_y;
		rbc->spring_stiffness_ang_z = con->spring_stiffness_ang_z;

		rbc->motor_lin_max_impulse = con->motor_lin_max_impulse;
		rbc->motor_lin_target_velocity = con->motor_lin_target_velocity;
		rbc->motor_ang_max_impulse = con->motor_ang_max_impulse;
		rbc->motor_ang_target_velocity = con->motor_ang_target_velocity;
	}
	else
	{
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
	}

	/* flag cache as outdated */
	BKE_rigidbody_cache_reset(rbw);

	/* return this object */
	return rbc;
}

/* Add rigid body constraint to the specified object */
RigidBodyShardCon *BKE_rigidbody_create_shard_constraint(Scene *scene, short type, bool reset)
{
	RigidBodyShardCon *rbc;
	RigidBodyWorld *rbw = scene->rigidbody_world;

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
		BKE_rigidbody_cache_reset(rbw);

	/* return this object */
	return rbc;
}

/* ************************************** */
/* Utilities API */

/* Get RigidBody world for the given scene, creating one if needed
 *
 * \param scene Scene to find active Rigid Body world for
 */
RigidBodyWorld *BKE_rigidbody_get_world(Scene *scene)
{
	/* sanity check */
	if (scene == NULL)
		return NULL;

	return scene->rigidbody_world;
}

void BKE_rigidbody_remove_shard_con(Scene *scene, RigidBodyShardCon *con)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
		if (rbw && rbw->shared->physics_world && con && con->physics_constraint) {
			RB_dworld_remove_constraint(rbw->shared->physics_world, con->physics_constraint);
		RB_constraint_delete(con->physics_constraint);
		con->physics_constraint = NULL;
	}
}

void BKE_rigidbody_remove_shard(Scene *scene, MeshIsland *mi)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	int i = 0;

	/* rbw can be NULL directly after linking / appending objects without their original scenes
	 * if an attempt to refracture is done then, this would crash here with null pointer access */
	if (mi->rigidbody != NULL && rbw != NULL) {

		RigidBodyShardCon *con;

		for (i = 0; i < mi->participating_constraint_count; i++) {
			con = mi->participating_constraints[i];
			BKE_rigidbody_remove_shard_con(scene, con);
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

		/* this SHOULD be the correct global index, mark with NULL as 'dirty' BEFORE deleting */
		/* need to check whether we didnt create the rigidbody world manually already, prior to fracture, in this
		 * case cache_index_map might be not initialized ! checking numbodies here, they should be 0 in a fresh
		 * rigidbody world */
#if 0
		if ((rbw->cache_index_map != NULL) && (rbw->numbodies > 0) && mi->linear_index < rbw->numbodies) {
			//mi->rigidbody = NULL;
			rbw->cache_index_map[mi->linear_index] = NULL;
		}
#endif

		if (rbw->cache_index_map != NULL) {
			MEM_freeN(rbw->cache_index_map);
			rbw->cache_index_map = NULL;
		}

		//BKE_rigidbody_update_ob_array(rbw);
	}
}

static bool do_remove_modifier(RigidBodyWorld* rbw, ModifierData *md, Object *ob)
{
	RigidBodyShardCon *con;
	MeshIsland *mi;
	FractureModifierData *fmd;
	bool modFound = false;

	if (md->type == eModifierType_Fracture)
	{
		fmd = (FractureModifierData *)md;
		modFound = true;
		CollectionObject *go;

		for (con = fmd->shared->meshConstraints.first; con; con = con->next) {
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
				for (con = fmdi->shared->meshConstraints.first; con; con = con->next) {
					if (rbw && rbw->shared->physics_world && con->physics_constraint) {
						RB_dworld_remove_constraint(rbw->shared->physics_world, con->physics_constraint);
						RB_constraint_delete(con->physics_constraint);
						con->physics_constraint = NULL;
					}
				}
			}
		}


		for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
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

				if (rbw->cache_index_map != NULL) {
					MEM_freeN(rbw->cache_index_map);
					rbw->cache_index_map = NULL;
				}

				MEM_freeN(mi->rigidbody);
				mi->rigidbody = NULL;
			}
		}
	}

	return modFound;
}

void BKE_rigidbody_remove_object(Main *bmain, Scene *scene, Object *ob)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	RigidBodyOb *rbo = ob->rigidbody_object;
	RigidBodyCon *rbc;
	ModifierData *md;
	int i;
	bool modFound = false;

	if (rbw) {
		for (md = ob->modifiers.first; md; md = md->next) {
			modFound = do_remove_modifier(rbw, md, ob);
		}

		if (!modFound) {
			/* remove from rigidbody world, free object won't do this */
					if (rbw->shared->physics_world && rbo->shared->physics_object)
					RB_dworld_remove_body(rbw->shared->physics_world, rbo->shared->physics_object);

			/* remove object from array */
			if (rbw && rbw->objects) {
				for (i = 0; i < rbw->numbodies; i++) {
					int index = rbw->cache_offset_map[i];
					if (rbw->objects[index] == ob) {
						rbw->objects[index] = NULL;
					}

					if (rbw->cache_index_map != NULL) {
						MEM_freeN(rbw->cache_index_map);
						rbw->cache_index_map = NULL;
					}
				}
			}

			/* remove object from rigid body constraints */
			if (rbw->constraints) {
				FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN(rbw->constraints, obt)
				{
					if (obt && obt->rigidbody_constraint) {
						rbc = obt->rigidbody_constraint;
						if (ELEM(ob, rbc->ob1, rbc->ob2)) {
							BKE_rigidbody_remove_constraint(scene, obt);
						}
					}
				}
				FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
			}
			BKE_collection_object_remove(bmain, rbw->group, ob, false);
		}

		/* remove object's settings */
		BKE_rigidbody_free_object(ob, rbw);

		/* flag cache as outdated */
		BKE_rigidbody_cache_reset(rbw);
	}
}

void BKE_rigidbody_remove_constraint(Scene *scene, Object *ob)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	RigidBodyCon *rbc = ob->rigidbody_constraint;

	/* remove from rigidbody world, free object won't do this */
	if (rbw && rbw->shared->physics_world && rbc->physics_constraint) {
		RB_dworld_remove_constraint(rbw->shared->physics_world, rbc->physics_constraint);
	}
	/* remove object's settings */
	BKE_rigidbody_free_constraint(ob);

	/* flag cache as outdated */
	BKE_rigidbody_cache_reset(rbw);
}

static int rigidbody_group_count_items(const ListBase *group, int *r_num_objects, int *r_num_shards)
{
	int num_gobjects = 0;
	ModifierData *md;
	FractureModifierData *rmd;
		CollectionObject *gob;

	if (r_num_objects == NULL || r_num_shards == NULL)
	{
		return num_gobjects;
	}

	*r_num_objects = 0;
	*r_num_shards = 0;

	for (gob = group->first; gob; gob = gob->next) {
		bool found_modifiers = false;
		for (md = gob->ob->modifiers.first; md; md = md->next) {
			if (md->type == eModifierType_Fracture) {
				rmd = (FractureModifierData *)md;
				if (BKE_rigidbody_modifier_active(rmd))
				{
					found_modifiers = true;
					if (rmd->shared->meshIslands.first != NULL)
					{
						*r_num_shards += BLI_listbase_count(&rmd->shared->meshIslands);
					}
				}
			}
		}
		if (found_modifiers == false && gob->ob->rigidbody_object) {
			(*r_num_objects)++;
		}
		num_gobjects++;
	}

	return num_gobjects;
}

static int object_sort_eval(const void *s1, const void *s2, void* context)
{
	Object **o1 = (Object**)s1;
	Object **o2 = (Object**)s2;

	FractureModifierData *fmd1 = (FractureModifierData*)modifiers_findByType(*o1, eModifierType_Fracture);
	FractureModifierData *fmd2 = (FractureModifierData*)modifiers_findByType(*o2, eModifierType_Fracture);


	if ((fmd1 && fmd1->dm_group && fmd1->use_constraint_group) &&
	   (fmd2 && fmd2->dm_group && fmd2->use_constraint_group))
	{
		return 0;
	}

	if ((fmd1 && fmd1->dm_group && fmd1->use_constraint_group) &&
	   (fmd2 && !(fmd2->dm_group && fmd2->use_constraint_group)))
	{
		return 1;
	}

	if ((fmd1 && fmd1->dm_group && fmd1->use_constraint_group) && !fmd2)
	{
		return 1;
	}

	if ((fmd1 && !(fmd1->dm_group && fmd1->use_constraint_group)) &&
	   (fmd2 && (fmd2->dm_group && fmd2->use_constraint_group)))
	{
		return -1;
	}

	if ((!fmd1) && (fmd2 && fmd2->dm_group && fmd2->use_constraint_group))
	{
		return -1;
	}

	return 0;
}

/* ************************************** */
/* Simulation Interface - Bullet */

/* Update object array and rigid body count so they're in sync with the rigid body group */
void BKE_rigidbody_update_ob_array(RigidBodyWorld *rbw, bool do_bake_correction)
{
	CollectionObject *go;
	FractureModifierData *rmd;
	MeshIsland *mi;
	int i, j = 0, l = 0, m = 0, n = 0, counter = 0;
	bool ismapped = false;
	RigidBodyOb ** tmp_index = NULL;
	int *tmp_offset = NULL;

	if (rbw->objects != NULL) {
		MEM_freeN(rbw->objects);
		rbw->objects = NULL;
	}

	if (rbw->cache_index_map != NULL) {
		MEM_freeN(rbw->cache_index_map);
		rbw->cache_index_map = NULL;
	}

	if (rbw->cache_offset_map != NULL) {
		MEM_freeN(rbw->cache_offset_map);
		rbw->cache_offset_map = NULL;
	}

	l = rigidbody_group_count_items(&rbw->group->gobject, &m, &n);

	rbw->numbodies = m + n;
	rbw->objects = MEM_mallocN(sizeof(Object *) * l, "objects");
	rbw->cache_index_map = MEM_mallocN(sizeof(RigidBodyOb *) * rbw->numbodies, "cache_index_map");
	rbw->cache_offset_map = MEM_mallocN(sizeof(int) * rbw->numbodies, "cache_offset_map");
	tmp_index =  MEM_mallocN(sizeof(RigidBodyOb *) * rbw->numbodies, "cache_index_map temp");
	tmp_offset = MEM_mallocN(sizeof(int) * rbw->numbodies, "cache_offset_map");

	printf("RigidbodyCount changed: %d\n", rbw->numbodies);

	//pre-sort rbw->objects... put such with fmd->dm_group and fmd->use_constraint_group after all without,
	//to enforce proper eval order

	for (go = rbw->group->gobject.first, i = 0; go; go = go->next, i++) {

		if (go->ob->rigidbody_object)
		{
			rbw->objects[i] = go->ob;
		}
	}

	BLI_qsort_r(rbw->objects, l, sizeof(Object *), object_sort_eval, NULL);

	//correct map if baked, it might be shifted
	for (i = 0; i < l; i++) {
		Object *ob = rbw->objects[i];
		printf("%s\n", ob->id.name + 2);

		rmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);
		if (rmd) {
			for (mi = rmd->shared->meshIslands.first, j = 0; mi; mi = mi->next) {
				//store original position of the object in the object array, to be able to rearrange it later so it matches the baked cache
				tmp_index[counter] = mi->rigidbody; /* map all shards of an object to this object index*/
				tmp_offset[counter] = i;
				mi->linear_index = counter;
				if (mi->rigidbody && !do_bake_correction) {
					//as we search by id now in the pointcache, we set the id here too
					mi->rigidbody->meshisland_index = counter;
					if (mi->object_index != -1)
					{
						mi->object_index = i;
					}
				}
				counter++;
				j++;
			}
			ismapped = true;
		}

		if (!ismapped && ob->rigidbody_object) {
			tmp_index[counter] = ob->rigidbody_object; /*1 object 1 index here (normal case)*/
			tmp_offset[counter] = i;
			if (ob->rigidbody_object && !do_bake_correction)
				ob->rigidbody_object->meshisland_index = counter;
			counter++;
		}

		ismapped = false;
	}

	//do shuffle for each rigidbody
	for (i = 0; i < rbw->numbodies; i++)
	{
		RigidBodyOb *rbo = tmp_index[i];
		int offset = tmp_offset[i];

		if (rbo->meshisland_index != i && do_bake_correction)
		{
			rbw->cache_index_map[rbo->meshisland_index] = rbo;
			rbw->cache_offset_map[rbo->meshisland_index] = offset;
			//printf("Shuffle %d -> %d\n", i, rbo->meshisland_index);
		}
		else {
			rbw->cache_index_map[i] = rbo;
			rbw->cache_offset_map[i] = offset;
		}
	}

	MEM_freeN(tmp_index);
	MEM_freeN(tmp_offset);
}

static void rigidbody_update_sim_world(Scene *scene, RigidBodyWorld *rbw, bool rebuild)
{
	float adj_gravity[3];
	bool skip_correction = rbw->flag & RBW_FLAG_REFRESH_MODIFIERS;

	/* adjust gravity to take effector weights into account */
	if (scene->physics_settings.flag & PHYS_GLOBAL_GRAVITY) {
		copy_v3_v3(adj_gravity, scene->physics_settings.gravity);
		mul_v3_fl(adj_gravity, rbw->effector_weights->global_gravity * rbw->effector_weights->weight[0]);
	}
	else {
		zero_v3(adj_gravity);
	}

	/* update gravity, since this RNA setting is not part of RigidBody settings */
	RB_dworld_set_gravity(rbw->shared->physics_world, adj_gravity);

	/* update object array in case there are changes */
	if (rebuild)
	{
		BKE_rigidbody_update_ob_array(rbw, (rbw->shared->pointcache->flag & PTCACHE_BAKED) && !skip_correction);
	}
}

void BKE_rigidbody_update_sim_ob(Scene *scene, RigidBodyWorld *rbw, Object *ob, RigidBodyOb *rbo, float centroid[3],
									MeshIsland *mi, float size[3], FractureModifierData *fmd, Depsgraph *depsgraph)
{
	float loc[3];
	float rot[4];
	float scale[3], centr[3];

	/* only update if rigid body exists */
	if (rbo->shared->physics_object == NULL)
		return;

	if (rbo->shape == RB_SHAPE_TRIMESH && rbo->flag & RBO_FLAG_USE_DEFORM) {
		Mesh *dm = NULL;

		if (rbo->mesh_source == RBO_MESH_DEFORM) {
			dm = ob->runtime.mesh_deform_eval;
		}
		else if (rbo->mesh_source == RBO_MESH_FINAL) {
			dm = ob->runtime.mesh_eval;
		}

		if (dm) {
			MVert *mvert = dm->mvert;
			int totvert = dm->totvert;
			BoundBox *bb = BKE_object_boundbox_get(ob);

			if (RB_shape_get_num_verts(rbo->shared->physics_shape) != totvert)
			{
				if (mi != NULL)
				{
					//fracture modifier case TODO, update mi->physicsmesh somehow and redraw
					rbo->flag |= RBO_FLAG_NEEDS_RESHAPE;
					BKE_rigidbody_shard_validate(rbw, mi, ob, false, false, size);
				}
				else
				{
					//regular rigidbody case
					rigidbody_validate_sim_shape(ob, true);
					RB_body_set_collision_shape(rbo->shared->physics_object, rbo->shared->physics_shape);
				}
			}
			else
			{
				RB_shape_trimesh_update(rbo->shared->physics_shape, (float *)mvert, totvert, sizeof(MVert), bb->vec[0],
										bb->vec[6]);
			}
		}
	}

	copy_v3_v3(centr, centroid);
	mat4_decompose(loc, rot, scale, ob->obmat);
	mul_v3_v3(scale, size);

	/* update scale for all objects */
	RB_body_set_scale(rbo->shared->physics_object, scale);
	/* compensate for embedded convex hull collision margin */
	if (!(rbo->flag & RBO_FLAG_USE_MARGIN) && rbo->shape == RB_SHAPE_CONVEXH) {
		RB_shape_set_margin(rbo->shared->physics_shape, RBO_GET_MARGIN(rbo) * MIN3(scale[0], scale[1], scale[2]));
	}

	/* make transformed objects temporarily kinmatic so that they can be moved by the user during simulation */
	if (ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ) {
		RB_body_set_kinematic_state(rbo->shared->physics_object, true);
		RB_body_set_mass(rbo->shared->physics_object, 0.0f);
	}

	/* update rigid body location and rotation for kinematic bodies */
	if ((rbo->flag & RBO_FLAG_KINEMATIC && rbo->force_thresh == 0.0f) || (ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ)) {
		if (((rbo->type == RBO_TYPE_ACTIVE || mi == NULL) && (rbo->flag & RBO_FLAG_KINEMATIC_REBUILD) == 0))
		{
			if (!(fmd && fmd->anim_mesh_ob && fmd->use_animated_mesh))
			{
				//override for externally animated rbs
				mul_v3_v3(centr, scale);
				mul_qt_v3(rot, centr);
				add_v3_v3(loc, centr);
				RB_body_activate(rbo->shared->physics_object);
				RB_body_set_loc_rot(rbo->shared->physics_object, loc, rot);
			}
		}
	}
	/* update influence of effectors - but don't do it on an effector (why not ?)*/
	/* only dynamic bodies need effector update */
	else if (rbo->type == RBO_TYPE_ACTIVE /* && ((ob->pd == NULL) || (ob->pd->forcefield == PFIELD_NULL))*/) {
		EffectorWeights *effector_weights = rbw->effector_weights;
		EffectedPoint epoint = {0};
		ListBase *effectors;

		/* get effectors present in the group specified by effector_weights */
		effectors = BKE_effectors_create(depsgraph, ob, NULL, effector_weights);
		if (effectors) {
			float eff_force[3] = {0.0f, 0.0f, 0.0f};
			float eff_loc[3], eff_vel[3];
			float thresh = rbo->force_thresh * rbo->force_thresh; /*use this to compare against squared length of vector */

			/* create dummy 'point' which represents last known position of object as result of sim */
			// XXX: this can create some inaccuracies with sim position, but is probably better than using unsimulated vals?
			RB_body_get_position(rbo->shared->physics_object, eff_loc);
			//mul_v3_v3(centr, scale);
			//add_v3_v3(eff_loc, centr);

			RB_body_get_linear_velocity(rbo->shared->physics_object, eff_vel);

			pd_point_from_loc(scene, eff_loc, eff_vel, 0, &epoint);

			/* calculate net force of effectors, and apply to sim object
			 *	- we use 'central force' since apply force requires a "relative position" which we don't have...
			 */
				BKE_effectors_apply(effectors, NULL, effector_weights, &epoint, eff_force, NULL);
			if ((rbo->flag & RBO_FLAG_KINEMATIC) && (thresh < len_squared_v3(eff_force)))
			{
				BKE_rigidbody_activate(rbo, rbw, mi, ob);
				RB_body_apply_central_force(rbo->shared->physics_object, eff_force);
			}
			else if (rbo->flag & RBO_FLAG_KINEMATIC)
			{
				if ((rbo->flag & RBO_FLAG_KINEMATIC_REBUILD) == 0)
				{   //XXXXX TODO, maybe this is wrong here
					/* do the same here as above, but here we needed the eff_force value to compare against threshold */

					if (!(fmd && fmd->anim_mesh_ob && fmd->use_animated_mesh))
					{	//same override as above
						mul_v3_v3(centr, scale);
						mul_qt_v3(rot, centr);
						add_v3_v3(loc, centr);
						RB_body_activate(rbo->shared->physics_object);
						RB_body_set_loc_rot(rbo->shared->physics_object, loc, rot);
					}
				}
			}
			else
			{
				if (G.f & G_DEBUG)
					printf("\tapplying force (%f,%f,%f) to '%s'\n", eff_force[0], eff_force[1], eff_force[2], ob->id.name + 2);
				/* activate object in case it is deactivated */
				if (!is_zero_v3(eff_force))
						RB_body_activate(rbo->shared->physics_object);
				RB_body_apply_central_force(rbo->shared->physics_object, eff_force);
			}
		}
		else if (G.f & G_DEBUG)
			printf("\tno forces to apply to '%s'\n", ob->id.name + 2);

		/* cleanup */
		 BKE_effectors_free(effectors);
	}
	/* NOTE: passive objects don't need to be updated since they don't move */

	/* NOTE: no other settings need to be explicitly updated here,
	 * since RNA setters take care of the rest :)
	 */
}

/* Updates and validates world, bodies and shapes.
 * < rebuild: rebuild entire simulation
 */
static void rigidbody_update_simulation(Scene *scene, RigidBodyWorld *rbw, bool rebuild, Depsgraph *depsgraph)
{

	CollectionObject *go;
	bool did_modifier = false;
	float centroid[3] = {0, 0, 0};
	float size[3] = {1.0f, 1.0f, 1.0f};
	float count = BLI_listbase_count(&rbw->group->gobject);
	int i = 0;

	/* update world */
	if (rebuild) {
		BKE_rigidbody_validate_sim_world(scene, rbw, true);
	}

	rigidbody_update_sim_world(scene, rbw, rebuild);

	/* update objects */
	for (i = 0; i < count; i++) {
		Object *ob = rbw->objects[i];

		if (ob && (ob->type == OB_MESH)) {
			did_modifier = BKE_rigidbody_modifier_update(scene, ob, rbw, rebuild, depsgraph);
		}

		if (!did_modifier && ob->type == OB_MESH) {
			/* validate that we've got valid object set up here... */
			RigidBodyOb *rbo = ob->rigidbody_object;
			/* update transformation matrix of the object so we don't get a frame of lag for simple animations */
			BKE_object_where_is_calc(depsgraph, scene, ob);

			if (rbw->flag & RBW_FLAG_REBUILD_CONSTRAINTS && rbo)
			{
				zero_v3(rbo->lin_vel);
				zero_v3(rbo->ang_vel);
			}

			if (rbo == NULL) {
				/* Since this object is included in the sim group but doesn't have
				 * rigid body settings (perhaps it was added manually), add!
				 *	- assume object to be active? That is the default for newly added settings...
				 */
				ob->rigidbody_object = BKE_rigidbody_create_object(scene, ob, RBO_TYPE_ACTIVE, NULL);
				rigidbody_validate_sim_object(rbw, ob, true, true);

				rbo = ob->rigidbody_object;
			}
			else {
				/* perform simulation data updates as tagged */
				/* refresh object... */
				if (rebuild) {
					/* World has been rebuilt so rebuild object */
					/* TODO(Sybren): rigidbody_validate_sim_object() can call rigidbody_validate_sim_shape(),
					 * but neither resets the RBO_FLAG_NEEDS_RESHAPE flag nor calls RB_body_set_collision_shape().
					 * This results in the collision shape being created twice, which is unnecessary. */
					rigidbody_validate_sim_object(rbw, ob, true, true);
				}
				else if (rbo->flag & RBO_FLAG_NEEDS_VALIDATE) {
					rigidbody_validate_sim_object(rbw, ob, false, true);
				}
				/* refresh shape... */
				if (rbo->flag & RBO_FLAG_NEEDS_RESHAPE) {
					/* mesh/shape data changed, so force shape refresh */
					rigidbody_validate_sim_shape(ob, true);
					/* now tell RB sim about it */
					// XXX: we assume that this can only get applied for active/passive shapes that will be included as rigidbodies
					 RB_body_set_collision_shape(rbo->shared->physics_object, rbo->shared->physics_shape);
				}
				rbo->flag &= ~(RBO_FLAG_NEEDS_VALIDATE | RBO_FLAG_NEEDS_RESHAPE);
			}

			/* update simulation object... */
			BKE_rigidbody_update_sim_ob(scene, rbw, ob, rbo, centroid, NULL, size, NULL, depsgraph);
		}
	}

	if (rbw->shared->physics_world && rbw->flag & RBW_FLAG_REBUILD_CONSTRAINTS)
	{
		double start = PIL_check_seconds_timer();
		RB_dworld_init_compounds(rbw->shared->physics_world);
		printf("Building compounds done, %g\n", PIL_check_seconds_timer() - start);
	}

	rbw->flag &= ~RBW_FLAG_REFRESH_MODIFIERS;

	/* update constraints */
	if (rbw->constraints == NULL) /* no constraints, move on */
		return;
	for (go = rbw->constraints->gobject.first; go; go = go->next) {
		Object *ob = go->ob;

		if (ob) {
			/* validate that we've got valid object set up here... */
			RigidBodyCon *rbc = ob->rigidbody_constraint;
			/* update transformation matrix of the object so we don't get a frame of lag for simple animations */
			BKE_object_where_is_calc(depsgraph, scene, ob);

			if (rbc == NULL) {
				/* Since this object is included in the group but doesn't have
				 * constraint settings (perhaps it was added manually), add!
				 */
				ob->rigidbody_constraint = BKE_rigidbody_create_constraint(scene, ob, RBC_TYPE_FIXED, NULL);
				rigidbody_validate_sim_constraint(scene, ob, true);

				rbc = ob->rigidbody_constraint;
			}
			else {
				/* perform simulation data updates as tagged */
				if (rebuild) {
					/* World has been rebuilt so rebuild constraint */
					rigidbody_validate_sim_constraint(scene, ob, true);
				}
				else if (rbc->flag & RBC_FLAG_NEEDS_VALIDATE) {
					rigidbody_validate_sim_constraint(scene, ob, false);
				}
				rbc->flag &= ~RBC_FLAG_NEEDS_VALIDATE;
			}

			if (rbc->physics_constraint)
			{
				RigidBodyOb *rbo1 = rbc->ob1->rigidbody_object;
				RigidBodyOb *rbo2 = rbc->ob2->rigidbody_object;

				if ((rbo1->force_thresh > 0 || rbo2->force_thresh > 0))
				{
					if (RB_constraint_get_applied_impulse(rbc->physics_constraint) >= rbo1->force_thresh + rbo2->force_thresh)
					{
						RB_constraint_set_enabled(rbc->physics_constraint, false);
						BKE_rigidbody_activate(rbo1, rbw, NULL, rbc->ob1);
						BKE_rigidbody_activate(rbo2, rbw, NULL, rbc->ob2);
					}
				}
			}
		}
	}
}

static ThreadMutex post_step_lock = BLI_MUTEX_INITIALIZER;
static void rigidbody_update_simulation_post_step(Depsgraph* depsgraph, RigidBodyWorld *rbw)
{
	CollectionObject *go;
	ModifierData *md;
	FractureModifierData *rmd;
	int modFound = false;
	RigidBodyOb *rbo;
	MeshIsland *mi;

	for (go = rbw->group->gobject.first; go; go = go->next) {

		Object *ob = go->ob;
		//handle fractured rigidbodies, maybe test for psys as well ?
		BLI_mutex_lock(&post_step_lock);
		for (md = ob->modifiers.first; md; md = md->next) {
			if (md->type == eModifierType_Fracture) {
				rmd = (FractureModifierData *)md;
				if (BKE_rigidbody_modifier_active(rmd)) {
					for (mi = rmd->shared->meshIslands.first; mi; mi = mi->next) {
						rbo = mi->rigidbody;
						if (!rbo) continue;
						/* reset kinematic state for transformed objects */
						if (rbo->shared->physics_object && ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ) {
							RB_body_set_kinematic_state(rbo->shared->physics_object, rbo->flag & RBO_FLAG_KINEMATIC ||
														rbo->flag & RBO_FLAG_DISABLED);
							RB_body_set_mass(rbo->shared->physics_object, RBO_GET_MASS(rbo));
							/* deactivate passive objects so they don't interfere with deactivation of active objects */
							if (rbo->type == RBO_TYPE_PASSIVE)
								RB_body_deactivate(rbo->shared->physics_object);
						}

						/* update stored velocities, can be set again after sim rebuild */
						if (rmd->fracture_mode == MOD_FRACTURE_DYNAMIC && rbo->shared->physics_object)
						{
							if (!(rbo->flag & RBO_FLAG_KINEMATIC))
							{
								RB_body_get_linear_velocity(rbo->shared->physics_object, rbo->lin_vel);
								RB_body_get_angular_velocity(rbo->shared->physics_object, rbo->ang_vel);
							}
						}
					}
					modFound = true;
					break;
				}
			}
		}
		BLI_mutex_unlock(&post_step_lock);

		/* handle regular rigidbodies */
		if (ob && ob->rigidbody_object && !modFound) {
			rbo = ob->rigidbody_object;
			/* reset kinematic state for transformed objects */
			if (rbo && (ob->flag & SELECT) && (G.moving & G_TRANSFORM_OBJ)) {
				RB_body_set_kinematic_state(rbo->shared->physics_object, rbo->flag & RBO_FLAG_KINEMATIC ||
											rbo->flag & RBO_FLAG_DISABLED);
				RB_body_set_mass(rbo->shared->physics_object, RBO_GET_MASS(rbo));
				/* deactivate passive objects so they don't interfere with deactivation of active objects */
				if (rbo->type == RBO_TYPE_PASSIVE)
					RB_body_deactivate(rbo->shared->physics_object);
			}

			if (!(rbo->flag & RBO_FLAG_KINEMATIC) && rbo->shared->physics_object)
			{
				RB_body_get_linear_velocity(rbo->shared->physics_object, rbo->lin_vel);
				RB_body_get_angular_velocity(rbo->shared->physics_object, rbo->ang_vel);
			}
		}
		modFound = false;
	}
}

bool BKE_rigidbody_check_sim_running(RigidBodyWorld *rbw, float ctime)
{
	return (rbw && (rbw->flag & RBW_FLAG_MUTED) == 0 && ctime > rbw->shared->pointcache->startframe);
}

/* Sync rigid body and object transformations */
static ThreadMutex modifier_lock = BLI_MUTEX_INITIALIZER;
void BKE_rigidbody_sync_transforms(Scene *scene, Object *ob, float ctime)
{
	RigidBodyWorld *rbw = scene->rigidbody_world; //take later other rbws into account too !
	RigidBodyOb *rbo = NULL;
	ModifierData *md;
	int modFound = false;

	if (rbw == NULL)
		return;

	BLI_mutex_lock(&modifier_lock);
	for (md = ob->modifiers.first; md; md = md->next) {
		modFound = BKE_rigidbody_modifier_sync(md, ob, scene, ctime);
		if (modFound)
			break;
	}
	BLI_mutex_unlock(&modifier_lock);

	if (!modFound)
	{
		rbo = ob->rigidbody_object;

		/* keep original transform for kinematic and passive objects */
		if (ELEM(NULL, rbw, rbo) || (((rbo->flag & RBO_FLAG_KINEMATIC) &&
			((rbo->flag & RBO_FLAG_KINEMATIC_REBUILD) == 0)) || rbo->type == RBO_TYPE_PASSIVE))
		{
			return;
		}

		/* use rigid body transform after cache start frame if objects is not being transformed */
		if (BKE_rigidbody_check_sim_running(rbw, ctime) && !(ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ))
		{
			float mat[4][4], size_mat[4][4], size[3];

			/* keep original transform when the simulation is muted */
			if (rbw->flag & RBW_FLAG_MUTED)
				return;

			normalize_qt(rbo->orn); // RB_TODO investigate why quaternion isn't normalized at this point
			quat_to_mat4(mat, rbo->orn);
			copy_v3_v3(mat[3], rbo->pos);

			mat4_to_size(size, ob->obmat);
			size_to_mat4(size_mat, size);
			mul_m4_m4m4(mat, mat, size_mat);

			copy_m4_m4(ob->obmat, mat);
		}
		/* otherwise set rigid body transform to current obmat */
		else {
			if (ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ)
				rbw->flag |= RBW_FLAG_OBJECT_CHANGED;

			mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);
		}
	}
}

static void do_reset_rigidbody(RigidBodyOb *rbo, Object *ob, MeshIsland* mi, float loc[3],
							  float rot[3], float quat[4], float rotAxis[3], float rotAngle)
{
	bool correct_delta = !(rbo->flag & RBO_FLAG_KINEMATIC || rbo->type == RBO_TYPE_PASSIVE);

	/* return rigid body and object to their initial states */
	copy_v3_v3(rbo->pos, ob->loc);
	if (mi)
		add_v3_v3(rbo->pos, mi->centroid);
	copy_v3_v3(ob->loc, loc);

	if (correct_delta) {
		add_v3_v3(rbo->pos, ob->dloc);
	}

	if (ob->rotmode > 0) {
		float qt[4];
		eulO_to_quat(qt, ob->rot, ob->rotmode);

		if (correct_delta) {
			float dquat[4];
			eulO_to_quat(dquat, ob->drot, ob->rotmode);

			mul_qt_qtqt(rbo->orn, dquat, qt);
		}
		else {
			copy_qt_qt(rbo->orn, qt);
		}

		if (mi)
			mul_qt_qtqt(rbo->orn, qt, mi->rot);

		copy_v3_v3(ob->rot, rot);
	}
	else if (ob->rotmode == ROT_MODE_AXISANGLE) {
		float qt[4];
		axis_angle_to_quat(qt, ob->rotAxis, ob->rotAngle);

		if (correct_delta) {
			float dquat[4];
			axis_angle_to_quat(dquat, ob->drotAxis, ob->drotAngle);

			mul_qt_qtqt(rbo->orn, dquat, qt);
		}
		else {
			copy_qt_qt(rbo->orn, qt);
		}

		if (mi)
			mul_qt_qtqt(rbo->orn, qt, mi->rot);

		copy_v3_v3(ob->rotAxis, rotAxis);
		ob->rotAngle = rotAngle;
	}
	else {
		if (correct_delta) {
			mul_qt_qtqt(rbo->orn, ob->dquat, ob->quat);
		}
		else {
			copy_qt_qt(rbo->orn, ob->quat);
		}

		if (mi)
			mul_qt_qtqt(rbo->orn, rbo->orn, mi->rot);

		copy_qt_qt(ob->quat, quat);
	}

	if (rbo->shared->physics_object) {
		/* allow passive objects to return to original transform */
		if (rbo->type == RBO_TYPE_PASSIVE)
			RB_body_set_kinematic_state(rbo->shared->physics_object, true);
		RB_body_set_loc_rot(rbo->shared->physics_object, rbo->pos, rbo->orn);
	}
}

/* Used when cancelling transforms - return rigidbody and object to initial states */
void BKE_rigidbody_aftertrans_update(Object *ob, float loc[3], float rot[3], float quat[4], float rotAxis[3], float rotAngle,
									 Depsgraph *depsgraph)
{
	RigidBodyOb *rbo;
	ModifierData *md;
	FractureModifierData *rmd;
	float imat[4][4];
	Scene* scene = DEG_get_input_scene(depsgraph); // or the eval one ?

	md = modifiers_findByType(ob, eModifierType_Fracture);
	if (md != NULL)
	{
		MeshIsland *mi;
		rmd = (FractureModifierData *)md;
		invert_m4_m4(imat, ob->obmat);
		for (mi = rmd->shared->meshIslands.first; mi; mi = mi->next)
		{
			rbo = mi->rigidbody;
			do_reset_rigidbody(rbo, ob, mi, loc, rot, quat, rotAxis, rotAngle);
			if (rbo->flag & RBO_FLAG_KINEMATIC)
			{
				BKE_rigidbody_passive_fake_parenting(rmd, ob, rbo, imat);
			}
			else
			{
				BKE_rigidbody_passive_hook(rmd, mi, ob, scene, depsgraph);
			}
		}

		//then update origmat
		copy_m4_m4(rmd->origmat, ob->obmat);
	}
	else {
		rbo = ob->rigidbody_object;
		do_reset_rigidbody(rbo, ob, NULL, loc, rot, quat, rotAxis, rotAngle);

		// RB_TODO update rigid body physics object's loc/rot for dynamic objects here as well (needs to be done outside bullet's update loop)
	}
}

void BKE_rigidbody_cache_reset(RigidBodyWorld *rbw)
{
	if (rbw) {
		rbw->shared->pointcache->flag |= PTCACHE_OUTDATED;
	}
}

/* ------------------ */

/* Rebuild rigid body world */
/* NOTE: this needs to be called before frame update to work correctly */
void BKE_rigidbody_rebuild_world(Depsgraph *depsgraph, Scene *scene, float ctime)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	PointCache *cache;
	PTCacheID pid;
	int startframe, endframe;
	int shards = 0, objects = 0;

	if ((int)ctime == -1)
	{
		if (rbw && rbw->group)
		{
			/*hack to be able to update the simulation data after loading from FM*/
			rigidbody_update_simulation(scene, rbw, true, depsgraph);
		}
		return;
	}

	BKE_ptcache_id_from_rigidbody(&pid, NULL, rbw);
	BKE_ptcache_id_time(&pid, scene, ctime, &startframe, &endframe, NULL);
	cache = rbw->shared->pointcache;

	/* flag cache as outdated if we don't have a world or number of objects in the simulation has changed */
	rigidbody_group_count_items(&rbw->group->gobject, &shards, &objects);
	if (rbw->shared->physics_world == NULL /*|| rbw->numbodies != (shards + objects)*/) {
		cache->flag |= PTCACHE_OUTDATED;
	}

	if (ctime == startframe + 1 && rbw->ltime == startframe) {
		if (cache->flag & PTCACHE_OUTDATED) {

			//if we destroy the cache, also reset dynamic data (if not baked)
			if (!(cache->flag & PTCACHE_BAKED))
			{
				//resetDynamic(rbw, true, false);
			}

			BKE_ptcache_id_reset(scene, &pid, PTCACHE_RESET_OUTDATED);
			rigidbody_update_simulation(scene, rbw, true, depsgraph);
			BKE_ptcache_validate(cache, (int)ctime);
			cache->last_exact = 0;
			cache->flag &= ~PTCACHE_REDO_NEEDED;
		}
	}
}

/* Run RigidBody simulation for the specified physics world */
void BKE_rigidbody_do_simulation(Depsgraph *depsgraph, Scene *scene, float ctime)
{
	float timestep;
	RigidBodyWorld *rbw = scene->rigidbody_world;
	PointCache *cache;
	PTCacheID pid;
	int startframe, endframe;
	bool was_changed = false;

	BKE_ptcache_id_from_rigidbody(&pid, NULL, rbw);
	BKE_ptcache_id_time(&pid, scene, ctime, &startframe, &endframe, NULL);
	cache = rbw->shared->pointcache;

	/*trigger dynamic update*/
	if (rbw->flag & RBW_FLAG_OBJECT_CHANGED)
	{
		rbw->flag &= ~RBW_FLAG_OBJECT_CHANGED;
		was_changed = rbw->flag & RBW_FLAG_REFRESH_MODIFIERS;
		rbw->flag &= ~RBW_FLAG_REFRESH_MODIFIERS;
		//if (!(cache->flag & PTCACHE_BAKED))
		{
			bool baked = cache->flag & PTCACHE_BAKED;
			bool from_cache = cache->last_exact > cache->simframe;
			if ((from_cache || baked) && ctime > startframe + 1)
			{
				/* dont mess with baked data */
				rigidbody_update_simulation(scene, rbw, true, depsgraph);
			}
		}
		//rbw->flag &= ~RBW_FLAG_REFRESH_MODIFIERS;
	}

	if (ctime <= startframe) {
		/* rebuild constraints */
		rbw->flag |= RBW_FLAG_REBUILD_CONSTRAINTS;

		rbw->ltime = startframe;
		if (rbw->flag & RBW_FLAG_OBJECT_CHANGED)
		{       /* flag modifier refresh at their next execution XXX TODO -> still used ? */
			rbw->flag |= RBW_FLAG_REFRESH_MODIFIERS;
			rbw->flag &= ~RBW_FLAG_OBJECT_CHANGED;
			rigidbody_update_simulation(scene, rbw, true, depsgraph);
		}
		return;
	}
	/* make sure we don't go out of cache frame range */
	else if (ctime > endframe) {
		ctime = endframe;
	}

	/* don't try to run the simulation if we don't have a world yet but allow reading baked cache */
	if (rbw->shared->physics_world == NULL && !(cache->flag & PTCACHE_BAKED))
		return;
	else if ((rbw->objects == NULL) || (rbw->cache_index_map == NULL))
	{
		if (!was_changed)
		{
			BKE_rigidbody_cache_reset(rbw);
			BKE_ptcache_id_reset(scene, &pid, PTCACHE_RESET_OUTDATED);
		}

		BKE_rigidbody_update_ob_array(rbw, cache->flag & PTCACHE_BAKED);
	}

	/* try to read from cache */
	// RB_TODO deal with interpolated, old and baked results
	bool can_simulate = (ctime == rbw->ltime + 1) && !(cache->flag & PTCACHE_BAKED);

	if (BKE_ptcache_read(&pid, ctime, can_simulate) == PTCACHE_READ_EXACT) {
		BKE_ptcache_validate(cache, (int)ctime);
		rbw->ltime = ctime;
		return;
	}

	if (!DEG_is_active(depsgraph)) {
		/* When the depsgraph is inactive we should neither write to the cache
		 * nor run the simulation. */
		return;
	}
	else if (rbw->ltime == startframe)
	{
		/*bool did_it = */BKE_restoreKinematic(rbw, false);
		//if (did_it)

		//make 1st run like later runs... hack...
		rigidbody_update_simulation(scene, rbw, true, depsgraph);
	}

	/* advance simulation, we can only step one frame forward */
	if (compare_ff_relative(ctime, rbw->ltime + 1, FLT_EPSILON, 64)) {
		/* write cache for first frame when on second frame */
		if (rbw->ltime == startframe && (cache->flag & PTCACHE_OUTDATED || cache->last_exact == 0)) {
			BKE_ptcache_write(&pid, startframe);
		}

		if (ctime >= startframe) {
			rbw->flag &= ~RBW_FLAG_REBUILD_CONSTRAINTS;
		}

		/* update and validate simulation */
		rigidbody_update_simulation(scene, rbw, false, depsgraph);

		/* calculate how much time elapsed since last step in seconds */
		timestep = 1.0f / (float)FPS * (ctime - rbw->ltime) * rbw->time_scale;
		/* step simulation by the requested timestep, steps per second are adjusted to take time scale into account */
		RB_dworld_step_simulation(rbw->shared->physics_world, timestep, INT_MAX,
								  1.0f / (float)rbw->steps_per_second * min_ff(rbw->time_scale, 1.0f));

		rigidbody_update_simulation_post_step(depsgraph, rbw);

		/* write cache for current frame */
		BKE_ptcache_validate(cache, (int)ctime);
		BKE_ptcache_write(&pid, (unsigned int)ctime);

		rbw->ltime = ctime;
	}
}
/* ************************************** */

#else  /* WITH_BULLET */

/* stubs */
#if defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

struct RigidBodyOb *BKE_rigidbody_copy_object(const Object *ob, const int flag) { return NULL; }
struct RigidBodyCon *BKE_rigidbody_copy_constraint(const Object *ob, const int flag) { return NULL; }
void BKE_rigidbody_validate_sim_world(Scene *scene, RigidBodyWorld *rbw, bool rebuild) {}
void BKE_rigidbody_calc_volume(Object *ob, float *r_vol) { if (r_vol) *r_vol = 0.0f; }
void BKE_rigidbody_calc_center_of_mass(Object *ob, float r_center[3]) { zero_v3(r_center); }
struct RigidBodyWorld *BKE_rigidbody_create_world(Scene *scene) { return NULL; }
struct RigidBodyWorld *BKE_rigidbody_world_copy(RigidBodyWorld *rbw, const int flag) { return NULL; }
void BKE_rigidbody_world_groups_relink(struct RigidBodyWorld *rbw) {}
void BKE_rigidbody_world_id_loop(struct RigidBodyWorld *rbw, RigidbodyWorldIDFunc func, void *userdata) {}
struct RigidBodyOb *BKE_rigidbody_create_object(Scene *scene, Object *ob, short type) { return NULL; }
struct RigidBodyCon *BKE_rigidbody_create_constraint(Scene *scene, Object *ob, short type) { return NULL; }
struct RigidBodyWorld *BKE_rigidbody_get_world(Scene *scene) { return NULL; }
void BKE_rigidbody_remove_object(struct Main *bmain, Scene *scene, Object *ob) {}
void BKE_rigidbody_remove_constraint(Scene *scene, Object *ob) {}
void BKE_rigidbody_sync_transforms(RigidBodyWorld *rbw, Object *ob, float ctime) {}
void BKE_rigidbody_aftertrans_update(Object *ob, float loc[3], float rot[3], float quat[4], float rotAxis[3], float rotAngle) {}
bool BKE_rigidbody_check_sim_running(RigidBodyWorld *rbw, float ctime) { return false; }
void BKE_rigidbody_cache_reset(RigidBodyWorld *rbw) {}
void BKE_rigidbody_rebuild_world(Depsgraph *depsgraph, Scene *scene, float ctime) {}
void BKE_rigidbody_do_simulation(Depsgraph *depsgraph, Scene *scene, float ctime) {}

#if defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic pop
#endif

#endif  /* WITH_BULLET */



/* -------------------- */
/* Depsgraph evaluation */

void BKE_rigidbody_rebuild_sim(Depsgraph *depsgraph,
							   Scene *scene)
{
	float ctime = DEG_get_ctime(depsgraph);
	DEG_debug_print_eval_time(depsgraph, __func__, scene->id.name, scene, ctime);
	/* rebuild sim data (i.e. after resetting to start of timeline) */
	if (BKE_scene_check_rigidbody_active(scene)) {
		BKE_rigidbody_rebuild_world(depsgraph, scene, ctime);
	}
}

void BKE_rigidbody_eval_simulation(Depsgraph *depsgraph,
								   Scene *scene)
{
	float ctime = DEG_get_ctime(depsgraph);
	DEG_debug_print_eval_time(depsgraph, __func__, scene->id.name, scene, ctime);

	/* evaluate rigidbody sim */
	if (!BKE_scene_check_rigidbody_active(scene)) {
		return;
	}
	BKE_rigidbody_do_simulation(depsgraph, scene, ctime);
}

void BKE_rigidbody_object_sync_transforms(Depsgraph *depsgraph,
										  Scene *scene,
										  Object *ob)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	float ctime = DEG_get_ctime(depsgraph);
	DEG_debug_print_eval_time(depsgraph, __func__, ob->id.name, ob, ctime);
	/* read values pushed into RBO from sim/cache... */
	BKE_rigidbody_sync_transforms(scene, ob, ctime);
}

