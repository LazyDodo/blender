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

#include "DEG_depsgraph_query.h"

#ifdef WITH_BULLET
#include "RBI_api.h"
#endif

#include "BKE_fracture.h"
#include "BKE_rigidbody.h"
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
static void check_fracture(rbContactPoint *cp, RigidBodyWorld *rbw, Object *obA, Object *obB);
static MeshIsland* findMeshIsland(FractureModifierData *fmd, int id);
static void test_deactivate_rigidbody(RigidBodyOb *rbo, MeshIsland *mi);
static float box_volume(float size[3]);


void BKE_rigidbody_activate(RigidBodyOb* rbo, RigidBodyWorld *UNUSED(rbw), MeshIsland *mi, Object *ob)
{
	RigidBodyShardCon *con;
	int i;

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

			different_cluster = ((con->mi1->particle_index != con->mi2->particle_index) ||
								((con->mi1->particle_index == -1) && (con->mi2->particle_index == -1)));

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
			for (con = rmd->shared->meshConstraints.first; con; con = con->next) {
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
			for (con = rmd->shared->meshConstraints.first; con; con = con->next) {
				if ((con->mi1 != NULL && con->mi1->rigidbody != NULL) &&
					(con->mi2 != NULL && con->mi2->rigidbody != NULL)) {
					sub_v3_v3v3(con_vec, con->mi1->centroid, con->mi2->centroid);
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
	if ((max_con_mass == 0) && (rmd->use_mass_dependent_thresholds)) {
		return;
	}

	if ((con->mi1 == NULL) || (con->mi2 == NULL)) {
		return;
	}

	max_thresh = thresh = rmd->breaking_threshold;
	if ((con->mi1->rigidbody != NULL) && (con->mi2->rigidbody != NULL)) {

		if (rmd->use_compounds)
		{
			float min_mass = MIN2(con->mi1->rigidbody->mass, con->mi2->rigidbody->mass);
			float max_mass = MAX2(con->mi1->rigidbody->mass, con->mi2->rigidbody->mass);

			thresh = ((min_mass + (rmd->mass_threshold_factor * max_mass)) / (min_mass + max_mass)) * max_thresh;
		}
		else
		{
			con_mass = con->mi1->rigidbody->mass + con->mi2->rigidbody->mass;
			if (rmd->use_mass_dependent_thresholds)
			{
				thresh = (con_mass / max_con_mass) * max_thresh;
			}
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
	BKE_mesh_boundbox_calc(dm, loc, size);
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

void BKE_rigidbody_calc_shard_mass(Object *ob, MeshIsland *mi, Mesh *orig_dm)
{
	Mesh *dm_ob = orig_dm, *dm_mi;
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
		dm_mi = mi->physics_mesh;
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

static void initNormals(struct MeshIsland *mi, Object *ob, FractureModifierData *fmd)
{
	/* hrm have to init Normals HERE, because we cant do this in readfile.c in case the file is loaded (have no access to the Object there) */
	if (mi->vertno == NULL && mi->vertices_cached != NULL) {
		KDTreeNearest n;
		int index = 0, i = 0;
		MVert* mvrt;

		Mesh *dm = ob->runtime.mesh_eval;
		if (dm == NULL) {
			dm = ob->data;
		}

		if (fmd->shared->nor_tree == NULL) {
		/* HRRRRRMMMM need to build the kdtree here as well if we start the sim after loading and not refreshing,
		 * again, no access to object.... */
			int totvert;
			KDTree *tree;
			MVert *mv, *mvert;

			mvert = dm->mvert;
			totvert = dm->totvert;
			tree = BLI_kdtree_new(totvert);

			for (i = 0, mv = mvert; i < totvert; i++, mv++) {
				BLI_kdtree_insert(tree, i, mv->co);
			}

			BLI_kdtree_balance(tree);
			fmd->shared->nor_tree = tree;
		}

		mi->vertno = MEM_callocN(sizeof(short) * 3 * mi->vertex_count, "mi->vertno");
		for (i = 0; i < mi->vertex_count; i++) {
			MVert *v = mi->vertices_cached[i];
			index = BLI_kdtree_find_nearest(fmd->shared->nor_tree, v->co, &n);
			mvrt = dm->mvert + index;
			mi->vertno[i * 3] = mvrt->no[0];
			mi->vertno[i * 3 + 1] = mvrt->no[1];
			mi->vertno[i * 3 + 2] = mvrt->no[2];
		}
	}
}

//should be done in the FM
void BKE_rigidbody_update_cell(struct MeshIsland *mi, Object *ob, float loc[3], float rot[4], FractureModifierData *rmd, int frame,
							   struct RigidBodyWorld *rbw)
{
	float startco[3], centr[3], size[3];
	short startno[3];
	int j, n = 0, x = 0;
	bool invalidData;

	/* hrm have to init Normals HERE, because we cant do this in readfile.c in case the file is loaded (have no access to the Object there)*/
	if (mi->vertno == NULL && rmd->fix_normals) {
		initNormals(mi, ob, rmd);
	}

	invalidData = (loc[0] == FLT_MIN) || (rot[0] == FLT_MIN);

	if (invalidData) {
		return;
	}

	//invert_m4_m4(ob->imat, ob->obmat);
	mat4_to_size(size, ob->obmat);

	if (rmd->fracture_mode != MOD_FRACTURE_DYNAMIC && frame >= mi->start_frame) {
		/*record only in prefracture case here, when you want to convert to keyframes*/
		n = frame - mi->start_frame + 1;
		x = frame - mi->start_frame;

		if (mi->locs == NULL || mi->rots == NULL)
		{
			float loca[3], rota[4], quat[4];
			mi->locs = MEM_mallocN(sizeof(float)*3, "mi->locs");
			mi->rots = MEM_mallocN(sizeof(float)*4, "mi->rots");
			mi->frame_count = 0;

			copy_v3_v3(loca, mi->centroid);
			mul_m4_v3(ob->obmat, loca);
			mat4_to_quat(quat, ob->obmat);

			copy_qt_qt(rota, mi->rot);
			mul_qt_qtqt(rota, quat, rota);

			mi->locs[0] = loca[0];
			mi->locs[1] = loca[1];
			mi->locs[2] = loca[2];

			mi->rots[0] = rota[0];
			mi->rots[1] = rota[1];
			mi->rots[2] = rota[2];
			mi->rots[3] = rota[3];
		}

		if (n > mi->frame_count) {
			mi->locs = MEM_reallocN(mi->locs, sizeof(float) * 3 * n);
			mi->rots = MEM_reallocN(mi->rots, sizeof(float) * 4 * n);

			mi->locs[x*3] = loc[0];
			mi->locs[x*3+1] = loc[1];
			mi->locs[x*3+2] = loc[2];

			mi->rots[x*4] = rot[0];
			mi->rots[x*4+1] = rot[1];
			mi->rots[x*4+2] = rot[2];
			mi->rots[x*4+3] = rot[3];
			mi->frame_count = n;
		}
	}

	for (j = 0; j < mi->vertex_count; j++) {
		struct MVert *vert;
		float fno[3];

		if (!mi->vertices_cached) {
			return;
		}

		vert = mi->vertices_cached[j];
		if (vert == NULL) break;
		if (vert->co == NULL) break;
		//if (rmd->refresh == true) break;

		startco[0] = mi->vertco[j * 3];
		startco[1] = mi->vertco[j * 3 + 1];
		startco[2] = mi->vertco[j * 3 + 2];

		if (rmd->fix_normals) {
			float irot[4], qrot[4];
			startno[0] = mi->vertno[j * 3];
			startno[1] = mi->vertno[j * 3 + 1];
			startno[2] = mi->vertno[j * 3 + 2];

			/*ignore global quaternion rotation here */
			normal_short_to_float_v3(fno, startno);
			mat4_to_quat(qrot, ob->obmat);
			invert_qt_qt(irot, qrot);
			mul_qt_v3(rot, fno);
			mul_qt_v3(irot, fno);
			normal_float_to_short_v3(vert->no, fno);
		}

		copy_v3_v3(vert->co, startco);
		mul_v3_v3(vert->co, size);
		mul_qt_v3(rot, vert->co);
		copy_v3_v3(centr, mi->centroid);
		mul_v3_v3(centr, size);
		mul_qt_v3(rot, centr);
		sub_v3_v3(vert->co, centr);
		add_v3_v3(vert->co, loc);
		mul_m4_v3(ob->imat, vert->co);
	}
}


/* --------------------- */

/* Create new physics sim collision shape for object and store it,
 * or remove the existing one first and replace...
 */
void BKE_rigidbody_validate_sim_shard_shape(MeshIsland *mi, Object *ob, short rebuild)
{
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

	INIT_MINMAX(min, max);
	if (!BKE_mesh_minmax(mi->physics_mesh, min, max)) {
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
			new_shape = BKE_rigidbody_get_shape_convexhull_from_mesh(mi->physics_mesh, hull_margin, &can_embed);
			if (!(rbo->flag & RBO_FLAG_USE_MARGIN))
				rbo->margin = (can_embed && has_volume) ? 0.04f : 0.0f;      /* RB_TODO ideally we shouldn't directly change the margin here */
			break;
		case RB_SHAPE_TRIMESH:
		{
			new_shape = BKE_rigidbody_get_shape_trimesh_from_mesh(ob, mi);
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
}

/* --------------------- */

/* Create physics sim representation of shard given RigidBody settings
 * < rebuild: even if an instance already exists, replace it
 */
void BKE_rigidbody_validate_sim_shard(RigidBodyWorld *rbw, MeshIsland *mi, Object *ob, FractureModifierData *fmd,
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

	/* at validation, reset frame count as well */

	/* make sure collision shape exists */
	/* FIXME we shouldn't always have to rebuild collision shapes when rebuilding objects, but it's needed for constraints to update correctly */
	if (rbo->shared->physics_shape == NULL || rebuild)
		BKE_rigidbody_validate_sim_shard_shape(mi, ob, true);

	if (rbo->shared->physics_object) {
		//if (rebuild == false /*|| mi->rigidbody->flag & RBO_FLAG_KINEMATIC_REBUILD*/)
		RB_dworld_remove_body(rbw->shared->physics_world, rbo->shared->physics_object);
	}

	if (!rbo->shared->physics_object || rebuild /*|| (fmd->use_animated_mesh && fmd->anim_mesh_ob)*/) {
		float size[3];

		/* remove rigid body if it already exists before creating a new one */
		if (rbo->shared->physics_object) {
			RB_body_delete(rbo->shared->physics_object);
		}

		copy_v3_v3(loc, rbo->pos);
		copy_qt_qt(rot, rbo->orn);
		copy_v3_v3(size, isize);

		mul_v3_v3(size, ob->size);
		rbo->shared->physics_object = RB_body_new(rbo->shared->physics_shape, loc, rot, fmd->use_compounds, fmd->impulse_dampening,
												  fmd->directional_factor, fmd->minimum_impulse, fmd->mass_threshold_factor, size);

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
		RB_dworld_add_body(rbw->shared->physics_world, rbo->shared->physics_object, rbo->col_groups, mi, ob, mi->linear_index);
	}

	rbo->flag &= ~RBO_FLAG_NEEDS_VALIDATE;
	//rbo->flag &= ~RBO_FLAG_KINEMATIC_REBUILD;
}

/* --------------------- */

MeshIsland* BKE_rigidbody_closest_meshisland_to_point(FractureModifierData* fmd, Object *ob, Object *ob2, Scene* scene,
													  RigidBodyCon *con)
{
	MeshIsland *mi, **mi_array = NULL;
	KDTree *tree;
	KDTreeNearest *n = NULL;
	int count = 0;
	int index = 0;
	float loc[3], min[3], max[3], vec[3] = {1, 1, 1};
	int i = 0, r = 0;

	count = BLI_listbase_count(&fmd->shared->meshIslands);
	tree = BLI_kdtree_new(count);
	mi_array = MEM_mallocN(sizeof(MeshIsland*) * count, "mi_array find_closest_meshisland");

	for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
		mul_v3_m4v3(loc, ob->obmat, mi->centroid);
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
			MeshIsland* mi2;
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

	if (fmd->fracture_mode == MOD_FRACTURE_EXTERNAL || (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && fmd->is_dynamic_external))
	{
		mul_v3_m4v3(loc, ob->obmat, rbc->pos);
		mat4_to_quat(rot, ob->obmat);
		mul_qt_qtqt(rot, rot, rbc->orn);
	}
	else
	{
		/* keep old constraint calculation for other fracture modes ! */
		/* do this for all constraints */
		/* location for fixed constraints doesnt matter, so keep old setting */
/*		if (rbc->type == RBC_TYPE_FIXED) {
			copy_v3_v3(rbc->pos, rbc->mi1->rigidbody->pos);
		}
		else*/ {
			/* else set location to center */
			add_v3_v3v3(rbc->pos, rbc->mi1->rigidbody->pos, rbc->mi2->rigidbody->pos);
			mul_v3_fl(rbc->pos, 0.5f);
		}

		copy_qt_qt(rbc->orn, rbc->mi1->rigidbody->orn);
		copy_v3_v3(loc, rbc->pos);
		copy_qt_qt(rot, rbc->orn);
	}

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

				if (fmd->fracture_mode == MOD_FRACTURE_EXTERNAL || (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && fmd->is_dynamic_external))
				{
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
				}
				else
				{
					/* no plastic mode available in other fracture modes */
					rigidbody_set_springs_active(rbc, true);
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
		//char id[64];
		//sprintf(id, "%d", rbc->id);
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

static bool do_activate(Object* ob, Object *ob2, MeshIsland *mi_compare, RigidBodyWorld *rbw, MeshIsland *mi_trigger, bool activate)
{
	FractureModifierData *fmd;
	bool valid = true;
	bool antiValid = ob2->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER;
	bool wouldActivate = false;
	MeshIsland *mi;

	fmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);
	valid = valid && (fmd != NULL);
	antiValid = antiValid && (fmd != NULL);

	valid = valid && (ob->rigidbody_object->flag & RBO_FLAG_USE_KINEMATIC_DEACTIVATION);
	valid = valid && ((ob2->rigidbody_object->flag & RBO_FLAG_IS_TRIGGER) || ((ob2->rigidbody_object->flag & RBO_FLAG_PROPAGATE_TRIGGER) &&
			((mi_trigger) && (mi_trigger->rigidbody->flag & RBO_FLAG_PROPAGATE_TRIGGER))));

	if (valid || antiValid)
	{
		for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next)
		{
			bool dissolve = ob->rigidbody_object->flag & RBO_FLAG_CONSTRAINT_DISSOLVE;

			bool same_cluster = ((mi->particle_index != -1) &&
								(mi->particle_index == mi_compare->particle_index));

			bool different_cluster = !same_cluster && dissolve;

			RigidBodyOb* rbo = mi->rigidbody;
			if ((((rbo->flag & RBO_FLAG_KINEMATIC) || different_cluster) &&
				 ((mi_compare == mi) || (same_cluster && !dissolve && fmd->activate_broken))) && valid)
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

static void fake_dynamic_collide(Object *ob1, Object *ob2, MeshIsland *mi1, MeshIsland *mi2, RigidBodyWorld *rbw)
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
				point.contact_body_indexA = mi1->linear_index;
			}
			else {
				copy_v3_v3(point.contact_pos_world_onA, ob1->loc);
				point.contact_body_indexA = -1;
			}

			if (mi2) {
				copy_v3_v3(point.contact_pos_world_onB, mi2->centroid);
				point.contact_body_indexB = mi2->linear_index;
			}
			else {
				copy_v3_v3(point.contact_pos_world_onB, ob2->loc);
				point.contact_body_indexB = -1;
			}

			check_fracture(&point, rbw, ob1, ob2);
		}
	}
}

static bool check_constraint_island(FractureModifierData* fmd, MeshIsland *mi1, MeshIsland *mi2)
{
	if (mi1 && mi2 && !fmd->use_compounds && (!fmd->use_constraint_collision || fmd->use_self_collision)) {

		float dist_sq = len_squared_v3v3(mi1->centroid, mi2->centroid);
		bool is_near = len_squared_v3v3(mi1->rigidbody->pos, mi2->rigidbody->pos) < dist_sq;
		bool same_island = mi1->constraint_index == mi2->constraint_index;
		bool same_near_self = same_island && fmd->use_self_collision && is_near;
		bool diff_clust_near_self = mi1->particle_index != mi2->particle_index && fmd->use_self_collision && is_near;
		bool regular_case = mi1 && mi2 && !fmd->use_constraint_collision && !same_island;

		if (fmd->use_self_collision)
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
		return !(fmd->dm_group && fmd->use_constraint_collision);
	}

	return true;
}

/* this allows partial object activation, only some shards will be activated, called from bullet(!) */
int BKE_rigidbody_filter_callback(void* scene, void* island1, void* island2, void *blenderOb1, void* blenderOb2, bool activate)
{
	//pass scene here, in order to hopefully get the original one from DEG
	MeshIsland* mi1, *mi2;
	Scene *sc = (Scene*)scene, *sc_orig;

	RigidBodyWorld *rbw;
	Object* ob1, *ob2;
	int ob_index1 = -1, ob_index2 = -1;
	bool validOb = true, check_activate = false;

	// oh man... the pleasures of CoW...
	sc_orig = DEG_get_original_id(sc);
	rbw = sc_orig->rigidbody_world;

	mi1 = (MeshIsland*)island1;
	mi2 = (MeshIsland*)island2;

	//MOOOOO
	ob1 = /*DEG_get_original_object(*/(Object*)blenderOb1;
	ob2 = /*DEG_get_original_object(*/(Object*)blenderOb2;

	FractureModifierData *fmd1 = (FractureModifierData*)modifiers_findByType(ob1, eModifierType_Fracture);
	FractureModifierData *fmd2 = (FractureModifierData*)modifiers_findByType(ob2, eModifierType_Fracture);

#if 0
	if ((fmd1 && fmd1->fracture_mode == MOD_FRACTURE_EXTERNAL) ||
	   (fmd2 && fmd2->fracture_mode == MOD_FRACTURE_EXTERNAL))
	{
		/*external doesnt need triggering, maybe the prefractured object (and dynamic ?)... later TODO */
		/* XXXX remove this in case of external, it interferes */
		ob1 = blenderOb1;
		ob2 = blenderOb2;
		return check_colgroup_ghost(ob1, ob2);
	}
#endif

	if (rbw == NULL)
	{
		/* just check for ghost flags here, do not activate anything */
		ob1 = blenderOb1;
		ob2 = blenderOb2;
		return check_colgroup_ghost(ob1, ob2);
	}

	/* cache offset map is a dull name for that... */
	if (mi1 != NULL && rbw->shared->cache_offset_map)
	{
		ob_index1 = rbw->shared->cache_offset_map[mi1->linear_index];
		ob1 = rbw->shared->objects[ob_index1];
	}
	else
	{
		ob1 = blenderOb1;
		ob_index1 = -1;
	}

	if (mi2 != NULL && rbw->shared->cache_offset_map)
	{
		ob_index2 = rbw->shared->cache_offset_map[mi2->linear_index];
		ob2 = rbw->shared->objects[ob_index2];
	}
	else
	{
		ob2 = blenderOb2;
		ob_index2 = -1;
	}

	if ((!ob1 && ob_index1 == -1) || (!ob2 && ob_index2 == -1))
		return false;

	if ((mi1 != NULL) && (mi2 != NULL) && ob_index1 != -1 && ob_index2 != -1) {
		validOb = (ob_index1 != ob_index2 && colgroup_check(ob1->rigidbody_object->col_groups, ob2->rigidbody_object->col_groups) &&
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


		if (ob1->rigidbody_object->flag & RBO_FLAG_USE_KINEMATIC_DEACTIVATION)
		{
			bool override = activate || (ob2->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER);
			check_activate = do_activate(ob1, ob2, mi1, rbw, mi2, override);
		}

		if (ob2->rigidbody_object->flag & RBO_FLAG_USE_KINEMATIC_DEACTIVATION)
		{
			bool override = activate || (ob1->rigidbody_object->flag & RBO_FLAG_ANTI_TRIGGER);
			check_activate = do_activate(ob2, ob1, mi2, rbw, mi1, override);
		}
	}

	//if ghost is involved, and dynafrac trigger is enabled, try to call check_fracture manually here, without forces and with centroid as contact point
	fake_dynamic_collide(ob1, ob2, mi1, mi2, rbw);
	fake_dynamic_collide(ob2, ob1, mi2, mi1, rbw);

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

static Shard* findShard(FractureModifierData *fmd, int id)
{
	Shard *t = fmd->shared->frac_mesh->shard_map.first;
	Shard *s = NULL;

	while (t)
	{
		if (t->shard_id == id && t->flag & SHARD_INTACT)
		{
			//printf("FOUND: %d\n", id);
			s = t;
			break;
		}
		t = t->next;
	}

	return s;
}

static MeshIsland* findMeshIsland(FractureModifierData *fmd, int id)
{
	MeshIsland *mi = fmd->shared->meshIslands.first;

	while (mi)
	{
		if (mi->rigidbody->meshisland_index == id)
		{
			return mi;
		}
		mi = mi->next;
	}

	return NULL;
}

static bool check_shard_size(FractureModifierData *fmd, int id)
{
	FractureID *fid;
	float size = fmd->dynamic_min_size, diff[3];
	Shard *s = NULL;

	s = findShard(fmd, id);

	if (s == NULL)
	{
		return false;
	}

	BKE_shard_calc_minmax(s);

	sub_v3_v3v3(diff, s->max, s->min);
	if (diff[max_axis_v3(diff)] < size)
	{
		return false;
	}

	for (fid = fmd->shared->fracture_ids.first; fid; fid = fid->next)
	{
		if (fid->shardID == id)
		{
			return false;
		}
	}

	printf("FRACTURE : %d\n", id);

	return true;
}

static bool check_constraints(FractureModifierData *fmd, MeshIsland *mi, RigidBodyWorld *rbw) {
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



static void check_fracture(rbContactPoint* cp, RigidBodyWorld *rbw, Object *obA, Object *obB)
{
	int linear_index1, linear_index2;
	Object* ob1 = NULL, *ob2 = NULL;
	int ob_index1, ob_index2;
	FractureModifierData *fmd1, *fmd2;
	float force;

	//printf("checkfracture\n");

	if (cp == NULL)
		return;

	force = cp->contact_force;

	linear_index1 = cp->contact_body_indexA;
	linear_index2 = cp->contact_body_indexB;

	if (rbw == NULL)
	{
		return;
	}

	if (linear_index2 > -1 && linear_index2 < rbw->shared->numbodies)
	{
		ob_index2 = rbw->shared->cache_offset_map[linear_index2];
		ob2 = rbw->shared->objects[ob_index2];
	}
	else if (obB) {
		ob2 = obB;
	}

	if (linear_index1 > -1 && linear_index1 < rbw->shared->numbodies)
	{
		ob_index1 = rbw->shared->cache_offset_map[linear_index1];
		ob1 = rbw->shared->objects[ob_index1];
		fmd1 = (FractureModifierData*)modifiers_findByType(ob1, eModifierType_Fracture);

		if (fmd1 && fmd1->fracture_mode == MOD_FRACTURE_DYNAMIC)
		{
			if (fmd1->shared->current_shard_entry && fmd1->shared->current_shard_entry->is_new)
			{
				RigidBodyOb *rbo = rbw->shared->cache_index_map[linear_index1];
				int id = rbo->meshisland_index;
				MeshIsland* mi = findMeshIsland(fmd1, id);
				Shard *s = findShard(fmd1, mi->id);
				if (mi->rigidbody->mass > 0) {
					force = force / mi->rigidbody->mass;
				}

				//printf("FORCE1:%f\n",force);
				bool canbreak = (force > fmd1->dynamic_force) || ((fmd1->limit_impact || obB) && can_break(ob2, ob1, fmd1->limit_impact));

				if (canbreak && check_constraints(fmd1, mi, rbw))
				{
					if (s) {
						float size[3] =  {1.0f, 1.0f, 1.0f};

						if (ob1 == ob2 || (ob2 && ob2->rigidbody_object && ob1->rigidbody_object->type == RBO_TYPE_PASSIVE)) {
							//todo calculate shard...
							size[0] = size[1] = size[2] = -1.0f;
						}
						else if (ob2) {
							BKE_object_dimensions_get(ob2, size);
						}

						copy_v3_v3(s->impact_loc, cp->contact_pos_world_onA);
						copy_v3_v3(s->impact_size, size);
					}
					/*only fracture on new entries, this is necessary because after loading a file
					 *the pointcache thinks it is empty and a fracture is attempted ! */
					if (check_shard_size(fmd1, mi->id))
					{
						FractureID* fid1 = MEM_mallocN(sizeof(FractureID), "contact_callback_fractureid1");
						fid1->shardID = mi->id;
						BLI_addtail(&fmd1->shared->fracture_ids, fid1);
						fmd1->update_dynamic = true;
					}
				}
			}
		}
	}
	else if (obA) {
		ob1 = obA;
	}

	if (linear_index2 > -1 && linear_index2 < rbw->shared->numbodies)
	{
		//ob_index2 = rbw->cache_offset_map[linear_index2];
		//ob2 = rbw->objects[ob_index2];
		fmd2 = (FractureModifierData*)modifiers_findByType(ob2, eModifierType_Fracture);

		if (fmd2 && fmd2->fracture_mode == MOD_FRACTURE_DYNAMIC)
		{
			if (fmd2->shared->current_shard_entry && fmd2->shared->current_shard_entry->is_new)
			{
				RigidBodyOb *rbo = rbw->shared->cache_index_map[linear_index2];
				int id = rbo->meshisland_index;
				MeshIsland* mi = findMeshIsland(fmd2, id);
				Shard *s = findShard(fmd2, mi->id);

				if (mi->rigidbody->mass > 0) {
					force = force / mi->rigidbody->mass;
				}

				//printf("FORCE2:%f\n",force);
				bool canbreak = (force > fmd2->dynamic_force) || ((fmd2->limit_impact || obA) && can_break(ob1, ob2, fmd2->limit_impact));

				if (canbreak && check_constraints(fmd2, mi, rbw))
				{
					if (s) {
						float size[3] = {1.0f, 1.0f, 1.0f};

						if (ob1 == ob2 || (ob1 && ob1->rigidbody_object && ob1->rigidbody_object->type == RBO_TYPE_PASSIVE)) {
							size[0] = size[1] = size[2] = -1.0f;
						}
						else if (ob1) {
							BKE_object_dimensions_get(ob1, size);
						}
						copy_v3_v3(s->impact_loc, cp->contact_pos_world_onB);
						copy_v3_v3(s->impact_size, size);
					}

					if (check_shard_size(fmd2, mi->id))
					{
						FractureID* fid2 = MEM_mallocN(sizeof(FractureID), "contact_callback_fractureid2");
						fid2->shardID = mi->id;
						BLI_addtail(&fmd2->shared->fracture_ids, fid2);
						fmd2->update_dynamic = true;
					}
				}
			}
		}
	}

	cp = NULL;
}

void BKE_rigidbody_contact_callback(rbContactPoint* cp, void* world)
{
	Scene* scene = DEG_get_original_id(world);
	RigidBodyWorld *rbw = scene->rigidbody_world;
	check_fracture(cp, rbw, NULL, NULL);
}

void BKE_rigidbody_id_callback(void *world, void* island, int* objectId, int* islandId)
{
	MeshIsland *mi = (MeshIsland*)island;
	RigidBodyWorld *rbw = (RigidBodyWorld*)world;

	*objectId = -1;
	*islandId = -1;

	if (mi)
	{
		*objectId = rbw->shared->cache_offset_map[mi->linear_index];
		*islandId = mi->id;
	}
}

static void rigidbody_passive_fake_hook(MeshIsland *mi, float co[3])
{
	//no reshape necessary as vertcount didnt change, but update rbo->pos / orn ? according to change of 1st vertex
	//fake hook system
	if (mi->rigidbody->type == RBO_TYPE_PASSIVE &&
		mi->rigidbody->shared->physics_object && !(mi->rigidbody->flag & RBO_FLAG_KINEMATIC))
	{
		float oldloc[3], loc[3], diff[3], pos[3];
		oldloc[0] = mi->vertco[0];
		oldloc[1] = mi->vertco[1];
		oldloc[2] = mi->vertco[2];

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

void BKE_rigidbody_shard_validate(RigidBodyWorld *rbw, MeshIsland *mi, Object *ob, FractureModifierData* fmd,
                                  int rebuild, int transfer_speed, float size[3])
{
	if (mi == NULL || mi->rigidbody == NULL) {
		return;
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

static void activateCluster(MeshIsland *mi, int particle_index, RigidBodyWorld *rbw, Object *ob) {
	RigidBodyShardCon *con;
	int i = 0;
	for (i = 0; i < mi->participating_constraint_count; i++)
	{
		con = mi->participating_constraints[i];
		if (con->physics_constraint && con->mi1->particle_index == particle_index) /*|| particle_index == -1)*/
		{
			if (con->mi1->rigidbody->flag & RBO_FLAG_KINEMATIC) {
				BKE_rigidbody_activate(con->mi1->rigidbody, rbw, con->mi1, ob);
				//activateCluster(con->mi1, particle_index, rbw, ob);
			}
		}

		if (con->physics_constraint && con->mi2->particle_index == particle_index) /*|| particle_index == -1)*/
		{
			if (con->mi2->rigidbody->flag & RBO_FLAG_KINEMATIC) {
				BKE_rigidbody_activate(con->mi2->rigidbody, rbw, con->mi2, ob);
				//activateCluster(con->mi2, particle_index, rbw, ob);
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

static void handle_breaking_percentage(FractureModifierData* fmd, Object *ob, MeshIsland *mi, RigidBodyWorld *rbw, int breaking_percentage)
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
				if (con->mi1->particle_index != con->mi2->particle_index)
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
				if (con && con->mi1->particle_index != con->mi2->particle_index) {
					if (fmd->use_breaking)
					{
						if (con->physics_constraint) {

							RB_constraint_set_enabled(con->physics_constraint, false);
							/*if (con->mi1->rigidbody->flag & RBO_FLAG_KINEMATIC ||
								con->mi2->rigidbody->flag & RBO_FLAG_KINEMATIC ) */
							if (fmd->activate_broken)
							{
								BKE_rigidbody_activate(con->mi1->rigidbody, rbw, con->mi1, ob);
								activateCluster(con->mi1, con->mi1->particle_index, rbw, ob);

								BKE_rigidbody_activate(con->mi2->rigidbody, rbw, con->mi2, ob);
								activateCluster(con->mi2, con->mi2->particle_index, rbw, ob);
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
				if (con && fmd->use_breaking)
				{
					if (con->physics_constraint) {
						RB_constraint_set_enabled(con->physics_constraint, false);
						if (fmd->activate_broken)
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

static void test_deactivate_rigidbody(RigidBodyOb *rbo, MeshIsland* mi)
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

			if (mi != NULL)
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
	//rbsc->flag |= RBC_FLAG_DISABLE_COLLISIONS;
	BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, true);
	//set_constraint_index(fmd, rbsc);

	//thresh = RB_constraint_get_breaking_threshold(rbsc->physics_constraint);
	RB_constraint_set_breaking_threshold(rbsc->physics_constraint, thresh * weakening);

	RB_body_deactivate(rbsc->mi1->rigidbody->shared->physics_object);
	RB_body_deactivate(rbsc->mi2->rigidbody->shared->physics_object);
}

static void handle_deform_angle(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
								float anglediff, float weight, float deform_angle)
{
	if ((fmd->deform_angle > 0 || (fmd->deform_angle_weighted && weight > 0)) &&
		(anglediff > deform_angle))
	{
		/* if we have cluster breaking angle, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_deform_angle > 0 && rbsc->mi1->particle_index == rbsc->mi2->particle_index) ||
			 fmd->cluster_deform_angle == 0)
		{
			deform_constraint(fmd, ob, rbsc, rbw);
		}
	}

	if ((fmd->cluster_deform_angle > 0) && (rbsc->mi1->particle_index != rbsc->mi2->particle_index)
		&& anglediff > fmd->cluster_deform_angle)
	{
		deform_constraint(fmd, ob, rbsc, rbw);
	}
}

static void handle_deform_dist(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
								float distdiff, float weight, float deform_dist)
{
	if ((fmd->deform_distance > 0 || (fmd->deform_angle_weighted && weight > 0)) &&
		(distdiff > deform_dist))
	{
		/* if we have cluster breaking angle, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_deform_distance > 0 && rbsc->mi1->particle_index == rbsc->mi2->particle_index) ||
			 fmd->cluster_deform_distance == 0)
		{
			deform_constraint(fmd, ob, rbsc, rbw);
		}
	}

	if ((fmd->cluster_deform_distance > 0) && (rbsc->mi1->particle_index != rbsc->mi2->particle_index)
		&& distdiff > fmd->cluster_deform_distance)
	{
		deform_constraint(fmd, ob, rbsc, rbw);
	}
}


static void handle_breaking_angle(FractureModifierData *fmd, Object *ob, RigidBodyShardCon *rbsc, RigidBodyWorld *rbw,
								  float anglediff, float weight, float breaking_angle)
{
	if ((fmd->breaking_angle > 0 || (fmd->breaking_angle_weighted && weight > 0)) &&
		(anglediff > breaking_angle))
	{
		/* if we have cluster breaking angle, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_breaking_angle > 0 && rbsc->mi1->particle_index == rbsc->mi2->particle_index &&
			 rbsc->mi1->particle_index != -1) || fmd->cluster_breaking_angle == 0)
		{
			if (fmd->use_breaking)
			{
				//break constraint
				if (rbsc->physics_constraint) {
					RB_constraint_set_enabled(rbsc->physics_constraint, false);
					if (fmd->activate_broken)
					{
						BKE_rigidbody_activate(rbsc->mi1->rigidbody, rbw, rbsc->mi1, ob);
						BKE_rigidbody_activate(rbsc->mi2->rigidbody, rbw, rbsc->mi2, ob);
					}
				}
			}
		}
	}

	if ((fmd->cluster_breaking_angle > 0) && (rbsc->mi1->particle_index != rbsc->mi2->particle_index)
		&& anglediff > fmd->cluster_breaking_angle)
	{
		if (fmd->use_breaking)
		{
			if (rbsc->physics_constraint) {
				RB_constraint_set_enabled(rbsc->physics_constraint, false);
				if (fmd->activate_broken)
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
	if ((fmd->breaking_distance > 0 || (fmd->breaking_distance_weighted && weight > 0)) &&
		(distdiff > breaking_distance))
	{
		/* if we have cluster breaking distance, then only treat equal cluster indexes like the default, else all */
		if ((fmd->cluster_breaking_distance > 0 && rbsc->mi1->particle_index == rbsc->mi2->particle_index &&
			 rbsc->mi1->particle_index != -1) || fmd->cluster_breaking_distance == 0)
		{
			if (fmd->use_breaking)
			{
				if (rbsc->physics_constraint) {
					RB_constraint_set_enabled(rbsc->physics_constraint, false);
					if (fmd->activate_broken)
					{
						BKE_rigidbody_activate(rbsc->mi1->rigidbody, rbw, rbsc->mi1, ob);
						BKE_rigidbody_activate(rbsc->mi2->rigidbody, rbw, rbsc->mi2, ob);
					}
				}
			}
		}
	}

	if ((fmd->cluster_breaking_distance > 0) && (rbsc->mi1->particle_index != rbsc->mi2->particle_index)
		&& distdiff > fmd->cluster_breaking_distance)
	{
		if (fmd->use_breaking)
		{
			if (rbsc->physics_constraint) {
				RB_constraint_set_enabled(rbsc->physics_constraint, false);
				if (fmd->activate_broken)
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
	bool exceededAngle = false, exceededDist = false, regularBroken = false;

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

	if (exceededDist || exceededAngle) //|| regularBroken)
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
	if ((exceededDist || exceededAngle) /*&& !regularBroken*/)
	{
		if (rbsc->type == RBC_TYPE_6DOF_SPRING && rbsc->flag & RBC_FLAG_PLASTIC_ACTIVE)
		{
			if (rbsc->physics_constraint)
			{
				rigidbody_set_springs_active(rbsc, false);
				RB_constraint_set_enabled(rbsc->physics_constraint, rbsc->flag & RBC_FLAG_ENABLED);
				//rbsc->flag &= ~RBC_FLAG_PLASTIC_ACTIVE;
			}
		}
	}
}

static void handle_regular_breaking(FractureModifierData *fmd, Object *ob, RigidBodyWorld *rbw, RigidBodyShardCon *rbsc, float max_con_mass, bool rebuild)
{
	float weight = MIN2(rbsc->mi1->thresh_weight, rbsc->mi2->thresh_weight);
	float breaking_angle = fmd->breaking_angle_weighted ? fmd->breaking_angle * weight : fmd->breaking_angle;
	float breaking_distance = fmd->breaking_distance_weighted ? fmd->breaking_distance * weight : fmd->breaking_distance;
	float deform_angle = fmd->deform_angle_weighted ? fmd->deform_angle * weight : fmd->deform_angle;
	float deform_distance = fmd->deform_distance_weighted ? fmd->deform_distance * weight : fmd->deform_distance;
	float dist, angle, distdiff, anglediff;

	if ((fmd->use_mass_dependent_thresholds || fmd->use_compounds /*|| fmd->mass_threshold_factor > 0.0f*/)) {
		BKE_rigidbody_calc_threshold(max_con_mass, fmd, rbsc);
	}

	if ((((fmd->breaking_angle) > 0) || (fmd->breaking_angle_weighted && weight > 0) ||
		(fmd->breaking_distance > 0) || (fmd->breaking_distance_weighted && weight > 0) ||
		 (fmd->cluster_breaking_angle > 0 || (fmd->cluster_breaking_distance > 0))) /*&& !rebuild*/)
	{
		calc_dist_angle(rbsc, &dist, &angle, false);
		anglediff = fabs(angle - rbsc->start_angle);
		distdiff = fabs(dist - rbsc->start_dist);

		/* handle breaking */
		handle_breaking_angle(fmd, ob, rbsc, rbw, anglediff, weight, breaking_angle);
		handle_breaking_distance(fmd, ob, rbsc, rbw, distdiff, weight, breaking_distance);
	}

	if ((((fmd->deform_angle) > 0) || (fmd->deform_angle_weighted && weight > 0) ||
		(fmd->deform_distance > 0) || (fmd->deform_distance_weighted && weight > 0) ||
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
			if ((rbsc->mi1->particle_index != -1) && (rbsc->mi1->particle_index == rbsc->mi2->particle_index)) {
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
	if ((iterations > 0) && (fmd->fracture_mode != MOD_FRACTURE_EXTERNAL)) {
		rbsc->flag |= RBC_FLAG_OVERRIDE_SOLVER_ITERATIONS;
		rbsc->num_solver_iterations = iterations;
	}
}

bool BKE_rigidbody_modifier_update(Scene* scene, Object* ob, RigidBodyWorld *rbw, bool rebuild, Depsgraph *depsgraph)
{
	MeshIsland *mi;
	RigidBodyShardCon *rbsc;
	short laststeps = rbw->steps_per_second;
	float lastscale = rbw->time_scale;
	int i = 0;
	FractureModifierData *fmd = NULL;

	fmd = (FractureModifierData*) modifiers_findByType(DEG_get_original_object(ob), eModifierType_Fracture);

	if (BKE_rigidbody_modifier_active(fmd)) {
		float max_con_mass = 0;
		bool is_empty = BLI_listbase_is_empty(&fmd->shared->meshIslands);
		int count = 0, brokencount = 0, plastic = 0;
		float frame = 0;
		float size[3] = {1.0f, 1.0f, 1.0f};
		float bbsize[3];
		float locbb[3];

		//hacky check for ob->derivedFinal validity
		/*if (ob->derivedFinal && ob->derivedFinal->getNumLoopTri(ob->derivedFinal) > 0)
		{
			DM_mesh_boundbox(ob->derivedFinal, locbb, bbsize);
		}
		else
		{
			BKE_mesh_boundbox_calc((Mesh*)ob->data, locbb, bbsize);
		}*/


		if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
		{
			int fr = (int)BKE_scene_frame_get(scene);
			if (BKE_fracture_dynamic_lookup_mesh_state(fmd, fr, true, scene))
			{
				BKE_rigidbody_update_ob_array(rbw,
				                              rbw->shared->pointcache->flag & PTCACHE_BAKED);
			}
		}
		//else
		{
			if (rebuild || is_zero_m4(fmd->passive_parent_mat))
			{
				copy_m4_m4(fmd->passive_parent_mat, ob->obmat);
				//print_m4("Passivemat: \n", fmd->passive_parent_mat);
			}

			//print_m4("Obmat: \n", ob->obmat);
			//print_m4("Passivemat: \n", fmd->passive_parent_mat);
			//printf("WHERE IS CALC\n");
		}

		BKE_object_where_is_calc(depsgraph, scene, ob);
		fmd->constraint_island_count = 1;

		if ((ob->rigidbody_object && (ob->rigidbody_object->flag & RBO_FLAG_KINEMATIC) //&&
			 /*fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED*/)) {

			if (fmd->use_animated_mesh && fmd->anim_mesh_ob)
			{
				BKE_fracture_animated_loc_rot(fmd, ob, false, depsgraph);
			}
		}

		for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
			if (mi->rigidbody == NULL) {
				continue;
			}
			else {  /* as usual, but for each shard now, and no constraints*/
				/* perform simulation data updates as tagged */
				/* refresh object... */
				int do_rebuild = rebuild;

				BKE_mesh_boundbox_calc(mi->physics_mesh, locbb, bbsize);

				if ((rbw->flag & RBW_FLAG_REBUILD_CONSTRAINTS) && fmd->fracture_mode != MOD_FRACTURE_DYNAMIC)
				{
					//reset speeds
					//printf("ZEROIZING speed (shard)\n");
					//zero_v3(mi->rigidbody->lin_vel);
					//zero_v3(mi->rigidbody->ang_vel);
					mi->rigidbody->flag |= RBO_FLAG_NEEDS_VALIDATE;
					mi->rigidbody->flag &= ~RBO_FLAG_PROPAGATE_TRIGGER;
				}

				BKE_rigidbody_passive_hook(fmd, mi, ob, scene, depsgraph);

				if (fmd->use_breaking && fmd->fracture_mode != MOD_FRACTURE_EXTERNAL)
				{
					float weight = mi->thresh_weight;
					int breaking_percentage = fmd->breaking_percentage_weighted ? (fmd->breaking_percentage * weight) :
																				  fmd->breaking_percentage;

					if (fmd->breaking_percentage > 0 || (fmd->breaking_percentage_weighted && weight > 0) ||
						(fmd->cluster_breaking_percentage > 0))
					{
						handle_breaking_percentage(fmd, ob, mi, rbw, breaking_percentage);
					}
				}

				BKE_rigidbody_shard_validate(rbw, is_empty ? NULL : mi, ob, fmd, do_rebuild,
											 fmd->fracture_mode == MOD_FRACTURE_DYNAMIC, bbsize);

				mi->constraint_index = mi->id;

			}

			/* update simulation object... */
			if (fmd->fracture_mode == MOD_FRACTURE_EXTERNAL)
			{
				Shard *s = BLI_findlink(&fmd->shared->frac_mesh->shard_map, mi->id);
				if (s)
					copy_v3_v3(size, s->impact_size);
			}

			BKE_rigidbody_update_sim_ob(scene, rbw, ob, mi->rigidbody, mi->centroid, mi, size, fmd, depsgraph);
		}

		if (fmd->use_mass_dependent_thresholds) {
			max_con_mass = BKE_rigidbody_calc_max_con_mass(ob);
		}

		frame = BKE_scene_frame_get(scene);

		for (rbsc = fmd->shared->meshConstraints.first; rbsc; rbsc = rbsc->next) {

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
				if (!fmd->is_dynamic_external) {
					RB_constraint_set_enabled(rbsc->physics_constraint, rbsc->flag & RBC_FLAG_ENABLED);
					rbsc->flag |= RBC_FLAG_NEEDS_VALIDATE;
				}

				if (((fmd->fracture_mode == MOD_FRACTURE_EXTERNAL) || ((fmd->fracture_mode == MOD_FRACTURE_DYNAMIC) && fmd->is_dynamic_external))
					&& (rbsc->type == RBC_TYPE_6DOF_SPRING))
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

			if (rebuild) { // || rbsc->mi1->rigidbody->flag & RBO_FLAG_KINEMATIC_REBUILD ||
				//rbsc->mi2->rigidbody->flag & RBO_FLAG_KINEMATIC_REBUILD) {
				/* World has been rebuilt so rebuild constraint */
				BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, true);
				BKE_rigidbody_start_dist_angle(rbsc, fmd->fracture_mode == MOD_FRACTURE_EXTERNAL ||
											   (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && fmd->is_dynamic_external), true);
				//TODO ensure evaluation on transform change too
			}

			else if (rbsc->flag & RBC_FLAG_NEEDS_VALIDATE || fmd->fracture_mode == MOD_FRACTURE_DYNAMIC) {
				BKE_rigidbody_validate_sim_shard_constraint(rbw, fmd, ob, rbsc, false);
				//if (fmd->fracture_mode == MOD_FRACTURE_EXTERNAL)
				//	BKE_rigidbody_start_dist_angle(rbsc, true);
			}

			if (fmd->fracture_mode != MOD_FRACTURE_EXTERNAL && !rebuild)
			{
				handle_regular_breaking(fmd, ob, rbw, rbsc, max_con_mass, rebuild);
			}

			if (((fmd->fracture_mode == MOD_FRACTURE_EXTERNAL) || (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && fmd->is_dynamic_external)) &&
				(rbsc->flag & RBC_FLAG_USE_BREAKING) && !rebuild)
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

		printf("Constraints: Frame %d , Total %d,  Intact %d,  Broken %d, Plastic %d\n", (int)frame, count, count-brokencount, brokencount, plastic);
		return true;
	}
	else
	{
		return false;
	}
}

void BKE_rigidbody_passive_hook(FractureModifierData *fmd, MeshIsland *mi, Object* ob, Scene* scene, Depsgraph *depsgraph)
{
	RigidBodyOb *rbo = mi->rigidbody;

	if (rbo->type == RBO_TYPE_PASSIVE && !(rbo->flag & RBO_FLAG_KINEMATIC))
	{
		Mesh *dm = fmd->shared->visible_mesh_cached;
		ModifierData *md;
		bool found = false;

		if (dm)
		{
			MVert *mv = dm->mvert + mi->vertex_indices[0];
			int totvert = dm->totvert;
			float acc[3];
			copy_v3_v3(acc, mv->co);

			for (md = ob->modifiers.first; md; md = md->next)
			{
				if (md->type == eModifierType_Fracture)
				{
					if ((FractureModifierData*)md == fmd)
					{
						found = true;
					}
				}

				//only eval following hookmodifiers, based on our derivedmesh
				if (md->type == eModifierType_Hook && found)
				{
					float (*vertexCos)[3], old[3], diff[3];
					const ModifierTypeInfo *mti = modifierType_getInfo(md->type);
					HookModifierData *hmd = (HookModifierData*)md;
					ModifierEvalContext mctx = {.depsgraph = depsgraph, .object = ob};

					//skip hook modifiers which were just added and arent valid yet
					if (!hmd->object)
						continue;

					BKE_object_where_is_calc(depsgraph, scene, hmd->object);

					vertexCos = BKE_mesh_vertexCos_get(dm, &totvert);
					copy_v3_v3(old, vertexCos[mi->vertex_indices[0]]);

					mti->deformVerts(md, &mctx, dm, vertexCos, totvert);

					sub_v3_v3v3(diff, vertexCos[mi->vertex_indices[0]], old);
					add_v3_v3(acc, diff);
					MEM_freeN(vertexCos);
				}
			}

			rigidbody_passive_fake_hook(mi, acc);
		}
	}
}

bool BKE_rigidbody_modifier_sync(ModifierData *md, Object *ob, Scene *scene, float ctime)
{
	RigidBodyWorld *rbw = scene->rigidbody_world;
	bool modFound = false;
	FractureModifierData *fmd = NULL;
	MeshIsland *mi;
	bool exploOK = false;
	RigidBodyOb *rbo;
	float size[3] = {1, 1, 1};
	float centr[3];
	float imat[4][4];


	if (md->type == eModifierType_Fracture) {

		fmd = (FractureModifierData *)md;
		bool mode = fmd->fracture_mode == MOD_FRACTURE_EXTERNAL;

		exploOK = !fmd->valid_mesh || (fmd->valid_mesh && fmd->shared->frac_mesh && fmd->shared->dm) ||
		          mode || fmd->is_dynamic_external;

		if (BKE_rigidbody_modifier_active(fmd) && exploOK) {
			modFound = true;

			if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC && !(ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ))
			{
				int frame = (int)ctime;

				if (BKE_fracture_dynamic_lookup_mesh_state(fmd, frame, true, scene))
				{
					BKE_rigidbody_update_ob_array(rbw,
					                              rbw->shared->pointcache->flag & PTCACHE_BAKED);
				}
			}

			invert_m4_m4(imat, fmd->passive_parent_mat);
			invert_m4_m4(ob->imat, ob->obmat);

			for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {

				rbo = mi->rigidbody;
				if (!rbo || !ob->rigidbody_object) {
					continue;
				}

				if (ob->rigidbody_object->type == RBO_TYPE_ACTIVE) {
#if 0
					if (rbo->type == RBO_TYPE_PASSIVE)
					{
						printf("PASSIVE: %s\n", mi->name);
						if (rbo->flag & RBO_FLAG_KINEMATIC) {
							printf("KINEMATIC: %s\n", mi->name);
						}
					}
#endif
					BKE_rigidbody_passive_fake_parenting(fmd, ob, rbo, imat);
				}

				/* use rigid body transform after cache start frame if objects is not being transformed */
				if (BKE_rigidbody_check_sim_running(rbw, ctime) && !(ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ)) {

					/* keep original transform when the simulation is muted */
					if (rbw->flag & RBW_FLAG_MUTED)
						break;
				}
				/* otherwise set rigid body transform to current obmat*/
				else {

					mat4_to_loc_quat(rbo->pos, rbo->orn, ob->obmat);
					mat4_to_size(size, ob->obmat);
					copy_v3_v3(centr, mi->centroid);
					mul_v3_v3(centr, size);
					mul_qt_v3(rbo->orn, centr);
					add_v3_v3(rbo->pos, centr);

					if (mode)
					{
						mul_qt_qtqt(rbo->orn, rbo->orn, mi->rot);
					}

					zero_v3(rbo->lin_vel);
					zero_v3(rbo->ang_vel);
				}

				if ((ob->rigidbody_object->type == RBO_TYPE_ACTIVE) && (rbo->type == RBO_TYPE_ACTIVE || rbo->flag & RBO_FLAG_KINEMATIC)) {

					float quat[4];

					if (mode)
					{
						float iquat[4];
						invert_qt_qt(iquat, mi->rot);
						mul_qt_qtqt(quat, rbo->orn, iquat);
					}
					else {
						copy_qt_qt(quat, rbo->orn);
					}

					BKE_rigidbody_update_cell(mi, ob, rbo->pos, quat, fmd, (int)ctime, rbw);
				}
			}

			copy_m4_m4(fmd->passive_parent_mat, ob->obmat);
			//print_m4("Passivemat2: \n", fmd->passive_parent_mat);

			if ((ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ) /* || (mode && rbw) */||
				((ob->rigidbody_object) && (ob->rigidbody_object->flag & RBO_FLAG_KINEMATIC)))
			{
				/* update "original" matrix */
				copy_m4_m4(fmd->origmat, ob->obmat);
				if (ob->flag & SELECT && G.moving & G_TRANSFORM_OBJ && rbw) {
					RigidBodyShardCon *con;

					rbw->flag |= RBW_FLAG_OBJECT_CHANGED;
					BKE_rigidbody_cache_reset(rbw);
					/* re-enable all constraints as well */
					for (con = fmd->shared->meshConstraints.first; con; con = con->next) {
						//con->flag |= RBC_FLAG_ENABLED;
						//con->flag |= RBC_FLAG_NEEDS_VALIDATE;
						if (con->physics_constraint)
							RB_constraint_set_enabled(con->physics_constraint, con->flag & RBC_FLAG_ENABLED);
					}
				}
			}

			if (!is_zero_m4(fmd->origmat) && rbw && !(rbw->flag & RBW_FLAG_OBJECT_CHANGED))
			{
				if (fmd->fracture_mode != MOD_FRACTURE_EXTERNAL /*&& !fmd->is_dynamic_external*/)
				{
					copy_m4_m4(ob->obmat, fmd->origmat);
					zero_m4(fmd->origmat);
				}
			}

			return modFound;
		}
	}

	return modFound;
}

void BKE_rigidbody_passive_fake_parenting(FractureModifierData *fmd, Object *ob, RigidBodyOb *rbo, float imat[4][4])
{
	if (rbo->type == RBO_TYPE_PASSIVE && rbo->shared->physics_object)
	{
		//fake parenting, move all passive rbos together with original object in FM case
		float quat[4];
		//float imat[4][4];

		//first get rid of old obmat (=passive_parent_mat) -> do outside loop, expensive function due to profiler
		//invert_m4_m4(imat, fmd->passive_parent_mat);
		mat4_to_quat(quat, imat);
		mul_m4_v3(imat, rbo->pos);
		mul_qt_qtqt(rbo->orn, quat, rbo->orn);

		//then apply new obmat
		mat4_to_quat(quat, ob->obmat);
		mul_m4_v3(ob->obmat, rbo->pos);
		mul_qt_qtqt(rbo->orn, quat, rbo->orn);

		RB_body_set_kinematic_state(rbo->shared->physics_object, true);
		RB_body_set_loc_rot(rbo->shared->physics_object, rbo->pos, rbo->orn);
	}
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
			triggered = go->ob->rigidbody_object->flag & RBO_FLAG_USE_KINEMATIC_DEACTIVATION;

			FractureModifierData *fmd = (FractureModifierData*)modifiers_findByType(go->ob, eModifierType_Fracture);
			if (fmd && triggered)
			{
				MeshIsland* mi;
				for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next)
				{
					if (mi->rigidbody)
					{
						mi->rigidbody->flag &= ~RBO_FLAG_KINEMATIC_REBUILD;
						if (kinematic)
						{
							if(!fmd->use_animated_mesh ||
							  (fmd->use_animated_mesh && (override_bind || mi->rigidbody->flag & RBO_FLAG_KINEMATIC_BOUND)))
							{
								if (fmd->use_animated_mesh && override_bind)
								{
									mi->rigidbody->flag |= RBO_FLAG_KINEMATIC_BOUND;
								}

								mi->rigidbody->flag |= RBO_FLAG_KINEMATIC;
							}
						}
						else
						{
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

#if 0
static ThreadMutex reset_lock = BLI_MUTEX_INITIALIZER;
static void resetDynamic(RigidBodyWorld *rbw, bool do_reset_always, bool override_bind, Depsgraph *depsgraph)
{
	CollectionObject *go;
	if (!rbw->group)
		return;

	for (go = rbw->group->gobject.first; go; go = go->next)
	{
		Object *ob = go->ob;
		FractureModifierData *fmd = (FractureModifierData*)modifiers_findByType(ob, eModifierType_Fracture);

		if ((fmd && fmd->fracture_mode != MOD_FRACTURE_DYNAMIC && !(rbw->pointcache->flag & PTCACHE_BAKED)))
		{
			//also purge distortion cache here too (called from cache reset
			fmd->distortion_cached = false;
			fmd->refresh_autohide = true;
		}

		if (fmd && fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
		{
			//Scene *scene = fmd->modifier.scene;
			MeshIsland *mi;

			if (do_reset_always && (!fmd->use_animated_mesh || (fmd->use_animated_mesh && override_bind)))
			{
				ModifierData *md;
				Mesh *dm = NULL;
				bool found = false;

				if (ob->type == OB_MESH)
				{
					dm = ob->data;
				}

				if (!dm)
				{
					return;
				}

				BLI_mutex_lock(&reset_lock);
				fmd->refresh = true;
				fmd->reset_shards = true;
				fmd->last_frame = INT_MAX;

				//sigh, apply all modifiers before fracture
				for (md = ob->modifiers.first; md; md = md->next)
				{
					const ModifierTypeInfo *mti = modifierType_getInfo(md->type);

					if (!found)
					{
						if (mti->deformVerts && (md->mode & (eModifierMode_Realtime | eModifierMode_Render)))
						{
							float (*vertexCos)[3];
							int totvert;

							BKE_mesh_vertexCos_get(dm, &totvert);
							mti->deformVerts(md, ob, dm, vertexCos, totvert, 0);
							CDDM_apply_vert_coords(dm, vertexCos);
							MEM_freeN(vertexCos);
						}

						if (mti->applyModifier && (md->mode & (eModifierMode_Realtime | eModifierMode_Render)))
						{
							DerivedMesh *ndm;

							if (md == (ModifierData*)fmd) {
								BLI_mutex_unlock(&reset_lock);
							}

							ndm = mti->applyModifier(md, ob, dm, 0);
							if (ndm != dm)
							{
								dm->needsFree = 1;
								dm->release(dm);
							}
							dm = ndm;
						}
					}

					if (md == (ModifierData*)fmd)
					{
						found = true;
						break;
					}
				}
				//BLI_mutex_unlock(&reset_lock);

				//DAG_id_tag_update(go->ob, OB_RECALC_ALL);
				//WM_main_add_notifier(NC_OBJECT | ND_MODIFIER, go->ob);
				//WM_main_add_notifier(NC_OBJECT | ND_TRANSFORM, go->ob);
			}

			for (mi = fmd->meshIslands.first; mi; mi = mi->next)
			{
				zero_v3(mi->rigidbody->lin_vel);
				zero_v3(mi->rigidbody->ang_vel);
			}
		}
	}
}

#endif
#endif
