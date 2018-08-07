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
 * Copyright (C) 2014 by Martin Felke.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/modifiers/intern/MOD_fracture.c
 *  \ingroup modifiers
 */

#include "MEM_guardedalloc.h"
#include "BKE_fracture.h"
#include "BKE_library_query.h"
#include "BKE_mesh.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "BLI_ghash.h"

#include "DNA_group_types.h"
#include "DNA_fracture_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

#include "MOD_util.h"

static void initData(ModifierData *md)
{
	FractureModifierData *fmd = (FractureModifierData *) md;

	fmd->fracture_all = false;
	fmd->cluster_count = 0;
	fmd->extra_group = NULL;
	fmd->frac_algorithm = MOD_FRACTURE_BOOLEAN;
	fmd->point_source = MOD_FRACTURE_UNIFORM;
	fmd->shard_id = -1;
	fmd->shard_count = 10;
	fmd->percentage = 100;

	fmd->visible_mesh = NULL;
	fmd->visible_mesh_cached = NULL;
	fmd->refresh = false;
    //zero_m4(fmd->origmat);
	fmd->breaking_threshold = 10.0f;
	fmd->use_constraints = false;
	fmd->contact_dist = 1.0f;
	fmd->use_mass_dependent_thresholds = false;
	fmd->explo_shared = false;
	fmd->constraint_limit = 50;
	fmd->breaking_distance = 0;
	fmd->breaking_angle = 0;
	fmd->breaking_percentage = 0;     /* disable by default*/
	fmd->max_vol = 0;
	fmd->refresh_constraints = false;

	fmd->cluster_breaking_threshold = 1000.0f;
	fmd->solver_iterations_override = 0;
	fmd->cluster_solver_iterations_override = 0;
	fmd->shards_to_islands = false;
	fmd->execute_threaded = false;
	fmd->nor_tree = NULL;
	fmd->fix_normals = false;
	fmd->auto_execute = false;
	fmd->face_pairs = NULL;
	fmd->autohide_dist = 0.0f;
	fmd->automerge_dist = 0.0f;

	fmd->breaking_percentage_weighted = false;
	fmd->breaking_angle_weighted = false;
	fmd->breaking_distance_weighted = false;

	/* XXX needed because of messy particle cache, shows incorrect positions when start/end on frame 1
	 * default use case is with this flag being enabled, disable at own risk */
	fmd->use_particle_birth_coordinates = true;
	fmd->splinter_length = 1.0f;
	fmd->nor_range = 1.0f;

	fmd->cluster_breaking_angle = 0;
	fmd->cluster_breaking_distance = 0;
	fmd->cluster_breaking_percentage = 0;

	/* used for advanced fracture settings now, XXX needs rename perhaps*/
	fmd->use_experimental = 0;
	fmd->use_breaking = true;
	fmd->use_smooth = false;

	fmd->fractal_cuts = 1;
	fmd->fractal_amount = 1.0f;
	fmd->physics_mesh_scale = 0.75f;
	fmd->fractal_iterations = 5;

	fmd->cluster_group = NULL;
	fmd->cutter_group = NULL;

	fmd->grease_decimate = 100.0f;
	fmd->grease_offset = 0.5f;
	fmd->use_greasepencil_edges = true;

	fmd->cutter_axis = MOD_FRACTURE_CUTTER_Z;
	fmd->cluster_constraint_type = RBC_TYPE_FIXED;
	fmd->constraint_type = RBC_TYPE_FIXED;
	fmd->vert_index_map = NULL;
	fmd->constraint_target = MOD_FRACTURE_CENTROID;
	fmd->vertex_island_map = NULL;

	fmd->meshIslands.first = NULL;
	fmd->meshIslands.last = NULL;
	fmd->meshConstraints.first = NULL;
	fmd->meshConstraints.last = NULL;

	fmd->fracture_mode = MOD_FRACTURE_PREFRACTURED;
    fmd->last_frame = 0; //INT_MIN
	fmd->dynamic_force = 10.0f;
	fmd->update_dynamic = false;
	fmd->limit_impact = false;
	fmd->reset_shards = false;

	fmd->use_compounds = false;
	fmd->impulse_dampening = 0.05f;
	fmd->directional_factor = 0.0f;
	fmd->minimum_impulse = 0.1f;
	fmd->mass_threshold_factor = 0.0f;

	fmd->autohide_filter_group = NULL;
	fmd->constraint_count = 0;

	fmd->boolean_double_threshold = 1e-6f;
	fmd->dynamic_percentage = 0.0f;
	fmd->dynamic_new_constraints = MOD_FRACTURE_ALL_DYNAMIC_CONSTRAINTS;
	fmd->dynamic_min_size = 1.0f;
	fmd->keep_cutter_shards = MOD_FRACTURE_KEEP_BOTH;
	fmd->use_constraint_collision = false;
	fmd->inner_crease = 0.0f;
	fmd->is_dynamic_external = false;

	fmd->mat_ofs_difference = 0;
	fmd->mat_ofs_intersect = 0;

	fmd->orthogonality_factor = 0.0f;
	fmd->keep_distort = false;
	fmd->do_merge = false;

	fmd->pack_storage.first = NULL;
	fmd->pack_storage.last = NULL;

	fmd->deform_weakening = 0.0f;
	fmd->distortion_cached = false;

	fmd->grid_resolution[0] = 10;
	fmd->grid_resolution[1] = 10;
	fmd->grid_resolution[2] = 10;

	fmd->use_centroids = false;
	fmd->use_vertices = false;
	fmd->use_self_collision = false;

	fmd->min_acceleration = 0.0f;
	fmd->max_acceleration = 1.0f;
	fmd->acceleration_fade = 0.85f;

	fmd->use_animated_mesh = false;
	fmd->anim_mesh_ob = NULL;
	fmd->anim_bind = NULL;
	fmd->anim_bind_len = 0;
	fmd->anim_mesh_rot = false;
	fmd->anim_bind_limit = 0.0f;
//	zero_v3(fmd->grid_offset);
//	zero_v3(fmd->grid_spacing);
	fmd->use_constraint_group = false;
	fmd->activate_broken = false;
}

static void freeData(ModifierData *md)
{
	FractureModifierData *fmd = (FractureModifierData *) md;

    if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
    {
        //BKE_fracture_dynamic_free(fmd, true, true);
    }
    else
    {
        //BKE_fracture_free(fmd, false, true);
    }
}


//XXX todo, simplify to copy generic stuff, maybe take shards over even, but re-init the meshisland verts as in packing system
static void copyData(ModifierData *md, ModifierData *target)
{
	FractureModifierData *rmd  = (FractureModifierData *)md;
	FractureModifierData *trmd = (FractureModifierData *)target;

	/*todo -> copy fracture stuff as well, and dont forget readfile / writefile...*/
    //zero_m4(trmd->origmat);

	/* vgroups  XXX TODO non ascii strings ?*/
/*	strncpy(trmd->thresh_defgrp_name, rmd->thresh_defgrp_name, strlen(rmd->thresh_defgrp_name));
	strncpy(trmd->ground_defgrp_name, rmd->ground_defgrp_name, strlen(rmd->ground_defgrp_name));
	strncpy(trmd->inner_defgrp_name, rmd->inner_defgrp_name, strlen(rmd->inner_defgrp_name));
    strncpy(trmd->acceleration_defgrp_name, rmd->acceleration_defgrp_name, strlen(rmd->acceleration_defgrp_name));*/

	trmd->visible_mesh = NULL;
	trmd->visible_mesh_cached = NULL;
	trmd->meshIslands.first = NULL;
	trmd->meshIslands.last = NULL;
	trmd->meshConstraints.first = NULL;
	trmd->meshConstraints.last = NULL;
	trmd->face_pairs = NULL;
	trmd->vert_index_map = NULL;
	trmd->vertex_island_map = NULL;

	trmd->breaking_threshold = rmd->breaking_threshold;
	trmd->use_constraints = rmd->use_constraints;
	trmd->contact_dist = rmd->contact_dist;
	trmd->use_mass_dependent_thresholds = rmd->use_mass_dependent_thresholds;
	trmd->explo_shared = rmd->explo_shared;

	trmd->refresh = true; //ensure valid data ?
	trmd->constraint_limit = rmd->constraint_limit;
	trmd->breaking_angle = rmd->breaking_angle;
	trmd->breaking_distance = rmd->breaking_distance;
	trmd->breaking_percentage = rmd->breaking_percentage;
	trmd->use_experimental = rmd->use_experimental;
	trmd->refresh_constraints = false;

	trmd->cluster_count = rmd->cluster_count;
	trmd->cluster_breaking_threshold = rmd->cluster_breaking_threshold;
	trmd->solver_iterations_override = rmd->solver_iterations_override;
	trmd->shards_to_islands = rmd->shards_to_islands;

	trmd->shard_count = rmd->shard_count;
	trmd->frac_algorithm = rmd->frac_algorithm;

	trmd->auto_execute = rmd->auto_execute;
	trmd->autohide_dist = rmd->autohide_dist;
	trmd->automerge_dist = rmd->automerge_dist;

	trmd->solver_iterations_override = rmd->solver_iterations_override;

	trmd->breaking_angle_weighted = rmd->breaking_angle_weighted;
	trmd->breaking_distance_weighted = rmd->breaking_distance_weighted;
	trmd->breaking_percentage_weighted = rmd->breaking_percentage_weighted;

	trmd->execute_threaded = rmd->execute_threaded;
	trmd->point_seed = rmd->point_seed;
	trmd->point_source = rmd->point_source;

	/*id refs ?*/
	trmd->inner_material = rmd->inner_material;
	trmd->extra_group = rmd->extra_group;

	/* sub object group  XXX Do we keep this ?*/
	trmd->dm_group = rmd->dm_group;

	trmd->cluster_group = rmd->cluster_group;
	trmd->cutter_group = rmd->cutter_group;

	trmd->use_particle_birth_coordinates = rmd->use_particle_birth_coordinates;
	trmd->splinter_length = rmd->splinter_length;
	trmd->splinter_axis = rmd->splinter_axis;
	trmd->cluster_solver_iterations_override = rmd->cluster_solver_iterations_override;

	trmd->cluster_breaking_angle = rmd->cluster_breaking_angle;
	trmd->cluster_breaking_distance = rmd->cluster_breaking_distance;
	trmd->cluster_breaking_percentage = rmd->cluster_breaking_percentage;

	trmd->use_breaking = rmd->use_breaking;
	trmd->use_smooth = rmd->use_smooth;
	trmd->fractal_cuts = rmd->fractal_cuts;
	trmd->fractal_amount = rmd->fractal_amount;

	trmd->grease_decimate = rmd->grease_decimate;
	trmd->grease_offset = rmd->grease_offset;
	trmd->use_greasepencil_edges = rmd->use_greasepencil_edges;
	trmd->cutter_axis = rmd->cutter_axis;

	trmd->cluster_constraint_type = rmd->cluster_constraint_type;
	trmd->constraint_target = rmd->constraint_target;

	trmd->fracture_mode = rmd->fracture_mode;
	trmd->last_frame = rmd->last_frame;
	trmd->dynamic_force = rmd->dynamic_force;
	trmd->update_dynamic = false;
	trmd->reset_shards = false;

	trmd->impulse_dampening = rmd->impulse_dampening;
	trmd->directional_factor = rmd->directional_factor;
	trmd->minimum_impulse = rmd->minimum_impulse;
	trmd->mass_threshold_factor = rmd->mass_threshold_factor;
	trmd->use_compounds = rmd->use_compounds;

	trmd->autohide_filter_group = rmd->autohide_filter_group;

	trmd->boolean_solver = rmd->boolean_solver;
	trmd->boolean_double_threshold = rmd->boolean_double_threshold;

	trmd->dynamic_percentage = rmd->dynamic_percentage;
	trmd->dynamic_new_constraints = rmd->dynamic_new_constraints;
	trmd->dynamic_min_size = rmd->dynamic_min_size;
	trmd->keep_cutter_shards = rmd->keep_cutter_shards;
	trmd->use_constraint_collision = rmd->use_constraint_collision;
	trmd->inner_crease = rmd->inner_crease;
	trmd->is_dynamic_external = rmd->is_dynamic_external;

	trmd->mat_ofs_difference = rmd->mat_ofs_difference;
	trmd->mat_ofs_intersect = rmd->mat_ofs_intersect;
	trmd->keep_distort = rmd->keep_distort;
	trmd->do_merge = rmd->do_merge;
	trmd->orthogonality_factor = rmd->orthogonality_factor;
	trmd->constraint_type = rmd->constraint_type;

	trmd->deform_angle = rmd->deform_angle;
	trmd->deform_angle_weighted = rmd->deform_angle_weighted;
	trmd->cluster_deform_angle = rmd->cluster_deform_angle;

	trmd->deform_distance = rmd->deform_distance;
	trmd->deform_distance_weighted = rmd->deform_distance_weighted;
	trmd->cluster_deform_distance = rmd->cluster_deform_distance;
	trmd->deform_weakening = rmd->deform_weakening;
	trmd->distortion_cached = rmd->distortion_cached;

//	copy_v3_v3_int(trmd->grid_resolution, rmd->grid_resolution);
	trmd->use_centroids = rmd->use_centroids;
	trmd->use_vertices = rmd->use_vertices;
	trmd->use_self_collision = rmd->use_self_collision;

	trmd->min_acceleration = rmd->min_acceleration;
	trmd->max_acceleration = rmd->max_acceleration;
	trmd->acceleration_fade = rmd->acceleration_fade;

	trmd->use_animated_mesh = rmd->use_animated_mesh;
	trmd->anim_mesh_ob = rmd->anim_mesh_ob;
	trmd->anim_bind_len = 0; //rmd->anim_bind_len;
	trmd->anim_bind = NULL;
	//memcpy(trmd->anim_bind, rmd->anim_bind, sizeof(int) * trmd->anim_bind_len);
	trmd->anim_mesh_rot = rmd->anim_mesh_rot;
	trmd->anim_bind_limit = rmd->anim_bind_limit;
//	copy_v3_v3(trmd->grid_offset, rmd->grid_offset);
//	copy_v3_v3(trmd->grid_spacing, rmd->grid_spacing);

	trmd->use_constraint_group = rmd->use_constraint_group;
	trmd->activate_broken = rmd->activate_broken;
}

static bool dependsOnTime(ModifierData *UNUSED(md))
{
	return true;
}

static bool dependsOnNormals(ModifierData *UNUSED(md))
{
	return true;
}

static void foreachIDLink(ModifierData *md, Object *ob,
                          IDWalkFunc walk, void *userData)
{
    FractureModifierData *fmd = (FractureModifierData *) md;

	walk(userData, ob, (ID **)&fmd->inner_material, IDWALK_CB_NOP);
	walk(userData, ob, (ID **)&fmd->extra_group, IDWALK_CB_NOP);
	walk(userData, ob, (ID **)&fmd->dm_group, IDWALK_CB_NOP);
	walk(userData, ob, (ID **)&fmd->cluster_group, IDWALK_CB_NOP);
	walk(userData, ob, (ID **)&fmd->cutter_group, IDWALK_CB_NOP);
	walk(userData, ob, (ID **)&fmd->autohide_filter_group, IDWALK_CB_NOP);
	walk(userData, ob, (ID **)&fmd->anim_mesh_ob, IDWALK_CB_NOP);
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *UNUSED(md))
{
	CustomDataMask dataMask = 0;
    dataMask |= CD_MASK_MDEFORMVERT | CD_MASK_MLOOPUV | CD_MASK_CREASE | CD_MASK_BWEIGHT | CD_MASK_MEDGE;
	return dataMask;
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
    FractureModifierData *fmd = (FractureModifierData *)md;
    if (fmd->anim_mesh_ob != NULL) {
        DEG_add_object_relation(ctx->node, fmd->anim_mesh_ob, DEG_OB_COMP_TRANSFORM, "Fracture Modifier Anim");
        DEG_add_object_relation(ctx->node, fmd->anim_mesh_ob, DEG_OB_COMP_GEOMETRY, "Fracture Modifier Anim");
    }

    if (fmd->extra_group) {
        CollectionObject *go;
        for (go = fmd->extra_group->gobject.first; go; go = go->next) {
            if (go->ob)
            {
                DEG_add_object_relation(ctx->node, go->ob, DEG_OB_COMP_TRANSFORM, "Fracture Modifier Extra");
                DEG_add_object_relation(ctx->node, go->ob, DEG_OB_COMP_GEOMETRY, "Fracture Modifier Extra");
            }
        }
    }

    if (fmd->autohide_filter_group) {
        CollectionObject *go;
        for (go = fmd->autohide_filter_group->gobject.first; go; go = go->next) {
            if (go->ob)
            {
                DEG_add_object_relation(ctx->node, go->ob, DEG_OB_COMP_TRANSFORM, "Fracture Modifier Autohide Filter");
                DEG_add_object_relation(ctx->node, go->ob, DEG_OB_COMP_GEOMETRY, "Fracture Modifier Autohide Filter");
            }
        }
    }

    if (fmd->cutter_group) {
        CollectionObject *go;
        for (go = fmd->cutter_group->gobject.first; go; go = go->next) {
            if (go->ob)
            {

                DEG_add_object_relation(ctx->node, go->ob, DEG_OB_COMP_TRANSFORM, "Fracture Modifier Cutter");
                DEG_add_object_relation(ctx->node, go->ob, DEG_OB_COMP_GEOMETRY, "Fracture Modifier Cutter");
            }
        }
    }


    /* We need own transformation as well. */
    DEG_add_object_relation(ctx->node, ctx->object, DEG_OB_COMP_TRANSFORM, "Fracture Modifier");
}

static void foreachObjectLink(
    ModifierData *md, Object *ob,
    ObjectWalkFunc walk, void *userData)
{
	FractureModifierData *fmd = (FractureModifierData *) md;

	if (fmd->anim_mesh_ob)
		walk(userData, ob, &fmd->anim_mesh_ob, IDWALK_CB_NOP);

	if (fmd->extra_group) {
        CollectionObject *go;
		for (go = fmd->extra_group->gobject.first; go; go = go->next) {
			if (go->ob) {
				walk(userData, ob, &go->ob, IDWALK_CB_NOP);
			}
		}
	}

	if (fmd->cutter_group) {
        CollectionObject *go;
		for (go = fmd->cutter_group->gobject.first; go; go = go->next) {
			if (go->ob) {
				walk(userData, ob, &go->ob, IDWALK_CB_NOP);
			}
		}
	}

	if (fmd->autohide_filter_group) {
        CollectionObject *go;
		for (go = fmd->autohide_filter_group->gobject.first; go; go = go->next) {
			if (go->ob) {
				walk(userData, ob, &go->ob, IDWALK_CB_NOP);
			}
		}
	}

	if (fmd->cutter_group) {
        CollectionObject *go;
		for (go = fmd->cutter_group->gobject.first; go; go = go->next) {
			if (go->ob)
			{
				walk(userData, ob, &go->ob, IDWALK_CB_NOP);
			}
		}
	}
}

static Mesh *applyModifier(ModifierData *md, const ModifierEvalContext *ctx, Mesh *derivedData)
{
	FractureModifierData *fmd = (FractureModifierData *) md;
    Mesh *pack_dm = BKE_fracture_mesh_from_packdata(fmd, derivedData);
    Mesh *final_dm = derivedData;
    Object* ob = ctx->object;
    Scene* scene = DEG_get_input_scene(ctx->depsgraph); //or evaluated ?

	if (fmd->fracture_mode == MOD_FRACTURE_PREFRACTURED)
	{
		bool init = false;

		//just track the frames for resetting automerge data when jumping
        int frame = (int)BKE_scene_frame_get(scene);

		//deactivate multiple settings for now, not working properly XXX TODO (also deactivated in RNA and python)
        final_dm = BKE_fracture_prefractured_apply(fmd, ob, pack_dm, ctx->depsgraph);

		if (init)
			fmd->shard_count = 10;

        fmd->last_frame = frame;
	}
	else if (fmd->fracture_mode == MOD_FRACTURE_DYNAMIC)
	{
		if (ob->rigidbody_object == NULL) {
			//initialize FM here once
			fmd->refresh = true;
		}

        final_dm = BKE_fracture_dynamic_apply(fmd, ob, pack_dm);
	}
	else if (fmd->fracture_mode == MOD_FRACTURE_EXTERNAL)
	{
        final_dm = BKE_fracture_external_apply(fmd, ob, pack_dm, derivedData, scene);
	}

	if (final_dm != derivedData)
	{
		//dont forget to create customdatalayers for crease and bevel weights (else they wont be drawn in editmode)
		final_dm->cd_flag |= (ME_CDFLAG_EDGE_CREASE | ME_CDFLAG_VERT_BWEIGHT | ME_CDFLAG_EDGE_BWEIGHT);
	}

	if (pack_dm != derivedData)
	{
        BKE_mesh_free(pack_dm);
		pack_dm = NULL;
	}

	return final_dm;
}

ModifierTypeInfo modifierType_Fracture = {
    /* name */              "Fracture",
    /* structName */        "FractureModifierData",
    /* structSize */        sizeof(FractureModifierData),
    /* type */              eModifierTypeType_Constructive,
    /* flags */             eModifierTypeFlag_AcceptsMesh |
                            eModifierTypeFlag_UsesPointCache |
                            eModifierTypeFlag_Single,

    /* copyData */          copyData,

    /* deformVerts_DM */    NULL,
    /* deformMatrices_DM */ NULL,
    /* deformVertsEM_DM */  NULL,
    /* deformMatricesEM_DM*/NULL,
    /* applyModifier_DM */  NULL,
    /* applyModifierEM_DM */NULL,

    /* deformVerts */       NULL,
    /* deformMatrices */    NULL,
    /* deformVertsEM */     NULL,
    /* deformMatricesEM */  NULL,
    /* applyModifier */     applyModifier,
    /* applyModifierEM */   NULL,

    /* initData */          initData,
    /* requiredDataMask */  requiredDataMask,
    /* freeData */          freeData,

    /* isDisabled */        NULL,
    /* updateDepsgraph */   updateDepsgraph,
    /* dependsOnTime */     dependsOnTime,
    /* dependsOnNormals */  dependsOnNormals,
    /* foreachObjectLink */ foreachObjectLink,
    /* foreachIDLink */     foreachIDLink,
    /* foreachTexLink */    NULL,
};

