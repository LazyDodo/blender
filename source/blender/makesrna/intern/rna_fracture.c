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
 * Copyright (C) 2017 by Martin Felke.
 * All rights reserved.
 *
 * The Original Code is: all of this file
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/makesrna/intern/rna_fracture.c
 *  \ingroup RNA
 */

#include <stdlib.h>

#include "DNA_mesh_types.h"
#include "DNA_meta_types.h"
#include "DNA_fracture_types.h"
#include "DNA_modifier_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_collection_types.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BKE_rigidbody.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "rna_internal.h"

#include "WM_api.h"
#include "WM_types.h"

#ifdef RNA_RUNTIME

#include "BKE_fracture.h"
#include "BKE_modifier.h"
#include "DEG_depsgraph.h"

#ifdef WITH_BULLET
#  include "RBI_api.h"
#endif

#define FM_FLAG_SET(dest, value, flag) { \
	if (value) \
		dest |= flag; \
	else \
		dest &= ~flag; \
}

static char *rna_Modifier_path(PointerRNA *ptr)
{
	ModifierData *md = ptr->data;
	char name_esc[sizeof(md->name) * 2];

	BLI_strescape(name_esc, md->name, sizeof(name_esc));
	return BLI_sprintfN("modifiers[\"%s\"]", name_esc);
}

static StructRNA *rna_Modifier_refine(struct PointerRNA *ptr)
{
	ModifierData *md = (ModifierData *)ptr->data;

	switch ((ModifierType)md->type) {
		case eModifierType_Subsurf:
			return &RNA_SubsurfModifier;
		case eModifierType_Lattice:
			return &RNA_LatticeModifier;
		case eModifierType_Curve:
			return &RNA_CurveModifier;
		case eModifierType_Build:
			return &RNA_BuildModifier;
		case eModifierType_Mirror:
			return &RNA_MirrorModifier;
		case eModifierType_Decimate:
			return &RNA_DecimateModifier;
		case eModifierType_Wave:
			return &RNA_WaveModifier;
		case eModifierType_Armature:
			return &RNA_ArmatureModifier;
		case eModifierType_Hook:
			return &RNA_HookModifier;
		case eModifierType_Softbody:
			return &RNA_SoftBodyModifier;
		case eModifierType_Boolean:
			return &RNA_BooleanModifier;
		case eModifierType_Array:
			return &RNA_ArrayModifier;
		case eModifierType_EdgeSplit:
			return &RNA_EdgeSplitModifier;
		case eModifierType_Displace:
			return &RNA_DisplaceModifier;
		case eModifierType_UVProject:
			return &RNA_UVProjectModifier;
		case eModifierType_Smooth:
			return &RNA_SmoothModifier;
		case eModifierType_Cast:
			return &RNA_CastModifier;
		case eModifierType_MeshDeform:
			return &RNA_MeshDeformModifier;
		case eModifierType_ParticleSystem:
			return &RNA_ParticleSystemModifier;
		case eModifierType_ParticleInstance:
			return &RNA_ParticleInstanceModifier;
		case eModifierType_Explode:
			return &RNA_ExplodeModifier;
		case eModifierType_Cloth:
			return &RNA_ClothModifier;
		case eModifierType_Collision:
			return &RNA_CollisionModifier;
		case eModifierType_Bevel:
			return &RNA_BevelModifier;
		case eModifierType_Shrinkwrap:
			return &RNA_ShrinkwrapModifier;
		case eModifierType_Fluidsim:
			return &RNA_FluidSimulationModifier;
		case eModifierType_Mask:
			return &RNA_MaskModifier;
		case eModifierType_SimpleDeform:
			return &RNA_SimpleDeformModifier;
		case eModifierType_Multires:
			return &RNA_MultiresModifier;
		case eModifierType_Surface:
			return &RNA_SurfaceModifier;
		case eModifierType_Smoke:
			return &RNA_SmokeModifier;
		case eModifierType_Solidify:
			return &RNA_SolidifyModifier;
		case eModifierType_Screw:
			return &RNA_ScrewModifier;
		case eModifierType_Ocean:
			return &RNA_OceanModifier;
		case eModifierType_Warp:
			return &RNA_WarpModifier;
		case eModifierType_WeightVGEdit:
			return &RNA_VertexWeightEditModifier;
		case eModifierType_WeightVGMix:
			return &RNA_VertexWeightMixModifier;
		case eModifierType_WeightVGProximity:
			return &RNA_VertexWeightProximityModifier;
		case eModifierType_DynamicPaint:
			return &RNA_DynamicPaintModifier;
		case eModifierType_Remesh:
			return &RNA_RemeshModifier;
		case eModifierType_Skin:
			return &RNA_SkinModifier;
		case eModifierType_LaplacianSmooth:
			return &RNA_LaplacianSmoothModifier;
		case eModifierType_Triangulate:
			return &RNA_TriangulateModifier;
		case eModifierType_UVWarp:
			return &RNA_UVWarpModifier;
		case eModifierType_MeshCache:
			return &RNA_MeshCacheModifier;
		case eModifierType_LaplacianDeform:
			return &RNA_LaplacianDeformModifier;
		case eModifierType_Wireframe:
			return &RNA_WireframeModifier;
		case eModifierType_DataTransfer:
			return &RNA_DataTransferModifier;
		case eModifierType_NormalEdit:
			return &RNA_NormalEditModifier;
		case eModifierType_CorrectiveSmooth:
			return &RNA_CorrectiveSmoothModifier;
		case eModifierType_MeshSequenceCache:
			return &RNA_MeshSequenceCacheModifier;
		case eModifierType_SurfaceDeform:
			return &RNA_SurfaceDeformModifier;
		case eModifierType_Fracture:
			return &RNA_FractureModifier;
		case eModifierType_WeightedNormal:
			return &RNA_WeightedNormalModifier;
		/* Default */
		case eModifierType_None:
		case eModifierType_ShapeKey:
		case NUM_MODIFIER_TYPES:
			return &RNA_Modifier;
	}

	return &RNA_Modifier;
}

/* Vertex Groups */

#define RNA_MOD_VGROUP_NAME_SET(_type, _prop)                                               \
static void rna_##_type##Modifier_##_prop##_set(PointerRNA *ptr, const char *value)         \
{                                                                                           \
	_type##ModifierData *tmd = (_type##ModifierData *)ptr->data;                            \
	rna_object_vgroup_name_set(ptr, value, tmd->_prop, sizeof(tmd->_prop));                 \
}

RNA_MOD_VGROUP_NAME_SET(Fracture, thresh_defgrp_name);
RNA_MOD_VGROUP_NAME_SET(Fracture, passive_defgrp_name);
RNA_MOD_VGROUP_NAME_SET(Fracture, inner_defgrp_name);

#undef RNA_MOD_VGROUP_NAME_SET

/* UV layers */

#define RNA_MOD_UVLAYER_NAME_SET(_type, _prop)                                              \
static void rna_##_type##Modifier_##_prop##_set(PointerRNA *ptr, const char *value)         \
{                                                                                           \
	_type##ModifierData *tmd = (_type##ModifierData *)ptr->data;                            \
	rna_object_uvlayer_name_set(ptr, value, tmd->_prop, sizeof(tmd->_prop));                \
}

RNA_MOD_UVLAYER_NAME_SET(Fracture, uvlayer_name);

#undef RNA_MOD_UVLAYER_NAME_SET


static void rna_FractureModifier_threshold_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->breaking_threshold = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_contact_dist_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->contact_dist = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_use_constraints_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData *)ptr->data;
	FM_FLAG_SET(rmd->flag, value, MOD_FRACTURE_USE_CONSTRAINTS);
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_use_constraint_collision_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData *)ptr->data;
	FM_FLAG_SET(rmd->flag, value, MOD_FRACTURE_USE_CONSTRAINT_COLLISION);
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}


static void rna_FractureModifier_mass_dependent_thresholds_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData *)ptr->data;
	FM_FLAG_SET(rmd->flag, value, MOD_FRACTURE_USE_MASS_DEP_THRESHOLDS);
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_constraint_limit_set(PointerRNA *ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->constraint_limit = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_breaking_percentage_set(PointerRNA *ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->breaking_percentage = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_breaking_angle_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->breaking_angle = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_breaking_distance_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->breaking_distance = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_deform_angle_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->deform_angle = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_deform_distance_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->deform_distance = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_deform_angle_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_deform_angle = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_deform_distance_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_deform_distance = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_threshold_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_breaking_threshold = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_deform_weakening_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->deform_weakening = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_solver_iterations_override_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->solver_iterations_override = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_solver_iterations_override_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_solver_iterations_override = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_autohide_dist_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->autohide_dist = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_AUTOHIDE);
	rmd->distortion_cached = false;
}

static void rna_FractureModifier_automerge_dist_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->automerge_dist = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_AUTOHIDE);
	rmd->distortion_cached = false;
}

static void rna_FractureModifier_cluster_breaking_angle_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_breaking_angle = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_breaking_distance_set(PointerRNA *ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_breaking_distance = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_breaking_percentage_set(PointerRNA *ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_breaking_percentage = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_use_breaking_set(PointerRNA *ptr, bool value)
{
	RigidBodyShardCon* rbsc;
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	FM_FLAG_SET(rmd->flag, value, MOD_FRACTURE_USE_BREAKING);

	for (rbsc = rmd->shared->constraints.first; rbsc; rbsc = rbsc->next)
	{
		if (value == true){
			rbsc->flag |= RBC_FLAG_USE_BREAKING;
		}
		else {
			rbsc->flag &= ~RBC_FLAG_USE_BREAKING;
		}

		rbsc->flag |= RBC_FLAG_NEEDS_VALIDATE;
	}
}

static void rna_FractureModifier_constraint_type_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->constraint_type = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_constraint_type_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cluster_constraint_type = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_constraint_target_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->constraint_target = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_frac_algorithm_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->frac_algorithm = value;
}

static void rna_FractureModifier_point_source_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->point_source = value;
}

static void rna_FractureModifier_point_seed_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->point_seed = value;
}

static void rna_FractureModifier_percentage_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->percentage = value;
}


static void rna_FractureModifier_extra_group_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->extra_group = value.data;
}

static void rna_FractureModifier_inner_material_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->inner_material = value.data;
}

static void rna_FractureModifier_splinter_length_set(PointerRNA* ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->splinter_length = value;
}

static void rna_FractureModifier_splinter_axis_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->splinter_axis = value;
}

static void rna_FractureModifier_normal_search_radius_set(PointerRNA* ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->normal_search_radius = value;
}


static void rna_FractureModifier_fractal_cuts_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->fractal_cuts = value;
}

static void rna_FractureModifier_fractal_amount_set(PointerRNA* ptr, float value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->fractal_amount = value;
}


static void rna_FractureModifier_fractal_iterations_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->fractal_iterations = value;
}

static void rna_FractureModifier_cutter_group_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->cutter_group = value.data;
}

static void rna_FractureModifier_pack_group_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->pack_group = value.data;
}


static void rna_FractureModifier_autohide_filter_group_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->autohide_filter_group = value.data;
}

static void rna_FractureModifier_do_merge_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData *)ptr->data;
	rmd->perform_merge = value;
}

static void rna_FractureModifier_cluster_count_set(PointerRNA* ptr, int value)
{
	FractureModifierData *rmd = (FractureModifierData *)ptr->data;
	rmd->cluster_count = value;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_cluster_group_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData *)ptr->data;
	rmd->cluster_group = value.data;
	FM_FLAG_SET(rmd->shared->flag, true, MOD_FRACTURE_REFRESH_CONSTRAINTS);
}

static void rna_FractureModifier_anim_mesh_ob_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->anim_mesh_ob = value.data;
}

static void rna_FractureModifier_dupli_ob_set(PointerRNA* ptr, PointerRNA value)
{
	FractureModifierData *rmd = (FractureModifierData*)ptr->data;
	rmd->dupli_ob = value.data;
	rmd->shared->flag |= MOD_FRACTURE_REFRESH;
}


static void rna_Modifier_update(Main *UNUSED(bmain), Scene *scene, PointerRNA *ptr)
{
	FractureModifierData *fmd = (FractureModifierData*)ptr->data;
	bool dupli = (fmd->flag & MOD_FRACTURE_USE_DUPLI) && fmd->dupli_ob;

	BKE_rigidbody_cache_reset(scene);

	if (dupli) {
		fmd->shared->flag |= MOD_FRACTURE_REFRESH;
	}

	DEG_id_tag_update(ptr->id.data, OB_RECALC_DATA | OB_RECALC_OB | OB_RECALC_TIME |
									DEG_TAG_COPY_ON_WRITE | DEG_TAG_BASE_FLAGS_UPDATE);
	DEG_id_tag_update(&scene->id, DEG_TAG_COPY_ON_WRITE | DEG_TAG_BASE_FLAGS_UPDATE);

	WM_main_add_notifier(NC_OBJECT | ND_MODIFIER, ptr->id.data);
	WM_main_add_notifier(NC_OBJECT | ND_POINTCACHE, NULL);
	WM_main_add_notifier(NC_OBJECT | ND_TRANSFORM, NULL);
	WM_main_add_notifier(NC_OBJECT | ND_PARENT, NULL);
	WM_main_add_notifier(NC_SCENE | ND_FRAME, scene);
	WM_main_add_notifier(NC_WINDOW, NULL);
	WM_main_add_notifier(NC_MATERIAL | ND_SHADING, NULL);
}

#endif

void RNA_def_fracture(BlenderRNA *brna)
{
	StructRNA *srna, *subrna;
	PropertyRNA *prop;

	int noteflag = 0;

	static EnumPropertyItem prop_fracture_algorithm[] = {
		{MOD_FRACTURE_BISECT_FAST, "BISECT_FAST", 0, "Fast Bisect", "Use a faster but more inaccurate bisection algorithm, also creates uglier shards."},
		{MOD_FRACTURE_BISECT_FAST_FILL, "BISECT_FAST_FILL", 0, "Fast Bisect + Fill ", "Use the faster but different bisection algorithm and fill cut faces"},
		{MOD_FRACTURE_BOOLEAN, "BOOLEAN", 0, "Voronoi + Boolean", "Use voronoi and boolean intersection as fracture algorithm"},
		{MOD_FRACTURE_BISECT_FILL, "BISECT_FILL", 0, "Voronoi + Bisect + Fill", "Use voronoi and mesh bisect as fracture algorithm, fill cut faces"},
		{MOD_FRACTURE_BISECT, "BISECT", 0, "Voronoi + Bisect", "Use voronoi and mesh bisect as fracture algorithm, don't fill cut faces"},
		{MOD_FRACTURE_BOOLEAN_FRACTAL, "BOOLEAN_FRACTAL", 0, "Voronoi + Fractal Boolean", "Use voronoi and boolean intersection with fractally subdivided cells" },
		{0, NULL, 0, NULL, NULL}
	};

	static EnumPropertyItem prop_point_source_items[] = {
		{MOD_FRACTURE_OWN_PARTICLES, "OWN_PARTICLES", 0, "Own Particles", "Use own particles as point cloud"},
		{MOD_FRACTURE_OWN_VERTS, "OWN_VERTS", 0, "Own Vertices", "Use own vertices as point cloud"},
		{MOD_FRACTURE_EXTRA_PARTICLES, "EXTRA_PARTICLES", 0, "Extra Particles", "Use particles of group objects as point cloud"},
		{MOD_FRACTURE_EXTRA_VERTS, "EXTRA_VERTS", 0, "Extra Vertices", "Use vertices of group objects as point cloud"},
		{MOD_FRACTURE_CUSTOM, "CUSTOM", 0, "Custom", "Use custom cutter object group"},
		//{MOD_FRACTURE_GREASEPENCIL, "GREASE_PENCIL", 0, "Grease Pencil", "Use grease pencil points as point cloud"},
		{MOD_FRACTURE_UNIFORM, "UNIFORM", 0, "Uniform", "Use a random uniform pointcloud generated over the bounding box"},
		{MOD_FRACTURE_GRID, "GRID", 0, "Grid", "Use a grid pointcloud generated over the bounding box"},
		{0, NULL, 0, NULL, NULL}
	};

	static EnumPropertyItem prop_splinter_axises[] = {
		{MOD_FRACTURE_SPLINTER_X, "SPLINTER_X", 0, "Splinter X", "Splinters in X Direction"},
		{MOD_FRACTURE_SPLINTER_Y, "SPLINTER_Y", 0, "Splinter Y", "Splinters in Y Direction"},
		{MOD_FRACTURE_SPLINTER_Z, "SPLINTER_Z", 0, "Splinter Z", "Splinters in Z Direction"},
		{0, NULL, 0, NULL, NULL}
	};

	static EnumPropertyItem prop_constraint_targets[] = {
		{MOD_FRACTURE_CENTROID, "CENTROID", 0, "Centroid", "Build constraints based on distances between centroids"},
		{MOD_FRACTURE_VERTEX, "VERTEX", 0, "Vertex", "Build constraints based on distances between vertices (use lower values here)"},
		{0, NULL, 0, NULL, NULL}
	};

	static EnumPropertyItem prop_keep_cutter_shards[] = {
		{MOD_FRACTURE_KEEP_BOTH, "KEEP_BOTH", 0, "Both", "Keep both shards(intersect and difference)"},
		{MOD_FRACTURE_KEEP_INTERSECT, "KEEP_INTERSECT", 0, "Intersect Only", "Keep intersected shards only"},
		{MOD_FRACTURE_KEEP_DIFFERENCE, "KEEP_DIFFERENCE", 0, "Difference Only", "Keep difference shards only"},
		{0, NULL, 0, NULL, NULL}
	};

	static EnumPropertyItem prop_dynamic_constraints[] = {
		{MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS, "NO_CONSTRAINTS", 0, "None", "Build no new constraints"},
		{MOD_FRACTURE_MIXED_DYNAMIC_CONSTRAINTS, "MIXED_CONSTRAINTS", 0, "Mixed", "Build constraints between new and old shards"},
		{MOD_FRACTURE_ALL_DYNAMIC_CONSTRAINTS, "ALL_CONSTRAINTS", 0, "All", "Build all new constraints"},
		{0, NULL, 0, NULL, NULL}
	};

	srna = RNA_def_struct(brna, "FractureModifier", "Modifier");
	RNA_def_struct_ui_text(srna, "Fracture Modifier", "Add a fracture container to this object");
	RNA_def_struct_sdna(srna, "FractureModifierData");
	RNA_def_struct_ui_icon(srna, ICON_MOD_EXPLODE);

	prop = RNA_def_property(srna, "cluster_count", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "cluster_count");
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_cluster_count_set", NULL);
	RNA_def_property_range(prop, 0, 100000);
	RNA_def_property_ui_text(prop, "Cluster Count", "Amount of clusters built from existing shards, 0 for none");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	//simulation stuff...
	prop = RNA_def_property(srna, "breaking_threshold", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "breaking_threshold");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_threshold_set", NULL);
	RNA_def_property_ui_text(prop, "Inner Breaking threshold", "Threshold to break constraints between shards in the same object");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.0001f, 5);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_constraints", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_CONSTRAINTS);
	RNA_def_property_boolean_funcs(prop, NULL, "rna_FractureModifier_use_constraints_set");
	RNA_def_property_ui_text(prop, "Use Constraints", "Create constraints between all shards");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "contact_dist", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "contact_dist");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_contact_dist_set", NULL);
	RNA_def_property_ui_text(prop, "Search Radius", "Limit search radius up to which two mesh shards are being connected, 0 for entire boundingbox");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_mass_dependent_thresholds", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_MASS_DEP_THRESHOLDS);
	RNA_def_property_boolean_funcs(prop, NULL, "rna_FractureModifier_mass_dependent_thresholds_set");
	RNA_def_property_ui_text(prop, "Use Mass Dependent Thresholds", "Match the breaking threshold according to the masses of the constrained shards");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "constraint_limit", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "constraint_limit");
	RNA_def_property_range(prop, 0, INT_MAX);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_constraint_limit_set", NULL);
	RNA_def_property_ui_text(prop, "Search Limit", "Maximum number of surrounding shards being taken into account per shard during constraint creation, 0 for unlimited");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "breaking_percentage", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "breaking_percentage");
	RNA_def_property_range(prop, 0, 100);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_breaking_percentage_set", NULL);
	RNA_def_property_ui_text(prop, "Breaking Percentage", "Percentage of broken constraints per island which leads to breaking of all others");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "breaking_angle", PROP_FLOAT, PROP_ANGLE);
	RNA_def_property_float_sdna(prop, NULL, "breaking_angle");
	RNA_def_property_range(prop, 0, DEG2RADF(360.0));
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_breaking_angle_set", NULL);
	RNA_def_property_ui_text(prop, "Breaking Angle", "Angle in degrees above which constraint should break");
	RNA_def_property_ui_range(prop, 0.0f, DEG2RADF(360.0), 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "breaking_distance", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "breaking_distance");
	RNA_def_property_range(prop, 0, FLT_MAX);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_breaking_distance_set", NULL);
	RNA_def_property_ui_text(prop, "Breaking Distance", "Distance above which constraint should break");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_breaking_threshold", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "cluster_breaking_threshold");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_cluster_threshold_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Breaking threshold", "Threshold to break constraints INSIDE a cluster of shards");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.01f, 5);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "solver_iterations_override", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "solver_iterations_override");
	RNA_def_property_range(prop, 0, INT_MAX);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_solver_iterations_override_set", NULL);
	RNA_def_property_ui_text(prop, "Solver Iterations Override", "Override the world constraint solver iteration value with this value, 0 means no override");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "frac_algorithm", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, prop_fracture_algorithm);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FractureModifier_frac_algorithm_set", NULL);
	RNA_def_property_ui_text(prop, "Fracture Algorithm", "Select type of fracture algorithm");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "shard_count", PROP_INT, PROP_NONE);
	RNA_def_property_range(prop, 1, 100000);
	RNA_def_property_int_default(prop, 10);
	RNA_def_property_ui_text(prop, "Shard Count", "How many sub-shards should be generated from the current shard");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "point_source", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, prop_point_source_items);
	RNA_def_property_flag(prop, PROP_ENUM_FLAG);
	RNA_def_property_enum_default(prop, MOD_FRACTURE_UNIFORM);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FractureModifier_point_source_set", NULL);
	RNA_def_property_ui_text(prop, "Point Source", "Source of point cloud");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "point_seed", PROP_INT, PROP_NONE);
	RNA_def_property_range(prop, 0, 100000);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_point_seed_set", NULL);
	RNA_def_property_ui_text(prop, "Seed", "Seed for uniform pointcloud");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "percentage", PROP_INT, PROP_NONE);
	RNA_def_property_range(prop, 0, 100);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_percentage_set", NULL);
	RNA_def_property_ui_text(prop, "Percentage", "Percentage of the sum of points of all selected pointsources to actually use for fracture");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "extra_group", PROP_POINTER, PROP_NONE);
    RNA_def_property_struct_type(prop, "Collection");
	RNA_def_property_ui_text(prop, "Extra Group",
	                         "Depending on whether extra particles or extra vertices is chosen, particles or vertices of objects from this group serve as point source. Both at the same time is not possible.");
    RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_extra_group_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "split_islands", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_SPLIT_TO_ISLANDS);
	RNA_def_property_ui_text(prop, "Split to Islands", "Split each piece to separate islands");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "thresh_vertex_group", PROP_STRING, PROP_NONE);
	RNA_def_property_string_sdna(prop, NULL, "thresh_defgrp_name");
	RNA_def_property_ui_text(prop, "Threshold Vertex Group", "Vertex group name for defining weighted thresholds on different mesh parts");
	RNA_def_property_string_funcs(prop, NULL, NULL, "rna_FractureModifier_thresh_defgrp_name_set");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "passive_vertex_group", PROP_STRING, PROP_NONE);
	RNA_def_property_string_sdna(prop, NULL, "passive_defgrp_name");
	RNA_def_property_ui_text(prop, "Passive Vertex Group", "Vertex group name for defining passive mesh parts (will remain static during rigidbody simulation");
	RNA_def_property_string_funcs(prop, NULL, NULL, "rna_FractureModifier_passive_defgrp_name_set");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "fix_normals", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_FIX_NORMALS);
	RNA_def_property_ui_text(prop, "Fix normals (WIP)", "Fix normals of fractured smooth objects, to let cracks nearly disappear");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "inner_material", PROP_POINTER, PROP_NONE);
	RNA_def_property_ui_text(prop, "Inner Material", "Material which will be applied to all inner faces generated by fracturing");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_inner_material_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "inner_vertex_group", PROP_STRING, PROP_NONE);
	RNA_def_property_string_sdna(prop, NULL, "inner_defgrp_name");
	RNA_def_property_ui_text(prop, "Inner Vertex Group", "Vertex group name for defining inner vertices (will contain vertices of inner faces (Boolean, Bisect + Fill only) ");
	RNA_def_property_string_funcs(prop, NULL, NULL, "rna_FractureModifier_inner_defgrp_name_set");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "auto_execute", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_AUTOEXECUTE);
	RNA_def_property_ui_text(prop, "Auto Execute", "Automatic execution of fracturing, CAUTION: this can be slow and buggy");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "pack_group", PROP_POINTER, PROP_NONE);
	RNA_def_property_struct_type(prop, "Collection");
	RNA_def_property_ui_text(prop, "Sub Object Group", "The meshes of each member object will be aggregated into a combined mesh. This will override the base mesh of the fractured object using this group.");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_pack_group_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "autohide_dist", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "autohide_dist");
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_autohide_dist_set", NULL);
	RNA_def_property_range(prop, 0.0f, 100.0f);
	RNA_def_property_ui_text(prop, "Autohide Distance", "Distance between faces below which both faces should be hidden");
	RNA_def_property_ui_range(prop, 0.0f, 100.0f, 0.01f, 5);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "automerge_dist", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "automerge_dist");
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_automerge_dist_set", NULL);
	RNA_def_property_range(prop, 0.0f, 100.0f);
	RNA_def_property_ui_range(prop, 0.0f, 100.0f, 0.01f, 5);
	RNA_def_property_ui_text(prop, "Automerge Distance",
 "Distance between faces below which vertices of both faces should be merged; (costs performance, use with smooth objects and fix normals to better hide cracks)");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "breaking_percentage_weighted", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_BREAKING_PERCENTAGE_WEIGHTED);
	RNA_def_property_ui_text(prop, "Weighted Percentage", "Modify breaking percentage by threshold weights");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "breaking_angle_weighted", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_BREAKING_ANGLE_WEIGHTED);
	RNA_def_property_ui_text(prop, "Weighted Angle", "Modify breaking angle by threshold weights");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "breaking_distance_weighted", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_BREAKING_DISTANCE_WEIGHTED);
	RNA_def_property_ui_text(prop, "Weighted Distance", "Modify breaking distance by threshold weights");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_particle_birth_coordinates", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_PARTICLE_BIRTH_COORDS);
	RNA_def_property_ui_text(prop, "Use Particle Birth Coordinates", "Use birth or simulated state particle coordinates");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "splinter_axis", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, prop_splinter_axises);
	RNA_def_property_flag(prop, PROP_ENUM_FLAG);
	RNA_def_property_enum_default(prop, MOD_FRACTURE_SPLINTER_Z);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FractureModifier_splinter_axis_set", NULL);
	RNA_def_property_ui_text(prop, "Splinter Axis", "Global direction of splinters");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "splinter_length", PROP_FLOAT, PROP_NONE);
	RNA_def_property_range(prop, 1.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "Splinter length", "Length of splinters");
	RNA_def_property_ui_range(prop, 1.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_splinter_length_set", NULL);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_solver_iterations_override", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "cluster_solver_iterations_override");
	RNA_def_property_range(prop, 0, INT_MAX);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_cluster_solver_iterations_override_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Solver Iterations Override", "Override the world constraint solver iteration value for INSIDE clusters with this value, 0 means no override");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "normal_search_radius", PROP_FLOAT, PROP_NONE);
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "Normal Search Radius", "Radius in which to search for valid normals");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_normal_search_radius_set", NULL);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_breaking_percentage", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "cluster_breaking_percentage");
	RNA_def_property_range(prop, 0, 100);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_cluster_breaking_percentage_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Breaking Percentage", "Percentage of broken constraints per cluster which leads to breaking of all others");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_breaking_angle", PROP_FLOAT, PROP_ANGLE);
	RNA_def_property_float_sdna(prop, NULL, "cluster_breaking_angle");
	RNA_def_property_range(prop, 0, DEG2RADF(360.0));
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_cluster_breaking_angle_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Breaking Angle", "Angle in degrees above which constraint between different clusters should break");
	RNA_def_property_ui_range(prop, 0.0f, DEG2RADF(360.0), 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_breaking_distance", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "cluster_breaking_distance");
	RNA_def_property_range(prop, 0, FLT_MAX);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_cluster_breaking_distance_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Breaking Distance", "Distance above which constraint between different clusters should break");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	/*Breakable constraints*/
	prop = RNA_def_property(srna, "use_breaking", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_BREAKING);
	RNA_def_property_boolean_funcs(prop, NULL, "rna_FractureModifier_use_breaking_set");
	RNA_def_property_ui_text(prop, "Breakable",
	                         "Constraints can be broken if it receives an impulse above the threshold");
	RNA_def_property_update(prop, /*NC_OBJECT | ND_POINTCACHE*/ 0, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_smooth", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_SMOOTH);
	RNA_def_property_ui_text(prop, "Smooth Inner Faces", "Set Inner Faces to Smooth Shading (needs refracture)");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "fractal_cuts", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "fractal_cuts");
	RNA_def_property_range(prop, 1, 10);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_fractal_cuts_set", NULL);
	RNA_def_property_ui_text(prop, "Fractal Grid Cuts", "Number of fractal cuts on each cell");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "fractal_amount", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "fractal_amount");
	RNA_def_property_range(prop, 0, 20);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_fractal_amount_set", NULL);
	RNA_def_property_ui_text(prop, "Fractal Displacement", "Amount of fractal displacement on each cell");
	RNA_def_property_ui_range(prop, 0.0f, 2.0f, 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "fractal_iterations", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "fractal_iterations");
	RNA_def_property_range(prop, 1, 10);
	RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_fractal_iterations_set", NULL);
	RNA_def_property_ui_text(prop, "Fractal Iterations", "Number of times the number of cuts will be made to the grid, with the given fractal amount");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_group", PROP_POINTER, PROP_NONE);
	RNA_def_property_pointer_sdna(prop, NULL, "cluster_group");
    RNA_def_property_struct_type(prop, "Collection");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_cluster_group_set", NULL, NULL);
	RNA_def_property_ui_text(prop, "Cluster Group", "Centroids of objects in this group determine where cluster centers will be");
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cutter_group", PROP_POINTER, PROP_NONE);
    RNA_def_property_struct_type(prop, "Collection");
	RNA_def_property_ui_text(prop, "Cutter Group", "A set of objects to make boolean cuts against");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_cutter_group_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "constraint_type", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "constraint_type");
	RNA_def_property_enum_items(prop, rna_enum_rigidbody_constraint_type_items);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FractureModifier_constraint_type_set", NULL);
	RNA_def_property_ui_text(prop, "Constraint Type", "Type of Rigid Body Constraint between shards and inside clusters");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_constraint_type", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "cluster_constraint_type");
	RNA_def_property_enum_items(prop, rna_enum_rigidbody_constraint_type_items);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FractureModifier_cluster_constraint_type_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Constraint Type", "Type of Rigid Body Constraint between clusters");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "constraint_target", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, prop_constraint_targets);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FractureModifier_constraint_target_set", NULL);
	RNA_def_property_enum_default(prop, MOD_FRACTURE_CENTROID);
	RNA_def_property_ui_text(prop, "Search Method", "Method to search constraints among surrounding shards");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "dynamic_force", PROP_FLOAT, PROP_NONE);
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "Dynamic force threshold", "Only break dynamically when force is above this threshold");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 3);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "limit_impact", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_LIMIT_IMPACT);
	RNA_def_property_ui_text(prop, "Limit Impact", "Activates only shards within the impact object size approximately");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "uv_layer", PROP_STRING, PROP_NONE);
	RNA_def_property_string_sdna(prop, NULL, "uvlayer_name");
	RNA_def_property_ui_text(prop, "Inner UV Map", "Name of UV map for inner faces");
	RNA_def_property_string_funcs(prop, NULL, NULL, "rna_FractureModifier_uvlayer_name_set");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "autohide_filter_group", PROP_POINTER, PROP_NONE);
    RNA_def_property_struct_type(prop, "Collection");
	RNA_def_property_ui_text(prop, "Autohide Filter Group",
	                         "Shards within the radius of each group member object will be excluded from autohide, which selectively reveals cracks ");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_autohide_filter_group_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "keep_cutter_shards", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, prop_keep_cutter_shards);
	RNA_def_property_enum_default(prop, MOD_FRACTURE_KEEP_BOTH);
	RNA_def_property_ui_text(prop, "Cutter Mode", "Determines which shards to keep from cutting process");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "boolean_double_threshold", PROP_FLOAT, PROP_DISTANCE);
	RNA_def_property_float_sdna(prop, NULL, "boolean_double_threshold");
	RNA_def_property_range(prop, 0, 1.0f);
	RNA_def_property_ui_range(prop, 0, 1, 0.0001, 6);
	RNA_def_property_ui_text(prop, "Overlap Threshold",  "Threshold for checking overlapping geometry");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "dynamic_percentage", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "dynamic_percentage");
	RNA_def_property_range(prop, 0, 100);
	//RNA_def_property_int_funcs(prop, NULL, "rna_FractureModifier_breaking_percentage_set", NULL);
	RNA_def_property_ui_text(prop, "Constraint Percentage", "Percentage of broken constraints per island which allows dynamic fracturing of this island");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "dynamic_new_constraints", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_items(prop, prop_dynamic_constraints);
	RNA_def_property_enum_default(prop, MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS);
	RNA_def_property_ui_text(prop, "New Constraints", "Which constraints are created while dynamically fracturing");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "dynamic_min_size", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "dynamic_min_size");
	RNA_def_property_range(prop, 0.001f, 10000.0f);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Minimum Size",  "Minimum shard size in blenderunits");
	RNA_def_property_ui_range(prop, 0.001f, 10000.0f, 0.1f, 2);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_constraint_collision", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_CONSTRAINT_COLLISION);
	RNA_def_property_boolean_funcs(prop, NULL, "rna_FractureModifier_use_constraint_collision_set");
	RNA_def_property_ui_text(prop, "Constrained Collision", "Let constrained shards collide with each other");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "inner_crease", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "inner_crease");
	RNA_def_property_range(prop, 0.0f, 1.0f);
	RNA_def_property_float_default(prop, 0.0f);
	RNA_def_property_ui_text(prop, "Inner Crease",  "Crease at edges of inner faces");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1f, 2);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "material_offset_intersect", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "mat_ofs_intersect");
	RNA_def_property_range(prop, 0, SHRT_MAX);
	RNA_def_property_ui_text(prop, "Intersect Material Offset", "Offset material index of intersected shards");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "material_offset_difference", PROP_INT, PROP_NONE);
	RNA_def_property_int_sdna(prop, NULL, "mat_ofs_difference");
	RNA_def_property_range(prop, 0, SHRT_MAX);
	RNA_def_property_ui_text(prop, "Difference Material Offset", "Offset material index of difference shards");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "rectangular_alignment", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "rectangular_alignment");
	RNA_def_property_range(prop, 0.0f, 1.0f);
	RNA_def_property_float_default(prop, 0.0f);
	RNA_def_property_ui_text(prop, "Rectangular Alignment",
	                         "1 means only orthogonal cuts, move down to 0 to get more diagonal-ish cuts");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.1f, 2);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "keep_distort", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_KEEP_DISTORT);
	RNA_def_property_ui_text(prop, "Keep Distortion", "Whether or not to make the distortion on torn merged shards persistent.");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "perform_merge", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_funcs(prop, NULL, "rna_FractureModifier_do_merge_set");
	RNA_def_property_ui_text(prop, "Perform Merge", "Whether or not to actually weld the prepared automerge geometry.");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "deform_angle", PROP_FLOAT, PROP_ANGLE);
	RNA_def_property_float_sdna(prop, NULL, "deform_angle");
	RNA_def_property_range(prop, 0, DEG2RADF(360.0));
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_deform_angle_set", NULL);
	RNA_def_property_ui_text(prop, "Deforming Angle", "Angle in degrees above which constraint should keep its deform");
	RNA_def_property_ui_range(prop, 0.0f, DEG2RADF(360.0), 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "deform_distance", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "deform_distance");
	RNA_def_property_range(prop, 0, FLT_MAX);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_deform_distance_set", NULL);
	RNA_def_property_ui_text(prop, "Deforming Distance", "Distance above which constraint should keep its deform");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_deform_angle", PROP_FLOAT, PROP_ANGLE);
	RNA_def_property_float_sdna(prop, NULL, "cluster_deform_angle");
	RNA_def_property_range(prop, 0, DEG2RADF(360.0));
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_cluster_deform_angle_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Deforming Angle", "Angle in degrees above which constraint between different clusters should keep its deform");
	RNA_def_property_ui_range(prop, 0.0f, DEG2RADF(360.0), 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "cluster_deform_distance", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "cluster_deform_distance");
	RNA_def_property_range(prop, 0, FLT_MAX);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_cluster_deform_distance_set", NULL);
	RNA_def_property_ui_text(prop, "Cluster Deforming Distance", "Distance above which constraint between different clusters should keep its deform");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "deform_angle_weighted", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_DEFORM_ANGLE_WEIGHTED);
	RNA_def_property_ui_text(prop, "Weighted Deforming Angle", "Modify deform angle by threshold weights");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "deform_distance_weighted", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_DEFORM_DISTANCE_WEIGHTED);
	RNA_def_property_ui_text(prop, "Weighted Deforming Distance", "Modify deform distance by threshold weights");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "deform_weakening", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "deform_weakening");
	RNA_def_property_range(prop, 0.0f, 1.0f);
	RNA_def_property_float_funcs(prop, NULL, "rna_FractureModifier_deform_weakening_set", NULL);
	RNA_def_property_ui_text(prop, "Deform Weakening Factor",
	                         "Multiplies the breaking threshold in each iteration with 1.0 - factor in order to weaken it at deform, 0 to disable");
	RNA_def_property_ui_range(prop, 0.0f, 1.0f, 0.0001f, 6);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "grid_resolution", PROP_INT, PROP_XYZ);
	RNA_def_property_int_sdna(prop, NULL, "grid_resolution");
	RNA_def_property_range(prop, 1, 10000);
	RNA_def_property_array(prop, 3);
	RNA_def_property_int_default(prop, 10);
	RNA_def_property_ui_text(prop, "Grid Resolution", "How many grid cells per Bounding Box Axis");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_centroids", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_CENTROIDS);
	RNA_def_property_ui_text(prop, "Use Centroids", "Only output the meshisland centroids as vertices");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_vertices", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_VERTICES);
	RNA_def_property_ui_text(prop, "Use Vertices", "Only output the meshisland vertices");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_self_collision", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_SELF_COLLISION);
	RNA_def_property_ui_text(prop, "Self Collision", "Allow collisions between constraint islands");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_animated_mesh", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_ANIMATED_MESH);
	RNA_def_property_ui_text(prop, "Use Animated Mesh", "Allow moving original vertices to move the shards as well");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "animated_mesh_input", PROP_POINTER, PROP_NONE);
	RNA_def_property_ui_text(prop, "Animated Mesh", "Input for moving original vertices to move the shards as well");
	RNA_def_property_pointer_sdna(prop, NULL, "anim_mesh_ob");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_anim_mesh_ob_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_flag(prop, PROP_ID_SELF_CHECK);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_animated_mesh_rotation", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_ANIMATED_MESH_ROTATION);
	RNA_def_property_ui_text(prop, "Use Rotation", "Allow moving original vertices to rotate the shards as well");
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "animated_mesh_limit", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "anim_bind_limit");
	RNA_def_property_range(prop, 0, FLT_MAX);
	RNA_def_property_ui_text(prop, "Limit", "Maximal distance between shards and verts to perform a bind between");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 4);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "grid_offset", PROP_FLOAT, PROP_XYZ);
	RNA_def_property_float_sdna(prop, NULL, "grid_offset");
	RNA_def_property_range(prop, 0.0f, 1.0f);
	RNA_def_property_array(prop, 3);
	RNA_def_property_float_default(prop, 0.0f);
	RNA_def_property_ui_text(prop, "Grid Offset", "How far odd rows are relatively offset compared to even ones");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "grid_spacing", PROP_FLOAT, PROP_XYZ);
	RNA_def_property_float_sdna(prop, NULL, "grid_spacing");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_array(prop, 3);
	RNA_def_property_float_default(prop, 0.0f);
	RNA_def_property_ui_text(prop, "Grid Spacing", "How much space inbetween the bricks, in each direction");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_constraint_group", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_GROUP_CONSTRAINTS_ONLY);
	//RNA_def_property_boolean_funcs(prop, NULL, "rna_FractureModifier_use_constraint_group_set");
	RNA_def_property_ui_text(prop, "Constraints Only", "Only manage the external constraints in this Fracture Modifier");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "activate_broken", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_ACTIVATE_BROKEN);
	RNA_def_property_ui_text(prop, "Activate Broken", "Activate both shards or all elements of the cluster if a constraint breaks");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	subrna = RNA_def_struct(brna, "FractureModifierShared", NULL);
	RNA_def_struct_sdna(subrna, "FractureModifierData_Shared");

	prop = RNA_def_property(srna, "shared", PROP_POINTER, PROP_NONE);
	RNA_def_property_pointer_sdna(prop, NULL, "shared");
	RNA_def_property_struct_type(prop, "FractureModifierShared");

	prop = RNA_def_property(srna, "use_dynamic", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_DYNAMIC);
	RNA_def_property_ui_text(prop, "Dynamic", "Dynamically fracture shards during the sim");
	//RNA_def_property_update(prop, noteflag, "rna_Modifier_update"); //keep cache

	prop = RNA_def_property(srna, "dynamic_activation_size", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "dynamic_activation_size");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_float_default(prop, 1.0f);
	RNA_def_property_ui_text(prop, "Activation Size",  "Activation shard size in blenderunits");
	RNA_def_property_ui_range(prop, 0.001f, FLT_MAX, 0.1f, 2);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "contact_size", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "contact_size");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_float_default(prop, 0.0f);
	RNA_def_property_ui_text(prop, "Search Size", "Limit size of shards being connected above this value, 0 disables.");
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 0.1f, 2);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "dynamic_shard_count", PROP_INT, PROP_NONE);
	RNA_def_property_range(prop, 1, 100000);
	RNA_def_property_int_default(prop, 10);
	RNA_def_property_ui_text(prop, "Dynamic Shard Count",
	                         "How many sub-shards should be generated from the current shard dynamically");
	//RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	//RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "use_dupli", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "flag", MOD_FRACTURE_USE_DUPLI);
	RNA_def_property_ui_text(prop, "Use Duplis", "Take Dupli objects as shards");
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	prop = RNA_def_property(srna, "dupli_input", PROP_POINTER, PROP_NONE);
	RNA_def_property_ui_text(prop, "Dupli Object", "Input object for duplis to use");
	RNA_def_property_pointer_sdna(prop, NULL, "dupli_ob");
	RNA_def_property_pointer_funcs(prop, NULL, "rna_FractureModifier_dupli_ob_set", NULL, NULL);
	RNA_def_property_flag(prop, PROP_EDITABLE);
	RNA_def_property_flag(prop, PROP_ID_SELF_CHECK);
	RNA_def_property_clear_flag(prop, PROP_ANIMATABLE);
	RNA_def_property_update(prop, noteflag, "rna_Modifier_update");

	RNA_api_fracture(brna, subrna);
}
