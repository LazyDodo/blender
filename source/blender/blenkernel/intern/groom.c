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
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/gpencil.c
 *  \ingroup bke
 */

 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_string_utils.h"

#include "BLT_translation.h"

#include "DNA_groom_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_animsys.h"
#include "BKE_customdata.h"
#include "BKE_global.h"
#include "BKE_groom.h"
#include "BKE_hair.h"
#include "BKE_bvhutils.h"
#include "BKE_library.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_mesh_sample.h"
#include "BKE_object.h"
#include "BKE_object_facemap.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "bmesh.h"


void BKE_groom_init(Groom *groom)
{
	BLI_assert(MEMCMP_STRUCT_OFS_IS_ZERO(groom, id));
	
	groom->bb = BKE_boundbox_alloc_unit();
	
	groom->curve_res = 12;
	
	groom->hair_system = BKE_hair_new();
	groom->hair_draw_settings = BKE_hair_draw_settings_new();
}

void *BKE_groom_add(Main *bmain, const char *name)
{
	Groom *groom = BKE_libblock_alloc(bmain, ID_GM, name, 0);

	BKE_groom_init(groom);

	return groom;
}

void BKE_groom_bundle_curve_cache_clear(GroomBundle *bundle)
{
	if (bundle->curvecache)
	{
		MEM_freeN(bundle->curvecache);
		bundle->curvecache = NULL;
		bundle->curvesize = 0;
		bundle->totcurvecache = 0;
	}
}

static void groom_bundles_free(ListBase *bundles)
{
	for (GroomBundle *bundle = bundles->first; bundle; bundle = bundle->next)
	{
		BKE_groom_bundle_curve_cache_clear(bundle);
		
		if (bundle->sections)
		{
			MEM_freeN(bundle->sections);
		}
		if (bundle->verts)
		{
			MEM_freeN(bundle->verts);
		}
		if (bundle->scalp_region)
		{
			MEM_freeN(bundle->scalp_region);
		}
		if (bundle->guides)
		{
			MEM_freeN(bundle->guides);
		}
		if (bundle->guide_shape_weights)
		{
			MEM_freeN(bundle->guide_shape_weights);
		}
	}
	BLI_freelistN(bundles);
}

/** Free (or release) any data used by this groom (does not free the groom itself). */
void BKE_groom_free(Groom *groom)
{
	BKE_groom_batch_cache_free(groom);
	
	if (groom->editgroom)
	{
		EditGroom *edit = groom->editgroom;
		
		groom_bundles_free(&edit->bundles);
		
		MEM_freeN(edit);
		groom->editgroom = NULL;
	}
	
	MEM_SAFE_FREE(groom->bb);
	
	if (groom->hair_system)
	{
		BKE_hair_free(groom->hair_system);
	}
	if (groom->hair_draw_settings)
	{
		BKE_hair_draw_settings_free(groom->hair_draw_settings);
	}
	
	groom_bundles_free(&groom->bundles);
	
	MEM_SAFE_FREE(groom->mat);
	
	BKE_animdata_free(&groom->id, false);
}

/**
 * Only copy internal data of Groom ID from source to already allocated/initialized destination.
 * You probably never want to use that directly, use id_copy or BKE_id_copy_ex for typical needs.
 *
 * WARNING! This function will not handle ID user count!
 *
 * \param flag  Copying options (see BKE_library.h's LIB_ID_COPY_... flags for more).
 */
void BKE_groom_copy_data(Main *UNUSED(bmain), Groom *groom_dst, const Groom *groom_src, const int UNUSED(flag))
{
	groom_dst->bb = MEM_dupallocN(groom_src->bb);
	
	BLI_duplicatelist(&groom_dst->bundles, &groom_src->bundles);
	for (GroomBundle *bundle = groom_dst->bundles.first; bundle; bundle = bundle->next)
	{
		if (bundle->curvecache)
		{
			bundle->curvecache = MEM_dupallocN(bundle->curvecache);
		}
		if (bundle->sections)
		{
			bundle->sections = MEM_dupallocN(bundle->sections);
		}
		if (bundle->verts)
		{
			bundle->verts = MEM_dupallocN(bundle->verts);
		}
		if (bundle->scalp_region)
		{
			bundle->scalp_region = MEM_dupallocN(bundle->scalp_region);
		}
		if (bundle->guides)
		{
			bundle->guides = MEM_dupallocN(bundle->guides);
		}
		if (bundle->guide_shape_weights)
		{
			bundle->guide_shape_weights = MEM_dupallocN(bundle->guide_shape_weights);
		}
	}
	
	groom_dst->editgroom = NULL;
	
	if (groom_dst->hair_system)
	{
		groom_dst->hair_system = BKE_hair_copy(groom_dst->hair_system);
	}
	if (groom_dst->hair_draw_settings)
	{
		groom_dst->hair_draw_settings = BKE_hair_draw_settings_copy(groom_dst->hair_draw_settings);
	}
	
	groom_dst->mat = MEM_dupallocN(groom_src->mat);
}

Groom *BKE_groom_copy(Main *bmain, const Groom *groom)
{
	Groom *groom_copy;
	BKE_id_copy_ex(bmain, &groom->id, (ID **)&groom_copy, 0, false);
	return groom_copy;
}

void BKE_groom_make_local(Main *bmain, Groom *groom, const bool lib_local)
{
	BKE_id_make_local_generic(bmain, &groom->id, true, lib_local);
}


bool BKE_groom_minmax(Groom *groom, float min[3], float max[3])
{
	bool result = false;
	for (GroomBundle *bundle = groom->bundles.first; bundle; bundle = bundle->next)
	{
		if (bundle->curvecache)
		{
			const int totcurvecache = (bundle->numshapeverts + 1) * bundle->curvesize;
			GroomCurveCache *cc = bundle->curvecache;
			for (int i = 0; i < totcurvecache; ++i, ++cc)
			{
				minmax_v3v3_v3(min, max, cc->co);
			}
			
			result = true;
		}
	}
	return result;
}

BoundBox *BKE_groom_boundbox_get(Object *ob)
{
	BLI_assert(ob->type == OB_GROOM);
	Groom *groom = ob->data;

	if (ob->bb)
		return ob->bb;

	if (groom->bb == NULL || (groom->bb->flag & BOUNDBOX_DIRTY)) {
		BKE_groom_boundbox_calc(groom);
	}

	return groom->bb;
}

void BKE_groom_boundbox_calc(Groom *groom)
{
	if (groom->bb == NULL)
	{
		groom->bb = MEM_callocN(sizeof(BoundBox), "boundbox");
	}

	float min[3], max[3];
	INIT_MINMAX(min, max);
	if (!BKE_groom_minmax(groom, min, max)) {
		min[0] = min[1] = min[2] = -1.0f;
		max[0] = max[1] = max[2] = 1.0f;
	}

	BKE_boundbox_init_from_minmax(groom->bb, min, max);
	groom->bb->flag &= ~BOUNDBOX_DIRTY;
}


/* === Scalp regions === */

void BKE_groom_bind_scalp_regions(Groom *groom, bool force_rebind)
{
	if (groom->editgroom)
	{
		for (GroomBundle *bundle = groom->editgroom->bundles.first; bundle; bundle = bundle->next)
		{
			BKE_groom_bundle_bind(groom, bundle, force_rebind);
		}
	}
	else
	{
		for (GroomBundle *bundle = groom->bundles.first; bundle; bundle = bundle->next)
		{
			BKE_groom_bundle_bind(groom, bundle, force_rebind);
		}
	}
}

/* Returns the transform at the root of the bundle */
static bool groom_get_bundle_transform_on_scalp(const GroomBundle *bundle, const Mesh *scalp, float r_loc[3], float r_rot[3][3])
{
	const int numshapeverts = bundle->numshapeverts;
	if (numshapeverts == 0 || bundle->scalp_region == NULL)
	{
		return false;
	}
	
	/* last sample is the center position */
	const MeshSample *center_sample = &bundle->scalp_region[numshapeverts];
	
	if (BKE_mesh_sample_eval(scalp, center_sample, r_loc, r_rot[2], r_rot[1]))
	{
		cross_v3_v3v3(r_rot[0], r_rot[1], r_rot[2]);
		return true;
	}
	else
	{
		zero_v3(r_loc);
		unit_m3(r_rot);
		return false;
	}
}

static bool groom_shape_rebuild(GroomBundle *bundle, Object *scalp_ob)
{
	BLI_assert(bundle->scalp_region != NULL);
	BLI_assert(scalp_ob->type == OB_MESH);
	const int numshapeverts = bundle->numshapeverts;
	
	bool result = true;
	float (*shape)[2] = MEM_mallocN(sizeof(*shape) * numshapeverts, "groom section shape");
	
	Mesh *me = scalp_ob->data;
	
	float center_loc[3];
	float center_mat[3][3];
	if (!groom_get_bundle_transform_on_scalp(bundle, me, center_loc, center_mat))
	{
		result = false;
		goto cleanup;
	}
	
	MeshSample *sample = bundle->scalp_region;
	GroomSectionVertex *vert0 = bundle->verts;
	for (int i = 0; i < numshapeverts; ++i, ++sample, ++vert0)
	{
		/* 3D position of the shape vertex origin on the mesh */
		float co[3], nor[3], tang[3];
		if (!BKE_mesh_sample_eval(me, sample, co, nor, tang))
		{
			result = false;
			goto cleanup;
		}
		/* Get relative offset from the center */
		sub_v3_v3(co, center_loc);
		/* Convert mesh surface positions to 2D shape
		 * by projecting onto the normal plane
		 */
		shape[i][0] = dot_v3v3(co, center_mat[0]);
		shape[i][1] = dot_v3v3(co, center_mat[1]);
	}
	
	bundle->totverts = numshapeverts * bundle->totsections;
	bundle->verts = MEM_reallocN_id(bundle->verts, sizeof(*bundle->verts) * bundle->totverts, "groom bundle vertices");
	/* Set the shape for all sections */
	GroomSectionVertex *vert = bundle->verts;
	for (int i = 0; i < bundle->totsections; ++i)
	{
		for (int j = 0; j < numshapeverts; ++j, ++vert)
		{
			copy_v2_v2(vert->co, shape[j]);
			vert->flag = 0;
		}
	}
	
cleanup:
	MEM_freeN(shape);
	
	return result;
}

static BMesh *groom_create_scalp_bmesh(Mesh *me)
{
	const BMAllocTemplate allocsize = BMALLOC_TEMPLATE_FROM_ME(me);
	
	BMesh *bm = BM_mesh_create(&allocsize, &((struct BMeshCreateParams){
	        .use_toolflags = true,
	         }));
	
	BM_mesh_bm_from_me(bm, me, (&(struct BMeshFromMeshParams){
	        .calc_face_normal = true, .use_shapekey = false,
	        }));
	
	return bm;
}

static bool groom_bundle_region_from_mesh_fmap(GroomBundle *bundle, Object *scalp_ob)
{
	BLI_assert(scalp_ob->type == OB_MESH);
	
	BKE_groom_bundle_curve_cache_clear(bundle);
	
	Mesh *me = scalp_ob->data;
	const int scalp_fmap_nr = BKE_object_facemap_name_index(scalp_ob, bundle->scalp_facemap_name);
	const int cd_fmap_offset = CustomData_get_offset(&me->pdata, CD_FACEMAP);
	if (scalp_fmap_nr < 0 || cd_fmap_offset < 0)
	{
		return false;
	}
	
	BMesh *bm = groom_create_scalp_bmesh(me);
	bool result = true;
	
	/* Tag faces in the face map for the BMO walker */
	{
		BMFace *f;
		BMIter iter;
		BM_ITER_MESH(f, &iter, bm, BM_FACES_OF_MESH)
		{
			int fmap = BM_ELEM_CD_GET_INT(f, cd_fmap_offset);
			BM_elem_flag_set(f, BM_ELEM_TAG, fmap == scalp_fmap_nr);
		}
	}
	
	BMOperator op;
	BMO_op_initf(bm, &op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
	             "face_island_boundary faces=%hf", BM_ELEM_TAG);
	
	BMO_op_exec(bm, &op);
	if (BMO_error_occurred(bm))
	{
		result = false;
		goto finalize;
	}
	
	const int numshapeverts = bundle->numshapeverts = BMO_slot_buffer_count(op.slots_out, "boundary");
	bundle->scalp_region = MEM_callocN(sizeof(*bundle->scalp_region) * (numshapeverts + 1), "groom bundle scalp region");
	
	float center_co[3]; /* average vertex location for placing the center */
	{
		BMLoop *l;
		BMOIter oiter;
		MeshSample *sample = bundle->scalp_region;
		zero_v3(center_co);
		BMO_ITER (l, &oiter, op.slots_out, "boundary", BM_LOOP)
		{
			sample->orig_poly = BM_elem_index_get(l->f);
			sample->orig_loops[0] = BM_elem_index_get(l);
			sample->orig_verts[0] = BM_elem_index_get(l->v);
			sample->orig_weights[0] = 1.0f;
			BLI_assert(BKE_mesh_sample_is_valid(sample));
			
			add_v3_v3(center_co, l->v->co);
			
			++sample;
		}
		if (numshapeverts > 0)
		{
			mul_v3_fl(center_co, 1.0f / (float)(numshapeverts));
		}
	}
	
	{
		/* BVH tree for binding the region center location */
		BVHTreeFromMesh bvhtree;
		BKE_bvhtree_from_mesh_get(&bvhtree, me, BVHTREE_FROM_LOOPTRI, 2);
		if (bvhtree.tree != NULL) {
			BVHTreeNearest nearest;
			nearest.index = -1;
			nearest.dist_sq = FLT_MAX;
			
			BLI_bvhtree_find_nearest(bvhtree.tree, center_co, &nearest, bvhtree.nearest_callback, &bvhtree);
			if (nearest.index >= 0)
			{
				/* last sample is the center position */
				MeshSample *center_sample = &bundle->scalp_region[numshapeverts];
				BKE_mesh_sample_weights_from_loc(center_sample, me, nearest.index, nearest.co);
				BLI_assert(BKE_mesh_sample_is_valid(center_sample));
			}
		}
		else
		{
			result = false;
		}
	
		free_bvhtree_from_mesh(&bvhtree);
	}
	
finalize:
	if (result == true)
	{
		groom_shape_rebuild(bundle, scalp_ob);
	}
	else
	{
		if (bundle->scalp_region)
		{
			MEM_freeN(bundle->scalp_region);
			bundle->scalp_region = NULL;
		}
	}
	
	BMO_op_finish(bm, &op);
	BM_mesh_free(bm);
	
	return result;
}

bool BKE_groom_bundle_bind(Groom *groom, GroomBundle *bundle, bool force_rebind)
{
	if (bundle->scalp_region && !force_rebind)
	{
		return true;
	}
	
	BKE_groom_bundle_unbind(bundle);
	if (!groom->scalp_object)
	{
		return false;
	}
	if (!BKE_object_facemap_find_name(groom->scalp_object, bundle->scalp_facemap_name))
	{
		return false;
	}
	
	if (groom->scalp_object->type == OB_MESH)
	{
		groom_bundle_region_from_mesh_fmap(bundle, groom->scalp_object);
	}
	
	return (bundle->scalp_region != NULL);
}

void BKE_groom_bundle_unbind(GroomBundle *bundle)
{
	if (bundle->scalp_region)
	{
		MEM_freeN(bundle->scalp_region);
		bundle->scalp_region = NULL;
	}
}


/* === Constraints === */

/* Apply constraints on groom geometry */
void BKE_groom_apply_constraints(Groom *groom, Mesh *scalp)
{
	ListBase *bundles = (groom->editgroom ? &groom->editgroom->bundles : &groom->bundles);
	for (GroomBundle *bundle = bundles->first; bundle; bundle = bundle->next)
	{
		if (bundle->totsections > 0)
		{
			GroomSection *section = &bundle->sections[0];
			
			if (scalp)
			{
				/* For bound regions the bundle should be attached to the scalp */
				groom_get_bundle_transform_on_scalp(bundle, scalp, section->center, section->mat);
			}
			else
			{
				if (bundle->totsections > 1)
				{
					/* align to the first segment */
					float dir[3];
					sub_v3_v3v3(dir, (section+1)->center, section->center);
					normalize_v3(dir);
					
					const float dir_prev[3] = {0.0f, 1.0f, 0.0f};
					rotation_between_vecs_to_mat3(section->mat, dir_prev, dir);
				}
				else
				{
					unit_m3(section->mat);
				}
			}
		}
	}
}


/* === Hair System === */

/* Set loop weights for all faces covered by the bundle region */
static bool groom_add_bundle_loop_weights(const Groom *groom, const GroomBundle *bundle,
                                          float *loop_weights)
{
	BLI_assert(groom->scalp_object && groom->scalp_object->type == OB_MESH);
	const Mesh *scalp = groom->scalp_object->data;
	
	const int fmap_nr = BKE_object_facemap_name_index(groom->scalp_object, bundle->scalp_facemap_name);
	if (fmap_nr < 0)
	{
		return false;
	}
	
	if (!CustomData_has_layer(&scalp->pdata, CD_FACEMAP))
	{
		return false;
	}
	
	const int *map = CustomData_get_layer(&scalp->pdata, CD_FACEMAP);
	if (!map)
	{
		return false;
	}
	
	/* Find all polys of this bundle's face map */
	const MPoly *mp = scalp->mpoly;
	for (int i = 0; i < scalp->totpoly; ++i, ++mp)
	{
		int fmap_value = map[i];
		if (fmap_value == fmap_nr)
		{
			/* Include poly by setting the weight for all its loops */
			for (int j = 0; j < mp->totloop; ++j)
			{
				loop_weights[mp->loopstart + j] = 1.0f;
			}
		}
	}
	
	return true;
}

/* Set loop weights for all faces covered */
static void groom_add_all_loop_weights(const Groom *groom,
                                       float *loop_weights)
{
	for (GroomBundle *bundle = groom->bundles.first; bundle; bundle = bundle->next)
	{
		groom_add_bundle_loop_weights(groom, bundle, loop_weights);
	}
}

/* Distribute points on the scalp to use as guide curve origins */
static void groom_generate_guides(
        Groom *groom,
        Mesh *scalp,
        unsigned int seed)
{
	const size_t loop_weights_size = sizeof(float) * scalp->totloop;
	float *loop_weights = MEM_mallocN(loop_weights_size, "groom scalp loop weights");
	
	int i = 0;
	GroomBundle *bundle = groom->bundles.first;
	for (; bundle; bundle = bundle->next, ++i)
	{
		bundle->guides = MEM_reallocN_id(bundle->guides, sizeof(GroomHairGuide) * bundle->guides_count, "groom bundle hair guides");
		
		/* Mask for the scalp region */
		memset(loop_weights, 0, loop_weights_size);
		groom_add_bundle_loop_weights(groom, bundle, loop_weights);
		
		{
			const int count = bundle->guides_count;
			unsigned int region_seed = BLI_ghashutil_combine_hash(seed, BLI_ghashutil_uinthash(i));
			float scalp_area = BKE_hair_calc_surface_area(scalp);
			float density = BKE_hair_calc_density_from_count(scalp_area, count);
			float min_distance = BKE_hair_calc_min_distance_from_density(density);
			MeshSampleGenerator *gen = BKE_mesh_sample_gen_surface_poissondisk(
			                               region_seed,
			                               min_distance,
			                               count,
			                               loop_weights);
			
			BKE_mesh_sample_generator_bind(gen, scalp);
			
			static const bool use_threads = false;
			bundle->totguides = BKE_mesh_sample_generate_batch_ex(
			                        gen,
			                        &bundle->guides->root,
			                        sizeof(GroomHairGuide),
			                        count,
			                        use_threads);
			
			BKE_mesh_sample_free_generator(gen);
		}
		
		/* Calculate weights for interpolating the guide between shape vertices */
		{
			bundle->guide_shape_weights = MEM_reallocN_id(bundle->guide_shape_weights, sizeof(float) * bundle->totguides * bundle->numshapeverts, "groom guide shape weights");
			
			/* Use first section as shape for computing weights */
			const int shapesize = bundle->numshapeverts;
			float (*shape)[2] = MEM_mallocN(sizeof(float) * 2 * shapesize, "bundle shape verts");
			/* Expecting at least one section */
			BLI_assert(bundle->totsections >= 1 && bundle->totverts >= shapesize);
			/* Need a dense array for the interp_weights_poly_v2 function */
			for (int j = 0; j < shapesize; ++j)
			{
				copy_v2_v2(shape[j], bundle->verts[j].co);
			}
			
			float *w = bundle->guide_shape_weights;
			for (int j = 0; j < bundle->totguides; ++j)
			{
				/* Define interpolation point by projecting the guide root
				 * onto the bundle base plane.
				 */
				float p[3], n[3], t[3];
				BKE_mesh_sample_eval(scalp, &bundle->guides[j].root, p, n, t);
				sub_v3_v3(p, bundle->sections[0].center);
				mul_transposed_m3_v3(bundle->sections[0].mat, p);
				
				/* Calculate interpolation weights */
				interp_weights_poly_v2(w, shape, shapesize, p);
				
				w += bundle->numshapeverts;
			}
			
			MEM_freeN(shape);
		}
	}
	
	MEM_freeN(loop_weights);
}

void BKE_groom_hair_distribute(Groom *groom, unsigned int seed, int hair_count)
{
	struct HairSystem *hsys = groom->hair_system;
	
	BLI_assert(groom->scalp_object && groom->scalp_object->type == OB_MESH);
	Mesh *scalp = groom->scalp_object->data;
	
	{
		/* Per-loop weights for limiting follicles to covered faces */
		float *loop_weights = MEM_callocN(sizeof(float) * scalp->totloop, "groom scalp loop weights");
		groom_add_all_loop_weights(groom, loop_weights);
		
		unsigned int hair_seed = BLI_ghashutil_combine_hash(seed, BLI_ghashutil_strhash("groom hair follicles"));
		BKE_hair_generate_follicles_ex(hsys, scalp, hair_seed, hair_count, loop_weights);
		
		MEM_freeN(loop_weights);
	}
	
	{
		unsigned int guides_seed = BLI_ghashutil_combine_hash(seed, BLI_ghashutil_strhash("groom guide curves"));
		groom_generate_guides(groom, scalp, guides_seed);
	}
}

void BKE_groom_hair_update_guide_curves(Groom *groom)
{
	struct HairSystem *hsys = groom->hair_system;
	const ListBase *bundles = groom->editgroom ? &groom->editgroom->bundles : &groom->bundles;
	
	/* Count guides for all regions combined */
	int totguides = 0;
	for (const GroomBundle *bundle = bundles->first; bundle; bundle = bundle->next)
	{
		totguides += bundle->totguides;
	}
	
	/* First declare all guide curves and lengths */
	BKE_hair_guide_curves_begin(hsys, totguides);
	for (const GroomBundle *bundle = bundles->first; bundle; bundle = bundle->next)
	{
		const int curvesize = bundle->curvesize;
		for (int i = 0; i < bundle->totguides; ++i)
		{
			BKE_hair_set_guide_curve(hsys, i, &bundle->guides[i].root, curvesize);
		}
	}
	BKE_hair_guide_curves_end(hsys);
	
	int idx = 0;
	for (const GroomBundle *bundle = bundles->first; bundle; bundle = bundle->next)
	{
		const int shapesize = bundle->numshapeverts;
		const int curvesize = bundle->curvesize;
		const float *weights = bundle->guide_shape_weights;
		float co[3];
		for (int guide_idx = 0; guide_idx < bundle->totguides; ++guide_idx)
		{
			for (int j = 0; j < curvesize; ++j)
			{
				/* Compute barycentric guide location from shape */
				zero_v3(co);
				for (int k = 0; k < shapesize; ++k)
				{
					madd_v3_v3fl(co, bundle->curvecache[k * curvesize + j].co, weights[k]);
				}
				
				BKE_hair_set_guide_vertex(hsys, idx, 0, co);
				++idx;
			}
			
			weights += shapesize;
		}
	}
	
	if (groom->scalp_object)
	{
		BLI_assert(groom->scalp_object->type == OB_MESH);
		Mesh *scalp = groom->scalp_object->data;
		BKE_hair_bind_follicles(hsys, scalp);
	}
}


/* === Curve cache === */

/* forward differencing method for cubic polynomial eval */
static void groom_forward_diff_cubic(float a, float b, float c, float d, float *p, int it, int stride)
{
	float f = (float)it;
	a *= 1.0f / (f*f*f);
	b *= 1.0f / (f*f);
	c *= 1.0f / (f);
	
	float q0 = d;
	float q1 = a + b + c;
	float q2 = 6 * a + 2 * b;
	float q3 = 6 * a;

	for (int i = 0; i <= it; i++) {
		*p = q0;
		p = POINTER_OFFSET(p, stride);
		q0 += q1;
		q1 += q2;
		q2 += q3;
	}
}

/* cubic bspline section eval */
static void groom_eval_curve_cache_section(GroomCurveCache *cache, int curve_res, const float *co0, const float *co1, const float *co2, const float *co3)
{
	float a, b, c, d;
	for (int k = 0; k < 3; ++k)
	{
		/* define tangents from segment direction */
		float n1, n2;
		
		if (co0)
		{
			n1 = 0.5f * (co2[k] - co0[k]);
		}
		else
		{
			n1 = co2[k] - co1[k];
		}
		
		if (co3)
		{
			n2 = 0.5f * (co3[k] - co1[k]);
		}
		else
		{
			n2 = co2[k] - co1[k];
		}
		
		/* Hermite spline interpolation */
		a = 2.0f * (co1[k] - co2[k]) + n1 + n2;
		b = 3.0f * (co2[k] - co1[k]) - 2.0f * n1 - n2;
		c = n1;
		d = co1[k];
		
		groom_forward_diff_cubic(a, b, c, d, cache->co + k, curve_res, sizeof(*cache));
	}
}

static void groom_eval_center_curve_section(
        GroomBundle *bundle,
        int curve_res)
{
	BLI_assert(bundle->totsections >= 2);
	BLI_assert(curve_res >= 1);
	
	/* last cache curve is the center curve */
	GroomCurveCache *cache = bundle->curvecache + bundle->curvesize * bundle->numshapeverts;
	for (int i = 0; i < bundle->totsections-1; ++i, cache += curve_res)
	{
		const GroomSection *section = &bundle->sections[i];
		const float *co0 = (i > 0) ? section[-1].center : NULL;
		const float *co1 = section[0].center;
		const float *co2 = section[1].center;
		const float *co3 = (i < bundle->totsections - 2) ? section[2].center : NULL;
		groom_eval_curve_cache_section(cache, curve_res, co0, co1, co2, co3);
	}
}

static void groom_eval_shape_vertex(const GroomBundle *bundle, const Mesh *scalp, int vertex_idx, int section_idx, float r_co[3])
{
	if (section_idx == 0 && bundle->scalp_region != NULL)
	{
		/* For bound regions use location on the scalp */
		const MeshSample *bound_loc = &bundle->scalp_region[vertex_idx];
		float nor[3], tang[3];
		BKE_mesh_sample_eval(scalp, bound_loc, r_co, nor, tang);
	}
	else
	{
		const GroomSection *section = &bundle->sections[section_idx];
		const GroomSectionVertex *vertex = &bundle->verts[vertex_idx + section_idx * bundle->numshapeverts];
		
		float tmp[3] = {0.0f, 0.0f, 0.0f};
		copy_v2_v2(tmp, vertex->co);
		mul_v3_m3v3(r_co, section->mat, tmp);
		add_v3_v3(r_co, section->center);
	}
}

static void groom_eval_shape_curves(
        GroomBundle *bundle,
        const Mesh *scalp,
        int curve_res)
{
	BLI_assert(bundle->totsections >= 2);
	BLI_assert(curve_res >= 1);
	
	for (int i = 0; i < bundle->numshapeverts; ++i)
	{
		GroomCurveCache *cache = bundle->curvecache + i * bundle->curvesize;
		for (int j = 0; j < bundle->totsections-1; ++j, cache += curve_res)
		{
			const float *co0 = NULL, *co1 = NULL, *co2 =NULL, *co3 = NULL;
			
			float vec0[3], vec1[3], vec2[3], vec3[3];
			if (j > 0)
			{
				groom_eval_shape_vertex(bundle, scalp, i, j-1, vec0);
				co0 = vec0;
			}
			{
				groom_eval_shape_vertex(bundle, scalp, i, j, vec1);
				co1 = vec1;
			}
			{
				groom_eval_shape_vertex(bundle, scalp, i, j+1, vec2);
				co2 = vec2;
			}
			if (j < bundle->totsections - 2)
			{
				groom_eval_shape_vertex(bundle, scalp, i, j+2, vec3);
				co3 = vec3;
			}
			
			groom_eval_curve_cache_section(cache, curve_res, co0, co1, co2, co3);
		}
	}
}

static void groom_eval_curve_step(float mat[3][3], const float mat_prev[3][3], const float co0[3], const float co1[3])
{
	float dir[3];
	sub_v3_v3v3(dir, co1, co0);
	normalize_v3(dir);
	
	float dir_prev[3];
	normalize_v3_v3(dir_prev, mat_prev[2]);
	float rot[3][3];
	rotation_between_vecs_to_mat3(rot, dir_prev, dir);
	
	mul_m3_m3m3(mat, rot, mat_prev);
}

/* Computes rotation matrices for all but the first segment of a bundle */
static void groom_eval_section_mats(GroomBundle *bundle, int curve_res)
{
	const int curvesize = bundle->curvesize;
	const int numshapeverts = bundle->numshapeverts;
	
	GroomSection *section = bundle->sections;
	/* last curve cache is center curve */
	const GroomCurveCache *cache = bundle->curvecache + bundle->curvesize * numshapeverts;
	
	float mat[3][3];
	/* initialize matrix */
	copy_m3_m3(mat, section->mat);
	++cache;
	++section;
	
	for (int i = 1; i < curvesize - 1; ++i, ++cache)
	{
		/* align interior points to average of prev and next segment */
		groom_eval_curve_step(mat, mat, cache[-1].co, cache[+1].co);
		
		if (i % curve_res == 0)
		{
			/* set section matrix */
			copy_m3_m3(section->mat, mat);
			++section;
		}
	}
	
	/* align to last segment */
	groom_eval_curve_step(mat, mat, cache[-1].co, cache[0].co);
	/* last section is not handled in the loop */
	copy_m3_m3(section->mat, mat);
}

void BKE_groom_curve_cache_update(Groom *groom, const Mesh *scalp)
{
	ListBase *bundles = (groom->editgroom ? &groom->editgroom->bundles : &groom->bundles);
	
	for (GroomBundle *bundle = bundles->first; bundle; bundle = bundle->next)
	{
		const int totsections = bundle->totsections;
		const int numshapeverts = bundle->numshapeverts;
		const int curve_res = groom->curve_res;
		if (totsections == 0)
		{
			BKE_groom_bundle_curve_cache_clear(bundle);
			
			/* nothing to do */
			continue;
		}
		
		bundle->curvesize = (totsections-1) * curve_res + 1;
		bundle->totcurvecache = bundle->curvesize * (numshapeverts + 1);
		bundle->curvecache = MEM_reallocN_id(bundle->curvecache, sizeof(GroomCurveCache) * bundle->totcurvecache, "groom bundle curve cache");
		
		if (totsections == 1)
		{
			/* degenerate case */
			copy_v3_v3(bundle->curvecache[numshapeverts].co, bundle->sections[0].center);
			
			for (int i = 0; i < numshapeverts; ++i)
			{
				copy_v2_v2(bundle->curvecache[i].co, bundle->verts[i].co);
				bundle->curvecache[i].co[2] = 0.0f;
			}
			
			continue;
		}
		
		/* Calculate center curve */
		groom_eval_center_curve_section(bundle, curve_res);
		
		/* Calculate coordinate frames for sections */
		groom_eval_section_mats(bundle, curve_res);
		
		/* Calculate shape curves */
		groom_eval_shape_curves(bundle, scalp, curve_res);
	}
}

void BKE_groom_curve_cache_clear(Groom *groom)
{
	for (GroomBundle *bundle = groom->bundles.first; bundle; bundle = bundle->next)
	{
		BKE_groom_bundle_curve_cache_clear(bundle);
	}
	if (groom->editgroom)
	{
		for (GroomBundle *bundle = groom->editgroom->bundles.first; bundle; bundle = bundle->next)
		{
			BKE_groom_bundle_curve_cache_clear(bundle);
		}
	}
}


/* === Depsgraph evaluation === */

void BKE_groom_eval_geometry(const struct Depsgraph *depsgraph, Groom *groom)
{
	if (G.debug & G_DEBUG_DEPSGRAPH) {
		printf("%s on %s\n", __func__, groom->id.name);
	}
	
	/* Apply constraints from scalp mesh to the groom geometry */
	Mesh *scalp = groom->scalp_object ? groom->scalp_object->data : NULL;
	Mesh *scalp_eval = (Mesh *)DEG_get_evaluated_id(depsgraph, &scalp->id);
	BKE_groom_apply_constraints(groom, scalp_eval);
	
	/* calculate curves for interpolating shapes */
	BKE_groom_curve_cache_update(groom, scalp_eval);
	
	/* generate actual guide curves for hair */
	BKE_groom_hair_update_guide_curves(groom);
	
	if (groom->bb == NULL || (groom->bb->flag & BOUNDBOX_DIRTY)) {
		BKE_groom_boundbox_calc(groom);
	}
}


/* === Draw Cache === */

void (*BKE_groom_batch_cache_dirty_cb)(Groom* groom, int mode) = NULL;
void (*BKE_groom_batch_cache_free_cb)(Groom* groom) = NULL;

void BKE_groom_batch_cache_dirty(Groom* groom, int mode)
{
	if (groom->batch_cache)
	{
		BKE_groom_batch_cache_dirty_cb(groom, mode);
	}
}

void BKE_groom_batch_cache_free(Groom *groom)
{
	if (groom->batch_cache)
	{
		BKE_groom_batch_cache_free_cb(groom);
	}
}

/* === Utility functions === */

struct Mesh* BKE_groom_get_scalp(struct Groom *groom)
{
	return groom->scalp_object ? groom->scalp_object->data : NULL;
}
