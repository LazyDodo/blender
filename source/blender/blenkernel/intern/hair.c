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

/** \file blender/blenkernel/intern/hair.c
 *  \ingroup bke
 */

#include <limits.h>
#include <string.h>

#include "MEM_guardedalloc.h"

#include "BLI_kdtree.h"
#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_sort.h"
#include "BLI_string_utf8.h"
#include "BLI_string_utils.h"

#include "DNA_hair_types.h"
#include "DNA_object_types.h"

#include "BKE_DerivedMesh.h"
#include "BKE_cdderivedmesh.h"
#include "BKE_hair.h"
#include "BKE_library.h"
#include "BKE_mesh.h"
#include "BKE_mesh_sample.h"

#include "BLT_translation.h"

HairSystem* BKE_hair_new(void)
{
	HairSystem *hair = MEM_callocN(sizeof(HairSystem), "hair system");
	
	hair->pattern = MEM_callocN(sizeof(HairPattern), "hair pattern");
	
	hair->material_index = 1;
	
	return hair;
}

HairSystem* BKE_hair_copy(HairSystem *hsys)
{
	HairSystem *nhsys = MEM_dupallocN(hsys);
	
	if (hsys->pattern)
	{
		nhsys->pattern = MEM_dupallocN(hsys->pattern);
		nhsys->pattern->follicles = MEM_dupallocN(hsys->pattern->follicles);
	}
	
	if (hsys->guides.curves)
	{
		nhsys->guides.curves = MEM_dupallocN(hsys->guides.curves);
	}
	if (hsys->guides.verts)
	{
		nhsys->guides.verts = MEM_dupallocN(hsys->guides.verts);
	}
	
	nhsys->draw_batch_cache = NULL;
	nhsys->draw_texture_cache = NULL;
	
	return nhsys;
}

void BKE_hair_free(HairSystem *hsys)
{
	BKE_hair_batch_cache_free(hsys);
	
	if (hsys->guides.curves)
	{
		MEM_freeN(hsys->guides.curves);
	}
	if (hsys->guides.verts)
	{
		MEM_freeN(hsys->guides.verts);
	}
	
	if (hsys->pattern)
	{
		if (hsys->pattern->follicles)
		{
			MEM_freeN(hsys->pattern->follicles);
		}
		MEM_freeN(hsys->pattern);
	}
	
	MEM_freeN(hsys);
}

/* Calculate surface area of a scalp mesh */
float BKE_hair_calc_surface_area(struct DerivedMesh *scalp)
{
	BLI_assert(scalp != NULL);
	
	int numpolys = scalp->getNumPolys(scalp);
	MPoly *mpolys = scalp->getPolyArray(scalp);
	MLoop *mloops = scalp->getLoopArray(scalp);
	MVert *mverts = scalp->getVertArray(scalp);

	float area = 0.0f;
	for (int i = 0; i < numpolys; ++i)
	{
		area += BKE_mesh_calc_poly_area(&mpolys[i], mloops + mpolys[i].loopstart, mverts);
	}
	return area;
}

/* Calculate a density value based on surface area and sample count */
float BKE_hair_calc_density_from_count(float area, int count)
{
	return area > 0.0f ? count / area : 0.0f;
}

/* Calculate maximum sample count based on surface area and density */
int BKE_hair_calc_max_count_from_density(float area, float density)
{
	return (int)(density * area);
}

/* Calculate a density value based on a minimum distance */
float BKE_hair_calc_density_from_min_distance(float min_distance)
{
	// max. circle packing density (sans pi factor): 1 / (2 * sqrt(3))
	static const float max_factor = 0.288675135;
	
	return min_distance > 0.0f ? max_factor / (min_distance * min_distance) : 0.0f;
}

/* Calculate a minimum distance based on density */
float BKE_hair_calc_min_distance_from_density(float density)
{
	// max. circle packing density (sans pi factor): 1 / (2 * sqrt(3))
	static const float max_factor = 0.288675135;
	
	return density > 0.0f ? sqrt(max_factor / density) : 0.0f;
}

/* Distribute hair follicles on a scalp mesh */
void BKE_hair_generate_follicles(
        HairSystem* hsys,
        struct DerivedMesh *scalp,
        unsigned int seed,
        int count)
{
	HairPattern *pattern = hsys->pattern;
	
	// Limit max_count to theoretical limit based on area
	float scalp_area = BKE_hair_calc_surface_area(scalp);
	float density = BKE_hair_calc_density_from_count(scalp_area, count);
	float min_distance = BKE_hair_calc_min_distance_from_density(density);
	
	if (pattern->follicles)
	{
		MEM_freeN(pattern->follicles);
	}
	pattern->follicles = MEM_callocN(sizeof(HairFollicle) * count, "hair follicles");
	
	{
		MeshSampleGenerator *gen = BKE_mesh_sample_gen_surface_poissondisk(seed, min_distance, count, NULL, NULL);
		
		BKE_mesh_sample_generator_bind(gen, scalp);
		
		static const bool use_threads = false;
		pattern->num_follicles = BKE_mesh_sample_generate_batch_ex(
		            gen,
		            &pattern->follicles->mesh_sample,
		            sizeof(HairFollicle),
		            count,
		            use_threads);
		
		BKE_mesh_sample_free_generator(gen);
	}
	
	hsys->flag |= HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING;
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

/* ================================= */

void BKE_hair_guide_curves_begin(HairSystem *hsys, int totcurves)
{
	if (totcurves != hsys->guides.totcurves)
	{
		hsys->guides.curves = MEM_reallocN(hsys->guides.curves, sizeof(HairGuideCurve) * totcurves);
		hsys->guides.totcurves = totcurves;

		hsys->flag |= HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING;
		BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
	}
}

void BKE_hair_set_guide_curve(HairSystem *hsys, int index, const MeshSample *mesh_sample, int numverts)
{
	BLI_assert(index <= hsys->guides.totcurves);
	
	HairGuideCurve *curve = &hsys->guides.curves[index];
	memcpy(&curve->mesh_sample, mesh_sample, sizeof(MeshSample));
	curve->numverts = numverts;
	
	hsys->flag |= HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING;
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

/* Calculate vertex start indices on all curves based on length.
 * Returns the total number of vertices.
 */
static int hair_guide_calc_vertstart(HairSystem *hsys)
{
	/* Recalculate vertex count and start offsets in curves */
	int vertstart = 0;
	for (int i = 0; i < hsys->guides.totcurves; ++i)
	{
		hsys->guides.curves[i].vertstart = vertstart;
		vertstart += hsys->guides.curves[i].numverts;
	}
	
	return vertstart;
}

void BKE_hair_guide_curves_end(HairSystem *hsys)
{
	const int totverts = hair_guide_calc_vertstart(hsys);

	if (totverts != hsys->guides.totverts)
	{
		hsys->guides.verts = MEM_reallocN(hsys->guides.verts, sizeof(HairGuideVertex) * totverts);
		hsys->guides.totverts = totverts;

		BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
	}
}

void BKE_hair_set_guide_vertex(HairSystem *hsys, int index, int flag, const float co[3])
{
	BLI_assert(index <= hsys->guides.totverts);
	
	HairGuideVertex *vertex = &hsys->guides.verts[index];
	vertex->flag = flag;
	copy_v3_v3(vertex->co, co);
	
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

void BKE_hair_set_hair_guides(HairSystem *hsys, HairGuideData *guides)
{
	if (hsys->guides.curves)
	{
		MEM_freeN(hsys->guides.curves);
	}
	hsys->guides.curves = MEM_dupallocN(hsys->guides.curves);
	hsys->guides.totcurves = guides->totcurves;

	if (hsys->guides.verts)
	{
		MEM_freeN(hsys->guides.verts);
	}
	hsys->guides.verts = MEM_dupallocN(hsys->guides.verts);
	hsys->guides.totverts = guides->totverts;

#ifndef NDEBUG
	const int vertcount = hair_guide_calc_vertstart(hsys);
	BLI_assert(vertcount <= hsys->guides.totverts);
#endif

	hsys->flag |= HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING;
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

void BKE_hair_clear_guides(HairSystem *hsys)
{
	if (hsys->guides.curves)
	{
		MEM_freeN(hsys->guides.curves);
		hsys->guides.curves = NULL;
	}
	hsys->guides.totcurves = 0;

	if (hsys->guides.verts)
	{
		MEM_freeN(hsys->guides.verts);
		hsys->guides.verts = NULL;
	}
	hsys->guides.totverts = 0;

	hsys->flag &= ~HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING;
	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
}

/* ================================= */

BLI_INLINE void hair_fiber_verify_weights(HairFollicle *follicle)
{
	const float *w = follicle->parent_weight;
	
	BLI_assert(w[0] >= 0.0f && w[1] >= 0.0f && w[2] >= 0.0f && w[3] >= 0.0f);
	float sum = w[0] + w[1] + w[2] + w[3];
	float epsilon = 1.0e-2;
	BLI_assert(sum > 1.0f - epsilon && sum < 1.0f + epsilon);
	UNUSED_VARS(sum, epsilon);
	
	BLI_assert(w[0] >= w[1] && w[1] >= w[2] && w[2] >= w[3]);
}

static void hair_fiber_sort_weights(HairFollicle *follicle)
{
	unsigned int *idx = follicle->parent_index;
	float *w = follicle->parent_weight;

#define FIBERSWAP(a, b) \
	SWAP(unsigned int, idx[a], idx[b]); \
	SWAP(float, w[a], w[b]);

	for (int k = 0; k < 3; ++k) {
		int maxi = k;
		float maxw = w[k];
		for (int i = k+1; i < 4; ++i) {
			if (w[i] > maxw) {
				maxi = i;
				maxw = w[i];
			}
		}
		if (maxi != k)
			FIBERSWAP(k, maxi);
	}
	
#undef FIBERSWAP
}

static void hair_fiber_find_closest_strand(
        HairFollicle *follicle,
        const float loc[3],
        const KDTree *tree,
        const float (*strandloc)[3])
{
	/* Use the 3 closest strands for interpolation.
	 * Note that we have up to 4 possible weights, but we
	 * only look for a triangle with this method.
	 */
	KDTreeNearest nearest[3];
	const float *sloc[3] = {NULL};
	int k, found = BLI_kdtree_find_nearest_n(tree, loc, nearest, 3);
	for (k = 0; k < found; ++k) {
		follicle->parent_index[k] = (unsigned int)nearest[k].index;
		sloc[k] = strandloc[nearest[k].index];
	}
	for (; k < 4; ++k) {
		follicle->parent_index[k] = HAIR_STRAND_INDEX_NONE;
		follicle->parent_weight[k] = 0.0f;
	}
	
	/* calculate barycentric interpolation weights */
	if (found == 3) {
		float closest[3];
		closest_on_tri_to_point_v3(closest, loc, sloc[0], sloc[1], sloc[2]);
		
		float w[3];
		interp_weights_tri_v3(w, sloc[0], sloc[1], sloc[2], closest);
		copy_v3_v3(follicle->parent_weight, w);
		/* float precisions issues can cause slightly negative weights */
		CLAMP3(follicle->parent_weight, 0.0f, 1.0f);
	}
	else if (found == 2) {
		follicle->parent_weight[1] = line_point_factor_v3(loc, sloc[0], sloc[1]);
		follicle->parent_weight[0] = 1.0f - follicle->parent_weight[1];
		/* float precisions issues can cause slightly negative weights */
		CLAMP2(follicle->parent_weight, 0.0f, 1.0f);
	}
	else if (found == 1) {
		follicle->parent_weight[0] = 1.0f;
	}
	
	hair_fiber_sort_weights(follicle);
}

void BKE_hair_bind_follicles(HairSystem *hsys, DerivedMesh *scalp)
{
	if (!(hsys->flag & HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING))
	{
		return;
	}
	hsys->flag &= ~HAIR_SYSTEM_UPDATE_FOLLICLE_BINDING;
	
	HairPattern *pattern = hsys->pattern;
	const int num_strands = hsys->guides.totcurves;
	if (num_strands == 0 || !pattern)
		return;
	
	float (*strandloc)[3] = MEM_mallocN(sizeof(float) * 3 * num_strands, "strand locations");
	{
		for (int i = 0; i < num_strands; ++i) {
			float nor[3], tang[3];
			if (!BKE_mesh_sample_eval(scalp, &hsys->guides.curves[i].mesh_sample, strandloc[i], nor, tang)) {
				zero_v3(strandloc[i]);
			}
		}
	}
	
	KDTree *tree = BLI_kdtree_new(num_strands);
	for (int c = 0; c < num_strands; ++c) {
		BLI_kdtree_insert(tree, c, strandloc[c]);
	}
	BLI_kdtree_balance(tree);
	
	HairFollicle *follicle = pattern->follicles;
	for (int i = 0; i < pattern->num_follicles; ++i, ++follicle) {
		float loc[3], nor[3], tang[3];
		if (BKE_mesh_sample_eval(scalp, &follicle->mesh_sample, loc, nor, tang)) {
			hair_fiber_find_closest_strand(follicle, loc, tree, strandloc);
			hair_fiber_verify_weights(follicle);
		}
	}
	
	BLI_kdtree_free(tree);
	MEM_freeN(strandloc);
}

/* === Export === */

/* Returns number of vertices in a curve after subdivision */
BLI_INLINE int hair_get_strand_subdiv_length(int orig_length, int subdiv)
{
	return ((orig_length - 1) << subdiv) + 1;
}

/* Returns total number of vertices after subdivision */
BLI_INLINE int hair_get_strand_subdiv_numverts(int numstrands, int numverts, int subdiv)
{
	return ((numverts - numstrands) << subdiv) + numstrands;
}

/* Subdivide a curve */
static int hair_guide_subdivide(const HairGuideCurve* curve, const HairGuideVertex* verts, int subdiv, HairGuideVertex *r_verts)
{
	{
		/* Move vertex positions from the dense array to their initial configuration for subdivision. */
		const int step = (1 << subdiv);
		HairGuideVertex *dst = r_verts;
		for (int i = 0; i < curve->numverts; ++i) {
			copy_v3_v3(dst->co, verts[i].co);
			dst += step;
		}
	}
	
	/* Subdivide */
	for (int d = 0; d < subdiv; ++d) {
		const int num_edges = (curve->numverts - 1) << d;
		const int hstep = 1 << (subdiv - d - 1);
		const int step = 1 << (subdiv - d);
		
		/* Calculate edge points */
		{
			int index = 0;
			for (int k = 0; k < num_edges; ++k, index += step) {
				add_v3_v3v3(r_verts[index + hstep].co, r_verts[index].co, r_verts[index + step].co);
				mul_v3_fl(r_verts[index + hstep].co, 0.5f);
			}
		}
		
		/* Move original points */
		{
			int index = step;
			for (int k = 1; k < num_edges; ++k, index += step) {
				add_v3_v3v3(r_verts[index].co, r_verts[index - hstep].co, r_verts[index + hstep].co);
				mul_v3_fl(r_verts[index].co, 0.5f);
			}
		}
	}
	
	const int num_verts = ((curve->numverts - 1) << subdiv) + 1;
	return num_verts;
}

/* Calculate tangent and normal vector changes from one segment to the next */
static void hair_guide_transport_frame(const float co1[3], const float co2[3],
                                       float prev_tang[3], float prev_nor[3],
                                       float r_tang[3], float r_nor[3])
{
	/* segment direction */
	sub_v3_v3v3(r_tang, co2, co1);
	normalize_v3(r_tang);
	
	/* rotate the frame */
	float rot[3][3];
	rotation_between_vecs_to_mat3(rot, prev_tang, r_tang);
	mul_v3_m3v3(r_nor, rot, prev_nor);
	
	copy_v3_v3(prev_tang, r_tang);
	copy_v3_v3(prev_nor, r_nor);
}

/* Calculate tangent and normal vectors for all vertices on a curve */
static void hair_guide_calc_vectors(const HairGuideVertex* verts, int numverts, float rootmat[3][3],
                                    float (*r_tangents)[3], float (*r_normals)[3])
{
	BLI_assert(numverts >= 2);
	
	float prev_tang[3], prev_nor[3];
	
	copy_v3_v3(prev_tang, rootmat[2]);
	copy_v3_v3(prev_nor, rootmat[0]);
	
	hair_guide_transport_frame(
	        verts[0].co, verts[1].co,
	        prev_tang, prev_nor,
	        r_tangents[0], r_normals[0]);
	
	for (int i = 1; i < numverts - 1; ++i)
	{
		hair_guide_transport_frame(
		        verts[i-1].co, verts[i+1].co,
		        prev_tang, prev_nor,
		        r_tangents[i], r_normals[i]);
	}
	
	hair_guide_transport_frame(
	        verts[numverts-2].co, verts[numverts-1].co,
	        prev_tang, prev_nor,
	        r_tangents[numverts-1], r_normals[numverts-1]);
}

/* Create a new export cache.
 * This can be used to construct full fiber data for rendering.
 */

HairExportCache* BKE_hair_export_cache_new(void)
{
	HairExportCache *cache = MEM_callocN(sizeof(HairExportCache), "hair export cache");
	return cache;
}

/* Update an existing export cache when data is invalidated.
 * Returns flags for data that has been updated.
 */

int BKE_hair_export_cache_update(HairExportCache *cache, const HairSystem *hsys,
                                 int subdiv, DerivedMesh *scalp, int data)
{
	/* Check for missing data */
	data |= BKE_hair_export_cache_get_required_updates(cache);
	
	if (data & HAIR_EXPORT_GUIDE_CURVES)
	{
		/* Cache subdivided guide curves */
		const int totguidecurves = cache->totguidecurves = hsys->guides.totcurves;
		cache->guide_curves = MEM_reallocN_id(cache->guide_curves, sizeof(HairGuideCurve) * totguidecurves, "hair export guide curves");
		
		int totguideverts = 0;
		for (int i = 0; i < totguidecurves; ++i) {
			const HairGuideCurve *curve_orig = &hsys->guides.curves[i];
			HairGuideCurve *curve = &cache->guide_curves[i];
			
			memcpy(&curve->mesh_sample, &curve_orig->mesh_sample, sizeof(MeshSample));
			curve->numverts = hair_get_strand_subdiv_length(curve_orig->numverts, subdiv);
			curve->vertstart = totguideverts;
			
			totguideverts  += curve->numverts;
		}
		cache->totguideverts = totguideverts;
	}
	
	if (data & HAIR_EXPORT_GUIDE_VERTICES)
	{
		const int totguidecurves = cache->totguidecurves;
		const int totguideverts = cache->totguideverts;
		cache->guide_verts = MEM_reallocN_id(cache->guide_verts, sizeof(HairGuideVertex) * totguideverts, "hair export guide verts");
		cache->guide_tangents = MEM_reallocN_id(cache->guide_tangents, sizeof(float[3]) * totguideverts, "hair export guide tangents");
		cache->guide_normals = MEM_reallocN_id(cache->guide_normals, sizeof(float[3]) * totguideverts, "hair export guide normals");
		
		for (int i = 0; i < totguidecurves; ++i) {
			const HairGuideCurve *curve_orig = &hsys->guides.curves[i];
			const HairGuideVertex *verts_orig = &hsys->guides.verts[curve_orig->vertstart];
			const HairGuideCurve *curve = &cache->guide_curves[i];
			HairGuideVertex *verts = &cache->guide_verts[curve->vertstart];
			float (*tangents)[3] = &cache->guide_tangents[curve->vertstart];
			float (*normals)[3] = &cache->guide_normals[curve->vertstart];
			
			hair_guide_subdivide(curve_orig, verts_orig, subdiv, verts);
			
			{
				/* Root matrix for defining the initial normal direction */
				float rootpos[3];
				float rootmat[3][3];
				BKE_mesh_sample_eval(scalp, &curve->mesh_sample, rootpos, rootmat[2], rootmat[0]);
				cross_v3_v3v3(rootmat[1], rootmat[2], rootmat[0]);
				
				hair_guide_calc_vectors(verts, curve->numverts, rootmat, tangents, normals);
			}
		}
	}

	if (hsys->pattern)
	{
		if (data & HAIR_EXPORT_FOLLICLE_BINDING)
		{
			cache->follicles = hsys->pattern->follicles;
			cache->totfibercurves = hsys->pattern->num_follicles;
		}

		if (data & HAIR_EXPORT_FIBER_VERTEX_COUNTS)
		{
			/* Calculate the length of each fiber from the weighted average of its guide strands */
			const int totguidecurves = cache->totguidecurves;
			const int totfibercurves = cache->totfibercurves;
			
			cache->fiber_numverts = MEM_reallocN_id(cache->fiber_numverts, sizeof(int) * totfibercurves, "fiber numverts");
			cache->totfiberverts = 0;
			
			const HairFollicle *follicle = hsys->pattern->follicles;
			for (int i = 0; i < totfibercurves; ++i, ++follicle) {
				float fiblen = 0.0f;
				
				for (int k = 0; k < 4; ++k) {
					const int si = follicle->parent_index[k];
					const float sw = follicle->parent_weight[k];
					if (si == HAIR_STRAND_INDEX_NONE || sw == 0.0f) {
						break;
					}
					BLI_assert(si < totguidecurves);
					
					fiblen += (float)cache->guide_curves[si].numverts * sw;
				}
				
				/* Use rounded number of segments */
				const int numverts = (int)(fiblen + 0.5f);
				cache->fiber_numverts[i] = numverts;
				cache->totfiberverts += numverts;
			}
		}

		if (data & HAIR_EXPORT_FIBER_ROOT_POSITIONS)
		{
			const int totfibercurves = cache->totfibercurves;
			
			cache->fiber_root_position = MEM_reallocN_id(cache->fiber_root_position, sizeof(float[3]) * totfibercurves, "fiber root position");
			const HairFollicle *follicle = hsys->pattern->follicles;
			for (int i = 0; i < totfibercurves; ++i, ++follicle) {
				/* Cache fiber root position */
				float nor[3], tang[3];
				BKE_mesh_sample_eval(scalp, &follicle->mesh_sample, cache->fiber_root_position[i], nor, tang);
			}
		}
	}
	else
	{
		cache->follicles = NULL;
		cache->totfibercurves = 0;
		
		if (cache->fiber_numverts)
		{
			MEM_freeN(cache->fiber_numverts);
			cache->fiber_numverts = NULL;
		}
		if (cache->fiber_root_position)
		{
			MEM_freeN(cache->fiber_root_position);
			cache->fiber_root_position = NULL;
		}
	}
	
	return data;
}


/* Update an existing export cache when data is invalidated.
 * Returns flags for data that has been updated.
 * XXX Mesh-based version for Cycles export, until DerivedMesh->Mesh conversion is done.
 */

int BKE_hair_export_cache_update_mesh(HairExportCache *cache, const HairSystem *hsys,
                                      int subdiv, struct Mesh *scalp, int data)
{
	DerivedMesh *dm = CDDM_from_mesh(scalp);
	int result = BKE_hair_export_cache_update(cache, hsys, subdiv, dm, data);
	dm->release(dm);
	return result;
}

/* Free the given export cache */

void BKE_hair_export_cache_free(HairExportCache *cache)
{
	if (cache->fiber_numverts)
	{
		MEM_freeN(cache->fiber_numverts);
	}
	if (cache->guide_curves)
	{
		MEM_freeN(cache->guide_curves);
	}
	if (cache->guide_verts)
	{
		MEM_freeN(cache->guide_verts);
	}
	if (cache->guide_tangents)
	{
		MEM_freeN(cache->guide_tangents);
	}
	if (cache->guide_normals)
	{
		MEM_freeN(cache->guide_normals);
	}
	if (cache->fiber_root_position)
	{
		MEM_freeN(cache->fiber_root_position);
	}
	MEM_freeN(cache);
}

/* Returns flags for missing data parts */

int BKE_hair_export_cache_get_required_updates(const HairExportCache *cache)
{
	int data = 0;
	if (!cache->guide_curves)
	{
		data |= HAIR_EXPORT_GUIDE_CURVES;
	}
	if (!cache->guide_verts || !cache->guide_normals || !cache->guide_tangents)
	{
		data |= HAIR_EXPORT_GUIDE_VERTICES;
	}
	if (!cache->follicles)
	{
		data |= HAIR_EXPORT_FOLLICLE_BINDING;
	}
	if (!cache->fiber_root_position)
	{
		data |= HAIR_EXPORT_FIBER_ROOT_POSITIONS;
	}
	if (!cache->fiber_numverts)
	{
		data |= HAIR_EXPORT_FIBER_VERTEX_COUNTS;
	}
	return data;
}

/* Invalidate all data in a hair export cache */

void BKE_hair_export_cache_clear(HairExportCache *cache)
{
	/* Invalidate everything */
	BKE_hair_export_cache_invalidate(cache, HAIR_EXPORT_ALL);
}

/* Invalidate part of the data in a hair export cache */

void BKE_hair_export_cache_invalidate(HairExportCache *cache, int invalidate)
{
	if (invalidate & HAIR_EXPORT_GUIDE_CURVES)
	{
		if (cache->guide_curves)
		{
			MEM_freeN(cache->guide_curves);
			cache->guide_curves = 0;
		}
	}
	if (invalidate & HAIR_EXPORT_GUIDE_VERTICES)
	{
		if (cache->guide_verts)
		{
			MEM_freeN(cache->guide_verts);
			cache->guide_verts = NULL;
		}
		if (cache->guide_tangents)
		{
			MEM_freeN(cache->guide_tangents);
			cache->guide_tangents = NULL;
		}
		if (cache->guide_normals)
		{
			MEM_freeN(cache->guide_normals);
			cache->guide_tangents = NULL;
		}
	}
	if (invalidate & HAIR_EXPORT_FOLLICLE_BINDING)
	{
		cache->follicles = NULL;
	}
	if (invalidate & HAIR_EXPORT_FIBER_ROOT_POSITIONS)
	{
		if (cache->fiber_root_position)
		{
			MEM_freeN(cache->fiber_root_position);
			cache->fiber_root_position = NULL;
		}
	}
	if (invalidate & HAIR_EXPORT_FIBER_VERTEX_COUNTS)
	{
		if (cache->fiber_numverts)
		{
			MEM_freeN(cache->fiber_numverts);
			cache->fiber_numverts = NULL;
		}
	}
}
