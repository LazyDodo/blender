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

/** \file blender/blenkernel/intern/hair_draw.c
 *  \ingroup bke
 */

#include <string.h>

#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_kdtree.h"
#include "BLI_rand.h"

#include "DNA_hair_types.h"

#include "BKE_DerivedMesh.h"
#include "BKE_mesh_sample.h"
#include "BKE_hair.h"

/* === Draw Settings === */

HairDrawSettings* BKE_hair_draw_settings_new(void)
{
	HairDrawSettings *draw_settings = MEM_callocN(sizeof(HairDrawSettings), "hair draw settings");
	
	draw_settings->follicle_mode = HAIR_DRAW_FOLLICLE_NONE;
	
	return draw_settings;
}

HairDrawSettings* BKE_hair_draw_settings_copy(HairDrawSettings *draw_settings)
{
	HairDrawSettings *ndraw_settings = MEM_dupallocN(draw_settings);
	return ndraw_settings;
}

void BKE_hair_draw_settings_free(HairDrawSettings *draw_settings)
{
	MEM_freeN(draw_settings);
}

/* === Draw Cache === */

typedef struct HairFiberTextureBuffer {
	unsigned int parent_index[4];
	float parent_weight[4];
	float root_position[3];
	int pad;
} HairFiberTextureBuffer;
BLI_STATIC_ASSERT_ALIGN(HairFiberTextureBuffer, 8)

typedef struct HairStrandVertexTextureBuffer {
	float co[3];
	float nor[3];
	float tang[3];
	int pad;
} HairStrandVertexTextureBuffer;
BLI_STATIC_ASSERT_ALIGN(HairStrandVertexTextureBuffer, 8)

typedef struct HairStrandMapTextureBuffer {
	unsigned int vertex_start;
	unsigned int vertex_count;
} HairStrandMapTextureBuffer;
BLI_STATIC_ASSERT_ALIGN(HairStrandMapTextureBuffer, 8)

static void hair_strand_transport_frame(const float co1[3], const float co2[3],
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

static void hair_strand_calc_vectors(const HairGuideVertex* verts, int num_verts, float rootmat[3][3],
                                     HairStrandVertexTextureBuffer *strand)
{
	for (int i = 0; i < num_verts; ++i) {
		copy_v3_v3(strand[i].co, verts[i].co);
	}
	
	// Calculate tangent and normal vectors
	{
		BLI_assert(num_verts >= 2);
		
		float prev_tang[3], prev_nor[3];
		
		copy_v3_v3(prev_tang, rootmat[2]);
		copy_v3_v3(prev_nor, rootmat[0]);
		
		hair_strand_transport_frame(strand[0].co, strand[1].co,
		        prev_tang, prev_nor,
		        strand[0].tang, strand[0].nor);
		
		for (int i = 1; i < num_verts - 1; ++i)
		{
			hair_strand_transport_frame(strand[i-1].co, strand[i+1].co,
			        prev_tang, prev_nor,
			        strand[i].tang, strand[i].nor);
		}
		
		hair_strand_transport_frame(strand[num_verts-2].co, strand[num_verts-1].co,
		        prev_tang, prev_nor,
		        strand[num_verts-1].tang, strand[num_verts-1].nor);
	}
}

static void hair_get_strand_buffer(
        const HairExportCache *cache,
        DerivedMesh *scalp,
        HairStrandMapTextureBuffer *strand_map_buffer,
        HairStrandVertexTextureBuffer *strand_vertex_buffer)
{
	for (int i = 0; i < cache->totguidecurves; ++i) {
		const HairGuideCurve *curve = &cache->guide_curves[i];
		const HairGuideVertex *verts = &cache->guide_verts[curve->vertstart];
		HairStrandMapTextureBuffer *smap = &strand_map_buffer[i];
		HairStrandVertexTextureBuffer *svert = &strand_vertex_buffer[curve->vertstart];
		
		smap->vertex_start = curve->vertstart;
		smap->vertex_count = curve->numverts;
		
		{
			float pos[3];
			float matrix[3][3];
			BKE_mesh_sample_eval(scalp, &curve->mesh_sample, pos, matrix[2], matrix[0]);
			cross_v3_v3v3(matrix[1], matrix[2], matrix[0]);
			hair_strand_calc_vectors(verts, curve->numverts, matrix, svert);
		}
	}
}

static void hair_get_fiber_buffer(const HairExportCache *cache, DerivedMesh *scalp,
                                  HairFiberTextureBuffer *fiber_buf)
{
	const int totfibers = cache->totfibercurves;
	const HairFollicle *follicle = cache->follicles;
	HairFiberTextureBuffer *fb = fiber_buf;
	float nor[3], tang[3];
	for (int i = 0; i < totfibers; ++i, ++fb, ++follicle) {
		BKE_mesh_sample_eval(scalp, &follicle->mesh_sample, fb->root_position, nor, tang);
		for (int k = 0; k < 4; ++k)
		{
			fb->parent_index[k] = follicle->parent_index[k];
			fb->parent_weight[k] = follicle->parent_weight[k];
		}
	}
}

void BKE_hair_get_texture_buffer_size(
        const HairExportCache *cache,
        int *r_size,
        int *r_strand_map_start,
        int *r_strand_vertex_start,
        int *r_fiber_start)
{
	*r_strand_map_start = 0;
	*r_strand_vertex_start = *r_strand_map_start + cache->totguidecurves * sizeof(HairStrandMapTextureBuffer);
	*r_fiber_start = *r_strand_vertex_start + cache->totguideverts * sizeof(HairStrandVertexTextureBuffer);
	*r_size = *r_fiber_start + cache->totfibercurves * sizeof(HairFiberTextureBuffer);
}

void BKE_hair_get_texture_buffer(
        const HairExportCache *cache,
        DerivedMesh *scalp,
        void *buffer)
{
	int size, strand_map_start, strand_vertex_start, fiber_start;
	BKE_hair_get_texture_buffer_size(cache, &size, &strand_map_start, &strand_vertex_start, &fiber_start);
	
	if (scalp)
	{
		HairStrandMapTextureBuffer *strand_map = (HairStrandMapTextureBuffer*)((char*)buffer + strand_map_start);
		HairStrandVertexTextureBuffer *strand_verts = (HairStrandVertexTextureBuffer*)((char*)buffer + strand_vertex_start);
		HairFiberTextureBuffer *fibers = (HairFiberTextureBuffer*)((char*)buffer + fiber_start);
		
		hair_get_strand_buffer(
		            cache,
		            scalp,
		            strand_map,
		            strand_verts);
		hair_get_fiber_buffer(
		            cache,
		            scalp,
		            fibers);
	}
	else
	{
		memset(buffer, 0, size);
	}
}

void (*BKE_hair_batch_cache_dirty_cb)(HairSystem* hsys, int mode) = NULL;
void (*BKE_hair_batch_cache_free_cb)(HairSystem* hsys) = NULL;

void BKE_hair_batch_cache_dirty(HairSystem* hsys, int mode)
{
	if (hsys->draw_batch_cache) {
		BKE_hair_batch_cache_dirty_cb(hsys, mode);
	}
}

void BKE_hair_batch_cache_free(HairSystem* hsys)
{
	if (hsys->draw_batch_cache || hsys->draw_texture_cache) {
		BKE_hair_batch_cache_free_cb(hsys);
	}
}

/* === Render API === */

/* Calculate required size for render buffers. */
void BKE_hair_render_get_buffer_size(
        const HairExportCache* cache,
        int *r_totcurves,
        int *r_totverts)
{
	if (r_totcurves) *r_totcurves = cache->totfibercurves;
	if (r_totverts) *r_totverts = cache->totfiberverts;
}

/* Create render data in existing buffers.
 * Buffers must be large enough according to BKE_hair_get_render_buffer_size.
 */
void BKE_hair_render_fill_buffers(
        const HairExportCache* cache,
        int vertco_stride,
        int *r_curvestart,
        int *r_curvelen,
        float *r_vertco)
{
	
}
