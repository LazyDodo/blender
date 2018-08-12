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

#include "BKE_mesh_sample.h"
#include "BKE_hair.h"

/* === Draw Settings === */

HairDrawSettings* BKE_hair_draw_settings_new(void)
{
	HairDrawSettings *draw_settings = MEM_callocN(sizeof(HairDrawSettings), "hair draw settings");
	
	draw_settings->follicle_mode = HAIR_DRAW_FOLLICLE_POINTS;
	draw_settings->fiber_mode = HAIR_DRAW_FIBER_CURVES;
	draw_settings->shape_flag = HAIR_DRAW_CLOSE_TIP;
	draw_settings->shape = 0.0f;
	draw_settings->root_radius = 1.0f;
	draw_settings->tip_radius = 0.0f;
	draw_settings->radius_scale = 0.01f;
	
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
	float len;
} HairStrandVertexTextureBuffer;
BLI_STATIC_ASSERT_ALIGN(HairStrandVertexTextureBuffer, 8)

typedef struct HairStrandMapTextureBuffer {
	unsigned int vertex_start;
	unsigned int vertex_count;
	
	/* Shape attributes */
	float taper_length;         /* Distance at which final thickness is reached */
	float taper_thickness;      /* Relative thickness of the strand */
} HairStrandMapTextureBuffer;
BLI_STATIC_ASSERT_ALIGN(HairStrandMapTextureBuffer, 8)

static void hair_get_strand_buffer(
        const HairExportCache *cache,
        HairStrandMapTextureBuffer *strand_map_buffer,
        HairStrandVertexTextureBuffer *strand_vertex_buffer)
{
	for (int i = 0; i < cache->totcurves; ++i) {
		const HairFiberCurve *curve = &cache->fiber_curves[i];
		const HairFiberVertex *verts = &cache->fiber_verts[curve->vertstart];
		const float (*tangents)[3] = &cache->fiber_tangents[curve->vertstart];
		const float (*normals)[3] = &cache->fiber_normals[curve->vertstart];
		HairStrandMapTextureBuffer *smap = &strand_map_buffer[i];
		HairStrandVertexTextureBuffer *svert = &strand_vertex_buffer[curve->vertstart];
		
		smap->vertex_start = curve->vertstart;
		smap->vertex_count = curve->numverts;
		smap->taper_length = curve->taper_length;
		smap->taper_thickness = curve->taper_thickness;
		
		float len = 0.0f;
		for (int j = 0; j < curve->numverts; ++j)
		{
			copy_v3_v3(svert[j].co, verts[j].co);
			copy_v3_v3(svert[j].tang, tangents[j]);
			copy_v3_v3(svert[j].nor, normals[j]);
			
			if (j > 0)
			{
				len += len_v3v3(verts[j-1].co, verts[j].co);
			}
			svert[j].len = len;
		}
	}
}

static void hair_get_fiber_buffer(const HairExportCache *cache,
                                  HairFiberTextureBuffer *fiber_buf)
{
	const int totfibers = cache->totfollicles;
	const HairFollicle *follicle = cache->follicles;
	HairFiberTextureBuffer *fb = fiber_buf;
	for (int i = 0; i < totfibers; ++i, ++fb, ++follicle) {
		copy_v3_v3(fb->root_position, cache->follicle_root_position[i]);
		
		fb->parent_index[0] = follicle->curve;
		fb->parent_index[1] = HAIR_CURVE_INDEX_NONE;
		fb->parent_index[2] = HAIR_CURVE_INDEX_NONE;
		fb->parent_index[3] = HAIR_CURVE_INDEX_NONE;
		fb->parent_weight[0] = 1.0f;
		fb->parent_weight[1] = 0.0f;
		fb->parent_weight[2] = 0.0f;
		fb->parent_weight[3] = 0.0f;
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
	*r_strand_vertex_start = *r_strand_map_start + cache->totcurves * sizeof(HairStrandMapTextureBuffer);
	*r_fiber_start = *r_strand_vertex_start + cache->totverts * sizeof(HairStrandVertexTextureBuffer);
	*r_size = *r_fiber_start + cache->totfollicles * sizeof(HairFiberTextureBuffer);
}

void BKE_hair_get_texture_buffer(
        const HairExportCache *cache,
        void *buffer)
{
	int size, strand_map_start, strand_vertex_start, fiber_start;
	BKE_hair_get_texture_buffer_size(cache, &size, &strand_map_start, &strand_vertex_start, &fiber_start);
	
	HairStrandMapTextureBuffer *strand_map = (HairStrandMapTextureBuffer*)((char*)buffer + strand_map_start);
	HairStrandVertexTextureBuffer *strand_verts = (HairStrandVertexTextureBuffer*)((char*)buffer + strand_vertex_start);
	HairFiberTextureBuffer *fibers = (HairFiberTextureBuffer*)((char*)buffer + fiber_start);
	
	hair_get_strand_buffer(
	            cache,
	            strand_map,
	            strand_verts);
	hair_get_fiber_buffer(
	            cache,
	            fibers);
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

/* === Fiber Curve Interpolation === */

/* NOTE: Keep this code in sync with the GLSL version!
 * see common_hair_fibers_lib.glsl
 */

/* Subdivide a curve */
static int hair_curve_subdivide(
        const HairFiberCurve* curve,
        const HairFiberVertex* verts,
        int subdiv,
        int vertco_stride,
        float *r_vertco)
{
	{
		/* Move vertex positions from the dense array to their initial configuration for subdivision.
		 * Also add offset to ensure the curve starts on the scalp surface.
		 */
		const int step = (1 << subdiv) * vertco_stride;
		BLI_assert(curve->numverts > 0);
		
		float *dst = r_vertco;
		for (int i = 0; i < curve->numverts; ++i) {
			copy_v3_v3(dst, verts[i].co);
			dst = POINTER_OFFSET(dst, step);
		}
	}
	
	/* Subdivide */
	for (int d = 0; d < subdiv; ++d) {
		const int num_edges = (curve->numverts - 1) << d;
		const int hstep = (1 << (subdiv - d - 1)) * vertco_stride;
		const int step = (1 << (subdiv - d)) * vertco_stride;
		
		/* Calculate edge points */
		{
			float *p = r_vertco;
			for (int k = 0; k < num_edges; ++k) {
				float *ps = POINTER_OFFSET(p, step);
				float *ph = POINTER_OFFSET(p, hstep);
				add_v3_v3v3(ph, p, ps);
				mul_v3_fl(ph, 0.5f);
				p = ps;
			}
		}
		
		/* Move original points */
		{
			float *p = r_vertco;
			for (int k = 1; k < num_edges; ++k) {
				float *ps = POINTER_OFFSET(p, step);
				float *ph = POINTER_OFFSET(p, hstep);
				float *hp = POINTER_OFFSET(p, -hstep);
				add_v3_v3v3(p, hp, ph);
				mul_v3_fl(p, 0.5f);
				p = ps;
			}
		}
	}
	
	const int num_verts = ((curve->numverts - 1) << subdiv) + 1;
	return num_verts;
}

/* === Render API === */

/* Calculate required size for render buffers. */
void BKE_hair_render_get_buffer_size(
        const HairExportCache* cache,
        int subdiv,
        int *r_totcurves,
        int *r_totverts)
{
	*r_totcurves = cache->totfollicles;
	
	const int subdiv_factor = 1 << subdiv;
	for (int i = 0; i < cache->totfollicles; ++i)
	{
		const int numverts = cache->fiber_curves[cache->follicles[i].curve].numverts;
		*r_totverts = (numverts - 1) * subdiv_factor + 1;
	}
}

/* Create render data in existing buffers.
 * Buffers must be large enough according to BKE_hair_get_render_buffer_size.
 */
void BKE_hair_render_fill_buffers(
        const HairExportCache* cache,
        int subdiv,
        int vertco_stride,
        int *r_curvestart,
        int *r_curvelen,
        float *r_vertco)
{
	int vertstart = 0;
	float *vert = r_vertco;
	for (int i = 0; i < cache->totfollicles; ++i)
	{
		const HairFiberCurve *curve = &cache->fiber_curves[cache->follicles[i].curve];
		const HairFiberVertex *verts = &cache->fiber_verts[curve->vertstart];
		const int numverts = curve->numverts;
		r_curvestart[i] = vertstart;
		r_curvelen[i] = numverts;
		
		hair_curve_subdivide(curve, verts, subdiv, vertco_stride, vert);
		
		vertstart += numverts;
		vert = POINTER_OFFSET(vert, vertco_stride * numverts);
	}
}
