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

/** \file draw_cache_impl_groom.c
 *  \ingroup draw
 *
 * \brief Groom API for render engines
 */

#include <string.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"
#include "BLI_listbase.h"
#include "BLI_math.h"

#include "DNA_groom_types.h"
#include "DNA_scene_types.h"
#include "DNA_userdef_types.h"

#include "BKE_groom.h"
#include "BKE_mesh.h"

#include "GPU_batch.h"

#include "draw_cache_impl.h"  /* own include */

#define SELECT   1

static void groom_batch_cache_clear(Groom *groom);

/* ---------------------------------------------------------------------- */
/* Groom Gwn_Batch Cache */

typedef struct GroomBatchCache {
	Gwn_VertBuf *pos;
	Gwn_IndexBuf *edges;
	Gwn_IndexBuf *faces;

	Gwn_Batch *all_verts;
	Gwn_Batch *all_edges;
	Gwn_Batch *all_triangles;
	Gwn_Batch **shaded_triangles;
	int mat_len;

	Gwn_Batch *overlay_verts;

	/* settings to determine if cache is invalid */
	bool is_dirty;

	bool is_editmode;
} GroomBatchCache;

/* Gwn_Batch cache management. */

static bool groom_batch_cache_valid(Groom *groom)
{
	GroomBatchCache *cache = groom->batch_cache;

	if (cache == NULL) {
		return false;
	}

	if (cache->is_editmode != (groom->editgroom != NULL)) {
		return false;
	}

	if (cache->is_dirty == false) {
		return true;
	}
	else {
		if (cache->is_editmode) {
			return false;
		}
	}

	return true;
}

static void groom_batch_cache_init(Groom *groom)
{
	GroomBatchCache *cache = groom->batch_cache;

	if (!cache) {
		cache = groom->batch_cache = MEM_callocN(sizeof(*cache), __func__);
	}
	else {
		memset(cache, 0, sizeof(*cache));
	}

	cache->is_editmode = (groom->editgroom != NULL);

	cache->is_dirty = false;
}

static GroomBatchCache *groom_batch_cache_get(Groom *groom)
{
	if (!groom_batch_cache_valid(groom)) {
		groom_batch_cache_clear(groom);
		groom_batch_cache_init(groom);
	}
	return groom->batch_cache;
}

void DRW_groom_batch_cache_dirty(Groom *groom, int mode)
{
	GroomBatchCache *cache = groom->batch_cache;
	if (cache == NULL) {
		return;
	}
	switch (mode) {
		case BKE_GROOM_BATCH_DIRTY_ALL:
			cache->is_dirty = true;
			break;
		case BKE_GROOM_BATCH_DIRTY_SELECT:
			/* TODO Separate Flag vbo */
			GWN_BATCH_DISCARD_SAFE(cache->overlay_verts);
			break;
		default:
			BLI_assert(0);
	}
}

static void groom_batch_cache_clear(Groom *groom)
{
	GroomBatchCache *cache = groom->batch_cache;
	if (!cache) {
		return;
	}

	GWN_BATCH_DISCARD_SAFE(cache->all_verts);
	GWN_BATCH_DISCARD_SAFE(cache->all_edges);
	GWN_BATCH_DISCARD_SAFE(cache->all_triangles);
	GWN_BATCH_DISCARD_SAFE(cache->overlay_verts);
	for (int i = 0; i < cache->mat_len; ++i)
	{
		GWN_BATCH_DISCARD_SAFE(cache->shaded_triangles[i]);
	}
	MEM_SAFE_FREE(cache->shaded_triangles);
	cache->mat_len = 0;

	GWN_VERTBUF_DISCARD_SAFE(cache->pos);
	GWN_INDEXBUF_DISCARD_SAFE(cache->edges);
	GWN_INDEXBUF_DISCARD_SAFE(cache->faces);
}

void DRW_groom_batch_cache_free(Groom *groom)
{
	groom_batch_cache_clear(groom);
	MEM_SAFE_FREE(groom->batch_cache);
}

enum {
	VFLAG_VERTEX_SELECTED = 1 << 0,
	VFLAG_VERTEX_ACTIVE   = 1 << 1,
};

BLI_INLINE char make_vertex_flag(bool active, bool selected)
{
	char vflag = 0;
	if (active)
	{
		vflag |= VFLAG_VERTEX_ACTIVE;
	}
	if (selected)
	{
		vflag |= VFLAG_VERTEX_SELECTED;
	}
	return vflag;
}

/* Parts of the groom object to render */
typedef enum GroomRenderPart
{
	GM_RENDER_REGIONS           = (1 << 0),   /* Draw scalp regions */
	GM_RENDER_CENTER_CURVES     = (1 << 1),   /* Draw center curves of bundles */
	GM_RENDER_SECTIONS          = (1 << 2),   /* Draw section curves */
	GM_RENDER_HULL              = (1 << 3),   /* Draw hull faces */
	
	GM_RENDER_ALL       = GM_RENDER_REGIONS | GM_RENDER_CENTER_CURVES | GM_RENDER_SECTIONS | GM_RENDER_HULL,
} GroomRenderPart;

static int groom_count_verts_regions(const ListBase *regions)
{
	// TODO
	UNUSED_VARS(regions);
	return 0;
}

static int groom_count_verts_center_curves(const ListBase *regions, bool use_curve_cache)
{
	int vert_len = 0;
	for (GroomRegion *region = regions->first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		if (use_curve_cache)
		{
			vert_len += bundle->curvesize;
		}
		else
		{
			vert_len += bundle->totsections;
		}
	}
	return vert_len;
}

static int groom_count_verts_sections(const ListBase *regions)
{
	int vert_len = 0;
	for (GroomRegion *region = regions->first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		vert_len += bundle->totverts;
	}
	return vert_len;
}

static int groom_count_verts_hull(const ListBase *regions, bool use_curve_cache)
{
	int vert_len = 0;
	for (GroomRegion *region = regions->first; region; region = region->next)
	{
		GroomBundle *bundle = &region->bundle;
		if (use_curve_cache)
		{
			vert_len += bundle->curvesize * region->numverts;
		}
		else
		{
			vert_len += bundle->totverts;
		}
	}
	return vert_len;
}

static int groom_count_verts(const ListBase *regions, int parts, bool use_curve_cache)
{
	int vert_len = 0;
	
	if (parts & GM_RENDER_REGIONS)
	{
		vert_len += groom_count_verts_regions(regions);
	}
	if (parts & GM_RENDER_CENTER_CURVES)
	{
		vert_len += groom_count_verts_center_curves(regions, use_curve_cache);
	}
	if (parts & GM_RENDER_SECTIONS)
	{
		vert_len += groom_count_verts_sections(regions);
	}
	if (parts & GM_RENDER_HULL)
	{
		vert_len += groom_count_verts_hull(regions, use_curve_cache);
	}
	
	return vert_len;
}

static int groom_count_edges(const ListBase *regions, int parts, bool use_curve_cache)
{
	int edge_len = 0;
	
	if (parts & GM_RENDER_REGIONS)
	{
		// TODO
	}
	if (parts & GM_RENDER_CENTER_CURVES)
	{
		for (GroomRegion *region = regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			if (use_curve_cache)
			{
				if (bundle->curvesize > 0)
				{
					edge_len += bundle->curvesize - 1;
				}
			}
			else
			{
				if (bundle->totsections > 0)
				{
					edge_len += bundle->totsections - 1;
				}
			}
		}
	}
	if (parts & GM_RENDER_SECTIONS)
	{
		for (GroomRegion *region = regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			// Closed edge loop, 1 edge per vertex
			edge_len += bundle->totverts;
		}
	}
	
	return edge_len;
}

static int groom_count_faces(const ListBase *regions, int parts, bool use_curve_cache)
{
	int face_len = 0;
	
	if (parts & GM_RENDER_REGIONS)
	{
		// TODO
	}
	if (parts & GM_RENDER_SECTIONS)
	{
		for (GroomRegion *region = regions->first; region; region = region->next)
		{
			/* Polygon on each section, except the first */
			GroomBundle *bundle = &region->bundle;
			if (region->numverts > 2)
			{
				int numpolys = bundle->totsections - 1;
				face_len += poly_to_tri_count(numpolys, region->numverts * numpolys);
			}
		}
	}
	if (parts & GM_RENDER_HULL)
	{
		for (GroomRegion *region = regions->first; region; region = region->next)
		{
			/* Closed triangle strip around each section */
			GroomBundle *bundle = &region->bundle;
			if (use_curve_cache)
			{
				if (bundle->curvesize > 0)
				{
					face_len += (bundle->curvesize - 1) * region->numverts * 2;
				}
			}
			else
			{
				if (bundle->totsections > 0)
				{
					face_len += (bundle->totsections - 1) * region->numverts * 2;
				}
			}
		}
	}
	
	return face_len;
}

typedef struct GroomRenderData
{
	const struct ListBase *regions;
	int curve_res;
	
	int tri_len;                    /* Total mlooptri array length */
	struct MLoopTri *mlooptri;      /* Triangulation data for sections */
} GroomRenderData;

static GroomRenderData* groom_render_data_create(Groom *groom)
{
	GroomRenderData *rdata = MEM_callocN(sizeof(GroomRenderData), "groom render data");
	EditGroom *edit = groom->editgroom;
	const ListBase *regions = edit ? &edit->regions : &groom->regions;
	
	rdata->regions = regions;
	rdata->curve_res = groom->curve_res;
	
	{
		int totpoly = 0;
		int totvert = 0;
		for (GroomRegion *region = regions->first; region; region = region->next)
		{
			/* Polygon on each section */
			GroomBundle *bundle = &region->bundle;
			if (region->numverts > 2)
			{
				int numpolys = bundle->totsections;
				totpoly += numpolys;
				totvert += region->numverts * numpolys;
				
				/* Polygons are unconnected, no shared vertices,
				 * same vertex number for each section polygon.
				 */
				int section_tri_len = poly_to_tri_count(1, region->numverts);
				rdata->tri_len += section_tri_len * numpolys;
			}
		}
		
		/* Temporary mesh buffers for calculating tesselation */
		MPoly *mpolys = MEM_mallocN(sizeof(*mpolys) * totpoly, __func__);
		MVert *mverts = MEM_mallocN(sizeof(*mverts) * totvert, __func__);
		MLoop *mloops = MEM_mallocN(sizeof(*mloops) * totvert, __func__);

		{
			int ipoly = 0;
			int ivert = 0;
			for (GroomRegion *region = regions->first; region; region = region->next)
			{
				/* Polygon on each section, except the first */
				GroomBundle *bundle = &region->bundle;
				if (region->numverts > 2)
				{
					const GroomSectionVertex *vert = bundle->verts;
					for (int i = 0; i < bundle->totsections; ++i)
					{
						mpolys[ipoly].loopstart = ivert;
						mpolys[ipoly].totloop = region->numverts;
						
						for (int j = 0; j < region->numverts; ++j)
						{
							/* No need to calculate final 3D positions,
							 * because the tesselation in local 2D space is just the same.
							 */
							copy_v2_v2(mverts[ivert].co, vert->co);
							mverts[ivert].co[2] = 0.0f;
							mloops[ivert].v = ivert;
							
							++vert;
							++ivert;
						}
						
						++ipoly;
					}
				}
			}
		}
		
		rdata->mlooptri = MEM_mallocN(sizeof(*rdata->mlooptri) * rdata->tri_len, __func__);
		BKE_mesh_recalc_looptri(mloops, mpolys, mverts, totvert, totpoly, rdata->mlooptri);
		
		MEM_freeN(mverts);
		MEM_freeN(mpolys);
		MEM_freeN(mloops);
	}
	
	return rdata;
}

static void groom_render_data_free(GroomRenderData *rdata)
{
	if (rdata->mlooptri)
	{
		MEM_freeN(rdata->mlooptri);
	}
	
	MEM_freeN(rdata);
}

#define GM_ATTR_ID_UNUSED 0xFFFFFFFF

static void groom_get_verts(
        GroomRenderData *rdata,
        int parts,
        bool use_curve_cache,
        Gwn_VertBuf *vbo,
        uint id_pos,
        uint id_nor,
        uint id_flag)
{
	int vert_len = groom_count_verts(rdata->regions, parts, use_curve_cache);
	
	GWN_vertbuf_data_alloc(vbo, vert_len);
	
	uint idx = 0;
	if (parts & GM_RENDER_REGIONS)
	{
		// TODO
	}
	if (parts & GM_RENDER_CENTER_CURVES)
	{
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			if (use_curve_cache)
			{
				/* curvecache has numverts+1 curves, the last one is the center curve */
				GroomCurveCache *cache = bundle->curvecache + bundle->curvesize * region->numverts;
				for (int i = 0; i < bundle->curvesize; ++i, ++cache)
				{
					if (id_pos != GM_ATTR_ID_UNUSED)
					{
						GWN_vertbuf_attr_set(vbo, id_pos, idx, cache->co);
					}
					if (id_flag != GM_ATTR_ID_UNUSED)
					{
						char vflag = make_vertex_flag(false, false);
						GWN_vertbuf_attr_set(vbo, id_flag, idx, &vflag);
					}
					
					++idx;
				}
			}
			else
			{
				const bool active = region->flag & GM_REGION_SELECT;
				GroomSection *section = bundle->sections;
				for (int i = 0; i < bundle->totsections; ++i, ++section)
				{
					if (id_pos != GM_ATTR_ID_UNUSED)
					{
						GWN_vertbuf_attr_set(vbo, id_pos, idx, section->center);
					}
					if (id_flag != GM_ATTR_ID_UNUSED)
					{
						char vflag = make_vertex_flag(active, section->flag & GM_SECTION_SELECT);
						GWN_vertbuf_attr_set(vbo, id_flag, idx, &vflag);
					}
					
					++idx;
				}
			}
		}
	}
	if (parts & GM_RENDER_SECTIONS)
	{
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			GroomSection *section = bundle->sections;
			GroomSectionVertex *vertex = bundle->verts;
			for (int i = 0; i < bundle->totsections; ++i, ++section)
			{
				const bool active = (region->flag & GM_REGION_SELECT) && (section->flag & GM_SECTION_SELECT);
				
				for (int j = 0; j < region->numverts; ++j, ++vertex)
				{
					if (id_pos != GM_ATTR_ID_UNUSED)
					{
						float co[3] = {vertex->co[0], vertex->co[1], 0.0f};
						mul_m3_v3(section->mat, co);
						add_v3_v3(co, section->center);
						GWN_vertbuf_attr_set(vbo, id_pos, idx, co);
					}
					if (id_nor != GM_ATTR_ID_UNUSED)
					{
						GWN_vertbuf_attr_set(vbo, id_nor, idx, section->mat[2]);
					}
					if (id_flag != GM_ATTR_ID_UNUSED)
					{
						char vflag = make_vertex_flag(active, vertex->flag & GM_VERTEX_SELECT);
						GWN_vertbuf_attr_set(vbo, id_flag, idx, &vflag);
					}
					
					++idx;
				}
			}
		}
	}
	if (parts & GM_RENDER_HULL)
	{
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			if (use_curve_cache)
			{
				
				GroomCurveCache *cache = bundle->curvecache;
				for (int i = 0; i < region->numverts; ++i)
				{
					for (int j = 0; j < bundle->curvesize; ++j, ++cache)
					{
						if (id_pos != GM_ATTR_ID_UNUSED)
						{
							GWN_vertbuf_attr_set(vbo, id_pos, idx, cache->co);
						}
						if (id_flag != GM_ATTR_ID_UNUSED)
						{
							char vflag = make_vertex_flag(false, false);
							GWN_vertbuf_attr_set(vbo, id_flag, idx, &vflag);
						}
						
						++idx;
					}
				}
			}
			else
			{
				GroomSection *section = bundle->sections;
				GroomSectionVertex *vertex = bundle->verts;
				for (int i = 0; i < bundle->totsections; ++i, ++section)
				{
					for (int j = 0; j < region->numverts; ++j, ++vertex)
					{
						if (id_pos != GM_ATTR_ID_UNUSED)
						{
							float co[3] = {vertex->co[0], vertex->co[1], 0.0f};
							mul_m3_v3(section->mat, co);
							add_v3_v3(co, section->center);
							GWN_vertbuf_attr_set(vbo, id_pos, idx, co);
						}
						if (id_flag != GM_ATTR_ID_UNUSED)
						{
							char vflag = make_vertex_flag(false, false);
							GWN_vertbuf_attr_set(vbo, id_flag, idx, &vflag);
						}
						
						++idx;
					}
				}
			}
		}
	}
}

static void groom_get_edges(
        GroomRenderData *rdata,
        int parts,
        bool use_curve_cache,
        Gwn_IndexBuf **ibo)
{
	Gwn_IndexBufBuilder elb;
	
	int vert_len = groom_count_verts(rdata->regions, parts, use_curve_cache);
	int edge_len = groom_count_edges(rdata->regions, parts, use_curve_cache);
	
	GWN_indexbuf_init(&elb, GWN_PRIM_LINES, edge_len, vert_len);
	
	uint idx = 0;
	if (parts & GM_RENDER_REGIONS)
	{
		// TODO
		idx += groom_count_verts_regions(rdata->regions);
	}
	if (parts & GM_RENDER_CENTER_CURVES)
	{
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			if (use_curve_cache)
			{
				for (int i = 0; i < bundle->curvesize - 1; ++i)
				{
					GWN_indexbuf_add_line_verts(&elb, idx + i, idx + i + 1);
				}
				idx += bundle->curvesize;
			}
			else
			{
				GroomSection *section = bundle->sections;
				for (int i = 0; i < bundle->totsections - 1; ++i, ++section)
				{
					GWN_indexbuf_add_line_verts(&elb, idx + i, idx + i + 1);
				}
				idx += bundle->totsections;
			}
		}
	}
	if (parts & GM_RENDER_SECTIONS)
	{
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			const int numshapeverts = region->numverts;
			if (numshapeverts > 1)
			{
				/* Skip the first section */
				for (int i = 1; i < bundle->totsections; ++i)
				{
					uint idx0 = idx + i * numshapeverts;
					for (int j = 0; j < numshapeverts - 1; ++j)
					{
						GWN_indexbuf_add_line_verts(&elb, idx0 + j, idx0 + j + 1);
					}
					// close the loop
					GWN_indexbuf_add_line_verts(&elb, idx0 + (numshapeverts-1), idx0);
				}
				
				idx += bundle->totverts;
			}
		}
	}
	if (parts & GM_RENDER_HULL)
	{
		idx += groom_count_verts_hull(rdata->regions, use_curve_cache);
	}
	
	*ibo = GWN_indexbuf_build(&elb);
}

static void groom_get_faces(
        GroomRenderData *rdata,
        int parts,
        bool use_curve_cache,
        Gwn_IndexBuf **ibo)
{
	Gwn_IndexBufBuilder elb;
	
	int vert_len = groom_count_verts(rdata->regions, parts, use_curve_cache);
	int face_len = groom_count_faces(rdata->regions, parts, use_curve_cache);
	
	GWN_indexbuf_init_ex(&elb, GWN_PRIM_TRIS, face_len, vert_len, true);
	
	uint idx = 0;
	if (parts & GM_RENDER_REGIONS)
	{
		// TODO
		idx += groom_count_verts_regions(rdata->regions);
	}
	if (parts & GM_RENDER_CENTER_CURVES)
	{
		idx += groom_count_verts_center_curves(rdata->regions, use_curve_cache);
	}
	if (parts & GM_RENDER_SECTIONS)
	{
		const MLoopTri *mtri = rdata->mlooptri;
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
			GroomBundle *bundle = &region->bundle;
			const int numshapeverts = region->numverts;
			if (numshapeverts > 1)
			{
				int section_tri_len = poly_to_tri_count(1, region->numverts);
				/* Skip the first section */
				mtri += section_tri_len;
				for (int i = 1; i < bundle->totsections; ++i)
				{
					for (int j = 0; j < section_tri_len; ++j, ++mtri)
					{
						GWN_indexbuf_add_tri_verts(
						            &elb,
						            idx + mtri->tri[0],
						            idx + mtri->tri[1],
						            idx + mtri->tri[2]);
					}
				}
				
			}
		}
		idx += groom_count_verts_sections(rdata->regions);
	}
	if (parts & GM_RENDER_HULL)
	{
		for (GroomRegion *region = rdata->regions->first; region; region = region->next)
		{
#if 0
			/* Closed triangle strip around each section */
			GroomBundle *bundle = &region->bundle;
			if (use_curve_cache)
			{
				GroomCurveCache *cache = bundle->curvecache;
				const int numshapeverts = region->numverts;
				if (numshapeverts > 1)
				{
					for (int i = 0; i < numshapeverts; ++i)
					{
						for (int j = 0; j < bundle->curvesize - 1; ++j, ++cache)
						{
							GWN_indexbuf_add_tri_verts(&elb, idx + j, idx + curvesize + j, idx + j + 1);
						}
					}
					idx += bundle->curvesize * numshapeverts;
				}
			}
			else
			{
				GroomSection *section = bundle->sections;
				for (int i = 0; i < bundle->totsections - 1; ++i, ++section)
				{
					GWN_indexbuf_add_line_verts(&elb, idx + i, idx + i + 1);
				}
				idx += bundle->totsections;
			}
#endif
		}
	}
	
	*ibo = GWN_indexbuf_build(&elb);
}

/* Gwn_Batch cache usage. */
static Gwn_VertBuf *groom_batch_cache_get_pos(GroomRenderData *rdata, GroomBatchCache *cache, int parts)
{
	if (cache->pos == NULL) {
		static Gwn_VertFormat format = { 0 };
		static struct { uint pos, nor; } attr_id;
		
		GWN_vertformat_clear(&format);
		
		/* initialize vertex format */
		attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
		attr_id.nor = GWN_vertformat_attr_add(&format, "nor", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
		
		cache->pos = GWN_vertbuf_create_with_format(&format);
		
		groom_get_verts(rdata, parts, true, cache->pos, attr_id.pos, attr_id.nor, GM_ATTR_ID_UNUSED);
	}

	return cache->pos;
}

static Gwn_IndexBuf *groom_batch_cache_get_edges(GroomRenderData *rdata, GroomBatchCache *cache, int parts)
{
	if (cache->edges == NULL) {
		groom_get_edges(rdata, parts, true, &cache->edges);
	}

	return cache->edges;
}

static Gwn_IndexBuf *groom_batch_cache_get_faces(GroomRenderData *rdata, GroomBatchCache *cache, int parts)
{
	if (cache->faces == NULL) {
		groom_get_faces(rdata, parts, true, &cache->faces);
	}

	return cache->faces;
}

static void groom_batch_cache_create_overlay_batches(GroomRenderData *rdata, GroomBatchCache *cache, int parts)
{
	if (cache->overlay_verts == NULL) {
		static Gwn_VertFormat format = { 0 };
		static struct { uint pos, nor, data; } attr_id;
		if (format.attrib_ct == 0) {
			/* initialize vertex format */
			attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
			attr_id.nor = GWN_vertformat_attr_add(&format, "nor", GWN_COMP_F32, 3, GWN_FETCH_FLOAT);
			attr_id.data = GWN_vertformat_attr_add(&format, "data", GWN_COMP_U8, 1, GWN_FETCH_INT);
		}
		
		Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);
		
		groom_get_verts(rdata, parts, false, vbo, attr_id.pos, attr_id.nor, attr_id.data);
		
		cache->overlay_verts = GWN_batch_create_ex(GWN_PRIM_POINTS, vbo, NULL, GWN_BATCH_OWNS_VBO);
	}	
}

Gwn_Batch *DRW_groom_batch_cache_get_all_edges(Groom *groom)
{
	GroomBatchCache *cache = groom_batch_cache_get(groom);

	if (cache->all_edges == NULL) {
		GroomRenderData *rdata = groom_render_data_create(groom);
		
		cache->all_edges = GWN_batch_create(
		                       GWN_PRIM_LINES,
		                       groom_batch_cache_get_pos(rdata, cache, GM_RENDER_ALL),
		                       groom_batch_cache_get_edges(rdata, cache, GM_RENDER_ALL));
		
		groom_render_data_free(rdata);
	}

	return cache->all_edges;
}

Gwn_Batch *DRW_groom_batch_cache_get_all_triangles(Groom *groom)
{
	GroomBatchCache *cache = groom_batch_cache_get(groom);

	if (cache->all_triangles == NULL) {
		GroomRenderData *rdata = groom_render_data_create(groom);
		
		cache->all_triangles = GWN_batch_create(
		                       GWN_PRIM_TRIS,
		                       groom_batch_cache_get_pos(rdata, cache, GM_RENDER_ALL),
		                       groom_batch_cache_get_faces(rdata, cache, GM_RENDER_ALL));
		
		groom_render_data_free(rdata);
	}

	return cache->all_triangles;
}

Gwn_Batch **DRW_groom_batch_cache_get_surface_shaded(Groom *groom, struct GPUMaterial **UNUSED(gpumat_array), uint gpumat_array_len)
{
	GroomBatchCache *cache = groom_batch_cache_get(groom);

	if (cache->shaded_triangles == NULL) {
		cache->mat_len = gpumat_array_len;
		cache->shaded_triangles = MEM_callocN(sizeof(*cache->shaded_triangles) * cache->mat_len, __func__);
		cache->shaded_triangles[0] = DRW_groom_batch_cache_get_all_triangles(groom);
		for (int i = 1; i < cache->mat_len; ++i) {
			cache->shaded_triangles[i] = NULL;
		}
	}
	return cache->shaded_triangles;
}

Gwn_Batch *DRW_groom_batch_cache_get_all_verts(Groom *groom)
{
	GroomBatchCache *cache = groom_batch_cache_get(groom);

	if (cache->all_verts == NULL) {
		GroomRenderData *rdata = groom_render_data_create(groom);
		
		cache->all_verts = GWN_batch_create(
		                       GWN_PRIM_POINTS,
		                       groom_batch_cache_get_pos(rdata, cache, GM_RENDER_ALL),
		                       NULL);
		
		groom_render_data_free(rdata);
	}

	return cache->all_verts;
}

Gwn_Batch *DRW_groom_batch_cache_get_overlay_verts(Groom *groom, int mode)
{
	GroomBatchCache *cache = groom_batch_cache_get(groom);

	if (cache->overlay_verts == NULL) {
		GroomRenderData *rdata = groom_render_data_create(groom);
		
		GroomRenderPart parts = 0;
		switch ((GroomEditMode)mode)
		{
			case GM_EDIT_MODE_REGIONS: parts |= GM_RENDER_REGIONS; break;
			case GM_EDIT_MODE_CURVES: parts |= GM_RENDER_CENTER_CURVES; break;
			case GM_EDIT_MODE_SECTIONS: parts |= GM_RENDER_SECTIONS; break;
		}
		
		groom_batch_cache_create_overlay_batches(rdata, cache, parts);
		
		groom_render_data_free(rdata);
	}

	return cache->overlay_verts;
}
