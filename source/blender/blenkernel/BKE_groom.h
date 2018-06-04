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

#ifndef __BKE_GROOM_H__
#define __BKE_GROOM_H__

/** \file BKE_groom.h
 *  \ingroup bke
 */

struct BoundBox;
struct Depsgraph;
struct Groom;
struct GroomBundle;
struct GroomRegion;
struct Main;
struct Mesh;
struct Object;
struct Scene;

void BKE_groom_init(struct Groom *groom);
void *BKE_groom_add(struct Main *bmain, const char *name);

void BKE_groom_free(struct Groom *groom);
void BKE_groom_bundle_curve_cache_clear(struct GroomBundle *bundle);

void BKE_groom_copy_data(struct Main *bmain, struct Groom *groom_dst, const struct Groom *groom_src, const int flag);
struct Groom *BKE_groom_copy(struct Main *bmain, const struct Groom *groom);

void BKE_groom_make_local(struct Main *bmain, struct Groom *groom, const bool lib_local);

bool BKE_groom_minmax(struct Groom *groom, float min[3], float max[3]);
void BKE_groom_boundbox_calc(struct Groom *groom);
struct BoundBox *BKE_groom_boundbox_get(struct Object *ob);


/* === Curve cache === */

void BKE_groom_curve_cache_update(struct Groom *groom, const struct Mesh *scalp);
void BKE_groom_curve_cache_clear(struct Groom *groom);


/* === Scalp regions === */

/* Try to bind bundles to their scalp regions */
void BKE_groom_bind_scalp_regions(struct Groom *groom, bool force_rebind);

bool BKE_groom_region_bind(struct Groom *groom, struct GroomRegion *region, bool force_rebind);
void BKE_groom_region_unbind(struct GroomRegion *region);


/* === Constraints === */

/* Apply constraints on groom geometry */
void BKE_groom_apply_constraints(struct Groom *groom, struct Mesh *scalp);


/* === Hair System === */

/* Create follicles on the scalp surface for hair fiber rendering */
void BKE_groom_hair_distribute(struct Groom *groom, unsigned int seed, int hair_count);

/* Calculate guide curve shapes based on groom bundle deformation */
void BKE_groom_hair_update_guide_curves(struct Groom *groom);


/* === Depsgraph evaluation === */

void BKE_groom_eval_geometry(const struct Depsgraph *depsgraph, struct Groom *groom);


/* === Draw Cache === */

enum {
	BKE_GROOM_BATCH_DIRTY_ALL = 0,
	BKE_GROOM_BATCH_DIRTY_SELECT,
};
void BKE_groom_batch_cache_dirty(struct Groom *groom, int mode);
void BKE_groom_batch_cache_free(struct Groom *groom);

/* === Iterators === */

/* Utility class for iterating over groom elements */
typedef struct GroomIterator
{
	int isection;                       /* section index */
	struct GroomSection *section;       /* section data pointer */
	
	int ivertex;                        /* vertex index */
	int isectionvertex;                 /* vertex index for the inner loop */
	struct GroomSectionVertex *vertex;  /* vertex data pointer */
} GroomIterator;

#define GROOM_ITER_SECTIONS(iter, bundle) \
	for (iter.isection = 0, iter.section = (bundle)->sections; \
	     iter.isection < (bundle)->totsections; \
	     ++iter.isection, ++iter.section)

#define GROOM_ITER_SECTION_LOOPS(iter, sectionvar, vertexvar, bundle) \
	for (iter.isection = 0, iter.section = (bundle)->sections, iter.ivertex = 0, iter.vertex = (bundle)->verts; \
	     iter.isection < (bundle)->totsections; \
	     ++iter.isection, ++iter.section) \
		for (iter.isectionvertex = 0; \
		     iter.isectionvertex < (bundle)->numloopverts; \
		     ++iter.isectionvertex, ++iter.vertex)

/* === Utility functions === */
struct Mesh* BKE_groom_get_scalp(struct Groom *groom);

#endif /*  __BKE_GROOM_H__ */
