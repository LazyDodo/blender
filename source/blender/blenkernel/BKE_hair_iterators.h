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

#ifndef __BKE_HAIR_ITERATORS_H__
#define __BKE_HAIR_ITERATORS_H__

/** \file blender/blenkernel/BKE_hair_iterators.h
 *  \ingroup bke
 */

#include "BLI_utildefines.h"

#include "DNA_hair_types.h"

#include "BKE_hair.h"

/* === Hair Element Iterators === */

typedef struct HairIter__curves
{
	HairFiberCurve *curve;
	HairFiberCurve *last;
} HairIter__curves;

typedef struct HairIter__verts
{
	HairFiberVertex *vertex;
	HairFiberVertex *last;
} HairIter__verts;

typedef struct HairIter__follicles
{
	HairFollicle *follicle;
	HairFollicle *last;
} HairIter__follicles;

typedef struct HairIter__follicle_curves
{
	HairFollicle *follicle;
	HairFollicle *last;
	HairFiberCurve *curves_array;
	HairFiberCurve *curve;
} HairIter__follicle_curves;

typedef struct HairIter__curve_verts
{
	HairFiberVertex *vertex;
	HairFiberVertex *last;
} HairIter__curve_verts;

typedef struct HairIterator
{
	union {
		HairIter__curves curves;
		HairIter__verts verts;
		HairIter__follicles follicles;
		HairIter__follicle_curves follicle_curves;
		HairIter__curve_verts curve_verts;
	} data;
} HairIterator;

/* === Inline functions for element iterators === */

/* Curves */

BLI_INLINE HairFiberCurve* BKE_hair_iter__curves_init(HairIterator *iter, const HairCurveData *curve_data)
{
	iter->data.curves.last = curve_data->curves + curve_data->totcurves;
	iter->data.curves.curve = curve_data->curves;
	return iter->data.curves.curve;
}

BLI_INLINE bool BKE_hair_iter__curves_valid(const HairIterator *iter)
{
	return iter->data.curves.curve != iter->data.curves.last;
}

BLI_INLINE HairFiberCurve* BKE_hair_iter__curves_next(HairIterator *iter)
{
	++iter->data.curves.curve;
	return iter->data.curves.curve;
}

/* Vertices */

BLI_INLINE HairFiberVertex* BKE_hair_iter__verts_init(HairIterator *iter, const HairCurveData *curve_data)
{
	iter->data.verts.last = curve_data->verts + curve_data->totverts;
	iter->data.verts.vertex = curve_data->verts;
	return iter->data.verts.vertex;
}

BLI_INLINE bool BKE_hair_iter__verts_valid(const HairIterator *iter)
{
	return iter->data.verts.vertex != iter->data.verts.last;
}

BLI_INLINE HairFiberVertex* BKE_hair_iter__verts_next(HairIterator *iter)
{
	++iter->data.verts.vertex;
	return iter->data.verts.vertex;
}

/* Follicles */

BLI_INLINE HairFollicle* BKE_hair_iter__follicles_init(HairIterator *iter, const HairCurveData *curve_data)
{
	iter->data.follicles.last = curve_data->follicles + curve_data->totfollicles;
	iter->data.follicles.follicle = curve_data->follicles;
	return iter->data.follicles.follicle;
}

BLI_INLINE bool BKE_hair_iter__follicles_valid(const HairIterator *iter)
{
	return iter->data.follicles.follicle != iter->data.follicles.last;
}

BLI_INLINE HairFollicle* BKE_hair_iter__follicles_next(HairIterator *iter)
{
	++iter->data.follicles.follicle;
	return iter->data.follicles.follicle;
}

/* Follicle Curves */

BLI_INLINE HairFollicle* BKE_hair_iter__follicle_curves_init(HairIterator *iter, const HairCurveData *curve_data)
{
	HairIter__follicle_curves *intern = &iter->data.follicle_curves;

	intern->last = curve_data->follicles + curve_data->totfollicles;
	intern->curves_array = curve_data->curves;

	intern->follicle = curve_data->follicles;
	while (intern->follicle != intern->last) {
		if (intern->follicle->curve != HAIR_CURVE_INDEX_NONE) {
			intern->curve = &intern->curves_array[intern->follicle->curve];
			break;
		}
		++intern->follicle;
	}
	return intern->follicle;
}

BLI_INLINE bool BKE_hair_iter__follicle_curves_valid(const HairIterator *iter)
{
	const HairIter__follicle_curves *intern = &iter->data.follicle_curves;

	return intern->follicle != intern->last;
}

BLI_INLINE HairFollicle* BKE_hair_iter__follicle_curves_next(HairIterator *iter)
{
	HairIter__follicle_curves *intern = &iter->data.follicle_curves;

	++intern->follicle;
	while (intern->follicle != intern->last) {
		if (intern->follicle->curve != HAIR_CURVE_INDEX_NONE) {
			intern->curve = &intern->curves_array[intern->follicle->curve];
			break;
		}
		++intern->follicle;
	}
	return intern->follicle;
}

/* Curve Vertices */

BLI_INLINE HairFiberVertex* BKE_hair_iter__curve_verts_init(HairIterator *iter, const HairCurveData *curve_data, const HairFiberCurve *curve)
{
	HairIter__curve_verts *intern = &iter->data.curve_verts;

	intern->last = &curve_data->verts[curve->vertstart + curve->numverts];
	intern->vertex = &curve_data->verts[curve->vertstart];
	return intern->vertex;
}

BLI_INLINE bool BKE_hair_iter__curve_verts_valid(const HairIterator *iter)
{
	const HairIter__curve_verts *intern = &iter->data.curve_verts;

	return intern->vertex != intern->last;
}

BLI_INLINE HairFiberVertex* BKE_hair_iter__curve_verts_next(HairIterator *iter)
{
	HairIter__curve_verts *intern = &iter->data.curve_verts;

	++intern->vertex;
	return intern->vertex;
}

/* Loop macros */

#define BKE_HAIR_ITER_CURVES(ele, iter, curve_data) \
	for ((ele) = BKE_hair_iter__curves_init((iter), (curve_data)); \
	     BKE_hair_iter__curves_valid(iter); \
	     (ele) = BKE_hair_iter__curves_next(iter))

#define BKE_HAIR_ITER_CURVES_INDEX(ele, iter, curve_data, indexvar) \
	for ((ele) = BKE_hair_iter__curves_init((iter), (curve_data)), (indexvar = 0); \
	     BKE_hair_iter__curves_valid(iter); \
	     (ele) = BKE_hair_iter__curves_next(iter), ++(indexvar))

#define BKE_HAIR_ITER_VERTS(ele, iter, curve_data) \
	for ((ele) = BKE_hair_iter__verts_init((iter), (curve_data)); \
	     BKE_hair_iter__verts_valid(iter); \
	     (ele) = BKE_hair_iter__verts_next(iter))

#define BKE_HAIR_ITER_VERTS_INDEX(ele, iter, curve_data, indexvar) \
	for ((ele) = BKE_hair_iter__verts_init((iter), (curve_data)), (indexvar = 0); \
	     BKE_hair_iter__verts_valid(iter); \
	     (ele) = BKE_hair_iter__verts_next(iter), ++(indexvar))

#define BKE_HAIR_ITER_FOLLICLES(ele, iter, curve_data) \
	for ((ele) = BKE_hair_iter__follicles_init((iter), (curve_data)); \
	     BKE_hair_iter__follicles_valid(iter); \
	     (ele) = BKE_hair_iter__follicles_next(iter))

#define BKE_HAIR_ITER_FOLLICLES_INDEX(ele, iter, curve_data, indexvar) \
	for ((ele) = BKE_hair_iter__follicles_init((iter), (curve_data)), (indexvar = 0); \
	     BKE_hair_iter__follicles_valid(iter); \
	     (ele) = BKE_hair_iter__follicles_next(iter), ++(indexvar))

#define BKE_HAIR_ITER_FOLLICLE_CURVES(ele, curvevar, iter, curve_data) \
	for ((ele) = BKE_hair_iter__follicle_curves_init((iter), (curve_data)), \
	         (curvevar) = (iter)->data.follicle_curves.curve; \
	     BKE_hair_iter__follicle_curves_valid(iter); \
	     (ele) = BKE_hair_iter__follicle_curves_next((iter)), \
	         (curvevar) = (iter)->data.follicle_curves.curve)

#define BKE_HAIR_ITER_CURVE_VERTS(ele, iter, curve_data, curve) \
	for ((ele) = BKE_hair_iter__curve_verts_init((iter), (curve_data), (curve)); \
	     BKE_hair_iter__curve_verts_valid(iter); \
	     (ele) = BKE_hair_iter__curve_verts_next(iter))

#define BKE_HAIR_ITER_CURVE_VERTS_INDEX(ele, iter, curve_data, curve, indexvar) \
	for ((ele) = BKE_hair_iter__curve_verts_init((iter), (curve_data), (curve)), (indexvar = 0); \
	     BKE_hair_iter__curve_verts_valid(iter); \
	     (ele) = BKE_hair_iter__curve_verts_next(iter), ++(indexvar))

#endif
