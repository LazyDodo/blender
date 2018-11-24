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
 * Copyright (C) 2014 by Martin Felke.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/BKE_fracture_util.h
 *  \ingroup blenkernel
 *  \brief CSG operations
 */

#ifndef BKE_FRACTURE_UTIL_H
#define BKE_FRACTURE_UTIL_H

#include "DNA_fracture_types.h"

typedef struct BisectContext {
	bool clear_inner;
	bool clear_outer;
	bool use_fill;
	bool do_fast_bisect;
	bool use_smooth_inner;

	char uv_layer[64];
	float normal[3];
	float centroid[3];
	float obmat[4][4];
	short inner_material_index;
	struct KDTree *geometry_limitation_tree;

} BisectContext;

typedef struct BooleanContext {
	short inner_material_index;
	int operation; /*0 == intersection, 2 == difference*/

	//fractal stuff
	bool use_fractal;
	bool use_smooth_inner;
	int num_cuts;
	int num_iterations;
	float fractal_amount;
	float cutter_plane_matrix[4][4];
	float cutter_plane_radius;

	char uv_layer[64];
	float thresh;
} BooleanContext;

Mesh* BKE_fracture_mesh_boolean(Mesh* geometry, Mesh* shard, Object* obj, BooleanContext *ctx);
Mesh* BKE_fracture_mesh_bisect(Mesh* geometry, Shard *raw_shard, BisectContext* ctx);
void BKE_fracture_mesh_boolean_fractal(Mesh* geometry, Mesh **outputA, Mesh** outputB, Object *obj, BooleanContext *ctx);
void BKE_fracture_mesh_bisect_fast(Mesh* geometry, Mesh **outputA, Mesh** outputB, BisectContext *ctx);

#endif /* BKE_FRACTURE_UTIL_H*/
