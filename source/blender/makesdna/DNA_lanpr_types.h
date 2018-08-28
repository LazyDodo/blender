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
 * The Original Code is Copyright (C) 2010 Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __DNA_LANPR_TYPES_H__
#define __DNA_LANPR_TYPES_H__

/** \file DNA_lanpr_types.h
 *  \ingroup DNA
 */

#include "DNA_listBase.h"
#include "DNA_ID.h"
#include "DNA_group_types.h"

struct Object;
struct Material;
struct Collection;

#define LANPR_MASTER_MODE_SOFTWARE     0
#define LANPR_MASTER_MODE_DPIX         1
#define LANPR_MASTER_MODE_SNAKE        2

#define LANPR_POST_PROCESSING_DISABLED 0
#define LANPR_POST_PROCESSING_ENABLED  1

#define LANPR_USE_DIFFERENT_TAPER      0
#define LANPR_USE_SAME_TAPER           1

#define LANPR_DISABLE_TIP_EXTEND       0
#define LANPR_ENABLE_TIP_EXTEND        1

#define LANPR_NORMAL_DONT_CARE         0
#define LANPR_NORMAL_DIRECTIONAL       1

#define LANPR_COMPONENT_MODE_ALL          0
#define LANPR_COMPONENT_MODE_OBJECT       1
#define LANPR_COMPONENT_MODE_MATERIAL     2
#define LANPR_COMPONENT_MODE_COLLECTION   3
//#define LANPR_COMPONENT_MODE_REST         4 // use _ALL

#define LANPR_COMPONENT_INCLUSIVE         0
#define LANPR_COMPONENT_EXCLUSIVE         1

#define LANPR_COMPONENT_LOGIG_OR          0
#define LANPR_COMPONENT_LOGIC_AND         1

struct DRWShadingGroup;

typedef struct LANPR_LineLayerComponent {
	struct LANPR_LineLayerComponent *next, *prev;

	struct Object     *object_select;
	struct Material   *material_select;
	struct Collection *collection_select;

	int component_mode;
	int what;

}LANPR_LineLayerComponent;

typedef struct LANPR_LineLayer {
	struct LANPR_LineLayer *next, *prev;

	int type;

	int use_qi_range;

	int qi_begin;
	int qi_end;        /* these are for QI Range thing... just occlusion levels */

	int enable_contour;
	int enable_crease;
	int enable_edge_mark;
	int enable_material_seperate;
	int enable_intersection;

	float thickness;    /* default is for contour */
	float thickness_crease;
	float thickness_material;
	float thickness_edge_mark;
	float thickness_intersection;

	float color[4];    /* default is for contour */
	float crease_color[4];
	float material_color[4];
	float edge_mark_color[4];
	float intersection_color[4];

	int normal_mode;
	int normal_effect_inverse;
    float normal_ramp_begin;
    float normal_ramp_end;
    float normal_thickness_begin;
    float normal_thickness_end;
	struct Object *normal_control_object;

	int logic_mode;      /* for component evaluation */

	int a;

	ListBase components;

	// still need legacy mode
	// should use runtime pointer
	struct DRWShadingGroup *shgrp;
	struct GPUBatch  *batch;

}LANPR_LineLayer;

#endif