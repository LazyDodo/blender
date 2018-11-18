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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor(s): Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file DNA_hair_types.h
 *  \ingroup DNA
 */

#ifndef __DNA_HAIR_TYPES_H__
#define __DNA_HAIR_TYPES_H__

#include "DNA_defs.h"
#include "DNA_listBase.h"
#include "DNA_ID.h"
#include "DNA_meshdata_types.h"
#include "DNA_customdata_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum HairFollicleFlag
{
	HAIR_FOLLICLE_SELECT            = (1 << 0),
} HairFollicleFlag;

/* Root point (follicle) of a hair on a surface */
typedef struct HairFollicle {
	int flag;
	unsigned int curve;         /* Index of the curve used by the fiber */
	MeshSample mesh_sample;     /* Sample on the scalp mesh for the root vertex */
} HairFollicle;

/* Collection of hair roots on a surface */
typedef struct HairPattern {
	struct HairFollicle *follicles;
	int num_follicles;
	int pad;

	/* Per follicle custom data */
	struct CustomData fdata;
} HairPattern;

typedef struct HairFiberCurve {
	int vertstart;                      /* Offset in the vertex array where the curve starts */
	int numverts;                       /* Number of vertices in the curve */

	/* Shape */
	float taper_length;                 /* Distance at which final thickness is reached */
	float taper_thickness;              /* Relative thickness of the strand */
} HairFiberCurve;

typedef struct HairFiberVertex {
	int flag;
	float co[3];
} HairFiberVertex;

typedef enum HairFiberVertexFlag
{
	HAIR_VERTEX_SELECT          = (1 << 0),
} HairFiberVertexFlag;

/* Hair curve data */
typedef struct HairCurveData
{
	/* Curves for shaping hair fibers */
	struct HairFiberCurve *curves;
	/* Control vertices on curves */
	struct HairFiberVertex *verts;
	/* Number of curves */
	int totcurves;
	/* Number of curve vertices */
	int totverts;

	/* Per vertex custom data */
	struct CustomData vdata;
	/* Per curve custom data */
	struct CustomData cdata;
} HairCurveData;

/* Editable hair data */
typedef struct EditHair {
	/* Set of hair follicles on the scalp mesh */
	struct HairPattern *pattern;

	/* Curve data */
	HairCurveData curve_data;
} EditHair;

/* not saved in file! */
typedef struct HairSystem_Runtime {
	int64_t cd_dirty_follicles;

	void *draw_batch_cache;     /* Data buffers for drawing */
} HairSystem_Runtime;

typedef struct HairSystem {
	ID id;                      /* Hair data is a datablock */
	struct AnimData *adt;       /* Animation data - for animating settings */

	struct BoundBox *bb;

	struct EditHair *edithair;

	int flag;
	int pad;

	/* Set of hair follicles on the scalp mesh */
	struct HairPattern *pattern;

	/* Curve data */
	HairCurveData curve_data;

	struct Material **mat;      /* Material slots */
	short totcol;               /* Number of material slots */
	short pad2;
	int pad3;

	struct HairDrawSettings *draw_settings;

	HairSystem_Runtime runtime;
} HairSystem;

// typedef enum eHairSystemFlag
// {
// } eHairSystemFlag;

typedef struct HairDrawSettings
{
	short follicle_mode;
	short fiber_mode;
	short shape_flag;
	short pad;
	
	float shape;
	float root_radius;
	float tip_radius;
	float radius_scale;
} HairDrawSettings;

typedef enum eHairDrawFollicleMode
{
	HAIR_DRAW_FOLLICLE_NONE     = 0,
	HAIR_DRAW_FOLLICLE_POINTS   = 1,
} eHairDrawFollicleMode;

typedef enum eHairDrawFiberMode
{
	HAIR_DRAW_FIBER_NONE        = 0,
	HAIR_DRAW_FIBER_CURVES      = 1,
} eHairDrawFiberMode;

typedef enum eHairDrawShapeFlag {
	HAIR_DRAW_CLOSE_TIP         = (1<<0),
} eHairDrawShapeFlag;

#ifdef __cplusplus
}
#endif

#endif /* __DNA_HAIR_TYPES_H__ */
