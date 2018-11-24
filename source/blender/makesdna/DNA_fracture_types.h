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
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Martin Felke
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file DNA_fracture_types.h
 *  \ingroup DNA
 */
 
#ifndef DNA_FRACTURE_TYPES_H
#define DNA_FRACTURE_TYPES_H

#include "BLI_utildefines.h"
#include "DNA_mesh_types.h"
#include "DNA_object_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct DerivedMesh;
struct KDTree;
struct Shard;
struct RNG;

enum {
	MOD_RIGIDBODY_CENTROIDS = 0,
	MOD_RIGIDBODY_VERTICES = 1,
};

enum {
	SHARD_INTACT   = 1 << 0,
	SHARD_FRACTURED = 1 << 1,
	SHARD_SKIP = 1 << 2,
	SHARD_DELETE = 1 << 3,
};

/* Fracture Modifier */
enum {
	MOD_FRACTURE_BISECT_FAST      = (1 << 0),
	MOD_FRACTURE_BISECT_FAST_FILL = (1 << 1),
	MOD_FRACTURE_BOOLEAN          = (1 << 2),
	MOD_FRACTURE_BISECT_FILL      = (1 << 3),
	MOD_FRACTURE_BISECT           = (1 << 4),
	MOD_FRACTURE_BOOLEAN_FRACTAL  = (1 << 5),
};

enum {
	MOD_FRACTURE_OWN_VERTS       = (1 << 0),
	MOD_FRACTURE_OWN_PARTICLES   = (1 << 1),
	MOD_FRACTURE_EXTRA_VERTS     = (1 << 2),
	MOD_FRACTURE_EXTRA_PARTICLES = (1 << 3),
	//MOD_FRACTURE_GREASEPENCIL    = (1 << 4),
	MOD_FRACTURE_UNIFORM         = (1 << 5),
	MOD_FRACTURE_GRID            = (1 << 6),
	MOD_FRACTURE_CUSTOM          = (1 << 7),
};

enum {
	MOD_FRACTURE_SPLINTER_X      = (1 << 0),
	MOD_FRACTURE_SPLINTER_Y      = (1 << 1),
	MOD_FRACTURE_SPLINTER_Z      = (1 << 2),
};

enum {
	MOD_FRACTURE_CUTTER_X      = (1 << 0),
	MOD_FRACTURE_CUTTER_Y      = (1 << 1),
	MOD_FRACTURE_CUTTER_Z      = (1 << 2),
};

enum {
	MOD_FRACTURE_CENTROID      = (1 << 0),
	MOD_FRACTURE_VERTEX        = (1 << 1),
};

enum {
	MOD_FRACTURE_PREFRACTURED      = (1 << 0),
	MOD_FRACTURE_DYNAMIC           = (1 << 1),
	MOD_FRACTURE_EXTERNAL          = (1 << 2),
};

enum {
	MOD_FRACTURE_KEEP_BOTH          = (1 << 0),
	MOD_FRACTURE_KEEP_INTERSECT      = (1 << 1),
	MOD_FRACTURE_KEEP_DIFFERENCE     = (1 << 2),
};

enum {
	MOD_FRACTURE_NO_DYNAMIC_CONSTRAINTS     = (1 << 0),
	MOD_FRACTURE_MIXED_DYNAMIC_CONSTRAINTS  = (1 << 1),
	MOD_FRACTURE_ALL_DYNAMIC_CONSTRAINTS     = (1 << 2),
};

/* FractureModifierShared->flag */
enum {
	MOD_FRACTURE_REFRESH             = (1 << 0),
	MOD_FRACTURE_REFRESH_DYNAMIC     = (1 << 1),
	MOD_FRACTURE_REFRESH_CONSTRAINTS = (1 << 2),
	MOD_FRACTURE_REFRESH_AUTOHIDE    = (1 << 3),
};

typedef struct FractureQueueEntry {
	struct FractureQueueEntry *next, *prev;
	struct Shard *mi;
} FractureQueueEntry;

typedef struct AnimBind {
	int v;
	int v1;
	int v2;
	int mi;
	int poly;
	float offset[3];
	float no[3];
	float quat[4];
} AnimBind;

typedef struct FractureModifierData_Shared {

	/* shards, constraints, shared verts are the components of the FM RAM cache, actually */
	/* all shards of this Fracture Modifier, including invisible ones */
	/* during simulation, while iterating over this list, we filter out valid shards by their frame range. */
	ListBase shards;

	/* all constraints of this Fracture Modifier, make sure to only have constraints between visible shards
	 * else there will be simulation problems like lockups, crashes etc */
	ListBase constraints;

	/* used for storing shared vertices for automerge */
	ListBase automerge_shared_verts;

	/*volatile storage of shards being "hit" or fractured currently, needs to be cleaned up after usage! */
	ListBase dynamic_fracture_queue;

	/* stores the pre-autohide stage assembled mesh, useful for packgroup */
	struct Mesh *mesh_cached;

	/* runtime only data */
	/* store original vertices here (coords), to find them later and reuse their normals */
	struct KDTree *vertex_normals_tree;

	/* index pairs of faces of different shards, which are adjacent to each other */
	struct GHash *adjacent_face_pairs;

	/*used for autoconversion of former objects to clusters,marks object membership of each vert*/
	struct GHash *vertex_index_map;

	/* used for constraint building based on vertex proximity, temporary runtime data */
	struct GHash *vertex_island_map;

	/* used to collect materials from objects to be packed, temporary data */
	struct GHash *material_index_map;

	/*used to collect vertexgroups from objects to be packed, temporary data */
	struct GHash *defgrp_index_map;

	/* bound animation data, stores which shard is bound to which "guiding" vertices */
	struct AnimBind *anim_bind;

	/* for optimization, only refracture islands which really changed */
	struct KDTree *last_island_tree;

	/* islands from previous fracture process */
	struct Shard **last_islands;

	/* sigh, random number generator... */
	struct RNG *rng;

	float splinter_matrix[4][4];

	/* length of bound animation data */
	int anim_bind_len;

	/* number of of previous fracture islands */
	int last_islands_count;

	/*DANGER...
	 * what happens if the new compound object has more materials than fit into 1 short ? shouldnt happen but can...		  so reserve an int here better */

	/* start indexes in material and vertexgroup array (?) */
	int matstart;
	int defstart;

	/* refresh flags */
	int flag;

	/* markers of dynamic cache, in case the start / end is changed, the dynamic cache needs to be reallocated */
	int last_cache_start;
	int last_cache_end;

	char pad[4];

} FractureModifierData_Shared;

typedef struct Shard {
	struct Shard *next, *prev;

	/* geometry data */
	struct Mesh *mesh;

	/* rigidbody simulation data */
	struct RigidBodyOb *rigidbody;

	/* constraints this shard is participating in, for quicker access */
	struct RigidBodyShardCon **participating_constraints;

	/* array of neighbor indexes, especially useful for voronoi cell like fracture,
	 * note: might be incorrect after boolean cuts if the geometry changes,
	 * so cant rely for all purposes on it, like finding adjacent face pairs */
	int *neighbors;

	/* motion data cache, each shard manages its own motion data */
	/* length is endframe - startframe + 1 */
	float *locs;
	float *rots;
	float *vels;
	float *aves; //angular velocities

	char name[66]; /* MAX_ID_NAME */
	char pad1[2];

	/* valid start and end frames, between those the shard is considered valid and visible */
	int startframe;
	int endframe;
	/* length of participating constraints array */
	int participating_constraint_count;

	/* per FM an unique identifier number for each shard */
	int id;

	/* location and rotation of this shard after fracture.
	 * Actual Motion data is kept in rigidbody's loc, rot
	 * need this also for constraints probably */

	float loc[3];
	float rot[4];

	/* threshold and passive weights */
	float thresh_weight;
	float passive_weight;

	/* cluster membership index */
	int cluster_index;

	/* index of which constraint this shard belongs to */
	int constraint_index;

	/* to which object in the rigidbody collection this shard belongs */
	int object_index;

	/*store number of used materials here, from the original object*/
	int totcol;

	/*store number of used vertexgroups here, from the original object*/
	int totdef;

	/* whether this shard is already has been processed for this actual fracture event */
	int fractured;

	/* corners of the shard's bounding box */
	float min[3], max[3];

	/* number of adjacent shards */
	int neighbor_count;

	/*store raw, unprocessed centroid here
	(might change when mesh shape changes via boolean / bisect) */
	float raw_centroid[3];

	/* same with the raw cell volume */
	float raw_volume;

	/* last impact location on this shard */
	float impact_loc[3];

	/* size of impact area (bbox-like) */
	float impact_size[3];

	char pad[4];
} Shard;

//TODO
typedef struct SharedVertGroup {
	struct SharedVertGroup *next, *prev;
	float rest_co[3]; /* original coordinate after fracture */
	float delta[3];
	int index, excession_frame;
	int exceeded, deltas_set, moved;
	char pad[4];
	ListBase verts;
} SharedVertGroup;

typedef struct SharedVert {
	struct SharedVert *next, *prev;
	float rest_co[3];
	float delta[3];
	int index, excession_frame;
	int exceeded, deltas_set, moved;
	char pad[4];
} SharedVert;

#ifdef __cplusplus
}
#endif

#endif /* DNA_FRACTURE_TYPES_H */
