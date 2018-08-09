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
	MOD_FRACTURE_GREASEPENCIL    = (1 << 4),
	MOD_FRACTURE_UNIFORM         = (1 << 5),
	MOD_FRACTURE_GRID            = (1 << 6),
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

typedef struct ShardSequence {
	struct ShardSequence *next, *prev;
	struct FracMesh *frac_mesh;
	int frame;
	int is_new;
} ShardSequence;

typedef struct MeshIslandSequence {
	struct MeshIslandSequence *next, *prev;
	struct Mesh *visible_dm;
	ListBase meshIslands;
	int frame;
	int is_new;
} MeshIslandSequence;

typedef struct FractureID {
	struct FractureID *next, *prev;
	int shardID;
	char pad[4];
} FractureID;

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

	struct FracMesh *frac_mesh; /* store only the current fracmesh here first, later maybe an entire history...*/
	struct Mesh *dm; /*untransformed assembled shard mesh */

	struct BMesh *visible_mesh; /* for halving / existing islands */
	struct Mesh *visible_mesh_cached; /*transformed output mesh*/

	ListBase meshIslands, meshConstraints;
	ListBase islandShards; /* the additional islands after splitting the shards, so count(meshislands) = count(shards + islandshards) */

	//dynamic stuff... heavy TODO
	ListBase shard_sequence; /* used as mesh cache / history for dynamic fracturing, for shards (necessary for conversion to DM) */
	ListBase meshIsland_sequence; /* used as mesh cache / history for dynamic fracturing, for meshIslands (necessary for loc/rot "pointcache") */
	ShardSequence *current_shard_entry; /*volatile storage of current shard entry, so we dont have to search in the list */
	MeshIslandSequence *current_mi_entry; /*analogous to current shard entry */
	ListBase fracture_ids; /*volatile storage of shards being "hit" or fractured currently, needs to be cleaned up after usage! */
	ListBase shared_verts; /* used for storing shared vertices for automerge */
	ListBase pack_storage; /*used to store packed geometry when switching modes */

	//runtime only data, can be shared too
	struct KDTree *nor_tree; /* store original vertices here (coords), to find them later and reuse their normals */
	struct GHash *face_pairs;
	struct GHash *vert_index_map; /*used for autoconversion of former objects to clusters, marks object membership of each vert*/
	struct GHash *vertex_island_map; /* used for constraint building based on vertex proximity, temporary data */
	struct GHash *material_index_map; /* used to collect materials from objects to be packed, temporary data */
	struct GHash *defgrp_index_map; /*used to collect vertexgroups from objects to be packed, temporary data */
	struct AnimBind *anim_bind; /* bound animation data */

	int anim_bind_len;

	/*DANGER... what happens if the new compound object has more materials than fit into 1 short ? shouldnt happen but can...*/
	/*so reserve an int here better */
	int matstart;
	int defstart;

	int refresh;
	int refresh_constraints;
	int refresh_autohide;
	int reset_shards;

	char pad[4];

} FractureModifierData_Shared;

typedef struct Shard {
	struct Shard *next, *prev;
	struct MVert *mvert;
	struct MPoly *mpoly;
	struct MLoop *mloop;
	struct MEdge *medge;

	struct CustomData vertData;
	struct CustomData polyData;
	struct CustomData loopData;
	struct CustomData edgeData;

	int totvert, totpoly, totloop, totedge;
	int pad[2];
	
	int *cluster_colors;
	float min[3], max[3];
	float centroid[3];  /* centroid of shard, calculated during fracture */
	float raw_centroid[3];  /*store raw, unprocessed centroid here (might change when mesh shape changes via boolean / bisect) */
	int *neighbor_ids;  /* neighbors of me... might be necessary for easier compounding or fracture, dont need to iterate over all */
	int shard_id;       /* the identifier */
	int neighbor_count; /* counts of neighbor islands */
	int parent_id;      /* the shard from which this shard originates, we keep all shards in the shardmap */
	int flag;           /* flag for fracture state (INTACT, FRACTURED)*/
	int setting_id;     /* to which setting this shard belongs, -1 for none or no settings available*/
	float raw_volume;
	float impact_loc[3]; /* last impact location on this shard */
	float impact_size[3]; /* size of impact area (simplified) */
	//char pad2[4];
} Shard;

typedef struct MeshIsland {
	struct MeshIsland *next, *prev;
	struct BMVert **vertices;
	struct MVert **vertices_cached;
	float *vertco;
	short *vertno;
	struct Mesh *physics_mesh;
	struct RigidBodyOb *rigidbody;
	int *neighbor_ids;
	int *vertex_indices;
	struct BoundBox *bb;
	struct RigidBodyShardCon **participating_constraints;
	float *locs;
	float *rots;

	char name[66]; /* MAX_ID_NAME */
	char pad1[2];

	int start_frame;
	int frame_count;
	int participating_constraint_count;
	int vertex_count, id, neighbor_count;
	float centroid[3], start_co[3];
	float rot[4]; /*hrm, need this for constraints probably */
	float thresh_weight, ground_weight;
	int linear_index;  /* index in rigidbody world */
	int particle_index;
	int constraint_index;
	int object_index;
	int totcol; /*store number of used materials here, from the original object*/
	int totdef; /*store number of used vertexgroups here, from the original object*/
	char pad[4];
} MeshIsland;

typedef struct FracMesh {
	struct KDTree *last_shard_tree;
	struct Shard **last_shards;
	ListBase shard_map;     /* groups mesh elements to islands, generated by fracture itself */
	int shard_count;        /* how many islands we have */
	short cancel;           /* whether the process is cancelled (from the job, ugly, but this way we dont need the entire modifier) */
	short running;          /* whether the process is currently in progress, so the modifier wont be touched from the main thread */
	int progress_counter;   /* counts progress */
	int last_expected_shards;
} FracMesh;

typedef struct SharedVertGroup {
	struct SharedVertGroup *next, *prev;
	float rest_co[3];
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
