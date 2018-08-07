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

/** \file blender/blenkernel/BKE_fracture.h
 *  \ingroup blenkernel
 */

#ifndef BKE_FRACTURE_H
#define BKE_FRACTURE_H

#include "BLI_sys_types.h"
#include "BKE_scene.h"

struct FracMesh;
struct Shard;
struct Mesh;

struct RigidBodyWorld;
struct FractureModifierData;
struct Object;
struct Group;
struct MeshIsland;
struct RigidBodyShardCon;
struct GHash;

struct BoundBox;
struct MVert;
struct MPoly;
struct MLoop;
struct MEdge;

struct BMesh;
struct CustomData;
struct Scene;
struct Main;
struct KDTree;

typedef int ShardID;

typedef struct FracPoint {
	float co[3];
	float offset[3];
} FracPoint;

typedef struct FracPointCloud {
	struct FracPoint *points;   /* just a bunch of positions in space*/
	int totpoints; /* number of positions */
} FracPointCloud;


void BKE_fracture_dynamic_free(struct FractureModifierData *fmd,
                               bool do_free_sequence, bool do_free_rigidbody, struct Scene* scene);

struct Mesh* BKE_fracture_prefractured_apply(struct FractureModifierData *fmd, struct Object *ob, struct Mesh *mesh,
                                             struct Depsgraph *depsgraph);

struct Mesh* BKE_fracture_dynamic_apply(struct FractureModifierData *fmd, struct Object *ob, struct Mesh *mesh, struct Scene *scene);
struct Mesh* BKE_fracture_external_apply(struct FractureModifierData *fmd, struct Object* ob, struct Mesh* mesh,
                                         struct Mesh *derivedData, struct Scene *scene);

struct Mesh* BKE_fracture_mesh_from_packdata(struct FractureModifierData *fmd, struct Mesh *derivedData);

void BKE_fracture_initialize(struct FractureModifierData* fmd, struct Object *ob,
                                     struct Mesh* mesh, struct Depsgraph *depsgraph);

void BKE_fracture_autohide_refresh(struct FractureModifierData* fmd, struct Object *ob);
void BKE_fracture_automerge_refresh(struct FractureModifierData* fmd);
struct Mesh *BKE_fracture_result_mesh(struct FractureModifierData* fmd, struct Mesh *dm, struct Object* ob, bool validMesh,
                                      struct Scene* scene);

FracPointCloud BKE_fracture_points_get(struct Depsgraph *depsgraph, struct FractureModifierData *emd,
                                       struct Object *ob, struct Mesh *fracmesh, ShardID id);

void BKE_fracture_face_calc_center_mean(struct Mesh *dm, struct MPoly *mp, float r_cent[3]);

struct Shard* BKE_fracture_shard_find(struct ListBase *shards, ShardID id);

void BKE_fracture_face_pairs(struct FractureModifierData *fmd, struct Mesh *dm, struct Object *ob);
void BKE_fracture_shared_vert_groups(struct FractureModifierData* fmd, struct Mesh* dm, struct ListBase *shared_verts);
void BKE_fracture_shared_verts_free(struct ListBase* lb);

struct Mesh *BKE_fracture_autohide_do(struct FractureModifierData *fmd, struct Mesh *dm, struct Object *ob, struct Scene* sc);

struct FracMesh* BKE_fracture_fracmesh_copy(struct FracMesh* fm);
void BKE_fracture_simulation_free(struct FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, struct Scene *scene);
void BKE_fracture_meshislands_free(struct FractureModifierData* fmd, struct ListBase* meshIslands, bool do_free_rigidbody,
                                   struct Scene* scene);

void BKE_fracture_free(struct FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, struct Scene *scene);

void BKE_fracture_do(struct FractureModifierData *fmd, ShardID id, struct Object *obj, struct Mesh *dm,
                     struct Depsgraph *depsgraph, struct Main *bmain);

void BKE_fracture_normal_find(struct Mesh *dm, struct KDTree *tree, float co[3], short no[3], short rno[3], float range);
void BKE_fracture_physics_mesh_normals_fix(struct FractureModifierData *fmd, struct Shard* s, struct MeshIsland* mi, int i,
                                           struct Mesh* orig_dm);

void BKE_fracture_collect_layers(struct Shard* s, struct Mesh *dm, int vertstart, int polystart, int loopstart, int edgestart);


void BKE_fracture_animated_loc_rot(struct FractureModifierData *fmd, struct Object *ob, bool do_bind, struct Depsgraph *depsgraph);


void BKE_fracture_constraints_refresh(struct FractureModifierData *fmd, struct Object *ob, struct Scene* scene);


/* external mode */
struct MeshIsland* BKE_fracture_mesh_island_add(struct Main* bmain, struct FractureModifierData *fmd, struct Object* own,
                                                struct Object *target, struct Scene *scene);

void BKE_fracture_mesh_island_remove(struct FractureModifierData *fmd, struct MeshIsland *mi, struct Scene* scene);
void BKE_fracture_mesh_island_remove_all(struct FractureModifierData *fmd, struct Scene* scene);

void BKE_fracture_mesh_constraint_remove(struct FractureModifierData *fmd, struct RigidBodyShardCon* con, struct Scene *scene);
void BKE_fracture_constraints_free(struct FractureModifierData *fmd, struct Scene *scene);

int BKE_fracture_update_visual_mesh(struct FractureModifierData *fmd, struct Object *ob, bool do_custom_data);

struct RigidBodyShardCon *BKE_fracture_mesh_constraint_create(struct Scene *scene, struct FractureModifierData *fmd,
                                                     struct MeshIsland *mi1, struct MeshIsland *mi2, short con_type);

void BKE_fracture_mesh_constraint_remove_all(struct FractureModifierData *fmd, struct Scene *scene);


/* Shard handling */
struct Shard *BKE_fracture_shard_create(struct MVert *mvert, struct MPoly *mpoly, struct MLoop *mloop,
                                        struct MEdge *medge, int totvert, int totpoly, int totloop, int totedge, bool copy);

struct Shard* BKE_fracture_shard_copy(struct Shard *s);
struct Mesh* BKE_fracture_shard_to_mesh(struct Shard *s, bool doCustomData);

float BKE_shard_calc_minmax(struct Shard *shard);
void BKE_fracture_shard_free(struct Shard *s, bool doCustomData);

bool BKE_fracture_shard_center_centroid_area(struct Shard *shard, float cent[3]);
void BKE_fracture_custom_data_mesh_to_shard(struct Shard *s, struct Mesh *dm);


void BKE_fracture_fracmesh_free(struct FracMesh *fm, bool doCustomData);
struct FracMesh *BKE_fracture_fracmesh_create(void);

void BKE_bm_mesh_hflag_flush_vert(struct BMesh *bm, const char hflag);

void BKE_fracture_constraint_create(struct Scene *scene, struct FractureModifierData *fmd,
                                                     struct MeshIsland *mi1, struct MeshIsland *mi2, short con_type, float thresh);

bool BKE_fracture_dynamic_lookup_mesh_state(struct FractureModifierData *fmd, int frame, int do_lookup, struct Scene* scene);

struct MDeformVert* BKE_fracture_shards_to_islands(struct FractureModifierData* fmd, struct Object* ob, struct Mesh *orig_dm,
                                                   struct Scene *scene);

void BKE_fracture_fill_vgroup(struct FractureModifierData *rmd, struct Mesh *dm, struct MDeformVert *dvert, struct Object *ob,
                              struct Mesh *old_cached);

void BKE_fracture_do_halving(struct FractureModifierData *fmd, struct Object* ob, struct Mesh *dm, struct Mesh *orig_dm,
                             bool is_prehalving, ShardID id, struct Scene* scene);

struct Mesh* BKE_fracture_assemble_mesh_from_shards(struct FractureModifierData *fmd, bool doCustomData, bool use_packed);

void BKE_fracture_modifier_free(struct FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, struct Scene *scene);

void BKE_fracture_mesh_island_free(struct FractureModifierData *rmd, struct MeshIsland *mi, bool remove_rigidbody, struct Scene* scene);

short BKE_fracture_collect_materials(struct Main* bmain, struct Object* o, struct Object* ob, int matstart, struct GHash** mat_index_map);

struct Mesh *BKE_fracture_prefractured_do(struct FractureModifierData *fmd, struct Object *ob, struct Mesh *dm,
                                          struct Mesh *orig_dm, char names [][66], int count, struct Scene* scene);

#endif /* BKE_FRACTURE_H */
