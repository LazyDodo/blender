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
#include "BKE_customdata.h"

struct Mesh;

struct RigidBodyWorld;
struct FractureModifierData;
struct ModifierData;
struct Object;
struct Group;
struct Shard;
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


static const CustomDataMask CD_MASK_ISLAND =
		CD_MASK_MDEFORMVERT |
		CD_MASK_PROP_FLT | CD_MASK_PROP_INT | CD_MASK_PROP_STR | CD_MASK_MDISPS |
		CD_MASK_MLOOPUV | CD_MASK_MLOOPCOL |
		CD_MASK_RECAST | CD_MASK_PAINT_MASK |
		CD_MASK_GRID_PAINT_MASK | CD_MASK_MVERT_SKIN | CD_MASK_FREESTYLE_EDGE | CD_MASK_FREESTYLE_FACE |
		CD_MASK_CUSTOMLOOPNORMAL | CD_MASK_FACEMAP;

void BKE_fracture_points(struct FractureModifierData *fmd, struct Object* obj, struct Shard *mi,
                         struct Depsgraph* depsgraph, struct Main *bmain, struct Scene *scene, bool is_initial);

void BKE_fracture_shard_by_points(struct FractureModifierData *fmd, struct FracPointCloud *pointcloud,
                                  struct Object *ob, struct Shard* mi, struct Scene *scene, struct Main *bmain);
void BKE_fracture_autohide_refresh(struct FractureModifierData* fmd, struct Object *ob, struct Mesh *me_assembled);
void BKE_fracture_automerge_refresh(struct FractureModifierData* fmd, struct Mesh *me_assembled);

FracPointCloud BKE_fracture_points_get(struct Depsgraph *depsgraph, struct FractureModifierData *emd,
                                       struct Object *ob, struct Shard *mi);

void BKE_fracture_face_calc_center_mean(struct Mesh *dm, struct MPoly *mp, float r_cent[3]);

void BKE_fracture_face_pairs(struct FractureModifierData *fmd, struct Mesh *dm, struct Object *ob);
void BKE_fracture_shared_vert_groups(struct FractureModifierData* fmd, struct Mesh* dm, struct ListBase *shared_verts);
void BKE_fracture_shared_verts_free(struct ListBase* lb);

struct Mesh *BKE_fracture_autohide_do(struct FractureModifierData *fmd, struct Mesh *dm, struct Object *ob, struct Scene* sc);

void BKE_fracture_meshislands_free(struct FractureModifierData* fmd, struct Scene* scene);

void BKE_fracture_free(struct FractureModifierData *fmd, struct Scene *scene);

void BKE_fracture_do(struct FractureModifierData *fmd, struct Shard *mi, struct Object *obj,
                     struct Depsgraph *depsgraph, struct Main *bmain, struct Scene *scene, bool is_init);


void BKE_fracture_animated_loc_rot(struct FractureModifierData *fmd, struct Object *ob, bool do_bind, struct Depsgraph *depsgraph);

void BKE_fracture_constraints_refresh(struct FractureModifierData *fmd, struct Object *ob, struct Scene* scene);


/* external mode */
struct Shard* BKE_fracture_mesh_island_add(struct Main* bmain, struct FractureModifierData *fmd, struct Object* own,
                                                struct Object *target, struct Scene *scene, int id);

void BKE_fracture_mesh_island_remove(struct FractureModifierData *fmd, struct Shard *mi, struct Scene* scene);
void BKE_fracture_mesh_island_remove_all(struct FractureModifierData *fmd, struct Scene* scene);

void BKE_fracture_mesh_constraint_remove(struct FractureModifierData *fmd, struct RigidBodyShardCon* con, struct Scene *scene);
void BKE_fracture_constraints_free(struct FractureModifierData *fmd, struct RigidBodyWorld *rbw);

struct RigidBodyShardCon *BKE_fracture_mesh_constraint_create(struct Scene *scene, struct FractureModifierData *fmd,
                                                     struct Shard *mi1, struct Shard *mi2, short con_type);

void BKE_fracture_mesh_constraint_remove_all(struct FractureModifierData *fmd, struct Scene *scene);


bool BKE_fracture_mesh_center_centroid_area(struct Mesh *shard, float cent[3]);

void BKE_bm_mesh_hflag_flush_vert(struct BMesh *bm, const char hflag);

void BKE_fracture_constraint_create(struct Scene *scene, struct FractureModifierData *fmd,
                                                     struct Shard *mi1, struct Shard *mi2, short con_type, float thresh);

void BKE_fracture_split_islands(struct FractureModifierData *fmd, struct Object* ob, struct Mesh *me, struct Mesh ***temp_islands,
int *count);

struct Mesh* BKE_fracture_assemble_mesh_from_islands(struct FractureModifierData *fmd, struct Object *ob, float ctime);

void BKE_fracture_modifier_free(struct FractureModifierData *fmd, struct Scene *scene);

void BKE_fracture_mesh_island_free(struct FractureModifierData *fmd, struct Shard *mi, struct Scene* scene);

short BKE_fracture_collect_materials(struct Main* bmain, struct Object* o, struct Object* ob, int matstart, struct GHash** mat_index_map);

struct Mesh* BKE_fracture_mesh_copy(struct Mesh* source, struct Object* ob);
struct BMesh* BKE_fracture_mesh_to_bmesh(struct Mesh* me);
struct Mesh* BKE_fracture_bmesh_to_mesh(struct BMesh* bm);

void BKE_update_velocity_layer(struct FractureModifierData *fmd, struct Mesh *dm);
bool BKE_rigidbody_remove_modifier(struct RigidBodyWorld* rbw, struct ModifierData *md, struct Object *ob);
void BKE_fracture_external_constraints_setup(struct FractureModifierData *fmd, struct Scene *scene, struct Object *ob, int frame);

struct Mesh* BKE_fracture_apply(struct FractureModifierData *fmd, struct Object *ob, struct Mesh *me, struct Depsgraph* depsgraph);

struct Shard *BKE_fracture_mesh_island_create(struct Mesh* me, struct Scene *scene, struct Object *ob, int frame);
void BKE_fracture_mesh_boundbox_calc(struct Mesh *me, float r_loc[], float r_size[]);
void BKE_fracture_mesh_free(struct Mesh *me);

void BKE_fracture_meshisland_check_realloc_cache(struct FractureModifierData *fmd, struct RigidBodyWorld *rbw, struct Shard* mi, int frame);

bool BKE_fracture_meshisland_check_frame(struct FractureModifierData *fmd, struct Shard* mi, int frame);
void BKE_fracture_dynamic_do(struct FractureModifierData *fmd, struct Object* ob, struct Scene* scene,
                             struct Depsgraph* depsgraph, struct Main* bmain);

void BKE_fracture_clear_cache(struct FractureModifierData* fmd, struct Scene *scene);
void BKE_fracture_meshisland_vertexgroups_do(struct FractureModifierData *fmd, struct Object *ob, struct Shard* mi);
void BKE_fracture_meshislands_pack(struct FractureModifierData *fmd, struct Object* obj, struct Main* bmain, struct Scene* scene);

void BKE_fracture_postprocess_meshisland(struct FractureModifierData *fmd, struct Object* ob, struct Shard* mi,
                                         struct Mesh ***temp_meshs, int count, struct Scene* scene, int frame,
                                         struct Shard **shards);

void BKE_fracture_meshisland_normals_fix(struct FractureModifierData *fmd, struct Shard* mi, struct Mesh* orig_me);

void BKE_fracture_copy_customdata(struct CustomData* src, struct CustomData* dst, CustomDataMask mask, int src_ofs,
                                  int dst_ofs, int copyelem, int totelem);

bool BKE_fracture_check_valid_shard(struct FractureModifierData *fmd, struct Shard *mi, struct Scene *scene);
void BKE_fracture_duplis_to_shards(struct FractureModifierData *fmd, struct Object *ob, struct Scene *scene,
                                   struct Depsgraph *depsgraph, struct Main *bmain, int frame);

bool BKE_fracture_handle_initial_shards(struct FractureModifierData* fmd, struct Object* ob,
                                        struct Depsgraph *depsgraph, struct Main* bmain,
                                        struct Scene* scene, int frame);

void BKE_fracture_meshislands_connect(struct Scene* scene, struct FractureModifierData* fmd,
                                      struct Shard* mi, struct Shard* mi2, short con_type, float thresh);

#endif /* BKE_FRACTURE_H */
