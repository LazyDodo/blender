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
 * The Original Code is Copyright (C) 2005 by the Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Daniel Dunbar
 *                 Ton Roosendaal,
 *                 Ben Batt,
 *                 Brecht Van Lommel,
 *                 Campbell Barton
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_build.c
 *  \ingroup modifiers
 */


#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"
#include "BLI_rand.h"
#include "BLI_math_vector.h"
#include "BLI_ghash.h"

#include "DNA_meshdata_types.h"
#include "DNA_mesh_types.h"
#include "DNA_object_types.h"

#include "DEG_depsgraph_query.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_particle.h"
#include "BKE_scene.h"



#ifdef _OPENMP
#  include "BKE_mesh.h"  /* BKE_MESH_OMP_LIMIT */
#endif

static void initData(ModifierData *md)
{
	BuildModifierData *bmd = (BuildModifierData *) md;

	bmd->start = 1.0;
	bmd->length = 100.0;
}

static void copyData(ModifierData *md, ModifierData *target)
{
#if 0
	BuildModifierData *bmd = (BuildModifierData *) md;
	BuildModifierData *tbmd = (BuildModifierData *) target;
#endif
	modifier_copyData_generic(md, target);
}

static bool dependsOnTime(ModifierData *UNUSED(md))
{
	return true;
}

static Mesh *applyModifier(ModifierData *md, struct Depsgraph *depsgraph,
                           Object *UNUSED(ob), struct Mesh *mesh,
                           ModifierApplyFlag UNUSED(flag))
{
	Mesh *result; /* TODO(sybren): remove and replace with direct modification of mesh? */
	BuildModifierData *bmd = (BuildModifierData *) md;
	unsigned int i, j, k;
	unsigned int numFaces_dst, numEdges_dst, numLoops_dst = 0;
	unsigned int *vertMap, *edgeMap, *faceMap;
	float frac;
	MPoly *mpoly_dst;
	MLoop *ml_dst, *ml_src /*, *mloop_dst */;
	GHashIterator gh_iter;
	/* maps vert indices in old mesh to indices in new mesh */
	GHash *vertHash = BLI_ghash_int_new("build ve apply gh");
	/* maps edge indices in old mesh to indices in new mesh */
	GHash *edgeHash = BLI_ghash_int_new("build ed apply gh");

	const unsigned int numVert_src = mesh->totvert;
	const unsigned int numEdge_src = mesh->totedge;
	const unsigned int numPoly_src = mesh->totpoly;
	MPoly *mpoly_src = mesh->mpoly;
	MLoop *mloop_src = mesh->mloop;
	MEdge *medge_src = mesh->medge;
	MVert *mvert_src = mesh->mvert;

	vertMap = MEM_malloc_arrayN(numVert_src, sizeof(*vertMap), "build modifier vertMap");
	edgeMap = MEM_malloc_arrayN(numEdge_src, sizeof(*edgeMap), "build modifier edgeMap");
	faceMap = MEM_malloc_arrayN(numPoly_src, sizeof(*faceMap), "build modifier faceMap");

	range_vn_u(vertMap, numVert_src, 0);
	range_vn_u(edgeMap, numEdge_src, 0);
	range_vn_u(faceMap, numPoly_src, 0);

	struct Scene *scene = DEG_get_input_scene(depsgraph);
	frac = (BKE_scene_frame_get(scene) - bmd->start) / bmd->length;
	CLAMP(frac, 0.0f, 1.0f);
	if (bmd->flag & MOD_BUILD_FLAG_REVERSE) {
		frac = 1.0f - frac;
	}

	numFaces_dst = numPoly_src * frac;
	numEdges_dst = numEdge_src * frac;

	BKE_mesh_tessface_ensure(mesh);

	/* if there's at least one face, build based on faces */
	if (numFaces_dst) {
		MPoly *mpoly, *mp;
		MLoop *ml, *mloop;
		unsigned int new_idx;
		
		if (bmd->flag & MOD_BUILD_FLAG_RANDOMIZE) {
			BLI_array_randomize(faceMap, sizeof(*faceMap),
			                    numPoly_src, bmd->seed);
		}

		/* get the set of all vert indices that will be in the final mesh,
		 * mapped to the new indices
		 */
		mpoly = mpoly_src;
		mloop = mloop_src;
		new_idx = 0;
		for (i = 0; i < numFaces_dst; i++) {
			mp = mpoly + faceMap[i];
			ml = mloop + mp->loopstart;

			for (j = 0; j < mp->totloop; j++, ml++) {
				void **val_p;
				if (!BLI_ghash_ensure_p(vertHash, SET_UINT_IN_POINTER(ml->v), &val_p)) {
					*val_p = SET_UINT_IN_POINTER(new_idx);
					new_idx++;
				}
			}

			numLoops_dst += mp->totloop;
		}
		BLI_assert(new_idx == BLI_ghash_len(vertHash));

		/* get the set of edges that will be in the new mesh (i.e. all edges
		 * that have both verts in the new mesh)
		 */
		new_idx = 0;
		for (i = 0; i < numEdge_src; i++) {
			MEdge *me = medge_src + i;

			if (BLI_ghash_haskey(vertHash, SET_UINT_IN_POINTER(me->v1)) &&
			    BLI_ghash_haskey(vertHash, SET_UINT_IN_POINTER(me->v2)))
			{
				BLI_ghash_insert(edgeHash, SET_UINT_IN_POINTER(i), SET_UINT_IN_POINTER(new_idx));
				new_idx++;
			}
		}
		BLI_assert(new_idx == BLI_ghash_len(edgeHash));
	}
	else if (numEdges_dst) {
		MEdge *medge, *me;
		unsigned int new_idx;

		if (bmd->flag & MOD_BUILD_FLAG_RANDOMIZE)
			BLI_array_randomize(edgeMap, sizeof(*edgeMap),
			                    numEdge_src, bmd->seed);

		/* get the set of all vert indices that will be in the final mesh,
		 * mapped to the new indices
		 */
		medge = medge_src;
		new_idx = 0;
		BLI_assert(BLI_ghash_len(vertHash) == 0);
		for (i = 0; i < numEdges_dst; i++) {
			void **val_p;
			me = medge + edgeMap[i];

			if (!BLI_ghash_ensure_p(vertHash, SET_UINT_IN_POINTER(me->v1), &val_p)) {
				*val_p = SET_UINT_IN_POINTER(new_idx);
				new_idx++;
			}
			if (!BLI_ghash_ensure_p(vertHash, SET_UINT_IN_POINTER(me->v2), &val_p)) {
				*val_p = SET_UINT_IN_POINTER(new_idx);
				new_idx++;
			}
		}
		BLI_assert(new_idx == BLI_ghash_len(vertHash));

		/* get the set of edges that will be in the new mesh */
		for (i = 0; i < numEdges_dst; i++) {
			j = BLI_ghash_len(edgeHash);

			BLI_ghash_insert(edgeHash,  SET_INT_IN_POINTER(edgeMap[i]), SET_INT_IN_POINTER(j));
		}
	}
	else {
		int numVerts = numVert_src * frac;

		if (bmd->flag & MOD_BUILD_FLAG_RANDOMIZE) {
			BLI_array_randomize(vertMap, sizeof(*vertMap),
			                    numVert_src, bmd->seed);
		}

		/* get the set of all vert indices that will be in the final mesh,
		 * mapped to the new indices
		 */
		for (i = 0; i < numVerts; i++) {
			BLI_ghash_insert(vertHash, SET_INT_IN_POINTER(vertMap[i]), SET_INT_IN_POINTER(i));
		}
	}

	/* now we know the number of verts, edges and faces, we can create
	 * the mesh
	 */
	printf("Creating new mesh from %s=%p\n", mesh->id.name, mesh);
	result = BKE_mesh_from_template(mesh, BLI_ghash_len(vertHash), BLI_ghash_len(edgeHash),
	                                0, numLoops_dst, numFaces_dst);
	printf("    : %s=%p\n", result->id.name, result);

	/* copy the vertices across */
	GHASH_ITER (gh_iter, vertHash) {
		MVert source;
		MVert *dest;
		unsigned int oldIndex = GET_UINT_FROM_POINTER(BLI_ghashIterator_getKey(&gh_iter));
		unsigned int newIndex = GET_UINT_FROM_POINTER(BLI_ghashIterator_getValue(&gh_iter));

		source = mvert_src[oldIndex];
		dest = &result->mvert[newIndex];

		CustomData_copy_data(&mesh->vdata, &result->vdata, oldIndex, newIndex, 1);
		*dest = source;
	}

	/* copy the edges across, remapping indices */
	GHASH_ITER (gh_iter, edgeHash) {
		MEdge source;
		MEdge *dest;
		unsigned int oldIndex = GET_UINT_FROM_POINTER(BLI_ghashIterator_getKey(&gh_iter));
		unsigned int newIndex = GET_UINT_FROM_POINTER(BLI_ghashIterator_getValue(&gh_iter));
		BLI_assert(newIndex < BLI_ghash_len(edgeHash));

		source = medge_src[oldIndex];
		dest = &result->medge[newIndex];
		
		CustomData_copy_data(&mesh->edata, &result->edata, oldIndex, newIndex, 1);
		*dest = source;
		dest->v1 = GET_UINT_FROM_POINTER(BLI_ghash_lookup(vertHash, SET_UINT_IN_POINTER(source.v1)));
		dest->v2 = GET_UINT_FROM_POINTER(BLI_ghash_lookup(vertHash, SET_UINT_IN_POINTER(source.v2)));
	}

	mpoly_dst = result->mpoly;
	ml_dst = result->mloop;
	
	/* copy the faces across, remapping indices */
	k = 0;
	for (i = 0; i < numFaces_dst; i++) {
		MPoly *source;
		MPoly *dest;

		source = mpoly_src + faceMap[i];
		dest = mpoly_dst + i;
		CustomData_copy_data(&mesh->pdata, &result->pdata, faceMap[i], i, 1);

		*dest = *source;
		dest->loopstart = k;
		CustomData_copy_data(&mesh->ldata, &result->ldata, source->loopstart, dest->loopstart, dest->totloop);

		ml_src = mloop_src + source->loopstart;
		for (j = 0; j < source->totloop; j++, k++, ml_src++, ml_dst++) {
			ml_dst->v = GET_INT_FROM_POINTER(BLI_ghash_lookup(vertHash, SET_INT_IN_POINTER(ml_src->v)));
			ml_dst->e = GET_INT_FROM_POINTER(BLI_ghash_lookup(edgeHash, SET_INT_IN_POINTER(ml_src->e)));
		}
	}

	BLI_ghash_free(vertHash, NULL, NULL);
	BLI_ghash_free(edgeHash, NULL, NULL);
	
	MEM_freeN(vertMap);
	MEM_freeN(edgeMap);
	MEM_freeN(faceMap);

	BKE_mesh_tessface_ensure(result);

	/* TODO(sybren): also copy flags & tags? */
	return result;
}


ModifierTypeInfo modifierType_Build = {
	/* name */              "Build",
	/* structName */        "BuildModifierData",
	/* structSize */        sizeof(BuildModifierData),
	/* type */              eModifierTypeType_Nonconstructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_AcceptsCVs,
	/* copyData */          copyData,

	/* deformVerts_DM */    NULL,
	/* deformMatrices_DM */ NULL,
	/* deformVertsEM_DM */  NULL,
	/* deformMatricesEM_DM*/NULL,
	/* applyModifier_DM */  NULL,
	/* applyModifierEM_DM */NULL,

	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,

	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepsgraph */   NULL,
	/* dependsOnTime */     dependsOnTime,
	/* dependsOnNormals */  NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
