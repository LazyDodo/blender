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

/** \file blender/modifiers/intern/MOD_boolean.c
 *  \ingroup modifiers
 */

// #ifdef DEBUG_TIME

#include <stdio.h>

#include "DNA_object_types.h"

#include "BLI_utildefines.h"
#include "BLI_math_matrix.h"

#include "BKE_library_query.h"
#include "BKE_modifier.h"

#include "MOD_util.h"

#include "BLI_alloca.h"
#include "BLI_math_geom.h"

#include "BKE_global.h"  /* only to check G.debug */
#include "BKE_library.h"
#include "BKE_material.h"
#include "BKE_mesh.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "DEG_depsgraph_query.h"

#include "MEM_guardedalloc.h"
#include "BKE_boolean.h"

#include "bmesh.h"
#include "bmesh_tools.h"
#include "tools/bmesh_intersect.h"

#ifdef DEBUG_TIME
#  include "PIL_time.h"
#  include "PIL_time_utildefines.h"
#endif

static void initData(ModifierData *md)
{
	BooleanModifierData *bmd = (BooleanModifierData *)md;

	bmd->double_threshold = 1e-6f;
}

static bool isDisabled(const struct Scene *UNUSED(scene), ModifierData *md, bool UNUSED(useRenderParams))
{
	BooleanModifierData *bmd = (BooleanModifierData *) md;

	return !bmd->object;
}

static void foreachObjectLink(
        ModifierData *md, Object *ob,
        ObjectWalkFunc walk, void *userData)
{
	BooleanModifierData *bmd = (BooleanModifierData *) md;

	walk(userData, ob, &bmd->object, IDWALK_CB_NOP);
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
	BooleanModifierData *bmd = (BooleanModifierData *)md;
	if (bmd->object != NULL) {
		DEG_add_object_relation(ctx->node, bmd->object, DEG_OB_COMP_TRANSFORM, "Boolean Modifier");
		DEG_add_object_relation(ctx->node, bmd->object, DEG_OB_COMP_GEOMETRY, "Boolean Modifier");
	}
	/* We need own transformation as well. */
	DEG_add_object_relation(ctx->node, ctx->object, DEG_OB_COMP_TRANSFORM, "Boolean Modifier");
}


static Mesh *applyModifier(ModifierData *md, const ModifierEvalContext *ctx, Mesh *mesh)
{
	BooleanModifierData *bmd = (BooleanModifierData *) md;
	Mesh *result = mesh;

	Mesh *mesh_other;
	bool mesh_other_free;

	if (bmd->object == NULL) {
		return result;
	}

    Object *other = DEG_get_evaluated_object(ctx->depsgraph, bmd->object);
    mesh_other = BKE_modifier_get_evaluated_mesh_from_evaluated_object(other, &mesh_other_free);

    result = BKE_boolean_operation(mesh, ctx->object, mesh_other, bmd->object, bmd->operation,
                                      bmd->double_threshold, bmd);

    /* if new mesh returned, return it; otherwise there was
     * an error, so delete the modifier object */
    if (result == NULL) {
        modifier_setError(md, "Cannot execute boolean operation");
    }

	if (mesh_other != NULL && mesh_other_free) {
		BKE_id_free(NULL, mesh_other);
	}

	return result;
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *UNUSED(md))
{
	CustomDataMask dataMask = CD_MASK_MTFACE | CD_MASK_MEDGE;

	dataMask |= CD_MASK_MDEFORMVERT;

	return dataMask;
}

ModifierTypeInfo modifierType_Boolean = {
	/* name */              "Boolean",
	/* structName */        "BooleanModifierData",
	/* structSize */        sizeof(BooleanModifierData),
	/* type */              eModifierTypeType_Nonconstructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_UsesPointCache,

	/* copyData */          modifier_copyData_generic,

	/* deformVerts_DM */    NULL,
	/* deformMatrices_DM */ NULL,
	/* deformVertsEM_DM */  NULL,
	/* deformMatricesEM_DM*/NULL,
	/* applyModifier_DM */  NULL,

	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,

	/* initData */          initData,
	/* requiredDataMask */  requiredDataMask,
	/* freeData */          NULL,
	/* isDisabled */        isDisabled,
	/* updateDepsgraph */   updateDepsgraph,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  NULL,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
