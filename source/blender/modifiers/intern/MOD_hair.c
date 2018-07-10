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

/** \file blender/modifiers/intern/MOD_hair.c
 *  \ingroup modifiers
 */

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"
#include "BLI_listbase.h"

#include "DNA_object_types.h"
#include "DNA_hair_types.h"

#include "BKE_hair.h"
#include "BKE_library.h"
#include "BKE_library_query.h"
#include "BKE_modifier.h"

#include "DEG_depsgraph_build.h"

#include "MOD_util.h"


static void initData(ModifierData *md)
{
	HairModifierData *hmd = (HairModifierData *) md;
	
	hmd->hair_system = BKE_hair_new();
	
	hmd->flag |= 0;
	
	hmd->follicle_count = 100000;
	
	hmd->draw_settings = BKE_hair_draw_settings_new();
}

static void copyData(const ModifierData *md, ModifierData *target, int flag)
{
	const HairModifierData *hmd = (HairModifierData *) md;
	HairModifierData *tfmd = (HairModifierData *) target;

    modifier_copyData_generic(md, target, flag);
	
	if (hmd->hair_system) {
		tfmd->hair_system = BKE_hair_copy(hmd->hair_system);
	}
	if (hmd->draw_settings)
	{
		tfmd->draw_settings = BKE_hair_draw_settings_copy(hmd->draw_settings);
	}
}

static void freeData(ModifierData *md)
{
	HairModifierData *hmd = (HairModifierData *) md;
	
	if (hmd->hair_system) {
		BKE_hair_free(hmd->hair_system);
	}
	if (hmd->draw_settings)
	{
		BKE_hair_draw_settings_free(hmd->draw_settings);
	}
	for (HairModifierFiberCurve *curve = hmd->fiber_curves.first; curve; curve = curve->next)
	{
		if (curve->verts)
		{
			MEM_freeN(curve->verts);
		}
	}
	BLI_freelistN(&hmd->fiber_curves);
}

static struct Mesh *applyModifier(ModifierData *md, const ModifierEvalContext *ctx,
                                  struct Mesh *mesh)
{
	HairModifierData *hmd = (HairModifierData *) md;
	
	UNUSED_VARS(hmd, ctx);
	
	return mesh;
}

static void foreachObjectLink(
        ModifierData *md,
        Object *ob,
        ObjectWalkFunc walk,
        void *userData)
{
	HairModifierData *hmd = (HairModifierData *) md;
	UNUSED_VARS(ob, walk, userData, hmd);
}

static void foreachIDLink(
        ModifierData *md,
        Object *ob,
        IDWalkFunc walk,
        void *userData)
{
	HairModifierData *hmd = (HairModifierData *) md;
	UNUSED_VARS(hmd);
	
	foreachObjectLink(md, ob, (ObjectWalkFunc)walk, userData);
}

ModifierTypeInfo modifierType_Hair = {
	/* name */              "Hair",
	/* structName */        "HairModifierData",
	/* structSize */        sizeof(HairModifierData),
	/* type */              eModifierTypeType_NonGeometrical,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsEditmode,

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
	/* freeData */          freeData,
	/* isDisabled */        NULL,
	/* updateDepsgraph */   NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     foreachIDLink,
	/* foreachTexLink */    NULL,
};
