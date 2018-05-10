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
 * The Original Code is Copyright (C) 2017, Blender Foundation
 * This is a new part of Blender
 *
 * Contributor(s): Antonio Vazquez
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_gpencilcolor.c
 *  \ingroup modifiers
 */

#include <stdio.h>

#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_gpencil_types.h"

#include "BLI_blenlib.h"
#include "BLI_ghash.h"
#include "BLI_math_color.h"
#include "BLI_math_vector.h"
#include "BLI_utildefines.h"

#include "BKE_global.h"
#include "BKE_context.h"
#include "BKE_gpencil.h"
#include "BKE_main.h"
#include "BKE_material.h"

#include "DEG_depsgraph.h"

#include "MOD_modifiertypes.h"
#include "MOD_gpencil_util.h"

static void initData(ModifierData *md)
{
	GpencilColorModifierData *gpmd = (GpencilColorModifierData *)md;
	gpmd->pass_index = 0;
	ARRAY_SET_ITEMS(gpmd->hsv, 1.0f, 1.0f, 1.0f);
	gpmd->layername[0] = '\0';
	gpmd->flag |= GP_COLOR_CREATE_COLORS;
}

static void copyData(const ModifierData *md, ModifierData *target)
{
	modifier_copyData_generic(md, target);
}

/* color correction strokes */
static void deformStroke(ModifierData *md, Depsgraph *UNUSED(depsgraph),
                         Object *ob, bGPDlayer *gpl, bGPDstroke *gps)
{

	GpencilColorModifierData *mmd = (GpencilColorModifierData *)md;
	float hsv[3], factor[3];

	if (!is_stroke_affected_by_modifier(ob,
	        mmd->layername, mmd->pass_index, 1, gpl, gps,
	        mmd->flag & GP_COLOR_INVERSE_LAYER, mmd->flag & GP_COLOR_INVERSE_PASS))
	{
		return;
	}
	
	copy_v3_v3(factor, mmd->hsv);
	add_v3_fl(factor, -1.0f);

	rgb_to_hsv_v(gps->tmp_rgb, hsv);
	add_v3_v3(hsv, factor);
	CLAMP3(hsv, 0.0f, 1.0f);
	hsv_to_rgb_v(hsv, gps->tmp_rgb);

	rgb_to_hsv_v(gps->tmp_fill, hsv);
	add_v3_v3(hsv, factor);
	CLAMP3(hsv, 0.0f, 1.0f);
	hsv_to_rgb_v(hsv, gps->tmp_fill);
}

static void bakeModifierGP(const bContext *C, Depsgraph *depsgraph,
                           ModifierData *md, Object *ob)
{
	GpencilColorModifierData *mmd = (GpencilColorModifierData *)md;
	Main *bmain = CTX_data_main(C);
	bGPdata *gpd = ob->data;
	
	GHash *gh_color = BLI_ghash_str_new("GP_Color modifier");
	for (bGPDlayer *gpl = gpd->layers.first; gpl; gpl = gpl->next) {
		for (bGPDframe *gpf = gpl->frames.first; gpf; gpf = gpf->next) {
			for (bGPDstroke *gps = gpf->strokes.first; gps; gps = gps->next) {

				Material *mat = give_current_material(ob, gps->mat_nr + 1);
				if (mat == NULL)
					continue;
				GpencilColorData *gpcolor = mat->gpcolor;
				/* skip stroke if it doesn't have color info */
				if (ELEM(NULL, gpcolor))
					continue;

				copy_v4_v4(gps->tmp_rgb, gpcolor->rgb);
				copy_v4_v4(gps->tmp_fill, gpcolor->fill);

				/* look for color */
				if (mmd->flag & GP_TINT_CREATE_COLORS) {
					Material *newmat = (Material *)BLI_ghash_lookup(gh_color, mat->id.name);
					if (newmat == NULL) {
						BKE_object_material_slot_add(ob);
						newmat = BKE_material_copy(bmain, mat);
						assign_material(ob, newmat, ob->totcol, BKE_MAT_ASSIGN_EXISTING);

						copy_v4_v4(newmat->gpcolor->rgb, gps->tmp_rgb);
						copy_v4_v4(newmat->gpcolor->fill, gps->tmp_fill);

						BLI_ghash_insert(gh_color, mat->id.name, newmat);
					}
					/* reasign color index */
					int idx = BKE_object_material_slot_find_index(ob, newmat);
					gps->mat_nr = idx - 1;
				}
				else {
					/* reuse existing color */
					copy_v4_v4(gpcolor->rgb, gps->tmp_rgb);
					copy_v4_v4(gpcolor->fill, gps->tmp_fill);
				}

				deformStroke(md, depsgraph, ob, gpl, gps);
			}
		}
	}
	/* free hash buffers */
	if (gh_color) {
		BLI_ghash_free(gh_color, NULL, NULL);
		gh_color = NULL;
	}
}

ModifierTypeInfo modifierType_Gpencil_Color = {
	/* name */              "Hue/Saturation",
	/* structName */        "GpencilColorModifierData",
	/* structSize */        sizeof(GpencilColorModifierData),
	/* type */              eModifierTypeType_Gpencil,
	/* flags */             eModifierTypeFlag_GpencilMod | eModifierTypeFlag_SupportsEditmode,

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
	/* applyModifier */     NULL,
	/* applyModifierEM */   NULL,

	/* deformStroke */      deformStroke,
	/* generateStrokes */   NULL,
	/* bakeModifierGP */    bakeModifierGP,

	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepsgraph */   NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
