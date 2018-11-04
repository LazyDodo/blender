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
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/hair/edithair_test.c
 *  \ingroup edhair
 */

#include "MEM_guardedalloc.h"

#include "BLI_math.h"

#include "DNA_hair_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_context.h"
#include "BKE_hair.h"
#include "BKE_mesh_sample.h"

#include "DEG_depsgraph.h"

#include "ED_hair.h"
#include "ED_screen.h"
#include "ED_view3d.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "WM_api.h"
#include "WM_types.h"

#include "UI_resources.h"

#include "BLT_translation.h"

#include "hair_intern.h"  /* own include */

/* -------------------------------------------------------------------- */
/** \name Add test hair
 * Test operator for quickly adding hair
 * \{ */


/* Distribute hair follicles on a scalp mesh.
 * Optional per-loop weights control follicle density on the scalp.
 */
static void hair_generate_follicles_ex(
        HairPattern *pattern,
        struct Mesh *scalp,
        unsigned int seed,
        int count,
        const float *loop_weights)
{
	// Limit max_count to theoretical limit based on area
	float scalp_area = BKE_hair_calc_surface_area(scalp);
	float density = BKE_hair_calc_density_from_count(scalp_area, count);
	float min_distance = BKE_hair_calc_min_distance_from_density(density);

	if (pattern->follicles)
	{
		MEM_freeN(pattern->follicles);
	}
	pattern->follicles = MEM_callocN(sizeof(HairFollicle) * count, "hair follicles");

	{
		MeshSampleGenerator *gen = BKE_mesh_sample_gen_surface_poissondisk(seed, min_distance, count, loop_weights);

		BKE_mesh_sample_generator_bind(gen, scalp);

		static const bool use_threads = false;
		pattern->num_follicles = BKE_mesh_sample_generate_batch_ex(
		            gen,
		            &pattern->follicles->mesh_sample,
		            sizeof(HairFollicle),
		            count,
		            use_threads);

		BKE_mesh_sample_free_generator(gen);
	}

	for (int i = 0; i < pattern->num_follicles; ++i) {
		HairFollicle *follicle = &pattern->follicles[i];

		follicle->curve = HAIR_CURVE_INDEX_NONE;
	}
}

static int add_test_hair_exec(bContext *C, wmOperator *op)
{
	struct Depsgraph *depsgraph = CTX_data_depsgraph(C);
	Object *obedit = CTX_data_edit_object(C);;
	HairSystem *hsys = obedit->data;
	EditHair *edit = hsys->edithair;
	struct Mesh *scalp = BKE_hair_get_scalp(depsgraph, obedit, hsys);
	if (!scalp)
	{
		return OPERATOR_CANCELLED;
	}

	const int seed = RNA_int_get(op->ptr, "seed");
	const int count = RNA_int_get(op->ptr, "count");

	hair_generate_follicles_ex(edit->pattern, scalp, seed, count, NULL);

	BKE_hair_batch_cache_dirty(hsys, BKE_HAIR_BATCH_DIRTY_ALL);
	DEG_id_tag_update(obedit->data, DEG_TAG_SELECT_UPDATE);
	WM_event_add_notifier(C, NC_GEOM|ND_DRAW, obedit);

	return OPERATOR_FINISHED;
}

void HAIR_OT_add_test_hair(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Add Test Hair";
	ot->description = "Test hair distribution";
	ot->idname = "HAIR_OT_add_test_hair";

	/* api callbacks */
	ot->invoke = WM_operator_props_popup_confirm;
	ot->exec = add_test_hair_exec;
	ot->poll = ED_operator_edithair;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	RNA_def_int(ot->srna, "seed", 0, 0, 1000000, "Seed", "Random seed value", 0, 1000000);
	RNA_def_int(ot->srna, "count", 100, 1, 1000000, "Count", "Number of hairs to add", 1, 10000);
}

/** \} */
