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

/** \file blender/editors/hair/edithair_brush.c
 *  \ingroup edhair
 */

#include "MEM_guardedalloc.h"

#include "BLI_math.h"

#include "DNA_hair_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_context.h"
#include "BKE_customdata.h"
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

// typedef struct HairOperator {
// } HairOperator;

/************************* Brush edit operator ********************/

typedef struct BrushEdit {
	Scene *scene;
	ViewLayer *view_layer;
	Object *ob;
	EditHair *edit;

	bool first;
	float lastmouse[2];
	float zfac;

	ViewContext vc;
} BrushEdit;

static bool brush_edit_init(bContext *C, wmOperator *op)
{
	Scene *scene = CTX_data_scene(C);
	ViewLayer *view_layer = CTX_data_view_layer(C);
	ARegion *ar = CTX_wm_region(C);
	const HairEditSettings *hset = &scene->toolsettings->hair_edit;
	Object *ob = CTX_data_active_object(C);
	HairSystem *hsys = ob->data;
	EditHair *edit = hsys->edithair;
	
	if (hset->brushtype < 0)
		return false;

	BrushEdit *bedit;
	bedit = MEM_callocN(sizeof(BrushEdit), "BrushEdit");
	op->customdata = bedit;

	bedit->scene = scene;
	bedit->view_layer = view_layer;
	bedit->ob = ob;
	bedit->edit = edit;
	bedit->first = true;

	/* set the 'distance factor' for grabbing (used in comb etc) */
	float min[3], max[3];
	INIT_MINMAX(min, max);
	// TODO implement selection mechanism for affected hair points/scalp area
	BKE_hair_minmax(hsys, min, max);
	mid_v3_v3v3(min, min, max);
	bedit->zfac = ED_view3d_calc_zfac(ar->regiondata, min, NULL);

	ED_hair_init_view3d(C, &bedit->vc);

	return true;
}

static void brush_edit_exit(wmOperator *op)
{
	BrushEdit *bedit = op->customdata;

	MEM_freeN(bedit);
}

static void brush_edit_apply(bContext *C, wmOperator *op, PointerRNA *itemptr)
{
	Depsgraph *depsgraph = CTX_data_depsgraph(C);
	BrushEdit *bedit = op->customdata;

	float mouse[2];
	RNA_float_get_array(itemptr, "mouse", mouse);
	bool flip = RNA_boolean_get(itemptr, "pen_flip");

	if (bedit->first) {
		bedit->lastmouse[0] = mouse[0];
		bedit->lastmouse[1] = mouse[1];
	}

#if 0
	BrushEdit *bedit = op->customdata;
	Depsgraph *depsgraph = CTX_data_depsgraph(C);
	Scene *scene = bedit->scene;
	Object *ob = bedit->ob;
	PTCacheEdit *edit = bedit->edit;
	ParticleEditSettings *pset = PE_settings(scene);
	ParticleSystemModifierData *psmd_eval = edit->psmd_eval;
	ParticleBrushData *brush = &pset->brush[pset->brushtype];
	ARegion *ar = CTX_wm_region(C);
	float vec[3], mousef[2];
	int mval[2];
	int flip, mouse[2], removed = 0, added = 0, selected = 0, tot_steps = 1, step = 1;
	float dx, dy, dmax;
	int lock_root = pset->flag & PE_LOCK_FIRST;

	if (!PE_start_edit(edit))
		return;

	RNA_float_get_array(itemptr, "mouse", mousef);
	mouse[0] = mousef[0];
	mouse[1] = mousef[1];
	flip = RNA_boolean_get(itemptr, "pen_flip");

	if (bedit->first) {
		bedit->lastmouse[0] = mouse[0];
		bedit->lastmouse[1] = mouse[1];
	}

	dx = mouse[0] - bedit->lastmouse[0];
	dy = mouse[1] - bedit->lastmouse[1];

	mval[0] = mouse[0];
	mval[1] = mouse[1];


	/* disable locking temporatily for disconnected hair */
	if (edit->psys && edit->psys->flag & PSYS_GLOBAL_HAIR)
		pset->flag &= ~PE_LOCK_FIRST;

	if (((pset->brushtype == PE_BRUSH_ADD) ?
	     (sqrtf(dx * dx + dy * dy) > pset->brush[PE_BRUSH_ADD].step) : (dx != 0 || dy != 0)) || bedit->first)
	{
		PEData data = bedit->data;
		data.context = C; // TODO(mai): why isnt this set in bedit->data?

		view3d_operator_needs_opengl(C);
		selected = (short)count_selected_keys(scene, edit);

		dmax = max_ff(fabsf(dx), fabsf(dy));
		tot_steps = dmax / (0.2f * pe_brush_size_get(scene, brush)) + 1;

		dx /= (float)tot_steps;
		dy /= (float)tot_steps;

		for (step = 1; step <= tot_steps; step++) {
			mval[0] = bedit->lastmouse[0] + step * dx;
			mval[1] = bedit->lastmouse[1] + step * dy;

			switch (pset->brushtype) {
				case PE_BRUSH_COMB:
				{
					const float mval_f[2] = {dx, dy};
					data.mval = mval;
					data.rad = pe_brush_size_get(scene, brush);

					data.combfac = (brush->strength - 0.5f) * 2.0f;
					if (data.combfac < 0.0f)
						data.combfac = 1.0f - 9.0f * data.combfac;
					else
						data.combfac = 1.0f - data.combfac;

					invert_m4_m4(ob->imat, ob->obmat);

					ED_view3d_win_to_delta(ar, mval_f, vec, bedit->zfac);
					data.dvec = vec;

					foreach_mouse_hit_key(&data, brush_comb, selected);
					break;
				}
				case PE_BRUSH_CUT:
				{
					if (edit->psys && edit->pathcache) {
						data.mval = mval;
						data.rad = pe_brush_size_get(scene, brush);
						data.cutfac = brush->strength;

						if (selected)
							foreach_selected_point(&data, brush_cut);
						else
							foreach_point(&data, brush_cut);

						removed = remove_tagged_particles(ob, edit->psys, pe_x_mirror(ob));
						if (pset->flag & PE_KEEP_LENGTHS)
							recalc_lengths(edit);
					}
					else
						removed = 0;

					break;
				}
				case PE_BRUSH_LENGTH:
				{
					data.mval = mval;

					data.rad = pe_brush_size_get(scene, brush);
					data.growfac = brush->strength / 50.0f;

					if (brush->invert ^ flip)
						data.growfac = 1.0f - data.growfac;
					else
						data.growfac = 1.0f + data.growfac;

					foreach_mouse_hit_point(&data, brush_length, selected);

					if (pset->flag & PE_KEEP_LENGTHS)
						recalc_lengths(edit);
					break;
				}
				case PE_BRUSH_PUFF:
				{
					if (edit->psys) {
						data.mesh = psmd_eval->mesh_final;
						data.mval = mval;
						data.rad = pe_brush_size_get(scene, brush);
						data.select = selected;

						data.pufffac = (brush->strength - 0.5f) * 2.0f;
						if (data.pufffac < 0.0f)
							data.pufffac = 1.0f - 9.0f * data.pufffac;
						else
							data.pufffac = 1.0f - data.pufffac;

						data.invert = (brush->invert ^ flip);
						invert_m4_m4(ob->imat, ob->obmat);

						foreach_mouse_hit_point(&data, brush_puff, selected);
					}
					break;
				}
				case PE_BRUSH_ADD:
				{
					if (edit->psys && edit->psys->part->from == PART_FROM_FACE) {
						data.mval = mval;

						added = brush_add(C, &data, brush->count);

						if (pset->flag & PE_KEEP_LENGTHS)
							recalc_lengths(edit);
					}
					else
						added = 0;
					break;
				}
				case PE_BRUSH_SMOOTH:
				{
					data.mval = mval;
					data.rad = pe_brush_size_get(scene, brush);

					data.vec[0] = data.vec[1] = data.vec[2] = 0.0f;
					data.tot = 0;

					data.smoothfac = brush->strength;

					invert_m4_m4(ob->imat, ob->obmat);

					foreach_mouse_hit_key(&data, brush_smooth_get, selected);

					if (data.tot) {
						mul_v3_fl(data.vec, 1.0f / (float)data.tot);
						foreach_mouse_hit_key(&data, brush_smooth_do, selected);
					}

					break;
				}
				case PE_BRUSH_WEIGHT:
				{
					if (edit->psys) {
						data.mesh = psmd_eval->mesh_final;
						data.mval = mval;
						data.rad = pe_brush_size_get(scene, brush);

						data.weightfac = brush->strength; /* note that this will never be zero */

						foreach_mouse_hit_key(&data, BKE_brush_weight_get, selected);
					}

					break;
				}
			}
			if ((pset->flag & PE_KEEP_LENGTHS) == 0)
				recalc_lengths(edit);

			if (ELEM(pset->brushtype, PE_BRUSH_ADD, PE_BRUSH_CUT) && (added || removed)) {
				if (pset->brushtype == PE_BRUSH_ADD && pe_x_mirror(ob))
					PE_mirror_x(scene, ob, 1);

				update_world_cos(depsgraph, ob, edit);
				psys_free_path_cache(NULL, edit);
				DEG_id_tag_update(&ob->id, ID_RECALC_GEOMETRY);
			}
			else {
				PE_update_object(depsgraph, scene, ob, 1);
			}
		}

		if (edit->psys) {
			WM_event_add_notifier(C, NC_OBJECT | ND_PARTICLE | NA_EDITED, ob);
			BKE_particle_batch_cache_dirty_tag(edit->psys, BKE_PARTICLE_BATCH_DIRTY_ALL);
			DEG_id_tag_update(&ob->id, ID_RECALC_SELECT);
		}
		else {
			DEG_id_tag_update(&ob->id, ID_RECALC_GEOMETRY);
			WM_event_add_notifier(C, NC_OBJECT | ND_MODIFIER, ob);
		}

		bedit->lastmouse[0] = mouse[0];
		bedit->lastmouse[1] = mouse[1];
		bedit->first = 0;
	}

	pset->flag |= lock_root;
#endif
}

static int brush_edit_exec(bContext *C, wmOperator *op)
{
	if (!brush_edit_init(C, op))
		return OPERATOR_CANCELLED;

	RNA_BEGIN(op->ptr, itemptr, "stroke")
	{
		brush_edit_apply(C, op, &itemptr);
	}
	RNA_END;

	brush_edit_exit(op);

	return OPERATOR_FINISHED;
}

static void brush_edit_apply_event(bContext *C, wmOperator *op, const wmEvent *event)
{
	PointerRNA itemptr;
	float mouse[2];

	VECCOPY2D(mouse, event->mval);

	/* fill in stroke */
	RNA_collection_add(op->ptr, "stroke", &itemptr);

	RNA_float_set_array(&itemptr, "mouse", mouse);
	RNA_boolean_set(&itemptr, "pen_flip", event->shift != false); // XXX hardcoded

	/* apply */
	brush_edit_apply(C, op, &itemptr);
}

static int brush_edit_invoke(bContext *C, wmOperator *op, const wmEvent *event)
{
	if (!brush_edit_init(C, op))
		return OPERATOR_CANCELLED;

	brush_edit_apply_event(C, op, event);

	WM_event_add_modal_handler(C, op);

	return OPERATOR_RUNNING_MODAL;
}

static int brush_edit_modal(bContext *C, wmOperator *op, const wmEvent *event)
{
	switch (event->type) {
		case LEFTMOUSE:
		case MIDDLEMOUSE:
		case RIGHTMOUSE: // XXX hardcoded
			if (event->val == KM_RELEASE) {
				brush_edit_exit(op);
				return OPERATOR_FINISHED;
			}
			break;
		case MOUSEMOVE:
			brush_edit_apply_event(C, op, event);
			break;
	}

	return OPERATOR_RUNNING_MODAL;
}

static void brush_edit_cancel(bContext *UNUSED(C), wmOperator *op)
{
	brush_edit_exit(op);
}

void HAIR_OT_brush_edit(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Brush Edit";
	ot->idname = "HAIR_OT_brush_edit";
	ot->description = "Apply a brush stroke to the hair";

	/* api callbacks */
	ot->exec = brush_edit_exec;
	ot->invoke = brush_edit_invoke;
	ot->modal = brush_edit_modal;
	ot->cancel = brush_edit_cancel;
	ot->poll = ED_hair_poll_view3d;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO | OPTYPE_BLOCKING;

	/* properties */
	PropertyRNA *prop;
	prop = RNA_def_collection_runtime(ot->srna, "stroke", &RNA_OperatorStrokeElement, "Stroke", "");
	RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);
}
