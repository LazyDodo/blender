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
 * The Original Code is Copyright (C) 2018, Blender Foundation
 * This is a new part of Blender
 *
 * Contributor(s): Antonio Vazquez
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/gpencil_modifiers/intern/MOD_gpencilarmature.c
 *  \ingroup modifiers
 */

#include <stdio.h>

#include "DNA_armature_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_gpencil_modifier_types.h"
#include "DNA_modifier_types.h"
#include "BLI_math.h"

#include "BLI_listbase.h"
#include "BLI_task.h"
#include "BLI_utildefines.h"

#include "BKE_action.h"
#include "BKE_armature.h"
#include "BKE_context.h"
#include "BKE_colortools.h"
#include "BKE_deform.h"
#include "BKE_gpencil.h"
#include "BKE_gpencil_modifier.h"
#include "BKE_modifier.h"
#include "BKE_library_query.h"
#include "BKE_scene.h"
#include "BKE_main.h"
#include "BKE_layer.h"

#include "MEM_guardedalloc.h"

#include "MOD_gpencil_util.h"
#include "MOD_gpencil_modifiertypes.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"
#include "DEG_depsgraph_query.h"

typedef struct bPoseChanDeform {
	Mat4     *b_bone_mats;
	DualQuat *dual_quat;
	DualQuat *b_bone_dual_quats;
} bPoseChanDeform;

typedef struct ArmatureBBoneDefmatsData {
	bPoseChanDeform *pdef_info_array;
	DualQuat *dualquats;
	bool use_quaternion;
} ArmatureBBoneDefmatsData;

static void initData(GpencilModifierData *md)
{
	ArmatureGpencilModifierData *gpmd = (ArmatureGpencilModifierData *)md;
	gpmd->object = NULL;
}

static void copyData(const GpencilModifierData *md, GpencilModifierData *target)
{
	BKE_gpencil_modifier_copyData_generic(md, target);
}

static void b_bone_deform(bPoseChanDeform *pdef_info, Bone *bone, float co[3], DualQuat *dq, float defmat[3][3])
{
	Mat4 *b_bone = pdef_info->b_bone_mats;
	float(*mat)[4] = b_bone[0].mat;
	float segment, y;
	int a;

	/* need to transform co back to bonespace, only need y */
	y = mat[0][1] * co[0] + mat[1][1] * co[1] + mat[2][1] * co[2] + mat[3][1];

	/* now calculate which of the b_bones are deforming this */
	segment = bone->length / ((float)bone->segments);
	a = (int)(y / segment);

	/* note; by clamping it extends deform at endpoints, goes best with
	 * straight joints in restpos. */
	CLAMP(a, 0, bone->segments - 1);

	if (dq) {
		copy_dq_dq(dq, &(pdef_info->b_bone_dual_quats)[a]);
	}
	else {
		mul_m4_v3(b_bone[a + 1].mat, co);

		if (defmat) {
			copy_m3_m4(defmat, b_bone[a + 1].mat);
		}
	}
}

static void pchan_deform_mat_add(bPoseChannel *pchan, float weight, float bbonemat[3][3], float mat[3][3])
{
	float wmat[3][3];

	if (pchan->bone->segments > 1)
		copy_m3_m3(wmat, bbonemat);
	else
		copy_m3_m4(wmat, pchan->chan_mat);

	mul_m3_fl(wmat, weight);
	add_m3_m3m3(mat, mat, wmat);
}

static void pchan_bone_deform(bPoseChannel *pchan, bPoseChanDeform *pdef_info, float weight, float vec[3], DualQuat *dq,
	float mat[3][3], const float co[3], float *contrib)
{
	float cop[3], bbonemat[3][3];
	DualQuat bbonedq;

	if (!weight)
		return;

	copy_v3_v3(cop, co);

	if (vec) {
		if (pchan->bone->segments > 1)
			/* applies on cop and bbonemat */
			b_bone_deform(pdef_info, pchan->bone, cop, NULL, (mat) ? bbonemat : NULL);
		else
			mul_m4_v3(pchan->chan_mat, cop);

		vec[0] += (cop[0] - co[0]) * weight;
		vec[1] += (cop[1] - co[1]) * weight;
		vec[2] += (cop[2] - co[2]) * weight;

		if (mat)
			pchan_deform_mat_add(pchan, weight, bbonemat, mat);
	}
	else {
		if (pchan->bone->segments > 1) {
			b_bone_deform(pdef_info, pchan->bone, cop, &bbonedq, NULL);
			add_weighted_dq_dq(dq, &bbonedq, weight);
		}
		else
			add_weighted_dq_dq(dq, pdef_info->dual_quat, weight);
	}

	(*contrib) += weight;
}

static float dist_bone_deform(bPoseChannel *pchan, bPoseChanDeform *pdef_info, float vec[3], DualQuat *dq,
	float mat[3][3], const float co[3])
{
	Bone *bone = pchan->bone;
	float fac, contrib = 0.0;
	float cop[3], bbonemat[3][3];
	DualQuat bbonedq;

	if (bone == NULL)
		return 0.0f;

	copy_v3_v3(cop, co);

	fac = distfactor_to_bone(cop, bone->arm_head, bone->arm_tail, bone->rad_head, bone->rad_tail, bone->dist);

	if (fac > 0.0f) {
		fac *= bone->weight;
		contrib = fac;
		if (contrib > 0.0f) {
			if (vec) {
				if (bone->segments > 1)
					/* applies on cop and bbonemat */
					b_bone_deform(pdef_info, bone, cop, NULL, (mat) ? bbonemat : NULL);
				else
					mul_m4_v3(pchan->chan_mat, cop);

				/* Make this a delta from the base position */
				sub_v3_v3(cop, co);
				madd_v3_v3fl(vec, cop, fac);

				if (mat)
					pchan_deform_mat_add(pchan, fac, bbonemat, mat);
			}
			else {
				if (bone->segments > 1) {
					b_bone_deform(pdef_info, bone, cop, &bbonedq, NULL);
					add_weighted_dq_dq(dq, &bbonedq, fac);
				}
				else
					add_weighted_dq_dq(dq, pdef_info->dual_quat, fac);
			}
		}
	}

	return contrib;
}

static void pchan_b_bone_defmats(bPoseChannel *pchan, bPoseChanDeform *pdef_info, const bool use_quaternion)
{
	Bone *bone = pchan->bone;
	Mat4 b_bone[MAX_BBONE_SUBDIV], b_bone_rest[MAX_BBONE_SUBDIV];
	Mat4 *b_bone_mats;
	DualQuat *b_bone_dual_quats = NULL;
	int a;

	b_bone_spline_setup(pchan, 0, b_bone);
	b_bone_spline_setup(pchan, 1, b_bone_rest);

	/* allocate b_bone matrices and dual quats */
	b_bone_mats = MEM_mallocN((1 + bone->segments) * sizeof(Mat4), "BBone defmats");
	pdef_info->b_bone_mats = b_bone_mats;

	if (use_quaternion) {
		b_bone_dual_quats = MEM_mallocN((bone->segments) * sizeof(DualQuat), "BBone dqs");
		pdef_info->b_bone_dual_quats = b_bone_dual_quats;
	}

	/* first matrix is the inverse arm_mat, to bring points in local bone space
	 * for finding out which segment it belongs to */
	invert_m4_m4(b_bone_mats[0].mat, bone->arm_mat);

	/* then we make the b_bone_mats:
	 * - first transform to local bone space
	 * - translate over the curve to the bbone mat space
	 * - transform with b_bone matrix
	 * - transform back into global space */

	for (a = 0; a < bone->segments; a++) {
		float tmat[4][4];

		invert_m4_m4(tmat, b_bone_rest[a].mat);
		mul_m4_series(b_bone_mats[a + 1].mat, pchan->chan_mat, bone->arm_mat, b_bone[a].mat, tmat, b_bone_mats[0].mat);

		if (use_quaternion)
			mat4_to_dquat(&b_bone_dual_quats[a], bone->arm_mat, b_bone_mats[a + 1].mat);
	}
}

static void gpencil_armature_bbone_defmats_cb(void *userdata, Link *iter, int index)
{
	ArmatureBBoneDefmatsData *data = userdata;
	bPoseChannel *pchan = (bPoseChannel *)iter;

	if (!(pchan->bone->flag & BONE_NO_DEFORM)) {
		bPoseChanDeform *pdef_info = &data->pdef_info_array[index];
		const bool use_quaternion = data->use_quaternion;

		if (pchan->bone->segments > 1) {
			pchan_b_bone_defmats(pchan, pdef_info, use_quaternion);
		}

		if (use_quaternion) {
			pdef_info->dual_quat = &data->dualquats[index];
			mat4_to_dquat(pdef_info->dual_quat, pchan->bone->arm_mat, pchan->chan_mat);
		}
	}
}

static void gpencil_armature_deform_verts(Object *armOb, Object *target, bGPDstroke *gps, float(*vertexCos)[3],
	float(*defMats)[3][3], int numVerts, int deformflag,
	float(*prevCos)[3], const char *defgrp_name)
{
	bPoseChanDeform *pdef_info_array;
	bPoseChanDeform *pdef_info = NULL;
	bArmature *arm = armOb->data;
	bPoseChannel *pchan, **defnrToPC = NULL;
	int *defnrToPCIndex = NULL;
	MDeformVert *dverts = NULL;
	bDeformGroup *dg;
	DualQuat *dualquats = NULL;
	float obinv[4][4], premat[4][4], postmat[4][4];
	const bool use_envelope = (deformflag & ARM_DEF_ENVELOPE) != 0;
	const bool use_quaternion = (deformflag & ARM_DEF_QUATERNION) != 0;
	const bool invert_vgroup = (deformflag & ARM_DEF_INVERT_VGROUP) != 0;
	int defbase_tot = 0;       /* safety for vertexgroup index overflow */
	int i, target_totvert = 0; /* safety for vertexgroup overflow */
	int armature_def_nr;
	int totchan;

	/* in editmode, or not an armature */
	if (arm->edbo || (armOb->pose == NULL)) {
		return;
	}

	if ((armOb->pose->flag & POSE_RECALC) != 0) {
		printf("ERROR! Trying to evaluate influence of armature '%s' which needs Pose recalc!", armOb->id.name);
		BLI_assert(0);
	}

	invert_m4_m4(obinv, target->obmat);
	copy_m4_m4(premat, target->obmat);
	mul_m4_m4m4(postmat, obinv, armOb->obmat);
	invert_m4_m4(premat, postmat);

	/* initialize B_bone matrices and dual quaternions */
	totchan = BLI_listbase_count(&armOb->pose->chanbase);

	if (use_quaternion) {
		dualquats = MEM_callocN(sizeof(DualQuat) * totchan, "dualquats");
	}

	pdef_info_array = MEM_callocN(sizeof(bPoseChanDeform) * totchan, "bPoseChanDeform");

	ArmatureBBoneDefmatsData data = {
		.pdef_info_array = pdef_info_array,.dualquats = dualquats,.use_quaternion = use_quaternion
	};
	BLI_task_parallel_listbase(&armOb->pose->chanbase, &data,
		gpencil_armature_bbone_defmats_cb, totchan > 512);

	/* get the def_nr for the overall armature vertex group if present */
	//creo que esto sobra
	armature_def_nr = defgroup_name_index(target, defgrp_name);

	defbase_tot = BLI_listbase_count(&target->defbase);
	dverts = gps->dvert;
	if (dverts)
		target_totvert = gps->totpoints;

	/* get a vertex-deform-index to posechannel array */
	defnrToPC = MEM_callocN(sizeof(*defnrToPC) * defbase_tot, "defnrToBone");
	defnrToPCIndex = MEM_callocN(sizeof(*defnrToPCIndex) * defbase_tot, "defnrToIndex");

	GHash *idx_hash = BLI_ghash_ptr_new("pose channel index by name");
	int pchan_index = 0;
	for (pchan = armOb->pose->chanbase.first; pchan != NULL; pchan = pchan->next, ++pchan_index) {
		BLI_ghash_insert(idx_hash, pchan, SET_INT_IN_POINTER(pchan_index));
	}
	for (i = 0, dg = target->defbase.first; dg; i++, dg = dg->next) {
		defnrToPC[i] = BKE_pose_channel_find_name(armOb->pose, dg->name);
		/* exclude non-deforming bones */
		if (defnrToPC[i]) {
			if (defnrToPC[i]->bone->flag & BONE_NO_DEFORM) {
				defnrToPC[i] = NULL;
			}
			else {
				defnrToPCIndex[i] = GET_INT_FROM_POINTER(BLI_ghash_lookup(idx_hash, defnrToPC[i]));
			}
		}
	}
	BLI_ghash_free(idx_hash, NULL, NULL);

	for (i = 0; i < numVerts; i++) {
		MDeformVert *dvert;
		DualQuat sumdq, *dq = NULL;
		float *co, dco[3];
		float sumvec[3], summat[3][3];
		float *vec = NULL, (*smat)[3] = NULL;
		float contrib = 0.0f;
		float armature_weight = 1.0f; /* default to 1 if no overall def group */
		float prevco_weight = 1.0f;   /* weight for optional cached vertexcos */

		if (use_quaternion) {
			memset(&sumdq, 0, sizeof(DualQuat));
			dq = &sumdq;
		}
		else {
			sumvec[0] = sumvec[1] = sumvec[2] = 0.0f;
			vec = sumvec;

			if (defMats) {
				zero_m3(summat);
				smat = summat;
			}
		}

		if (armature_def_nr != -1) {
			BLI_assert(i < gps->totpoints);
			dvert = gps->dvert + i;
		}

		if (armature_def_nr != -1 && dvert) {
			armature_weight = defvert_find_weight(dvert, armature_def_nr);

			if (invert_vgroup)
				armature_weight = 1.0f - armature_weight;

			/* hackish: the blending factor can be used for blending with prevCos too */
			if (prevCos) {
				prevco_weight = armature_weight;
				armature_weight = 1.0f;
			}
		}

		/* check if there's any  point in calculating for this vert */
		if (armature_weight == 0.0f)
			continue;

		/* get the coord we work on */
		co = prevCos ? prevCos[i] : vertexCos[i];

		/* Apply the object's matrix */
		mul_m4_v3(premat, co);

		if (dvert && dvert->totweight) { 
			MDeformWeight *dw = dvert->dw;
			int deformed = 0;
			unsigned int j;

			for (j = dvert->totweight; j != 0; j--, dw++) {
				const int index = dw->def_nr;
				if (index >= 0 && index < defbase_tot && (pchan = defnrToPC[index])) {
					float weight = dw->weight;
					Bone *bone = pchan->bone;
					pdef_info = pdef_info_array + defnrToPCIndex[index];

					deformed = 1;

					if (bone && bone->flag & BONE_MULT_VG_ENV) {
						weight *= distfactor_to_bone(co, bone->arm_head, bone->arm_tail,
							bone->rad_head, bone->rad_tail, bone->dist);
					}
					pchan_bone_deform(pchan, pdef_info, weight, vec, dq, smat, co, &contrib);
				}
			}
			/* if there are vertexgroups but not groups with bones
			 * (like for softbody groups) */
			if (deformed == 0 && use_envelope) {
				pdef_info = pdef_info_array;
				for (pchan = armOb->pose->chanbase.first; pchan; pchan = pchan->next, pdef_info++) {
					if (!(pchan->bone->flag & BONE_NO_DEFORM))
						contrib += dist_bone_deform(pchan, pdef_info, vec, dq, smat, co);
				}
			}
		}

		/* actually should be EPSILON? weight values and contrib can be like 10e-39 small */
		if (contrib > 0.0001f) {
			if (use_quaternion) {
				normalize_dq(dq, contrib);

				if (armature_weight != 1.0f) {
					copy_v3_v3(dco, co);
					mul_v3m3_dq(dco, (defMats) ? summat : NULL, dq);
					sub_v3_v3(dco, co);
					mul_v3_fl(dco, armature_weight);
					add_v3_v3(co, dco);
				}
				else
					mul_v3m3_dq(co, (defMats) ? summat : NULL, dq);

				smat = summat;
			}
			else {
				mul_v3_fl(vec, armature_weight / contrib);
				add_v3_v3v3(co, vec, co);
			}

			if (defMats) {
				float pre[3][3], post[3][3], tmpmat[3][3];

				copy_m3_m4(pre, premat);
				copy_m3_m4(post, postmat);
				copy_m3_m3(tmpmat, defMats[i]);

				if (!use_quaternion) /* quaternion already is scale corrected */
					mul_m3_fl(smat, armature_weight / contrib);

				mul_m3_series(defMats[i], post, smat, pre, tmpmat);
			}
		}

		/* always, check above code */
		mul_m4_v3(postmat, co);

		/* interpolate with previous modifier position using weight group */
		if (prevCos) {
			float mw = 1.0f - prevco_weight;
			vertexCos[i][0] = prevco_weight * vertexCos[i][0] + mw * co[0];
			vertexCos[i][1] = prevco_weight * vertexCos[i][1] + mw * co[1];
			vertexCos[i][2] = prevco_weight * vertexCos[i][2] + mw * co[2];
		}
	}

	if (dualquats)
		MEM_freeN(dualquats);
	if (defnrToPC)
		MEM_freeN(defnrToPC);
	if (defnrToPCIndex)
		MEM_freeN(defnrToPCIndex);

	/* free B_bone matrices */
	pdef_info = pdef_info_array;
	for (pchan = armOb->pose->chanbase.first; pchan; pchan = pchan->next, pdef_info++) {
		if (pdef_info->b_bone_mats)
			MEM_freeN(pdef_info->b_bone_mats);
		if (pdef_info->b_bone_dual_quats)
			MEM_freeN(pdef_info->b_bone_dual_quats);
	}

	MEM_freeN(pdef_info_array);
}

/* deform stroke */
static void deformStroke(
        GpencilModifierData *md, Depsgraph *UNUSED(depsgraph),
        Object *ob, bGPDlayer *gpl, bGPDstroke *gps)
{
	ArmatureGpencilModifierData *mmd = (ArmatureGpencilModifierData *)md;
	if (!mmd->object) {
		return;
	}

	//int vindex = defgroup_name_index(ob, mmd->vgname);
	//float weight = 1.0f;

	//bPoseChannel *pchan = BKE_pose_channel_find_name(mmd->object->pose, mmd->subtarget);
	//float dmat[4][4];

	/* loop points and apply deform */
	for (int i = 0; i < gps->totpoints; i++) {
		bGPDspoint *pt = &gps->points[i];
		MDeformVert *dvert = &gps->dvert[i];

		//weight = get_modifier_point_weight(dvert, 0, vindex);
		//if (weight < 0) {
		//	continue;
		//}
	//		gp_armature_co_apply(&tData, weight, pt);
	}
}

static void bakeModifier(
        Main *bmain, Depsgraph *depsgraph,
        GpencilModifierData *md, Object *ob)
{
	ArmatureGpencilModifierData *mmd = (ArmatureGpencilModifierData *)md;
	Scene *scene = DEG_get_evaluated_scene(depsgraph);
	bGPdata *gpd = ob->data;
	int oldframe = (int)DEG_get_ctime(depsgraph);

	if (mmd->object == NULL)
		return;

	for (bGPDlayer *gpl = gpd->layers.first; gpl; gpl = gpl->next) {
		for (bGPDframe *gpf = gpl->frames.first; gpf; gpf = gpf->next) {
			/* apply armature effects on this frame
			 * NOTE: this assumes that we don't want armature animation on non-keyframed frames
			 */
			CFRA = gpf->framenum;
			BKE_scene_graph_update_for_newframe(depsgraph, bmain);

			/* compute armature effects on this frame */
			for (bGPDstroke *gps = gpf->strokes.first; gps; gps = gps->next) {
				deformStroke(md, depsgraph, ob, gpl, gps);
			}
		}
	}

	/* return frame state and DB to original state */
	CFRA = oldframe;
	BKE_scene_graph_update_for_newframe(depsgraph, bmain);
}

static bool isDisabled(GpencilModifierData *md, int UNUSED(userRenderParams))
{
	ArmatureGpencilModifierData *mmd = (ArmatureGpencilModifierData *)md;

	return !mmd->object;
}

static void updateDepsgraph(GpencilModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
	ArmatureGpencilModifierData *lmd = (ArmatureGpencilModifierData *)md;
	if (lmd->object != NULL) {
		DEG_add_object_relation(ctx->node, lmd->object, DEG_OB_COMP_EVAL_POSE, "Armature Modifier");
		DEG_add_object_relation(ctx->node, lmd->object, DEG_OB_COMP_TRANSFORM, "Armature Modifier");
	}
	DEG_add_object_relation(ctx->node, ctx->object, DEG_OB_COMP_TRANSFORM, "Armature Modifier");
}

static void foreachObjectLink(
        GpencilModifierData *md, Object *ob,
        ObjectWalkFunc walk, void *userData)
{
	ArmatureGpencilModifierData *mmd = (ArmatureGpencilModifierData *)md;

	walk(userData, ob, &mmd->object, IDWALK_CB_NOP);
}

GpencilModifierTypeInfo modifierType_Gpencil_Armature = {
	/* name */              "Armature",
	/* structName */        "ArmatureGpencilModifierData",
	/* structSize */        sizeof(ArmatureGpencilModifierData),
	/* type */              eGpencilModifierTypeType_Gpencil,
	/* flags */             0,

	/* copyData */          copyData,

	/* deformStroke */      deformStroke,
	/* generateStrokes */   NULL,
	/* bakeModifier */    bakeModifier,

	/* initData */          initData,
	/* freeData */          NULL,
	/* isDisabled */        isDisabled,
	/* updateDepsgraph */   updateDepsgraph,
	/* dependsOnTime */     NULL,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
