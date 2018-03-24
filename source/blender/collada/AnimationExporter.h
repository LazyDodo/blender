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
 * Contributor(s): Chingiz Dyussenov, Arystanbek Dyussenov, Jan Diederich, Tod Liverseed.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __BC_ANIMATION_EXPORTER_H__
#define __BC_ANIMATION_EXPORTER_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
extern "C"
{
#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_anim_types.h"
#include "DNA_action_types.h"
#include "DNA_curve_types.h"
#include "DNA_lamp_types.h"
#include "DNA_camera_types.h"
#include "DNA_armature_types.h"
#include "DNA_material_types.h"
#include "DNA_constraint_types.h"
#include "DNA_scene_types.h"

#include "BLI_math.h"
#include "BLI_string.h"
#include "BLI_listbase.h"
#include "BLI_utildefines.h"

#include "BKE_DerivedMesh.h"
#include "BKE_fcurve.h"
#include "BKE_animsys.h"
#include "BKE_scene.h"
#include "BKE_action.h" // pose functions
#include "BKE_armature.h"
#include "BKE_object.h"
#include "BKE_constraint.h"
#include "BIK_api.h"
#include "ED_object.h"
}

#include "MEM_guardedalloc.h"

#include "RNA_access.h"

#include "COLLADASWSource.h"
#include "COLLADASWInstanceGeometry.h"
#include "COLLADASWInputList.h"
#include "COLLADASWPrimitves.h"
#include "COLLADASWVertices.h"
#include "COLLADASWLibraryAnimations.h"
#include "COLLADASWParamTemplate.h"
#include "COLLADASWParamBase.h"
#include "COLLADASWSampler.h"
#include "COLLADASWConstants.h"
#include "COLLADASWBaseInputElement.h"

#include "EffectExporter.h"
#include "BCAnimationCurveContainer.h"
#include "collada_internal.h"

#include "IK_solver.h"

#include <vector>
#include <map>
#include <algorithm> // std::find

class AnimationExporter: COLLADASW::LibraryAnimations
{
private:
	bContext *mContext;
	const ExportSettings *export_settings;
	EvaluationContext * eval_ctx;
	Scene *scene;
	COLLADASW::StreamWriter *sw;
	//AnimationCurveCache cache;

	std::vector<std::vector<std::string>> anim_meta;

public:

	AnimationExporter(bContext *C, EvaluationContext *eval_ctx, COLLADASW::StreamWriter *sw, const ExportSettings *export_settings):
		mContext(C),
		eval_ctx(eval_ctx),
		COLLADASW::LibraryAnimations(sw),
		export_settings(export_settings)
	{
		this->sw = sw;
	}

	bool exportAnimations(Scene *sce);

	// Main entry point into Animation export (called for each exported object)
	void exportObjectAnimation(Object *ob, BCAnimationSampler &sampler);

protected:

	void export_curve_animation_set(
		Object *ob,
		BCAnimationSampler &sampler);

	void export_curve_animation(
		Object *ob,
		const BCAnimationCurve &curve);

	void export_collada_curve_animation(
		std::string id,
		std::string name,
		std::string target,
		std::string axis,
		const BCAnimationCurve &curve);

	void export_matrix_animation_set(
		Object *ob,
		BCAnimationSampler &sampler);

	void export_matrix_animation (
		Object *ob,
		BCFrames &frames,
		BCFrameSampleMap &outmats,
		BCAnimationSampler &sampler);

	void export_bone_animation_recursive(
		Object *ob_arm, 
		Bone *bone, 
		BCAnimationSampler &sampler);

	void export_bone_animation(
		Object *ob,
		Bone *bone,
		BCFrames &frames,
		BCFrameSampleMap &outmats);

	void export_collada_matrix_animation(
		std::string id,
		std::string name,
		std::string target,
		BCFrames &frames,
		BCFrameSampleMap &outmats);

	/* Helper functions */
	void openAnimationWithClip(std::string id, std::string name);
	bool open_animation_container(bool has_container, Object *ob);
	void close_animation_container(bool has_container);

	std::string create_source_from_values(COLLADASW::InputSemantic::Semantics semantic, std::vector<float> &values, bool is_rot, const std::string& anim_id, const std::string axis_name);
	std::string create_4x4_source_from_values(
		BCFrameSampleMap &cache,
		const std::string& anim_id);

	std::string create_linear_interpolation_source(int tot, const std::string& anim_id);

	std::string get_semantic_suffix(COLLADASW::InputSemantic::Semantics semantic);

	void add_source_parameters(COLLADASW::SourceBase::ParameterNameList& param,
		COLLADASW::InputSemantic::Semantics semantic,
		bool is_rot,
		const std::string axis,
		bool transform);

	float convert_angle(float angle);

	// Now we get into the mess ...

	void export_morph_animation(
		Object *ob, 
		BCAnimationSampler &sampler);

	bool is_bone_deform_group(Bone * bone);

	int get_source_values(BezTriple *bezt, COLLADASW::InputSemantic::Semantics semantic, bool is_angle, float *values);
	int get_source_values(const BCAnimationCurve &curve, float sample_frame, COLLADASW::InputSemantic::Semantics semantic, bool is_angle, float *values);

	void get_eul_source_for_quat(std::vector<float> &cache, Object *ob );

	void create_keyframed_animation(
		Object *ob, 
		FCurve *fcu, 
		BC_animation_transform_type mt_type,
		bool is_param, 
		BCAnimationSampler &sampler, 
		Material *ma = NULL);

	std::string create_tangent_from_curve(COLLADASW::InputSemantic::Semantics semantic, const BCAnimationCurve &curve, std::vector<float>frames, const std::string& anim_id, const std::string axis_name);
	std::string create_lens_source_from_fcurve(Camera *cam, COLLADASW::InputSemantic::Semantics semantic, const BCAnimationCurve &curve, const std::string& anim_id);
	std::string create_interpolation_source(const BCAnimationCurve &curve, const std::string& anim_id, std::string axis_name, bool *has_tangents);
	
	// for rotation, axis name is always appended and the value of append_axis is ignored
	std::string get_transform_sid(const std::string path, const std::string axis_name);
	std::string get_param_sid(BC_animation_transform_type mt_type, const std::string axis_name);
	std::string get_param_sid(const std::string path, const std::string axis_name);
	BC_animation_transform_type get_transform_type(const std::string path);

	// enable fcurves driving a specific bone, disable all the rest
	// if bone_name = NULL enable all fcurves
	void enable_fcurves(bAction *act, char *bone_name);
	
	inline const std::string extract_transform_name(const std::string path)
	{
		return bc_string_after(path, '.');
	}

	std::string get_subchannel(std::string channel, int id);

};

#endif
