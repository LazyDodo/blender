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

#include "GeometryExporter.h"
#include "AnimationExporter.h"
#include "AnimationClipExporter.h"
#include "BCAnimationCurveContainer.h"
#include "MaterialExporter.h"
#include "collada_utils.h"



std::map<std::string, std::vector<std::string>> BC_CHANNEL_NAME_FROM_TYPE = {
	{ "color"         , { "R", "G", "B" } },
	{ "specular_color", { "R", "G", "B" } },
	{ "diffuse_color",  { "R", "G", "B" } },
	{ "alpha",          { "R", "G", "B" } },
	{ "scale",          { "X", "Y", "Z" } },
	{ "location",       { "X", "Y", "Z" } },
	{ "rotation_euler", { "X", "Y", "Z" } }
};

std::string EMPTY_STRING;

std::string AnimationExporter::get_subchannel(std::string channel, int id)
{
	std::map<std::string, std::vector<std::string>>::const_iterator it;
	it = BC_CHANNEL_NAME_FROM_TYPE.find(channel);
	if (it == BC_CHANNEL_NAME_FROM_TYPE.end())
		return "";

	const std::vector<std::string> &subchannel = it->second;
	if (id >= subchannel.size())
		return "";
	return subchannel[id];
}

bool AnimationExporter::open_animation_container(bool has_container, Object *ob)
{
	if (!has_container) {
		char anim_id[200];
		sprintf(anim_id, "action_container-%s", translate_id(id_name(ob)).c_str());
		openAnimation(anim_id, id_name(ob));
	}
	return true;
}

void AnimationExporter::openAnimationWithClip(std::string action_id, std::string action_name)
{
	std::vector<std::string> anim_meta_entry;
	anim_meta_entry.push_back(translate_id(action_id));
	anim_meta_entry.push_back(action_name);
	anim_meta.push_back(anim_meta_entry);

	openAnimation(translate_id(action_id), action_name);
}

void AnimationExporter::close_animation_container(bool has_container)
{
	if (has_container)
		closeAnimation();
}

bool AnimationExporter::exportAnimations(Scene *sce)
{
	bool has_anim_data = BCAnimationSampler::has_animations(sce, this->export_settings->export_set);
	if (has_anim_data) {
		BCAnimationSampler sampler(mContext);

		this->scene = sce;

		std::set<Object *> animated_subset;
		BCAnimationSampler::get_animated_subset(animated_subset, this->export_settings->export_set);

		LinkNode *node;
		for (node = this->export_settings->export_set; node; node = node->next) {
			Object *ob = (Object *)node->link;
			if (animated_subset.find(ob) != animated_subset.end()) {
				//curve_container.addObject(ob);
				sampler.add_object(ob);
			}
		}

		sampler.sample_scene(scene,
			export_settings->sampling_rate,
			/*keyframe_at_end = */ true,
			export_settings->open_sim,
			export_settings->keep_keyframes,
			export_settings->export_animation_type
			);

		openLibrary();
		
		std::set<Object *>::iterator it;
		for (it = animated_subset.begin(); it != animated_subset.end(); ++it) {
			Object *ob = *it;
			exportObjectAnimation(ob, sampler);
		}

		closeLibrary();

#if 0
		/* TODO: If all actions shall be exported, we need to call the
		 * AnimationClipExporter which will figure out which actions
		 * need to be exported for which objects
		 */ 
		if (this->export_settings->include_all_actions) {
			AnimationClipExporter ace(eval_ctx, sw, export_settings, anim_meta);
			ace.exportAnimationClips(sce);
		}
#endif
	}
	return has_anim_data;
}

/* called for each exported object */
void AnimationExporter::exportObjectAnimation(Object *ob, BCAnimationSampler &sampler)
{
	bool container_is_open = false;

	//Transform animations (trans, rot, scale)
	container_is_open = open_animation_container(container_is_open, ob);

	/* Now take care of the Object Animations
	 * Note: For Armatures the skeletal animation has already been exported (see above)
	 * However Armatures also can have Object animation.
	 */
	if (this->export_settings->export_transformation_type == BC_TRANSFORMATION_TYPE_MATRIX) {
		export_matrix_animation_set(ob, sampler); // actually one matrix curve
	}
	else {
		export_curve_animation_set(ob, sampler); // each curve might have different frames
	}

	if (ob->type == OB_ARMATURE) {
		container_is_open = exportArmatureAnimation(ob, sampler, container_is_open);
	}

	close_animation_container(container_is_open);
}

bool AnimationExporter::exportArmatureAnimation(Object *ob, BCAnimationSampler &sampler, bool has_container)
{
	/* TODO: This needs to be handled by extra profiles, postponed for now
	* export_morph_animation(ob);
	*/

	if (ob->type == OB_ARMATURE) {
		/* Export skeletal animation (if any) */
		bArmature *arm = (bArmature *)ob->data;
		for (Bone *root_bone = (Bone *)arm->bonebase.first; root_bone; root_bone = root_bone->next)
			export_bone_animation_recursive(ob, root_bone, sampler);
	}

	return has_container;
}

/*
 * Export all animation FCurves of an Object.
 *
 * Note: This uses the keyframes as sample points,
 * and exports "baked keyframes" while keeping the tangent information
 * of the FCurves intact. This works for simple cases, but breaks
 * especially when negative scales are involved in the animation.
 * And when parent inverse matrices are involved (when exporting
 * object hierarchies)
 *
 */
void AnimationExporter::export_curve_animation_set(Object *ob, BCAnimationSampler &sampler)
{
	BCFrameSampleMap samples;
	BCAnimationCurveMap curves;

	sampler.get_curves(curves, ob);
	bool is_flat = sampler.get_samples(samples, ob);

	if (is_flat)
		return;

	BCAnimationCurveMap::iterator it;
	for (it = curves.begin(); it != curves.end(); ++it) {
		BCAnimationCurve &curve = it->second;
		if (curve.get_channel_target() == "rotation_quaternion") {
			/*
			   Can not export Quaternion animation in Collada as far as i know)
			   Maybe automatically convert to euler rotation?
			   Discard for now.
			*/
			continue;
		}

		sampler.add_value_set(curve, samples, this->export_settings->export_animation_type);  // prepare curve
		if (curve.is_flat())
			continue;

		export_curve_animation(ob, curve);
	}
}

void AnimationExporter::export_matrix_animation_set(Object *ob, BCAnimationSampler &sampler)
{
	std::vector<float> frames;
	sampler.get_frame_set(frames, ob);
	if (frames.size() > 0) {
		BCMatrixSampleMap samples;
		bool is_flat = sampler.get_samples(samples, ob);
		if (!is_flat) {
			export_matrix_animation(ob, frames, samples, sampler); // there is just one curve to export here
		}
	}
}

void AnimationExporter::export_matrix_animation(
	Object *ob, 
	BCFrames &frames, 
	BCMatrixSampleMap &samples,
	BCAnimationSampler &sampler)
{
	bAction *action = bc_getSceneObjectAction(ob);
	std::string name = id_name(ob);
	std::string action_name = (action == nullptr) ? name+"-action" : id_name(action);
	std::string channel_type = "transform";
	std::string axis = "";
	std::string id = bc_get_action_id(action_name, name, channel_type, axis);

	std::string target = translate_id(name) + '/' + channel_type;

	export_collada_matrix_animation( id, name, target, frames, samples);
}

//write bone animations in transform matrix sources
void AnimationExporter::export_bone_animation_recursive(Object *ob, Bone *bone, BCAnimationSampler &sampler)
{
	std::vector<float> frames;
	sampler.get_frame_set(frames, ob, bone);
	
	if (frames.size()) {
		BCMatrixSampleMap samples;
		bool is_flat = sampler.get_samples(samples, ob, bone);
		if (!is_flat) {
			export_bone_animation(ob, bone, frames, samples);
		}
	}

	for (Bone *child = (Bone *)bone->childbase.first; child; child = child->next)
		export_bone_animation_recursive(ob, child, sampler);
}

void AnimationExporter::export_morph_animation(Object *ob, BCAnimationSampler &sampler)
{
	FCurve *fcu;
	Key *key = BKE_key_from_object(ob);
	if (!key) return;

	if (key->adt && key->adt->action) {
		fcu = (FCurve *)key->adt->action->curves.first;

		while (fcu) {
			BC_animation_transform_type tm_type = get_transform_type(fcu->rna_path);
			create_keyframed_animation(ob, fcu, tm_type, true, sampler);
			fcu = fcu->next;
		}
	}

}

/* Euler sources from quternion sources 
 * Important: We assume the object has a scene action.
 * If it has not, then Blender will die
*/
void AnimationExporter::get_eul_source_for_quat(std::vector<float> &values, Object *ob)
{
	bAction *action = bc_getSceneObjectAction(ob);
	FCurve *fcu = (FCurve *)action->curves.first;
	const int keys = fcu->totvert;
	std::vector<std::vector<float>> quats;
	quats.resize(keys);
	for (int i = 0; i < keys; i++)
		quats[i].resize(4);

	int curve_count = 0;
	while (fcu) {
		std::string transformName = extract_transform_name(fcu->rna_path);

		if (transformName == "rotation_quaternion" ) {
			curve_count += 1;
			for (int i = 0; i < fcu->totvert; i++) {
				std::vector<float> &quat = quats[i];
				quat[fcu->array_index] = fcu->bezt[i].vec[1][1];
			}
			if (curve_count == 4)
				break; /* Quaternion curves can not use more the 4 FCurves!*/
		}
		fcu = fcu->next;
	}

	float feul[3];
	for (int i = 0; i < keys; i++) {
		std::vector<float> &quat = quats[i];
		quat_to_eul(feul, &quat[0]);

		for (int k = 0; k < 3; k++)
			values.push_back(feul[k]);
	}
}

void AnimationExporter::create_keyframed_animation(
	Object *ob, 
	FCurve *fcu, 
	BC_animation_transform_type tm_type,
	bool is_param, 
	BCAnimationSampler &sampler, 
	Material *ma)
{
	BCAnimationCurve curve(BC_ANIMATION_CURVE_TYPE_MATERIAL, fcu);
	export_curve_animation(ob, curve);
}

/* convert f-curves to animation curves and write
* Important: We assume the object has a scene action.
* If it has not, then Blender will die!
*/
void AnimationExporter::export_curve_animation(Object *ob, const BCAnimationCurve &curve)
{
	std::string channel = curve.get_channel_target();
	BC_animation_curve_type channel_type = curve.get_channel_type();
	int array_index = curve.get_array_index();
	std::string axis = get_subchannel(channel, array_index); // RGB or XYZ

	std::string action_name;
	bAction *action = bc_getSceneObjectAction(ob);
	action_name = (action) ? id_name(action) : "constraint_anim";

	const std::string curve_name = curve.get_animation_name(ob);
	std::string id = bc_get_action_id(action_name, curve_name, channel, axis, ".");
	
	std::string target = translate_id(curve_name);

	if (channel_type == BC_ANIMATION_CURVE_TYPE_LIGHT || channel_type == BC_ANIMATION_CURVE_TYPE_CAMERA) {
		target += "/" + get_param_sid(channel, axis);
	}
	else if (channel_type == BC_ANIMATION_CURVE_TYPE_MATERIAL) {
		int material_index = curve.get_tag();
		Material *ma = give_current_material(ob, material_index + 1);
		if (ma) {
			target = id_name(ma) + "-effect/common/" + curve.get_sid(axis);
		}
	}

	else {
		target += "/" + get_transform_sid(channel, axis);
	}

	export_collada_curve_animation(id, curve_name, target, axis, curve);
}

void AnimationExporter::export_bone_animation(Object *ob, Bone *bone, BCFrames &frames, BCMatrixSampleMap &samples)
{
	bAction* action = bc_getSceneObjectAction(ob);
	std::string bone_name(bone->name);
	std::string name = id_name(ob);
	std::string id = bc_get_action_id(id_name(action), name, bone_name, "pose_matrix");
	std::string target = translate_id(id_name(ob) + "_" + bone_name) + "/transform";

	export_collada_matrix_animation(id, name, target, frames, samples);
}

bool AnimationExporter::is_bone_deform_group(Bone *bone)
{
	bool is_def;
	//Check if current bone is deform
	if ((bone->flag & BONE_NO_DEFORM) == 0) return true;
	//Check child bones
	else {
		for (Bone *child = (Bone *)bone->childbase.first; child; child = child->next) {
			//loop through all the children until deform bone is found, and then return
			is_def = is_bone_deform_group(child);
			if (is_def) return true;
		}
	}
	//no deform bone found in children also
	return false;
}

void AnimationExporter::export_collada_curve_animation(
	std::string id,
	std::string name,
	std::string target,
	std::string axis,
	const BCAnimationCurve &curve)
{
	BCFrames frames;
	BCValues values;
	curve.get_sampled_frames(frames);
	curve.get_sampled_values(values);

	fprintf(stdout, "Export animation curve %s (%d control points)\n", id.c_str(), int(frames.size()));
	openAnimation(id, name);
	std::string intangent_id;
	std::string outtangent_id;
	bool has_tangents = false;
	bool is_rot = curve.is_rot();

	std::string input_id = create_source_from_values(COLLADASW::InputSemantic::INPUT, frames, false, id, axis);
	std::string output_id = create_source_from_values(COLLADASW::InputSemantic::OUTPUT, values, is_rot, id, axis);

	std::string interpolation_id;
	if (this->export_settings->keep_smooth_curves)
		interpolation_id = create_interpolation_source(curve, id, axis, &has_tangents);
	else
		interpolation_id = create_linear_interpolation_source(frames.size(), id);

	if (has_tangents) {
		intangent_id = create_tangent_from_curve(COLLADASW::InputSemantic::IN_TANGENT, curve, frames, id, axis);
		outtangent_id = create_tangent_from_curve(COLLADASW::InputSemantic::OUT_TANGENT, curve, frames, id, axis);
	}

	std::string sampler_id = std::string(id) + SAMPLER_ID_SUFFIX;
	COLLADASW::LibraryAnimations::Sampler sampler(sw, sampler_id);

	sampler.addInput(COLLADASW::InputSemantic::INPUT, COLLADABU::URI(EMPTY_STRING, input_id));
	sampler.addInput(COLLADASW::InputSemantic::OUTPUT, COLLADABU::URI(EMPTY_STRING, output_id));
	sampler.addInput(COLLADASW::InputSemantic::INTERPOLATION, COLLADABU::URI(EMPTY_STRING, interpolation_id));

	if (has_tangents) {
		sampler.addInput(COLLADASW::InputSemantic::IN_TANGENT, COLLADABU::URI(EMPTY_STRING, intangent_id));
		sampler.addInput(COLLADASW::InputSemantic::OUT_TANGENT, COLLADABU::URI(EMPTY_STRING, outtangent_id));
	}

	addSampler(sampler);
	addChannel(COLLADABU::URI(EMPTY_STRING, sampler_id), target);

	closeAnimation();
}

void AnimationExporter::export_collada_matrix_animation(std::string id, std::string name, std::string target, BCFrames &frames, BCMatrixSampleMap &samples)
{
	fprintf(stdout, "Export animation matrix %s (%d control points)\n", id.c_str(), int(frames.size()));

	openAnimationWithClip(id, name);

	std::string input_id = create_source_from_values(COLLADASW::InputSemantic::INPUT, frames, false, id, "");
	std::string output_id = create_4x4_source_from_values(samples, id);
	std::string interpolation_id = create_linear_interpolation_source(frames.size(), id);

	std::string sampler_id = std::string(id) + SAMPLER_ID_SUFFIX;
	COLLADASW::LibraryAnimations::Sampler sampler(sw, sampler_id);

	sampler.addInput(COLLADASW::InputSemantic::INPUT, COLLADABU::URI(EMPTY_STRING, input_id));
	sampler.addInput(COLLADASW::InputSemantic::OUTPUT, COLLADABU::URI(EMPTY_STRING, output_id));
	sampler.addInput(COLLADASW::InputSemantic::INTERPOLATION, COLLADABU::URI(EMPTY_STRING, interpolation_id));

	// Matrix animation has no tangents

	addSampler(sampler);
	addChannel(COLLADABU::URI(EMPTY_STRING, sampler_id), target);

	closeAnimation();
}


float AnimationExporter::convert_angle(float angle)
{
	return COLLADABU::Math::Utils::radToDegF(angle);
}

std::string AnimationExporter::get_semantic_suffix(COLLADASW::InputSemantic::Semantics semantic)
{
	switch (semantic) {
		case COLLADASW::InputSemantic::INPUT:
			return INPUT_SOURCE_ID_SUFFIX;
		case COLLADASW::InputSemantic::OUTPUT:
			return OUTPUT_SOURCE_ID_SUFFIX;
		case COLLADASW::InputSemantic::INTERPOLATION:
			return INTERPOLATION_SOURCE_ID_SUFFIX;
		case COLLADASW::InputSemantic::IN_TANGENT:
			return INTANGENT_SOURCE_ID_SUFFIX;
		case COLLADASW::InputSemantic::OUT_TANGENT:
			return OUTTANGENT_SOURCE_ID_SUFFIX;
		default:
			break;
	}
	return "";
}

void AnimationExporter::add_source_parameters(COLLADASW::SourceBase::ParameterNameList& param,
	COLLADASW::InputSemantic::Semantics semantic,
	bool is_rot, 
	const std::string axis, 
	bool transform)
{
	switch (semantic) {
		case COLLADASW::InputSemantic::INPUT:
			param.push_back("TIME");
			break;
		case COLLADASW::InputSemantic::OUTPUT:
			if (is_rot) {
				param.push_back("ANGLE");
			}
			else {
				if (axis != "") {
					param.push_back(axis);
				}
				else
				if (transform) {
					param.push_back("TRANSFORM");
				}
				else {     //assumes if axis isn't specified all axises are added
					param.push_back("X");
					param.push_back("Y");
					param.push_back("Z");
				}
			}
			break;
		case COLLADASW::InputSemantic::IN_TANGENT:
		case COLLADASW::InputSemantic::OUT_TANGENT:
			param.push_back("X");
			param.push_back("Y");
			break;
		default:
			break;
	}
}

/*
Experimental stuff: try to find a tangent that lies between the closest fcurve points.
if the sample_point hits a function point, then the tangents of the function point are preserved.
Cool eh?
*/
int AnimationExporter::get_source_values(const BCAnimationCurve &curve, float val, COLLADASW::InputSemantic::Semantics semantic, bool is_angle, float *values)
{
	const FCurve *fcu = curve.get_fcurve();

	int lower_index = curve.closest_index_below(val);
	int upper_index = curve.closest_index_above(val, lower_index);
	if (lower_index==upper_index)
		return get_source_values(&fcu->bezt[lower_index], semantic, is_angle, values);

	/* This actually does not work */
	float lower[2];
	float upper[2];

	lower[0] = fcu->bezt[lower_index].vec[1][0]; // ctime
	lower[1] = fcu->bezt[lower_index].vec[1][1]; // value
	upper[0] = fcu->bezt[upper_index].vec[1][0]; // ctime
	upper[1] = fcu->bezt[upper_index].vec[1][1]; // value

	float fraction = (val-lower[1]) / (upper[1]-lower[1]);

	float lower_values[3];
	float upper_values[3];

	int length = get_source_values(&fcu->bezt[lower_index], semantic, is_angle, lower_values);
	             get_source_values(&fcu->bezt[upper_index], semantic, is_angle, upper_values);

	for (int i = 0; i < 2; i++) {
		values[i] = lower_values[i] + (upper_values[i] - lower_values[i]) * fraction;
	}

	return length;
}

int AnimationExporter::get_source_values(BezTriple *bezt, COLLADASW::InputSemantic::Semantics semantic, bool is_angle, float *values)
{
	int length;
	switch (semantic) {
		case COLLADASW::InputSemantic::INPUT:
			length = 1;
			values[0] = FRA2TIME(bezt->vec[1][0]);
			break;
		case COLLADASW::InputSemantic::OUTPUT:
			length = 1;
			if (is_angle) {
				values[0] = RAD2DEGF(bezt->vec[1][1]);
			}
			else {
				values[0] = bezt->vec[1][1];
			}
			break;

		case COLLADASW::InputSemantic::IN_TANGENT:
			length = 2;
			values[0] = FRA2TIME(bezt->vec[0][0]);
			if (bezt->ipo != BEZT_IPO_BEZ) {
				// We're in a mixed interpolation scenario, set zero as it's irrelevant but value might contain unused data
				values[0] = 0;
				values[1] = 0;
			}
			else if (is_angle) {
				values[1] = RAD2DEGF(bezt->vec[0][1]);
			}
			else {
				values[1] = bezt->vec[0][1];
			}
			break;

		case COLLADASW::InputSemantic::OUT_TANGENT:
			length = 2;
			values[0] = FRA2TIME(bezt->vec[2][0]);
			if (bezt->ipo != BEZT_IPO_BEZ) {
				// We're in a mixed interpolation scenario, set zero as it's irrelevant but value might contain unused data
				values[0] = 0;
				values[1] = 0;
			}
			else if (is_angle) {
				values[1] = RAD2DEGF(bezt->vec[2][1]);
			}
			else {
				values[1] = bezt->vec[2][1];
			}
			break;
		default:
			length = 0;
			break;
	}
	return length;
}

std::string AnimationExporter::create_tangent_from_curve(COLLADASW::InputSemantic::Semantics semantic, const BCAnimationCurve &curve, std::vector<float>frames, const std::string& anim_id, std::string axis_name)
{
	std::string channel = curve.get_channel_target();

	const std::string source_id = anim_id + get_semantic_suffix(semantic);

	bool is_angle = (bc_startswith(channel, "rotation") || channel == "spot_size");
	bool is_euler = (channel == "rotation_euler");

	COLLADASW::FloatSourceF source(mSW);
	source.setId(source_id);
	source.setArrayId(source_id + ARRAY_ID_SUFFIX);
	source.setAccessorCount(curve.size());
	source.setAccessorStride(2);

	COLLADASW::SourceBase::ParameterNameList &param = source.getParameterNameList();
	add_source_parameters(param, semantic, is_angle, axis_name, false);

	source.prepareToAppendValues();

	std::vector<float> values;
	curve.get_sampled_values(values);

	const FCurve *fcu = curve.get_fcurve(); // need this to get the original tangents

	for (unsigned int frame_index = 0; frame_index < values.size(); frame_index++) {
		float sampled_val = values[frame_index];

		if (is_angle) {
			sampled_val = RAD2DEGF(sampled_val);
		}

		float vals[3]; // be careful!
		int length = get_source_values(curve, frames[frame_index], semantic, is_angle, vals);
		float offset = 0;
		float bases[3];
		int len = get_source_values(curve, frames[frame_index], COLLADASW::InputSemantic::OUTPUT, is_angle, bases);
		sampled_val += vals[1] - bases[0];

		source.appendValues(vals[0]);
		source.appendValues(sampled_val);

	}
	source.finish();
	return source_id;
}

/*
 * Similar to create_source_from_fcurve, but adds conversion of lens
 * animation data from focal length to FOV.
 */
std::string AnimationExporter::create_lens_source_from_fcurve(Camera *cam, COLLADASW::InputSemantic::Semantics semantic, const BCAnimationCurve &curve, const std::string& anim_id)
{
	std::string source_id = anim_id + get_semantic_suffix(semantic);

	COLLADASW::FloatSourceF source(mSW);
	source.setId(source_id);
	source.setArrayId(source_id + ARRAY_ID_SUFFIX);
	source.setAccessorCount(curve.size());

	source.setAccessorStride(1);

	COLLADASW::SourceBase::ParameterNameList &param = source.getParameterNameList();
	add_source_parameters(param, semantic, false, "", false);

	source.prepareToAppendValues();

	const FCurve *fcu = curve.get_fcurve();
	for (unsigned int i = 0; i < curve.size(); i++) {
		float values[3]; // be careful!
		int length = get_source_values(&fcu->bezt[i], semantic, false, values);
		for (int j = 0; j < length; j++)
		{
			float val = RAD2DEGF(focallength_to_fov(values[j], cam->sensor_x));
			source.appendValues(val);
		}
	}

	source.finish();

	return source_id;
}

std::string AnimationExporter::create_source_from_values(COLLADASW::InputSemantic::Semantics semantic, std::vector<float> &values, bool is_rot, const std::string& anim_id, const std::string axis_name)
{
	/* T can be float, int or double */

	int stride = 1;
	int entry_count = values.size() / stride;
	std::string source_id = anim_id + get_semantic_suffix(semantic);

	COLLADASW::FloatSourceF source(mSW);
	source.setId(source_id);
	source.setArrayId(source_id + ARRAY_ID_SUFFIX);
	source.setAccessorCount(entry_count);
	source.setAccessorStride(stride);

	COLLADASW::SourceBase::ParameterNameList &param = source.getParameterNameList();
	add_source_parameters(param, semantic, is_rot, axis_name, false);

	source.prepareToAppendValues();

	for (int i = 0; i < entry_count; i++) {
		float val = values[i];
		if (semantic == COLLADASW::InputSemantic::INPUT)
			val = FRA2TIME(val);
		if (is_rot)
			val = RAD2DEGF(val);
		source.appendValues(val);
	}

	source.finish();

	return source_id;
}

/*
 * Create a collada matrix source for a set of samples
*/
std::string AnimationExporter::create_4x4_source_from_values(BCMatrixSampleMap &samples, const std::string &anim_id)
{
	COLLADASW::InputSemantic::Semantics semantic = COLLADASW::InputSemantic::OUTPUT;
	std::string source_id = anim_id + get_semantic_suffix(semantic);

	COLLADASW::Float4x4Source source(mSW);
	source.setId(source_id);
	source.setArrayId(source_id + ARRAY_ID_SUFFIX);
	source.setAccessorCount(samples.size());
	source.setAccessorStride(16);

	COLLADASW::SourceBase::ParameterNameList &param = source.getParameterNameList();
	add_source_parameters(param, semantic, false, "", true);

	source.prepareToAppendValues();

	BCMatrixSampleMap::iterator it;
	int j = 0;
	int precision = (this->export_settings->limit_precision) ? 6 : -1; // could be made configurable
	for (it = samples.begin(); it != samples.end(); it++) {
		const BCMatrix *sample = it->second;
		double daemat[4][4];
		sample->get_matrix(daemat, true, precision);
		source.appendValues(daemat);
	}

	source.finish();
	return source_id;
}

std::string AnimationExporter::create_interpolation_source(const BCAnimationCurve &curve, 
	const std::string& anim_id, 
	const std::string axis,
	bool *has_tangents)
{
	std::string source_id = anim_id + get_semantic_suffix(COLLADASW::InputSemantic::INTERPOLATION);

	COLLADASW::NameSource source(mSW);
	source.setId(source_id);
	source.setArrayId(source_id + ARRAY_ID_SUFFIX);
	source.setAccessorCount(curve.size());
	source.setAccessorStride(1);

	COLLADASW::SourceBase::ParameterNameList &param = source.getParameterNameList();
	param.push_back("INTERPOLATION");

	source.prepareToAppendValues();

	*has_tangents = false;

	const FCurve *fcu = curve.get_fcurve();
	std::vector<float>frames;
	curve.get_sampled_frames(frames);

	for (unsigned int i = 0; i < curve.size(); i++) {
		float frame = frames[i];
		int ipo = curve.get_ipo(frame);
		if (ipo == BEZT_IPO_BEZ) {
			source.appendValues(BEZIER_NAME);
			*has_tangents = true;
		}
		else if (ipo == BEZT_IPO_CONST) {
			source.appendValues(STEP_NAME);
		}
		else { // BEZT_IPO_LIN
			source.appendValues(LINEAR_NAME);
		}
	}
	// unsupported? -- HERMITE, CARDINAL, BSPLINE, NURBS

	source.finish();

	return source_id;
}

std::string AnimationExporter::create_linear_interpolation_source(int tot, const std::string& anim_id)
{
	std::string source_id = anim_id + get_semantic_suffix(COLLADASW::InputSemantic::INTERPOLATION);

	COLLADASW::NameSource source(mSW);
	source.setId(source_id);
	source.setArrayId(source_id + ARRAY_ID_SUFFIX);
	source.setAccessorCount(tot);
	source.setAccessorStride(1);

	COLLADASW::SourceBase::ParameterNameList &param = source.getParameterNameList();
	param.push_back("INTERPOLATION");

	source.prepareToAppendValues();

	for (int i = 0; i < tot; i++) {
		source.appendValues(LINEAR_NAME);
	}

	source.finish();

	return source_id;
}

extern std::map<BC_animation_transform_type, std::string> BC_ANIMATION_NAME_FROM_TYPE;
extern std::map<std::string, BC_animation_transform_type> BC_ANIMATION_TYPE_FROM_NAME;

BC_animation_transform_type AnimationExporter::get_transform_type(std::string path)
{
	BC_animation_transform_type tm_type;
	// when given rna_path, overwrite tm_type from it
	std::string name = extract_transform_name(path);
	std::map<std::string, BC_animation_transform_type>::iterator type_it = BC_ANIMATION_TYPE_FROM_NAME.find(name);
	tm_type = (type_it != BC_ANIMATION_TYPE_FROM_NAME.end()) ? type_it->second : BC_ANIMATION_TYPE_UNKNOWN;

	return tm_type;
}

std::string AnimationExporter::get_param_sid(BC_animation_transform_type tm_type, const std::string axis_name)
{
	std::string tm_name;
	std::map<BC_animation_transform_type, std::string>::iterator name_it = BC_ANIMATION_NAME_FROM_TYPE.find(tm_type);
	tm_name = name_it->second;

	return tm_name;
}

std::string AnimationExporter::get_param_sid(std::string path, const std::string axis_name)
{
	BC_animation_transform_type tm_type = get_transform_type(path);
	return get_param_sid(tm_type, axis_name);
}

/*
* Assign sid of the animated parameter or transform for rotation,
* axis name is always appended and the value of append_axis is ignored
*/
std::string AnimationExporter::get_transform_sid(const std::string channel, const std::string axis_name)
{
	std::string tm_name;
	BC_animation_transform_type tm_type = get_transform_type(channel);

	if (
		tm_type == BC_ANIMATION_TYPE_ROTATION_EULER ||
		tm_type == BC_ANIMATION_TYPE_ROTATION_QUAT
		)
	{
		tm_type = BC_ANIMATION_TYPE_ROTATION;
	}

	bool is_angle = (
		tm_type == BC_ANIMATION_TYPE_ROTATION ||
		tm_type == BC_ANIMATION_TYPE_ROTATION_EULER ||
		tm_type == BC_ANIMATION_TYPE_ROTATION_QUAT);

	std::map<BC_animation_transform_type, std::string>::iterator name_it = BC_ANIMATION_NAME_FROM_TYPE.find(tm_type);
	tm_name = name_it->second;

	if (tm_name.size()) {
		if (is_angle)
			return tm_name + std::string(axis_name) + ".ANGLE";
		else
			if (axis_name[0])
				return tm_name + "." + std::string(axis_name);
			else
				return tm_name;
	}

	return tm_name;
}

/*
 * enable fcurves driving a specific bone, disable all the rest
 * if bone_name = NULL enable all fcurves
 */
void AnimationExporter::enable_fcurves(bAction *act, char *bone_name)
{
	FCurve *fcu;
	char prefix[200];

	if (bone_name)
		BLI_snprintf(prefix, sizeof(prefix), "pose.bones[\"%s\"]", bone_name);

	for (fcu = (FCurve *)act->curves.first; fcu; fcu = fcu->next) {
		if (bone_name) {
			if (STREQLEN(fcu->rna_path, prefix, strlen(prefix)))
				fcu->flag &= ~FCURVE_DISABLED;
			else
				fcu->flag |= FCURVE_DISABLED;
		}
		else {
			fcu->flag &= ~FCURVE_DISABLED;
		}
	}
}
