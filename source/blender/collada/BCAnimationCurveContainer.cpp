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

#include <vector>
#include <map>
#include <algorithm> // std::find

#include "ExportSettings.h"
#include "BCAnimationCurve.h"
#include "BCAnimationCurveContainer.h"

extern "C" {
#include "BKE_action.h"
#include "BKE_constraint.h"
#include "BKE_material.h"
#include "BLI_listbase.h"
#include "DNA_scene_types.h"
#include "DNA_constraint_types.h"
#include "ED_object.h"
}

static std::string EMPTY_STRING;
static BCAnimationCurveMap BCEmptyAnimationCurves;
static const BCMatrix UNIT_MATRIX;

/* 
    BCAnimationSampler
*/

BCAnimationSampler::BCAnimationSampler(bContext *C):
	mContext(C)
{
	// nothing todo here
}

void BCAnimationSampler::addObject(Object *ob)
{
	BCFrameSet &keyframes = objects[ob];
	get_keyframes(ob, keyframes);
}

void BCAnimationSampler::enable_fcurves(bAction *act, char *bone_name)
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

bool BCAnimationSampler::is_animated_by_constraint(Object *ob, ListBase *conlist, std::set<Object *> &animated_objects)
{
	bConstraint *con;
	for (con = (bConstraint *)conlist->first; con; con = con->next) {
		ListBase targets = { NULL, NULL };

		const bConstraintTypeInfo *cti = BKE_constraint_typeinfo_get(con);

		if (!bc_validateConstraints(con))
			continue;

		if (cti && cti->get_constraint_targets) {
			bConstraintTarget *ct;
			Object *obtar;
			cti->get_constraint_targets(con, &targets);
			for (ct = (bConstraintTarget *)targets.first; ct; ct = ct->next) {
				obtar = ct->tar;
				if (obtar) {
					if (animated_objects.find(obtar) != animated_objects.end())
						return true;
				}
			}
		}
	}
	return false;
}

void BCAnimationSampler::find_depending_animated(std::set<Object *> &animated_objects, std::set<Object *> &candidates)
{
	bool found_more;
	do {
		found_more = false;
		std::set<Object *>::iterator it;
		for (it = candidates.begin(); it != candidates.end(); ++it) {
			Object *cob = *it;
			ListBase *conlist = get_active_constraints(cob);
			if (is_animated_by_constraint(cob, conlist, animated_objects)) {
				animated_objects.insert(cob);
				candidates.erase(cob);
				found_more = true;
				break;
			}
		}
	} while (found_more && candidates.size() > 0);
}

void BCAnimationSampler::get_animated_subset(std::set<Object *> &animated_objects, LinkNode *export_set)
{
	/*
	Check if this object is animated. That is: Check if it has its own action, or

	- Check if it has constraints to other objects
	- at least one of the other objects is animated as well
	*/

	animated_objects.clear();
	std::set<Object *> static_objects;
	std::set<Object *> candidates;

	LinkNode *node;
	for (node = export_set; node; node = node->next) {
		Object *cob = (Object *)node->link;
		if (bc_getSceneObjectAction(cob))
			animated_objects.insert(cob);
		else {
			ListBase conlist = cob->constraints;
			if (conlist.first)
				candidates.insert(cob);
		}
	}
	find_depending_animated(animated_objects, candidates);
}
bool BCAnimationSampler::bone_matrix_local_get(Object *ob, Bone *bone, float(&mat)[4][4], bool for_opensim)
{

	/* Ok, lets be super cautious and check if the bone exists */
	bPose *pose = ob->pose;
	bPoseChannel *pchan = BKE_pose_channel_find_name(pose, bone->name);
	if (!pchan) {
		return false;
	}

	bAction *action = bc_getSceneObjectAction(ob);
	bPoseChannel *parchan = pchan->parent;

	enable_fcurves(action, bone->name);
	float ipar[4][4];

	if (bone->parent) {
		invert_m4_m4(ipar, parchan->pose_mat);
		mul_m4_m4m4(mat, ipar, pchan->pose_mat);
	}
	else
		copy_m4_m4(mat, pchan->pose_mat);

	/* OPEN_SIM_COMPATIBILITY
	* AFAIK animation to second life is via BVH, but no
	* reason to not have the collada-animation be correct
	*/
	if (for_opensim) {
		float temp[4][4];
		copy_m4_m4(temp, bone->arm_mat);
		temp[3][0] = temp[3][1] = temp[3][2] = 0.0f;
		invert_m4(temp);

		mul_m4_m4m4(mat, mat, temp);

		if (bone->parent) {
			copy_m4_m4(temp, bone->parent->arm_mat);
			temp[3][0] = temp[3][1] = temp[3][2] = 0.0f;

			mul_m4_m4m4(mat, temp, mat);
		}
	}
	enable_fcurves(action, NULL);
	return true;
}

static void get_sample_frames(BCFrameSet &sample_frames, int sampling_rate, bool keyframe_at_end, Scene *scene)
{
	sample_frames.clear();

	if (sampling_rate < 1)
		return; // no sample frames in this case

	float sfra = scene->r.sfra;
	float efra = scene->r.efra;

	int frame_index;
	for (frame_index = int(sfra); frame_index < efra; frame_index += sampling_rate) {
		sample_frames.insert(frame_index);
	}
	
	if (frame_index >= efra && keyframe_at_end)
	{
		sample_frames.insert(efra);
	}
}

static bool is_object_keyframe(Object *ob, int frame_index)
{
	return false;
}

static void add_keyframes_from(bAction *action, BCFrameSet &frameset)
{
	if (action) {
		FCurve *fcu = nullptr;
		for (fcu = (FCurve *)action->curves.first; fcu; fcu = fcu->next) {
			BezTriple  *bezt = fcu->bezt;
			for (int i = 0; i < fcu->totvert; bezt++, i++) {
				int frame_index = bezt->vec[1][0];
				frameset.insert(frame_index);
			}
		}
	}
}

/*
 * Collect all keyframes from all animation curves related to the object 
 * The bc_get... functions check for nullptr and correct object type
 * The add_keyframes_from() function checks for nullptr 
*/
void BCAnimationSampler::get_keyframes(Object *ob, BCFrameSet &frameset)
{
	frameset.clear();
	add_keyframes_from(bc_getSceneObjectAction(ob), frameset);
	add_keyframes_from(bc_getSceneCameraAction(ob), frameset);
	add_keyframes_from(bc_getSceneLampAction(ob), frameset);

	for (int a = 0; a < ob->totcol; a++) {
		Material *ma = give_current_material(ob, a + 1);
		add_keyframes_from(bc_getSceneMaterialAction(ma), frameset);
	}
}

void BCAnimationSampler::sample_scene(
	Scene *scene,
	int sampling_rate,
	int keyframe_at_end,
	bool for_opensim,
	bool keep_keyframes,
	BC_export_animation_type export_animation_type)
{
	BCFrameSet scene_sample_frames;
	get_sample_frames(scene_sample_frames, sampling_rate, keyframe_at_end, scene);
	BCFrameSet::iterator it;

	int startframe = scene->r.sfra;
	int endframe = scene->r.efra;

	for (int frame_index = startframe; frame_index <= endframe; ++frame_index) {
		/* Loop over all frames and decide for each frame if sampling is necessary */
		bool is_scene_sample_frame = false;
		bool needs_update = true;
		if (scene_sample_frames.find(frame_index) != scene_sample_frames.end()) {
			bc_update_scene(mContext, scene, frame_index);
			needs_update = false;
			is_scene_sample_frame = true;
		}

		bool needs_sampling = is_scene_sample_frame || keep_keyframes || export_animation_type == BC_ANIMATION_TYPE_KEYS;
		if (!needs_sampling) {
			continue;
		}

		BCAnimatedObjectMap::iterator obit;
		for (obit = objects.begin(); obit != objects.end(); ++obit) {
			Object *ob = obit->first;
			BCFrameSet &object_keyframes = obit->second;
			if (is_scene_sample_frame || object_keyframes.find(frame_index) != object_keyframes.end() ) {

				if (needs_update) {
					bc_update_scene(mContext, scene, frame_index);
					needs_update = false;
				}

				float mat[4][4];
				BKE_object_matrix_local_get(ob, mat);
				BCMatrix ob_mat(mat);
				sample_data.add(ob, ob_mat, frame_index);

				if (ob->type == OB_ARMATURE) {
					bPoseChannel *pchan;
					for (pchan = (bPoseChannel *)ob->pose->chanbase.first; pchan; pchan = pchan->next) {
						Bone *bone = pchan->bone;
						if (bone_matrix_local_get(ob, bone, mat, for_opensim)) {
							BCMatrix bone_mat(mat);
							sample_data.add(ob, bone, bone_mat, frame_index);
						}
					}
				}

				// TODO ...
				bAction *action;
				if (ob->type == OB_CAMERA) {
					action = bc_getSceneCameraAction(ob);
				}
				else if (ob->type == OB_LAMP) {
					action = bc_getSceneLampAction(ob);
				}
				if (action) {
					// TODO add curves to sample list...
				}

				for (int a = 0; a < ob->totcol; a++) {
					Material *ma = give_current_material(ob, a + 1);
					if (ma) {
						action = bc_getSceneMaterialAction(ma);
						if (action) {
							// TODO add curves to sample list...
						}
					}	    
				}
			}
		}
	}
}


bool BCAnimationSampler::is_flat_line(BCMatrixMap &values) const
{
	static float MIN_DISTANCE = 0.00001;

	if (values.size() < 2)
		return true; // need at least 2 entries to be not flat

	BCMatrixMap::iterator it;
	const BCMatrix *refmat = nullptr;
	for (it = values.begin(); it != values.end(); ++it) {
		const BCMatrix *matrix = it->second;

		if (refmat == nullptr) {
			refmat = matrix;
			continue;
		}

		if (!matrix->in_range(*refmat, MIN_DISTANCE))
			return false;
	}
	return true;
}

bool BCAnimationSampler::is_flat_line(std::vector<float> &values) const
{
	return BCAnimationCurve::is_flat_line(values);
}

void BCAnimationSampler::get_frame_set(BCFrames &frames, Object *ob)
{
	BCSampleKey key(ob);
	sample_data.get_frames(key, frames);
}

void BCAnimationSampler::get_frame_set(BCFrames &frames, Object *ob, Bone *bone)
{
	BCSampleKey key(ob, bone);
	sample_data.get_frames(key, frames);
}

void BCAnimationSampler::get_frame_set(BCFrames &frames, Object *ob, const BCAnimationCurve &curve)
{
	curve.get_frames(frames); // TODO: get frames from curve...
}

bool BCAnimationSampler::get_matrix_set(BCMatrixMap &matrices, Object *ob, Bone *bone)
{
	const BCSampleKey key(ob, bone);
	sample_data.get_matrices(key, matrices);
	return is_flat_line(matrices);
}

bool BCAnimationSampler::get_matrix_set(BCMatrixMap &matrices, Object *ob)
{
	const BCSampleKey key(ob);
	sample_data.get_matrices(key, matrices);
	return is_flat_line(matrices);
}

const float BCAnimationSampler::get_value(const BCMatrix &mat, const std::string &path, const int array_index) const
{
	std::string channel = bc_string_after(path, '.');

	if (channel == "location") {
		const float(&loc)[3] = mat.location();
		return loc[array_index];
	}
	else if (channel == "scale") {
		const float(&size)[3] = mat.scale();
		return size[array_index];
	}
	else if (channel == "rotation_euler") {
		const float(&rot)[3] = mat.rotation();
		return rot[array_index];
	}
	else
	{
		const float(&quat)[4] = mat.quat();
		return quat[array_index];
	}
}

/*
   Add sampled values to FCurve 
   If no FCurve exists, create a temporary FCurve;
   Note: The temporary FCureve will later be removed when the
   BCAnimationSampler is removed (by its destructor)
*/
void BCAnimationSampler::add_value_set(
	BCMatrixMap &matrices, 
	BCAnimationCurve &curve, 
	BC_export_animation_type animation_type)
{
	std::string rna_path = curve.get_rna_path();
	int array_index = curve.get_array_index();

	BCMatrixMap::iterator it;
	for (it = matrices.begin(); it != matrices.end(); ++it) {
		const int frame_index = it->first;
		if (animation_type == BC_ANIMATION_TYPE_SAMPLE || curve.is_keyframe(frame_index)) {
			const BCMatrix *matrix = it->second;
			float val = get_value(*matrix, rna_path, array_index);
			curve.add_value(val, frame_index);
		}
	}
	curve.remove_unused_keyframes();
	curve.calchandles();
}

/*
   We assume here that the curve is populated with
   the to be exported keyframes.
   I.E. add_value_set() must have been called before
*/
const bool BCAnimationSampler::get_value_set(BCValues &values, BCFrames &frames, BCAnimationCurve &curve)
{
	values.clear();
	curve.get_values(values);
	return is_flat_line(values);
}

void BCAnimationSampler::generate_transform(
	const std::string prep,
	const std::string path,
	const int index,
	const BC_animation_curve_type type, BCAnimationCurveMap &curves)
{
	std::string rna_path = prep + path;
	CurveKey key(rna_path, index);

	BCAnimationCurveMap::const_iterator it = curves.find(key);
	if (it == curves.end()) {
		curves[key].init(type, rna_path, index);
	}
}

void BCAnimationSampler::generate_transforms(
	const std::string prep,
	const BC_animation_curve_type type,
	BCAnimationCurveMap &curves)
{
	generate_transform(prep, "location", 0, type, curves);
	generate_transform(prep, "location", 1, type, curves);
	generate_transform(prep, "location", 2, type, curves);
	generate_transform(prep, "rotation_euler", 0, type, curves);
	generate_transform(prep, "rotation_euler", 1, type, curves);
	generate_transform(prep, "rotation_euler", 2, type, curves);
	generate_transform(prep, "scale", 0, type, curves);
	generate_transform(prep, "scale", 1, type, curves);
	generate_transform(prep, "scale", 2, type, curves);
}

void BCAnimationSampler::generate_transforms(Bone *bone, BCAnimationCurveMap &curves)
{
	std::string prep = "pose.bones[\"" + std::string(bone->name) + "\"].";
	generate_transforms(prep, BC_ANIMATION_CURVE_TYPE_BONE, curves);

	for (Bone *child = (Bone *)bone->childbase.first; child; child = child->next)
		generate_transforms(child, curves);
}

void BCAnimationSampler::get_curves(BCAnimationCurveMap &curves, Object *ob)
{
	BC_animation_curve_type curve_type = BC_ANIMATION_CURVE_TYPE_OBJECT;

	bAction *action = bc_getSceneObjectAction(ob);
	if (action) {
		FCurve *fcu = (FCurve *)action->curves.first;

		for (; fcu; fcu = fcu->next) {
			if (ob->type == OB_ARMATURE) {
				char *boneName = BLI_str_quoted_substrN(fcu->rna_path, "pose.bones[");
				if (boneName) {
					curve_type = BC_ANIMATION_CURVE_TYPE_BONE;
				}
			}

			/* Adding action curves */
			CurveKey key(fcu->rna_path, fcu->array_index);
			curves[key].init(curve_type, fcu);
		}
	}

	/* Add missing curves */
	curve_type = BC_ANIMATION_CURVE_TYPE_OBJECT;
	generate_transforms(EMPTY_STRING, curve_type, curves);
	if (ob->type == OB_ARMATURE) {
		bArmature *arm = (bArmature *)ob->data;
		for (Bone *root_bone = (Bone *)arm->bonebase.first; root_bone; root_bone = root_bone->next)
			generate_transforms(root_bone, curves);
	}


	action = NULL;
	if (ob->type == OB_CAMERA) {
		action = bc_getSceneCameraAction(ob);
		curve_type = BC_ANIMATION_CURVE_TYPE_CAMERA;
	}
	else if (ob->type == OB_LAMP) {
		action = bc_getSceneLampAction(ob);
		curve_type = BC_ANIMATION_CURVE_TYPE_LIGHT;
	}

	if (action) {
		/* Add lamp action or Camera action */
		FCurve *fcu = (FCurve *)action->curves.first;
		for (; fcu; fcu = fcu->next) {
			CurveKey key(fcu->rna_path, fcu->array_index);
			curves[key].init(curve_type, fcu);
		}
	}
}

/* ==================================================================== */

void BCSampleFrame::add(Object *ob, BCMatrix &matrix)
{
	sampleMap[BCSampleKey(ob)] = matrix;
}

/* Add a new Bone to this map with the given Matrix*/
void BCSampleFrame::add(Object *ob, Bone *bone, BCMatrix &matrix)
{
	sampleMap[BCSampleKey(ob, bone)] = matrix;
}

/* Get the matrix for the given key, returns Unity when the key does not exist */
const BCMatrix &BCSampleFrame::matrix(const BCSampleKey key) const
{
	BCSamplesMap::const_iterator it = sampleMap.find(key);
	if (it == sampleMap.end()) {
		return UNIT_MATRIX;
	}
	return it->second;
}

/* Get the matrix for the given Object, returns Unity when the Objewct is not sampled */
const BCMatrix &BCSampleFrame::matrix(Object *ob) const
{
	const BCSampleKey key(ob);
	return matrix(key);
}

/* Get the matrix for the given Bone, returns Unity when the Objewct is not sampled */
const BCMatrix &BCSampleFrame::matrix(Object *ob, Bone *bone) const
{
	const BCSampleKey key(ob, bone);
	return matrix(key);
}

/* Check if the key is in this BCSampleFrame */
const bool BCSampleFrame::contains(const BCSampleKey &key) const
{
	return sampleMap.find(key) != sampleMap.end();
}

/* Check if the Object is in this BCSampleFrame */
const bool BCSampleFrame::contains(Object *ob) const
{
	const BCSampleKey key(ob);
	return contains(key);
}

/* Check if the Bone is in this BCSampleFrame */
const bool BCSampleFrame::contains(Object *ob, Bone *bone) const
{
	const BCSampleKey key(ob, bone);
	return contains(key);
}

/* Return the BCSampleMap for this BCSampleFrame */
const BCSamplesMap &BCSampleFrame::get_samples() const
{
	return sampleMap;
}

/* ==================================================================== */

/* Add object for frame. Creates a new BCSampleFrame if it does not yet exist */
BCSampleFrame &BCSampleFrames::add(Object *ob, BCMatrix &matrix, int frame_index)
{
	BCSampleFrame &frame = sample_frames[frame_index];
	frame.add(ob, matrix);
	return frame;
}


/* Add object+bone for frame. Creates a new BCSampleFrame if it does not yet exist */
BCSampleFrame &BCSampleFrames::add(Object *ob, Bone *bone, BCMatrix &matrix, int frame_index)
{
	BCSampleFrame &frame = sample_frames[frame_index];
	frame.add(ob, bone, matrix);
	return frame;
}

/* ====================================================== */
/* Below are the getters which we need to export the data */
/* ====================================================== */

/* Return either the BCSampleFrame or nullptr if frame does not exist*/
BCSampleFrame * BCSampleFrames::get_frame(int frame_index)
{
	BCSampleFrameMap::iterator it = sample_frames.find(frame_index);
	BCSampleFrame *frame = (it == sample_frames.end()) ? nullptr : &it->second;
	return frame;
}

/* Return a list of all frames that need to be sampled */
const int BCSampleFrames::get_frames(std::vector<int> &frames) const
{
	frames.clear(); // safety;
	BCSampleFrameMap::const_iterator it;
	for (it = sample_frames.begin(); it != sample_frames.end(); ++it) {
		frames.push_back(it->first);
	}
	return frames.size();
}

const int BCSampleFrames::get_frames(BCSampleKey &key, std::vector<int> &frames) const
{
	frames.clear(); // safety;
	BCSampleFrameMap::const_iterator it;
	for (it = sample_frames.begin(); it != sample_frames.end(); ++it) {
		const BCSampleFrame &frame = it->second;
		if (frame.contains(key)) {
			frames.push_back(it->first);
		}
	}
	return frames.size();
}

const int BCSampleFrames::get_frames(BCSampleKey &key, BCFrames &frames) const
{
	frames.clear(); // safety;
	BCSampleFrameMap::const_iterator it;
	for (it = sample_frames.begin(); it != sample_frames.end(); ++it) {
		const BCSampleFrame &frame = it->second;
		if (frame.contains(key)) {
			frames.push_back(it->first);
		}
	}
	return frames.size();
}

const int BCSampleFrames::get_matrices(const BCSampleKey &key, BCMatrixMap &matrices) const
{
	matrices.clear(); // safety;
	BCSampleFrameMap::const_iterator it;
	for (it = sample_frames.begin(); it != sample_frames.end(); ++it) {
		const BCSampleFrame &frame = it->second;
		if (frame.contains(key)) {
			const BCMatrix &matrix = frame.matrix(key);
			matrices[it->first] = &matrix;
		}
	}
	return matrices.size();
}

/* For convenience */
const int BCSampleFrames::get_frames(Object *ob, std::vector<int> &frames) const
{
	BCSampleKey key(ob);
	return get_frames(key, frames);
}

const int BCSampleFrames::get_frames(Object *ob, Bone *bone, std::vector<int> &frames) const
{
	BCSampleKey key(ob, bone);
	return get_frames(key, frames);
}


const int BCSampleFrames::get_matrices(Object *ob, BCMatrixMap &matrices) const
{
	BCSampleKey key(ob);
	return get_matrices(key, matrices);
}

const int BCSampleFrames::get_matrices(Object *ob, Bone *bone, BCMatrixMap &matrices) const
{
	BCSampleKey key(ob, bone);
	return get_matrices(key, matrices);
}
