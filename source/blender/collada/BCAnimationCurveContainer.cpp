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
#include "BLI_listbase.h"
#include "DNA_scene_types.h"
}

static std::string EMPTY_STRING;
static BCAnimationCurves BCEmptyAnimationCurves;


BCAnimationCurveContainer::BCAnimationCurveContainer()
{
}

BCAnimationCurveContainer::~BCAnimationCurveContainer()
{
	clear_cache();
}

void BCAnimationCurveContainer::clear_cache()
{

}

void BCAnimationCurveContainer::clear_cache(Object *ob)
{

}

void BCAnimationCurveContainer::create_curves(Object *ob)
{

}

void BCAnimationCurveContainer::addObject(Object *ob)
{
	BCAnimationCurves curves;
	cached_objects[ob]=curves;
}

bool BCAnimationCurveContainer::bone_matrix_local_get(Object *ob, Bone *bone, float(&mat)[4][4], bool for_opensim)
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

void BCAnimationCurveContainer::sampleMain(Scene *scene,
	BC_export_transformation_type atm_type,
	bool for_opensim)
{
	std::map<int, std::vector<BCSamplePoint>>::iterator frame;
	for (frame = sample_frames.begin(); frame != sample_frames.end(); frame++) {
		int frame_index = frame->first;
		std::vector<BCSamplePoint> &sample_points = frame->second;

		bc_update_scene(scene, frame_index);

		for (int spi = 0; spi < sample_points.size(); spi++) {
			BCSamplePoint &point = sample_points[spi];
			Object *ob = point.get_object();
			float mat[4][4];

			if (ob->type == OB_ARMATURE) {
				/* For Armatures we need to check if this maybe is a pose sample point*/
				Bone *bone = point.get_bone();
				if (bone) {
					if (bone_matrix_local_get(ob, bone, mat, for_opensim)) {
						BCMatrix matrix(mat);
						point.set_matrix(matrix); // XXX maybe not needed
						BCAnimationCurve *curve = point.get_curve();
						if (curve) {
							curve->add_value(matrix, frame_index);
						}
					}
					continue;
				}
			}

			/* When this SamplePoint is not for a Bone, 
			 * then we just store the Object local matrix here
			 */

			BKE_object_matrix_local_get(ob, mat);
			BCMatrix matrix(mat);
			point.set_matrix(matrix); // XXX maybe not needed

			BCAnimationCurve *curve = point.get_curve();
			if (curve) {
				curve->add_value(matrix, frame_index);
			}
		}
	}
}

/*
* enable fcurves driving a specific bone, disable all the rest
* if bone_name = NULL enable all fcurves
*/
void BCAnimationCurveContainer::enable_fcurves(bAction *act, char *bone_name)
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

void BCAnimationCurveContainer::sampleScene(Scene *scene,
	BC_export_transformation_type atm_type,
	int sampling_rate,
	bool for_opensim,
	bool keyframe_at_end)
{
	if (sampling_rate > 0)
		create_sample_frames_generated(scene->r.sfra, scene->r.efra, sampling_rate, keyframe_at_end);
	else
		create_sample_frames_from_keyframes();

	sampleMain(scene, atm_type, for_opensim);
}

std::vector<BCSamplePoint> &BCAnimationCurveContainer::getFrameInfos(int frame_index)
{
	std::map<int, std::vector<BCSamplePoint>>::iterator frames = sample_frames.find(frame_index);
	if (frames == sample_frames.end()) {
		std::vector<BCSamplePoint> sample_points;
		sample_frames[frame_index] = sample_points;
	}
	return sample_frames[frame_index];
}


void BCAnimationCurveContainer::add_sample_point(BCSamplePoint &point, int frame_index)
{
	std::vector<BCSamplePoint> &frame_infos = getFrameInfos(frame_index);
	frame_infos.push_back(point);
}

void BCAnimationCurveContainer::setup_curves(Object *ob, BCAnimationCurves &curves) {
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

			BCAnimationCurve curve(fcu, ob, curve_type);
			curves.push_back(curve);
		}
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
		FCurve *fcu = (FCurve *)action->curves.first;
		for (; fcu; fcu = fcu->next) {
			BCAnimationCurve curve(fcu, ob, curve_type);
			curves.push_back(curve);
		}
	}
}
/*
* loop over all cached objects
*     loop over all fcurves
*         record all keyframes
* 
* The vector sample_frames finally contains a list of vectors
* where each vector contains a list of SamplePoints which
* need to be processed when evaluating the animation.
*/
void BCAnimationCurveContainer::create_sample_frames_from_keyframes()
{
	sample_frames.clear();
	BCAnimationObjectMap::iterator it;
	for (it = cached_objects.begin(); it != cached_objects.end(); it++) {

		Object *ob = it->first;
		BCAnimationCurves &curves = it->second;

		setup_curves(ob, curves);

		BCAnimationCurves::iterator fit;
		for (fit=curves.begin(); fit!= curves.end(); ++fit) {
			BCAnimationCurve &curve = *fit;
			const FCurve *fcu = curve.get_fcurve();
			for (unsigned int i = 0; i < fcu->totvert; i++) {
				int frame_index = int(fcu->bezt[i].vec[1][0]); // X value = ctime
				float frame_value = fcu->bezt[i].vec[1][1]; // Y original value
				curve.add_value(frame_value, frame_index); // create a placeholder

				// probably no longer needed:
				BCSamplePoint sample_point(ob, &curve, i);
				add_sample_point(sample_point, frame_index);
			}
		}
	}

}

/*
* loop over all cached objects
*     loop over active action using a stesize of sampling_rate
*         record all frames
*
* The vector sample_frames finally contains a list of vectors
* where each vector contains a list of SamplePoints which
* need to be processed when evaluating the animation.
* Note: The FCurves of the objects will not be used here.
*/
void BCAnimationCurveContainer::create_sample_frames_generated(float sfra, float efra, int sampling_rate, int keyframe_at_end)
{
	sample_frames.clear();
	BCAnimationObjectMap::iterator it;
	for (it = cached_objects.begin(); it != cached_objects.end(); it++) {

		Object *ob = it->first;
		BCAnimationCurves &curves = it->second;

		setup_curves(ob, curves);

		BCAnimationCurves::iterator fit;
		for (fit = curves.begin(); fit != curves.end(); ++fit) {
			BCAnimationCurve &curve = *fit;

			float f = sfra;
			do {
				int frame_index = int(f);
				float frame_value = 0; // don't know yet;
				curve.add_value(frame_value, frame_index); // create a placeholder

				// probably no longer needed:
				BCSamplePoint sample_point(ob, &curve, frame_index);
				add_sample_point(sample_point, frame_index);

#if 0
				//TODO: add more curves for objects and bones which are moved by constraints:

				if (ob && ob->type == OB_ARMATURE) {
					LISTBASE_FOREACH(bPoseChannel *, pchan, &ob->pose->chanbase) {
						BCSamplePoint bone_point(ob, pchan->bone);
						add_sample_point(bone_point, frame_index);
					}
				}
#endif

				if (f == efra)
					break;
				f += sampling_rate;
				if (f > efra)
					if (keyframe_at_end)
						f = efra; // make sure the last frame is always exported
					else
						break;
			} while (true);
		}
	}
}

void BCAnimationCurveContainer::get_frame_set(std::vector<float> &frames) const
{
	frames.clear();
	std::map<int, std::vector<BCSamplePoint>>::const_iterator it;
	for (it = sample_frames.begin(); it != sample_frames.end(); ++it) {
		int frame_index = it->first;
		frames.push_back(frame_index);
	}
}

void BCAnimationCurveContainer::get_frame_set(std::vector<float> &frames, Object *ob, Bone *bone) const
{
	frames.clear();
	std::map<int, std::vector<BCSamplePoint>>::const_iterator mit;
	for (mit = sample_frames.begin(); mit != sample_frames.end(); ++mit) {
		int frame_index = mit->first;
		const std::vector<BCSamplePoint> &points = mit->second;

		std::vector<BCSamplePoint>::const_iterator pit;
		for (pit = points.begin(); pit != points.end(); ++pit) {
			const BCSamplePoint &point = *pit;
			if (point.get_bone() == bone && point.get_object() == ob) {
				frames.push_back(frame_index);
				break;
			}
		}
		
	}
}

/*
Combine all sample points from all ColladaCurves defined for this object
*/
void BCAnimationCurveContainer::get_frame_set(std::vector<float> &frames, Object *ob) const
{
	std::set<float>collection;
	BCAnimationCurves::const_iterator it;

	const BCAnimationCurves &curves = get_curves(ob);
	for (it = curves.begin(); it != curves.end(); ++it) {
		const BCAnimationCurve &curve = *it;
		curve.get_frames(collection);
	}

	frames.clear();
	std::copy(collection.begin(), collection.end(), std::back_inserter(frames));
}


void BCAnimationCurveContainer::get_frame_set(std::vector<float> &frames, Object *ob, const BCAnimationCurve &curve) const
{
	frames.clear();
	curve.get_frames(frames);
}

bool BCAnimationCurveContainer::is_flat_line(std::vector<BCMatrix> &values) const
{
	static float MIN_DISTANCE = 0.00001;

	if (values.size() < 2)
		return true; // need at least 2 entries to be not flat

	BCMatrix &refmat = values[0];
	for (int index = 1; index < values.size(); index++) {
		BCMatrix &mat = values[index];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (!mat.in_range(refmat, MIN_DISTANCE))
					return false;
			}
		}
	}
	return true;
}

bool BCAnimationCurveContainer::is_flat_line(std::vector<float> &values) const
{
	static float MIN_DISTANCE = 0.00001;

	if (values.size() < 2)
		return true; // need at least 2 entries to be not flat

	float ref = values[0];
	for (int index = 1; index < values.size(); index++) {
		float val = values[index];
		if (fabs(val - ref) > MIN_DISTANCE) {
			return false;
		}
	}
	return true;
}

float BCAnimationCurveContainer::get_value(const BCMatrix &mat, std::string &path, int array_index) const
{
	std::string channel = bc_string_after(path, '.');

	if (channel == "location") {
		const float (&loc)[3] = mat.location();
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

const BCAnimationCurves &BCAnimationCurveContainer::get_curves(Object *ob) const
{
	BCAnimationObjectMap::const_iterator it = cached_objects.find(ob);
	if (it == cached_objects.end()) {
		return BCEmptyAnimationCurves; // empty vector returned as const -> not modifiable
	}
	return it->second;
}

bool BCAnimationCurveContainer::get_value_set(std::vector<float> &values, const BCAnimationCurve &curve) const
{
	values.clear();
	curve.get_values(values);
	return is_flat_line(values);

}

bool BCAnimationCurveContainer::get_matrix_set(std::vector<BCMatrix> &matrices, Object *ob, Bone *bone) const
{
	matrices.clear();

	std::map<int, std::vector<BCSamplePoint>>::const_iterator mit;
	for (mit = sample_frames.begin(); mit != sample_frames.end(); ++mit) {
		int frame_index = mit->first;
		const std::vector<BCSamplePoint> &points = mit->second;

		std::vector<BCSamplePoint>::const_iterator pit;
		for (pit = points.begin(); pit != points.end(); ++pit) {
			const BCSamplePoint &point = *pit;
			if (point.get_bone() == bone && point.get_object() == ob) {
				const BCMatrix &mat = point.get_matrix();
				matrices.push_back(mat);
				break;
			}
		}
	}

	return is_flat_line(matrices);
}

bool BCAnimationCurveContainer::get_matrix_set(std::vector<BCMatrix> &matrices, Object *ob) const
{
	return get_matrix_set(matrices, ob, NULL);
}

/* ========================================================================== */

BCSamplePoint::BCSamplePoint(Object *ob)
{
	this->curve = NULL;
	this->ob = ob;
	this->pose_bone = NULL;
}

BCSamplePoint::BCSamplePoint(Object *ob, BCAnimationCurve *curve, int index)
{
	this->curve = curve;
	this->ob = ob;
	this->pose_bone = NULL;

	/* Further elaborate on what this Fcurve is doing by checking
	 * its rna_path
     */

	if (ob && ob->type == OB_ARMATURE) {
		const FCurve *fcu = curve->get_fcurve();
		if (fcu) {
			char *boneName = BLI_str_quoted_substrN(fcu->rna_path, "pose.bones[");
			bPose *pose = ob->pose;
			if (boneName) {
				bPoseChannel *pchan = BKE_pose_channel_find_name(pose, boneName);
				this->pose_bone = pchan->bone;
			}
		}
	}
}

BCSamplePoint::BCSamplePoint(Object *ob, Bone *bone)
{
	this->curve = NULL;
	this->ob = ob;
	this->pose_bone = bone;
	//this->path = "pose.bones[\"" + id_name(bone) + "\"].matrix";
}

const BCMatrix &BCSamplePoint::get_matrix() const
{
	return matrix;
}

const BCMatrix &BCSamplePoint::set_matrix(BCMatrix &mat)
{
	this->matrix.set_matrix(mat);
	return this->matrix;
}

const BCMatrix &BCSamplePoint::set_matrix(float(&mat)[4][4])
{
	matrix.set_matrix(mat);
	return this->matrix;
}

const BCMatrix &BCSamplePoint::set_matrix(double(&mat)[4][4])
{
	matrix.set_matrix(mat);
	return this->matrix;
}

BCSamplePoint::~BCSamplePoint()
{
	int x = 0;
}

Object *BCSamplePoint::get_object() const
{
	return this->ob;
}

Bone *BCSamplePoint::get_bone() const
{
	return this->pose_bone;
}

BCAnimationCurve *BCSamplePoint::get_curve() const
{
	return curve;
}


const std::string BCSamplePoint::get_path() const
{
	if (!curve)
		return "";

	const FCurve *fcu = curve->get_fcurve();
	if (!fcu)
		return "";

	return std::string(fcu->rna_path);
}
