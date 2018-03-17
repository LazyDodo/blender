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

#ifndef __ANIMATION_CURVE_CONTAINER_H__
#define __ANIMATION_CURVE_CONTAINER_H__

#include "BCAnimationCurve.h"

class BCSamplePoint {

private:

	Object * ob = nullptr;
	Bone *pose_bone = nullptr;
	BCAnimationCurve *curve = nullptr;

	BCMatrix matrix; /* Local matrix, by default unit matrix, will be set when sampling */

public:

	BCSamplePoint(Object *ob);
	BCSamplePoint(Object *ob, BCAnimationCurve *curve, int index);
	BCSamplePoint(Object *ob, Bone *bone);
	~BCSamplePoint();

	Object *get_object() const;
	Bone *get_bone() const;
	BCAnimationCurve *get_curve() const;
	const BCMatrix &get_matrix() const;
	const std::string get_path() const;

	const BCMatrix &set_matrix(BCMatrix &matrix);
	const BCMatrix &set_matrix(float(&mat)[4][4]);
	const BCMatrix &set_matrix(double(&mat)[4][4]);

};


typedef std::map<CurveKey, BCAnimationCurve> BCAnimationCurves;
typedef std::map<Object *, BCAnimationCurves> BCAnimationObjectMap;
typedef std::map<int, std::vector<BCSamplePoint>> BCSamplePointMap;

class BCAnimationCurveContainer {
private:

	void create_curves(Object *ob);

	BCAnimationObjectMap animated_objects; // list of objects for animating
	BCSamplePointMap sample_frames; // list of frames where objects need to be sampled

	std::vector<BCSamplePoint> &getFrameInfos(int frame_index);
	void add_sample_point(BCSamplePoint &point, int frame_index);
	void enable_fcurves(bAction *act, char *bone_name);
	bool bone_matrix_local_get(Object *ob, Bone *bone, float (&mat)[4][4], bool for_opensim);
	bool is_flat_line(std::vector<BCMatrix> &matrices) const;
	bool is_flat_line(std::vector<float> &values) const;
	static bool is_animated_by_constraint(Object *ob, ListBase *conlist, std::set<Object *> &animated_objects);

	void generate_transform(Object *ob, std::string prep, std::string path, int index, BC_animation_curve_type type, BCAnimationCurves &curves);
	void generate_transforms(Object *ob, std::string prep, BC_animation_curve_type type, BCAnimationCurves &curves);
	void generate_transforms(Object *ob, Bone *bone, BCAnimationCurves &curves);

public:

	BCAnimationCurveContainer();

	void addObject(Object *obj);
	
	void sampleMain(Scene *scene, 
		BC_export_transformation_type atm_type,
		bool for_opensim);

	void sampleScene(Scene *scene, 
		BC_export_transformation_type atm_type,
		int sampling_rate,
		bool for_opensim,
		bool keyframe_at_end = true ); // generate keyframes for frames use timeline boundaries

	void setup_curves(Object *ob, BCAnimationCurves &curves);
	void create_sample_frames_from_keyframes();
	void create_sample_frames_generated(float sfra, float efra, int sampling_rate, int keyframe_at_end);
	void get_frame_set(std::vector<float> &frames) const;
	void get_frame_set(std::vector<float> &frames, Object *ob) const;
	void get_frame_set(std::vector<float> &frames, Object *ob, Bone *bone) const;
	void get_frame_set(std::vector<float> &frames, Object *ob, const BCAnimationCurve &curve) const;

	bool get_matrix_set(std::vector<BCMatrix> &matrices, Object *ob) const;
	bool get_matrix_set(std::vector<BCMatrix> &matrices, Object *ob, Bone *bone) const;

	/* Returns the stride (number of values per frame), returns negative stride if flat_line*/
	bool get_value_set(std::vector<float> &values, const BCAnimationCurve &curve) const;
	float get_value(const BCMatrix &mat, std::string &path, int array_index) const;

	const BCAnimationCurves &get_curves(Object *ob) const;
	static void get_animated_subset(std::set<Object *> &animated_objects, LinkNode *export_set);
	static void find_indirect_animated(std::set<Object *> &animated_objects, std::set<Object *> &candidates );

};


#endif