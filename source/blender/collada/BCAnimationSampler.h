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
* The Original Code is Copyright (C) 2008 Blender Foundation.
* All rights reserved.
*
* Contributor(s): Blender Foundation
*
* ***** END GPL LICENSE BLOCK *****
*/

#ifndef __BC_ANIMATION_CURVE_CONTAINER_H__
#define __BC_ANIMATION_CURVE_CONTAINER_H__

#include "BCAnimationCurve.h"
#include "BCSampleData.h"
#include "collada_utils.h"

extern "C" {
#include "BKE_action.h"
#include "BLI_math_rotation.h"
#include "DNA_action_types.h"
}

/* Collection of animation curves */
class BCAnimation {
public:
	BCFrameSet frame_set;
	BCAnimationCurveMap curve_map;
	Object *reference = NULL;

	~BCAnimation()
	{
		BCAnimationCurveMap::iterator it;
		for (it = curve_map.begin(); it != curve_map.end(); ++it) {
			delete it->second;
		}
		curve_map.clear();
	}
};

typedef std::map<Object *, BCAnimation> BCAnimationObjectMap;

class BCSampleFrame {

	/*
	Each frame on the timeline that needs to be sampled will have
	one BCSampleFrame where we collect sample information about all objects
	that need to be sampled for that frame.
	*/

private:
	BCSampleMap sampleMap;

public:

	~BCSampleFrame()
	{
		BCSampleMap::iterator it;
		for (it = sampleMap.begin(); it != sampleMap.end(); ++it) {
			BCSample *sample = it->second;
			delete sample;
		}
		sampleMap.clear();
	}

	BCSample &add(Object *ob);

	/* Following methods return NULL if object is not in the sampleMap*/
	const BCSample *get_sample(Object *ob) const;
	const BCMatrix *get_sample_matrix(Object *ob) const;
	const BCMatrix *get_sample_matrix(Object *ob, Bone *bone) const;

	const bool has_sample_for(Object *ob) const;
	const bool has_sample_for(Object *ob, Bone *bone) const;
};

typedef std::map<int, BCSampleFrame> BCSampleFrameMap;

class BCSampleFrameContainer {

	/*
	An Animation is made of multiple FCurve where each FCurve can have multiple
	keyframes. When we want to export the animation we
	first need to resample it fully to resolve things like:

	- animations by constraints
	- animations by drivers
	- resampled animations

	For this purpose we need to step through the entire animation and
	then sample each frame that contains at least one keyFrame or 
	sampleFrame. Then for each frame we have to store the transform
	information for all exported objects.

	This is HERE! The BCSampleFrames object stores a map of BCSampleFrame objects with the timline frame as key.
	*/

private:
	BCSampleFrameMap sample_frames;

public:

	~BCSampleFrameContainer()
	{
		int x = 0;
	}

	BCSample &add(Object *ob, int frame_index);
	BCSampleFrame * get_frame(int frame_index); // returns NULL if frame does not exist

	const int get_frames(std::vector<int> &frames) const;
	const int get_frames(Object *ob, BCFrames &frames) const;
	const int get_frames(Object *ob, Bone *bone, BCFrames &frames) const;

	const int get_samples(Object *ob, BCFrameSampleMap &samples) const;
	const int get_matrices(Object *ob, BCMatrixSampleMap &matrices) const;
	const int get_matrices(Object *ob, Bone *bone, BCMatrixSampleMap &bones) const;
};

class BCAnimationSampler {
private:
	bContext *mContext;
	BCSampleFrameContainer sample_data;
	BCAnimationObjectMap objects;

	void generate_transform(Object *ob, const BCCurveKey &key, BCAnimationCurveMap &curves);
	void generate_transforms(Object *ob, const std::string prep, const BC_animation_type type, BCAnimationCurveMap &curves);
	void generate_transforms(Object *ob, Bone *bone, BCAnimationCurveMap &curves);

	void initialize_curves(BCAnimationCurveMap &curves, Object *ob);
	void initialize_keyframes(BCFrameSet &frameset, Object *ob);

	void update_animation_curves(BCAnimation &animation, Object *ob, int frame_index);
	void check_property_is_animated(BCAnimation &animation, float *ref, float *val, std::string data_path, int length);

public:

	BCAnimationSampler(bContext *C);
	~BCAnimationSampler();

	void add_object(Object *ob);
	void add_objects(BCObjectSet &animated_subset);

	void sample_scene(Scene *scene,
		int sampling_rate,
		int keyframe_at_end,
		bool for_opensim,
		bool keep_keyframes,
		BC_export_animation_type export_animation_type);

	bool is_animated(BCMatrixSampleMap &values) const;

	BCAnimationCurveMap *get_curves(Object *ob);

	void get_object_frames(BCFrames &frames, Object *ob);
	bool get_object_samples(BCMatrixSampleMap &samples, Object *ob);
	void get_bone_frames(BCFrames &frames, Object *ob, Bone *bone);
	bool get_bone_samples(BCMatrixSampleMap &samples, Object *ob, Bone *bone);

	static void get_animated_from_export_set(std::set<Object *> &animated_objects, LinkNode &export_set);
	static void find_depending_animated(std::set<Object *> &animated_objects, std::set<Object *> &candidates);
	static bool is_animated_by_constraint(Object *ob, ListBase *conlist, std::set<Object *> &animated_objects);

	static bool has_animations(Scene *sce, LinkNode &node);
	static bool has_animations(Object *ob);

};

#endif
