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

extern "C" {
#include "BKE_action.h"
#include "DNA_action_types.h"
}

/* Usefull typedefs ============================================== */
typedef std::map<Object *, BCFrameSet> BCAnimatedObjectMap;
typedef std::map<int, const BCSample *> BCFrameSampleMap;
typedef std::map<CurveKey, BCAnimationCurve> BCAnimationCurveMap;
typedef std::map<Object *, BCAnimationCurveMap> BCAnimationObjectMap;

/* =============================================================== */


/* ============================================================== */
typedef std::map<BCSampleKey, BCSample *> BCSampleKeysMap;

class BCSampleFrame {

	/*
	Each frame on the timeline that needs to be sampled will have
	one BCSampleFrame where we collect all elements that need to be sampled.
	This can be Objects or Bones.
	Those elements are stored in a BCSampleMap which uses
	a BCSampleKey to identify the sampled item and a BCMatrix which contains
	the transform data for the item. Note that one item can have
	multiple FCurves. However all those FCurves can be feeded by the BCMatrix.

	TODO: expand this to Materials, Lights, Cameras...
	*/

private:
	BCSampleKeysMap sampleMap;

public:

	~BCSampleFrame()
	{
		BCSampleKeysMap::iterator it;
		for (it = sampleMap.begin(); it != sampleMap.end(); ++it) {
			BCSample *sample = it->second;
			delete sample;
		}
		sampleMap.clear();
	}

	/* Add a new Object to this map with the given Matrix*/
	/* Add a new Bone to this map with the given Matrix*/
	/* Add a new Material to this map with the given Matrix*/
	BCSample &add(Object *ob, BCMatrix &mat);
	BCSample &add(Object *ob, Bone *bone, BCMatrix &mat);

	/* Get the matrix for the given key, returns Unity when the key does not exist */
	/* Get the matrix for the given Object, returns Unity when the Objewct is not sampled */
	/* Get the matrix for the given Bone, returns Unity when the Object is not sampled */
	/* Get the matrix for the given Material, returns Unity when the Object is not sampled */
	const BCSample &get_sample(const BCSampleKey key) const;
	const BCSample &get_sample(Object *ob) const;
	const BCSample &get_sample(Object *ob, Bone *bone) const;
	const BCSample &get_sample(Object *ob, Material *mat) const;

	/* Check if the key is in this BCSampleFrame */
	/* Check if the Object is in this BCSampleFrame */
	/* Check if the Bone is in this BCSampleFrame */
	/* Check if the Material is in this BCSampleFrame */
	const bool contains(const BCSampleKey &key) const;
	const bool contains(Object *ob) const;
	const bool contains(Object *ob, Bone *bone) const;
	const bool contains(Object *ob, Material *mat) const;

	/* Return the BCSampleMap for this BCSampleFrame */
	const BCSampleKeysMap &get_samples() const;
};

typedef std::map<int, BCSampleFrame> BCSampleFrameMap;

class BCSampleFrames {

	/*
	An Animation is made of multiple keyframes or sample frames 
	within a timeline.  When we want to export the animation we
	first need to resample it fully to resolve things like:

	- animations by constraints
	- animations by drivers
	- resampled animations

	For this ourpose we need to step through the entire animation and
	then sample each frame that contains at least one keyFrame or 
	sampleFrame. Then for each frame we have to store the transform
	information for all exported objects.

	This is HERE! The BCSampleFrames class is a collector for all that information.
	The basic idea is:

	for each frame in the scene action:
	    for each object(and bone) that needs to be sampled in this frame:
			sample_frames.add(ob, matrix, frame_index)
			for each bone (when its an armature):
				sample_frames.add(object, bone, matrix, frame_index)
	
	Once the scene is fully sampled, we can get all export data from the
	BCSampleFrames instance without need to reevaluate things over and
	over again.
	*/

private:
	BCSampleFrameMap sample_frames;

public:

	~BCSampleFrames()
	{
		int x = 0;
	}

	/* Add object for frame. Creates a new BCSampleFrame if it does not yet exist */
	/* Add object+bone for frame. Creates a new BCSampleFrame if it does not yet exist */
	BCSample &add(Object *ob, BCMatrix &mat, int frame_index);
	BCSample &add(Object *ob, Bone *bone, BCMatrix &mat, int frame_index);

	/* ====================================================== */
	/* Below are the getters which we need to export the data */
	/* ====================================================== */

	/* Return either the BCSampleFrame or nullptr if frame does not exist*/
	BCSampleFrame * get_frame(int frame_index);

	/* Return a list of all frames that need to be sampled */
	const int get_frames(std::vector<int> &frames) const;
	const int get_frames(BCSampleKey &key, BCFrames &frames) const;
	const int get_matrices(const BCSampleKey &key, BCFrameSampleMap &matrices) const;
};

class BCAnimationSampler {
private:
	bContext *mContext;
	BCSampleFrames sample_data;
	BCAnimatedObjectMap objects;


	void generate_transform(
		const std::string prep,
		const std::string path,
		const int index,
		const BC_animation_curve_type type,
		BCAnimationCurveMap &curves);

	void generate_transforms(
		const std::string prep,
		const BC_animation_curve_type type,
		BCAnimationCurveMap &curves);

	void generate_transforms(Bone *bone, BCAnimationCurveMap &curves);

	/* Helper methods */
	static void enable_fcurves(bAction *act, char *bone_name);
	static bool bone_matrix_local_get(Object *ob,
		Bone *bone,
		float(&mat)[4][4],
		bool for_opensim);

public:

	BCAnimationSampler(bContext *C);
	~BCAnimationSampler()
	{
		int x = 0;
	}

	void add_object(Object *ob);

	void sample_scene(Scene *scene,
		int sampling_rate,
		int keyframe_at_end,
		bool for_opensim,
		bool keep_keyframes,
		BC_export_animation_type export_animation_type);

	bool is_flat_line(BCFrameSampleMap &values) const;
	bool is_flat_line(std::vector<float> &values) const;

	/* =========================================================================== */
	/* The Getters ... */
	/* =========================================================================== */

	void get_frame_set(BCFrames &frames, Object *ob);
	void get_frame_set(BCFrames &frames, Object *ob, Bone *bone);
	void get_frame_set(BCFrames &frames, Object *ob, const BCAnimationCurve &curve);
	bool get_samples(BCFrameSampleMap &samples, Object *ob, Bone *bone);
	bool get_samples(BCFrameSampleMap &samples, Object *ob);
	

	static void get_keyframes(Object *ob, BCFrameSet &frameset);
	static void get_animated_subset(std::set<Object *> &animated_objects, LinkNode *export_set);
	static void find_depending_animated(std::set<Object *> &animated_objects, std::set<Object *> &candidates);

	/* Possibly do not belong to this class ? */

	static bool is_animated_by_constraint(Object *ob, ListBase *conlist, std::set<Object *> &animated_objects);
	static bool has_animations(Scene *sce, LinkNode *node);
	static bool has_animations(Object *ob);

	void add_value_set(BCAnimationCurve &curve, BCFrameSampleMap &matrices, BC_export_animation_type animation_type);
	const bool get_value_set(BCValues &values, BCFrames &frames, BCAnimationCurve &curve);
	void get_curves(BCAnimationCurveMap &curves, Object *ob);

};

#endif
