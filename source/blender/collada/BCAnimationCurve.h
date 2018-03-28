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

#ifndef __BC_ANIMATION_CURVE_H__
#define __BC_ANIMATION_CURVE_H__

#include "collada_utils.h"
#include "BCSampleData.h"

extern "C"
{
#include "MEM_guardedalloc.h"
#include "BKE_fcurve.h"
#include "BKE_armature.h"
#include "ED_anim_api.h"
#include "ED_keyframing.h"
#include "ED_keyframes_edit.h"
}

typedef float(TangentPoint)[2];

typedef std::set<float> BCFrameSet;
typedef std::vector<float> BCFrames;
typedef std::vector<float> BCValues;
typedef std::vector<float> BCTimes;
typedef std::map<int, float> BCValueMap;

typedef enum BC_animation_curve_type {
	
	BC_ANIMATION_CURVE_TYPE_UNKNOWN = -1,

	BC_ANIMATION_CURVE_TYPE_OBJECT,
	BC_ANIMATION_CURVE_TYPE_BONE,
	BC_ANIMATION_CURVE_TYPE_CAMERA, // not needed, data in BC_ANIMATION_CURVE_TYPE_OBJECT
	BC_ANIMATION_CURVE_TYPE_MATERIAL,
	BC_ANIMATION_CURVE_TYPE_LIGHT // not needed, data in BC_ANIMATION_CURVE_TYPE_OBJECT

} BC_animation_curve_type;

class CurveKey {
private:
	std::string rna_path;
	int array_index;

public:

	CurveKey()
	{
		rna_path = "";
		array_index = -1;
	}

	CurveKey(const std::string rna_path, const int index)
	{
		this->init(rna_path, index);
	}

	void init(const std::string rna_path, const int index)
	{
		this->rna_path = rna_path;
		this->array_index = index;
	}

	const std::string path() const
	{
		return this->rna_path;
	}

	const int index() const
	{
		return this->array_index;
	}

	bool operator< (const CurveKey& key) const
	{
		if (rna_path == key.rna_path) {
			return key.array_index < this->array_index;
		}
		return ((rna_path < key.rna_path));
	}
};

typedef enum BCBeztDataType {
	BCBeztInTangent = 0,
	BCBeztOutTangent = 2,
	BCBeztValue = 1
} BCBeztDataType;

class BCBezTriple {
public:
	BezTriple & bezt;

	BCBezTriple(BezTriple bezt):
		bezt(bezt){}

	float get_frame()
	{
		return bezt.vec[1][0];
	}

	float get_value(BCBeztDataType type, float val[2] = nullptr)
	{
		if (val) {
			val[0] = bezt.vec[type][0];
			val[1] = bezt.vec[type][1];
		}
		return bezt.vec[type][1];
	}

};

class BCAnimationCurve {
private:
	BC_animation_curve_type type = BC_ANIMATION_CURVE_TYPE_UNKNOWN;
	BCValueMap samples;
	float min = 0;
	float max = 0;
	int tag = -1;

	CurveKey curve_key;
	bool curve_is_local_copy = false;
	FCurve *fcurve;

	float(*convert)(BCSample *); // conversion function for values

	void delete_fcurve(FCurve *fcu);
	FCurve *create_fcurve(int array_index, const char *rna_path);
	void create_bezt(float frame, float output);

public:
	BCAnimationCurve();
	BCAnimationCurve(const BCAnimationCurve &other);
	BCAnimationCurve(const BC_animation_curve_type type, FCurve *fcu);
	~BCAnimationCurve();

	const BC_animation_curve_type get_channel_type() const;
	const BC_animation_transform_type get_transform_type() const;
	const void set_transform_type(std::string path, int axis_index);

	const std::string get_channel_target() const;
	const std::string get_animation_name(Object *ob) const;
	const int get_array_index() const;
	const std::string get_rna_path() const;

	const int size() const;
	const int closest_index_above(const float sample_frame, const int start_at) const;
	const int closest_index_below(const float sample_frame) const;
	const int get_ipo(float sample_frame) const;
	Bone *get_bone(Object *ob) const;

	void calchandles();
	bool add_value(const BCSample &sample, int frame);
	void add_value(const float val, const int frame, bool modify_curve=false);
	void init(const BC_animation_curve_type type, const std::string rna_path, const int index);
	void init(const BC_animation_curve_type type, FCurve *fcu, int material_index=-1);
	void reset_values();
	FCurve *get_edit_fcurve();
	const FCurve *get_fcurve() const;

	float get_vertical_span()
	{
		return max - min;
	}

	void remove_unused_keyframes();
	void clean_handles();
	/*
	Return the frames of the sampled curve;
	Note: If the curve was not sampled, the
	returned vector is empty
	*/
	void get_sampled_frames(BCFrameSet &frames) const;

	/*
	Return the ctimes of the sampled curve;
	Note: If the curve was not sampled, the
	returned vector is empty
	*/
	void get_times(BCTimes &times, Scene *scene) const;
	const int get_tag() const;

	/*
	Return the ctimes of the sampled curve;
	Note: If the curve was not sampled, the
	returned vector is empty
	*/
	BCValueMap &get_value_map();
	void get_key_values(BCValues &values) const;
	void get_sampled_values(BCValues &values, bool fallback = true) const;
	void get_key_frames(BCFrames &frames) const;
	void get_sampled_frames(BCFrames &frames, bool fallback = true) const;
	bool is_flat();
	bool is_rot() const;
	bool is_keyframe(int frame);

	/* For convenience */
	static bool is_flat_line(BCValues &values);
	const bool is_transform_curve() const;
	static float get_time(BezTriple *bezt, Scene *scene);
	static float get_value(BezTriple *bezt, bool as_angle = false);
	static void get_in_tangent(BezTriple *bezt, Scene *scene, float point[2], bool as_angle = false);
	static void get_out_tangent(BezTriple *bezt, Scene *scene, float point[2], bool as_angle = false);
	static void get_tangent(BezTriple *bezt, Scene *scene, float point[2], bool as_angle, int index);

};

#endif
