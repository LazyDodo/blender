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

#ifndef __BCANIMATION_CURVE_H__
#define __BCANIMATION_CURVE_H__

#include "collada_utils.h"
#include "BCMatrix.h"

extern "C"
{
#include "MEM_guardedalloc.h"
#include "BKE_fcurve.h"
#include "BKE_armature.h"
#include "ED_keyframing.h"
}

typedef enum BC_animation_curve_type {
	BC_ANIMATION_CURVE_TYPE_UNKNOWN = -1,

	BC_ANIMATION_CURVE_TYPE_OBJECT,
	BC_ANIMATION_CURVE_TYPE_BONE,
	BC_ANIMATION_CURVE_TYPE_CAMERA,
	BC_ANIMATION_CURVE_TYPE_MATERIAL,
	BC_ANIMATION_CURVE_TYPE_LIGHT

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

	CurveKey(std::string rna_path, int array_index)
	{
		this->init(rna_path, array_index);
	}

	void init(std::string rna_path, int array_index)
	{
		this->rna_path = rna_path;
		this->array_index = array_index;
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

class BCAnimationCurve {
private:
	BC_animation_curve_type type = BC_ANIMATION_CURVE_TYPE_UNKNOWN;
	std::vector<int> sample_frames;
	CurveKey curve_key;
	bool curve_is_local_copy = false;
	FCurve *fcurve;

	void delete_fcurve(FCurve *fcu);
	FCurve *create_fcurve(int array_index, const char *rna_path);
	void create_bezt(float frame, float output);

public:
	BCAnimationCurve();
	~BCAnimationCurve();

	const BC_animation_curve_type get_channel_type() const;
	const std::string get_channel_target() const;
	const std::string get_animation_name(Object *ob) const;
	const int get_array_index() const;
	const std::string get_rna_path() const;

	const int size() const;
	const int closest_index_above(float sample_frame, int start_at) const;
	const int closest_index_below(float sample_frame) const;
	const int get_ipo(float sample_frame) const;
	Bone *get_bone(Object *ob) const;

	void add_value(const float val, const int frame);
	void init(Object *ob, BC_animation_curve_type type, std::string rna_path, int index);
	void init(Object *ob, BC_animation_curve_type type, FCurve *fcu);
	FCurve *get_edit_fcurve();
	const FCurve *get_fcurve() const;

	/*
	Pick the value from the matrix accoridng to the definition of the FCurve
	Note: This works only for "scale", "rotation", "rotation_euler" and "location"
	*/
	void add_value(BCMatrix mat, int frame);

	/*
	Return the frames of the sampled curve;
	Note: If the curve was not sampled, the
	returned vector is empty
	*/
	void get_frames(std::vector<float> &frames) const;
	void get_frames(std::set<float> &frames) const;

	/*
	Return the ctimes of the sampled curve;
	Note: If the curve was not sampled, the
	returned vector is empty
	*/
	void get_times(std::vector<float> &times, Scene *scene) const;

	/*
	Return the ctimes of the sampled curve;
	Note: If the curve was not sampled, the
	returned vector is empty
	*/
	void get_values(std::vector<float> &values) const;
	static bool is_flat_line(std::vector<float> &values);
};


#endif