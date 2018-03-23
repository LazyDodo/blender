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

#include "BCAnimationCurve.h"

std::map<BC_animation_transform_type, std::string> BC_ANIMATION_NAME_FROM_TYPE = {
{ BC_ANIMATION_TYPE_ROTATION, "rotation" },
{ BC_ANIMATION_TYPE_ROTATION_EULER, "rotation_euler" },
{ BC_ANIMATION_TYPE_ROTATION_QUAT, "rotation_quaternion" },
{ BC_ANIMATION_TYPE_SCALE, "scale" },
{ BC_ANIMATION_TYPE_LOCATION, "location" },

/* Materials */
{ BC_ANIMATION_TYPE_SPECULAR_COLOR, "specular" },
{ BC_ANIMATION_TYPE_DIFFUSE_COLOR, "diffuse" },
{ BC_ANIMATION_TYPE_IOR, "index_of_refraction" },
{ BC_ANIMATION_TYPE_SPECULAR_HARDNESS, "specular_hardness" },
{ BC_ANIMATION_TYPE_ALPHA, "alpha" },

/* Lamps */
{ BC_ANIMATION_TYPE_COLOR, "color" },
{ BC_ANIMATION_TYPE_FALL_OFF_ANGLE, "fall_off_angle" },
{ BC_ANIMATION_TYPE_FALL_OFF_EXPONENT, "fall_off_exponent" },
{ BC_ANIMATION_TYPE_BLENDER_DIST, "blender/blender_dist" },

/* Cameras */
{ BC_ANIMATION_TYPE_XFOV, "xfov" },
{ BC_ANIMATION_TYPE_XMAG, "xmag" },
{ BC_ANIMATION_TYPE_ZFAR, "zfar" },
{ BC_ANIMATION_TYPE_ZNEAR, "znear" },

{ BC_ANIMATION_TYPE_UNKNOWN, "" }
};

std::map<std::string, BC_animation_transform_type> BC_ANIMATION_TYPE_FROM_NAME = {
	{ "rotation", BC_ANIMATION_TYPE_ROTATION },
{ "rotation_euler", BC_ANIMATION_TYPE_ROTATION_EULER },
{ "rotation_quaternion", BC_ANIMATION_TYPE_ROTATION_QUAT },
{ "scale", BC_ANIMATION_TYPE_SCALE },
{ "location", BC_ANIMATION_TYPE_LOCATION },

/* Materials */
{ "specular_hardness", BC_ANIMATION_TYPE_SPECULAR_HARDNESS },
{ "specular_color", BC_ANIMATION_TYPE_SPECULAR_COLOR },
{ "diffuse_color", BC_ANIMATION_TYPE_DIFFUSE_COLOR },
{ "alpha", BC_ANIMATION_TYPE_ALPHA },
{ "ior", BC_ANIMATION_TYPE_IOR },

/* Lamps */
{ "color", BC_ANIMATION_TYPE_COLOR },
{ "fall_off_angle", BC_ANIMATION_TYPE_FALL_OFF_ANGLE },
{ "fall_off_exponent", BC_ANIMATION_TYPE_FALL_OFF_EXPONENT },
{ "blender/blender_dist", BC_ANIMATION_TYPE_BLENDER_DIST },
/* Lamp  RNA to animation type */
{ "spot_size", BC_ANIMATION_TYPE_FALL_OFF_ANGLE },
{ "spot_blend", BC_ANIMATION_TYPE_FALL_OFF_EXPONENT },
{ "distance", BC_ANIMATION_TYPE_BLENDER_DIST },

/* Cameras */
{ "xfov", BC_ANIMATION_TYPE_XFOV },
{ "xmag", BC_ANIMATION_TYPE_XMAG },
{ "zfar", BC_ANIMATION_TYPE_ZFAR },
{ "znear", BC_ANIMATION_TYPE_ZNEAR },
/* Camera RNA to animation type */
{ "lens", BC_ANIMATION_TYPE_XFOV },
{ "ortho_scale", BC_ANIMATION_TYPE_XMAG },
{ "clip_end", BC_ANIMATION_TYPE_ZFAR },
{ "clip_start", BC_ANIMATION_TYPE_ZNEAR }

};

void BCAnimationCurve::delete_fcurve(FCurve *fcu)
{
	free_fcurve(fcu);
}

FCurve *BCAnimationCurve::create_fcurve(int array_index, const char *rna_path)
{
	FCurve *fcu = (FCurve *)MEM_callocN(sizeof(FCurve), "FCurve");
	fcu->flag = (FCURVE_VISIBLE | FCURVE_AUTO_HANDLES | FCURVE_SELECTED);
	fcu->rna_path = BLI_strdupn(rna_path, strlen(rna_path));
	fcu->array_index = array_index;
	return fcu;
}

void BCAnimationCurve::create_bezt(float frame, float output)
{
	FCurve *fcu = get_edit_fcurve();
	BezTriple bez;
	memset(&bez, 0, sizeof(BezTriple));
	bez.vec[1][0] = frame;
	bez.vec[1][1] = output;
	bez.ipo = U.ipo_new; /* use default interpolation mode here... */
	bez.f1 = bez.f2 = bez.f3 = SELECT;
	bez.h1 = bez.h2 = HD_AUTO;
	insert_bezt_fcurve(fcu, &bez, INSERTKEY_NOFLAGS);
	calchandles_fcurve(fcu);
}

BCAnimationCurve::BCAnimationCurve()
{
	this->fcurve = nullptr;
	this->type = BC_ANIMATION_CURVE_TYPE_UNKNOWN;
	this->curve_is_local_copy = false;
}

BCAnimationCurve::BCAnimationCurve(const BC_animation_curve_type type, FCurve *fcu)
{
	init(type, fcu);
}

void BCAnimationCurve::init(const BC_animation_curve_type type, const std::string path, const int index)
{
	this->curve_key.init(path, index);
	this->fcurve = nullptr; // create_fcurve(index, path.c_str());
	this->type = type;
	this->curve_is_local_copy = false;
}

void BCAnimationCurve::init(const BC_animation_curve_type type, FCurve *fcu)
{
	this->curve_key.init(std::string(fcu->rna_path), fcu->array_index);
	this->fcurve = fcu;
	this->type = type;
	this->curve_is_local_copy = false; // make sure the curve is destroyed later;
}

BCAnimationCurve::~BCAnimationCurve()
{
	if (curve_is_local_copy && fcurve) {
		//fprintf(stderr, "removed fcurve %s\n", fcurve->rna_path);
		delete_fcurve(fcurve);
		this->fcurve = nullptr;
	}
}

const BC_animation_curve_type BCAnimationCurve::get_channel_type() const
{
	return type;
}

const BC_animation_transform_type BCAnimationCurve::get_transform_type() const
{
	BC_animation_transform_type tm_type;
	std::string target = get_channel_target();
	std::map<std::string, BC_animation_transform_type>::iterator it = BC_ANIMATION_TYPE_FROM_NAME.find(target);
	tm_type = (it != BC_ANIMATION_TYPE_FROM_NAME.end()) ? it->second : BC_ANIMATION_TYPE_UNKNOWN;
	return tm_type;
}

const std::string BCAnimationCurve::get_sid(const std::string axis_name) const
{
	std::string tm_name;
	BC_animation_transform_type tm_type = get_transform_type();
	std::map<BC_animation_transform_type, std::string>::iterator name_it = BC_ANIMATION_NAME_FROM_TYPE.find(tm_type);
	tm_name = name_it->second;

	if (axis_name != "")
		tm_name += '.' + axis_name;

	return tm_name;
}

const std::string BCAnimationCurve::get_channel_target() const
{
	const std::string path = curve_key.path();
	return bc_string_after(path, '.');
}

const std::string BCAnimationCurve::get_animation_name(Object *ob) const
{
	std::string name;

	switch (type) {
	case BC_ANIMATION_CURVE_TYPE_OBJECT:
		name = id_name(ob);
		break;

	case BC_ANIMATION_CURVE_TYPE_BONE:
		if (fcurve == NULL || fcurve->rna_path == NULL)
			name = "";
		else {
			const char *boneName = BLI_str_quoted_substrN(fcurve->rna_path, "pose.bones[");
			name = (boneName) ? std::string(boneName) : "";
		}
		break;

	case BC_ANIMATION_CURVE_TYPE_CAMERA:
		name = id_name(ob) + "-camera";
		break;

	case BC_ANIMATION_CURVE_TYPE_LIGHT:
		name = id_name(ob) + "-light";
		break;

	case BC_ANIMATION_CURVE_TYPE_MATERIAL:
		name = id_name(ob) + "-material";
		break;

	default:
		name = "";
	}
	return name;
}

const int BCAnimationCurve::get_array_index() const
{
	return curve_key.index();
}

const std::string BCAnimationCurve::get_rna_path() const
{
	return curve_key.path();
}

const int BCAnimationCurve::size() const
{
	return samples.size();
}

const int BCAnimationCurve::closest_index_above(const float sample_frame, const int start_at) const
{
	if (fcurve == NULL)
		return -1;

	const int cframe = fcurve->bezt[start_at].vec[1][0]; // inacurate!

	if (fabs(cframe - sample_frame) < 0.00001)
		return start_at;
	return (fcurve->totvert > start_at + 1) ? start_at + 1 : start_at;
}

const int BCAnimationCurve::closest_index_below(const float sample_frame) const
{
	if (fcurve == nullptr)
		return -1;

	float lower_frame = sample_frame;
	float upper_frame = sample_frame;
	int lower_index = 0;
	int upper_index = 0;

	for (int fcu_index = 0; fcu_index < fcurve->totvert; ++fcu_index) {
		upper_index = fcu_index;

		const int cframe = fcurve->bezt[fcu_index].vec[1][0]; // inacurate!
		if (cframe <= sample_frame) {
			lower_frame = cframe;
			lower_index = fcu_index;
		}
		if (cframe >= sample_frame) {
			upper_frame = cframe;
			break;
		}
	}

	if (lower_index == upper_index)
		return lower_index;

	const float fraction = float(sample_frame - lower_frame) / (upper_frame - lower_frame);
	return (fraction < 0.5) ? lower_index : upper_index;
}

const int BCAnimationCurve::get_ipo(float sample_frame) const
{
	const int index = closest_index_below(sample_frame);
	if (index < 0)
		return BEZT_IPO_BEZ;
	return fcurve->bezt[index].ipo;
}

const FCurve *BCAnimationCurve::get_fcurve() const
{
	return fcurve;
}

FCurve *BCAnimationCurve::get_edit_fcurve()
{
	if (!curve_is_local_copy) {

		if (fcurve) {
			fcurve = copy_fcurve(fcurve);
			//fprintf(stderr, "Copy to temporary fcurve %s (for editing)\n", fcurve->rna_path);
		}
		else {
			const int index = curve_key.index();
			const std::string &path = curve_key.path();
			fcurve = create_fcurve(index, path.c_str());
			//fprintf(stderr, "Create temporary fcurve %s (for editing)\n", fcurve->rna_path);
		}

		/* Replacing the pointer here is OK because the original value
		of FCurve was a const pointer into Blender territory. We do not
		touch that! We use the local copy to prepare for export.
		*/
		curve_is_local_copy = true;
	}
	return fcurve;
}

Bone *BCAnimationCurve::get_bone(Object *ob) const
{
	if (ob->type != OB_ARMATURE ||
		type != BC_ANIMATION_CURVE_TYPE_BONE ||
		!fcurve ||
		!fcurve->rna_path)
		return nullptr;

	const char *boneName = BLI_str_quoted_substrN(fcurve->rna_path, "pose.bones[");

	if (!boneName)
		return nullptr;

	Bone *bone = BKE_armature_find_bone_name((bArmature *)ob->data, boneName);
	return bone;

}

void BCAnimationCurve::calchandles()
{
	FCurve *fcu = this->get_edit_fcurve();
	if(fcu)
		calchandles_fcurve(fcu);
}

void BCAnimationCurve::remove_unused_keyframes()
{
	FCurve *fcu = get_edit_fcurve();
	if (fcu) {
		BezTriple  *bezt = fcu->bezt;
		for (int i = 0; i < fcu->totvert; bezt++, i++) {
			const int frame_index = bezt->vec[1][0];
			BCValueMap::const_iterator it = samples.find(frame_index);
			if (it == samples.end()) {
				fcu->bezt[i].f2 |= SELECT;
			}
			else {
				fcu->bezt[i].f2 &= ~SELECT;
			}
		}
		delete_fcurve_keys(fcu);
	}
}

void BCAnimationCurve::add_value(const float val, const int frame_index)
{
	FCurve *fcu = get_edit_fcurve();
	if (fcu) {
		const float eval = evaluate_fcurve(fcu, frame_index);

		/*
		* This is a bit tricky here. We actually only insert a keyframe into the FCurve
		* Preserving the current value. Then we add the frame index and the "true" value
		* into a separate value_map <frame, value>
		*
		* Reason: we need the Fcurve handles later when we want to export the values as a Bezier curve
		* XXX: This is actually not correct. But sometimes it works nicely. this needs more
		* experimenting
		*/

		insert_vert_fcurve(fcu, frame_index, eval, (eBezTriple_KeyframeType)BEZT_IPO_BEZ, INSERTKEY_NO_USERPREF);
		samples[frame_index] = val;

		if (samples.size() == 1)
			min = max = val;
		else {
			if (val < min)
				min = val;
			if (val > max)
				max = val;
		}
	}
}


/*
Pick the value from the matrix according to the definition of the FCurve
Note: This works only for "scale", "rotation", "rotation_euler" and "location"
*/
bool BCAnimationCurve::add_value(const BCSample &sample, int frame)
{
	std::string target = get_channel_target();
	const int array_index = curve_key.index();
	float val = 0;

	bool good = sample.get_value_for(target, array_index, &val);
	if (good) {
		add_value(val, frame);
	}
	return good;
}

/*
Return the frames of the sampled curve;
Note: If the curve was not sampled, the
returned vector is empty
*/

void BCAnimationCurve::get_key_frames(BCFrames &frames) const
{
	if (fcurve) {
		for (int i = 0; i < fcurve->totvert; i++) {
			const float val = fcurve->bezt[i].vec[1][0];
			frames.push_back(val);
		}
	}
}

void BCAnimationCurve::get_sampled_frames(BCFrames &frames, bool fallback) const
{
	frames.clear();
	if (samples.size() == 0 && fallback) {
		return get_key_frames(frames);
	}
	else if (samples.size() > 0) {
		BCValueMap::const_iterator it;
		for (it = samples.begin(); it != samples.end(); ++it) {
			//float val = evaluate_fcurve(fcurve, *it);
			const int val = it->first;
			frames.push_back(val);
		}
	}
}


void BCAnimationCurve::get_sampled_frames(BCFrameSet &frames) const
{
	BCValueMap::const_iterator it;
		for (it = samples.begin(); it != samples.end(); ++it) {
			const float val = it->first;
			frames.insert(val);
		}
}

/*
Return the ctimes of the sampled curve;
Note: If the curve was not sampled, the
returned vector is empty
*/
void BCAnimationCurve::get_times(BCTimes &times, Scene *scene) const
{
	BCValueMap::const_iterator it;
		for (it = samples.begin(); it != samples.end(); ++it) {
			const float time = FRA2TIME(it->first); // implicit use of scene
			times.push_back(time);
		}
}

/*
Return the ctimes of the sampled curve;
Note: If the curve was not sampled, the
returned vector is empty
*/
void BCAnimationCurve::get_key_values(BCValues &values) const
{
	if (fcurve) {
		for (int i = 0; i < fcurve->totvert; i++) {
			const float val = fcurve->bezt[i].vec[1][1];
			values.push_back(val);
		}
	}
}

void BCAnimationCurve::get_sampled_values(BCValues &values, bool fallback) const
{
	values.clear();
	if (samples.size() == 0 && fallback) {
		return get_key_values(values);
	}
	else if (samples.size() > 0) {
		BCValueMap::const_iterator it;
		for (it = samples.begin(); it != samples.end(); ++it) {
			//float val = evaluate_fcurve(fcurve, *it);
			const float val = it->second;
			values.push_back(val);
		}
	}
}


bool BCAnimationCurve::is_flat()
{
	static float MIN_DISTANCE = 0.00001;
	return fabs(max - min) < MIN_DISTANCE;
}

bool BCAnimationCurve::is_flat_line(BCValues &values)
{
	static float MIN_DISTANCE = 0.00001;

	if (values.size() < 2)
		return true; // need at least 2 entries to be not flat

	const float ref = values[0];
	for (int index = 1; index < values.size(); index++) {
		const float val = values[index];
		if (fabs(val - ref) > MIN_DISTANCE) {
			return false;
		}
	}
	return true;
}

bool BCAnimationCurve::is_rot() const
{
	return bc_startswith(get_channel_target(), "rotation");
}

bool BCAnimationCurve::is_keyframe(int frame) {
	if (this->fcurve == nullptr)
		return false;

	for (int i = 0; i < fcurve->totvert; ++i) {
		const int cframe = fcurve->bezt[i].vec[1][0];
		if (cframe == frame)
			return true;
		if (cframe > frame)
			break;
	}
	return false;
}

/* Needed for adding a BCAnimationCurve into a BCAnimationCurveSet */
inline bool operator< (const BCAnimationCurve& lhs, const BCAnimationCurve& rhs) {
	std::string lhtgt = lhs.get_channel_target();
	std::string rhtgt = rhs.get_channel_target();
	if (lhtgt == rhtgt)
	{
		const int lha = lhs.get_array_index();
		const int rha = rhs.get_array_index();
		return lha < rha;
	}
	else
		return lhtgt < rhtgt;
}
