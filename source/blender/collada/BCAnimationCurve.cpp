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

#include "BCAnimationCurve.h"

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
	insert_bezt_fcurve(fcu, &bez, 0);
	calchandles_fcurve(fcu);
}

BCAnimationCurve::BCAnimationCurve(FCurve *fcu, Object *ob, BC_animation_curve_type type) // = BC_ANIMATION_CURVE_TYPE_OBJECT )
{
	this->ob = ob;
	this->fcurve = fcu;
	this->type = type;
}

BCAnimationCurve::~BCAnimationCurve()
{
	if (curve_is_local_copy && fcurve) {
		fprintf(stderr, "removed fcurve %s\n", fcurve->rna_path);
		delete_fcurve(fcurve);
	}
}

const BC_animation_curve_type BCAnimationCurve::get_channel_type() const
{
	return type;
}

const std::string BCAnimationCurve::get_channel_target() const
{
	if (!fcurve || !fcurve->rna_path)
		return "";

	return bc_string_after(fcurve->rna_path, '.');
}

const std::string BCAnimationCurve::get_animation_name() const
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
			char *boneName = BLI_str_quoted_substrN(fcurve->rna_path, "pose.bones[");
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
	return fcurve->array_index;
}

const int BCAnimationCurve::size() const
{
	return export_values.size();
}

const int BCAnimationCurve::closest_index_above(float sample_frame, int start_at) const
{
	if (fcurve == NULL)
		return -1;

	int cframe = fcurve->bezt[start_at].vec[1][0]; // inacurate!

	if (fabs(cframe - sample_frame) < 0.00001)
		return start_at;
	return (fcurve->totvert > start_at + 1) ? start_at + 1 : start_at;
}

const int BCAnimationCurve::closest_index_below(float sample_frame) const
{
	if (fcurve == NULL)
		return -1;

	float lower_frame = sample_frame;
	float upper_frame = sample_frame;
	int lower_index = 0;
	int upper_index = 0;

	for (int fcu_index = 0; fcu_index < fcurve->totvert; ++fcu_index) {
		upper_index = fcu_index;

		int cframe = fcurve->bezt[fcu_index].vec[1][0]; // inacurate!
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

	float fraction = float(sample_frame - lower_frame) / (upper_frame - lower_frame);
	return (fraction < 0.5) ? lower_index : upper_index;
}

const int BCAnimationCurve::get_ipo(float sample_frame) const
{
	int index = closest_index_below(sample_frame);
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

		fcurve = copy_fcurve(fcurve);

		/* Replacing the pointer here is OK because the original value
		of FCurve was a const pointer into Blender territory. We do not
		touch that! We use the local copy to prepare for export.
		*/
		curve_is_local_copy = true;
		fprintf(stderr, "Copy fcurve %s (for editing)\n", fcurve->rna_path);
	}
	return fcurve;
}

Object *BCAnimationCurve::get_object() const
{
	return ob;
}

void BCAnimationCurve::add_value(const float val, const int frame)
{
	FCurve *fcu = get_edit_fcurve();
	insert_vert_fcurve(fcu, frame, val, BEZT_IPO_BEZ, INSERTKEY_NO_USERPREF);
	export_values[frame] = val;
}

/*
Pick the value from the matrix accoridng to the definition of the FCurve
Note: This works only for "scale", "rotation", "rotation_euler" and "location"
*/
void BCAnimationCurve::add_value(BCMatrix mat, int frame)
{
	std::string target = get_channel_target();
	int array_index = fcurve->array_index;
	float val;

	if (target == "location") {
		const float(&loc)[3] = mat.location();
		val = loc[array_index];
	}
	else if (target == "scale") {
		const float(&size)[3] = mat.scale();
		val = size[array_index];
	}
	else if (target == "rotation_euler") {
		const float(&rot)[3] = mat.rotation();
		val = rot[array_index];
	}
	else
	{
		const float(&quat)[4] = mat.quat();
		val = quat[array_index];
	}

	add_value(val, frame);
}

/*
Return the frames of the sampled curve;
Note: If the curve was not sampled, the
returned vector is empty
*/
void BCAnimationCurve::get_frames(std::vector<float> &frames) const
{
	std::map<int, float>::const_iterator it;
	for (it = export_values.begin(); it != export_values.end(); ++it) {
		frames.push_back(it->first);
	}
}

void BCAnimationCurve::get_frames(std::set<float> &frames) const
{
	std::map<int, float>::const_iterator it;
	for (it = export_values.begin(); it != export_values.end(); ++it) {
		frames.insert(it->first);
	}
}

/*
Return the ctimes of the sampled curve;
Note: If the curve was not sampled, the
returned vector is empty
*/
void BCAnimationCurve::get_times(std::vector<float> &times, Scene *scene) const
{
	std::map<int, float>::const_iterator it;
	for (it = export_values.begin(); it != export_values.end(); ++it) {
		float time = FRA2TIME(it->first); // implicit use of scene
		times.push_back(it->first);
	}
}

/*
Return the ctimes of the sampled curve;
Note: If the curve was not sampled, the
returned vector is empty
*/
void BCAnimationCurve::get_values(std::vector<float> &values) const
{
	std::map<int, float>::const_iterator it;
	for (it = export_values.begin(); it != export_values.end(); ++it) {
		values.push_back(it->second);
	}
}


/* Needed for adding a BCAnimationCurve into a BCAnimationCurveSet */
inline bool operator< (const BCAnimationCurve& lhs, const BCAnimationCurve& rhs) {
	std::string lhtgt = lhs.get_channel_target();
	std::string rhtgt = rhs.get_channel_target();
	if (lhtgt == rhtgt)
	{
		int lha = lhs.get_array_index();
		int rha = rhs.get_array_index();
		return lha < rha;
	}
	else
		return lhtgt < rhtgt;
}
