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

#ifndef __BC_SAMPLE_H__
#define __BC_SAMPLE_H__

#include <string>
#include <map>
#include <algorithm>

extern "C"
{
#include "DNA_object_types.h"
#include "DNA_armature_types.h"
#include "DNA_material_types.h"
}

/*
 * The list of currently supported animation types
 * TODO: Maybe this can be made more general.
*/
typedef enum BC_animation_transform_type {

	/* Translation channels */
	BC_ANIMATION_TYPE_ROTATION_EULER = 0,
	BC_ANIMATION_TYPE_ROTATION_QUAT = 1,
	BC_ANIMATION_TYPE_SCALE = 2,
	BC_ANIMATION_TYPE_LOCATION = 3,

	/* Material channels */
	BC_ANIMATION_TYPE_SPECULAR_HARDNESS = 4,
	BC_ANIMATION_TYPE_SPECULAR_COLOR = 5,
	BC_ANIMATION_TYPE_DIFFUSE_COLOR = 6,
	BC_ANIMATION_TYPE_ALPHA = 7,
	BC_ANIMATION_TYPE_IOR = 8,

	/* Lamp channels */
	BC_ANIMATION_TYPE_LIGHT_COLOR,
	BC_ANIMATION_TYPE_FALL_OFF_ANGLE,
	BC_ANIMATION_TYPE_FALL_OFF_EXPONENT,
	BC_ANIMATION_TYPE_BLENDER_DIST,

	/* Camera channels */
	BC_ANIMATION_TYPE_XFOV,
	BC_ANIMATION_TYPE_XMAG,
	BC_ANIMATION_TYPE_ZFAR,
	BC_ANIMATION_TYPE_ZNEAR,

	/* other */
	BC_ANIMATION_TYPE_ROTATION,
	BC_ANIMATION_TYPE_UNKNOWN = -1,

} BC_animation_transform_type;


class BCMaterial {
public:
	float specular_hardness;
	float specular_color[3];
	float diffuse_color[3];
	float alpha;
	float ior;
};

class BCLamp {
public:
	float light_color[3];
	float falloff_angle;
	float falloff_exponent;
	float blender_dist;
};

class BCCamera {
public:
	float xfov;
	float xmag;
	float zfar;
	float znear;
};


/*
* BCMatrix is a neat helper class to handle float matrix[4][4] items
* in a convenient way.
* Actually i could not figure out how to make BCMatrix an equivalent
* of a float matrix[4][4] That is why i added the Matrix typedef
* below. 
*
* As the code stands now, i must do this for example:
*
* float mat[4][4];
* BCMatrix obmat;
* copy_m4_m4(obmat.matrix, mat);
*
* I would prefer if i had only BCMatrix and do:
*
* float mat[4][4];
* BCMatrix obmat;
* copy_m4_m4(obmat, mat);
*
* Can this be achieved ?
*/

typedef float(Matrix)[4][4];
class BCMatrix {

	mutable float size[3];
	mutable float rot[3];
	mutable float loc[3];
	mutable float q[4];
	mutable bool decomposed = false;

	/* Private methods */
	void decompose() const;
	void unit();
	void copy(float(&r)[4][4], float(&a)[4][4]);

public:
	mutable float matrix[4][4];

	float(&location() const)[3];
	float(&rotation() const)[3];
	float(&scale() const)[3];
	float(&quat() const)[4];

	BCMatrix()
	{
		unit();
	}
	void get_matrix(double(&mat)[4][4], const bool transposed = false, const int precision = -1) const;
	const bool in_range(const BCMatrix &other, float distance) const;
};

typedef std::map<int, BCMaterial *> BCMaterialMap;
typedef std::map<Bone *, BCMatrix *> BCBoneMatrixMap;

class BCSample{
private:
	
	/* For Object Transformations */
	BCMatrix matrix;

	/* XXX: The following parts are exclusive, each BCSample has 
	   at most one of those filled with data. 
	   Maybe we can make a union here?
	*/
	BCMaterialMap material_map; /* For Material animation */
	BCBoneMatrixMap bone_matrix_map; /* For Armature animation */
	BCLamp lamp; /* For Lamp channels */
	BCCamera camera; /* For Camera channels */


public:
	BCSample(double(&mat)[4][4]);
	BCSample(float(&mat)[4][4]);
	~BCSample();

	void set_matrix(double(&mat)[4][4]);
	void set_matrix(float(&mat)[4][4]);
	void set_matrix(BCSample &other);
	void set_material(Material *ma);
	void set_bone(Bone *bone, Matrix &mat);
	const BCMatrix *get_sampled_matrix() const;
	const BCMatrix *get_sampled_matrix(Bone *bone) const;

	const bool set_vector(BC_animation_transform_type channel, float val[3]);
	const bool set_value(BC_animation_transform_type channel, const int array_index, float val);

	const float(&get_matrix() const)[4][4];
	const bool get_value(BC_animation_transform_type channel, const int array_index, float *val) const;
	const bool get_value(int ma_index, BC_animation_transform_type channel, const int array_index, float *val) const;

	/* Convenient helper functions */
	static void sanitize(float(&matrix)[4][4], int precision);
	static void transpose(float(&matrix)[4][4]);
	bool in_range(const BCSample &other, const float distance) const;

};

#endif