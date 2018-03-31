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

#include "BCSampleData.h"
#include "collada_utils.h"

BCSample::BCSample(Object *ob):
	matrix(ob)
{}

BCSample::~BCSample()
{
	int x = 0;
	BCMaterialMap::iterator it;
	for (it = material_map.begin(); it != material_map.end(); ++it) {
		delete it->second;
	}
	material_map.clear();
}

void BCSample::set_bone(Bone *bone, Matrix &mat)
{
	BCMatrix *matrix;
	BCBoneMatrixMap::const_iterator it = bone_matrix_map.find(bone);
	if (it == bone_matrix_map.end()) {
		matrix = new BCMatrix(mat);
		bone_matrix_map[bone] = matrix;
	}
	else {
		/*xxx: not sure if this is ever called (keep for now) */
		matrix = it->second;
		matrix->set_transform(mat);
	}
}

BCMatrix::BCMatrix(Matrix &mat)
{
	set_transform(mat);
}

BCMatrix::BCMatrix(Object *ob)
{
	set_transform(ob);
}

void BCMatrix::set_transform(Matrix &mat)
{
	copy_m4_m4(matrix, mat);
	mat4_decompose(this->loc, this->q, this->size, mat);
	quat_to_eul(this->rot, this->q);
}

void BCMatrix::set_transform(Object *ob)
{
	Matrix lmat;

	BKE_object_matrix_local_get(ob, lmat);
	copy_m4_m4(matrix, lmat);

	mat4_decompose(this->loc, this->q, this->size, lmat);
	quat_to_compatible_eul(this->rot, ob->rot, this->q);
}

const BCMatrix *BCSample::get_matrix(Bone *bone) const
{
	BCBoneMatrixMap::const_iterator it = bone_matrix_map.find(bone);
	if (it == bone_matrix_map.end()) {
		return nullptr;
	}
	return it->second;
}

const BCMatrix *BCSample::get_matrix() const
{
	return &matrix;
}

const BCLamp &BCSample::get_lamp() const
{
	return this->lamp;
}

void BCSample::set_lamp(Lamp *lamp)
{
	this->lamp.light_color[0] = lamp->r;
	this->lamp.light_color[1] = lamp->g;
	this->lamp.light_color[2] = lamp->b;
	this->lamp.falloff_angle = 0;    //TODO from where to get this
	this->lamp.falloff_exponent = 0; // TODO from where to get this
	this->lamp.blender_dist = lamp->dist;
}

const BCCamera &BCSample::get_camera() const
{
	return this->camera;
}

void BCSample::set_camera(Camera *camera)
{
	this->camera.lens = camera->lens;
	this->camera.xsensor = camera->sensor_x;
	this->camera.xfov = RAD2DEGF(focallength_to_fov(camera->lens, camera->sensor_x));
	this->camera.ysensor = camera->sensor_y;
	this->camera.xmag = camera->ortho_scale; // TODO: check this
	this->camera.zfar = camera->clipend;
	this->camera.znear = camera->clipsta;
}

const BCMaterial &BCSample::get_material(int index)
{
	const BCMaterial *mat = this->material_map[index];
	return *mat;
}

void BCSample::set_material(Material *ma)
{
	BCMaterial *material;
	BCMaterialMap::const_iterator it = material_map.find(ma->index);
	if (it == material_map.end()) {
		material = new BCMaterial();
		material_map[ma->index] = material;
	}
	else {
		material = it->second;
	}

	material->diffuse_color[0] = ma->r;
	material->diffuse_color[1] = ma->g;
	material->diffuse_color[2] = ma->b;
	material->specular_color[0] = ma->specr;
	material->specular_color[1] = ma->specg;
	material->specular_color[2] = ma->specb;
	material->alpha = ma->alpha;
	material->ior = ma->refrac;
}

/* Set single float vaules 
 * TODO: replace this by a std::map<rna_path, float>
 * so we can store values more freely
*/
const bool BCSample::set_value(BC_animation_transform_type channel, const int array_index, float val)
{
	switch (channel) {

	/* Light animation */
	case BC_ANIMATION_TYPE_LIGHT_FALLOFF_ANGLE:
		lamp.falloff_angle = val;
		break;
	case BC_ANIMATION_TYPE_LIGHT_FALLOFF_EXPONENT:
		lamp.falloff_exponent = val;
		break;
	case BC_ANIMATION_TYPE_LIGHT_BLENDER_DIST:
		lamp.blender_dist = val;
		break;

	/* Camera animation */
	case BC_ANIMATION_TYPE_LENS:
		camera.lens = val;
		break;
	case BC_ANIMATION_TYPE_XMAG:
		camera.xmag = val;
		break;
	case BC_ANIMATION_TYPE_ZFAR:
		camera.zfar = val;
		break;
	case BC_ANIMATION_TYPE_ZNEAR:
		camera.znear = val;
		break;

	default:
		return false;
	}

	return true;

}

/* Set vector values */
const bool BCSample::set_vector(BC_animation_transform_type channel, float val[3])
{
	float *vp;
	switch (channel) {
	/* Lamp animation */
	case BC_ANIMATION_TYPE_LIGHT_COLOR:
		vp = lamp.light_color;
		break;
	default:
		return false;
	}

	for (int i = 0; i < 3; ++i)
		vp[i] = val[i];

	return true;
}

/* Get channel value */
const bool BCSample::get_value(BC_animation_transform_type channel, const int array_index, float *val, int ma_index) const
{
	BCMaterialMap::const_iterator it = material_map.find(ma_index);
	if (it != material_map.end()) {
		BCMaterial *material = it->second;
		switch (channel) {
			/* Material animation*/
		case BC_ANIMATION_TYPE_SPECULAR_HARDNESS:
			*val = material->specular_hardness;
			break;
		case BC_ANIMATION_TYPE_SPECULAR_COLOR:
			*val = material->specular_color[array_index];
			break;
		case BC_ANIMATION_TYPE_DIFFUSE_COLOR:
			*val = material->diffuse_color[array_index];
			break;
		case BC_ANIMATION_TYPE_ALPHA:
			*val = material->alpha;
			break;
		case BC_ANIMATION_TYPE_IOR:
			*val = material->ior;
		default:
			*val = 0;
			return false;
		}
		return true;
	}
	return false;
}

/* Get channel value */
const bool BCSample::get_value(BC_animation_transform_type channel, const int array_index, float *val) const
{
	switch (channel) {

	/* Object animation */
	case BC_ANIMATION_TYPE_LOCATION:
		*val = matrix.location()[array_index];
		break;
	case BC_ANIMATION_TYPE_SCALE:
		*val = matrix.scale()[array_index];
		break;
	case BC_ANIMATION_TYPE_ROTATION:
	case BC_ANIMATION_TYPE_ROTATION_EULER:
		*val = matrix.rotation()[array_index];
		break;
	case BC_ANIMATION_TYPE_ROTATION_QUAT:
		*val = matrix.quat()[array_index];
		break;


	/* Lamp animation */
	case BC_ANIMATION_TYPE_LIGHT_COLOR:
		*val = lamp.light_color[array_index];
		break;
	case BC_ANIMATION_TYPE_LIGHT_FALLOFF_ANGLE:
		*val = lamp.falloff_angle;
		break;
	case BC_ANIMATION_TYPE_LIGHT_FALLOFF_EXPONENT:
		*val = lamp.falloff_exponent;
		break;
	case BC_ANIMATION_TYPE_LIGHT_BLENDER_DIST:
		*val = lamp.blender_dist;
		break;


	/* Camera animation */
	case BC_ANIMATION_TYPE_LENS:
		*val = camera.lens;
		break;
	case BC_ANIMATION_TYPE_XMAG:
		*val = camera.xmag;
		break;
	case BC_ANIMATION_TYPE_ZFAR:
		*val = camera.zfar;
		break;
	case BC_ANIMATION_TYPE_ZNEAR:
		*val = camera.znear;
		break;

	case BC_ANIMATION_TYPE_UNKNOWN:
	default:
		*val = 0;
		return false;
	}

	return true;
}

void BCMatrix::copy(Matrix &r, Matrix &a)
{
	/* destination comes first: */
	memcpy(r, a, sizeof(Matrix));
}

void BCMatrix::transpose(Matrix &mat)
{
	transpose_m4(mat);
}

void BCMatrix::sanitize(Matrix &mat, int precision)
{
	bc_sanitize_mat(mat, precision);
}

void BCMatrix::unit()
{
	unit_m4(matrix);
}

/*
We need double here because the OpenCollada API needs it.
precision = -1 indicates to not limit the precision
*/
void BCMatrix::get_matrix(double(&mat)[4][4], const bool transposed, const int precision) const
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			float val = (transposed) ? matrix[j][i] : matrix[i][j];
			if (precision >= 0)
				val = floor((val * pow(10, precision) + 0.5)) / pow(10, precision);
			mat[i][j] = val;
		}
}

const bool BCMatrix::in_range(const BCMatrix &other, float distance) const
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (fabs(other.matrix[i][j] - matrix[i][j]) > distance) {
				return false;
			}
		}
	}
	return true;
}

float(&BCMatrix::location() const)[3]
{
	return loc;
}

float(&BCMatrix::rotation() const)[3]
{
	return rot;
}

float(&BCMatrix::scale() const)[3]
{
	return size;
}

float(&BCMatrix::quat() const)[4]
{
	return q;
}


