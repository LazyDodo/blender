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

BCMatrix::BCMatrix()
{
	unit();
}

BCMatrix::BCMatrix(double(&mat)[4][4])
{
	set_matrix(mat);
}

BCMatrix::BCMatrix(float(&mat)[4][4])
{
	set_matrix(mat);
}

BCMatrix::~BCMatrix()
{
	int x = 0;
}

const bool BCMatrix::get_value_for(const std::string &target, const int array_index, float *val) const
{
	if (target == "location") {
		const float(&floc)[3] = location();
		*val = floc[array_index];
	}
	else if (target == "scale") {
		const float(&fsize)[3] = scale();
		*val = fsize[array_index];
	}
	else if (
		target == "rotation" ||
		target == "rotation_euler") {
		const float(&frot)[3] = rotation();
		*val = frot[array_index];
	}
	else if (
		target == "rotation_quaternion") {
		const float(&qt)[4] = quat();
		*val = qt[array_index];
	}
	else
	{
		*val = 0;
		return false;
	}
	return true;
}

void BCMatrix::copy(float(&r)[4][4], float(&a)[4][4])
{
	/* destination comes first: */
	memcpy(r, a, sizeof(float[4][4]));
}

void BCMatrix::transpose(float(&mat)[4][4])
{
	transpose_m4(mat);
}

void BCMatrix::sanitize(float(&matrix)[4][4], int precision)
{
	bc_sanitize_mat(matrix, precision);
}

void BCMatrix::unit()
{
	unit_m4(matrix);
}

void BCMatrix::set_matrix(double(&mat)[4][4])
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrix[i][j] = mat[i][j];
}

void BCMatrix::set_matrix(float(&mat)[4][4])
{
	copy_m4_m4(matrix, mat);
}

void BCMatrix::set_matrix(BCMatrix &other)
{
	set_matrix(other.matrix);
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

void BCMatrix::get_matrix(float(&mat)[4][4]) const
{
	// copy_m4_m4(mat, matrix);  // does not work because copy_m4_m4 does not declare 2nd parameter as const

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat[i][j] = matrix[i][j];
}

bool BCMatrix::in_range(const BCMatrix &other, float distance) const
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

void BCMatrix::decompose() const
{
	float mat[4][4];
	get_matrix(mat);
	mat4_decompose(loc, q, size, mat);
	quat_to_eul(rot, q);
	decomposed = true;
}

const float(&BCMatrix::location() const)[3]
{
	if (!decomposed)
	decompose();

return loc;
}

const float(&BCMatrix::rotation() const)[3]
{
	if (!decomposed)
	decompose();

return rot;

}

const float(&BCMatrix::scale() const)[3]
{
	if (!decomposed)
	decompose();

return size;
}

const float(&BCMatrix::quat() const)[4]
{
	if (!decomposed)
	decompose();

return q;
}

