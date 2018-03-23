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

#ifndef __BC_MATRIX_H__
#define __BC_MATRIX_H__

class BCSample {
private:
public:
};

class BCMatrix: public BCSample {
private:
	mutable bool decomposed = false;
	float matrix[4][4];
	mutable float size[3];
	mutable float rot[3];
	mutable float loc[3];
	mutable float q[4];

	void decompose() const;
	void unit();
	void copy(float(&r)[4][4], float(&a)[4][4]);

public:
	BCMatrix();
	~BCMatrix();
	BCMatrix(double(&mat)[4][4]);
	BCMatrix(float(&mat)[4][4]);
	void set_matrix(double(&mat)[4][4]);
	void set_matrix(float(&mat)[4][4]);
	void set_matrix(BCMatrix &other);
	void get_matrix(double(&mat)[4][4], const bool transposed = false, const int precision = -1) const;
	void get_matrix(float(&mat)[4][4]) const;
	static void sanitize(float(&mat)[4][4], int precision);
	static void transpose(float(&mat)[4][4]);
	bool in_range(BCMatrix other, float distance) const;

	const float(&location() const)[3];
	const float(&rotation() const)[3];
	const float(&scale() const)[3];
	const float(&quat() const)[4];

};

#endif