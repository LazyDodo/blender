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
 * Contributor(s): Martin Felke
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/modifiers/intern/MOD_meshcache_util.h
 *  \ingroup modifiers
 */

#ifndef MOD_BOOLEAN_UTIL_BMESH_H
#define MOD_BOOLEAN_UTIL_BMESH_H

struct Mesh *BKE_boolean_operation(struct Mesh *mesh, struct Object *ob,
                                   struct Mesh *mesh_other, struct Object *ob_other, int op_type,
                                   float double_threshold, struct BooleanModifierData *bmd);

#endif // MOD_BOOLEAN_UTIL_BMESH_H
