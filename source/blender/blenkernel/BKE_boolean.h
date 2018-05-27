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
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/BKE_boolean.h
 *  \ingroup blenkernel
 */


#ifndef __BKE_BOOLEAN_H__
#define __BKE_BOOLEAN_H__

struct Object;
struct DerivedMesh;
struct BooleanModifierData;

/* Performs a boolean between two mesh objects, it is assumed that both objects
 * are in fact mesh object. On success returns a DerivedMesh. On failure
 * returns NULL and reports an error. */

struct DerivedMesh *BKE_boolean_bmesh(struct DerivedMesh *dm, struct Object *ob,
                                               struct DerivedMesh *dm_other, struct Object *ob_other, int op_type,
                                               float double_threshold, struct BooleanModifierData *bmd);

#endif  /* BKE_BOOLEAN_H */
