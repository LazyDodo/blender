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
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file ED_hair.h
 *  \ingroup editors
 */

#ifndef __ED_HAIR_H__
#define __ED_HAIR_H__

struct bContext;
struct HairSystem;
struct Object;
struct wmOperator;
struct wmKeyConfig;
struct EditHair;

/* groom_ops.c */
void    ED_operatortypes_hair(void);
void    ED_operatormacros_hair(void);
void    ED_keymap_hair(struct wmKeyConfig *keyconf);

/* editgroom.c */
void    undo_push_hair(struct bContext *C, const char *name);

void ED_hair_edithair_load(struct Object *obedit);
void ED_hair_edithair_make(struct Object *obedit);
void ED_hair_edithair_free(struct Object *obedit);

int ED_hair_object_poll(struct bContext *C);

#endif /* __ED_HAIR_H__ */
