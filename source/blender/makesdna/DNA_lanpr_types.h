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
 * The Original Code is Copyright (C) 2010 Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __DNA_LANPR_TYPES_H__
#define __DNA_LANPR_TYPES_H__

/** \file DNA_npr_types.h
 *  \ingroup DNA
 */

//#include "DNA_listBase.h"
//#include "DNA_ID.h"

struct Object;
struct Material;
struct Collection;

typedef struct LANPR_LineStyleComponent{
    struct LANPR_LineStyle *next, *prev;

    struct Object     object_select;
    struct Material   material_select;
    struct Collection collection_select;

}LANPR_LineStyleComponent;

typedef struct LANPR_LineStyle{
	struct LANPR_LineStyle *next, *prev;

    int      qi_begin;
    int      qi_end;   /* these are for QI Range thing... just occlusion levels */

    float    thickness;
    float    line_thickness_crease;
	float    line_thickness_material;
	float    line_thickness_edge_mark;

    float    color[4];
    float    crease_color[4];
    float    material_color[4];
    float    edge_mark_color[4];

    int      logic_mode; /* for component evaluation */
    ListBase components;

}LANPR_LineStyle;




#endif