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

#define LANPR_POST_PROCESSING_DISABLED 0
#define LANPR_POST_PROCESSING_ENABLED  1

#define LANPR_USE_DIFFERENT_TAPER      0
#define LANPR_USE_SAME_TAPER           1

#define LANPR_DISABLE_TIP_EXTEND       0
#define LANPR_ENABLE_TIP_EXTEND        1

typedef struct LANPRLineStyle{
	struct LANPRLineStyle *next, *prev;

    float thickness;
    float color[4];
    
    int use_camera_distance;
    float camera_distance_influence;
    float camera_distance_exp;

} LANPRLineStyle;




#endif