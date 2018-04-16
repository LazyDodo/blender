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
 * The Original Code is Copyright (C) 2018 Blender Foundation.
 * All rights reserved.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/interface/interface_init_exit.h
 *  \ingroup edinterface
 */

#include "BLI_rect.h"

#include "DNA_screen_types.h"
#include "DNA_userdef_types.h"

#include "BKE_screen.h"

#include "MEM_guardedalloc.h"

#include "RNA_access.h"

#include "UI_interface.h"

#include "WM_api.h"
#include "WM_types.h"

#include "interface_intern.h"



void ui_init_button_group_types(void)
{
	WM_uibuttongrouptype_init();
	WM_uibuttongrouptype_add(UI_BGT_sortable_id_tabs());
}

void ui_exit_button_group_types(void)
{
	WM_uibuttongrouptype_free();
}
