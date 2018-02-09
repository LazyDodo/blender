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
 * Contributor(s): Brecht Van Lommel.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __BKE_VOLUME_H__
#define __BKE_VOLUME_H__

/** \file BKE_volume.h
 *  \ingroup bke
 *  \brief General operations for volumes.
 */

struct Main;
struct Volume;

void BKE_volume_init(struct Volume *volume);
void *BKE_volume_add(struct Main *bmain, const char *name);
void BKE_volume_copy_data(struct Main *bmain, struct Volume *volume_dst, const struct Volume *volume_src, const int flag);
struct Volume *BKE_volume_copy(struct Main *bmain, const struct Volume *volume);
void BKE_volume_make_local(struct Main *bmain, struct Volume *volume, const bool lib_local);
void BKE_volume_free(struct Volume *volume);

void BKE_volume_reload(struct Main *bmain, struct Volume *volume);

#endif
