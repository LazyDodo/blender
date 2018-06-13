/*
 * Copyright 2018, Blender Foundation.
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
 * Contributor(s): Blender Institute
 *
 */

/** \file eevee_lightcache.h
 *  \ingroup eevee
 */

#ifndef __EEVEE_LIGHTCACHE_H__
#define __EEVEE_LIGHTCACHE_H__

#include "BLI_sys_types.h"  /* for bool */

struct ViewLayer;
struct Scene;

void *EEVEE_lightcache_job_data_alloc(struct Main *bmain, struct ViewLayer *viewlayer, struct Scene *scene, bool run_as_job);
void EEVEE_lightcache_job_data_free(void *custom_data);
void EEVEE_lightcache_bake_job(void *custom_data, short *stop, short *do_update, float *progress);

#endif /* __EEVEE_LIGHTCACHE_H__ */