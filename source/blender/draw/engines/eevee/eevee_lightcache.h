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
struct EEVEE_LightCache;

/* Light Bake */
void *EEVEE_lightbake_job_data_alloc(struct Main *bmain, struct ViewLayer *viewlayer, struct Scene *scene, bool run_as_job);
void EEVEE_lightbake_job_data_free(void *custom_data);
void EEVEE_lightbake_update(void *custom_data);
void EEVEE_lightbake_job(void *custom_data, short *stop, short *do_update, float *progress);

/* Light Cache */
struct EEVEE_LightCache *EEVEE_lightcache_create(const int cube_count, const int cube_size, const int irr_size[3]);
void EEVEE_lightcache_free(struct EEVEE_LightCache *lcache);

#endif /* __EEVEE_LIGHTCACHE_H__ */