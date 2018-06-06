/*
 * Copyright 2016, Blender Foundation.
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

/** \file draw_hair.c
 *  \ingroup draw
 */

#include "BLI_utildefines.h"

#include "DNA_hair_types.h"
#include "DNA_scene_types.h"

#include "BKE_hair.h"

#include "DRW_render.h"

#include "GPU_extensions.h"
#include "GPU_texture.h"

#include "draw_common.h"

const char* DRW_hair_shader_defines(void)
{
	static char str[256];
	
	BLI_snprintf(str, sizeof(str), "#define HAIR_SHADER_FIBERS\n#define HAIR_SHADER_TEX_WIDTH %d\n",
	             GPU_max_texture_size());
	
	return str;
}

static DRWShadingGroup *drw_shgroup_create_hair_fibers_ex(
        Scene *UNUSED(scene), Object *object, HairSystem *hsys, struct Mesh *scalp,
        DRWPass *hair_pass,
        struct GPUMaterial *gpu_mat, GPUShader *gpu_shader)
{
	/* TODO */
	const int subdiv = 0;
	
	HairExportCache *hair_export = BKE_hair_export_cache_new();
	BKE_hair_export_cache_update(hair_export, hsys, subdiv, scalp, HAIR_EXPORT_ALL);
	
	const DRWHairFiberTextureBuffer *fiber_buffer = NULL;
	struct Gwn_Batch *hair_geom = DRW_cache_hair_get_fibers(hsys, hair_export, &fiber_buffer);
	
	BKE_hair_export_cache_free(hair_export);
	
	DRWShadingGroup *shgrp;
	if (gpu_mat) {
		shgrp = DRW_shgroup_material_create(gpu_mat, hair_pass);
	}
	else if (gpu_shader) {
		shgrp = DRW_shgroup_create(gpu_shader, hair_pass);
	}
	else {
		BLI_assert(0);
	}

	if (!hsys->draw_texture_cache) {
		hsys->draw_texture_cache = DRW_texture_create_2D(fiber_buffer->width, fiber_buffer->height,
		                                                 GPU_RG32F, 0, fiber_buffer->data);
	}
	GPUTexture **fibertex = (GPUTexture **)(&hsys->draw_texture_cache);
	
	DRW_shgroup_uniform_texture_ref(shgrp, "fiber_data", fibertex);
	DRW_shgroup_uniform_int(shgrp, "strand_map_start", &fiber_buffer->strand_map_start, 1);
	DRW_shgroup_uniform_int(shgrp, "strand_vertex_start", &fiber_buffer->strand_vertex_start, 1);
	DRW_shgroup_uniform_int(shgrp, "fiber_start", &fiber_buffer->fiber_start, 1);

	/* TODO(fclem): Until we have a better way to cull the hair and render with orco, bypass culling test. */
	DRW_shgroup_call_object_add_no_cull(shgrp, hair_geom, object);

	return shgrp;
}

DRWShadingGroup *DRW_shgroup_hair_fibers_create(
        Scene *scene, Object *object, HairSystem *hsys, struct Mesh *scalp,
        DRWPass *hair_pass,
        GPUShader *shader)
{
	return drw_shgroup_create_hair_fibers_ex(scene, object, hsys, scalp, hair_pass, NULL, shader);
}

DRWShadingGroup *DRW_shgroup_material_hair_fibers_create(
        Scene *scene, Object *object, HairSystem *hsys, struct Mesh *scalp,
        DRWPass *hair_pass,
        struct GPUMaterial *material)
{
	return drw_shgroup_create_hair_fibers_ex(scene, object, hsys, scalp, hair_pass, material, NULL);
}

void DRW_shgroup_hair(
        Object *ob,
        HairSystem *hsys,
        HairDrawSettings *draw_settings,
        struct Mesh *scalp,
        DRWShadingGroup *shgrp_verts,
        DRWShadingGroup *shgrp_edges)
{
	HairExportCache *hair_export = BKE_hair_export_cache_new();
	BKE_hair_export_cache_update(hair_export, hsys, 0, scalp, HAIR_EXPORT_GUIDE_CURVES | HAIR_EXPORT_GUIDE_VERTICES);
	
	switch (draw_settings->follicle_mode)
	{
		case HAIR_DRAW_FOLLICLE_NONE:
			break;
		case HAIR_DRAW_FOLLICLE_POINTS:
		{
			struct Gwn_Batch *geom = DRW_cache_hair_get_follicle_points(hsys, hair_export);
			DRW_shgroup_call_add(shgrp_verts, geom, ob->obmat);
			break;
		}
	}

	switch (draw_settings->guide_mode)
	{
		case HAIR_DRAW_GUIDE_NONE:
			break;
		case HAIR_DRAW_GUIDE_CURVES:
		{
			struct Gwn_Batch *geom = DRW_cache_hair_get_guide_curve_edges(hsys, hair_export);
			DRW_shgroup_call_add(shgrp_edges, geom, ob->obmat);
			break;
		}
	}
	
	BKE_hair_export_cache_free(hair_export);
}
