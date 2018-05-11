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

/** \file eevee_lightcache.c
 *  \ingroup draw_engine
 *
 * Eevee's indirect lighting cache.
 */

#include "DRW_render.h"

#include "BKE_global.h"

#include "DEG_depsgraph_build.h"
#include "DEG_depsgraph_query.h"

#include "DNA_lightprobe_types.h"

#include "PIL_time.h"

#include "eevee_lightcache.h"
#include "eevee_private.h"

typedef struct EEVEE_LightBake {
	Depsgraph *depsgraph;
	ViewLayer *view_layer;
	Scene *scene;
	struct Main *bmain;

	Object *probe;                   /* Current probe being rendered. */
	GPUTexture *rt_color;            /* Target cube color texture. */
	GPUTexture *rt_depth;            /* Target cube depth texture. */
	GPUFrameBuffer *rt_fb;           /* Target cube framebuffer. */

	/* Shared */
	int layer;                       /* Target layer to store the data to. */
	float samples_ct, invsamples_ct; /* Sample count for the convolution. */
	float lod_factor;                /* Sampling bias during convolution step. */
	float lod_max;                   /* Max cubemap LOD to sample when convolving. */
	int cube_count, grid_count;      /*  */

	/* Irradiance grid */
	int irr_cube_res;                /* Target cubemap at MIP 0. */
	int total_irradiance_samples;    /* Total for all grids */
	int bounce_curr, bounce_count;   /* The current light bounce being evaluated. */
	float vis_range, vis_blur;       /* Sample Visibility compression and bluring. */
	float vis_res;                   /* Resolution of the Visibility shadowmap. */
	GPUTexture *tex_double;          /* Result of previous light bounce. */

	/* Reflection probe */
	int ref_cube_res;                /* Target cubemap at MIP 0. */
	float probemat[6][4][4];         /* ViewProjection matrix for each cube face. */
	float texel_size, padding_size;  /* Texel and padding size for the final octahedral map. */
	float roughness;                 /* Roughness level of the current mipmap. */
} EEVEE_LightBake;

typedef struct EEVEE_LightCache {
	int flag;

	/* only a single cache for now */
	int cube_count, grid_count;      /* Number of probes that are ready */
	GPUTexture *grid_tex;
	GPUTexture *cube_tex;
	EEVEE_LightProbe *cube_data;
	EEVEE_LightGrid  *grid_data;
} EEVEE_LightCache;

/* EEVEE_LightCache->flag */
enum {
	LIGHTCACHE_BAKED            = (1 << 0),
	LIGHTCACHE_BAKING           = (1 << 1),
	LIGHTCACHE_CUBE_READY       = (1 << 2),
	LIGHTCACHE_GRID_READY       = (1 << 3),
	/* Update tagging */
	LIGHTCACHE_UPDATE_CUBE      = (1 << 4),
	LIGHTCACHE_UPDATE_GRID      = (1 << 5),
	LIGHTCACHE_UPDATE_WORLD     = (1 << 6),
};

void *EEVEE_lightcache_job_data_alloc(struct Main *bmain, struct ViewLayer *view_layer, struct Scene *scene)
{
	EEVEE_LightBake *lbake = MEM_callocN(sizeof(EEVEE_LightBake), "EEVEE_LightBake");

	lbake->view_layer = view_layer;
	lbake->scene = scene;
	lbake->bmain = bmain;
	lbake->depsgraph = DEG_graph_new(scene, view_layer, DAG_EVAL_RENDER);

    DEG_graph_relations_update(lbake->depsgraph, bmain, scene, view_layer);

	return lbake;
}

void EEVEE_lightcache_job_data_free(void *custom_data)
{
	EEVEE_LightBake *lbake = (EEVEE_LightBake *)custom_data;

	DEG_graph_free(lbake->depsgraph);

	MEM_freeN(lbake);
}

void EEVEE_lightcache_bake_job(void *custom_data, short *UNUSED(stop), short *UNUSED(do_update), float *UNUSED(progress))
{
	EEVEE_LightBake *lbake = (EEVEE_LightBake *)custom_data;
	EEVEE_ViewLayerData *sldata = EEVEE_view_layer_data_ensure_ex(lbake->view_layer);
	EEVEE_LightCache *lcache = EEVEE_lightcache_ensure(sldata);
	Depsgraph *depsgraph = lbake->depsgraph;

	int frame = 0; /* TODO make it user param. */
	DEG_evaluate_on_framechange(lbake->bmain, lbake->depsgraph, frame);

	/* First convert all lightprobes to tight UBO data from all lightprobes in the scene.
	 * This allows a large number of probe to be precomputed. */
	lcache->cube_count = lcache->grid_count = 0;
	lbake->cube_count = lbake->grid_count = 0;

	DEG_OBJECT_ITER_FOR_RENDER_ENGINE_BEGIN(depsgraph, ob, DEG_ITER_OBJECT_MODE_RENDER)
	{
		if (ob->type == OB_LIGHTPROBE) {
			LightProbe *prb = (LightProbe *)ob->data;

			if ((prb->type & LIGHTPROBE_TYPE_GRID) != 0) {
				lbake->grid_count += 1;
				lcache->grid_data = MEM_reallocN(lcache->grid_data, sizeof(EEVEE_LightGrid) * lbake->grid_count);
				// EEVEE_lightprobe_grid_data_from_object(ob, &lbake->grid_data[lbake->grid_count - 1]);
			}
			else if ((prb->type & LIGHTPROBE_TYPE_CUBE) != 0) {
				lbake->cube_count += 1;
				lcache->cube_data = MEM_reallocN(lcache->cube_data, sizeof(EEVEE_LightProbe) * lbake->cube_count);
				// EEVEE_lightprobe_cube_data_from_object(ob, &lbake->grid_data[lbake->cube_count - 1]);
			}
		}
	}
	DEG_OBJECT_ITER_FOR_RENDER_ENGINE_END;

#if 0
	/* Render irradiance */
	lbake->bounce_curr = 0;
	while (lbake->bounce_curr < lbake->bounce_count) {
		for (LinkData *link = ) {
			/* Ask for lower importance draw manager lock. */

			/* TODO: make DRW manager instanciable (and only lock on drawing) */

			/* Create passes */
			/* Iter through objects */
			/* Render one cubemap/irradiance sample. */
			if (*stop != 0) {
				goto cleanup;
			}
		}
		lbake->bounce_curr += 1;
	}

	/* Render reflections */
	for (prb in cube_probes) {
		/* Ask for lower importance draw manager lock. */

		/* Create passes */
		/* Iter through objects */
		/* Render one cubemap/irradiance sample. */
		if (*stop != 0) {
			goto cleanup;
		}
	}
#endif
}

EEVEE_LightCache *EEVEE_lightcache_ensure(EEVEE_ViewLayerData *sldata)
{
	if (sldata->light_cache == NULL) {
		sldata->light_cache = MEM_callocN(sizeof(EEVEE_LightCache), "EEVEE_LightCache");
		sldata->light_cache->cube_data = MEM_mallocN(0, "EEVEE Cube Data Cache");
		sldata->light_cache->grid_data = MEM_mallocN(0, "EEVEE Grid Data Cache");
	}

	return sldata->light_cache;
}

void EEVEE_lightcache_free(EEVEE_LightCache *lcache)
{
	MEM_freeN(lcache->cube_data);
	MEM_freeN(lcache->grid_data);
	MEM_freeN(lcache);
}
