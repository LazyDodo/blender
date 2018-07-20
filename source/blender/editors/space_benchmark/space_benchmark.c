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
 * Contributor(s): Campbell Barton
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/space_benchmark/space_benchmark.c
 *  \ingroup spbenchmark
 */

#include <string.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_utildefines.h"

#include "BKE_context.h"
#include "BKE_screen.h"

#include "ED_space_api.h"
#include "ED_screen.h"

#include "WM_api.h"
#include "WM_types.h"

/* ******************** default callbacks for benchmark space ***************** */

static SpaceLink *benchmark_new(const ScrArea *UNUSED(area), const Scene *UNUSED(scene))
{
	SpaceBenchmark *sbenchmark = MEM_callocN(sizeof(SpaceBenchmark), "initbenchmark");
	sbenchmark->spacetype = SPACE_BENCHMARK;

	// TODO sbenchmark->lheight =  14;

	/* main region */
	ARegion *ar = MEM_callocN(sizeof(ARegion), "main region for text");
	BLI_addtail(&sbenchmark->regionbase, ar);
	ar->regiontype = RGN_TYPE_WINDOW;

	return (SpaceLink *)sbenchmark;
}

/* not spacelink itself */
static void benchmark_free(SpaceLink *UNUSED(sl))
{
}

/* spacetype; init callback */
static void benchmark_init(struct wmWindowManager *UNUSED(wm), ScrArea *UNUSED(sa))
{
}

static SpaceLink *benchmark_duplicate(SpaceLink *sl)
{
	SpaceBenchmark *sbenchmarkn = MEM_dupallocN(sl);
	return (SpaceLink *)sbenchmarkn;
}

/* add handlers, stuff you only do once or on area/region changes */
static void benchmark_main_region_init(wmWindowManager *wm, ARegion *ar)
{
	ED_region_panels_init(wm, ar);
}

static void benchmark_main_region_draw(const bContext *C, ARegion *ar)
{
	ED_region_panels(C, ar);
}

static void benchmark_main_region_listener(
        wmWindow *UNUSED(win), ScrArea *UNUSED(sa), ARegion *UNUSED(ar),
        wmNotifier *UNUSED(wmn), const Scene *UNUSED(scene))
{
}

static void benchmark_keymap(struct wmKeyConfig *UNUSED(keyconf))
{
	// wmKeyMap *keymap = WM_keymap_find(keyconf, "Benchmark", SPACE_BENCHMARK, 0);
}

/* only called once, from space/spacetypes.c */
void ED_spacetype_benchmark(void)
{
	SpaceType *st = MEM_callocN(sizeof(SpaceType), "spacetype benchmark");
	ARegionType *art;

	st->spaceid = SPACE_BENCHMARK;
	strncpy(st->name, "Benchmark", BKE_ST_MAXNAME);

	st->new = benchmark_new;
	st->free = benchmark_free;
	st->init = benchmark_init;
	st->duplicate = benchmark_duplicate;
	st->keymap = benchmark_keymap;

	/* regions: main window */
	art = MEM_callocN(sizeof(ARegionType), "spacetype benchmark region");
	art->regionid = RGN_TYPE_WINDOW;
	art->keymapflag = ED_KEYMAP_UI | ED_KEYMAP_VIEW2D;

	art->init = benchmark_main_region_init;
	art->draw = benchmark_main_region_draw;
	art->listener = benchmark_main_region_listener;

	BLI_addhead(&st->regiontypes, art);

	BKE_spacetype_register(st);
}
