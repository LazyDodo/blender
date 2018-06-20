#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
#include "lanpr_all.h"
#include "DRW_render.h"
#include "BKE_object.h"
#include "DNA_mesh_types.h"
#include "DNA_camera_types.h"
#include "GPU_immediate.h"
#include "GPU_immediate_util.h"
#include "GPU_framebuffer.h"
#include "DNA_lanpr_types.h"
#include "DNA_meshdata_types.h"
#include "BKE_customdata.h"
#include "DEG_depsgraph_query.h"
#include "GPU_draw.h"

#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"
#include "bmesh.h"

#include "WM_types.h"
#include "WM_api.h"

#include <math.h>

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/

struct Object;

void lanpr_generate_geom_buffer(struct Object *ob){
    
}

static int lanpr_compute_feature_lines_exec(struct bContext *C, struct wmOperator *op){

	return OPERATOR_FINISHED;
}

static void lanpr_compute_feature_lines_cancel(struct bContext *C, struct wmOperator *op){

    return;
}


void SCENE_OT_lanpr_calculate_feature_lines(struct wmOperatorType* ot){

	/* identifiers */
	ot->name = "Calculate Feature Lines";
	ot->description = "LANPR calculates feature line in current scene";
	ot->idname = "SCENE_OT_lanpr_calculate";

	/* api callbacks */
	//ot->invoke = screen_render_invoke; /* why we need both invoke and exec? */
	//ot->modal = screen_render_modal;
	ot->cancel = lanpr_compute_feature_lines_cancel;
	ot->exec = lanpr_compute_feature_lines_exec;
}