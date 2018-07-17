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

#include <math.h>

int lanpr_GetLineBoundingAreas(LANPR_RenderBuffer *rb, LANPR_RenderLine *rl, int *RowBegin, int *RowEnd, int *ColBegin, int *ColEnd) ;

real lanpr_ThreePointCosine(LANPR_RenderVert* l,LANPR_RenderVert* c, LANPR_RenderVert* r){

}

LANPR_RenderLine* lanpr_GetConnectedRenderLine(LANPR_BoundingArea* ba, LANPR_RenderVert* rv, LANPR_RenderVert** NextV){
    nListItemPointer* lip;
    LANPR_RenderLine* nrl;
    real cosine;

    for(lip = ba->LinkedLines.pFirst; lip; lip=lip->Item.pNext){
        nrl = lip->p;
        if(rv == nrl->L) *NextV = nrl->R;
        elif(rv == nrl->R) *NextV = nrl->L;
        else continue;

        // always chain connected lines for now.
        // if(cosine whatever) continue;

        return nrl;
    }

    return 0;
}

int lanpr_GetNearByRenderLine(LANPR_BoundingArea* ba, LANPR_RenderLine* rl){
    
}

void lanpr_ChainFeatureLines_NO_THREAD(LANPR_RenderBuffer *rb){
    LANPR_RenderLineChain* rlc;
    LANPR_RenderLine* rl;

    for(rl = rb->AllRenderLines.pFirst; rl;rl=rl->Item.pNext){

        if(rl->Flags & LANPR_EDGE_FLAG_CHAIN_PICKED) continue;

		int r1, r2, c1, c2, row, col;
        LANPR_RenderLine* new_rl = rl;
		if (lanpr_GetLineBoundingAreas(rb, rl, &r1, &r2, &c1, &c2)) {
            for (row = r1; row != r2 + 1; row++) {
                for (col = c1, col != c2 + 1; col++) {
                    while(new_rl = lanpr_GetConnectedRenderLine(&rb->InitialBoundingAreas[row * 20 + col],new_rl)){
                        new_rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;

                    }
                    while(new_rl = lanpr_GetNearByRenderLine(&rb->InitialBoundingAreas[row * 20 + col],new_rl)){
                        new_rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;
                    }
                }
            }
        }
    }


}