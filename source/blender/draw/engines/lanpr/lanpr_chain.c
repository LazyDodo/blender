#include "DRW_engine.h"
#include "DRW_render.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
#include "BLI_math.h"
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

        if(nrl->Flags & LANPR_EDGE_FLAG_CHAIN_PICKED) continue;

        // always chain connected lines for now.
        // if(cosine whatever) continue;

        if(rv == nrl->L) *NextV = nrl->R;
        elif(rv == nrl->R) *NextV = nrl->L;
        else continue;

        return nrl;
    }

    return 0;
}

int lanpr_GetNearByRenderLine(LANPR_BoundingArea* ba, LANPR_RenderLine* rl){
    // hold on
}

LANPR_RenderLineChain* lanpr_create_render_line_chain(LANPR_RenderBuffer *rb){
    LANPR_RenderLineChain* rlc;
    rlc = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineChain));

    lstAppendItem(&rb->Chains,rlc);

    return rlc;
}

LANPR_RenderLineChainItem* lanpr_append_render_line_chain_point(LANPR_RenderBuffer *rb, LANPR_RenderLineChain* rlc, LANPR_RenderVert*rv, char type){
    LANPR_RenderLineChainItem* rlci;
    rlci = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineChainItem));

    rlci->rv = rv;
    rlci->LineType = type&LANPR_EDGE_FLAG_ALL_TYPE;
    lstAppendItem(&rlc->Chain,rlci);

    return rlci;
}

LANPR_RenderLineChainItem* lanpr_push_render_line_chain_point(LANPR_RenderBuffer *rb, LANPR_RenderLineChain* rlc, LANPR_RenderVert*rv, char type){
    LANPR_RenderLineChainItem* rlci;
    rlci = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineChainItem));

    rlci->rv = rv;
    rlci->LineType = type&LANPR_EDGE_FLAG_ALL_TYPE;
    lstPushItem(&rlc->Chain,rlci);

    return rlci;
}


void lanpr_reduce_render_line_chain_recursive(LANPR_RenderLineChain* rlc, LANPR_RenderLineChainItem* from, LANPR_RenderLineChainItem* to, float dist_threshold){
    LANPR_RenderLineChainItem* rlci,*next_rlci;
    float l[2],r[2],c[2];
    float max_dist=0;
    LANPR_RenderLineChainItem* max_rlci=0;

    copy_v2fl_v2db(l,from->rv->FrameBufferCoord);
    copy_v2fl_v2db(r,to->rv->FrameBufferCoord);

    // find the max distance item
    for(rlci = from; rlci!= to->Item.pNext; rlci=next_rlci){
        next_rlci = rlci->Item.pNext;

        copy_v2fl_v2db(c,rlci->rv->FrameBufferCoord);
        float dist = dist_to_line_segment_v2(c,l,r);

        if(dist>dist_threshold && dist>max_dist){ max_dist = dist; max_rlci=rlci; continue;}

        if(dist<=dist_threshold) lstRemoveItem(&rlc->Chain, rlci);
    }

    lanpr_reduce_render_line_chain_recursive(rlc, from, max_rlci, dist_threshold);
    lanpr_reduce_render_line_chain_recursive(rlc, max_rlci, to, dist_threshold);
}


#define LANPR_OTHER_RV(rl,rv)\
((rv) == (rl)->L?(rl)->R:(rl)->L) 

void lanpr_ChainFeatureLines_NO_THREAD(LANPR_RenderBuffer *rb, float dist_threshold){
    LANPR_RenderLineChain* rlc;
    LANPR_RenderLine* rl;

    for(rl = rb->AllRenderLines.pFirst; rl;rl=rl->Item.pNext){

        if(rl->Flags & LANPR_EDGE_FLAG_CHAIN_PICKED) continue;

        rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;

        rlc = lanpr_create_render_line_chain(rb);

		int r1, r2, c1, c2, row, col;
        LANPR_RenderLine* new_rl = rl;
        LANPR_RenderVert* new_rv;
		if (lanpr_GetLineBoundingAreas(rb, rl, &r1, &r2, &c1, &c2)) {
            for (row = r1; row != r2 + 1; row++) {
                for (col = c1, col != c2 + 1; col++) {

                    //grow left side
                    new_rv = rl->L;
                    lanpr_push_render_line_chain_point(rb,rlc,new_rv,rl->Flags);
                    while(new_rl = lanpr_GetConnectedRenderLine(&rb->InitialBoundingAreas[row * 20 + col],&new_rv)){
                        new_rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;
                        new_rv = LANPR_OTHER_RV(new_rl,new_rv);
                        
                    }

                    //grow right side
                    new_rv = rl->R;
                    lanpr_append_render_line_chain_point(rb,rlc,new_rv,rl->Flags);
                    while(new_rl = lanpr_GetConnectedRenderLine(&rb->InitialBoundingAreas[row * 20 + col],&new_rv)){
                        new_rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;
                        new_rv = LANPR_OTHER_RV(new_rl,new_rv);
                    }

                    //simplification
                    lanpr_reduce_render_line_chain_recursive(rlc, rlc->Chain.pFirst; rlc->Chain.pLast, dist_threshold);
                }
            }
        }
    }
}

int lanpr_CountChainVertices(LANPR_RenderLineChain* rlc){
    LANPR_RenderLineChainItem* rlci;
    int Count = 0;
    for(rlci = rlc->Chain.pFirst;rlci = rlci->Item.pNext){
        Count++;
    }
    return Count;
}

void lanpr_ChainGenerateDrawCommand(LANPR_RenderBuffer *rb){
    LANPR_RenderLineChain* rlc;
    LANPR_RenderLineChainItem* rlci;
    int vert_count;
    int i=0;
    float point[2];
    float last_point[2];
    float offset_accum=0;

    static Gwn_VertFormat format = { 0 };
	static struct { uint pos, offset, type, level; } attr_id;
	if (format.attr_len == 0) {
		attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
        attr_id.offset = GWN_vertformat_attr_add(&format, "offset", GWN_COMP_F32, 1, GWN_FETCH_FLOAT);
		attr_id.type = GWN_vertformat_attr_add(&format, "type", GWN_COMP_I32, 1, GWN_FETCH_FLOAT);
		attr_id.level = GWN_vertformat_attr_add(&format, "level", GWN_COMP_I32, 1, GWN_FETCH_INT);
	}

	Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);

    for(rlc = rb->Chains.pFirst; rlc;rlc=rlc->Item.pNext){
        vert_count += lanpr_CountChainVertices(rlc);
    }

    GWN_vertbuf_data_alloc(vbo, vert_count);

    for(rlc = rb->Chains.pFirst; rlc;rlc=rlc->Item.pNext){
        for(rlci = rlc->Chain.pFirst;rlci = rlci->Item.pNext){
            copy_v2fl_v2db(point,rlci->rv->FrameBufferCoord);
            
            GWN_vertbuf_attr_set(vbo, attr_id.pos, i, c);
            GWN_vertbuf_attr_set(vbo, attr_id.offset, i, &offset_accum);

            offset_accum += len_v2v2(point,last_point);
            copy_v2fl_v2db(last_point,point);

            i++;
        }
    }

    rb->ChainDrawBatch = GWN_batch_create_ex(GWN_PRIM_POINTS, vbo, 0, GWN_USAGE_DYNAMIC | GWN_BATCH_OWNS_VBO);

}