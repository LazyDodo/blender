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

int lanpr_get_line_bounding_areas(LANPR_RenderBuffer *rb, LANPR_RenderLine *rl, int *rowBegin, int *rowEnd, int *colBegin, int *colEnd);
LANPR_BoundingArea* lanpr_get_point_bounding_area(LANPR_RenderBuffer *rb, real x, real y);

#define LANPR_OTHER_RV(rl,rv) ((rv) == (rl)->L?(rl)->R:(rl)->L) 

LANPR_RenderLine* lanpr_get_connected_render_line(LANPR_BoundingArea* ba, LANPR_RenderVert* rv, LANPR_RenderVert** new_rv) {
	nListItemPointer* lip;
	LANPR_RenderLine* nrl;
	real cosine;

	for (lip = ba->LinkedLines.pFirst; lip; lip = lip->pNext) {
		nrl = lip->p;

		if ((!(nrl->Flags&LANPR_EDGE_FLAG_ALL_TYPE)) || (nrl->Flags & LANPR_EDGE_FLAG_CHAIN_PICKED)) continue;

		// always chain connected lines for now.
		// simplification will take care of the sharp points.
		// if(cosine whatever) continue;

		if (rv != nrl->L && rv != nrl->R) {
			if (nrl->Flags&LANPR_EDGE_FLAG_INTERSECTION) {
				if (rv->FrameBufferCoord[0] == nrl->L->FrameBufferCoord[0] && rv->FrameBufferCoord[1] == nrl->L->FrameBufferCoord[1]) {
					*new_rv = LANPR_OTHER_RV(nrl, nrl->L);
					return nrl;
				}elif(rv->FrameBufferCoord[0] == nrl->R->FrameBufferCoord[0] && rv->FrameBufferCoord[1] == nrl->R->FrameBufferCoord[1]){
					*new_rv = LANPR_OTHER_RV(nrl, nrl->R);
					return nrl;
				}
			}
			continue;
		}

		*new_rv = LANPR_OTHER_RV(nrl, rv);
        return nrl;
    }

    return 0;
}

int lanpr_get_nearby_render_line(LANPR_BoundingArea* ba, LANPR_RenderLine* rl){
    // hold on
	return 1;
}

LANPR_RenderLineChain* lanpr_create_render_line_chain(LANPR_RenderBuffer *rb){
    LANPR_RenderLineChain* rlc;
    rlc = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineChain));

    lstAppendItem(&rb->Chains,rlc);

    return rlc;
}

LANPR_RenderLineChainItem* lanpr_append_render_line_chain_point(LANPR_RenderBuffer *rb, LANPR_RenderLineChain* rlc, float x, float y, char type, int level){
    LANPR_RenderLineChainItem* rlci;
    rlci = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineChainItem));

    rlci->pos[0] = x;
    rlci->pos[1] = y;
    rlci->LineType = type&LANPR_EDGE_FLAG_ALL_TYPE;
	rlci->OcclusionLevel = level;
    lstAppendItem(&rlc->Chain,rlci);

    //printf("a %f,%f %d\n", x, y, level);

    return rlci;
}

LANPR_RenderLineChainItem* lanpr_push_render_line_chain_point(LANPR_RenderBuffer *rb, LANPR_RenderLineChain* rlc, float x, float y, char type, int level){
    LANPR_RenderLineChainItem* rlci;
    rlci = memStaticAquire(&rb->RenderDataPool, sizeof(LANPR_RenderLineChainItem));

    rlci->pos[0] = x;
    rlci->pos[1] = y;
    rlci->LineType = type&LANPR_EDGE_FLAG_ALL_TYPE;
	rlci->OcclusionLevel = level;
    lstPushItem(&rlc->Chain,rlci);

	//printf("p %f,%f %d\n", x, y, level);

    return rlci;
}

// refer to http://karthaus.nl/rdp/ for description
void lanpr_reduce_render_line_chain_recursive(LANPR_RenderLineChain* rlc, LANPR_RenderLineChainItem* from, LANPR_RenderLineChainItem* to, float dist_threshold){
    LANPR_RenderLineChainItem* rlci,*next_rlci;
    float l[2],r[2],c[2];
    float max_dist=0;
    LANPR_RenderLineChainItem* max_rlci=0;

    // find the max distance item
    for(rlci = from->Item.pNext; rlci!= to; rlci=next_rlci){
        next_rlci = rlci->Item.pNext;

        if(next_rlci && (next_rlci->OcclusionLevel!= rlci->OcclusionLevel || next_rlci->LineType!= rlci->LineType)) continue;

		float dist = dist_to_line_segment_v2(rlci->pos, from->pos, to->pos);
		if (dist>dist_threshold && dist>max_dist) { max_dist = dist; max_rlci = rlci; }
		//if (dist <= dist_threshold) lstRemoveItem(&rlc->Chain, (void*)rlci);
    }

	if (!max_rlci) {
		if (from->Item.pNext == to) return;
		for (rlci = from->Item.pNext; rlci != to; rlci = next_rlci) {
			next_rlci = rlci->Item.pNext;
			if (next_rlci && (next_rlci->OcclusionLevel != rlci->OcclusionLevel || next_rlci->LineType != rlci->LineType)) continue;
			lstRemoveItem(&rlc->Chain, (void*)rlci);
		}
	}else {
		if (from->Item.pNext != max_rlci) lanpr_reduce_render_line_chain_recursive(rlc, from, max_rlci, dist_threshold);
		if (to->Item.pPrev != max_rlci)   lanpr_reduce_render_line_chain_recursive(rlc, max_rlci, to, dist_threshold);
	}
}


void lanpr_NO_THREAD_chain_feature_lines(LANPR_RenderBuffer *rb, float dist_threshold){
    LANPR_RenderLineChain* rlc;
	LANPR_RenderLineChainItem* rlci;
    LANPR_RenderLine* rl;
    LANPR_BoundingArea* ba;
    LANPR_RenderLineSegment* rls;

    for(rl = rb->AllRenderLines.pFirst; rl;rl=rl->Item.pNext){

        if((!(rl->Flags&LANPR_EDGE_FLAG_ALL_TYPE)) || (rl->Flags & LANPR_EDGE_FLAG_CHAIN_PICKED)) continue;

        rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;

        rlc = lanpr_create_render_line_chain(rb);

		int r1, r2, c1, c2, row, col;
        LANPR_RenderLine* new_rl = rl;
        LANPR_RenderVert* new_rv;

        // step 1: grow left
        ba = lanpr_get_point_bounding_area(rb,rl->L->FrameBufferCoord[0], rl->L->FrameBufferCoord[1]);
        new_rv = rl->L;
		rls = rl->Segments.pFirst;
        lanpr_push_render_line_chain_point(rb,rlc,new_rv->FrameBufferCoord[0],new_rv->FrameBufferCoord[1],rl->Flags, rls->OcclusionLevel);
        while(ba &&(new_rl = lanpr_get_connected_render_line(ba,new_rv,&new_rv))){
            new_rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;

            int last_occlusion;
            
            if(new_rv==new_rl->L){
                for(rls = new_rl->Segments.pLast; rls;rls=rls->Item.pPrev){
                    float px,py;
                    px = tnsLinearItp(new_rl->L->FrameBufferCoord[0], new_rl->R->FrameBufferCoord[0], rls->at);
                    py = tnsLinearItp(new_rl->L->FrameBufferCoord[1], new_rl->R->FrameBufferCoord[1], rls->at);
                    lanpr_push_render_line_chain_point(rb,rlc,px,py, new_rl->Flags, rls->OcclusionLevel);
                }
            }elif(new_rv==new_rl->R){
                rls = new_rl->Segments.pFirst;
                last_occlusion = rls->OcclusionLevel;
                rls=rls->Item.pNext;
                for(rls; rls; rls= rls->Item.pNext){
                    float px,py;
                    px = tnsLinearItp(new_rl->L->FrameBufferCoord[0], new_rl->R->FrameBufferCoord[0], rls->at);
                    py = tnsLinearItp(new_rl->L->FrameBufferCoord[1], new_rl->R->FrameBufferCoord[1], rls->at);
                    lanpr_push_render_line_chain_point(rb,rlc,px,py, new_rl->Flags,last_occlusion);
                    last_occlusion = rls->OcclusionLevel;
                }
                lanpr_push_render_line_chain_point(rb,rlc, new_rl->R->FrameBufferCoord[0], new_rl->R->FrameBufferCoord[1], new_rl->Flags,last_occlusion);
            }
            ba = lanpr_get_point_bounding_area(rb,new_rv->FrameBufferCoord[0], new_rv->FrameBufferCoord[1]);
        }

        // step 2: this line
		rls = rl->Segments.pFirst;
        for(rls = rls->Item.pNext; rls;rls=rls->Item.pNext){
            float px,py;
            px = tnsLinearItp(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0], rls->at);
            py = tnsLinearItp(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1], rls->at);
            lanpr_append_render_line_chain_point(rb,rlc,px,py,rl->Flags, rls->OcclusionLevel);
        }
		lanpr_append_render_line_chain_point(rb, rlc, rl->R->FrameBufferCoord[0], rl->R->FrameBufferCoord[1], rl->Flags, 0);

        // step 3: grow right
        ba = lanpr_get_point_bounding_area(rb,rl->R->FrameBufferCoord[0], rl->R->FrameBufferCoord[1]);
        new_rv = rl->R;
        // below already done in step 2
        // lanpr_push_render_line_chain_point(rb,rlc,new_rv->FrameBufferCoord[0],new_rv->FrameBufferCoord[1],rl->Flags,0);
        while(ba && (new_rl = lanpr_get_connected_render_line(ba,new_rv,&new_rv))){
            new_rl->Flags |= LANPR_EDGE_FLAG_CHAIN_PICKED;

            int last_occlusion;

			// fix leading vertex type
			rlci = rlc->Chain.pLast;
			rlci->LineType = new_rl->Flags&LANPR_EDGE_FLAG_ALL_TYPE;
            
			if (new_rv == new_rl->L) {
				rls = new_rl->Segments.pLast;
				last_occlusion = rls->OcclusionLevel;
				rlci->OcclusionLevel = last_occlusion;
				rls = rls->Item.pPrev;
				if (rls) last_occlusion = rls->OcclusionLevel;
				for (rls = new_rl->Segments.pLast; rls; rls = rls->Item.pPrev) {
					float px, py;
					px = tnsLinearItp(new_rl->L->FrameBufferCoord[0], new_rl->R->FrameBufferCoord[0], rls->at);
					py = tnsLinearItp(new_rl->L->FrameBufferCoord[1], new_rl->R->FrameBufferCoord[1], rls->at);
					last_occlusion = rls->Item.pPrev ? ((LANPR_RenderLineSegment*)rls->Item.pPrev)->OcclusionLevel : 0;
					lanpr_append_render_line_chain_point(rb, rlc, px, py, new_rl->Flags, last_occlusion);
				}
			}elif(new_rv == new_rl->R) {
				rls = new_rl->Segments.pFirst;
				last_occlusion = rls->OcclusionLevel;
				rlci->OcclusionLevel = last_occlusion;
				rls = rls->Item.pNext;
				for (rls; rls; rls = rls->Item.pNext) {
					float px, py;
					px = tnsLinearItp(new_rl->L->FrameBufferCoord[0], new_rl->R->FrameBufferCoord[0], rls->at);
					py = tnsLinearItp(new_rl->L->FrameBufferCoord[1], new_rl->R->FrameBufferCoord[1], rls->at);
					lanpr_append_render_line_chain_point(rb, rlc, px, py, new_rl->Flags, rls->OcclusionLevel);
					//last_occlusion = rls->OcclusionLevel;
				}
				lanpr_append_render_line_chain_point(rb, rlc, new_rl->R->FrameBufferCoord[0], new_rl->R->FrameBufferCoord[1], new_rl->Flags, 0);
			}
			ba = lanpr_get_point_bounding_area(rb, new_rv->FrameBufferCoord[0], new_rv->FrameBufferCoord[1]);
        }

		//LANPR_RenderLineChainItem* rlci;
		//printf("line:\n");
		//for (rlci = rlc->Chain.pFirst; rlci; rlci = rlci->Item.pNext) {
		//	printf("  %f,%f %d\n", rlci->pos[0],rlci->pos[1], rlci->OcclusionLevel);
		//}
		//printf("--------\n");

        //lanpr_reduce_render_line_chain_recursive(rlc,rlc->Chain.pFirst, rlc->Chain.pLast, dist_threshold);
    }
}

int lanpr_count_chain(LANPR_RenderLineChain* rlc){
    LANPR_RenderLineChainItem* rlci;
    int Count = 0;
	for (rlci = rlc->Chain.pFirst; rlci; rlci = rlci->Item.pNext) {
        Count++;
    }
    return Count;
}

float lanpr_compute_chain_length(LANPR_RenderLineChain* rlc, float* lengths, int begin_index) {
	LANPR_RenderLineChainItem* rlci;
	int i=0;
	float offset_accum = 0;
	float dist;
	float last_point[2];

	rlci = rlc->Chain.pFirst;
	copy_v2_v2(last_point, rlci->pos);
	for (rlci = rlc->Chain.pFirst; rlci; rlci = rlci->Item.pNext) {
		dist = len_v2v2(rlci->pos, last_point);
		offset_accum += dist;
		lengths[begin_index + i] = offset_accum;
		copy_v2_v2(last_point, rlci->pos);
		i++;
	}
	return offset_accum;
}

int lanpr_get_gpu_line_type(LANPR_RenderLineChainItem* rlci) {
	switch (rlci->LineType) {
		case LANPR_EDGE_FLAG_CONTOUR:         return 0;
		case LANPR_EDGE_FLAG_CREASE:          return 1;
		case LANPR_EDGE_FLAG_MATERIAL:        return 2;
		case LANPR_EDGE_FLAG_EDGE_MARK:       return 3;
		case LANPR_EDGE_FLAG_INTERSECTION:    return 4;
		default: return 0;
	}
}

void lanpr_chain_generate_draw_command(LANPR_RenderBuffer *rb){
    LANPR_RenderLineChain* rlc;
    LANPR_RenderLineChainItem* rlci;
    int vert_count=0;
	int i = 0;
	int arg;
    float total_length;
	float* lengths;
	float length_target[2];

    static GPUVertFormat format = { 0 };
	static struct { uint pos, uvs, type, level; } attr_id;
	if (format.attr_len == 0) {
		attr_id.pos = GPU_vertformat_attr_add(&format, "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
        attr_id.uvs = GPU_vertformat_attr_add(&format, "uvs", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
		attr_id.type = GPU_vertformat_attr_add(&format, "type", GPU_COMP_I32, 1, GPU_FETCH_INT);
		attr_id.level = GPU_vertformat_attr_add(&format, "level", GPU_COMP_I32, 1, GPU_FETCH_INT);
	}

	GPUVertBuf *vbo = GPU_vertbuf_create_with_format(&format);

    for(rlc = rb->Chains.pFirst; rlc;rlc=rlc->Item.pNext){
		int count = lanpr_count_chain(rlc);
		//printf("seg contains %d verts\n", count);
		vert_count += count;
    }

    GPU_vertbuf_data_alloc(vbo, vert_count+1); // serve as end point's adj.

	lengths = MEM_callocN(sizeof(float)*vert_count, "chain lengths");

    GPUIndexBufBuilder elb;
	GPU_indexbuf_init_ex(&elb, GPU_PRIM_LINES_ADJ, vert_count*4, vert_count, true);// elem count will not exceed vert_count

    for(rlc = rb->Chains.pFirst; rlc; rlc=rlc->Item.pNext){

		total_length = lanpr_compute_chain_length(rlc, lengths, i);

		for (rlci = rlc->Chain.pFirst; rlci; rlci = rlci->Item.pNext) {

			length_target[0] = lengths[i];
			length_target[1] = total_length - lengths[i];

            GPU_vertbuf_attr_set(vbo, attr_id.pos, i, rlci->pos);
            GPU_vertbuf_attr_set(vbo, attr_id.uvs, i, length_target);

			arg = lanpr_get_gpu_line_type(rlci);
			GPU_vertbuf_attr_set(vbo, attr_id.type, i, &arg);

			arg = (int)rlci->OcclusionLevel;
			GPU_vertbuf_attr_set(vbo, attr_id.level, i, &arg);

			if (rlci == rlc->Chain.pLast) {
				if (rlci->Item.pPrev == rlc->Chain.pFirst) {
					length_target[1] = total_length;
					GPU_vertbuf_attr_set(vbo, attr_id.uvs, i, length_target);
				}
				i++; 
				continue; 
			}

			if (rlci == rlc->Chain.pFirst) {
				if (rlci->Item.pNext == rlc->Chain.pLast) GPU_indexbuf_add_line_adj_verts(&elb, vert_count, i, i + 1, vert_count);
				else GPU_indexbuf_add_line_adj_verts(&elb, vert_count, i, i + 1, i + 2);
			}
			else {
				if (rlci->Item.pNext == rlc->Chain.pLast) GPU_indexbuf_add_line_adj_verts(&elb, i-1, i, i + 1, vert_count);
				else GPU_indexbuf_add_line_adj_verts(&elb, i-1, i, i + 1, i + 2);
			}

			i++;
        }
    }
	//set end point flag value.
	length_target[0] = 3e30f;
	length_target[1] = 3e30f;
	GPU_vertbuf_attr_set(vbo, attr_id.pos, vert_count, length_target);

	MEM_freeN(lengths);

	if (rb->ChainDrawBatch) GPU_batch_discard(rb->ChainDrawBatch);
    rb->ChainDrawBatch = GPU_batch_create_ex(GPU_PRIM_LINES_ADJ, vbo, GPU_indexbuf_build(&elb), GPU_USAGE_DYNAMIC | GPU_BATCH_OWNS_VBO);

}