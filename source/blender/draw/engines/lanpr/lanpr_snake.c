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
#include "GPU_draw.h"

#include "GPU_batch.h"
#include "GPU_framebuffer.h"
#include "GPU_shader.h"
#include "GPU_uniformbuffer.h"
#include "GPU_viewport.h"
#include "bmesh.h"

extern struct LANPROneTimeInit OneTime;

int _TNS_ColOffsets[] = { -1,0,1,1,1,0,-1,-1 };
int _TNS_RowOffsets[] = { -1,-1,-1,0,1,1,1,0 };

int _TNS_Deviates[8][8] = {
	{ 0,1,2,3,4,3,2,1 },
	{ 1,0,1,2,3,4,3,2 },
	{ 2,1,0,1,2,3,4,3 },
	{ 3,2,1,0,1,2,3,4 },
	{ 4,3,2,1,0,1,2,3 },
	{ 3,4,3,2,1,0,1,2 },
	{ 2,3,4,3,2,1,0,1 },
	{ 1,2,3,4,3,2,1,0 }
};

#define TNS_CLAMP_TEXTURE_W(t,Col)\
	{if (Col >= t->width) Col = t->width - 1; if (Col < 0) Col = 0;}

#define TNS_CLAMP_TEXTURE_H(t,Row)\
	{if (Row >= t->height) Row = t->height - 1;if (Row < 0) Row = 0;}

#define TNS_CLAMP_TEXTURE_CONTINUE(t,Col,Row)\
	{if (Col >= t->width) continue; if (Col < 0) continue;\
     if (Row >= t->height) continue; if (Row < 0)continue; }


static LANPR_TextureSample* lanpr_any_uncovered_samples(LANPR_PrivateData* pd){
	return BLI_pophead(&pd->pending_samples);
}

int lanpr_direction_deviate(int From, int To) {
	return _TNS_Deviates[From - 1][To - 1];
}

int lanpr_detect_direction(LANPR_PrivateData* pd, int Col, int Row, int LastDirection) {
	int Deviate[9] = {100};
	int MinDeviate = 0;
	int i;
	LANPR_TextureSample* ts;

	for (i = 0; i < 8; i++) {
		TNS_CLAMP_TEXTURE_CONTINUE(pd, (_TNS_ColOffsets[i] + Col), (_TNS_RowOffsets[i] + Row));
		if (ts = pd->sample_table[(_TNS_ColOffsets[i] + Col) + (_TNS_RowOffsets[i] + Row) * pd->width]) {
			if (!LastDirection) return i + 1;
			Deviate[i+1] = lanpr_direction_deviate(i, LastDirection);
			if (!MinDeviate || Deviate[MinDeviate] > Deviate[i + 1]) MinDeviate = i + 1;
		}
	}

	return MinDeviate;
}

LANPR_LineStrip* lanpr_create_line_strip(LANPR_PrivateData* pd) {
	LANPR_LineStrip* ls = BLI_mempool_calloc(pd->mp_line_strip);
	return ls;
}
LANPR_LineStripPoint* lanpr_append_point(LANPR_PrivateData* pd, LANPR_LineStrip* ls, real X, real Y, real Z) {
	LANPR_LineStripPoint* lsp = BLI_mempool_calloc(pd->mp_line_strip_point);

	lsp->P[0] = X;
	lsp->P[1] = Y;
	lsp->P[2] = Z;

	BLI_addtail(&ls->points,lsp);

	ls->point_count++;

	return lsp;
}
LANPR_LineStripPoint* lanpr_push_point(LANPR_PrivateData* pd, LANPR_LineStrip* ls, real X, real Y, real Z) {
	LANPR_LineStripPoint* lsp = BLI_mempool_calloc(pd->mp_line_strip_point);

	lsp->P[0] = X;
	lsp->P[1] = Y;
	lsp->P[2] = Z;

	BLI_addhead(&ls->points, lsp);

	ls->point_count++;

	return lsp;
}

//void lanpr_remove_point(LANPR_LineStrip* ls,tnsLineStripPoint* lsp) {
//	lstRemoveItem(&ls->Points, lsp);
//	FreeMem(lsp);
//}

void lanpr_destroy_line_strip(LANPR_PrivateData* pd, LANPR_LineStrip* ls) {
	LANPR_LineStripPoint* lsp;
	while (lsp = BLI_pophead(&ls->points)) {
		BLI_mempool_free(pd->mp_line_strip_point, lsp);
	}
	BLI_mempool_free(pd->mp_line_strip, ls);
}

void lanpr_remove_sample(LANPR_PrivateData* pd, int Row, int Col) {
	LANPR_TextureSample* ts;
	ts = pd->sample_table[Row*pd->width + Col];
	pd->sample_table[Row*pd->width + Col] = 0;
	
	BLI_remlink(&pd->pending_samples, ts);
	ts->Item.prev=0;ts->Item.next=0;
	BLI_addtail(&pd->erased_samples, ts);
}

void lanpr_grow_snake_r(LANPR_PrivateData* pd, LANPR_LineStrip* ls, LANPR_LineStripPoint* ThisP, int Direction) {
	LANPR_LineStripPoint* NewP = ThisP,*p2;
	int Length = 5;
	int l = 0;
	int Deviate,Dir=Direction,NewDir;
	int AddPoint;
	int TX = NewP->P[0], TY = NewP->P[1];

	while (NewDir = lanpr_detect_direction(pd, TX, TY, Dir)) {
		AddPoint = 0;
		Deviate = lanpr_direction_deviate(NewDir, Dir);
		Dir = NewDir;

		l++;
		TX += _TNS_ColOffsets[NewDir-1];
		TY += _TNS_RowOffsets[NewDir-1];

		if (Deviate < 2) {
			lanpr_remove_sample(pd, TY, TX);
		}elif(Deviate < 3) {
			lanpr_remove_sample(pd, TY, TX);
			AddPoint = 1;
		}else {
			lanpr_remove_sample(pd, TY, TX);
			return;
		}

		if (AddPoint || l == Length) {
			p2 = lanpr_append_point(pd, ls, TX, TY, 0);
			NewP = p2;
			l = 0;
		}
	}
	if(TX!=ThisP->P[0] || TY!=ThisP->P[1])
		lanpr_append_point(pd, ls, TX, TY, 0);
}

void lanpr_grow_snake_l(LANPR_PrivateData* pd, LANPR_LineStrip* ls, LANPR_LineStripPoint* ThisP, int Direction) {
	LANPR_LineStripPoint* NewP = ThisP, *p2;
	int Length = 5;
	int l = 0;
	int Deviate, Dir = Direction, NewDir;
	int AddPoint;
	int TX = NewP->P[0], TY = NewP->P[1];

	while (NewDir = lanpr_detect_direction(pd, TX, TY, Dir)) {
		AddPoint = 0;
		Deviate = lanpr_direction_deviate(NewDir, Dir);
		Dir = NewDir;

		l++;
		TX += _TNS_ColOffsets[NewDir - 1];
		TY += _TNS_RowOffsets[NewDir - 1];

		if (Deviate < 2) {
			lanpr_remove_sample(pd, TY, TX);
		}elif(Deviate < 4) {
			lanpr_remove_sample(pd, TY, TX);
			AddPoint = 1;
		}
		else {
			lanpr_remove_sample(pd, TY, TX);
			return;
		}

		if (AddPoint || l == Length) {
			p2 = lanpr_push_point(pd, ls, TX, TY, 0);
			NewP = p2;
			l = 0;
		}
	}
	if(TX!=ThisP->P[0] || TY!=ThisP->P[1])
		lanpr_push_point(pd, ls, TX, TY, 0);
}

int lanpr_reverse_direction(int From) {
	From -= 4;
	if (From <= 0)From += 8;
	return From;
}



void lanpr_texture_to_ndc(int x,int y, int w,int h, float* x_ndc, float* y_ndc){
    *x_ndc = tnsLinearItp(-1,1,(float)x/(float)w);
	*y_ndc = tnsLinearItp(-1,1, (float)y/(float)h);
}

void lanpr_count_drawing_elements(LANPR_PrivateData* pd, int* vert_count, int* index_adjacent_count){
    int v_count=0;
	int e_count=0;
	LANPR_LineStrip* ls;
	for (ls = (LANPR_LineStrip *)(pd->line_strips.first); ls; ls = (LANPR_LineStrip *)(ls->Item.next)) {
		v_count += (ls->point_count);
		e_count += ((ls->point_count - 1) * 4);
	}
	*vert_count = v_count;
	*index_adjacent_count = e_count;
}

Gwn_Batch *lanpr_get_snake_batch(LANPR_PrivateData* pd){
	LANPR_LineStrip* ls;
	LANPR_LineStripPoint* lsp, *plsp;
	int i;
	//u32bit *Index_adjacent;
	float* Verts;
	float* Lengths;
	float TotalLength=0;
	int v_count,e_count;

	lanpr_count_drawing_elements(pd,&v_count,&e_count);

	//Index_adjacent = MEM_callocN(sizeof(unsigned int) * e_count, "Index_adjacent buffer pre alloc");
	Verts = MEM_callocN(sizeof(float) * v_count * 2, "Verts buffer pre alloc");
	Lengths = MEM_callocN(sizeof(float)* v_count * 2, "Length buffer pre alloc");

	Gwn_IndexBufBuilder elb;
	GWN_indexbuf_init_ex(&elb, GWN_PRIM_LINES_ADJ, e_count, v_count, true);

	int vert_offset=0;

	for (ls = (LANPR_LineStrip *)(pd->line_strips.first); ls; ls = (LANPR_LineStrip *)(ls->Item.next)) {
		for (i = 0; i < ls->point_count-1; i++) {
			int v1 = i+vert_offset-1;
			int v2 = i+vert_offset;
			int v3 = i+vert_offset+1;
			int v4 = i+vert_offset+2;
			if(v1<0) v1=0;
			if(v4>= v_count) v4 = v_count -1;
			GWN_indexbuf_add_line_adj_verts(&elb, v1,v2,v3,v4);
		}

		i = 0;
		float xf,yf;
		TotalLength=0;
		for (lsp = (LANPR_LineStripPoint *)(ls->points.first); lsp; lsp = (LANPR_LineStripPoint *)(lsp->Item.next)) {
			lanpr_texture_to_ndc(lsp->P[0],lsp->P[1],pd->width, pd->height, &xf,&yf);
			Verts[vert_offset*2 + i * 2 + 0] = xf;
			Verts[vert_offset*2 + i * 2 + 1] = yf;
			if (plsp = (LANPR_LineStripPoint *)(lsp->Item.prev)) {
				TotalLength += tMatDist2v(plsp->P, lsp->P);
				Lengths[(vert_offset + i) * 2] = TotalLength;
			}
			i++;
		}

		ls->total_length = TotalLength;
        i = 0;
		for (lsp = (LANPR_LineStripPoint *)(ls->points.first); lsp; lsp = (LANPR_LineStripPoint *)(lsp->Item.next)) {
			if (plsp = (LANPR_LineStripPoint *)(lsp->Item.prev)) {
				Lengths[(vert_offset + i) * 2 + 1] = ls->total_length - Lengths[(vert_offset + i) * 2];
			}
			i++;
		}

		vert_offset+=(ls->point_count);
	}

	static Gwn_VertFormat format = { 0 };
	static struct { uint pos, uvs; } attr_id;
	if (format.attrib_ct == 0) {
		attr_id.pos = GWN_vertformat_attr_add(&format, "pos", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
		attr_id.uvs = GWN_vertformat_attr_add(&format, "uvs", GWN_COMP_F32, 2, GWN_FETCH_FLOAT);
	}

	Gwn_VertBuf *vbo = GWN_vertbuf_create_with_format(&format);
	GWN_vertbuf_data_alloc(vbo, v_count);

	for (int i = 0; i < v_count; ++i) {
		GWN_vertbuf_attr_set(vbo, attr_id.pos, i, &Verts[i*2]);
		GWN_vertbuf_attr_set(vbo, attr_id.uvs, i, &Lengths[i*2]);
	}

	MEM_freeN(Verts);
	MEM_freeN(Lengths);
	
	return GWN_batch_create_ex(GWN_PRIM_LINES_ADJ, vbo, GWN_indexbuf_build(&elb), GWN_USAGE_STREAM);
}

void lanpr_snake_draw_scene(LANPR_TextureList* txl, LANPR_FramebufferList * fbl, LANPR_PassList *psl, LANPR_PrivateData *pd, SceneLANPR *lanpr){
    GPUFrameBufferBits clear_bits = GPU_COLOR_BIT;
    float clear_col[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	float clear_depth = 1.0f;
	uint clear_stencil = 0xFF;
    
    DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
	DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();

    const DRWContextState *draw_ctx = DRW_context_state_get();
	View3D *v3d = draw_ctx->v3d;
	RegionView3D *rv3d = draw_ctx->rv3d;
	Object *camera = (rv3d->persp == RV3D_CAMOB) ? v3d->camera : NULL;
    
    pd->znear = camera? ((Camera*)camera->data)->clipsta:0.1;
    pd->zfar = camera? ((Camera*)camera->data)->clipend:100;
    pd->normal_clamp =    lanpr->normal_clamp;
    pd->normal_strength = lanpr->normal_strength;
    pd->depth_clamp =     lanpr->depth_clamp;
    pd->depth_strength =  lanpr->depth_strength;

    GPU_framebuffer_bind(fbl->edge_intermediate);
    DRW_draw_pass(psl->edge_intermediate);

    if((!lanpr->enable_vector_trace) && (!lanpr->display_thinning_result)){
        GPU_framebuffer_bind(dfbl->default_fb);
        DRW_transform_to_display(txl->edge_intermediate);
        return;
    }

    if(lanpr->display_thinning_result || lanpr->enable_vector_trace){
        pd->stage = 0;
        GPU_framebuffer_bind(fbl->edge_thinning);
        clear_bits = GPU_COLOR_BIT;
        GPU_framebuffer_clear(fbl->edge_thinning, clear_bits, clear_col, clear_depth, clear_stencil);
        DRW_draw_pass(psl->edge_thinning);

        pd->stage = 1;
        GPU_framebuffer_bind(fbl->edge_intermediate);
        //GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
        DRW_draw_pass(psl->edge_thinning_2);

        pd->stage = 0;
        GPU_framebuffer_bind(fbl->edge_thinning);
        GPU_framebuffer_clear(fbl->edge_thinning, clear_bits, clear_col, clear_depth, clear_stencil);
        DRW_draw_pass(psl->edge_thinning);

        pd->stage = 1;
        GPU_framebuffer_bind(fbl->edge_intermediate);
        //GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, clear_col, clear_depth, clear_stencil);
        DRW_draw_pass(psl->edge_thinning_2);

        if(!lanpr->enable_vector_trace){
            GPU_framebuffer_bind(dfbl->default_fb);
            DRW_transform_to_display(txl->edge_intermediate);
            return;
        }
    }
    
    int texw = GPU_texture_width(txl->edge_intermediate) ,texh = GPU_texture_height(txl->edge_intermediate);;
	int tsize = texw*texh;
    int recreate=0;
    if(tsize != pd->width*pd->height) recreate =1;

    if(recreate){
        if(pd->line_result) MEM_freeN(pd->line_result);
        pd->line_result = MEM_callocN(sizeof(float) * tsize,"Texture readback buffer");

        if(pd->line_result_8bit) MEM_freeN(pd->line_result_8bit);
        pd->line_result_8bit = MEM_callocN(sizeof(unsigned char) * tsize,"Texture readback buffer 8bit");

        if(pd->sample_table) MEM_freeN(pd->sample_table);
        pd->sample_table = MEM_callocN(sizeof(void*) * tsize,"Texture readback buffer 8bit");

        pd->width = texw;
        pd->height = texh;
    }

    GPU_framebuffer_read_color(fbl->edge_intermediate,0,0,texw, texh,1,0, pd->line_result);

    float sample;
    int h, w;
    for (h = 0; h < texh; h++) {
        for (w = 0; w < texw; w++) {
            int index = h*texw + w;
            if ((sample = pd->line_result[index]) > 0.9) {
                pd->line_result_8bit[index] = 255;
                LANPR_TextureSample* ts = BLI_mempool_calloc(pd->mp_sample);
                BLI_addtail(&pd->pending_samples, ts);
                pd->sample_table[index] = ts;
                ts->X = w;
                ts->Y = h;
            }else{
                pd->sample_table[index] = 0;
            }
        }
    }

    LANPR_TextureSample *ts;
    LANPR_LineStrip* ls;
    LANPR_LineStripPoint* lsp;
    while(ts = lanpr_any_uncovered_samples(pd)){
        int Direction=0;
        LANPR_LineStripPoint tlsp = { 0 };

        tlsp.P[0] = ts->X;
        tlsp.P[1] = ts->Y;

        if (Direction = lanpr_detect_direction(pd, ts->X,ts->Y, Direction)) {
            BLI_addtail(&pd->line_strips, (ls = lanpr_create_line_strip(pd)));
            lsp = lanpr_append_point(pd, ls, ts->X, ts->Y, 0);
            lanpr_remove_sample(pd, ts->Y, ts->X);

            lanpr_grow_snake_r(pd, ls, lsp, Direction);

            lanpr_grow_snake_l(pd, ls, lsp, lanpr_reverse_direction(Direction));
        }

        //count++;
    }

    //GPU_framebuffer_bind()
    GPU_framebuffer_clear(fbl->edge_intermediate, clear_bits, lanpr->background_color, clear_depth, clear_stencil);

    float* tld = &lanpr->taper_left_distance, *tls = &lanpr->taper_left_strength,
        *trd = &lanpr->taper_right_distance, *trs = &lanpr->taper_right_strength;

    Gwn_Batch* snake_batch = lanpr_get_snake_batch(pd);

    psl->snake_pass = DRW_pass_create("Snake Visualization Pass", DRW_STATE_WRITE_COLOR);
    pd->snake_shgrp = DRW_shgroup_create(OneTime.snake_connection_shader, psl->snake_pass);
    DRW_shgroup_uniform_float(pd->snake_shgrp, "LineWidth", &lanpr->line_thickness, 1);
    DRW_shgroup_uniform_float(pd->snake_shgrp, "TaperLDist", tld, 1);
    DRW_shgroup_uniform_float(pd->snake_shgrp, "TaperLStrength", tls, 1);
    DRW_shgroup_uniform_float(pd->snake_shgrp, "TaperRDist", lanpr->use_same_taper?tld:trd, 1);
    DRW_shgroup_uniform_float(pd->snake_shgrp, "TaperRStrength", lanpr->use_same_taper?tls:trs, 1);
    DRW_shgroup_uniform_vec4(pd->snake_shgrp, "LineColor", lanpr->line_color, 1);

    DRW_shgroup_call_add(pd->snake_shgrp, snake_batch, NULL);
    DRW_draw_pass(psl->snake_pass);
    

    BLI_mempool_clear(pd->mp_sample);
    BLI_mempool_clear(pd->mp_line_strip);
    BLI_mempool_clear(pd->mp_line_strip_point);

    pd->pending_samples.first = pd->pending_samples.last = 0;
    pd->erased_samples.first = pd->erased_samples.last = 0;
    pd->line_strips.first = pd->line_strips.last = 0;
    //pd->line_strip_point.first = pd->line_strip_point.last = 0;


    GPU_framebuffer_bind(dfbl->default_fb);

    DRW_transform_to_display(txl->edge_intermediate);   
}