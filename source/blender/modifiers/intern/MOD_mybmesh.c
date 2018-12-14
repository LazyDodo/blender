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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_mybmesh.c
 *  \ingroup modifiers
 */

/* This code is based of the tessellation part of the paper
 * "Computing Smooth Surface Contours with Accurate Topology"
 * (Pierre Bénard, Aaron Hertzmann, Michael Kass).
 * Currently available at:
 * http://www.labri.fr/perso/pbenard/publications/contours.html
 *
 * The numbers in the comments refers to the chapters in the paper.
 */


#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_string.h"
#include "BLI_utildefines.h"
#include "BLI_buffer.h"
#include "BLI_alloca.h"
#include "BLI_ghash.h"
#include "BLI_gsqueue.h"
#include "BLI_memarena.h"
#include "BLI_rand.h"
#include "BLI_listbase.h"
#include "BLI_threads.h"

#include "BKE_library_query.h"
#include "BKE_modifier.h"
#include "BKE_mesh.h"
#include "BKE_mesh_mapping.h"
#include "BKE_deform.h"

#include "bmesh.h"
#include "bmesh_tools.h"

#include "MOD_util.h"
#include "MOD_modifiertypes.h"

#include "DEG_depsgraph_build.h"

//BMesh intern
#include "../../bmesh/intern/bmesh_private.h"

//TODO this modifier depends on OpenSubDiv. So if it's not compiled in, remove this modifier

#include "../../../../intern/opensubdiv/opensubdiv_capi.h"
#include "../../../../intern/opensubdiv/opensubdiv_converter_capi.h"
#include "../../../../intern/opensubdiv/opensubdiv_evaluator_capi.h"
#include "../../../../intern/opensubdiv/opensubdiv_topology_refiner_capi.h"

#include "../../blenkernel/intern/subdiv_converter.h"

#include "PIL_time.h"
#include "PIL_time_utildefines.h"

struct OpenSubdiv_Evaluator;

typedef struct {
	BMVert *vert; //Only used for shifted verts
	BMEdge *orig_edge;
	BMFace *orig_face;
	float u, v;
} Vert_buf;

typedef struct {
	BMEdge *cusp_e;
	BMFace *orig_face;
	float cusp_co[3];
	float cusp_no[3];

	float u, v;
} Cusp;

typedef struct {
	bool b_arr[3];
    bool kr_arr[3];
	float co_arr[3][3];
	float u_arr[3];
	float v_arr[3];
} Cusp_triang;

typedef struct {
	BMVert *vert;
	BMVert *C_vert;
	//Can we extend this radial vert?
	bool extendable;
	bool is_B;
	float radi_plane_no[3];
} Radi_vert;

typedef struct {
	BMFace *face;
	//Should be front or back facing?
	bool back_f;
} IncoFace;

typedef struct {
	BMesh *bm;
	BMesh *bm_orig;

	float cam_loc[3];

    GHash *vert_hash;

	BLI_Buffer *new_vert_buffer;
	BLI_Buffer *shifted_verts;
	BLI_Buffer *cusp_edges;
	BLI_Buffer *C_verts;
	BLI_Buffer *cusp_verts;
	BLI_Buffer *radi_vert_buffer;
	//Radial edge vert start idx
	int new_vert_idx;
	int radi_start_idx;

	struct OpenSubdiv_Evaluator *eval;
} MeshData;

//TODO for Kr look in subdiv.cpp in coutours source code (II)

//TODO dynamic arrays, use BLI_stack, BLI_buffer, BLI_mempool, BLI_memarena.

static void verts_to_limit(BMesh *bm, struct OpenSubdiv_Evaluator *eval){

	int i, j;

	BMIter iter_v, iter_f;
	BMVert *vert;
	BMFace *f;

	//TODO is it possible to only get non adjusted verts?
	//IE not moving a vert more than once.

	BM_ITER_MESH_INDEX (f, &iter_f, bm, BM_FACES_OF_MESH, i) {
			BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, j) {
				float u,v, du[3], dv[3];
				switch(j){
					case 1 :
						u = 1, v = 0;
						break;
					case 2 :
						u = v = 1;
						break;
					case 3 :
						u = 0, v = 1;
						break;
					default:
						u = v = 0;
						break;
				}
				eval->evaluateLimit(eval, i, u, v, vert->co, du, dv);
				//Adjust vert normal to the limit normal
				cross_v3_v3v3(vert->no, du, dv);
				normalize_v3(vert->no);
				//printf("j: %d\n",j);
			}
			//printf("i: %d\n",i);
			//printf("face i: %d\n", BM_elem_index_get(f));
	}

}

static Vert_buf* get_shift_vert( BMVert *vert, MeshData *m_d ){
	int vert_i;

	//check if vert is in the buffer
	for(vert_i = 0; vert_i < m_d->shifted_verts->count; vert_i++){
		Vert_buf *buf = &BLI_buffer_at(m_d->shifted_verts, Vert_buf, vert_i);
		if( vert == buf->vert ){
			return buf;
		}
	}
	return NULL;
}

static void add_shifted_vert( BMVert *vert, BMFace *orig_face, float uv[2], MeshData *m_d ){
	Vert_buf *buf = get_shift_vert( vert, m_d );

	if(buf != NULL){
		buf->orig_face = orig_face;
		buf->u = uv[0];
		buf->v = uv[1];
		return;
	}

	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int v_idx = BM_elem_index_get(vert);
	if( (v_idx + 1) > orig_verts){
		//If this is a non original vert, just update the data instead
		Vert_buf *v_buf = &BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);
		v_buf->orig_face = orig_face;
		v_buf->orig_edge = NULL;
		v_buf->u = uv[0];
		v_buf->v = uv[1];
	} else {
		Vert_buf new_buf;
		new_buf.orig_face = orig_face;
		new_buf.vert = vert;
		new_buf.u = uv[0];
		new_buf.v = uv[1];
		BLI_buffer_append(m_d->shifted_verts, Vert_buf, new_buf);
	}
}

static bool calc_if_B(const float cam_loc[3], const float P[3], const float du[3], const float dv[3]){
	//Is the point back facing?
	float nor[3], view_vec[3];

	cross_v3_v3v3(nor, du, dv);
	//TODO normalization is probably not needed
	normalize_v3(nor);
	sub_v3_v3v3(view_vec, cam_loc, P);

	return ( dot_v3v3(nor, view_vec) < 0);
}

static bool calc_if_B_nor(const float cam_loc[3], const float P[3], const float nor[3]){
	//Is the point back facing?
	float view_vec[3];

	sub_v3_v3v3(view_vec, cam_loc, P);

	return ( dot_v3v3(nor, view_vec) < 0);
}


static float get_facing_dir(const float cam_loc[3], const float P[3], const float du[3], const float dv[3]){
	//Get if point is back facing (-) or front facing (+). Zero if it's directly on the contour
	float nor[3], view_vec[3];

	cross_v3_v3v3(nor, du, dv);
	normalize_v3(nor);
	sub_v3_v3v3(view_vec, cam_loc, P);
	normalize_v3(view_vec);

	return dot_v3v3(nor, view_vec);
}


static float get_facing_dir_nor(const float cam_loc[3], const float P[3], const float nor[3]){
	//Get if point is back facing (-) or front facing (+). Zero if it's directly on the contour
	float view_vec[3];

	sub_v3_v3v3(view_vec, cam_loc, P);
	normalize_v3(view_vec);

	return dot_v3v3(nor, view_vec);
}


static BMVert* split_edge_and_move_nor(MeshData *m_d, BMesh *bm, BMEdge *edge, const float new_pos[3], const float new_no[3]){
	//Split edge one time and move the created vert to new_pos
	BMVert *vert;
    BMFace *face_arr[2];

	BMIter iter;
	BMFace *face;
	int i;

	//Save the connected faces for triangulation later
	BM_ITER_ELEM_INDEX (face, &iter, edge, BM_FACES_OF_EDGE, i){
		face_arr[i] = face;
	}
	//printf("Split edge!\n");

	//TODO perhaps use BM_edge_split instead?
    vert = bmesh_kernel_split_edge_make_vert(bm, edge->v1, edge, NULL);
	if( m_d != NULL ){
		BM_elem_index_set(vert, m_d->new_vert_idx++);
	}

	/*{
		MemArena *pf_arena;
		pf_arena = BLI_memarena_new(BLI_POLYFILL_ARENA_SIZE, __func__);
		LinkNode *faces_double = NULL;
		BM_face_triangulate(
				bm, face_arr[0],
				NULL, NULL,
				NULL, NULL,
				&faces_double,
				MOD_TRIANGULATE_QUAD_FIXED,
				0,false,
				pf_arena,
				NULL, NULL);

		BM_face_triangulate(
				bm, face_arr[1],
				NULL, NULL,
				NULL, NULL,
				&faces_double,
				MOD_TRIANGULATE_QUAD_FIXED,
				0,false,
				pf_arena,
				NULL, NULL);

		while (faces_double) {
			LinkNode *next = faces_double->next;
			BM_face_kill(bm, faces_double->link);
			MEM_freeN(faces_double);
			faces_double = next;
		}
		BLI_memarena_free(pf_arena);
	}*/

	//Triangulate the faces connected to the new vert
	{
		BMLoop *start;
		BMLoop *end;

		BMFace *new_f;

		for( int j = 0; j < i; j++ ){
			start = BM_face_vert_share_loop( face_arr[j], vert );
			end = (start->next)->next;

			//TODO maybe use bmesh_sfme instead?
			new_f = BM_face_split(bm, face_arr[j], start, end, NULL, NULL, false);

			BM_face_normal_update(face_arr[j]);
			BM_face_normal_update(new_f);
		}
	}

	copy_v3_v3(vert->co, new_pos);
	//Adjust vert normal to the limit normal
	copy_v3_v3(vert->no, new_no);
	normalize_v3(vert->no);

	return vert;
}

static BMVert *split_edge_and_move_vert(MeshData *m_d, BMesh *bm, BMEdge *edge, const float new_pos[3],
									const float du[3], const float dv[3]){

	float new_no[3];

	cross_v3_v3v3(new_no, du, dv);
	normalize_v3(new_no);
	return split_edge_and_move_nor(m_d, bm, edge, new_pos, new_no);
}

static bool get_uv_coord(BMVert *vert, BMFace *f, float *u, float *v){
	//Get U,V coords of a vertex
	int i;
	BMIter iter;
	BMVert *vert_iter;

	BM_ITER_ELEM_INDEX (vert_iter, &iter, f, BM_VERTS_OF_FACE, i) {
		if(vert == vert_iter){
			switch(i){
				case 1 :
					*u = 1, *v = 0;
					break;
				case 2 :
					*u = *v = 1;
					break;
				case 3 :
					*u = 0, *v = 1;
					break;
				default:
					*u = *v = 0;
					break;
			}
			return true;
		}
	}

	return false;
}

static bool is_vert_in_buffer(BMVert *v, BLI_Buffer *vert_buf){
	int vert_j;
	for(vert_j = 0; vert_j < vert_buf->count; vert_j++){
		BMVert *c_vert = BLI_buffer_at(vert_buf, BMVert*, vert_j);
		if( c_vert == v ){
			return true;
		}
	}
	return false;
}

static bool point_inside_v2(const float mat[3][3], const float point[2], BMFace *f){
	//TODO maybe add a sanity check to see if the face is not a quad or a triangle
	float (*mat_coords)[2] = BLI_array_alloca(mat_coords, f->len);
	BMVert *vert;
	BMIter iter_v;
	int vert_idx;

	BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, vert_idx) {
		mul_v2_m3v3(mat_coords[vert_idx], mat, vert->co);
	}

	if( f->len == 3 ){
		return isect_point_tri_v2(point, mat_coords[0], mat_coords[1], mat_coords[2]);
	}
	return isect_point_quad_v2(point, mat_coords[0], mat_coords[1], mat_coords[2], mat_coords[3]);
}

static bool point_inside(const float mat[3][3], const float point[3], BMFace *f){
	float mat_new_pos[2];
	mul_v2_m3v3(mat_new_pos, mat, point);

	return point_inside_v2(mat, mat_new_pos, f);
}

static void get_uv_point(BMFace *face, float uv[2], const float point_v2[2], const float mat[3][3] ){
	int vert_idx;
	float st[4][2];

	BMVert *v;
	BMIter iter_v;

	BM_ITER_ELEM_INDEX (v, &iter_v, face, BM_VERTS_OF_FACE, vert_idx) {
		mul_v2_m3v3(st[vert_idx], mat, v->co);
	}

	resolve_quad_uv_v2(uv, point_v2, st[0], st[1], st[2], st[3]);

    if( uv[0] > 1.0f ){
		uv[0] = 1.0f;
	} else if( uv[0] < 0.0f ){
		uv[0] = 0.0f;
	}

    if( uv[1] > 1.0f ){
		uv[1] = 1.0f;
	} else if( uv[1] < 0.0f ){
		uv[1] = 0.0f;
	}

}

static void get_st_point(BMFace *face, const float u, const float v, const float mat[3][3], float st[2]){
	int vert_idx;
	float face_st[4][3];
	float temp_st[3];

	BMVert *vert;
	BMIter iter_v;

	BM_ITER_ELEM_INDEX (vert, &iter_v, face, BM_VERTS_OF_FACE, vert_idx) {
		zero_v3(face_st[vert_idx]);
		mul_v2_m3v3(face_st[vert_idx], mat, vert->co);
	}

	interp_bilinear_quad_v3(face_st, u, v, temp_st);
	copy_v2_v2(st, temp_st);
}

static void get_vert_st_pos(MeshData *m_d, const float mat[3][3], BMVert *vert, float st[2]){
	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	Vert_buf* shift_vert = get_shift_vert( vert, m_d );
	if( shift_vert != NULL ){
		//This vert has been shifted
		get_st_point(shift_vert->orig_face, shift_vert->u, shift_vert->v, mat, st);
	} else {
		//Check if edge verts doesn't belong to orig_face
		int v_idx = BM_elem_index_get(vert);
		if( (v_idx + 1) > orig_verts){
			Vert_buf v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);

			get_st_point(v_buf.orig_face, v_buf.u, v_buf.v, mat, st);
		} else {
			BMVert *temp_v = BLI_ghash_lookup(m_d->vert_hash, vert);

			mul_v2_m3v3(st, mat, temp_v->co);
		}
	}
}

typedef struct FFBB_thread_data {
	//Shared data
	MeshData *m_d;
	BMEdge **edges;
	int orig_edges;
	float *step_arr;

	int cur_e;
	int tot_e;
	SpinLock spin;
} FFBB_thread_data;

static int FFBB_queue_next_e(FFBB_thread_data *queue)
{
	int edge_idx = -1;

	BLI_spin_lock(&queue->spin);
	if (queue->cur_e < queue->tot_e) {
		edge_idx = queue->cur_e;
		queue->cur_e++;
	}
	BLI_spin_unlock(&queue->spin);

	return edge_idx;
}

static void split_BB_FF_edges_thread(void *data_v) {
	//Split BB,FF edges if they have sign crossings
    FFBB_thread_data *th_data = (FFBB_thread_data *)data_v;
	MeshData *m_d = th_data->m_d;

	int i, face_index;
	BMIter iter_f;
	BMEdge *e;
	BMFace *f;
	BMVert *v1, *v2;
	float v1_u, v1_v, v2_u, v2_v;
	float v1_face_dir, v2_face_dir;
	bool is_B;

	while ((i = FFBB_queue_next_e(th_data)) >= 0) {
		Vert_buf v_buf;
		e = th_data->edges[i];

		v1_face_dir = get_facing_dir_nor(m_d->cam_loc, e->v1->co, e->v1->no);
		v2_face_dir = get_facing_dir_nor(m_d->cam_loc, e->v2->co, e->v2->no);

		if( (v1_face_dir < 0) != (v2_face_dir < 0) ){
			//This is not a FF or BB edge
			continue;
		}

		//Is this edge FF or BB?
		is_B = (v1_face_dir < 0);

		//Check if any of the verts facing are near flipping
		//It there is not, then it is highly unlikely that there is a zero crossing between them
		if( fabs(v1_face_dir) > 0.2f && fabs(v2_face_dir) > 0.2f ){
        	continue;
		}

		if( i < th_data->orig_edges ){
			//This edge exists on the original mesh
			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
			//why is table_ensure not fixing the assert?
			BMEdge *orig_e = BM_edge_at_index_find(m_d->bm_orig, i);

			//Get face connected to edge from orig mesh
			//TODO is it wise to use BM_ITER_ELEM here?
			BM_ITER_ELEM (f, &iter_f, orig_e, BM_FACES_OF_EDGE) {
				//Get first face
				break;
			}

			face_index = BM_elem_index_get(f);

			v1 = orig_e->v1;
			v2 = orig_e->v2;

			v_buf.orig_edge = orig_e;
			v_buf.orig_face = f;
		} else {
			BMVert *vert_arr[2];

			//This should be safe because the vert count is still the same as the original mesh.
			v1 = BLI_ghash_lookup(m_d->vert_hash, e->v1);
			v2 = BLI_ghash_lookup(m_d->vert_hash, e->v2);

			vert_arr[0] = v1;
			vert_arr[1] = v2;

			//TODO add checks if to hande if there is no face
			BLI_spin_lock(&th_data->spin);
			f = BM_face_exists_overlap(vert_arr, 2);
			face_index = BM_elem_index_get(f);
			BLI_spin_unlock(&th_data->spin);

			v_buf.orig_edge = NULL;
			v_buf.orig_face = f;
		}

		//TODO can I just check each vert once?
		get_uv_coord(v1, f, &v1_u, &v1_v);
		get_uv_coord(v2, f, &v2_u, &v2_v);
		{
			int j;
			float u, v;
			float P[3], du[3], dv[3];

			if( v1_u == v2_u ){
				u = v1_u;

				for(j=0; j < 10; j++){
					v = th_data->step_arr[j];
					m_d->eval->evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						BLI_spin_lock(&th_data->spin);
						split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						BLI_spin_unlock(&th_data->spin);
						break;
					}
				}
			} else if ( v1_v == v2_v ){
				v = v1_v;

				for(j=0; j < 10; j++){
					u = th_data->step_arr[j];
					m_d->eval->evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						BLI_spin_lock(&th_data->spin);
						split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						BLI_spin_unlock(&th_data->spin);
						break;
					}
				}
			} else {
				bool alt_diag;
				if((v1_u == 0 && v1_v == 0) || (v2_u == 0 && v2_v == 0)){
					alt_diag = false;
				} else {
					alt_diag = true;
				}
				for(j=0; j < 10; j++){
					if(alt_diag){
						u = 1.0f - th_data->step_arr[j];
					} else {
						u = th_data->step_arr[j];
					}

					v = th_data->step_arr[j];
					m_d->eval->evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						BLI_spin_lock(&th_data->spin);
						split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						BLI_spin_unlock(&th_data->spin);
						break;
					}
				}
			}
		}

	}

}

static void split_BB_FF_edges_thread_start(MeshData *m_d){
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);

	void    *edge_stack[BM_DEFAULT_ITER_STACK_SIZE];
	int      edge_len;
	BMEdge **edges = BM_iter_as_arrayN(m_d->bm, BM_EDGES_OF_MESH, NULL, &edge_len,
	                                   edge_stack, BM_DEFAULT_ITER_STACK_SIZE);

	//Do 10 samples but don't check end and start point
	float step = 1.0f/11.0f;
	float step_arr[] = { step*5.0f, step*6.0f, step*4.0f, step*7.0f, step*3.0f,
		step*8.0f, step*2.0f, step*9.0f, step*1.0f, step*10.0f };

	FFBB_thread_data th_data;
	th_data.m_d = m_d;
	th_data.edges = edges;
	th_data.orig_edges = orig_edges;
	th_data.step_arr = step_arr;

	th_data.cur_e = 0;
	th_data.tot_e = edge_len;

	ListBase threads;
	int tot_thread = BLI_system_thread_count();

	BLI_threadpool_init(&threads, split_BB_FF_edges_thread, tot_thread);

	BLI_spin_init(&th_data.spin);

	/* fill in threads handles */
	for (int i = 0; i < tot_thread; i++) {
		BLI_threadpool_insert(&threads, &th_data);
	}

	BLI_threadpool_end(&threads);

	BLI_spin_end(&th_data.spin);

	MEM_freeN(edges);
}

static void split_BB_FF_edges(MeshData *m_d) {
	//Split BB,FF edges if they have sign crossings

	int i, face_index;
	BMIter iter_f, iter_e;
	BMEdge *e;
	BMFace *f;
	BMVert *v1, *v2;
	float v1_u, v1_v, v2_u, v2_v;
	bool is_B;
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);
	int initial_edges = BM_mesh_elem_count(m_d->bm, BM_EDGE);

	//Do 10 samples but don't check end and start point
	float step = 1.0f/11.0f;
	float step_arr[] = { step*5.0f, step*6.0f, step*4.0f, step*7.0f, step*3.0f,
		step*8.0f, step*2.0f, step*9.0f, step*1.0f, step*10.0f };

	BM_ITER_MESH_INDEX (e, &iter_e, m_d->bm, BM_EDGES_OF_MESH, i) {
		Vert_buf v_buf;

		if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

		is_B = calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no);

		if( is_B != calc_if_B_nor(m_d->cam_loc, e->v2->co, e->v2->no) ){
			//This is not a FF or BB edge
			continue;
		}

		if( i < orig_edges ){
			//This edge exists on the original mesh
			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
			//why is table_ensure not fixing the assert?
			BMEdge *orig_e = BM_edge_at_index_find(m_d->bm_orig, i);

			//Get face connected to edge from orig mesh
			//TODO is it wise to use BM_ITER_ELEM here?
			BM_ITER_ELEM (f, &iter_f, orig_e, BM_FACES_OF_EDGE) {
				//Get first face
				break;
			}

			face_index = BM_elem_index_get(f);

			v1 = orig_e->v1;
			v2 = orig_e->v2;

			v_buf.orig_edge = orig_e;
			v_buf.orig_face = f;
		} else {
			BMVert *vert_arr[2];

			//This should be safe because the vert count is still the same as the original mesh.
			v1 = BLI_ghash_lookup(m_d->vert_hash, e->v1);
			v2 = BLI_ghash_lookup(m_d->vert_hash, e->v2);

			vert_arr[0] = v1;
			vert_arr[1] = v2;

			//TODO add checks if to hande if there is no face
			f = BM_face_exists_overlap(vert_arr, 2);
			face_index = BM_elem_index_get(f);

			v_buf.orig_edge = NULL;
			v_buf.orig_face = f;
		}

		//TODO can I just check each vert once?
		get_uv_coord(v1, f, &v1_u, &v1_v);
		get_uv_coord(v2, f, &v2_u, &v2_v);

		{
			int j;
			float u, v;
			float P[3], du[3], dv[3];

			if( v1_u == v2_u ){
				u = v1_u;

				for(j=0; j < 10; j++){
					v = step_arr[j];
					m_d->eval->evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						break;
					}
				}
			} else if ( v1_v == v2_v ){
				v = v1_v;

				for(j=0; j < 10; j++){
					u = step_arr[j];
					m_d->eval->evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						break;
					}
				}
			} else {
				bool alt_diag;
				if((v1_u == 0 && v1_v == 0) || (v2_u == 0 && v2_v == 0)){
					alt_diag = false;
				} else {
					alt_diag = true;
				}
				for(j=0; j < 10; j++){
					if(alt_diag){
						u = 1.0f - step_arr[j];
					} else {
						u = step_arr[j];
					}

					v = step_arr[j];
					m_d->eval->evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						break;
					}
				}
			}
		}

	}

}

static float get_k_r(struct OpenSubdiv_Evaluator *eval, int face_index, float u, float v, const float cam_loc[3]){
	float du[3], dv[3], dudu[3], dudv[3], dvdv[3], P[3], no[3], d1[3], d2[3];
	float k1, k2;
	float I[2][2], II[2][2];
	eval->evaluateLimit2(eval, face_index, u, v, P, du, dv, dudu, dudv, dvdv);

	cross_v3_v3v3(no, du, dv);
	normalize_v3(no);

	//http://en.wikipedia.org/wiki/Principal_curvature

	I[0][0] = dot_v3v3(du, du);
	I[0][1] = dot_v3v3(du, dv);
	I[1][0] = dot_v3v3(du, dv);
	I[1][1] = dot_v3v3(dv, dv);

	II[0][0] = dot_v3v3(dudu, no);
	II[0][1] = dot_v3v3(dudv, no);
	II[1][0] = dot_v3v3(dudv, no);
	II[1][1] = dot_v3v3(dvdv, no);
	{
		float S[2][2];
		float detI = determinant_m2(I[0][0], I[0][1], I[1][0], I[1][1]);

		if(fabsf(detI) < 1e-14){
			detI = 1e-14;
			//printf("detI near zero!!!\n");
		}

		S[0][0] = (II[0][1]*I[0][1] - II[0][0]*I[1][1]) / detI;
		S[0][1] = (II[1][1]*I[0][1] - II[0][1]*I[1][1]) / detI;
		S[1][0] = (II[0][0]*I[0][1] - II[0][1]*I[0][0]) / detI;
		S[1][1] = (II[0][1]*I[0][1] - II[1][1]*I[0][0]) / detI;

		{
			//TODO perhaps remove pdir2
			float pdir1[2], pdir2[2];
			float traceS = S[0][0] + S[1][1];
			float detS = determinant_m2(S[0][0], S[0][1], S[1][0], S[1][1]);
			float diff = traceS*traceS - 4.0f * detS;

			if(diff >= 0){
				float sqrtDiff = sqrtf(diff);
				k1 = 0.5f * (traceS + sqrtDiff);
				k2 = 0.5f * (traceS - sqrtDiff);
				if(fabsf(k1) < fabsf(k2)){
					float swap = k1;
					k1 = k2;
					k2 = swap;
				}
				if(fabsf(S[1][0]) > 1e-14){
					copy_v2_fl2(pdir1, k1 - S[1][1], S[1][0]);
					copy_v2_fl2(pdir2, k2 - S[1][1], S[1][0]);
				}else if (fabsf(S[0][1]) > 1e-14){
					copy_v2_fl2(pdir1, S[0][1], k1 - S[0][0]);
					copy_v2_fl2(pdir2, S[0][1], k2 - S[0][0]);
				}
				normalize_v2(pdir1);
				normalize_v2(pdir2);
			}else{
				//k1 = 0.0f;
				//k2 = 0.0f;
				//d1 = tanU;
				//d2 = tanV;
				printf("diff neg\n");
				return 0;
			}

			mul_v3_fl(du, pdir1[0]);
			mul_v3_fl(dv, pdir1[1]);
			add_v3_v3v3(d1, du, dv);
			normalize_v3(d1);
			cross_v3_v3v3(d2, no, d1);
			//d2 = d1 ^ limitNormal; //tanU * pdir2[0] + tanV * pdir2[1];
		}
	}
	{
		float view_vec[3], ndotv, sintheta, u2, v2, k_r;

		sub_v3_v3v3(view_vec, cam_loc, P);
		normalize_v3(view_vec);

		ndotv = dot_v3v3(no, view_vec);
		sintheta = 1.0f - ndotv*ndotv;
		u = dot_v3v3(view_vec, d1);
		v = dot_v3v3(view_vec, d2);
		u2 = u*u;
		v2 = v*v;
		k_r = (k1 * u2 + k2 * v2) / sintheta;
		return k_r;
	}
}

static void convert_uv_to_new_face(BMEdge *e, BMFace *old_f, BMFace *f, float *u, float *v){
	//convert the old u/v coords to the new face coord

	float v1_u, v1_v, v2_u, v2_v;
	float old_v1_u, old_v1_v, old_v2_u, old_v2_v;

	if(old_f == f){
		//Already have the correct uv coords
		return;
	}

	get_uv_coord(e->v1, f, &v1_u, &v1_v);
	get_uv_coord(e->v2, f, &v2_u, &v2_v);

	get_uv_coord(e->v1, old_f, &old_v1_u, &old_v1_v);
	get_uv_coord(e->v2, old_f, &old_v2_u, &old_v2_v);

	//Which axis are we moving along?
	if( *u == 0.0f || *u == 1.0f ){
		//Along v axis
		if( v1_u == v2_u ){
			//Still along the v axis in the new face
			if(v1_v != old_v1_v){
				*v = 1.0f - *v;
			}
			*u = v1_u;
		} else {
			//Changed axis to u
			if(v1_u != old_v1_v){
				*u = 1.0f - *v;
			} else {
				*u = *v;
			}
			*v = v1_v;
		}
	} else {
		//Along u axis
		if( v1_v == v2_v ){
			//Still along the u axis in the new face
			if(v1_u != old_v1_u){
				*u = 1.0f - *u;
			}
			*v = v1_v;
		} else {
			//Changed axis to v
			if(v1_v != old_v1_u){
				*v = 1.0f - *u;
			} else {
				*v = *u;
			}
			*u = v1_u;
		}
	}

}

static bool append_vert(BLI_Buffer *C_verts, BMVert *vert){
	//check if vert is already in the buffer
	for(int vert_i = 0; vert_i < C_verts->count; vert_i++){
		if( vert == BLI_buffer_at(C_verts, BMVert*, vert_i)){
			return false;
		}
	}
	BLI_buffer_append(C_verts, BMVert*, vert);
	return true;
}

static bool append_face(BLI_Buffer *faces, BMFace *face){
	//check if vert is already in the buffer
	for(int face_i = 0; face_i < faces->count; face_i++){
		if( face == BLI_buffer_at(faces, BMFace*, face_i)){
			return false;
		}
	}
	BLI_buffer_append(faces, BMFace*, face);
	return true;
}

static bool check_and_shift(BMVert *vert, const float new_loc[3], const float new_no[3], MeshData *m_d){
	//TODO add all shiftability checks from the paper

	float old_loc[3];

	{
		//Check if we will try to shift a pole.
		//If we shift it, it may prevent interplation later on
        //TODO remove this check as the "ab" multi face interpol method is removed
		/*
		BMVert *v = BLI_ghash_lookup(m_d->vert_hash, vert);
		if( v && BM_vert_edge_count(v) > 4 ){
			return false;
		}
		*/
	}
	copy_v3_v3( old_loc, vert->co );

	// Will the shift create folds?
	{

		BMFace* f;
		BMIter iter_f;
		float face_no[3];
		bool fold = false;

		copy_v3_v3( vert->co, new_loc );

		//Check if the new position lies inside any of the vert faces (potential fold if outside)
		BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {

			BM_face_calc_normal(f, face_no);
			// BM_face_point_inside_test is too inaccurate to use here as some overhangs are missed with it.
			if( dot_v3v3(new_no, face_no) < 0.0f ){
				fold = true;
				break;
			}
		}
		copy_v3_v3( vert->co, old_loc );
		if( fold ){
			//printf("Skipped shift vert!\n");
			return false;
		}
	}

	//Will this shift a cusp edge. It does, then it might lead to a CCC triangle.
	{
		BMEdge *edge;
		BMIter iter_e;

		BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
			for(int cusp_i = 0; cusp_i < m_d->cusp_edges->count; cusp_i++){
				Cusp cusp = BLI_buffer_at(m_d->cusp_edges, Cusp, cusp_i);
				if( edge == cusp.cusp_e ){
					return false;
				}
			}
		}
	}

	//Check if the shift might/will cause a CCC face
	{
		BLI_buffer_declare_static(BMVert*, c_verts, BLI_BUFFER_NOP, 32);
		BMFace* f;
		BMIter iter_f;
		int num_cross = 0; //number of zero crossing edges

		BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
			BMEdge* e;
			BMIter iter_e;
			BM_ITER_ELEM (e, &iter_e, f, BM_EDGES_OF_FACE) {
				if( e->v1 != vert && e->v2 != vert ){
					int vert_i;
					int edge_c_verts = 0;
					for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
						BMVert* C_vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);
						if( e->v1 == C_vert ) {
							edge_c_verts++;

							if( append_vert( &c_verts, e->v1 ) ){
								num_cross++;
							}
						}
						if( e->v2 == C_vert ){
							edge_c_verts++;

							if( append_vert( &c_verts, e->v2 ) ){
								num_cross++;
							}
						}
					}

					if ( edge_c_verts >= 2 ){
						//TODO if > 2 then we have added duplicate verts to C_verts, should not happen!
						BLI_buffer_free(&c_verts);
						return false;
					}

					if( edge_c_verts == 0 && calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no) != calc_if_B_nor(m_d->cam_loc, e->v2->co, e->v2->no) ){
						num_cross++;
					}

					if( num_cross > 2 ) {
						BLI_buffer_free(&c_verts);
						return false;
					}

				}
			}
		}
		BLI_buffer_free(&c_verts);
	}
	//Adjust vert normal to the limit normal
	copy_v3_v3(vert->no, new_no);
	copy_v3_v3(vert->co, new_loc);
	return true;
}

static void mult_face_search( BMFace *f, BMFace *f2, BMEdge *e, const float v1_uv[2], const float v2_uv[2], MeshData *m_d ){
	//Create a list of faces that should be used when searching for the split
	BMVert *vert;
	BMFace *face;
	BMIter iter_f, iter_v;
	//There must be some overlap between the faces connected to f and f2. Otherwise this might fail
	bool found_overlap = false;
	float edge_no[3], no[3];

	interp_v3_v3v3(edge_no, e->v1->no, e->v2->no, 0.5f);

	BLI_buffer_declare_static(BMFace*, faces, BLI_BUFFER_NOP, 32);

	//First face
	BM_ITER_ELEM (vert, &iter_v, f, BM_VERTS_OF_FACE) {
		BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
			BM_face_calc_normal(face, no);
			//Do not add faces that can create overlaps in the projection mapping
			if( dot_v3v3(no, edge_no) > 0.0f ){
				append_face(&faces, face);
			}
		}
	}

	//Keep track of what faces were added from f2
    int f2_start_idx = faces.count;

	//Second face
	BM_ITER_ELEM (vert, &iter_v, f2, BM_VERTS_OF_FACE) {
		BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
			bool dup_face = false;
			//check if face is already in the buffer
			for(int face_i = 0; face_i < faces.count; face_i++){
				if( face == BLI_buffer_at(&faces, BMFace*, face_i)){
					if( face_i < f2_start_idx ){
						//This face was added from f. We have overlap!
						found_overlap = true;
					}
					dup_face = true;
					break;
				}
			}
			if (!dup_face){
				BM_face_calc_normal(face, no);
				//Do not add faces that can create overlaps in the projection mapping
				if( dot_v3v3(no, edge_no) > 0.0f ){
					append_face(&faces, face);
				}
				BLI_buffer_append(&faces, BMFace*, face);
			}
		}
	}

	if( !found_overlap ){
		//We can't easily interpolate this edge, do not try to insert a new vertex here
		//TODO in theory, this should never happen. Should probably exit if it does.
		printf("Couldn't find any suitable interpolation face map!\n");
		return;
	}

	{

		//Now we can begin interpolating along the edge
		{
			float face_dir, uv_P[2], P[3], du[3], dv[3], new_no[3];
			float step = 0.5f;
			float step_len = 0.25f;
			int face_index;
			bool v1_face_dir = calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no);
			BMFace *cur_face;

			float mat[3][3];
			float start[2], end[2], cur_v2[2];

			axis_dominant_v3_to_m3(mat, edge_no);
			get_st_point(f, v1_uv[0], v1_uv[1], mat, start);
			get_st_point(f2, v2_uv[0], v2_uv[1], mat, end);

			for( int i = 0; i < 10; i++ ){
				interp_v2_v2v2(cur_v2, start, end, step);

				for(int face_i = 0; face_i < faces.count; face_i++){
						cur_face = BLI_buffer_at(&faces, BMFace*, face_i);
						if( point_inside_v2( mat, cur_v2, cur_face ) ){
							get_uv_point( cur_face, uv_P, cur_v2, mat );
							break;
						}
				}
				face_index = BM_elem_index_get(cur_face);
				m_d->eval->evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

				face_dir = get_facing_dir(m_d->cam_loc, P, du, dv);

				if( fabs(face_dir) < 1e-14 ){
					//We got lucky and found the zero crossing!
					printf("--->> got lucky\n");
					break;
				}

				if( (face_dir < 0) == v1_face_dir ){
					step += step_len;
				} else {
					step -= step_len;
				}
				step_len = step_len/2.0f;
			}

			cross_v3_v3v3(new_no, du, dv);
			normalize_v3(new_no);

			if( len_v3v3(P, e->v1->co) < BM_edge_calc_length(e) * 0.2f ){
				if(check_and_shift(e->v1, P, new_no, m_d) ){
					//Do not insert a new vert here, shift it instead
					append_vert(m_d->C_verts, e->v1);
					add_shifted_vert( e->v1, cur_face, uv_P, m_d );
					return;
				}
			} else if (len_v3v3(P, e->v2->co) < BM_edge_calc_length(e) * 0.2f ){
				if(check_and_shift(e->v2, P, new_no, m_d) ){
				//Do not insert a new vert here, shift it instead
				append_vert(m_d->C_verts, e->v2);
				add_shifted_vert( e->v2, cur_face, uv_P, m_d );
				return;
				}
			}

			{
				//Insert a new vert

				Vert_buf new_buf;

				BMVert *c_vert = split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
				append_vert(m_d->C_verts, c_vert);

				new_buf.orig_face = cur_face;
				new_buf.orig_edge = NULL;
				new_buf.u = uv_P[0];
				new_buf.v = uv_P[1];
				BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
			}
		}
	}
	BLI_buffer_free(&faces);

}

static bool bisect_search(const float v1_uv[2], const float v2_uv[2], BMEdge *e, BMFace *orig_face, float uv_result[2], MeshData *m_d ){
	//Search edge for sign crossing and split it!
	int i;
	float face_dir, uv_P[2], P[3], du[3], dv[3], new_no[3];
	float step = 0.5f;
	float step_len = 0.25f;
	bool v1_face_dir = calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no);
	int face_index = BM_elem_index_get(orig_face);

	for( i = 0; i < 10; i++){
		interp_v2_v2v2( uv_P, v1_uv, v2_uv, step);
		m_d->eval->evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);
		face_dir = get_facing_dir(m_d->cam_loc, P, du, dv);

		if( fabs(face_dir) < 1e-14 ){
			//We got lucky and found the zero crossing!
			printf("--->> got lucky\n");
			break;
		}

		if( (face_dir < 0) == v1_face_dir ){
			step += step_len;
		} else {
			step -= step_len;
		}
		step_len = step_len/2.0f;
	}

	cross_v3_v3v3(new_no, du, dv);
	normalize_v3(new_no);

	if( len_v3v3(P, e->v1->co) < BM_edge_calc_length(e) * 0.2f ){
		if( check_and_shift(e->v1, P, new_no, m_d) ){
			//Do not insert a new vert here, shift it instead
			append_vert(m_d->C_verts, e->v1);
			add_shifted_vert( e->v1, orig_face, uv_P, m_d );
			return false;
		}
	} else if (len_v3v3(P, e->v2->co) < BM_edge_calc_length(e) * 0.2f ){
		if( check_and_shift(e->v2, P, new_no, m_d) ){
			//Do not insert a new vert here, shift it instead
			append_vert(m_d->C_verts, e->v2);
			add_shifted_vert( e->v2, orig_face, uv_P, m_d );
			return false;
		}
	}

	copy_v2_v2(uv_result, uv_P);
	{
		BMVert *vert = split_edge_and_move_vert(m_d, m_d->bm, e, P, du, dv);
		append_vert(m_d->C_verts, vert);
	}
	return true;
}

static void search_edge( const int i, BMEdge *e, MeshData *m_d){

	float v1_u, v1_v, v2_u, v2_v;
	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);
	int v1_idx = BM_elem_index_get(e->v1);
	int v2_idx = BM_elem_index_get(e->v2);
	Vert_buf v_buf1, v_buf2;
	BMFace *f, *f2;
	BMEdge *orig_e = NULL;
	BMVert *v1 = NULL, *v2 = NULL;
	bool v1_has_face = false, v2_has_face = false, diff_faces = false;

	v1 = BLI_ghash_lookup(m_d->vert_hash, e->v1);

	if( v1 == NULL ){
		v_buf1 = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v1_idx - orig_verts);
		v1_u = v_buf1.u;
		v1_v = v_buf1.v;
		if( v_buf1.orig_edge == NULL ){
			v1_has_face = true;
		}
	}

	v2 = BLI_ghash_lookup(m_d->vert_hash, e->v2);

	if( v2 == NULL ){
		v_buf2 = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v2_idx - orig_verts);
		v2_u = v_buf2.u;
		v2_v = v_buf2.v;
		if( v_buf2.orig_edge == NULL ){
			v2_has_face = true;
		}
	}

	if( v1 && v2 ){

		if( i < orig_edges ){
			//this edge is on the original mesh
			BMIter iter_f;

			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
			//why is table_ensure not fixing the assert?
			orig_e = BM_edge_at_index_find(m_d->bm_orig, i);

			//Get face connected to edge from orig mesh
			//TODO is it wise to use BM_ITER_ELEM here?
			BM_ITER_ELEM (f, &iter_f, orig_e, BM_FACES_OF_EDGE) {
				//Get first face
				break;
			}
		} else {
			BMVert *vert_arr[] = {v1 ,v2};
			f = BM_face_exists_overlap(vert_arr, 2);
		}
		get_uv_coord(v1, f, &v1_u, &v1_v);
		get_uv_coord(v2, f, &v2_u, &v2_v);
	} else if ( v1 ){
		if( v2_has_face ){
			f = v_buf2.orig_face;
		} else {
			BMVert *vert_arr[3];

			vert_arr[0] = v_buf2.orig_edge->v1;
			vert_arr[1] = v_buf2.orig_edge->v2;
			vert_arr[2] = v1;
			//TODO check if get face fails
			if( v_buf2.orig_edge->v1 == v1 || v_buf2.orig_edge->v2 == v1 ){
				f = BM_face_exists_overlap(vert_arr, 2);
			} else {
				f = BM_face_exists_overlap(vert_arr, 3);
			}
			convert_uv_to_new_face( v_buf2.orig_edge, v_buf2.orig_face, f, &v2_u, &v2_v);
		}
		get_uv_coord(v1, f, &v1_u, &v1_v);
	} else if ( v2 ){
		if( v1_has_face ){
			f = v_buf1.orig_face;
		} else {
			BMVert *vert_arr[3];

			vert_arr[0] = v_buf1.orig_edge->v1;
			vert_arr[1] = v_buf1.orig_edge->v2;
			vert_arr[2] = v2;
			//TODO check if get face fails
			if( v_buf1.orig_edge->v1 == v2 || v_buf1.orig_edge->v2 == v2 ){
				f = BM_face_exists_overlap(vert_arr, 2);
			} else {
				f = BM_face_exists_overlap(vert_arr, 3);
			}
			convert_uv_to_new_face( v_buf1.orig_edge, v_buf1.orig_face, f, &v1_u, &v1_v);
		}
		get_uv_coord(v2, f, &v2_u, &v2_v);
	} else {
		if( v1_has_face || v2_has_face ){
			if( v1_has_face && v2_has_face ){
				if( v_buf1.orig_face != v_buf2.orig_face ){
					diff_faces = true;
					f2 = v_buf2.orig_face;
				}
				f = v_buf1.orig_face;
			} else if ( v1_has_face ) {
				f = v_buf1.orig_face;
				convert_uv_to_new_face( v_buf2.orig_edge, v_buf2.orig_face, f, &v2_u, &v2_v);
			} else {
				f = v_buf2.orig_face;
				convert_uv_to_new_face( v_buf1.orig_edge, v_buf1.orig_face, f, &v1_u, &v2_v);
			}
		} else {
			//No orig face. So this in on a orig edge. So just get the face from the v1 edge
			BMVert *vert_arr[] = {v_buf1.orig_edge->v1 ,v_buf1.orig_edge->v2};
			f = BM_face_exists_overlap(vert_arr, 2);
			convert_uv_to_new_face( v_buf1.orig_edge, v_buf1.orig_face, f, &v1_u, &v1_v);
			convert_uv_to_new_face( v_buf2.orig_edge, v_buf2.orig_face, f, &v2_u, &v2_v);
		}
	}

	{
		//TODO rewrite the above checks so that we do not need to do this check
		//Check if one of the verts has been shifted or not
		Vert_buf* vert1 = get_shift_vert( e->v1, m_d );
		Vert_buf* vert2 = get_shift_vert( e->v2, m_d );

		if( vert1 != NULL ){
			if(vert1->orig_face != f){
				if(diff_faces){
					f = vert1->orig_face;
				} else {
					diff_faces = true;
					f2 = f;
					f = vert1->orig_face;
				}
			}
			v1_u = vert1->u;
			v1_v = vert1->v;
		}

		if( vert2 != NULL ){
			if(diff_faces){
				if(vert2->orig_face != f2){
					f2 = vert2->orig_face;
				}
			} else if(vert2->orig_face != f) {
				diff_faces = true;
				f2 = vert2->orig_face;
			}

			v2_u = vert2->u;
			v2_v = vert2->v;
		}

		//Check if a shifted vert caused the verts to be on the same face
		if(diff_faces && f == f2){
			diff_faces = false;
		}
	}

	{
		Vert_buf new_buf;
		float uv_result[2];
		float v1_uv[2] = { v1_u, v1_v };
		float v2_uv[2] = { v2_u, v2_v };

		if(diff_faces) {
			//The edge spawns over multiple original edges, try to interpolate along this edge.
			//If it fails, do not insert any new verts here
			//printf("Mult face search\n");
			mult_face_search( f, f2, e, v1_uv, v2_uv, m_d );
			return;
		}

		if( (v1_u == 0 && v2_u == 0) || (v1_u == 1 && v2_u == 1) ||
			(v1_v == 0 && v2_v == 0) || (v1_v == 1 && v2_v == 1) )
		{
			//Along an original edge, save orig face for uv conversion
			new_buf.orig_face = f;
			if( v1 && v2 ){
				new_buf.orig_edge = orig_e;
			} else if ( v1 ){
				new_buf.orig_edge = v_buf2.orig_edge;
			} else {
				new_buf.orig_edge = v_buf1.orig_edge;
			}

		} else {
			new_buf.orig_face = f;
			new_buf.orig_edge = NULL;
		}

		if( bisect_search( v1_uv, v2_uv, e, f, uv_result, m_d) ){
			//if a new vert is inserted add it to the buffer
			new_buf.u = uv_result[0];
			new_buf.v = uv_result[1];
			BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
		}
	}
}

static void contour_insertion( MeshData *m_d ) {
	int i, cusp_i;
	BMEdge *e;
	BMIter iter_e;

	int initial_edges = BM_mesh_elem_count(m_d->bm, BM_EDGE);

	BM_ITER_MESH_INDEX (e, &iter_e, m_d->bm, BM_EDGES_OF_MESH, i) {
		if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

		{
			bool cont = false;
			//Check if this is a cusp edge
			for(cusp_i = 0; cusp_i < m_d->cusp_edges->count; cusp_i++){
				Cusp cusp = BLI_buffer_at(m_d->cusp_edges, Cusp, cusp_i);
				if(cusp.cusp_e == e){
					//Do not split this edge yet
					//printf("skipped cusp edge\n");
					cont = true;
					break;
				}
			}
			if( cont ){
				continue;
			}
		}

		if( calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no) == calc_if_B_nor(m_d->cam_loc, e->v2->co, e->v2->no) ){
			//This is not a FB or BF edge
			continue;
		}

		//Check if the edge already has a C vert
		{
			int vert_i;
			bool skip_edge = false;
			for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
				BMVert* C_vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);
				if( e->v1 == C_vert || e->v2 == C_vert){
					skip_edge = true;
				}
			}

			if(skip_edge) {
				continue;
			}
		}
		search_edge( i, e, m_d );
	}
}

static bool sign_cross(const bool bool_arr[3]) {
	int i;
	bool temp = bool_arr[0];
	for (i = 1; i < 3; i++) {
		if (temp != bool_arr[i]) {
			return true;
		}
	}
	return false;
}

static bool cusp_triangle(struct OpenSubdiv_Evaluator *eval, const float cam_loc[3], const int face_index, Cusp_triang *c_tri, Cusp *cusp){

		GSQueue *tri_que = BLI_gsqueue_new(sizeof(Cusp_triang));
        BLI_gsqueue_push_back(tri_que, c_tri);

		//Add this because it seems to get stuck sometimes because the triangles never becomes small enough
		//TODO maybe find a better end condition?
		int iteration = 0;

		while ( !BLI_gsqueue_is_empty(tri_que) ){
            Cusp_triang cur_tri;
			BLI_gsqueue_pop(tri_que, &cur_tri);

			if( !( sign_cross(cur_tri.b_arr) && sign_cross(cur_tri.kr_arr) ) ){
				continue;
			}

            iteration++;

			if( area_tri_v3(cur_tri.co_arr[0], cur_tri.co_arr[1], cur_tri.co_arr[2]) <= 1e-14 ||
				iteration > 1000){
				float cusp_co[3];
				float cusp_no[3];
				float du[3], dv[3];

				float uv1[3] = { cur_tri.u_arr[0], cur_tri.v_arr[0], 0 };
				float uv2[3] = { cur_tri.u_arr[1], cur_tri.v_arr[1], 0 };
				float uv3[3] = { cur_tri.u_arr[2], cur_tri.v_arr[2], 0 };

				float res_uv[3];

				mid_v3_v3v3v3(res_uv, uv1, uv2, uv3);
				eval->evaluateLimit(eval, face_index, res_uv[0], res_uv[1], cusp_co, du, dv);
				cross_v3_v3v3(cusp_no, du, dv);
				normalize_v3(cusp_no);

				copy_v3_v3(cusp->cusp_co, cusp_co);
				copy_v3_v3(cusp->cusp_no, cusp_no);
				cusp->u = res_uv[0];
				cusp->v = res_uv[1];
				BLI_gsqueue_free(tri_que);
				return true;
			}

            int best_edge = 0;

			{
				float edge_len = len_v3v3(cur_tri.co_arr[0], cur_tri.co_arr[1]);
				for( int i = 1; i < 3; i++ ){
                    float len = len_v3v3(cur_tri.co_arr[i], cur_tri.co_arr[(i+1)%3]);
					if( len > edge_len ){
						best_edge = i;
						edge_len = len;
					}
				}
			}

			float uv_1[] = {cur_tri.u_arr[best_edge], cur_tri.v_arr[best_edge]};
			float uv_2[] = {cur_tri.u_arr[(best_edge+1)%3], cur_tri.v_arr[(best_edge+1)%3]};
			float new_uv[2], new_co[3], du[3], dv[3];

			interp_v2_v2v2( new_uv, uv_1, uv_2, 0.5f);
			eval->evaluateLimit(eval, face_index, new_uv[0], new_uv[1], new_co, du, dv);
			bool new_b = calc_if_B(cam_loc, new_co, du, dv);
			bool new_kr = get_k_r(eval, face_index, new_uv[0], new_uv[1], cam_loc) > 0;

			for( int i = 0; i < 2; i++ ){
				Cusp_triang new_tri;
				copy_v3_v3((float*)new_tri.b_arr, (float*)cur_tri.b_arr);
				copy_v3_v3((float*)new_tri.kr_arr, (float*)cur_tri.kr_arr);
				copy_v3_v3(new_tri.u_arr, cur_tri.u_arr);
				copy_v3_v3(new_tri.v_arr, cur_tri.v_arr);
				copy_v3_v3(new_tri.co_arr[0], cur_tri.co_arr[0]);
				copy_v3_v3(new_tri.co_arr[1], cur_tri.co_arr[1]);
				copy_v3_v3(new_tri.co_arr[2], cur_tri.co_arr[2]);

				new_tri.b_arr[ (best_edge + i)%3 ] = new_b;
				new_tri.kr_arr[ (best_edge + i)%3 ] = new_kr;
				new_tri.u_arr[ (best_edge + i)%3 ] = new_uv[0];
				new_tri.v_arr[ (best_edge + i)%3 ] = new_uv[1];
				copy_v3_v3(new_tri.co_arr[ (best_edge + i)%3 ], new_co);

				BLI_gsqueue_push_back(tri_que, &new_tri);
			}

		}

		BLI_gsqueue_free(tri_que);
		return false;
}

static BMFace *get_orig_face(int orig_verts, BMVert *vert_arr_in[3], float u_arr[3], float v_arr[3], float co_arr[3][3], MeshData *m_d){
	int i;

	BMEdge *edge_arr[] = {NULL, NULL, NULL};
	BMFace *edge_face_arr[] = {NULL, NULL, NULL};
	BMFace *orig_face = NULL;
	BMVert *vert_arr[3] = {vert_arr_in[0], vert_arr_in[1], vert_arr_in[2]};

	//check if all verts are on the orignal mesh
	for(i = 0; i < 3; i++){
		int v_idx = BM_elem_index_get(vert_arr[i]);
        BMVert *temp_v;

		//Copy coords for later use
		copy_v3_v3(co_arr[i], vert_arr[i]->co);

        temp_v = BLI_ghash_lookup(m_d->vert_hash, vert_arr[i]);

		if( temp_v == NULL ){
			Vert_buf v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);
			u_arr[i] = v_buf.u;
			v_arr[i] = v_buf.v;
			if( v_buf.orig_edge == NULL ){
				orig_face = v_buf.orig_face;
				vert_arr[i] = NULL;
			} else {
				//Make sure we don't pick one of the verts that we already have
				int idx1 = BM_elem_index_get(v_buf.orig_edge->v1);
				BMVert *temp_v1 = vert_arr[ mod_i(i-1, 3) ];
				BMVert *temp_v2 = vert_arr[ mod_i(i+1, 3) ];

				if( ( temp_v1 != NULL && idx1 != BM_elem_index_get(temp_v1) ) &&
						( temp_v2 != NULL && idx1 != BM_elem_index_get(temp_v2) ) )
				{
					vert_arr[i] = v_buf.orig_edge->v1;
				} else {
					//If v1 is a duplicate then v2 has to be unique
					vert_arr[i] = v_buf.orig_edge->v2;
				}

				edge_arr[i] = v_buf.orig_edge;
				edge_face_arr[i] = v_buf.orig_face;
			}
		} else {
			vert_arr[i] = temp_v;
		}
	}

	if(orig_face == NULL){
		orig_face = BM_face_exists_overlap(vert_arr, 3);
	}

	for(i = 0; i < 3; i++){
		if(vert_arr[i] == NULL){
			continue;
		}

		if(edge_arr[i] != NULL){
			//Make use we have the correct uv coords
			convert_uv_to_new_face( edge_arr[i], edge_face_arr[i], orig_face, &u_arr[i], &v_arr[i]);
		} else {
			get_uv_coord(vert_arr[i], orig_face, &u_arr[i], &v_arr[i]);
		}
	}
	return orig_face;
}

static void cusp_detection( MeshData *m_d ){
	BMFace *f;
	BMIter iter_f;

	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);
	int initial_faces = BM_mesh_elem_count(m_d->bm, BM_FACE);
	int f_idx;

	BM_ITER_MESH_INDEX (f, &iter_f, m_d->bm, BM_FACES_OF_MESH, f_idx) {
		BMVert *vert;
		BMVert *vert_arr[3];
		BMIter iter_v;
		int vert_idx;
		bool first_vert, back_face, found_face, b_arr[3];
		first_vert = true;
		found_face = false;

		if( !(f_idx < initial_faces) ){
			//TODO perhaps insert every cusp edge after we iterated over all faces instead?
			//this might not be ok because will miss some faces that gets split
			break;
		}

		BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, vert_idx) {
			if(first_vert){
				first_vert = false;
				back_face = calc_if_B_nor(m_d->cam_loc, vert->co, vert->no);
				b_arr[vert_idx] = back_face;
			} else {
				//If one or more of the verts do not have the same facing, then we want to look for cusps
				bool temp = calc_if_B_nor(m_d->cam_loc, vert->co, vert->no);
				b_arr[vert_idx] = temp;
				if(temp != back_face){
					found_face = true;
				}
			}
			vert_arr[vert_idx] = vert;
		}

		if(!found_face){
			continue;
		}

		{
			//See if we are trying to insert a cusp on a face that already has a detected
			//cusp edge. Skip searching for cusps if this is the case.
			int edge_i;
			BMEdge *edge;
			BMIter iter_e;

			bool found_cusp_edge = false;

			BM_ITER_ELEM (edge, &iter_e, f, BM_EDGES_OF_FACE) {
				for(edge_i = 0; edge_i < m_d->cusp_edges->count; edge_i++){
					BMEdge* cusp_edge = BLI_buffer_at(m_d->cusp_edges, Cusp, edge_i).cusp_e;
					if( edge == cusp_edge ){
						found_cusp_edge = true;
						break;
					}
				}
				if( found_cusp_edge ){
					break;
				}
			}

			if( found_cusp_edge ){
				continue;
			}
		}

		//Find original mesh face + uv coords
		{
			float u_arr[3]; //array for u-coords (v1_u, v2_u ...)
			float v_arr[3];
			float co_arr[3][3];
			BMFace *orig_face = get_orig_face(orig_verts, vert_arr, u_arr, v_arr, co_arr, m_d);

			{
				int face_index = BM_elem_index_get(orig_face);
				//Check for k_r sign crossings
				float k_r1 = get_k_r(m_d->eval, face_index, u_arr[0], v_arr[0], m_d->cam_loc);
				float k_r2 = get_k_r(m_d->eval, face_index, u_arr[1], v_arr[1], m_d->cam_loc);
				float k_r3 = get_k_r(m_d->eval, face_index, u_arr[2], v_arr[2], m_d->cam_loc);
				bool k_r_crossing = false;
				if( (k_r1 > 0) != (k_r2 > 0) ){
					//k_r sign crossing!
					//printf("found k_r sign crossing\n");
					k_r_crossing = true;
				} else {
					//check last vert
					if( (k_r1 > 0) != (k_r3 > 0) ){
						//printf("found k_r sign crossing\n");
						k_r_crossing = true;
					}
				}

				if(k_r_crossing){
					Cusp cusp;
					Cusp_triang c_tri;
					cusp.orig_face = orig_face;

                    copy_v3_v3((float*)c_tri.b_arr, (float*)b_arr);
                    copy_v3_v3(c_tri.u_arr, u_arr);
                    copy_v3_v3(c_tri.v_arr, v_arr);
                    copy_v3_v3(c_tri.co_arr[0], co_arr[0]);
                    copy_v3_v3(c_tri.co_arr[1], co_arr[1]);
                    copy_v3_v3(c_tri.co_arr[2], co_arr[2]);
					c_tri.kr_arr[0] = k_r1 > 0;
					c_tri.kr_arr[1] = k_r2 > 0;
					c_tri.kr_arr[2] = k_r3 > 0;

					//Start looking for the cusp in the triangle
					if(cusp_triangle(m_d->eval, m_d->cam_loc, face_index, &c_tri, &cusp)){
						//We found a cusp!
						float uv_1[2], uv_2[2], uv_3[2];
						BMEdge *edge;
						BMVert *cusp_e_vert;

						//printf("Found a cusp point!\n");
						if(b_arr[0] == b_arr[1]){
							uv_1[0] = u_arr[0];
							uv_2[0] = u_arr[1];
							uv_3[0] = u_arr[2];

							uv_1[1] = v_arr[0];
							uv_2[1] = v_arr[1];
							uv_3[1] = v_arr[2];
							edge = BM_edge_exists( vert_arr[0], vert_arr[1] );
							cusp_e_vert = vert_arr[2];
						} else if(b_arr[0] == b_arr[2]){
							uv_1[0] = u_arr[0];
							uv_2[0] = u_arr[2];
							uv_3[0] = u_arr[1];

							uv_1[1] = v_arr[0];
							uv_2[1] = v_arr[2];
							uv_3[1] = v_arr[1];
							edge = BM_edge_exists( vert_arr[0], vert_arr[2] );
							cusp_e_vert = vert_arr[1];
						} else {
							uv_1[0] = u_arr[1];
							uv_2[0] = u_arr[2];
							uv_3[0] = u_arr[0];

							uv_1[1] = v_arr[1];
							uv_2[1] = v_arr[2];
							uv_3[1] = v_arr[0];
							edge = BM_edge_exists( vert_arr[1], vert_arr[2] );
							cusp_e_vert = vert_arr[0];
						}

						{
							float P[3], du[3], dv[3];
							float edge_uv[2];
							float cusp_uv[2] = {cusp.u, cusp.v};
							Vert_buf v_buf;
							int edge_idx = BM_elem_index_get(edge);
							int v1_idx = BM_elem_index_get(edge->v1);
							int v2_idx = BM_elem_index_get(edge->v2);

							if( isect_line_line_v2_point( uv_1, uv_2, uv_3, cusp_uv, edge_uv ) != ISECT_LINE_LINE_CROSS ){
								printf("Couldn't find intersection point to edge from cusp!\n");
								//TODO this is a big error so quit instead
								continue;
							}

							m_d->eval->evaluateLimit(m_d->eval, face_index, edge_uv[0], edge_uv[1], P, du, dv);

							float cusp_dist_to_edge1 = dist_to_line_v3(cusp.cusp_co, cusp_e_vert->co, edge->v1->co);
							float cusp_dist_to_edge2 = dist_to_line_v3(cusp.cusp_co, cusp_e_vert->co, edge->v2->co);
							//Check if we should use an existing edge (no new verts)
							if( cusp_dist_to_edge1 < BM_edge_calc_length(edge) * 0.2f ||
								cusp_dist_to_edge2 < BM_edge_calc_length(edge) * 0.2f ){

								BMVert *edge_vert;
								BMEdge *cusp_edge;
								BMIter iter_e;
								float new_no[3];

								//Which edge should we move? (Which is closest to the split point?)
								if( len_v3v3(P, edge->v1->co) < len_v3v3(P, edge->v2->co) ){
									edge_vert = edge->v1;
								} else {
									edge_vert = edge->v2;
								}

								BM_ITER_ELEM (cusp_edge, &iter_e, f, BM_EDGES_OF_FACE) {
									if( cusp_edge != edge && (cusp_edge->v1 == edge_vert || cusp_edge->v2 == edge_vert) ){
										//Found edge
										break;
									}
								}

								cross_v3_v3v3(new_no, du, dv);
								normalize_v3(new_no);

								//Can we shift this vertex?
								if( check_and_shift(edge_vert, P, new_no, m_d) ){
									cusp.cusp_e = cusp_edge;
									add_shifted_vert( edge_vert, orig_face, edge_uv, m_d );
									BLI_buffer_append(m_d->cusp_edges, Cusp, cusp);

									//printf("Used existing edge for cusp!\n");
									continue;
								}
							}

							if( (edge_idx < orig_edges) ){
								//Point on orig edge
								BMEdge *orig_e = BM_edge_at_index_find(m_d->bm_orig, edge_idx);
								v_buf.orig_edge = orig_e;
								v_buf.orig_face = orig_face;
							} else if( edge_uv[0] == 0 || edge_uv[0] == 1 || edge_uv[1] == 0 || edge_uv[1] == 1 ){
								if( (v1_idx + 1) > orig_verts){
									Vert_buf v_buf_old = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v1_idx - orig_verts);
									v_buf.orig_edge = v_buf_old.orig_edge;
									v_buf.orig_face = orig_face;
								} else if( (v2_idx + 1) > orig_verts) {
									Vert_buf v_buf_old = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v2_idx - orig_verts);
									v_buf.orig_edge = v_buf_old.orig_edge;
									v_buf.orig_face = orig_face;
								}
							} else {
								//On orig face
								v_buf.orig_edge = NULL;
								v_buf.orig_face = orig_face;
							}

							v_buf.u = edge_uv[0];
							v_buf.v = edge_uv[1];

							BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
							{
								int j;
								int cur_edge = BM_mesh_elem_count(m_d->bm, BM_EDGE);
								BMEdge *cusp_e = NULL;

								split_edge_and_move_vert(m_d, m_d->bm, edge, P, du, dv);

								//Get the cusp edge
								for(j = 0; j < 3; j++){
									cur_edge++;
									cusp_e = BM_edge_at_index_find( m_d->bm, cur_edge );
									if( cusp_e->v1 == cusp_e_vert || cusp_e->v2 == cusp_e_vert ){
										break;
									}
									cusp_e = NULL;
								}
								if(cusp_e == NULL){
									//printf("No cusp edge found!!!\n");
								} else {
									cusp.cusp_e = cusp_e;
									BLI_buffer_append(m_d->cusp_edges, Cusp, cusp);
								}
							}
						}
					}
				}
			}
		}
	}
}

static void cusp_insertion(MeshData *m_d){
	int cusp_i;

	for(cusp_i = 0; cusp_i < m_d->cusp_edges->count; cusp_i++){
		BMVert *vert;
		Vert_buf new_buf;
		Cusp cusp = BLI_buffer_at(m_d->cusp_edges, Cusp, cusp_i);

		if( len_v3v3(cusp.cusp_co, cusp.cusp_e->v1->co) < BM_edge_calc_length(cusp.cusp_e) * 0.2f ){
			if(  BM_vert_edge_count(cusp.cusp_e->v1) == 4 && check_and_shift(cusp.cusp_e->v1, cusp.cusp_co, cusp.cusp_no, m_d) ){
				float uv_P[2] = { cusp.u, cusp.v };
				append_vert(m_d->C_verts, cusp.cusp_e->v1);
				append_vert(m_d->cusp_verts, cusp.cusp_e->v1);
				add_shifted_vert( cusp.cusp_e->v1, cusp.orig_face, uv_P, m_d );
				continue;
			}
		} else if( len_v3v3(cusp.cusp_co, cusp.cusp_e->v2->co) < BM_edge_calc_length(cusp.cusp_e) * 0.2f ){
			if(  BM_vert_edge_count(cusp.cusp_e->v2) == 4 && check_and_shift(cusp.cusp_e->v2, cusp.cusp_co, cusp.cusp_no, m_d) ){
				float uv_P[2] = { cusp.u, cusp.v };
				append_vert(m_d->C_verts, cusp.cusp_e->v2);
				append_vert(m_d->cusp_verts, cusp.cusp_e->v2);
				add_shifted_vert( cusp.cusp_e->v2, cusp.orig_face, uv_P, m_d );
				continue;
			}
		}

		vert = split_edge_and_move_nor(m_d, m_d->bm, cusp.cusp_e, cusp.cusp_co, cusp.cusp_no);
		append_vert(m_d->C_verts, vert);
		append_vert(m_d->cusp_verts, vert);

		new_buf.orig_face = cusp.orig_face;
		new_buf.orig_edge = NULL;
		new_buf.u = cusp.u;
		new_buf.v = cusp.v;

		BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
	}
}

static bool poke_and_move(BMFace *f, const float new_pos[3], const float du[3], const float dv[3], Radi_vert *r_vert, MeshData *m_d){
	BMVert *vert;

	BMEdge *edge = NULL;
	bool rot_edge = false;
	float mat[3][3];
	float new_norm[3];

	cross_v3_v3v3(new_norm, du, dv);
	normalize_v3(new_norm);

	axis_dominant_v3_to_m3(mat, new_norm);

	// BM_face_point_inside_test is too inaccurate to use here as some overhangs are missed with it.
	if( !point_inside(mat, new_pos, f) ){
		BMIter iter_e;
		BMIter iter_f;
		BMEdge *e;
		BMFace *face;

		rot_edge = true;

		BM_ITER_ELEM (e, &iter_e, f, BM_EDGES_OF_FACE){
			BM_ITER_ELEM (face, &iter_f, e, BM_FACES_OF_EDGE){
				if( face != f ){

					if( point_inside(mat, new_pos, face) ){
						edge = e;
						break;
					}
				}
			}
		}
	}

	if( rot_edge ){
		if( edge == NULL || !BM_edge_rotate_check(edge) ){
			//Do not insert a radial edge here
			return false;
		}
		if( is_vert_in_buffer( edge->v1, m_d->C_verts) && is_vert_in_buffer( edge->v2, m_d->C_verts) ){
			return false;
		}
	}

	{
		//Borrowed from bmo_poke.c
        int i = 0;
		BMFace *f_new;
		BMLoop *l_iter, *l_first;
		/* only interpolate the central loop from the face once,
		 * then copy to all others in the fan */
		BMLoop *l_center_example;

		BMesh *bm = m_d->bm;

		vert = BM_vert_create(bm, new_pos, NULL, BM_CREATE_NOP);

        BM_elem_index_set(vert, m_d->new_vert_idx++);

		//Silence asserts in BM_loop_interp_from_face
		//TODO perhaps work around this in some other way?
        BM_face_normal_update(f);

		l_iter = l_first = BM_FACE_FIRST_LOOP(f);
		do {
			BMLoop *l_new;

			f_new = BM_face_create_quad_tri(bm, l_iter->v, l_iter->next->v, vert, NULL, f, BM_CREATE_NOP);
			l_new = BM_FACE_FIRST_LOOP(f_new);

            BM_face_normal_update(f_new);

			if (i == 0) {
				l_center_example = l_new->prev;
				BM_loop_interp_from_face(bm, l_center_example, f, true, false);
			}
			else {
				BM_elem_attrs_copy(bm, bm, l_center_example, l_new->prev);
			}

			/* Copy Loop Data */
			BM_elem_attrs_copy(bm, bm, l_iter, l_new);
			BM_elem_attrs_copy(bm, bm, l_iter->next, l_new->next);

		} while ((void)i++, (l_iter = l_iter->next) != l_first);

		/* Kill Face */
		BM_face_kill(bm, f);
	}
	//Adjust vert normal to the limit normal
	copy_v3_v3(vert->no, new_norm);

    r_vert->vert = vert;

	if( rot_edge ){
		BM_edge_rotate(m_d->bm, edge, true, 0);
		//printf("rotated edge!\n");
	}

	return true;
}

static void mult_radi_search( BLI_Buffer *diff_f, const float mat[3][3], const float edge1_mid[2], const float edge2_mid[2],
							const float val_1, const float val_2, bool is_B,
							const float rad_plane_no[3], BMVert *C_vert, BMFace *poke_face, MeshData *m_d ){
	//Try to find a vert that is connected to both faces
	BMVert *vert;
	BMFace *face;
	BMIter iter_f, iter_v;
	float poke_face_no[3], no[3];
	zero_v3(poke_face_no);

	BLI_buffer_declare_static(BMFace*, search_faces, BLI_BUFFER_NOP, 32);

	BM_ITER_ELEM (vert, &iter_v, poke_face, BM_VERTS_OF_FACE) {
		//Take the mean normal as face normal. It should represent the surface normal of the surface better
		add_v3_v3(poke_face_no, vert->no);
	}
	mul_v3_fl( poke_face_no, 1.0f / 3.0f );

	for(int f_idx = 0; f_idx < diff_f->count; f_idx++){
		BMFace *f = BLI_buffer_at(diff_f, BMFace*, f_idx);
		BM_ITER_ELEM (vert, &iter_v, f, BM_VERTS_OF_FACE) {
			BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
				//Do not add faces that can create overlaps in the projection mapping
				BM_face_calc_normal(face, no);
				if( dot_v3v3(no, poke_face_no) > 0.0f ){
					append_face(&search_faces, face);
				}
			}
		}
	}

	//Find the faces for our three points
	{
		float uvs[2][2];
		BMFace *face_ids[2];
		float rad_dir[2];
		rad_dir[0] = signf(val_1);
		rad_dir[1] = signf(val_2);

		for (int face_idx = 0; face_idx < search_faces.count; face_idx++) {
			int found_points = 0;
			for (int i = 0; i < 2; i++) {
				float point[2];

				switch(i){
					case 1 :
						copy_v2_v2(point, edge2_mid);
						break;
					default:
						copy_v2_v2(point, edge1_mid);
						break;
				}

				BMFace *f = BLI_buffer_at(&search_faces, BMFace*, face_idx);

				if( point_inside_v2(mat, point, f) ){
					get_uv_point(f, uvs[i], point, mat);

					face_ids[i] = f;
					found_points++;
				}
			}
			if (found_points == 2){
				//We have found all points inside the same face!
				break;
			}
		}

		{
			float search_val, uv_P[2], P[3], du[3], dv[3], temp[3];
			float step = 0.5f;
			float step_len = 0.25f;
			int face_index;
			BMFace *orig_face;
			Vert_buf v_buf;
			/*
			   print_v3("cent", cent);
			   print_v3("edge1_mid", edge1_mid);
			   print_v3("edge2_mid", edge2_mid);
			   print_v3("rad_dir", rad_dir);
			   print_v2("UV_cent", uvs[0]);
			   print_v2("UV_edge1", uvs[1]);
			   print_v2("UV_edge2", uvs[2]);
			   */
			if( face_ids[0] == face_ids[1] ){
				//We can work in pure uv space
				//printf("UV space\n");
				orig_face = face_ids[0];
				face_index = BM_elem_index_get(face_ids[0]);
				for(int i = 0; i < 10; i++ ){
					interp_v2_v2v2( uv_P, uvs[0], uvs[1], step);
					m_d->eval->evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

					sub_v3_v3v3(temp, P, C_vert->co);
					search_val = dot_v3v3(rad_plane_no, temp);

					if( fabs(search_val) < 1e-14 ){
						//We got lucky and found the zero crossing!
						printf("got lucky\n");
						break;
					}

					search_val = signf(search_val);

					if( signf(search_val) == rad_dir[0] ){
						step += step_len;
					} else {
						step -= step_len;
					}

					step_len = step_len/2.0f;
				}
			} else {
				//Work in st space
				float cur_st[3];

				//printf("St space\n");

				for(int i = 0; i < 10; i++ ){
					interp_v2_v2v2(cur_st, edge1_mid, edge2_mid, step);

					for (int face_idx = 0; face_idx < search_faces.count; face_idx++) {
						BMFace *f = BLI_buffer_at(&search_faces, BMFace*, face_idx);
						if( point_inside_v2(mat, cur_st, f) ){
							get_uv_point(f, uv_P, cur_st, mat);

							orig_face = f;
							face_index = BM_elem_index_get(f);
							m_d->eval->evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

							break;
						}
					}

					sub_v3_v3v3(temp, P, C_vert->co);
					search_val = dot_v3v3(rad_plane_no, temp);

					if( fabs(search_val) < 1e-14 ){
						//We got lucky and found the zero crossing!
						printf("got lucky\n");
						break;
					}

					search_val = signf(search_val);

					if( search_val == rad_dir[0] ){
						step += step_len;
					} else {
						step -= step_len;
					}

					step_len = step_len/2.0f;
				}
			}

			v_buf.orig_edge = NULL;
			v_buf.orig_face = orig_face;
			v_buf.u = uv_P[0];
			v_buf.v = uv_P[1];
			Radi_vert r_vert;
			if( poke_and_move(poke_face, P, du, dv, &r_vert, m_d) ){

				r_vert.extendable = true;
				copy_v3_v3(r_vert.radi_plane_no, rad_plane_no);
				r_vert.C_vert = C_vert;
				r_vert.is_B = is_B;

				BLI_buffer_append(m_d->radi_vert_buffer, Radi_vert, r_vert);
				BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
			}
		}

	}

	BLI_buffer_free(&search_faces);
}

static void radial_insertion( MeshData *m_d ){

	int vert_i, vert_j;
	BMFace *f;
	BMIter iter_f;
	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int initial_verts = BM_mesh_elem_count(m_d->bm, BM_VERT);

	for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
		int vert_idx, CC_idx, CC2_idx;
		float co_arr[3][3];
		BMVert *vert_arr[3];
		BMIter iter_v;

		BMVert *vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);

		int face_i;
		int face_count = BM_vert_face_count(vert);
		BMFace **face_arr = BLI_array_alloca(face_arr, face_count);
        /*
		if( BM_elem_index_get(vert) != 318 ){
			continue;
		}
        */
		BM_ITER_ELEM_INDEX (f, &iter_f, vert, BM_FACES_OF_VERT, face_i) {
			face_arr[face_i] = f;
		}

		for( face_i = 0; face_i < face_count; face_i++){
			BMVert *cur_vert;
			CC2_idx = -1;
			bool skip_face = false;

			f = face_arr[face_i];

			BM_ITER_ELEM_INDEX (cur_vert, &iter_v, f, BM_VERTS_OF_FACE, vert_idx) {
				if( BM_elem_index_get(cur_vert) > (initial_verts - 1) ){
					// This is a face we already worked on
					skip_face = true;
					break;
				}

				//Do not attempt to insert radial edges on the "inside" region of the cusp
				//That is, if the opposite edge as a cusp vert, do not try to insert a radial edge here.
				/*
				if (cur_vert != vert && is_vert_in_buffer(cur_vert, m_d->cusp_verts)){
					skip_face = true;
					break;
				}
				*/

				for(vert_j = 0; vert_j < m_d->C_verts->count; vert_j++){
					BMVert *vert2 = BLI_buffer_at(m_d->C_verts, BMVert*, vert_j);

					if( cur_vert == vert){
						CC_idx = vert_idx;
						break;
					} else if ( cur_vert == vert2 ){
						CC2_idx = vert_idx;
						break;
					}
				}

				vert_arr[vert_idx] = cur_vert;
				copy_v3_v3(co_arr[vert_idx], cur_vert->co);
			}

			if( skip_face ){
				continue;
			}

			{
				float val_1, val_2;
				float temp[3];

				float cam_vec[3], rad_plane_no[3];

				sub_v3_v3v3(cam_vec, m_d->cam_loc, co_arr[CC_idx]);
				cross_v3_v3v3(rad_plane_no, vert_arr[CC_idx]->no, cam_vec);
				normalize_v3(rad_plane_no);

				{
					//check if the radial plane intersects the vert's opposite edge
					sub_v3_v3v3(temp, co_arr[ mod_i(CC_idx-1, 3) ], co_arr[CC_idx]);

					val_1 = dot_v3v3(rad_plane_no, temp);

					sub_v3_v3v3(temp, co_arr[ mod_i(CC_idx+1, 3) ], co_arr[CC_idx]);

					val_2 = dot_v3v3(rad_plane_no, temp);

					if( signf(val_1) == signf(val_2) ){
						//Do not insert a radial edge here
						continue;
					}
				}

				if( CC2_idx != -1 ){
					//This face has an CC edge
					//Do the radial planes intersect?
					float rad_plane_no2[3];
					float val2_1, val2_2;

					sub_v3_v3v3(cam_vec, m_d->cam_loc, co_arr[CC2_idx]);
					cross_v3_v3v3(rad_plane_no2, vert_arr[CC2_idx]->no, cam_vec);
					normalize_v3(rad_plane_no2);

					//check if the radial plane intersects the vert's opposite edge
					sub_v3_v3v3(temp, co_arr[ mod_i(CC2_idx-1, 3) ], co_arr[CC2_idx]);

					val2_1 = dot_v3v3(rad_plane_no2, temp);

					sub_v3_v3v3(temp, co_arr[ mod_i(CC2_idx+1, 3) ], co_arr[CC2_idx]);

					val2_2 = dot_v3v3(rad_plane_no2, temp);

					if( signf(val2_1) != signf(val2_2) ){
						//TODO Implement this edge case
						printf("Radial intersect!\n");
					}
				}

				{
					float sts[3][2]; //2d coords in non uv space
					float mat[3][3];

					axis_dominant_v3_to_m3(mat, vert->no);
					BLI_buffer_declare_static(BMFace*, faces, BLI_BUFFER_NOP, 32);
					//Check if the triangle has been shifted so we can't use the original face for UV coords
					for(int i = 0; i < 3; i++){
						Vert_buf* shift_vert = get_shift_vert( vert_arr[i], m_d );
						if( shift_vert != NULL ){
							//This vert has been shifted
							get_st_point(shift_vert->orig_face, shift_vert->u, shift_vert->v, mat, sts[i]);
							append_face(&faces, shift_vert->orig_face);
						} else {
							//Check if edge verts doesn't belong to orig_face
							int v_idx = BM_elem_index_get(vert_arr[i]);
							if( (v_idx + 1) > orig_verts){
								Vert_buf v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);

								get_st_point(v_buf.orig_face, v_buf.u, v_buf.v, mat, sts[i]);

								if( v_buf.orig_edge != NULL ){
									//Add the other faces for searching
									BMIter iter;
									BMFace *face;

									BM_ITER_ELEM (face, &iter, v_buf.orig_edge, BM_FACES_OF_EDGE){
										BLI_buffer_append(&faces, BMFace*, face);
										append_face(&faces, face);
									}
								} else {
									append_face(&faces, v_buf.orig_face);
								}
							} else {
								BMIter iter;
								BMFace *face;
								BMVert *temp_v = BLI_ghash_lookup(m_d->vert_hash, vert_arr[i]);

								mul_v2_m3v3(sts[i], mat, temp_v->co);

								BM_ITER_ELEM (face, &iter, temp_v, BM_FACES_OF_VERT) {
									append_face(&faces, face);
								}

							}
						}
					}

					//This search will spawn multiple faces.
					float edge1_mid[2];
					float edge2_mid[2];
					//We need to figure out which facing the radial edge is supposed to have so we can
					//try to fix it later if the inserted radial vert has the wrong facing
					bool is_B;
					BMVert *b_vert;

					if( mod_i(CC_idx+1,3) == CC2_idx ){
						b_vert = vert_arr[ mod_i(CC_idx-1,3) ];
					} else {
						b_vert = vert_arr[ mod_i(CC_idx+1,3) ];
					}

					is_B = calc_if_B_nor(m_d->cam_loc, b_vert->co, b_vert->no);

					interp_v2_v2v2(edge1_mid, sts[CC_idx], sts[ mod_i(CC_idx-1, 3) ], 0.5f);
					interp_v2_v2v2(edge2_mid, sts[CC_idx], sts[ mod_i(CC_idx+1, 3) ], 0.5f);

					//printf("Diff faces\n");
					mult_radi_search(&faces, mat, edge1_mid, edge2_mid, val_1, val_2, is_B, rad_plane_no, vert_arr[CC_idx], f, m_d);
					BLI_buffer_free(&faces);
				}

			}
		}
	}
}

static bool radial_C_vert(BMVert *v, MeshData *m_d){


	if( !(BM_elem_index_get(v) < m_d->radi_start_idx) ){
		//This is a radial vert
		return true;
	}

	if( is_vert_in_buffer( v, m_d->C_verts) ){
		return true;
	}

	return false;
}

static void radial_flip( MeshData *m_d ){

	bool done = false;

	int iter = 0;

	while( !done ){

	int flips = 0;
	int vert_i;
	for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
		BMEdge *e;
		BMIter iter_e;
		int edge_idx;
		int edge_i;
		BMVert *vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);
		int edge_count = BM_vert_edge_count(vert);
		BMEdge **edge_arr = BLI_array_alloca(edge_arr, edge_count);

        if( is_vert_in_buffer(vert, m_d->cusp_verts) ){
			//Do not flip cusp edges
			//continue;
		}

		BM_ITER_ELEM_INDEX (e, &iter_e, vert, BM_EDGES_OF_VERT, edge_idx) {
			edge_arr[edge_idx] = e;
		}

		for( edge_i = 0; edge_i < edge_count; edge_i++){
			BMVert *edge_vert;
			e = edge_arr[edge_i];

			if(e->v1 != vert){
				edge_vert = e->v1;
			} else {
				edge_vert = e->v2;
			}

			if( radial_C_vert( edge_vert, m_d ) ){
				//This is a radial or CC edge, do not try to flip it.
				continue;
			}

			//See if it's possible to rotate the edge at all
			// IE check for mesh border edges etc
			if( !BM_edge_rotate_check(e) ){
				//Not possible, bail
				//printf("Couldn't rotate edge!\n");
				continue;
			}

			{
				//Check if we can just do a simple rotation that doesn't create any bad faces
				BMLoop *l1, *l2;
				BM_edge_calc_rotate(e, true, &l1, &l2);

				// new_v1 = l1->v;
				// new_v2 = l2->v;

				if( BM_elem_index_get(l1->v) < m_d->radi_start_idx &&
					BM_elem_index_get(l2->v) < m_d->radi_start_idx ){
					//The flip will not increase the number of standard radial triangles
					//Do not flip it
					continue;
				}
				//Check if the flip creates any folds
				{
					float mat[3][3];
					float mat_coords[3][2], mat_new_pos[2];

					axis_dominant_v3_to_m3(mat, edge_vert->no);
					mul_v2_m3v3(mat_new_pos, mat, edge_vert->co);

					mul_v2_m3v3(mat_coords[0], mat, vert->co);
					mul_v2_m3v3(mat_coords[1], mat, l1->v->co);
					mul_v2_m3v3(mat_coords[2], mat, l2->v->co);

					if( !isect_point_tri_v2(mat_new_pos, mat_coords[0], mat_coords[1], mat_coords[2]) ){

						//We can simply rotate it!
						BM_edge_rotate(m_d->bm, e, true, 0);
						flips++;
						continue;
					}
				}

				//printf("Try to dissolve vert!\n");

				{
					int edge_count2 = BM_vert_edge_count(edge_vert);
					BMEdge *cur_e = e;
					BMEdge **edge_arr2 = BLI_array_alloca(edge_arr2, edge_count2);
					edge_idx = 0;

					BMEdge *rad1_edge = BM_edge_exists(l1->v, edge_vert);
					BMEdge *rad2_edge = BM_edge_exists(l2->v, edge_vert);

					BMLoop *first_loop = BM_face_vert_share_loop( cur_e->l->f, edge_vert);
					BMLoop *cur_loop = first_loop;

					if( edge_count2 - 3 < 1 ){
						continue;
					}

					if( rad1_edge == NULL || rad2_edge == NULL){
						printf("Couldn't find the edges connected to the radial vert! (dissolve vert)\n");
						continue;
					}

					if( cur_loop == NULL ){
						printf("Couldn't find face loop!\n");
					}

					while (((cur_loop = BM_vert_step_fan_loop(cur_loop, &cur_e)) != first_loop) && (cur_loop != NULL)) {
						if(cur_e == e || cur_e == rad1_edge || cur_e == rad2_edge){
							continue;
						}

						edge_arr2[edge_idx] = cur_e;
						edge_idx++;

					}

					if(cur_loop == NULL){
						continue;
					}

					for( edge_idx = 0; edge_idx < edge_count2 - 3; edge_idx++){
						BMLoop *loop1, *loop2;
						BM_edge_calc_rotate(edge_arr2[edge_idx], true, &loop1, &loop2);

						if( BM_edge_rotate_check_degenerate(edge_arr2[edge_idx], loop1, loop2) ){
							BM_edge_rotate(m_d->bm, edge_arr2[edge_idx], true, 0);
						} else {
							//Try to rotate from the other side instead
							//printf("Try from other side!\n");
							break;
						}
					}

					if( edge_idx != edge_count2 - 3 ){
						int op_idx = edge_count2 -4;
						bool failed_rotate = false;

						for(; edge_idx <= op_idx; op_idx--){

							BMLoop *loop1, *loop2;
							BM_edge_calc_rotate(edge_arr2[op_idx], true, &loop1, &loop2);

							if( edge_idx == op_idx ){
								//This is the last edge that is going to be rotated
								// check_degenerate is too strict in this case (in most cases the face area will be near zero)
								//TODO check for folds
								if(	BM_edge_rotate(m_d->bm, edge_arr2[op_idx], true, 0) != NULL ){
									continue;
								} else {
									failed_rotate = true;
									break;
								}
							}

							if( BM_edge_rotate_check_degenerate(edge_arr2[op_idx], loop1, loop2) ){
								BM_edge_rotate(m_d->bm, edge_arr2[op_idx], true, 0);
							} else {
								failed_rotate = true;
								break;
							}
						}

						if( failed_rotate ){
							//printf("Failed to flip and delete in radial edge flip!\n");
							continue;
						}

					}
				}

				if( BM_disk_dissolve(m_d->bm, edge_vert) ){
					flips++;
					//printf("Dissolved vert!\n");
				}

			}

		}
	}

	if(flips == 0 || iter > 5){
		done = true;
	}
	iter++;
	}
}

static int radial_extention( MeshData *m_d ){
	int exten = 0;

	for(int vert_i = 0; vert_i < m_d->radi_vert_buffer->count; vert_i++){
		Radi_vert r_vert = BLI_buffer_at(m_d->radi_vert_buffer, Radi_vert, vert_i);
		BMFace *face;
		BMIter iter;

		if( !(r_vert.extendable) ){
			continue;
		}

		bool b_f = r_vert.is_B;
		int prev_inco_faces = 0;
		float prev_face_diff = 0;
		BMEdge *flip_edge = NULL;
		bool flipped_edge = false;
		float cent_f[3], no[3];
		BM_ITER_ELEM (face, &iter, r_vert.vert, BM_FACES_OF_VERT) {
			BM_face_calc_center_median(face, cent_f); /* yiming:    suppose should be this fix */
			BM_face_calc_normal(face, no);

			float face_dir = get_facing_dir_nor(m_d->cam_loc, cent_f, no);
			if( b_f != (face_dir < 0) ){
				prev_face_diff += fabs(face_dir);
				prev_inco_faces++;

				BMIter iter_e;
				BMEdge *edge;

				BM_ITER_ELEM (edge, &iter_e, face, BM_EDGES_OF_FACE) {
					//Don't flip any edge directly connected to this radi vert
                    if( edge->v1 == r_vert.vert || edge->v2 == r_vert.vert ){
						continue;
					}
					float plane[4];
					plane_from_point_normal_v3(plane, r_vert.C_vert->co, r_vert.radi_plane_no);

                    //Make sure that the radi plane will cut the edge.
					//IE the edge points lies on opposite sides of the plane
					if( (dist_signed_to_plane_v3( edge->v1->co, plane ) < 0) ==
						(dist_signed_to_plane_v3( edge->v2->co, plane ) < 0) ){
						continue;
					}

					if( !BM_edge_rotate_check(edge) ){
						continue;
					}
					//Check so we don't try to flip any contour/radial edges
					//There has to be at least one C vert for it to be a true radial edge
					BLI_Buffer *cv;
					cv = m_d->C_verts;
					BMLoop *l1, *l2;
					BM_edge_calc_rotate(edge, true, &l1, &l2);
					//TODO actually check if we flip a radial edge!
					if( !is_vert_in_buffer(edge->v1, cv) && !is_vert_in_buffer(edge->v2, cv) &&
						!is_vert_in_buffer(l1->v, cv) && !is_vert_in_buffer(l2->v, cv)){
						float lambda;
						float temp[3];
						sub_v3_v3v3(temp, edge->v2->co, edge->v1->co);

						//Does the radial plane intersect the opposite edge?
						if( isect_ray_plane_v3(edge->v1->co, temp, plane, &lambda, true) ){
							flip_edge = edge;
							break;
						}
					}
				}
			}
		}

		if( prev_inco_faces == 0 ){
			continue;
		}

		if( flip_edge != NULL ){
			BMLoop *loop1, *loop2;
			BM_edge_calc_rotate(flip_edge, true, &loop1, &loop2);
			if( BM_edge_rotate_check_degenerate(flip_edge, loop1, loop2) ){
				flip_edge = BM_edge_rotate(m_d->bm, flip_edge, true, 0);
				flipped_edge = true;
			}
		} else {
			continue;
		}

		//Begin extenting the radi edge
		{
			float mat[3][3];
			float cur_pos_v2[2], c_pos_v2[2], r_pos_v2[2];
			float old_pos[3], old_no[3];
			float best_pos[3], best_no[3], best_diff_facing;

			int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
			int idx = BM_elem_index_get(r_vert.vert);
			BMFace *cur_face;
			Vert_buf v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, idx - orig_verts);

			bool found_better_pos = false;

            cur_face = v_buf.orig_face;

			copy_v3_v3( old_pos, r_vert.vert->co );
			copy_v3_v3( old_no, r_vert.vert->no );

			best_diff_facing = prev_face_diff;

			axis_dominant_v3_to_m3(mat, r_vert.vert->no);

			get_st_point(v_buf.orig_face, v_buf.u, v_buf.v, mat, r_pos_v2);

			//Get C_vert st pos
			get_vert_st_pos(m_d, mat, r_vert.C_vert, c_pos_v2);

			for( int i=1; i < 11; i++ ){
				float t = 1.0f + (float)i/10.0f;
				float P[3], du[3], dv[3];
				float uv_P[2];

				interp_v2_v2v2(cur_pos_v2, c_pos_v2, r_pos_v2, t);

				if( !point_inside_v2( mat, cur_pos_v2, cur_face ) ){
					BMFace *f;
					BMVert *v;
					BMIter iter_v, iter_f;
					bool found_face = false;

					BM_ITER_ELEM (v, &iter_v, cur_face, BM_VERTS_OF_FACE) {
						BM_ITER_ELEM (f, &iter_f, v, BM_FACES_OF_VERT) {
							if( point_inside_v2( mat, cur_pos_v2, f ) ){
								cur_face = f;
								found_face = true;
								break;
							}
						}
						if( found_face ){
							break;
						}
					}
					if( !found_face ){
						continue;
					}
				}

				get_uv_point( cur_face, uv_P, cur_pos_v2, mat );
				m_d->eval->evaluateLimit(m_d->eval, BM_elem_index_get(cur_face), uv_P[0], uv_P[1], P, du, dv);

				cross_v3_v3v3(no, du, dv);
				normalize_v3(no);

				copy_v3_v3(r_vert.vert->co, P);
				copy_v3_v3(r_vert.vert->no, no);

                //Did the nr of consistent triangles increase?
				{
                    int new_inco_faces = 0;
					float new_diff_facing = 0;

					BM_ITER_ELEM (face, &iter, r_vert.vert, BM_FACES_OF_VERT) {
						BM_face_calc_center_median(face, cent_f);
						BM_face_calc_normal(face, no);

						float face_dir = get_facing_dir_nor(m_d->cam_loc, cent_f, no);
						if( b_f != (face_dir < 0) ){
							new_diff_facing += fabs(face_dir);
							new_inco_faces++;
						}
					}

					if( new_inco_faces < prev_inco_faces ){
						found_better_pos = true;
						copy_v3_v3(best_pos, P);
						copy_v3_v3(best_no, no);
                        prev_inco_faces = new_inco_faces;
						best_diff_facing = new_diff_facing;
					} else if ( new_inco_faces == prev_inco_faces && new_diff_facing < best_diff_facing ){
						found_better_pos = true;
						copy_v3_v3(best_pos, P);
						copy_v3_v3(best_no, no);
						best_diff_facing = new_diff_facing;
					}
				}
			}

			if(found_better_pos){
				copy_v3_v3(r_vert.vert->co, best_pos);
				copy_v3_v3(r_vert.vert->no, best_no);
				exten++;
				continue;
			}

			//Move back to original position and flip back edge
			copy_v3_v3( r_vert.vert->co, old_pos );
			copy_v3_v3( r_vert.vert->no, old_no );

			if( flipped_edge ){
				BM_edge_rotate(m_d->bm, flip_edge, false, 0);
			}

		}
	}
	return exten;
}

static void null_opti_edge(MeshData *m_d, BMEdge *e, bool back_f, BLI_Buffer *inco_faces){
	BMFace *f;
	BMIter iter;
	BM_ITER_ELEM (f, &iter, e, BM_FACES_OF_EDGE) {
		float no[3], P[3];
		BM_face_calc_normal(f, no);
		BM_face_calc_center_median(f, P);
		bool found_face = false;
        bool face_good = (back_f == calc_if_B_nor(m_d->cam_loc, P, no));

		for(int i = 0; i < inco_faces->count; i++){
			IncoFace *inface = &BLI_buffer_at(inco_faces, IncoFace, i);
			if( inface->face != NULL && inface->face == f ){
				found_face = true;
				if( face_good ){
					inface->face = NULL;
				}
				break;
			}
		}
		if( !found_face && !face_good ){
			IncoFace new_inface;
			new_inface.face = f;
			new_inface.back_f = back_f;
			BLI_buffer_append(inco_faces, IncoFace, new_inface);
		}
	}
}

static void null_opti_vert(MeshData *m_d, BMVert *v, bool back_f, BLI_Buffer *inco_faces){
	BMFace *f;
	BMIter iter;
	BM_ITER_ELEM (f, &iter, v, BM_FACES_OF_VERT) {
		float no[3], P[3];
		BM_face_calc_normal(f, no);
		BM_face_calc_center_median(f, P);
		bool found_face = false;
        bool face_good = (back_f == calc_if_B_nor(m_d->cam_loc, P, no));

		for(int i = 0; i < inco_faces->count; i++){
			IncoFace *inface = &BLI_buffer_at(inco_faces, IncoFace, i);
			if( inface->face != NULL && inface->face == f ){
				found_face = true;
				if( face_good ){
					inface->face = NULL;
				}
				break;
			}
		}
		if( !found_face && !face_good ){
			IncoFace new_inface;
			new_inface.face = f;
			new_inface.back_f = back_f;
			BLI_buffer_append(inco_faces, IncoFace, new_inface);
		}
	}
}

//TODO change name, this doesn't create a fan copy anymore...
static void create_fan_copy(BMesh *bm_copy, BMFace *input_face, GHash *vhash, GHash *ehash, GHash *fhash){
	BMVert *v, *input_face_v;
	BMEdge *e;
	BMFace *f;
	BMIter iter_f, iter_e, iter_v, iter_input_f;

	BM_ITER_ELEM (input_face_v, &iter_input_f, input_face, BM_VERTS_OF_FACE) {

		BM_ITER_ELEM (f, &iter_f, input_face_v, BM_FACES_OF_VERT) {
			//Add verts to vhash
			BM_ITER_ELEM (v, &iter_v, f, BM_VERTS_OF_FACE) {
				if (BLI_ghash_lookup(vhash, v) == NULL){
					BMVert *new_vert;

					new_vert = BM_vert_create(bm_copy, v->co, NULL, BM_CREATE_SKIP_CD);
					BLI_ghash_insert(vhash, v, new_vert);
				}
			}
			//Add edges to ehash
			BM_ITER_ELEM (e, &iter_e, f, BM_EDGES_OF_FACE) {
				if (BLI_ghash_lookup(ehash, e) == NULL){
					BMEdge *new_edge;
					BMVert *v1, *v2;

					/* Lookup v1 and v2 */
					v1 = BLI_ghash_lookup(vhash, e->v1);
					v2 = BLI_ghash_lookup(vhash, e->v2);

					/* Create a new edge */
					new_edge = BM_edge_create(bm_copy, v1, v2, NULL, BM_CREATE_SKIP_CD);
					BLI_ghash_insert(ehash, e, new_edge);
				}
			}

			if (BLI_ghash_lookup(fhash, f) == NULL){
				//Create faces
				BMVert **vtar = BLI_array_alloca(vtar, f->len);
				BMEdge **edar = BLI_array_alloca(edar, f->len);
                BMFace *new_face;
				BMLoop *l_iter_src, *l_first_src;

				l_first_src = BM_FACE_FIRST_LOOP(f);

				/* lookup edge and vert order */
				l_iter_src = l_first_src;
				int i = 0;
				do {
					vtar[i] = BLI_ghash_lookup(vhash, l_iter_src->v);
					edar[i] = BLI_ghash_lookup(ehash, l_iter_src->e);
					i++;
				} while ((l_iter_src = l_iter_src->next) != l_first_src);

				new_face = BM_face_create(bm_copy, vtar, edar, f->len, f, BM_CREATE_SKIP_CD);

				BLI_ghash_insert(fhash, f, new_face);
			}
		}
	}
}

static int opti_vertex_wiggle( MeshData *m_d, BLI_Buffer *inco_faces ){
	int face_i;
	int wiggled_verts = 0;

	for(face_i = 0; face_i < inco_faces->count; face_i++){
		IncoFace *inface = &BLI_buffer_at(inco_faces, IncoFace, face_i);

		BMVert *vert;
		BMIter iter_v;

		if( inface->face == NULL ){
			//Already fixed this edge
			continue;
		}

		BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
			//Do not try to wiggle C verts
			if (is_vert_in_buffer( vert, m_d->C_verts)){
				continue;
			}

			{
				int face_count = BM_vert_face_count(vert);
				int face_idx, vert_idx;
				int nr_inco_faces = 0;
				float (*store_2d)[3][2] = BLI_array_alloca(store_2d, face_count);
				float *face_area = BLI_array_alloca(face_area, face_count);
				bool *radial_face = BLI_array_alloca(radial_face, face_count);
				float tot_face_area = 0;
				float tot_diff_facing = 0;
				float mat[3][3];

				bool done = false;

				axis_dominant_v3_to_m3(mat, vert->no);

				BMFace *f;
				BMVert *face_vert;
				BMIter iter_f, iter_f_v;
				BM_ITER_ELEM_INDEX (f, &iter_f, vert, BM_FACES_OF_VERT, face_idx) {
					radial_face[face_idx] = false;

					BM_ITER_ELEM_INDEX (face_vert, &iter_f_v, f, BM_VERTS_OF_FACE, vert_idx) {
						get_vert_st_pos(m_d, mat, face_vert, store_2d[face_idx][vert_idx]);
						if( is_vert_in_buffer(face_vert, m_d->C_verts) ) {
							radial_face[face_idx] = true;
						}
					}

					float no[3];
					float P[3];
					BM_face_calc_normal(f, no);
					BM_face_calc_center_median(f, P);

					float face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
					if( inface->back_f != (face_dir < 0) ){
						tot_diff_facing += fabs(face_dir);
						nr_inco_faces++;
					}
					face_area[face_idx] = area_tri_v2( store_2d[face_idx][0], store_2d[face_idx][1], store_2d[face_idx][2]);
					tot_face_area += face_area[face_idx];
				}

				BMEdge *edge;
				BMIter iter_e;
				float cent[3];
				zero_v3(cent);
				BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
					add_v3_v3(cent, BM_edge_other_vert(edge, vert)->co);
				}

				mul_v3_fl( cent, 1.0f / (float)BM_vert_edge_count(vert) );
				{
					float old_pos[3], best_pos[3];
					float best_diff_facing = tot_diff_facing;
					float best_dist = len_v3v3(cent, vert->co);
					int best_inco_faces = nr_inco_faces;

					BMFace *best_face;
					float best_uv[2];

					RNG *rng = BLI_rng_new(0);

					copy_v3_v3( old_pos, vert->co );
					copy_v3_v3( best_pos, vert->co );

					BLI_buffer_declare_static(BMFace*, orig_faces, BLI_BUFFER_NOP, 32);
					BMVert *orig_v = BLI_ghash_lookup(m_d->vert_hash, vert);
					if( orig_v != NULL ){
						BM_ITER_ELEM (f, &iter_f, orig_v, BM_FACES_OF_VERT) {
							append_face(&orig_faces, f);
						}
					} else {
						//This is a new vert, look the orig face it belongs to
						int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
						int v_idx = BM_elem_index_get(vert);
						Vert_buf v_buf;

						v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);

						BM_ITER_ELEM (face_vert, &iter_f_v, v_buf.orig_face, BM_VERTS_OF_FACE) {
							BM_ITER_ELEM (f, &iter_f, face_vert, BM_FACES_OF_VERT) {
								append_face(&orig_faces, f);
							}
						}
					}
					for( face_idx = 0; face_idx < face_count; face_idx++ ){
						int samples = (face_area[face_idx] / tot_face_area) * 100;

						for(int i=0; i < samples; i++ ){
							float cur_v2[2];
							float uv_P[2];
							BMFace *uv_face;

							// TODO check if the new point lies inside any of the new mesh faces
							BLI_rng_get_tri_sample_float_v2(rng, store_2d[face_idx][0], store_2d[face_idx][1], store_2d[face_idx][2], cur_v2);

							for (int f_idx = 0; f_idx < orig_faces.count; f_idx++) {
								uv_face = BLI_buffer_at(&orig_faces, BMFace*, f_idx);
								if( point_inside_v2( mat, cur_v2, uv_face ) ){
									float P[3], du[3], dv[3], no[3];

									get_uv_point( uv_face, uv_P, cur_v2, mat );
									m_d->eval->evaluateLimit(m_d->eval, BM_elem_index_get(uv_face), uv_P[0], uv_P[1], P, du, dv);

									cross_v3_v3v3(no, du, dv);
									normalize_v3(no);

									copy_v3_v3(vert->co, P);
									copy_v3_v3(vert->no, no);
									//No need to iterate over the remaining faces
									break;
								}
							}

							int new_inco_faces = 0;
							int f_idx;
							float new_dist = len_v3v3(cent, vert->co);
							float new_diff_facing = 0;
							bool skip_pos = false;

							BM_ITER_ELEM_INDEX (f, &iter_f, vert, BM_FACES_OF_VERT, f_idx) {
								float no[3];
								float P[3];
								BM_face_calc_normal(f, no);
								BM_face_calc_center_median(f, P);

								if( dot_v3v3(no, vert->no) < 0.0f ){
									//Punish flipped faces
									skip_pos = true;
									break;
								}

								float face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
								if( inface->back_f != (face_dir < 0) ){
									new_diff_facing += fabs(face_dir);
									new_inco_faces++;
									if( radial_face[f_idx] ){
										//Don't allow new pos to create inco radi faces
										skip_pos = true;
									}
								}
							}

							if( skip_pos ){
								continue;
							}

							if( new_inco_faces < best_inco_faces ){
								best_inco_faces = new_inco_faces;
								best_dist = new_dist;
								best_diff_facing = new_diff_facing;
								copy_v3_v3( best_pos, vert->co );
								done = true;

								best_face = uv_face;
								copy_v2_v2(best_uv, uv_P);
							} else if (new_inco_faces == 0 && new_dist < best_dist){
								best_dist = new_dist;
								copy_v3_v3( best_pos, vert->co );
								done = true;

								best_face = uv_face;
								copy_v2_v2(best_uv, uv_P);
							} else if (new_inco_faces == best_inco_faces && new_diff_facing < best_diff_facing){
								best_diff_facing = new_diff_facing;
								copy_v3_v3( best_pos, vert->co );
								done = true;

								best_face = uv_face;
								copy_v2_v2(best_uv, uv_P);
							}

						}
					}

					BLI_buffer_free(&orig_faces);
					BLI_rng_free(rng);
					copy_v3_v3(vert->co, best_pos);

					if( done ){
						add_shifted_vert( vert, best_face, best_uv, m_d );
						null_opti_vert(m_d, vert, inface->back_f, inco_faces);
						wiggled_verts++;
						//printf("Vertex wiggle\n");
						break;
					} else {
						//printf("Bad Vertex wiggle\n");
					}
				}


			}
		}
	}
	return wiggled_verts;
}

static void optimization( MeshData *m_d ){

	// 1. Radial edge extension
	{
		//How many radial edges did we extend this iteration?
		int exten = 0;
		do {
			exten = radial_extention( m_d );
			//printf("exten: %d\n", exten);
		} while (exten > 0);
	}

	BLI_buffer_declare_static(IncoFace, inco_faces, BLI_BUFFER_NOP, 32);
	//Find and save all inconsistent faces before we begin with the other optimization steps
	{
		BMVert *vert;
		for(int vert_i = 0; vert_i < m_d->radi_vert_buffer->count; vert_i++){
			Radi_vert r_vert = BLI_buffer_at(m_d->radi_vert_buffer, Radi_vert, vert_i);
			vert = r_vert.vert;

			{
				BMFace *face;
				BMIter iter_f;
				bool b_f = r_vert.is_B;
				float P[3], no[3];

				BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
					// and only check each face once
					// taken from BM_face_exists_overlap for marks
					if(BM_ELEM_API_FLAG_TEST(face, _FLAG_OVERLAP) != 0){
						//Already added this face to inco_faces
						continue;
					}

					{
						// This shouldn't be needed, but in case we manage to have inconsistent
						// faces that borders our contour line, don't mark it for adjustment.
						BMVert *v;
						BMIter iter;
						bool found_c_vert = false;

						BM_ITER_ELEM (v, &iter, face, BM_VERTS_OF_FACE) {
							if( is_vert_in_buffer(v, m_d->C_verts) ) {
								found_c_vert = true;
								break;
							}
						}

						if( found_c_vert ) {
							continue;
						}

					}
					BM_face_calc_center_median(face, P);
					BM_face_calc_normal(face, no);

					if( b_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
						IncoFace inface;
						inface.face = face;
						inface.back_f = b_f;
						BM_ELEM_API_FLAG_ENABLE(face, _FLAG_OVERLAP);
						BLI_buffer_append(&inco_faces, IncoFace, inface);
					}

				}
			}

		}

		//Clear _OVERLAP flag
		for(int face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);
			BM_ELEM_API_FLAG_DISABLE(inface->face, _FLAG_OVERLAP);
		}
	}

	// 2. Edge flipping
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMEdge *edge;
			BMIter iter_e;
			BMEdge *best_edge = NULL;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (edge, &iter_e, inface->face, BM_EDGES_OF_FACE) {

				if( !BM_edge_rotate_check(edge) ){
					continue;
				}

				BLI_Buffer *cv;
				cv = m_d->C_verts;
				BMLoop *l1, *l2;
				BM_edge_calc_rotate(edge, true, &l1, &l2);
				if( !is_vert_in_buffer(edge->v1, cv) && !is_vert_in_buffer(edge->v2, cv) &&
					!is_vert_in_buffer(l1->v, cv) && !is_vert_in_buffer(l2->v, cv)){
					//This is not a radial triangle edge, see if we can flip it

					if( !BM_edge_rotate_check_degenerate(edge, l1, l2) ){
						continue;
					}

					if( !BM_edge_rotate_check_beauty(edge, l1, l2) ){
						continue;
					}

					//Calculate nr of info faces of egde
					int nr_inco_faces = 0;
					float P[3], no[3];
					float face_dir;
					float cur_diff_facing = 0;
					BMFace *face;
					BMIter iter_f;
					BM_ITER_ELEM (face, &iter_f, edge, BM_FACES_OF_EDGE) {
						BM_face_calc_normal(face, no);
						BM_face_calc_center_median(face, P);

						//Calc facing of face
						face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);

						if( inface->back_f != (face_dir < 0) ){
							cur_diff_facing += fabs(face_dir);
							nr_inco_faces++;
						}
					}

					{
						int new_inco_faces = 0;
						float new_diff_facing = 0;

						//Taken from BM_edge_rotate_check_degenerate

						//Check if the flip will cause a potential fold
						float ed_dir_new[3], ed_dir_v1_new[3], ed_dir_v2_new[3];
						BMVert *v1, *v2;
						BMVert *v1_old, *v2_old;
						BMVert *v1_alt, *v2_alt;

						BM_edge_ordered_verts(edge, &v1_old, &v2_old);

						v1 = l1->v;
						v2 = l2->v;

						/* get the next vert along */
						v1_alt = BM_face_other_vert_loop(l1->f, v1_old, v1)->v;
						v2_alt = BM_face_other_vert_loop(l2->f, v2_old, v2)->v;

						sub_v3_v3v3(ed_dir_new, v1->co, v2->co);
						normalize_v3(ed_dir_new);

						sub_v3_v3v3(ed_dir_v1_new, v1->co, v1_alt->co);
						sub_v3_v3v3(ed_dir_v2_new, v2->co, v2_alt->co);
						normalize_v3(ed_dir_v1_new);
						normalize_v3(ed_dir_v2_new);

						//First flipped face
						cross_v3_v3v3(no, ed_dir_new, ed_dir_v1_new);
						normalize_v3(no);

						//Calc center mean of new face
						zero_v3(P);
						add_v3_v3( P, v1->co );
						add_v3_v3( P, l1->v->co );
						add_v3_v3( P, l2->v->co );

						mul_v3_fl( P, 1.0f / 3.0f );

						//Facing of first new flip face
						face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
						if( inface->back_f != (face_dir < 0) ){
							//This is not a good flip!
							//printf("Opti flip, first face not good\n");
							new_diff_facing += fabs(face_dir);
							new_inco_faces++;
						}

						//Second flipped face
						cross_v3_v3v3(no, ed_dir_new, ed_dir_v2_new);
						normalize_v3(no);

						//Calc center mean of new face
						zero_v3(P);
						add_v3_v3( P, v2->co );
						add_v3_v3( P, l1->v->co );
						add_v3_v3( P, l2->v->co );

						mul_v3_fl( P, 1.0f / 3.0f );

						//Facing of second new flip face
						face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
						if( inface->back_f != (face_dir < 0) ){
							//This is not a good flip!
							//printf("Opti flip, second face not good\n");
							new_diff_facing += fabs(face_dir);
							new_inco_faces++;
						}
						if (new_inco_faces < nr_inco_faces){
							best_edge = edge;
						} else if ( new_inco_faces == nr_inco_faces && new_diff_facing < cur_diff_facing ){
							best_edge = edge;
						}
					}
				}

			}
			if (best_edge != NULL){
				//printf("Opti filped an edge!\n");

				//When we rotate, the connected inco faces might be lost.
				//Remove them before flipping and add the new inco faces later (if any)
				BMFace *face;
				BMIter iter_f;
				BM_ITER_ELEM (face, &iter_f, best_edge, BM_FACES_OF_EDGE) {
					for(int i = 0; i < inco_faces.count; i++){
						IncoFace *edge_inface = &BLI_buffer_at(&inco_faces, IncoFace, i);
						if( edge_inface->face != NULL && edge_inface->face == face ){
							edge_inface->face = NULL;
						}
					}
				}

				best_edge = BM_edge_rotate(m_d->bm, best_edge, true, 0);

				null_opti_edge(m_d, best_edge, inface->back_f, &inco_faces);
			}

		}
	}

	// 3. Vertex wiggling in paramter space
	int wiggled_verts = 0;
	do {
		wiggled_verts = opti_vertex_wiggle(m_d, &inco_faces);
	} while( wiggled_verts > 0 );

	// 4. Edge Splitting
	{
		int face_i;
		int initial_inco_len = inco_faces.count; //Don't try to split inco faces added in this step

		for(face_i = 0; face_i < initial_inco_len; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BMesh *bm_fan_copy;
			bm_fan_copy = BM_mesh_create(&bm_mesh_allocsize_default, &((struct BMeshCreateParams){0}));

			GHash *vhash = BLI_ghash_ptr_new("opti edge split vhash");
			GHash *ehash = BLI_ghash_ptr_new("opti edge split ehash");
			GHash *fhash = BLI_ghash_ptr_new("opti face split fhash");

			create_fan_copy(bm_fan_copy, inface->face, vhash, ehash, fhash);

			//Calculate edge split positions
			Cusp split_data[3][10];

			int edge_idx;
			BMEdge *edge;
			BMIter iter_e;
			BMFace *f, *face_vert;
			BMIter iter_f, iter_f_v;
			BM_ITER_ELEM_INDEX (edge, &iter_e, inface->face, BM_EDGES_OF_FACE, edge_idx) {
				int c_radi_verts = 0;
				for(int vert_i = 0; vert_i < m_d->radi_vert_buffer->count; vert_i++){
					Radi_vert r_vert = BLI_buffer_at(m_d->radi_vert_buffer, Radi_vert, vert_i);
					if (r_vert.vert == edge->v1 || r_vert.vert == edge->v2){
						c_radi_verts++;
					}
				}
				//Do not try to split CC edges
				c_radi_verts += is_vert_in_buffer( edge->v1, m_d->C_verts) + is_vert_in_buffer( edge->v2, m_d->C_verts);

				if (c_radi_verts == 2){
					for( int i = 0; i < 10; i++ ){
						split_data[edge_idx][i].orig_face = NULL;
					}
					continue;
				}

				//Sample all edges at 10 points. These 10 points are the potential edge split points.
				//Do 10 samples but don't check end and start point
				float step = 1.0f/11.0f;
				float step_arr[] = { step*5.0f, step*6.0f, step*4.0f, step*7.0f, step*3.0f,
					step*8.0f, step*2.0f, step*9.0f, step*1.0f, step*10.0f };

				float mat[3][3];
				float edge_no[3];
				float start[2], end[2], cur_v2[2];

				interp_v3_v3v3(edge_no, edge->v1->no, edge->v2->no, 0.5f);

				axis_dominant_v3_to_m3(mat, edge_no);
				get_vert_st_pos(m_d, mat, edge->v1, start);
				get_vert_st_pos(m_d, mat, edge->v2, end);

				BLI_buffer_declare_static(BMFace*, orig_faces, BLI_BUFFER_NOP, 32);

				BMVert *orig_v = BLI_ghash_lookup(m_d->vert_hash, edge->v1);

				if (orig_v == NULL){
					orig_v = BLI_ghash_lookup(m_d->vert_hash, edge->v2);
				}

				if (orig_v != NULL){
					BM_ITER_ELEM (f, &iter_f, orig_v, BM_FACES_OF_VERT) {
						append_face(&orig_faces, f);
					}
				} else {
					//This is a new vert, look the orig face it belongs to
					int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
					int v_idx = BM_elem_index_get(edge->v1);
					Vert_buf v_buf;

					v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);

					BM_ITER_ELEM (face_vert, &iter_f_v, v_buf.orig_face, BM_VERTS_OF_FACE) {
						BM_ITER_ELEM (f, &iter_f, face_vert, BM_FACES_OF_VERT) {
							append_face(&orig_faces, f);
						}
					}
				}

				for( int i = 0; i < 10; i++ ){
					interp_v2_v2v2(cur_v2, start, end, step_arr[i]);

					split_data[edge_idx][i].orig_face = NULL;
					for (int f_idx = 0; f_idx < orig_faces.count; f_idx++) {
						f = BLI_buffer_at(&orig_faces, BMFace*, f_idx);
						if( point_inside_v2( mat, cur_v2, f ) ){
							float P[3], du[3], dv[3];
							float uv_P[2];

							get_uv_point( f, uv_P, cur_v2, mat );
							m_d->eval->evaluateLimit(m_d->eval, BM_elem_index_get(f), uv_P[0], uv_P[1], P, du, dv);

							copy_v3_v3(split_data[edge_idx][i].cusp_co, P);
							cross_v3_v3v3(split_data[edge_idx][i].cusp_no, du, dv);
							normalize_v3(split_data[edge_idx][i].cusp_no);

							split_data[edge_idx][i].u = uv_P[0];
							split_data[edge_idx][i].v = uv_P[1];

							split_data[edge_idx][i].orig_face = f;
							//No need to iterate over the remaining faces
							break;
						}
					}
				}
				BLI_buffer_free(&orig_faces);
			}

			//Calculate number of inco faces in the vicinity
			int nr_inco_faces = 0;
			BM_ITER_MESH (f, &iter_f, bm_fan_copy, BM_FACES_OF_MESH){
				float no[3];
				float P[3];
				BM_face_calc_normal(f, no);
				BM_face_calc_center_median(f, P);

				float face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
				if( inface->back_f != (face_dir < 0) ){
					nr_inco_faces++;
				}
			}
			{
				BMFace *copy_f = BLI_ghash_lookup(fhash, inface->face);

				int new_inco_faces;
				float new_diff_facing;
				int best_inco_faces = nr_inco_faces;
				int best_edge;
				Cusp best_edge_split;
				float best_dist = INFINITY;
				bool done = false;

				BM_ITER_ELEM_INDEX (edge, &iter_e, copy_f, BM_EDGES_OF_FACE, edge_idx) {
					BMesh *bm_temp = BM_mesh_copy(bm_fan_copy);
					BMEdge *temp_e = BM_edge_at_index_find(bm_temp, BM_elem_index_get(edge));
					BMVert *split_vert = NULL;

					for(int j = 0; j < 10; j++){
						if (split_data[edge_idx][j].orig_face == NULL){
							continue;
						}
						bool fold = false;
						new_diff_facing = 0;
						new_inco_faces = 0;

						if( split_vert == NULL){
							split_vert = split_edge_and_move_nor(NULL, bm_temp, temp_e, split_data[edge_idx][j].cusp_co, split_data[edge_idx][j].cusp_no);
						} else {
							copy_v3_v3(split_vert->co, split_data[edge_idx][j].cusp_co);
							copy_v3_v3(split_vert->no, split_data[edge_idx][j].cusp_no);
						}

						BM_ITER_MESH (f, &iter_f, bm_temp, BM_FACES_OF_MESH){
							float no[3];
							float P[3];
							BM_face_calc_normal(f, no);
							BM_face_calc_center_median(f, P);

							float face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
							if( inface->back_f != (face_dir < 0) ){
								new_diff_facing += fabs(face_dir);
								if( dot_v3v3(no, split_vert->no) < 0.0f ){
									//Punish flipped faces
									fold = true;
								}
								new_inco_faces++;
							}
						}

						if( fold ){
                        	continue;
						}

						BMEdge *tmp_e;
						BMIter tmp_iter;
						float edge_max = 0;
						float edge_min = INFINITY;

						BM_ITER_ELEM (tmp_e, &tmp_iter, split_vert, BM_EDGES_OF_VERT) {
							float edge_len = BM_edge_calc_length(tmp_e);
							if (edge_len > edge_max){
								edge_max = edge_len;
							}
							if (edge_len < edge_min){
								edge_min = edge_len;
							}
						}

						float new_dist = edge_max - edge_min;

						if( new_inco_faces < best_inco_faces ){
							best_edge = edge_idx;
							best_inco_faces = new_inco_faces;
							best_dist = new_dist;
							best_edge_split = split_data[edge_idx][j];
							done = true;
						} else if ( new_inco_faces == best_inco_faces && new_dist < best_dist){
							best_edge = edge_idx;
							best_inco_faces = new_inco_faces;
							best_dist = new_dist;
							best_edge_split = split_data[edge_idx][j];
							done = true;
						}
					}
					BM_mesh_free(bm_temp);
				}

				if( done ){
					BMVert *split_vert;
					BM_ITER_ELEM_INDEX (edge, &iter_e, inface->face, BM_EDGES_OF_FACE, edge_idx) {
						if( edge_idx == best_edge ){
							Vert_buf new_buf;
							//print_v3("best split pos", best_edge_split_pos);
							split_vert = split_edge_and_move_nor(m_d, m_d->bm, edge, best_edge_split.cusp_co, best_edge_split.cusp_no);

							new_buf.orig_face = best_edge_split.orig_face;
							new_buf.orig_edge = NULL;
							new_buf.u = best_edge_split.u;
							new_buf.v = best_edge_split.v;

							BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
							break;
						}
					}
					null_opti_vert(m_d, split_vert, inface->back_f, &inco_faces);
					//printf("Edge wiggle\n");
				} else {
					//printf("Bad edge wiggle\n");
				}
			}

			BLI_ghash_free(vhash, NULL, NULL);
			BLI_ghash_free(ehash, NULL, NULL);
			BLI_ghash_free(fhash, NULL, NULL);
			BM_mesh_free(bm_fan_copy);
		}
	}

	//Try to wiggle again after split
	do {
		wiggled_verts = opti_vertex_wiggle(m_d, &inco_faces);
	} while( wiggled_verts > 0 );

	// 2.a (Not in the paper) Vertex dissolve
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				//Do not try to wiggle C verts
				if (is_vert_in_buffer( vert, m_d->C_verts)){
					continue;
				}
				if( BM_elem_index_get(vert) < m_d->radi_start_idx ){
					//Not a radial vert, see if we can dissolve it to improve the consistency

					if( BM_vert_edge_count(vert) == 3 && BM_vert_face_count(vert) == 3 ){
						float vec1[3], vec2[3], P[3], no[3];
						BMLoop *l1, *l2;
						BMVert *v1, *v2;
						BM_edge_calc_rotate(vert->e, true, &l1, &l2);

						BM_edge_ordered_verts(vert->e, &v1, &v2);

						if( vert == v1 ){
							sub_v3_v3v3(vec1, v2->co, l1->v->co);
							sub_v3_v3v3(vec2, v2->co, l2->v->co);

							cross_v3_v3v3(no, vec2, vec1);
							normalize_v3(no);

							//Calc center mean of new face
							zero_v3(P);
							add_v3_v3( P, v2->co );
							add_v3_v3( P, l1->v->co );
							add_v3_v3( P, l2->v->co );

							mul_v3_fl( P, 1.0f / 3.0f );
						} else {
							sub_v3_v3v3(vec1, v1->co, l1->v->co);
							sub_v3_v3v3(vec2, v1->co, l2->v->co);

							cross_v3_v3v3(no, vec1, vec2);
							normalize_v3(no);

							//Calc center mean of new face
							zero_v3(P);
							add_v3_v3( P, v1->co );
							add_v3_v3( P, l1->v->co );
							add_v3_v3( P, l2->v->co );

							mul_v3_fl( P, 1.0f / 3.0f );
						}

						if(inface->back_f == calc_if_B_nor(m_d->cam_loc, P, no)){
							//printf("Opti dissolve\n");
							//TODO remove this or only when debug

							if(!BM_disk_dissolve(m_d->bm, vert)){
								printf("Failed to opti dissolve\n");
								continue;
							}
							inface->face = NULL;
							//Done with this face
							break;
						}

					}

				}
			}
		}
	}
	// 2.b (Not in the paper) Smooth vertex position
	// TODO perhaps move this to before wiggling in normal direction (IE after step 4)
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				if( !is_vert_in_buffer( vert, m_d->C_verts) ){
					//not a radial vert, try to smooth the vertex pos and see if the consistency improves

					float old_pos[3], co[3], co2[3];
					int i = 0;
					bool done = true;

					BMEdge *edge;
					BMIter iter_e;

					zero_v3(co);
					copy_v3_v3(old_pos, vert->co);


					BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
						copy_v3_v3(co2, BM_edge_other_vert(edge, vert)->co);
						add_v3_v3v3(co, co, co2);
						i += 1;
					}

					mul_v3_fl(co, 1.0f / (float)i);
					mid_v3_v3v3(co, co, vert->co);

					copy_v3_v3(vert->co, co);

					{
						BMFace *face;
						BMIter iter_f;

						BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
							float no[3];
							float P[3];
							BM_face_calc_normal(face, no);
							BM_face_calc_center_median(face, P);

							if( inface->back_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
								//Bad vertex move
								//printf("Bad vert smooth\n");
								//Move the vertex back to it's original position

								copy_v3_v3(vert->co, old_pos);
								done = false;
								break;
							}

						}

						if( done ){
							//Good vert smooth
							//printf("vert smooth\n");
							null_opti_vert(m_d, vert, inface->back_f, &inco_faces);
							break;
						}
					}

				}
			}
		}
	}
	// 5. Vertex wiggling in normal direction
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				//Do not try to wiggle C verts
				if (is_vert_in_buffer( vert, m_d->C_verts)){
					continue;
				}

				float tot_diff_facing = 0;
				int nr_inco_faces = 0;

				BMFace *f;
				BMIter iter_f;
				BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
					float no[3];
					float P[3];
					BM_face_calc_normal(f, no);
					BM_face_calc_center_median(f, P);

					float face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
					if( inface->back_f != (face_dir < 0) ){
						tot_diff_facing += fabs(face_dir);
						nr_inco_faces++;
					}
				}

				BMEdge *edge;
				BMIter iter_e;
				float old_pos[3];
				float new_pos[3];
				float len = 0.0f;
				int i = 0;
				bool done = false;

				copy_v3_v3(old_pos, vert->co);

				BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
					len += BM_edge_calc_length(edge);
					i++;
				}

				len = len / (float)i;
				len = len * 0.2f;

				for( i = 0; i < 100; i++ ){
					float co[3];
					float cur_len;

					float new_diff_facing = 0;
					int new_inco_faces = 0;

					copy_v3_v3(co, vert->no);

					if( i % 2 ){
						cur_len = len * (float)(i) / 100.0f;
					} else {
						cur_len = len * (float)(i+1) / -100.0f;
					}

					mul_v3_fl(co, cur_len);

					add_v3_v3v3(vert->co, old_pos, co);

					BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
						float no[3];
						float P[3];
						BM_face_calc_normal(f, no);
						BM_face_calc_center_median(f, P);

						float face_dir = get_facing_dir_nor(m_d->cam_loc, P, no);
						if( inface->back_f != (face_dir < 0) ){
							new_diff_facing += fabs(face_dir);
							new_inco_faces++;
						}
					}

					if( new_inco_faces < nr_inco_faces ){
						nr_inco_faces = new_inco_faces;
						tot_diff_facing = new_diff_facing;
						copy_v3_v3(new_pos, vert->co);
						done = true;
					} else if (new_inco_faces == nr_inco_faces && new_diff_facing < tot_diff_facing ){
						tot_diff_facing = new_diff_facing;
						copy_v3_v3(new_pos, vert->co);
						done = true;
					}

				}

				if( done ){
					copy_v3_v3(vert->co, new_pos);
					null_opti_vert(m_d, vert, inface->back_f, &inco_faces);
					//printf("Opti normal wiggle\n");
					if(inface->face == NULL){
						break;
					}
				} else {
					copy_v3_v3(vert->co, old_pos);
				}
			}
		}
	}

	//Debug color
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}
			inface->face->mat_nr = 4;
		}
	}

	//Cleanup
	BLI_buffer_free(&inco_faces);
}

static struct OpenSubdiv_Evaluator *create_osd_eval(BMesh *bm, Mesh *mesh){

	struct OpenSubdiv_Evaluator *osd_evaluator;
	int subdiv_levels = 1; //Need at least one subdiv level to compute the limit surface

	OpenSubdiv_Converter converter;
	OpenSubdiv_TopologyRefiner *topology_refiner;
    OpenSubdiv_TopologyRefinerSettings osd_settings;
    SubdivSettings sub_settings;

	osd_settings.level = subdiv_levels;
	osd_settings.is_adaptive = true;

    sub_settings.is_simple = false; //Use catmark subdiv
	sub_settings.fvar_linear_interpolation = SUBDIV_FVAR_LINEAR_INTERPOLATION_ALL;

	BKE_subdiv_converter_init_for_mesh(&converter, &sub_settings, mesh);

	topology_refiner = openSubdiv_createTopologyRefinerFromConverter(&converter, &osd_settings);

	//Free user_data from converter as it is not needed after creation
	BKE_subdiv_converter_free(&converter);

	osd_evaluator = openSubdiv_createEvaluatorFromTopologyRefiner(topology_refiner);

	if (osd_evaluator == NULL) {
		printf("OpenSubdiv initialization failed, should not happen.");
		return NULL;
	}

	BMIter iter_v;
	BMVert *vert;
	int no_of_verts = BM_mesh_elem_count(bm, BM_VERT);
	int j;
	float *vert_array = BLI_array_alloca(vert_array, 3 * no_of_verts);

	BM_ITER_MESH_INDEX (vert, &iter_v, bm, BM_VERTS_OF_MESH, j) {
		vert_array[3*j] = vert->co[0];
		vert_array[3*j + 1] = vert->co[1];
		vert_array[3*j + 2] = vert->co[2];
	}

	osd_evaluator->setCoarsePositions(osd_evaluator,
			vert_array,
			0,
			no_of_verts);

	osd_evaluator->refine(osd_evaluator);
	return osd_evaluator;
}

static void recalc_face_normals(BMesh *bm){
	BMIter iter;
	BMFace *f;

	BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH){
		BM_face_normal_update(f);
	}
}

static void debug_colorize(BMesh *bm, const float cam_loc[3]){
	BMIter iter;
	BMFace *f;
	float P[3];

	BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH){
		if( f->mat_nr == 4 ){
			continue;
		}
		BM_face_calc_center_median(f, P);
		if( calc_if_B_nor(cam_loc, P, f->no) ){
			f->mat_nr = 1;
		} else {
			f->mat_nr = 0;
		}
	}
}

static void debug_colorize_radi( MeshData *m_d ){
	BMIter iter, iter_e;
	BMFace *f;
	BMEdge *e;
	float P[3];

	int vert_i;
	for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
		BMVert *vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);

		BM_ITER_ELEM (e, &iter_e, vert, BM_EDGES_OF_VERT) {
			BMVert *edge_vert;

			if(e->v1 != vert){
				edge_vert = e->v1;
			} else {
				edge_vert = e->v2;
			}

			if( !(BM_elem_index_get(edge_vert) < m_d->radi_start_idx) ){
				//This is a radial/CC edge vert.
				BM_ITER_ELEM (f, &iter, e, BM_FACES_OF_EDGE){
					if( f->mat_nr == 4 ){
						continue;
					}
					BM_face_calc_center_median(f, P);
					if( calc_if_B_nor(m_d->cam_loc, P, f->no) ){
						f->mat_nr = 3;
					} else {
						f->mat_nr = 2;
					}
				}
			}
		}
	}
}

static void create_vert_mapping(MeshData *m_d){
    BMVert *new, *old; //Key, value
	BMIter new_iter, old_iter;

	unsigned int tot_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	m_d->vert_hash = BLI_ghash_int_new_ex("vert map", tot_verts);

	BM_CHECK_TYPE_ELEM_ASSIGN(new) = BM_iter_new(&new_iter, m_d->bm, BM_VERTS_OF_MESH, NULL);
	BM_CHECK_TYPE_ELEM_ASSIGN(old) = BM_iter_new(&old_iter, m_d->bm_orig, BM_VERTS_OF_MESH, NULL);

	//When this function is called, bm is just a copy of bm_orig. So we can now map the original verts together.
	for(int i = 0; i < tot_verts; i++){

		BLI_ghash_insert(m_d->vert_hash, new, old);

		BM_CHECK_TYPE_ELEM_ASSIGN(new) = BM_iter_step(&new_iter);
		BM_CHECK_TYPE_ELEM_ASSIGN(old) = BM_iter_step(&old_iter);
	}
}

/* bmesh only function */
static Mesh *mybmesh_do(Mesh *mesh, MyBMeshModifierData *mmd, float cam_loc[3])
{

	Mesh *result;
	BMesh *bm_orig, *bm;
	bool quad_mesh = true;

	struct OpenSubdiv_Evaluator *osd_eval;

	//TODO do not calc normals as we overwrite them later
	bm = BKE_mesh_to_bmesh_ex(
	         mesh,
	         &((struct BMeshCreateParams){0}),
	         &((struct BMeshFromMeshParams){.calc_face_normal = true,}));


	TIMEIT_START(quad_check);

	//TODO use this to check if we need to subdivide the mesh to get a quad mesh.
	{
		BMIter iter;
		BMFace *f, *f_next;
		int sides = 4;
		/* use the mutable iterator so we can remove data as its looped over */
		BM_ITER_MESH_MUTABLE (f, f_next, &iter, bm, BM_FACES_OF_MESH) {
			if (f->len != sides) {
				quad_mesh = false;
				break;
			}
		}
	}

	if (!quad_mesh){
		if( mmd->camera_ob != NULL ){
			debug_colorize(bm, cam_loc);
		}

		result = BKE_mesh_from_bmesh_nomain(bm, &((struct BMeshToMeshParams){0}));

		BM_mesh_free(bm);

		result->runtime.cd_dirty_vert |= CD_MASK_NORMAL;

		return result;
	}
	TIMEIT_END(quad_check);

	TIMEIT_START(osd_create);

	if (mmd->osd_eval == NULL){
		osd_eval = create_osd_eval(bm, mesh);
		mmd->osd_eval = osd_eval;
	} else {
		//TODO make a smarter cache system. Or even better, resuse the OSD data from a subsurf modifier!
		osd_eval = mmd->osd_eval;
	}

	TIMEIT_END(osd_create);

	// (6.1) Initialization
	TIMEIT_START(vert_limit);
	verts_to_limit(bm, osd_eval);
	TIMEIT_END(vert_limit);
	//Keep a copy of the quad mesh
	bm_orig = BM_mesh_copy(bm);

	if (mmd->flag & MOD_MYBMESH_TRIANG) {
		//TODO check if shortest diagonal is better
		BM_mesh_triangulate(bm, MOD_TRIANGULATE_QUAD_FIXED, MOD_TRIANGULATE_NGON_BEAUTY, false, NULL, NULL, NULL);
	}

	if( mmd->camera_ob == NULL){
		//Can't proceed without camera obj
		result = BKE_mesh_from_bmesh_nomain(bm, &((struct BMeshToMeshParams){0}));
		BM_mesh_free(bm);
		BM_mesh_free(bm_orig);
		return result;
	}
	{
		BLI_buffer_declare_static(Vert_buf, new_vert_buffer, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(Vert_buf, shifted_verts, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(Cusp, cusp_edges, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(BMVert*, C_verts, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(BMVert*, cusp_verts, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(Radi_vert, radi_vert_buffer, BLI_BUFFER_NOP, 32);

		MeshData mesh_data;

		mesh_data.bm = bm;
		mesh_data.bm_orig = bm_orig;
		mesh_data.new_vert_idx = BM_mesh_elem_count(bm, BM_VERT);

		copy_v3_v3( mesh_data.cam_loc, cam_loc );

		mesh_data.new_vert_buffer = &new_vert_buffer;
		mesh_data.shifted_verts = &shifted_verts;
		mesh_data.cusp_edges = &cusp_edges;
		mesh_data.C_verts = &C_verts;
		mesh_data.cusp_verts = &cusp_verts;
		mesh_data.radi_vert_buffer = &radi_vert_buffer;
		mesh_data.eval = osd_eval;

        create_vert_mapping(&mesh_data);

		if (mmd->flag & MOD_MYBMESH_FF_SPLIT) {
			TIMEIT_START(split_bb_ff);
			split_BB_FF_edges_thread_start(&mesh_data);
			//split_BB_FF_edges(&mesh_data);
			TIMEIT_END(split_bb_ff);
			//printf("New verts: %d\n", new_vert_buffer.count);
		}
		// (6.2) Contour Insertion

		if (mmd->flag & MOD_MYBMESH_CUSP_D) {
			TIMEIT_START(cusp_detect);
			cusp_detection(&mesh_data);
			TIMEIT_END(cusp_detect);
		}

		if (mmd->flag & MOD_MYBMESH_FB_SPLIT) {
			TIMEIT_START(contour_insert);
			contour_insertion(&mesh_data);
			TIMEIT_END(contour_insert);
		}
		if (mmd->flag & MOD_MYBMESH_CUSP_I) {
			TIMEIT_START(cusp_insert);
			cusp_insertion(&mesh_data);
			TIMEIT_END(cusp_insert);
		}

		// (6.3) Radialization

		mesh_data.radi_start_idx = BM_mesh_elem_count(bm, BM_VERT);

		if (mmd->flag & MOD_MYBMESH_RAD_I){
			TIMEIT_START(radi_insert);
			radial_insertion(&mesh_data);
			TIMEIT_END(radi_insert);
		}

		if (mmd->flag & MOD_MYBMESH_RAD_FLIP){
			TIMEIT_START(radi_flip);
			radial_flip(&mesh_data);
			TIMEIT_END(radi_flip);
		}

		//Recalculate normals
		recalc_face_normals(bm);

		// (6.4) Optimization
		if (mmd->flag & MOD_MYBMESH_OPTI){
			TIMEIT_START(opti);
			optimization(&mesh_data);
			TIMEIT_END(opti);
			//Recalculate normals
			recalc_face_normals(bm);
		}

		TIMEIT_START(debug_color);
		debug_colorize(bm, cam_loc);
		debug_colorize_radi(&mesh_data);
		TIMEIT_END(debug_color);
		BLI_ghash_free(mesh_data.vert_hash, NULL, NULL);
		BLI_buffer_free(&new_vert_buffer);
		BLI_buffer_free(&shifted_verts);
		BLI_buffer_free(&cusp_edges);
		BLI_buffer_free(&C_verts);
		BLI_buffer_free(&cusp_verts);
		BLI_buffer_free(&radi_vert_buffer);
	}
	result = BKE_mesh_from_bmesh_nomain(bm, &((struct BMeshToMeshParams){0}));

	BM_mesh_free(bm);
	BM_mesh_free(bm_orig);

	result->runtime.cd_dirty_vert |= CD_MASK_NORMAL;

	return result;
}

/* MyBMesh */
static void initData(ModifierData *md)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *) md;

	mmd->camera_ob = NULL;
	mmd->osd_eval = NULL;
}

static void freeData(ModifierData *md)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *) md;
	if (mmd->osd_eval != NULL){
	openSubdiv_deleteEvaluator(mmd->osd_eval);
	}
}

static Mesh *applyModifier(ModifierData *md, const ModifierEvalContext *ctx,
								Mesh *mesh)
{
	Mesh *result;
	float cam_loc[3];
	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;

	if(mmd->camera_ob){
		float trans_mat[4][4];

		//Create tranformation matrix to get relative coordinates of the camera obj.
		invert_m4_m4(trans_mat, ctx->object->obmat);

		copy_v3_v3(cam_loc, mmd->camera_ob->loc);
		/*
		printf("Cam loc:\n");
		printf("1: %f\n", cam_loc[0]);
		printf("2: %f\n", cam_loc[1]);
		printf("3: %f\n", cam_loc[2]);
		*/
		//convert camera origin from world coord to the modifier obj local coords
		mul_m4_v3(trans_mat, cam_loc);
		/*
		printf("Cam loc 2:\n");
		printf("1: %f\n", cam_loc[0]);
		printf("2: %f\n", cam_loc[1]);
		printf("3: %f\n", cam_loc[2]);
		*/
	}

	if (!(result = mybmesh_do(mesh, mmd, cam_loc))) {
		return mesh;
	}

	return result;
}

static void foreachObjectLink(
        ModifierData *md, Object *ob,
        ObjectWalkFunc walk, void *userData)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;

	walk(userData, ob, &mmd->camera_ob, IDWALK_NOP);
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;
	if (mmd->camera_ob != NULL) {
		DEG_add_object_relation(ctx->node, mmd->camera_ob, DEG_OB_COMP_TRANSFORM, "MyBmesh Modifier");
	}
	DEG_add_object_relation(ctx->node, ctx->object, DEG_OB_COMP_TRANSFORM, "MyBmesh Modifier");
}

static bool dependsOnNormals(ModifierData *UNUSED(md))
{
	//TODO it does depend on normals. return true here?
	return false;
}

ModifierTypeInfo modifierType_MyBMesh = {
	/* name */              "MyBMesh",
	/* structName */        "MyBMeshModifierData",
	/* structSize */        sizeof(MyBMeshModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	/* flags cont */        eModifierTypeFlag_SupportsEditmode |
	/* flags cont */        eModifierTypeFlag_EnableInEditmode,

	/* copyData */          modifier_copyData_generic,

	/* deformVerts_DM */    NULL,
	/* deformMatrices_DM */ NULL,
	/* deformVertsEM_DM */  NULL,
	/* deformMatricesEM_DM*/NULL,
	/* applyModifier_DM */  NULL,

	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,

	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          freeData,
	/* isDisabled */        NULL,
	/* updateDepsgraph */   updateDepsgraph,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  dependsOnNormals,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
