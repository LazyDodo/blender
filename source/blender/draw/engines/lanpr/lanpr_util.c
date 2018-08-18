/*

   Ported from NUL4.0

   Author(s):WuYiming - xp8110@outlook.com

 */
#define _CRT_SEQURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include "lanpr_util.h"
#include "lanpr_all.h"
#include <math.h>


//===================================================================[slt]


void list_handle_empty(ListBase *h) {
	h->first = h->last = 0;
}

void *lst_get_top(ListBase *Handle){
	return Handle->first;
};

int list_remove_segment(ListBase *Handle, nListItem *Begin, nListItem *End) {
	if (!Begin->pPrev)
		Handle->first = End->pNext;
	else
		((nListItem *)Begin->pPrev)->pNext = End->pNext;

	if (!End->pNext)
		Handle->last = Begin->pPrev;
	else
		((nListItem *)End->pNext)->pPrev = Begin->pPrev;

	End->pNext = Begin->pPrev = 0;

	return 1;
};
void list_insert_item_before(ListBase *Handle, nListItem *toIns, nListItem *pivot){
	if (!pivot) { BLI_addhead(Handle, toIns); return; }

	if (pivot->pPrev) {
		((nListItem *)pivot->pPrev)->pNext = toIns;
		toIns->pPrev = pivot->pPrev;
	}
	else {
		Handle->first = toIns;
	}

	toIns->pNext = pivot;
	pivot->pPrev = toIns;
};
void list_insert_item_after(ListBase *Handle, nListItem *toIns, nListItem *pivot) {
	if (!pivot) { BLI_addtail(Handle, toIns); return; }

	if (pivot->pNext) {
		((nListItem *)pivot->pNext)->pPrev = toIns;
		toIns->pNext = pivot->pNext;
	}
	else {
		Handle->last = toIns;
	}

	toIns->pPrev = pivot;
	pivot->pNext = toIns;
}

void *list_append_pointer_only(ListBase *h, void *p) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = CreateNew(nListItemPointer);
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_append_pointer_sized_only(ListBase *h, void *p, int size) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = calloc(1, size);
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_push_pointer_only(ListBase *h, void *p) {
	nListItemPointer *lip = 0;
	if (!h) return 0;
	lip = CreateNew(nListItemPointer);
	lip->p = p;
	BLI_addhead(h, lip);
	return lip;
}
void *list_push_pointer_sized_only(ListBase *h, void *p, int size) {
	nListItemPointer *lip = 0;
	if (!h) return 0;
	lip = calloc(1, size);
	lip->p = p;
	BLI_addhead(h, lip);
	return lip;
}

void *list_append_pointer(ListBase *h, void *p) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = memAquireOnly(sizeof(nListItemPointer));
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_append_pointer_sized(ListBase *h, void *p, int size) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = memAquireOnly(size);
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_push_pointer(ListBase *h, void *p) {
	nListItemPointer *lip = 0;
	if (!h) return 0;
	lip = memAquireOnly(sizeof(nListItemPointer));
	lip->p = p;
	BLI_addhead(h, lip);
	return lip;
}
void *list_push_pointer_sized(ListBase *h, void *p, int size) {
	nListItemPointer *lip = 0;
	if (!h) return 0;
	lip = memAquireOnly(size);
	lip->p = p;
	BLI_addhead(h, lip);
	return lip;
}

void *list_append_pointer_static(ListBase *h, nStaticMemoryPool *smp, void *p) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = mem_static_aquire(smp, sizeof(nListItemPointer));
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_append_pointer_static_sized(ListBase *h, nStaticMemoryPool *smp, void *p, int size) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = mem_static_aquire(smp, size);
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_push_pointer_static(ListBase *h, nStaticMemoryPool *smp, void *p) {
	nListItemPointer *lip = 0;
	if (!h) return 0;
	lip = mem_static_aquire(smp, sizeof(nListItemPointer));
	lip->p = p;
	BLI_addhead(h, lip);
	return lip;
}
void *list_push_pointer_static_sized(ListBase *h, nStaticMemoryPool *smp, void *p, int size) {
	nListItemPointer *lip = 0;
	if (!h) return 0;
	lip = mem_static_aquire(smp, size);
	lip->p = p;
	BLI_addhead(h, lip);
	return lip;
}

void *list_pop_pointer_only(ListBase *h) {
	nListItemPointer *lip;
	void *rev = 0;
	if (!h) return 0;
	lip = BLI_pophead(h); BLI_remlink(h, h->first);
	rev = lip ? lip->p : 0;
	FreeMem(lip);
	return rev;
}
void list_remove_pointer_item_only(ListBase *h, nListItemPointer *lip) {
	BLI_remlink(h, (void *)lip);
	FreeMem(lip);
}
void list_remove_pointer_only(ListBase *h, void *p) {
	nListItemPointer *i;
	for (i = h->first; i; i = i->pNext) {
		if (i->p == p) {
			list_remove_pointer_item(h, i);
			break;
		}
	}
}
void list_clear_pointer_only(ListBase *h) {
	while (h && h->first) {
		list_pop_pointer(h);
	}
}
void list_generate_pointer_list_only(ListBase *from1, ListBase *from2, ListBase *to) {
	nListItemPointer *lip = from2 ? from2->last : 0;

	while (lip) {
		list_push_pointer(to, lip->p);
		lip = lip->pPrev;
	}

	lip = from1 ? from1->last : 0;

	while (lip) {
		list_push_pointer(to, lip->p);
		lip = lip->pPrev;
	}
}

void *list_pop_pointer(ListBase *h) {
	nListItemPointer *lip;
	void *rev = 0;
	if (!h) return 0;
	lip = BLI_pophead(h);
	rev = lip ? lip->p : 0;
	mem_free(lip);
	return rev;
}
void list_remove_pointer_item(ListBase *h, nListItemPointer *lip) {
	BLI_remlink(h, (void *)lip);
	mem_free(lip);
}
void list_remove_pointer(ListBase *h, void *p) {
	nListItemPointer *i;
	for (i = h->first; i; i = i->pNext) {
		if (i->p == p) {
			list_remove_pointer_item(h, i);
			break;
		}
	}
}
void list_clear_pointer(ListBase *h) {
	nListItemPointer *i;
	while (h && h->first) {
		list_pop_pointer(h);
	}
}
void list_generate_pointer_list(ListBase *from1, ListBase *from2, ListBase *to) {
	nListItemPointer *lip = from2 ? from2->last : 0;

	while (lip) {
		list_push_pointer(to, lip->p);
		lip = lip->pPrev;
	}

	lip = from1 ? from1->last : 0;

	while (lip) {
		list_push_pointer(to, lip->p);
		lip = lip->pPrev;
	}
}


void *list_append_pointer_static_pool(nStaticMemoryPool *mph, ListBase *h, void *p) {
	nListItemPointer *lip;
	if (!h) return 0;
	lip = mem_static_aquire(mph, sizeof(nListItemPointer));
	lip->p = p;
	BLI_addtail(h, lip);
	return lip;
}
void *list_pop_pointer_no_free(ListBase *h) {
	nListItemPointer *lip;
	void *rev = 0;
	if (!h) return 0;
	lip = BLI_pophead(h);
	rev = lip ? lip->p : 0;
	return rev;
}
void list_remove_pointer_item_no_free(ListBase *h, nListItemPointer *lip) {
	BLI_remlink(h, (void *)lip);
}


void list_copy_handle(ListBase *target, ListBase *src){
	target->first = src->first;
	target->last = src->last;
};
void list_clear_handle(ListBase *h){
	h->first = 0;
	h->last = 0;
}
void list_clear_prev_next(nListItem *li){
	li->pNext = 0;
	li->pPrev = 0;
}

void list_move_up(ListBase *h, nListItem *li) {
	void *pprev = li->pPrev ? ((nListItem *)li->pPrev)->pPrev : 0;
	if (!h || !li)
		return;
	if (li == h->first) return;
	else {
		if (li == h->last) h->last = li->pPrev;
		((nListItem *)li->pPrev)->pNext = li->pNext;
		((nListItem *)li->pPrev)->pPrev = li;
		if (li->pNext) ((nListItem *)li->pNext)->pPrev = li->pPrev;
		li->pNext = li->pPrev;
		li->pPrev = pprev;
		if (pprev) ((nListItem *)pprev)->pNext = li;
	}
	if (!li->pPrev) h->first = li;
}
void list_move_down(ListBase *h, nListItem *li) {
	void *ppnext = li->pNext ? ((nListItem *)li->pNext)->pNext : 0;
	if (!h || !li)
		return;
	if (li == h->last) return;
	else {
		if (li == h->first) h->first = li->pNext;
		((nListItem *)li->pNext)->pPrev = li->pPrev;
		((nListItem *)li->pNext)->pNext = li;
		if (li->pPrev) ((nListItem *)li->pPrev)->pNext = li->pNext;
		li->pPrev = li->pNext;
		li->pNext = ppnext;
		if (ppnext) ((nListItem *)ppnext)->pPrev = li;
	}
	if (!li->pNext) h->last = li;
}

void mem_init_pool(nMemoryPool *mph, int NodeSize) {
	int Count = 4096;

	if (mph->NodeSize) return;

	mph->NodeSize = NodeSize;
	mph->CountPerPool = Count;

	return;
}
void mem_init_pool_small(nMemoryPool *mph, int NodeSize) {
	int Count = 16;

	if (mph->NodeSize) return;

	mph->NodeSize = NodeSize;
	mph->CountPerPool = Count;

	return;
}
inline nMemoryPoolPart *mem_new_pool_part(nMemoryPool *mph) {

	if (!mph->NodeSize) return 0;

	int RealNodeSize = mph->NodeSize + sizeof(nMemoryPoolNode);
	int NodeCount = mph->CountPerPool;
	int TotalSize = sizeof(nMemoryPoolPart) + NodeCount * RealNodeSize;
	int i;
	nMemoryPoolPart *mp = calloc(1, TotalSize);

	void *BeginMem = ((BYTE *)mp) + sizeof(nMemoryPoolPart);

	mp->PoolRoot = mph;

	nMemoryPoolNode *mpn;

	mp->FreeMemoryNodes.first = mp->FreeMemoryNodes.last = BeginMem;

	mpn = (void *)((BYTE *)BeginMem);
	mpn->InPool = mp;

	for (i = 1; i < NodeCount; i++) {
		mpn = (void *)(((BYTE *)BeginMem) + RealNodeSize * i);
		mpn->InPool = mp;
		BLI_addtail(&mp->FreeMemoryNodes, mpn);
	}
	BLI_addhead(&mph->Pools, mp);

	return mp;
}
void mem_free(void *Data) {
	if (!Data) return;
	nMemoryPoolNode *mpn = (void *)(((BYTE *)Data) - sizeof(nMemoryPoolNode));
	nMemoryPoolPart *mp = mpn->InPool;
	nMemoryPool *mph = mp->PoolRoot;
	BLI_remlink(&mp->MemoryNodes, (void *)mpn);
	BLI_addtail(&mp->FreeMemoryNodes, (void *)mpn);
	memset(Data, 0, mph->NodeSize);
	if (!mp->MemoryNodes.first) {
		BLI_remlink(&mph->Pools, (void *)mp);
		FreeMem(mp);
	}
	//if (!mph->Pools.first) {
	//	mph->CountPerPool = 0;
	//	mph->NodeSize = 0;
	//}
}
void mem_destroy_pool(nMemoryPool *Handle) {
	nMemoryPool *mp;
	while ((mp = BLI_pophead(&Handle->Pools))) {
		FreeMem(mp);
	}
}

nStaticMemoryPoolNode *mem_new_static_pool(nStaticMemoryPool *smp) {
	nStaticMemoryPoolNode *smpn = MEM_callocN(NUL_MEMORY_POOL_128MB, "mempool");
	smpn->UsedByte = sizeof(nStaticMemoryPoolNode);
	BLI_addhead(&smp->Pools, smpn);
	return smpn;
}
void *mem_static_aquire(nStaticMemoryPool *smp, int size) {
	nStaticMemoryPoolNode *smpn = smp->Pools.first;
	void *ret;

	if (!smpn || (smpn->UsedByte + size) > NUL_MEMORY_POOL_128MB)
		smpn = mem_new_static_pool(smp);

	ret = ((BYTE *)smpn) + smpn->UsedByte;

	smpn->UsedByte += size;

	return ret;
}
void *mem_static_aquire_thread(nStaticMemoryPool *smp, int size) {
	nStaticMemoryPoolNode *smpn = smp->Pools.first;
	void *ret;

	BLI_spin_lock(&smp->csMem);

	if (!smpn || (smpn->UsedByte + size) > NUL_MEMORY_POOL_128MB)
		smpn = mem_new_static_pool(smp);

	ret = ((BYTE *)smpn) + smpn->UsedByte;

	smpn->UsedByte += size;

	BLI_spin_unlock(&smp->csMem);

	return ret;
}
void *mem_static_destroy(nStaticMemoryPool *smp) {
	nStaticMemoryPoolNode *smpn;
	void *ret = 0;

	while (smpn = BLI_pophead(&smp->Pools)) {
		FreeMem(smpn);
	}

	smp->EachSize = 0;

	return ret;
}

//=======================================================================[str]


void tmat_load_identity_44d(tnsMatrix44d m) {
	memset(m, 0, sizeof(tnsMatrix44d));
	m[0] = 1.0f;
	m[5] = 1.0f;
	m[10] = 1.0f;
	m[15] = 1.0f;
};

real tmat_dist_idv2(real x1, real y1, real x2, real y2) {
	real x = x2 - x1, y = y2 - y1;
	return sqrt((x * x + y * y));
}
real tmat_dist_3dv(tnsVector3d l, tnsVector3d r) {
	real x = r[0] - l[0];
	real y = r[1] - l[1];
	real z = r[2] - l[2];
	return sqrt(x * x + y * y + z * z);
}
real tmat_dist_2dv(tnsVector2d l, tnsVector2d r) {
	real x = r[0] - l[0];
	real y = r[1] - l[1];
	return sqrt(x * x + y * y);
}

real tmat_length_3d(tnsVector3d l) {
	return (sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]));
}
real tmat_length_2d(tnsVector3d l) {
	return (sqrt(l[0] * l[0] + l[1] * l[1]));
}
void tmat_normalize_3d(tnsVector3d result, tnsVector3d l) {
	real r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
	result[2] = l[2] / r;
}
void tmat_normalize_3f(tnsVector3f result, tnsVector3f l) {
	float r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
	result[2] = l[2] / r;
}
void tmat_normalize_2d(tnsVector2d result, tnsVector2d l) {
	real r = sqrt(l[0] * l[0] + l[1] * l[1]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
}
void tmat_normalize_self_3d(tnsVector3d result) {
	real r = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
	result[0] /= r;
	result[1] /= r;
	result[2] /= r;
}
real tmat_dot_3d(tnsVector3d l, tnsVector3d r, int normalize) {
	tnsVector3d ln, rn;
	if (normalize) {
		tmat_normalize_3d(ln, l); tmat_normalize_3d(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1] + ln[2] * rn[2]);
	}
	return (l[0] * r[0] + l[1] * r[1] + l[2] * r[2]);
}
real tmat_dot_3df(tnsVector3d l, tnsVector3f r, int normalize) {
	tnsVector3d ln; tnsVector3f rn;
	if (normalize) {
		tmat_normalize_3d(ln, l); tmat_normalize_3f(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1] + ln[2] * rn[2]);
	}
	return (l[0] * r[0] + l[1] * r[1] + l[2] * r[2]);
}
real tmat_dot_2d(tnsVector2d l, tnsVector2d r, int normalize) {
	tnsVector3d ln, rn;
	if (normalize) {
		tmat_normalize_2d(ln, l); tmat_normalize_2d(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1]);
	}
	return (l[0] * r[0] + l[1] * r[1]);
}
real tmat_vector_cross_3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[1] * r[2] - l[2] * r[1];
	result[1] = l[2] * r[0] - l[0] * r[2];
	result[2] = l[0] * r[1] - l[1] * r[0];
	return tmat_length_3d(result);
}
void tmat_vector_cross_only_3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[1] * r[2] - l[2] * r[1];
	result[1] = l[2] * r[0] - l[0] * r[2];
	result[2] = l[0] * r[1] - l[1] * r[0];
}
real tmat_angle_rad_3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference) {
	if (PositiveReference) {
		tnsVector3d res;
		tmat_vector_cross_3d(res, from, to);
		if (tmat_dot_3d(res, PositiveReference, 1) > 0)
			return acosf(tmat_dot_3d(from, to, 1));
		else
			return -acosf(tmat_dot_3d(from, to, 1));
	}
	return acosf(tmat_dot_3d(from, to, 1));
}
void tmat_apply_rotation_33d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	result[0] = mat[0] * v[0] + mat[1] * v[1] + mat[2] * v[2];
	result[1] = mat[3] * v[0] + mat[4] * v[1] + mat[5] * v[2];
	result[2] = mat[6] * v[0] + mat[7] * v[1] + mat[8] * v[2];
}
void tmat_apply_rotation_43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	result[0] = mat[0] * v[0] + mat[1] * v[1] + mat[2] * v[2];
	result[1] = mat[4] * v[0] + mat[5] * v[1] + mat[6] * v[2];
	result[2] = mat[8] * v[0] + mat[9] * v[1] + mat[10] * v[2];
}
void tmat_apply_transform_43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	w = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
	result[0] /= w;
	result[1] /= w;
	result[2] /= w;
}
void tmat_apply_normal_transform_43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
}
void tmat_apply_normal_transform_43df(tnsVector3d result, tnsMatrix44d mat, tnsVector3f v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
}
void tmat_apply_transform_44d(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
}
void tmat_apply_transform_43dfND(tnsVector4d result, tnsMatrix44d mat, tnsVector3f v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
}
void tmat_apply_transform_43df(tnsVector4d result, tnsMatrix44d mat, tnsVector3f v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	real w = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
	//result[0] /= w;
	//result[1] /= w;
	//result[2] /= w;
}
void tmat_apply_transform_44dTrue(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * v[3];
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * v[3];
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * v[3];
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * v[3];
}

void tmat_remove_translation_44d(tnsMatrix44d result, tnsMatrix44d mat) {
	tmat_load_identity_44d(result);
	result[0] = mat[0];
	result[1] = mat[1];
	result[2] = mat[2];
	result[4] = mat[4];
	result[5] = mat[5];
	result[6] = mat[6];
	result[8] = mat[8];
	result[9] = mat[9];
	result[10] = mat[10];
}
void tmat_clear_translation_44d(tnsMatrix44d mat) {
	mat[3] = 0;
	mat[7] = 0;
	mat[11] = 0;
}


void tmat_extract_xyz_euler_44d(tnsMatrix44d mat, real *x_result, real *y_result, real *z_result) {
	real xRot, yRot, zRot;

	if (mat[2] < 1) {
		if (mat[2] > -1) {
			yRot = asin(mat[2]);
			xRot = atan2(-mat[6], mat[10]);
			zRot = atan2(-mat[1], mat[0]);
		}
		else {
			yRot = -TNS_PI / 2;
			xRot = -atan2(-mat[4], mat[5]);
			zRot = 0;
		}
	}
	else {
		yRot = TNS_PI / 2;
		xRot = atan2(-mat[4], mat[5]);
		zRot = 0;
	}

	(*x_result) = -xRot;
	(*y_result) = -yRot;
	(*z_result) = -zRot;
}
void tmat_extract_location_44d(tnsMatrix44d mat, real *x_result, real *y_result, real *z_result) {
	*x_result = mat[12];
	*y_result = mat[13];
	*z_result = mat[14];
}
void tmat_extract_uniform_scale_44d(tnsMatrix44d mat, real *result) {
	tnsVector3d v = { mat[0], mat[1], mat[2] };
	*result = tmat_length_3d(v);
}



#define L(row, col)  l[(col << 2) + row]
#define R(row, col)  r[(col << 2) + row]
#define P(row, col)  result[(col << 2) + row]

void tmat_print_matrix_44d(tnsMatrix44d l) {
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			printf("%.5f ", L(i, j));
		}
		printf("\n");
	}
}

void tmat_obmat_to_16d(float obmat[4][4], tnsMatrix44d out) {
	out[0] = obmat[0][0];
	out[1] = obmat[0][1];
	out[2] = obmat[0][2];
	out[3] = obmat[0][3];
	out[4] = obmat[1][0];
	out[5] = obmat[1][1];
	out[6] = obmat[1][2];
	out[7] = obmat[1][3];
	out[8]  = obmat[2][0];
	out[9]  = obmat[2][1];
	out[10] = obmat[2][2];
	out[11] = obmat[2][3];
	out[12] = obmat[3][0];
	out[13] = obmat[3][1];
	out[14] = obmat[3][2];
	out[15] = obmat[3][3];
}

void tmat_copy_matrix_44d(tnsMatrix44d from, tnsMatrix44d to) {
	memcpy(to, from, sizeof(tnsMatrix44d));
}
void tmat_multiply_44d(tnsMatrix44d result, tnsMatrix44d l, tnsMatrix44d r) {
	int i;
	for (i = 0; i < 4; i++) {
		real ai0 = L(i, 0), ai1 = L(i, 1), ai2 = L(i, 2), ai3 = L(i, 3);
		P(i, 0) = ai0 * R(0, 0) + ai1 * R(1, 0) + ai2 * R(2, 0) + ai3 * R(3, 0);
		P(i, 1) = ai0 * R(0, 1) + ai1 * R(1, 1) + ai2 * R(2, 1) + ai3 * R(3, 1);
		P(i, 2) = ai0 * R(0, 2) + ai1 * R(1, 2) + ai2 * R(2, 2) + ai3 * R(3, 2);
		P(i, 3) = ai0 * R(0, 3) + ai1 * R(1, 3) + ai2 * R(2, 3) + ai3 * R(3, 3);
	}
};
void tmat_inverse_44d(tnsMatrix44d inverse, tnsMatrix44d mat) {
	int i, j, k;
	double temp;
	tnsMatrix44d tempmat;
	real max;
	int maxj;

	tmat_load_identity_44d(inverse);

	tmat_copy_matrix_44d(mat, tempmat);

	for (i = 0; i < 4; i++) {
		/* Look for row with max pivot */
		max = fabsf(tempmat[i * 5]);
		maxj = i;
		for (j = i + 1; j < 4; j++) {
			if (fabsf(tempmat[j * 4 + i]) > max) {
				max = fabsf(tempmat[j * 4 + i]);
				maxj = j;
			}
		}

		/* Swap rows if necessary */
		if (maxj != i) {
			for (k = 0; k < 4; k++) {
				real t;
				t = tempmat[i * 4 + k];
				tempmat[i * 4 + k] = tempmat[maxj * 4 + k];
				tempmat[maxj * 4 + k] = t;

				t = inverse[i * 4 + k];
				inverse[i * 4 + k] = inverse[maxj * 4 + k];
				inverse[maxj * 4 + k] = t;
			}
		}

		//if (UNLIKELY(tempmat[i][i] == 0.0f)) {
		//	return false;  /* No non-zero pivot */
		//}

		temp = (double)tempmat[i * 5];
		for (k = 0; k < 4; k++) {
			tempmat[i * 4 + k] = (real)((double)tempmat[i * 4 + k] / temp);
			inverse[i * 4 + k] = (real)((double)inverse[i * 4 + k] / temp);
		}
		for (j = 0; j < 4; j++) {
			if (j != i) {
				temp = tempmat[j * 4 + i];
				for (k = 0; k < 4; k++) {
					tempmat[j * 4 + k] -= (real)((double)tempmat[i * 4 + k] * temp);
					inverse[j * 4 + k] -= (real)((double)inverse[i * 4 + k] * temp);
				}
			}
		}
	}
}
void tmat_make_translation_matrix_44d(tnsMatrix44d mTrans, real x, real y, real z) {
	tmat_load_identity_44d(mTrans);
	mTrans[12] = x;
	mTrans[13] = y;
	mTrans[14] = z;
}
void tmat_make_perspective_matrix_44d(tnsMatrix44d mProjection, real fFov_rad, real fAspect, real zMin, real zMax) {
	real yMax;
	real yMin;
	real xMin;
	real xMax;

	if (fAspect < 1) {
		yMax = zMin * tanf(fFov_rad * 0.5f);
		yMin = -yMax;
		xMin = yMin * fAspect;
		xMax = -xMin;
	}else {
		xMax = zMin * tanf(fFov_rad * 0.5f);
		xMin = -xMax;
		yMin = xMin / fAspect;
		yMax = -yMin;
	}

	tmat_load_identity_44d(mProjection);

	mProjection[0] = (2.0f * zMin) / (xMax - xMin);
	mProjection[5] = (2.0f * zMin) / (yMax - yMin);
	mProjection[8] = (xMax + xMin) / (xMax - xMin);
	mProjection[9] = (yMax + yMin) / (yMax - yMin);
	mProjection[10] = -((zMax + zMin) / (zMax - zMin));
	mProjection[11] = -1.0f;
	mProjection[14] = -((2.0f * (zMax * zMin)) / (zMax - zMin));
	mProjection[15] = 0.0f;
}
void tmat_make_z_tracking_matrix_44d(tnsMatrix44d mat, tnsVector3d this, tnsVector3d that, tnsVector3d up) {
	tnsVector4d fwd, l, t, rt;
	fwd[3] = l[3] = t[3] = rt[3] = 1;
	t[0] = up[0];
	t[1] = up[1];
	t[2] = up[2];
	fwd[0] = that[0] - this[0];
	fwd[1] = that[1] - this[1];
	fwd[2] = that[2] - this[2];

	tmat_load_identity_44d(mat);

	tmat_vector_cross_3d(l, t, fwd);
	tmat_vector_cross_3d(rt, fwd, l);

	tmat_normalize_self_3d(l);
	tmat_normalize_self_3d(rt);
	tmat_normalize_self_3d(fwd);

	mat[0] = l[0];
	mat[1] = l[1];
	mat[2] = l[2];

	mat[4] = rt[0];
	mat[5] = rt[1];
	mat[6] = rt[2];

	mat[8] = fwd[0];
	mat[9] = fwd[1];
	mat[10] = fwd[2];
}
void tmat_make_z_tracking_delta_matrix_44d(tnsMatrix44d mat, tnsVector3d delta, tnsVector3d up) {
	tnsVector4d fwd, l, t, rt;
	fwd[3] = l[3] = t[3] = rt[3] = 1;
	t[0] = up[0];
	t[1] = up[1];
	t[2] = up[2];
	fwd[0] = delta[0];
	fwd[1] = delta[1];
	fwd[2] = delta[2];

	tmat_load_identity_44d(mat);

	tmat_vector_cross_3d(l, t, fwd);
	tmat_vector_cross_3d(rt, fwd, l);

	tmat_normalize_self_3d(l);
	tmat_normalize_self_3d(rt);
	tmat_normalize_self_3d(fwd);

	mat[0] = l[0];
	mat[1] = l[1];
	mat[2] = l[2];

	mat[4] = rt[0];
	mat[5] = rt[1];
	mat[6] = rt[2];

	mat[8] = fwd[0];
	mat[9] = fwd[1];
	mat[10] = fwd[2];
}
void tmat_make_ortho_matrix_44d(tnsMatrix44d mProjection, real xMin, real xMax, real yMin, real yMax, real zMin, real zMax) {
	tmat_load_identity_44d(mProjection);

	mProjection[0] = 2.0f / (xMax - xMin);
	mProjection[5] = 2.0f / (yMax - yMin);
	mProjection[10] = -2.0f / (zMax - zMin);
	mProjection[12] = -((xMax + xMin) / (xMax - xMin));
	mProjection[13] = -((yMax + yMin) / (yMax - yMin));
	mProjection[14] = -((zMax + zMin) / (zMax - zMin));
	mProjection[15] = 1.0f;
}
void tmat_make_rotation_matrix_44d(tnsMatrix44d m, real angle_rad, real x, real y, real z)
{
	real mag, s, c;
	real xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

	s = (real)sin(angle_rad);
	c = (real)cos(angle_rad);

	mag = (real)sqrt(x * x + y * y + z * z);

	// Identity matrix
	if (mag == 0.0f) {
		tmat_load_identity_44d(m);
		return;
	}

	// Rotation matrix is normalized
	x /= mag;
	y /= mag;
	z /= mag;

#define M(row, col)  m[col * 4 + row]

	xx = x * x;
	yy = y * y;
	zz = z * z;
	xy = x * y;
	yz = y * z;
	zx = z * x;
	xs = x * s;
	ys = y * s;
	zs = z * s;
	one_c = 1.0f - c;

	M(0, 0) = (one_c * xx) + c;
	M(0, 1) = (one_c * xy) + zs;
	M(0, 2) = (one_c * zx) + ys;
	M(0, 3) = 0.0f;

	M(1, 0) = (one_c * xy) - zs;
	M(1, 1) = (one_c * yy) + c;
	M(1, 2) = (one_c * yz) + xs;
	M(1, 3) = 0.0f;

	M(2, 0) = (one_c * zx) + ys;
	M(2, 1) = (one_c * yz) - xs;
	M(2, 2) = (one_c * zz) + c;
	M(2, 3) = 0.0f;

	M(3, 0) = 0.0f;
	M(3, 1) = 0.0f;
	M(3, 2) = 0.0f;
	M(3, 3) = 1.0f;

#undef M
}
void tmat_make_rotation_x_matrix_44d(tnsMatrix44d m, real angle_rad) {
	tmat_load_identity_44d(m);
	m[5] = cos(angle_rad);
	m[6] = sin(angle_rad);
	m[9] = -sin(angle_rad);
	m[10] = cos(angle_rad);
}
void tmat_make_rotation_y_matrix_44d(tnsMatrix44d m, real angle_rad) {
	tmat_load_identity_44d(m);
	m[0] = cos(angle_rad);
	m[2] = -sin(angle_rad);
	m[8] = sin(angle_rad);
	m[10] = cos(angle_rad);
}
void tmat_make_rotation_z_matrix_44d(tnsMatrix44d m, real angle_rad) {
	tmat_load_identity_44d(m);
	m[0] = cos(angle_rad);
	m[1] = sin(angle_rad);
	m[4] = -sin(angle_rad);
	m[5] = cos(angle_rad);
}
void tmat_make_scale_matrix_44d(tnsMatrix44d m, real x, real y, real z) {
	tmat_load_identity_44d(m);
	m[0] = x;
	m[5] = y;
	m[10] = z;
}
void tmat_make_viewport_matrix_44d(tnsMatrix44d m, real w, real h, real Far, real Near) {
	tmat_load_identity_44d(m);
	m[0] = w / 2;
	m[5] = h / 2;
	m[10] = (Far - Near) / 2;
	m[12] = w / 2;
	m[13] = h / 2;
	m[14] = (Far + Near) / 2;
	m[15] = 1;
	//m[0] = 2/w;
	//m[5] = 2/h;
	//m[10] = 1;
	//m[12] = 2/w;
	//m[13] = 2/h;
	//m[14] = 1;
	//m[15] = 1;
}


real lanpr_LinearInterpolate(real L, real R, real T) {
	return tnsLinearItp(L, R, T);
}
void lanpr_LinearInterpolate2dv(real *L, real *R, real T, real *Result) {
	Result[0] = tnsLinearItp(L[0], R[0], T);
	Result[1] = tnsLinearItp(L[1], R[1], T);
}
void lanpr_LinearInterpolate3dv(real *L, real *R, real T, real *Result) {
	Result[0] = tnsLinearItp(L[0], R[0], T);
	Result[1] = tnsLinearItp(L[1], R[1], T);
	Result[2] = tnsLinearItp(L[2], R[2], T);
}

