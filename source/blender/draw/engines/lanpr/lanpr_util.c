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

nSafeStringCollection SSC;



void* nutCalloc(int size, int num) {
	return calloc(size, num);
}

void nutFreeMem(void** ptr){
	//free_total+=1;
	if (!*ptr) return;
	free(*ptr);
		*ptr = 0;
}

int nutFloatCompare(real l, real r){
	return (l>r - 0.00005 && l<r + 0.00005);
}

int nutSameAddress(void* l, void* r){
	return (l == r);
}


//===================================================================[slt]


void lstEmptyDirect(nListHandle* h) {
	h->pFirst = h->pLast = 0;
}


void lstPushSingle(void** Head, nListSingle* Item) {
	Item->pNext = *Head;
	*Head = Item;
}
void* lstPopSingle(void** Head, nListSingle* Item) {
	*Head = ((nListSingle*)(*Head))->pNext;
	Item->pNext = 0;
	return *Head;
}

int   lstHaveItemInList(nListHandle* Handle){
	if (Handle->pFirst)
		return 1;
	return 0;
};

void lstAppendItem2(nListHandle* Handle, void* Item) {

	nListItem2* li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pFirst)
		Handle->pFirst = Item;

	if (Handle->pLast)
		((nListItem2*)Handle->pLast)->pNext = li;

	li->pPrev = Handle->pLast;
	li->pNext = 0;
	Handle->pLast = li;

};
void lstPushItem2(nListHandle* Handle, void* Item) {

	nListItem2* li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pLast)
		Handle->pLast = Item;

	li->pNext = Handle->pFirst;

	if (Handle->pFirst)
		((nListItem2*)Handle->pFirst)->pPrev = Item;

	Handle->pFirst = li;

};
void* lstPopItem2(nListHandle* Handle) {
	void* popitem;
	nListItem2* next;
	if (!Handle->pFirst)
		return 0;

	popitem = Handle->pFirst;

	next = ((nListItem2*)Handle->pFirst)->pNext;
	if (!next) {
		Handle->pFirst = 0;
		Handle->pLast = 0;
	}
	else {
		Handle->pFirst = next;
		if (next)
			next->pPrev = 0;
	};

	return popitem;
};

void lstAppendItem3(nListHandle* Handle, void* Item) {

	nListItem3* li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pFirst)
		Handle->pFirst = Item;

	if (Handle->pLast)
		((nListItem3*)Handle->pLast)->pNext = li;

	li->pPrev = Handle->pLast;
	li->pNext = 0;
	Handle->pLast = li;

};
void lstPushItem3(nListHandle* Handle, void* Item) {

	nListItem3* li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pLast)
		Handle->pLast = Item;

	li->pNext = Handle->pFirst;

	if (Handle->pFirst)
		((nListItem3*)Handle->pFirst)->pPrev = Item;

	Handle->pFirst = li;

};
void* lstPopItem3(nListHandle* Handle) {
	void* popitem;
	nListItem3* next;
	if (!Handle->pFirst)
		return 0;

	popitem = Handle->pFirst;

	next = ((nListItem3*)Handle->pFirst)->pNext;
	if (!next) {
		Handle->pFirst = 0;
		Handle->pLast = 0;
	}
	else {
		Handle->pFirst = next;
		if (next)
			next->pPrev = 0;
	};

	return popitem;
};


void* lstGetTop(nListHandle* Handle){
	return Handle->pFirst;
};
int lstRemoveItem2(nListHandle* Handle, nListItem2* li) {
	if (!li->pPrev)
		Handle->pFirst = li->pNext;
	else
		((nListItem2*)li->pPrev)->pNext = li->pNext;

	if (!li->pNext)
		Handle->pLast = li->pPrev;
	else
		((nListItem2*)li->pNext)->pPrev = li->pPrev;

	li->pNext = li->pPrev = 0;
};
int lstRemoveItem3(nListHandle* Handle, nListItem2* li) {
	if (!li->pPrev)
		Handle->pFirst = li->pNext;
	else
		((nListItem3*)li->pPrev)->pNext = li->pNext;

	if (!li->pNext)
		Handle->pLast = li->pPrev;
	else
		((nListItem3*)li->pNext)->pPrev = li->pPrev;

	li->pNext = li->pPrev = 0;
};

int lstRemoveSegment(nListHandle* Handle, nListItem* Begin, nListItem* End) {
	if (!Begin->pPrev)
		Handle->pFirst = End->pNext;
	else
		((nListItem*)Begin->pPrev)->pNext = End->pNext;

	if (!End->pNext)
		Handle->pLast = Begin->pPrev;
	else
		((nListItem*)End->pNext)->pPrev = Begin->pPrev;

	End->pNext = Begin->pPrev = 0;
};
void lstInsertItemBefore(nListHandle* Handle, nListItem* toIns, nListItem* pivot){
	if (!pivot) { lstPushItem(Handle, toIns); return; }

	if (pivot->pPrev){
		((nListItem*)pivot->pPrev)->pNext = toIns;
		toIns->pPrev = pivot->pPrev;
	}
	else {
		Handle->pFirst = toIns;
	}

	toIns->pNext = pivot;
	pivot->pPrev = toIns;
};
void lstInsertItemAfter(nListHandle* Handle, nListItem* toIns, nListItem* pivot) {
	if (!pivot) { lstAppendItem(Handle, toIns); return; }

	if (pivot->pNext) {
		((nListItem*)pivot->pNext)->pPrev = toIns;
		toIns->pNext = pivot->pNext;
	}
	else {
		Handle->pLast = toIns;
	}

	toIns->pPrev = pivot;
	pivot->pNext = toIns;
}
void lstInsertSegmentBefore(nListHandle* Handle, nListItem* Begin, nListItem* End, nListItem* pivot) {

	if (pivot->pPrev) {
		((nListItem*)pivot->pPrev)->pNext = Begin;
		Begin->pPrev = pivot->pPrev;
	}
	else {
		Handle->pFirst = Begin;
	}

	End->pNext = pivot;
	pivot->pPrev = End;
};
void lstInsertSegmentAfter(nListHandle* Handle, nListItem* Begin, nListItem* End, nListItem* pivot) {
	if (pivot->pNext) {
		((nListItem*)pivot->pNext)->pPrev = End;
		End->pNext = pivot->pNext;
	}else {
		Handle->pLast = End;
	}

	Begin->pPrev = pivot;
	pivot->pNext = Begin;
}


void* lstAppendPointerOnly(nListHandle* h, void* p) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = CreateNew(nListItemPointer);
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstAppendPointerSizedOnly(nListHandle* h, void* p, int size) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = calloc(1, size);
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstPushPointerOnly(nListHandle* h, void* p) {
	nListItemPointer* lip = 0;
	if (!h) return 0;
	lip = CreateNew(nListItemPointer);
	lip->p = p;
	lstPushItem(h, lip);
	return lip;
}
void* lstPushPointerSizedOnly(nListHandle* h, void* p, int size) {
	nListItemPointer* lip = 0;
	if (!h) return 0;
	lip = calloc(1, size);
	lip->p = p;
	lstPushItem(h, lip);
	return lip;
}

void* lstAppendPointer(nListHandle* h, void* p) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = memAquireOnly(sizeof(nListItemPointer));
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstAppendPointerSized(nListHandle* h, void* p,int size) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = memAquireOnly(size);
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstPushPointer(nListHandle* h, void* p) {
	nListItemPointer* lip=0;
	if (!h) return 0;
	lip = memAquireOnly(sizeof(nListItemPointer));
	lip->p = p;
	lstPushItem(h, lip);
	return lip;
}
void* lstPushPointerSized(nListHandle* h, void* p,int size) {
	nListItemPointer* lip = 0;
	if (!h) return 0;
	lip = memAquireOnly(size);
	lip->p = p;
	lstPushItem(h, lip);
	return lip;
}

void* lstAppendPointerStatic(nListHandle* h, nStaticMemoryPool* smp, void* p) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = memStaticAquire(smp,sizeof(nListItemPointer));
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstAppendPointerStaticSized(nListHandle* h, nStaticMemoryPool* smp, void* p, int size) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = memStaticAquire(smp, size);
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstPushPointerStatic(nListHandle* h, nStaticMemoryPool* smp, void* p) {
	nListItemPointer* lip = 0;
	if (!h) return 0;
	lip = memStaticAquire(smp, sizeof(nListItemPointer));
	lip->p = p;
	lstPushItem(h, lip);
	return lip;
}
void* lstPushPointerStaticSized(nListHandle* h, nStaticMemoryPool* smp, void* p, int size) {
	nListItemPointer* lip = 0;
	if (!h) return 0;
	lip = memStaticAquire(smp, size);
	lip->p = p;
	lstPushItem(h, lip);
	return lip;
}

void* lstPopPointerOnly(nListHandle* h) {
	nListItemPointer* lip;
	void* rev = 0;
	if (!h) return 0;
	lip = lstPopItem(h);
	rev = lip ? lip->p : 0;
	FreeMem(lip);
	return rev;
}
void lstRemovePointerItemOnly(nListHandle* h, nListItemPointer* lip) {
	lstRemoveItem(h, (void*)lip);
	FreeMem(lip);
}
void lstRemovePointerOnly(nListHandle* h, void* p) {
	nListItemPointer* i;
	for (i = h->pFirst; i; i = i->pNext) {
		if (i->p == p) {
			lstRemovePointerItem(h, i);
			break;
		}
	}
}
void lstClearPointerOnly(nListHandle* h) {
	nListItemPointer* i;
	while (h&&h->pFirst) {
		lstPopPointer(h);
	}
}
void lstGeneratePointerListOnly(nListHandle* from1, nListHandle* from2, nListHandle* to) {
	nListItemPointer* lip = from2 ? from2->pLast : 0;

	while (lip) {
		lstPushPointer(to, lip->p);
		lip = lip->pPrev;
	}

	lip = from1 ? from1->pLast : 0;

	while (lip) {
		lstPushPointer(to, lip->p);
		lip = lip->pPrev;
	}
}

void* lstPopPointer(nListHandle* h) {
	nListItemPointer* lip;
	void* rev=0;
	if (!h) return 0;
	lip = lstPopItem(h);
	rev = lip ? lip->p : 0;
	memFree(lip);
	return rev;
}
void lstRemovePointerItem(nListHandle* h, nListItemPointer* lip) {
	lstRemoveItem(h, (void*)lip);
	memFree(lip);
}
void lstRemovePointer(nListHandle* h, void* p) {
	nListItemPointer* i;
	for (i = h->pFirst; i; i = i->pNext) {
		if (i->p == p) {
			lstRemovePointerItem(h, i);
			break;
		}
	}
}
void lstClearPointer(nListHandle* h) {
	nListItemPointer* i;
	while (h&&h->pFirst) {
		lstPopPointer(h);
	}
}
void lstGeneratePointerList(nListHandle* from1, nListHandle* from2, nListHandle* to) {
	nListItemPointer* lip = from2?from2->pLast:0;

	while (lip) {
		lstPushPointer(to, lip->p);
		lip = lip->pPrev;
	}

	lip = from1?from1->pLast:0;

	while (lip) {
		lstPushPointer(to, lip->p);
		lip = lip->pPrev;
	}
}


void* lstAppendPointerStaticPool(nStaticMemoryPool* mph, nListHandle* h, void* p) {
	nListItemPointer* lip;
	if (!h) return 0;
	lip = memStaticAquire(mph, sizeof(nListItemPointer));
	lip->p = p;
	lstAppendItem(h, lip);
	return lip;
}
void* lstPopPointerNoFree(nListHandle* h) {
	nListItemPointer* lip;
	void* rev = 0;
	if (!h) return 0;
	lip = lstPopItem(h);
	rev = lip ? lip->p : 0;
	return rev;
}
void lstRemovePointerItemNoFree(nListHandle* h, nListItemPointer* lip) {
	lstRemoveItem(h, (void*)lip);
}


void lstCopyHandle(nListHandle* target, nListHandle* src){
	target->pFirst = src->pFirst;
	target->pLast = src->pLast;
};
void lstClearHandle(nListHandle* h){
	h->pFirst = 0;
	h->pLast = 0;
}
void lstClearPrevNext(nListItem* li){
	li->pNext = 0;
	li->pPrev = 0;
}

void lstMoveUp(nListHandle* h, nListItem* li) {
	void* pprev = li->pPrev ? ((nListItem*)li->pPrev)->pPrev : 0;
	if (!h || !li)
		return;
	if (li == h->pFirst) return;
	else {
		if (li == h->pLast) h->pLast = li->pPrev;
		((nListItem*)li->pPrev)->pNext = li->pNext;
		((nListItem*)li->pPrev)->pPrev = li;
		if(li->pNext)((nListItem*)li->pNext)->pPrev = li->pPrev;
		li->pNext = li->pPrev;
		li->pPrev = pprev;
		if (pprev) ((nListItem*)pprev)->pNext = li;
	}
	if (!li->pPrev) h->pFirst = li;
}
void lstMoveDown(nListHandle* h, nListItem* li) {
	void* ppnext = li->pNext ? ((nListItem*)li->pNext)->pNext : 0;
	if (!h || !li)
		return;
	if (li == h->pLast) return;
	else {
		if (li == h->pFirst) h->pFirst = li->pNext;
		((nListItem*)li->pNext)->pPrev = li->pPrev;
		((nListItem*)li->pNext)->pNext= li;
		if (li->pPrev)((nListItem*)li->pPrev)->pNext = li->pNext;
		li->pPrev = li->pNext;
		li->pNext = ppnext;
		if (ppnext) ((nListItem*)ppnext)->pPrev = li;
	}
	if (!li->pNext) h->pLast = li;
}

void lstForAllItemsDo(nListDoFunc func, nListHandle* hList){
	nListItem* it = hList->pFirst;
	for (; it; it = it->pNext){
		func(it);
	}
};
void lstForAllItemsDoLNRR(nListNonRecursiveDoFunc func, nListHandle* hList){
	nListItem* it = hList->pFirst;
	for (; it; it = it->pNext){
		func(0,it,0);
	}
};
void lstForAllItemsDo_DirectFree(nListDoFunc func, nListHandle* hList){
	nListItem* it;
	while (it = lstPopItem(hList)){
		if (func) func(it);
		FreeMem(it);
	}
};
void lstForAllItemsDo_arg_ptr(nListDoFuncArgp func, nListHandle* hList, void* arg){
	nListItem* it = hList->pFirst;
	for (; it; it = it->pNext){
		func(it, arg);
	};
};
void lstForAllItemsDo_NonRecursive_Root(nListHandle* FirstHandle, nListNonRecursiveDoFunc func, int bFreeItem, void* custom_data, nListCustomDataRemover remover){
	nListItem* li = 0, *NextLi;
	nListNonRecursiveRoot root = { 0 };
	nListNonRecursiveItem* nrItem = CreateNew(nListNonRecursiveItem);

	nrItem->bFreeList = bFreeItem;
	nrItem->func = func;
	nrItem->CustomData = custom_data;
	nrItem->remover = remover;
	lstCopyHandle(&nrItem->handle, FirstHandle);

	lstAppendItem(&root.NSItems, nrItem);

	while (lstHaveItemInList(&root.NSItems)){
		nrItem = lstPopItem(&root.NSItems);

		for (li = nrItem->handle.pFirst; li/*!=nrItem->handle.pLast*/; li = NextLi){
			if (nrItem->func) nrItem->func(&root, li, custom_data);
			NextLi = li->pNext;
			if (nrItem->bFreeList){
				nListItem* fli = li;
				FreeMem(fli);
			}
			if (li == nrItem->handle.pLast) break;
		}
		if (nrItem->remover) nrItem->remover(nrItem->CustomData);
		FreeMem(nrItem);
	}
};
void lstAddNonRecursiveListHandle(nListNonRecursiveRoot* root, nListHandle* newHandle, nListNonRecursiveDoFunc nrFunc, int bFreeList, void* custom_data, nListCustomDataRemover remover){
	nListNonRecursiveItem* nrItem = CreateNew(nListNonRecursiveItem);

	nrItem->bFreeList = bFreeList;
	nrItem->func = nrFunc;
	nrItem->CustomData = custom_data;
	nrItem->remover = remover;
	lstCopyHandle(&nrItem->handle, newHandle);

	lstAppendItem(&root->NSItems, nrItem);
};
void lstCopy_NonRecursive_Root(nListHandle* FromHandle, nListHandle* ToHandle, int SizeEachNode, nListNonRecursiveCopyFunc func, void* custom_data, nListCustomDataRemover remover){
	nListItem* li = 0, *tli = 0;
	nListNonRecursiveRoot root = { 0 };
	nListNonRecursiveItem* nrItem = CreateNew(nListNonRecursiveItem);
	nListItem* NextLi;

	nrItem->CopyFunc = func;
	lstCopyHandle(&nrItem->handle, FromHandle);
	nrItem->ToHandle = ToHandle;//Pointer
	lstClearHandle(ToHandle);
	nrItem->CustomData = custom_data;
	nrItem->remover = remover;
	nrItem->SizeEachNode = SizeEachNode;

	lstAppendItem(&root.NSItems, nrItem);

	while (lstHaveItemInList(&root.NSItems)){
		nrItem = lstPopItem(&root.NSItems);
		if (nrItem->CopyFunc){
			for (li = nrItem->handle.pFirst; li; li = li->pNext){
				tli = CreateNew_Size(nrItem->SizeEachNode);

				nrItem->CopyFunc(&root, li, tli, nrItem->CustomData);

				lstClearPrevNext(tli);
				lstAppendItem(nrItem->ToHandle, tli);
			}
			if (nrItem->remover) nrItem->remover(nrItem->CustomData);
		}else if (nrItem->func){
			for (li = nrItem->handle.pFirst; li/*!=nrItem->handle.pLast*/; li = NextLi){
				if (nrItem->func) nrItem->func(&root, li, custom_data);
				NextLi = li->pNext;
				if (nrItem->bFreeList){
					nListItem* fli = li;
					FreeMem(fli);
				}
				if (li == nrItem->handle.pLast) break;
			}
			if (nrItem->remover) nrItem->remover(nrItem->CustomData);
		}
		FreeMem(nrItem);
	}
};
void lstAddNonRecursiveListCopier(nListNonRecursiveRoot* root, nListHandle* oldHandle, nListHandle* newHandle, int sizeEach, nListNonRecursiveCopyFunc nrCpyFunc, void* custom_data, nListCustomDataRemover remover){
	nListNonRecursiveItem* nrItem = CreateNew(nListNonRecursiveItem);

	nrItem->CopyFunc = nrCpyFunc;
	lstCopyHandle(&nrItem->handle, oldHandle);
	nrItem->ToHandle = newHandle;
	nrItem->CustomData = custom_data;
	nrItem->remover = remover;
	nrItem->SizeEachNode = sizeEach;

	lstAppendItem(&root->NSItems, nrItem);
};

void* lstFindItem(void* CmpData, nCompareFunc func, nListHandle* hList){
	nListItem* it;

	if (!CmpData || !hList)
		return 0;

	it = hList->pFirst;
	for (; it; it = it->pNext){
		if (func(it, CmpData))
			return it;
	};
	return 0;
};
void lstCombineLists(nListHandle* dest, nListHandle* src){
	if ((!dest) || (!src))
		return;

	if ((!dest->pFirst) && (!dest->pLast)){
		dest->pFirst = src->pFirst;
		dest->pLast = src->pLast;
	}
	else{
		if (src->pLast){
			((nListItem*)src->pFirst)->pPrev = dest->pLast;
			((nListItem*)dest->pLast)->pNext = src->pFirst;
			dest->pLast = src->pLast;
		}
	}

	src->pFirst = 0;
	src->pLast = 0;
}
void lstDestroyList(nListHandle* hlst){
	nListItem* li, *nextli;
	for (li = hlst->pFirst; li; li = nextli){
		nextli = li->pNext;
		memFree(li);
	}
}
void lstDestroyListA(nListHandle* hlst) {
	nListItem* li, *nextli;
	for (li = hlst->pFirst; li; li = nextli) {
		nextli = li->pNext;
		FreeMem(li);
	}
}
void lstDestroyList_User(nListHandle* hlst, nListDoFunc func){
	nListItem* it = hlst->pFirst;
	for (; it; it = it->pNext){
		func(it);
		FreeMem(it);
	}
};
void lstCopyList(nListHandle* hOldlst, nListHandle* hNewList, int SizeEachNode, nCopyListFunc func){
	nListItem* li, *nextli, *newli;
	for (li = hOldlst->pFirst; li; li = nextli){
		newli = (nListItem*)CreateNew_Size(SizeEachNode);
		func(li, newli);
		lstAppendItem(hNewList, newli);
		nextli = li->pNext;
	}
}

void* lstReMatch(nListHandle* SearchHandle, nListHandle* CurrentHandle, void* ItemToFind){
	nListItem* sl = 0, *rl = 0;

	if (!SearchHandle || !CurrentHandle || !ItemToFind) return 0;

	sl = SearchHandle->pFirst; rl = CurrentHandle->pFirst;

	while (sl && rl){
		if (ItemToFind == sl){
			return rl;
		}
		else{
			sl = sl->pNext;
			rl = rl->pNext;
		}
	}
	return 0;
}
//void* lstReMatchEx(nListHandle* SearchHandle, nListHandle* CurrentHandle, void* ItemToFind, MatcherFunc func){
//	nListItem* sl = 0, *rl = 0;
//
//	if (!SearchHandle || !CurrentHandle || !ItemToFind) return 0;
//
//	sl = SearchHandle->pFirst; rl = CurrentHandle->pFirst;
//
//	while (sl && rl){
//		if (func(ItemToFind, sl)){
//			return rl;
//		}
//		else{
//			sl = sl->pNext;
//			rl = rl->pNext;
//		}
//	}
//	return 0;
//}

void lstAddElement(nListHandle* hlst, void* ext){
	nElementListItem* eli = CreateNew(nElementListItem);
	eli->Ext = ext;
	lstAppendItem(hlst, eli);
}
void lstDestroyElementList(nListHandle* hlst){
	nElementListItem* eli,*NextEli;
	for (eli = hlst->pFirst; eli; eli = NextEli){
		lstRemoveItem(hlst, (void*)eli);
		NextEli = eli->Item.pNext;
		FreeMem(eli);
	}
}


unsigned char hsh256DoHashCSTR(char* buckle){
	int i, len = 0;
	unsigned char rev = 0;

	if (buckle) len = strlen(buckle);

	for (i = 0; i<len; i++){
		rev = rev * 31 + buckle[i];
	}

	return rev;
}

void hsh256InsertItemCSTR(nHash256* hash, nListItem* li, char* buckle){
	int a = hsh256DoHashCSTR(buckle);
	lstAppendItem(&hash->Entries[a], li);
};
void hsh256InsertItem(nHash256* hash, nListItem* li, unsigned char buckle) {
	lstAppendItem(&hash->Entries[buckle], li);
};
void hsh65536InsertItem(nHash65536* hash, nListItem* li, long buckle) {
	lstAppendItem(&hash->Entries[(unsigned short)((buckle >> 10))], li);
	//hsh256InsertItem(&hash->HashHandles[(unsigned char)((buckle >> 8) / 8)], li, (unsigned char)(buckle/8));
	//printf("%d %d\n", (unsigned char)(buckle >> 5), (unsigned char)(buckle >> 6));
};


nListItem* hsh256FindItemCSTR(nHash256* hash, nCompareFunc func, char* buckle){
	unsigned char hsh;

	hsh = hsh256DoHashCSTR(buckle);

	//if(hash->Entries[hsh].pFirst == hash->Entries[hsh].pLast)
	//    return hash->Entries[hsh].pFirst;

	return lstFindItem(buckle, func, &hash->Entries[hsh]);
}

void memInitPool(nMemoryPool* mph, int NodeSize) {
	int Count = 4096;

	if (mph->NodeSize) return;

	mph->NodeSize = NodeSize;
	mph->CountPerPool = Count;

	return;
}
void memInitPoolSmall(nMemoryPool* mph, int NodeSize) {
	int Count = 16;

	if (mph->NodeSize) return;

	mph->NodeSize = NodeSize;
	mph->CountPerPool = Count;

	return;
}
inline nMemoryPoolPart* memNewPoolPart(nMemoryPool* mph) {

	if (!mph->NodeSize)return 0;

	int RealNodeSize = mph->NodeSize + sizeof(nMemoryPoolNode);
	int NodeCount = mph->CountPerPool;
	int TotalSize = sizeof(nMemoryPoolPart) + NodeCount * RealNodeSize;
	int i;
	nMemoryPoolPart* mp = calloc(1, TotalSize);

	void* BeginMem = ((BYTE*)mp) + sizeof(nMemoryPoolPart);

	mp->PoolRoot = mph;

	nMemoryPoolNode* mpn;

	mp->FreeMemoryNodes.pFirst = mp->FreeMemoryNodes.pLast = BeginMem;

	mpn = (void*)((BYTE*)BeginMem);
	mpn->InPool = mp;

	for (i = 1; i < NodeCount; i++) {
		mpn = (void*)(((BYTE*)BeginMem) + RealNodeSize*i);
		mpn->InPool = mp;
		lstAppendItem(&mp->FreeMemoryNodes, mpn);
	}
	lstPushItem(&mph->Pools, mp);

	return mp;
}

nDBInst* memGetDBInst(void* mem) {
	nMemoryPoolNode* mpn = (void*)(((BYTE*)mem) - sizeof(nMemoryPoolNode));
	return mpn->DBInst;
}
void memAssignDBInst(void* mem, nDBInst* DBInst) {
	nMemoryPoolNode* mpn = (void*)(((BYTE*)mem) - sizeof(nMemoryPoolNode));
	mpn->DBInst = DBInst;
}
void memFree(void* Data) {
	if (!Data)return;
	nMemoryPoolNode* mpn = (void*)(((BYTE*)Data) - sizeof(nMemoryPoolNode));
	nMemoryPoolPart* mp = mpn->InPool;
	nMemoryPool* mph = mp->PoolRoot;
	lstRemoveItem(&mp->MemoryNodes, (void*)mpn);
	lstAppendItem(&mp->FreeMemoryNodes, (void*)mpn);
	memset(Data, 0, mph->NodeSize);
	if (!mp->MemoryNodes.pFirst) {
		lstRemoveItem(&mph->Pools, (void*)mp);
		FreeMem(mp);
	}
	//if (!mph->Pools.pFirst) {
	//	mph->CountPerPool = 0;
	//	mph->NodeSize = 0;
	//}
}
void memDestroyPool(nMemoryPool* Handle) {
	nMemoryPool* mp;
	while ((mp = lstPopItem(&Handle->Pools))) {
		FreeMem(mp);
	}
}

nStaticMemoryPoolNode* memNewStaticPool(nStaticMemoryPool* smp) {
	nStaticMemoryPoolNode* smpn = calloc(1, NUL_MEMORY_POOL_128MB);
	smpn->UsedByte = sizeof(nStaticMemoryPoolNode);
	lstPushItem(&smp->Pools, smpn);
	return smpn;
}
void* memStaticAquire(nStaticMemoryPool*smp, int size) {
	nStaticMemoryPoolNode* smpn = smp->Pools.pFirst;
	void* ret;

	if (!smpn || (smpn->UsedByte + size) > NUL_MEMORY_POOL_128MB)
		smpn = memNewStaticPool(smp);

	ret = ((BYTE*)smpn) + smpn->UsedByte;

	smpn->UsedByte += size;

	return ret;
}
void* memStaticAquireThread(nStaticMemoryPool*smp, int size) {
	nStaticMemoryPoolNode* smpn = smp->Pools.pFirst;
	void* ret;

	//EnterCriticalSection(&smp->csMem);

	if (!smpn || (smpn->UsedByte + size) > NUL_MEMORY_POOL_128MB)
		smpn = memNewStaticPool(smp);

	ret = ((BYTE*)smpn) + smpn->UsedByte;

	smpn->UsedByte += size;

	//LeaveCriticalSection(&smp->csMem);

	return ret;
}
void* memStaticDestroy(nStaticMemoryPool*smp) {
	nStaticMemoryPoolNode* smpn;
	void* ret;

	while (smpn = lstPopItem(&smp->Pools)) {
		FreeMem(smpn);
	}

	smp->EachSize = 0;

	return ret;
}

//=======================================================================[str]

char* strGetNextString(char** pivot,char* NextMark){
	int lenth = 0;
	char* countP = *pivot;
	char* result = 0;
	int FloatArg = 0;
	int i;

	if (**pivot == '\0')
		return 0;

	if (*NextMark == '~') FloatArg = 1;

	//   container@identifier=window  container#window  contianer%
	while (!lenth) {
		for (countP; *countP != '.' && *(*pivot) != '\0' && *countP != '\0'
			&& *countP != '@'
			&& *countP != '='
			&& *countP != '#'
			&& *countP != '$';
			countP++) {
			lenth++;
		}
		if (lenth || (*countP)==0) break;
		(*pivot)++;
		countP++;
	}

	*NextMark = (*pivot)[lenth];
	if (!(*NextMark))*NextMark = '.';

	if (lenth){
		result = CreateNewBuffer(char, lenth + 1);

		for (i = 0; i<lenth; i++)
			result[i] = (*pivot)[i];

		result[lenth] = '\0';

		if ((*pivot)[lenth] == '\0') *pivot = &((*pivot)[lenth]);
		else (*pivot) += lenth + 1;

		return result;
	}
	else{
		return 0;
	}
};
int  strGetStringTerminateBy(char* content, char terminator, char* Out){
	int Ofst = 0;
	int Skip = 0;
	int i = 0;

	if ((!content) || (*content == '\0'))
		return 0;

	for (Ofst; content[Ofst] != terminator && content[Ofst] != '\0'; Ofst++){
		Out[i] = content[Ofst];
		i++;
	}
	Out[i] = 0;

	return i;
};
char*  strGetNewStringTerminateBy_PivotOver(char* content, char terminator,char** NewPivot,int IgnoreSpace) {
	int Ofst = 0;
	int Skip = 0;
	int i = 0;
	char* NewString;

	if (!content || *content == '\0')
		return 0;

	if(IgnoreSpace)	for (i; content[i] == ' '; i++) {
		;
	}

	for (Ofst; content[Ofst] != terminator && content[Ofst] != '\0'; Ofst++) {
		//Out[i] = content[Ofst];
		//i++;
		;
	}

	NewString = CreateNewBuffer(char, Ofst+1-i);

	memcpy(NewString, &content[i], sizeof(char)*(Ofst-i));

	NewString[Ofst-i] = '\0';

	*NewPivot = &content[Ofst + 1];

	return NewString;
};

int strHeadOfStringMatch(char* Str, char* SubStr) {
	int len = strlen(SubStr);
	int i = 0;
	for (i; i < len; i++) {
		if (Str[i] != SubStr[i]) return 0;
	}
	return 1;
}
int strSkipSegmet(char** pivot, char* content) {
	if (!pivot || !(*pivot) || !(*(*pivot)) || !content) return 0;

	if (strHeadOfStringMatch(*pivot, content)) {
		(*pivot) += strlen(content);
		return 1;
	}
	return 0;
}
char* strgetLastSegmentSeperateBy(char* Content, char Seperator) {
	char* p = Content;
	char* pn = Content;
	while (1) {
		while (*pn != Seperator) {
			if (!(*pn)) return p;
			pn++;
		}
		pn++;
		p = pn;
	}
}
void strDiscardLastSegmentSeperateBy(char* Content, char Seperator) {
	char* p = Content;
	char* pn = Content;
	while (1) {
		while (*pn != Seperator) {
			if (!(*pn)) {
				*p = 0;
				return;
			}
			pn++;
		}
		p = pn;
		pn++;
	}
}
void strDiscardSameBeginningSeperatedBy(char* s1, char* s2, char** Result1, char** Result2, char Seperator) {
	int i = 0;
	int p = 0;
	while (s1[i] == s2[i]) {
		i++;
		if (s1[i] == Seperator) p = i;
		if (!s1[i]) { p = i; break; }
		if (!s2[i]) { p = i; break; }
	}
	*Result1 = &s1[p];
	*Result2 = &s2[p];
}
int strCountSegmentSeperateBy(char* Content, char Seperator) {
	char* p = Content;
	char* pn = Content;
	int c = Content[0]?(Content[0]==Seperator?0:1):0;
	while (1) {
		while (*pn != Seperator) {
			if (!(*pn)) {
				if ((*p) == Seperator) c--;

				return c;
			}
			p = pn;
			pn++;
		}
		c++;
		pn++;
	}

	return c;
}

void strMakeDifferentName(char* Target) {
	char* p = strgetLastSegmentSeperateBy(Target, '.');
	int Temp;
	if (!sscanf(p, "%d", &Temp)) {
		int l = strlen(p);
		if (p[l - 1] != '.')	strcat(p, ".");
		strPrintIntAfter(Target, 0, 001);
	}else {
		sprintf(p, "%d", Temp + 1);
	};
}


void strReplaceCharacter(char* Str, char Find, char Replace) {
	char *p = Str;
	if (!p) return;
	while(*p){
		if (*p == Find) *p = Replace;
		p++;
	}
}
void strToUpperCase(char* Str) {
	char *p = Str;
	if (!p) return;
	while (*p) {
		if (*p>='a' && *p<='z') *p += 'A'-'a';
		p++;
	}
}
void strToLowerCase(char* Str) {
	char *p = Str;
	if (!p) return;
	while (*p) {
		if (*p >= 'A' && *p <= 'A') *p -= 'A'-'a';
		p++;
	}
}

nStringSplitor* strSplitPath(char* path){
	nStringPart*   sp;
	nStringSplitor* ss;
	char*          pivot = path;
	char*          temp_result;
	char           Type = '.';
	char           NextType='.';

	if (!path || !path[0])
		return 0;

	ss = memAquireOnly(sizeof(nStringSplitor));

	while (temp_result = strGetNextString(&pivot,&NextType)){
		if (*temp_result != '\0'){
			sp = memAquireOnly(sizeof(nStringPart));
			sp->Content = temp_result;
			lstAppendItem(&ss->parts, sp);
			ss->NumberParts += 1;
			if (NextType == '$') sp->Type = '$'; else sp->Type = Type;
			if (sp->Type == '=') {
				if (sp->Content[0] >= '0' && sp->Content[0] <= 9) {
					sscanf(sp->Content, "%d", &sp->IntValue);
				}
			}
			if (NextType == '$') NextType = '.';
			Type = NextType;
		}
	}

	if (ss->NumberParts == 0){
		strDestroyStringSplitor(&ss);
		return 0;
	}

	return ss;
};

void DF_ClearStingParts(nStringPart* sp){
	FreeMem(sp->Content);
};

int strDestroyStringSplitor(nStringSplitor** ss){
	if (!(*ss))
		return 0;

	lstForAllItemsDo(DF_ClearStingParts, &(*ss)->parts);
	lstDestroyList(&(*ss)->parts);
	memFree(*ss);
	*ss = 0;

	return 1;
}

char buff[128];
int strMakeInstructions(nStringSplitor** result, char* content){
	nStringPart*   sp;
	nStringSplitor* ss = *result;
	char*          pivot = content;
	char*          temp_result;

	if (!content || !content[0])
		return 0;

	if(!ss) ss = *result = memAquireOnly(sizeof(nStringSplitor));

	while (temp_result = strGetNewStringTerminateBy_PivotOver(pivot,'=',&pivot, 0)) {
		if (*temp_result != '\0') {
			sp = memAquireOnly(sizeof(nStringPart));
			sp->Content = temp_result;
			lstAppendItem(&ss->parts, sp);
			ss->NumberParts += 1;
		}
		temp_result = strGetNewStringTerminateBy_PivotOver(pivot, ';', &pivot, 0);
		if (!temp_result) break;
		if (*temp_result != '\0') {
			sp = memAquireOnly(sizeof(nStringPart));
			sp->Content = temp_result;
			lstAppendItem(&ss->parts, sp);
			ss->NumberParts += 1;
			if (temp_result[0] >= '0' && temp_result[0] <= '9' || temp_result[0]=='-') {
				sscanf(temp_result, "%d", &sp->IntValue);
				sscanf(temp_result, "%lf", &sp->FloatValue);
			}
		}
	}

	if (ss->NumberParts == 0) {
		strDestroyStringSplitor(&ss);
		return 0;
	}

	return 1;
}
nStringPart* strGetArgument(nStringSplitor* ss, char* content) {
	nStringPart* sp;
	if (!ss) return 0;
	for (sp = ss->parts.pFirst; sp; sp = sp->Item.pNext?((nListItem*)sp->Item.pNext)->pNext:0) {
		if (strIsTheSame(content, sp->Content))
			return sp->Item.pNext;
	}
	return 0;
}
char* strGetArgumentString(nStringSplitor* ss,char* content) {
	nStringPart* sp;
	if (!ss) return 0;
	for (sp = ss->parts.pFirst; sp; sp = sp->Item.pNext?((nListItem*)sp->Item.pNext)->pNext:0) {
		if (strIsTheSame(content, sp->Content) )
			return sp->Item.pNext?((nStringPart*)sp->Item.pNext)->Content:0;
	}
	return 0;
}
int strArgumentMatch(nStringSplitor* ss, char* id, char* value) {
	nStringPart* sp;
	if (!ss) return 0;
	for (sp = ss->parts.pFirst; sp; sp = sp->Item.pNext ? ((nListItem*)sp->Item.pNext)->pNext : 0) {
		if (strIsTheSame(id, sp->Content))
			return (strIsTheSame(((nStringPart*)sp->Item.pNext)->Content,value));
	}
	return 0;
}
int strGetIntSimple(char* content) {
	int a;
	sscanf(content, "%d", &a);
	return a;
}
real strGetFloatSimple(char* content) {
	real a;
	sscanf(content, "%lf", &a);
	return a;
}


void strConvInt_CString(int src, char* dest, int lenth){
	sprintf_s(dest, lenth, "%d", src);
};
void strConvFloat_CString(real src, char* dest, int lenth){
	sprintf_s(dest, lenth, "%lf", src);
};

void strCopyFull(char* dest, char* src){
	if (src && dest)
		strcpy(dest, src);
}
void strCopySized(char* dest, int LenthLim, char* src){
	if (src && dest)
		strcpy_s(dest, LenthLim + 1, src);
}
#define strAppend strcat_s
void strPrintFloatAfter(char* dest, int LenthLim, int bits, real data){
	char temp[64];
	sprintf_s(temp, LenthLim, "%.*lf", bits, data);
	strcat_s(dest, LenthLim, temp);
}
void strPrintIntAfter(char* dest, int LenthLim, int data){
	char temp[64];
	sprintf(&temp[0], "%d", data);
	strcat(dest, temp);
}
int  strIsTheSame(char* src, char*dest) {
	return (src && !strcmp(src, dest));
}


void strSafeDestroy(nSafeString** ss) {
	if (!*ss)
		return;
	lstRemoveItem(&SSC.SafeStrings, (void*)*ss);
	memFree((*ss)->Ptr);
	memFree(*ss);
}
void strSafeSet(nSafeString** ss, char* Content) {
	int len;
	if (*ss) strSafeDestroy(ss);
	if (!Content) return;
	len = strlen(Content);
	if (len < 1) return;
	(*ss) = memAquireOnly(sizeof(nSafeString));
	(*ss)->Ptr = memAquireOnly(sizeof(char)*(len + 1));
	strCopyFull((*ss)->Ptr, Content);
	lstAppendItem(&SSC.SafeStrings, *ss);
}


// =================================================


void tMatLoadIdentity44d(tnsMatrix44d m) {
	memset(m, 0, sizeof(tnsMatrix44d));
	m[0] = 1.0f;
	m[5] = 1.0f;
	m[10] = 1.0f;
	m[15] = 1.0f;
};

real tMatDistIdv2(real x1, real y1, real x2, real y2) {
	real x = x2 - x1, y = y2 - y1;
	return sqrt((x*x + y*y));
}
real tMatDist3dv(tnsVector3d l, tnsVector3d r) {
	real x = r[0] - l[0];
	real y = r[1] - l[1];
	real z = r[2] - l[2];
	return sqrt(x*x + y*y + z*z);
}
real tMatDist2dv(tnsVector2d l, tnsVector2d r) {
	real x = r[0] - l[0];
	real y = r[1] - l[1];
	return sqrt(x*x + y*y);
}

real tMatLength3d(tnsVector3d l) {
	return (sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]));
}
real tMatLength2d(tnsVector3d l) {
	return (sqrt(l[0] * l[0] + l[1] * l[1]));
}
void tMatNormalize3d(tnsVector3d result, tnsVector3d l) {
	real r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
	result[2] = l[2] / r;
}
void tMatNormalize2d(tnsVector2d result, tnsVector2d l) {
	real r = sqrt(l[0] * l[0] + l[1] * l[1]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
}
void tMatNormalizeSelf3d(tnsVector3d result) {
	real r = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
	result[0] /= r;
	result[1] /= r;
	result[2] /= r;
}
real tMatDot3d(tnsVector3d l, tnsVector3d r, int normalize) {
	tnsVector3d ln, rn;
	if (normalize) {
		tMatNormalize3d(ln, l); tMatNormalize3d(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1] + ln[2] * rn[2]);
	}
	return (l[0] * r[0] + l[1] * r[1] + l[2] * r[2]);
}
real tMatDot2d(tnsVector2d l, tnsVector2d r, int normalize) {
	tnsVector3d ln, rn;
	if (normalize) {
		tMatNormalize2d(ln, l); tMatNormalize2d(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1]);
	}
	return (l[0] * r[0] + l[1] * r[1]);
}
real tMatVectorCross3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[1] * r[2] - l[2] * r[1];
	result[1] = l[2] * r[0] - l[0] * r[2];
	result[2] = l[0] * r[1] - l[1] * r[0];
	return tMatLength3d(result);
}
void tMatVectorCrossOnly3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[1] * r[2] - l[2] * r[1];
	result[1] = l[2] * r[0] - l[0] * r[2];
	result[2] = l[0] * r[1] - l[1] * r[0];
}
real tMatAngleRad3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference) {
	if (PositiveReference) {
		tnsVector3d res;
		tMatVectorCross3d(res, from, to);
		if (tMatDot3d(res, PositiveReference, 1)>0)
			return acosf(tMatDot3d(from, to, 1));
		else
			return -acosf(tMatDot3d(from, to, 1));
	}
	return acosf(tMatDot3d(from, to, 1));
}
void tMatApplyRotation33d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	result[0] = mat[0] * v[0] + mat[1] * v[1] + mat[2] * v[2];
	result[1] = mat[3] * v[0] + mat[4] * v[1] + mat[5] * v[2];
	result[2] = mat[6] * v[0] + mat[7] * v[1] + mat[8] * v[2];
}
void tMatApplyRotation43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	result[0] = mat[0] * v[0] + mat[1] * v[1] + mat[2] * v[2];
	result[1] = mat[4] * v[0] + mat[5] * v[1] + mat[6] * v[2];
	result[2] = mat[8] * v[0] + mat[9] * v[1] + mat[10] * v[2];
}
void tMatApplyTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	w = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
	result[0] /= w;
	result[1] /= w;
	result[2] /= w;
}
void tMatApplyNormalTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
}
void tMatApplyTransform44d(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
}
void tMatApplyTransform43df(tnsVector4d result, tnsMatrix44d mat, tnsVector3f v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
}
void tMatApplyTransform44dTrue(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * v[3];
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * v[3];
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * v[3];
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * v[3];
}

void tMatRemoveTranslation44d(tnsMatrix44d result, tnsMatrix44d mat) {
	tMatLoadIdentity44d(result);
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
void tMatClearTranslation44d(tnsMatrix44d mat) {
	mat[3] = 0;
	mat[7] = 0;
	mat[11] = 0;
}


void tMatExtractXYZEuler44d(tnsMatrix44d mat, real* x_result, real* y_result, real* z_result) {
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
void tMatExtractLocation44d(tnsMatrix44d mat, real* x_result, real* y_result, real* z_result) {
	*x_result = mat[12];
	*y_result = mat[13];
	*z_result = mat[14];
}
void tMatExtractUniformScale44d(tnsMatrix44d mat, real* result) {
	tnsVector3d v = { mat[0],mat[1],mat[2] };
	*result = tMatLength3d(v);
}



#define L(row,col)  l[(col<<2)+row]
#define R(row,col)  r[(col<<2)+row]
#define P(row,col)  result[(col<<2)+row]

void tMatPrintMatrix44d(tnsMatrix44d l) {
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			printf("%.5f ", L(i, j));
		}
		printf("\n");
	}
}

void tMatCopyMatrix44d(tnsMatrix44d from, tnsMatrix44d to) {
	memcpy(to, from, sizeof(tnsMatrix44d));
}
void tMatMultiply44d(tnsMatrix44d result, tnsMatrix44d l, tnsMatrix44d r) {
	int i;
	for (i = 0; i < 4; i++) {
		real ai0 = L(i, 0), ai1 = L(i, 1), ai2 = L(i, 2), ai3 = L(i, 3);
		P(i, 0) = ai0 * R(0, 0) + ai1 * R(1, 0) + ai2 * R(2, 0) + ai3 * R(3, 0);
		P(i, 1) = ai0 * R(0, 1) + ai1 * R(1, 1) + ai2 * R(2, 1) + ai3 * R(3, 1);
		P(i, 2) = ai0 * R(0, 2) + ai1 * R(1, 2) + ai2 * R(2, 2) + ai3 * R(3, 2);
		P(i, 3) = ai0 * R(0, 3) + ai1 * R(1, 3) + ai2 * R(2, 3) + ai3 * R(3, 3);
	}
};
void tMatInverse44d(tnsMatrix44d inverse, tnsMatrix44d mat) {
	int i, j, k;
	double temp;
	tnsMatrix44d tempmat;
	real max;
	int maxj;

	tMatLoadIdentity44d(inverse);

	tMatCopyMatrix44d(mat, tempmat);

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
void tMatMakeTranslationMatrix44d(tnsMatrix44d mTrans, real x, real y, real z) {
	tMatLoadIdentity44d(mTrans);
	mTrans[12] = x;
	mTrans[13] = y;
	mTrans[14] = z;
}
void tMatMakePerspectiveMatrix44d(tnsMatrix44d mProjection, real fFov_rad, real fAspect, real zMin, real zMax) {
	real yMax = zMin * tanf(fFov_rad * 0.5f);
	real yMin = -yMax;
	real xMin = yMin * fAspect;
	real xMax = -xMin;

	tMatLoadIdentity44d(mProjection);

	mProjection[0] = (2.0f * zMin) / (xMax - xMin);
	mProjection[5] = (2.0f * zMin) / (yMax - yMin);
	mProjection[8] = (xMax + xMin) / (xMax - xMin);
	mProjection[9] = (yMax + yMin) / (yMax - yMin);
	mProjection[10] = -((zMax + zMin) / (zMax - zMin));
	mProjection[11] = -1.0f;
	mProjection[14] = -((2.0f * (zMax*zMin)) / (zMax - zMin));
	mProjection[15] = 0.0f;
}
void tMatMakeZTrackingMatrix44d(tnsMatrix44d mat, tnsVector3d this, tnsVector3d that, tnsVector3d up) {
	tnsVector4d fwd, l, t, rt;
	fwd[3] = l[3] = t[3] = rt[3] = 1;
	t[0] = up[0];
	t[1] = up[1];
	t[2] = up[2];
	fwd[0] = that[0] - this[0];
	fwd[1] = that[1] - this[1];
	fwd[2] = that[2] - this[2];

	tMatLoadIdentity44d(mat);

	tMatVectorCross3d(l, t, fwd);
	tMatVectorCross3d(rt, fwd, l);

	tMatNormalizeSelf3d(l);
	tMatNormalizeSelf3d(rt);
	tMatNormalizeSelf3d(fwd);

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
void tMatMakeZTrackingMatrixDelta44d(tnsMatrix44d mat, tnsVector3d delta, tnsVector3d up) {
	tnsVector4d fwd, l, t, rt;
	fwd[3] = l[3] = t[3] = rt[3] = 1;
	t[0] = up[0];
	t[1] = up[1];
	t[2] = up[2];
	fwd[0] = delta[0];
	fwd[1] = delta[1];
	fwd[2] = delta[2];

	tMatLoadIdentity44d(mat);

	tMatVectorCross3d(l, t, fwd);
	tMatVectorCross3d(rt, fwd, l);

	tMatNormalizeSelf3d(l);
	tMatNormalizeSelf3d(rt);
	tMatNormalizeSelf3d(fwd);

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
void tMatMakeOrthographicMatrix44d(tnsMatrix44d mProjection, real xMin, real xMax, real yMin, real yMax, real zMin, real zMax) {
	tMatLoadIdentity44d(mProjection);

	mProjection[0] = 2.0f / (xMax - xMin);
	mProjection[5] = 2.0f / (yMax - yMin);
	mProjection[10] = -2.0f / (zMax - zMin);
	mProjection[12] = -((xMax + xMin) / (xMax - xMin));
	mProjection[13] = -((yMax + yMin) / (yMax - yMin));
	mProjection[14] = -((zMax + zMin) / (zMax - zMin));
	mProjection[15] = 1.0f;
}
void tMatMakeRotationMatrix44d(tnsMatrix44d m, real angle_rad, real x, real y, real z)
{
	real mag, s, c;
	real xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

	s = (real)sin(angle_rad);
	c = (real)cos(angle_rad);

	mag = (real)sqrt(x*x + y*y + z*z);

	// Identity matrix
	if (mag == 0.0f) {
		tMatLoadIdentity44d(m);
		return;
	}

	// Rotation matrix is normalized
	x /= mag;
	y /= mag;
	z /= mag;

#define M(row,col)  m[col*4+row]

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
void tMatMakeRotationXMatrix44d(tnsMatrix44d m, real angle_rad) {
	tMatLoadIdentity44d(m);
	m[5] = cos(angle_rad);
	m[6] = sin(angle_rad);
	m[9] = -sin(angle_rad);
	m[10] = cos(angle_rad);
}
void tMatMakeRotationYMatrix44d(tnsMatrix44d m, real angle_rad) {
	tMatLoadIdentity44d(m);
	m[0] = cos(angle_rad);
	m[2] = -sin(angle_rad);
	m[8] = sin(angle_rad);
	m[10] = cos(angle_rad);
}
void tMatMakeRotationZMatrix44d(tnsMatrix44d m, real angle_rad) {
	tMatLoadIdentity44d(m);
	m[0] = cos(angle_rad);
	m[1] = sin(angle_rad);
	m[4] = -sin(angle_rad);
	m[5] = cos(angle_rad);
}
void tMatMakeScaleMatrix44d(tnsMatrix44d m, real x, real y, real z) {
	tMatLoadIdentity44d(m);
	m[0] = x;
	m[5] = y;
	m[10] = z;
}
void tMatMakeViewportMatrix44d(tnsMatrix44d m, real w, real h, real Far, real Near) {
	tMatLoadIdentity44d(m);
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


real tnsLinearInterpolate(real L, real R, real T) {
	return tnsLinearItp(L, R, T);
}
void tnsLinearInterpolate2dv(real* L, real* R, real T, real* Result) {
	Result[0] = tnsLinearItp(L[0], R[0], T);
	Result[1] = tnsLinearItp(L[1], R[1], T);
}
void tnsLinearInterpolate3dv(real* L, real* R, real T, real* Result) {
	Result[0] = tnsLinearItp(L[0], R[0], T);
	Result[1] = tnsLinearItp(L[1], R[1], T);
	Result[2] = tnsLinearItp(L[2], R[2], T);
}

