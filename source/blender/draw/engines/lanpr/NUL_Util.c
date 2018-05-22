/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/
#define _CRT_SEQURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "NUL_Util.h"
#include "NUL_Interface.h"

nSafeStringCollection SSC;

extern NUL MAIN;


struct tm* nulGetFullTime() {
	time_t t = time(0);
	return localtime(&t);
}

void nulRecordTime(nTimeRecorder* tr) {
	//GetSystemTime(&tr->Time);
	//time(&tr->At);
	tr->At = clock();
}
real nulTimeElapsedSecondsf(nTimeRecorder* End, nTimeRecorder* Begin) {
	return (real)(End->At - Begin->At) / CLOCKS_PER_SEC;
}
int nulTimeElapsedMilliseconds(nTimeRecorder* End, nTimeRecorder* Begin) {
	return(
		End->Time.wMilliseconds - Begin->Time.wMilliseconds +
		(End->Time.wSecond - Begin->Time.wSecond) * 1000 +
		(End->Time.wMinute - Begin->Time.wMinute) * 60000 +
		(End->Time.wHour - Begin->Time.wHour) * 3600000
		);
}

void nulSetAuthorInfo(char* Name,char* CopyrightString){
	strSafeSet(&MAIN.Author.Name, Name);
	strSafeSet(&MAIN.Author.CopyrightString, CopyrightString);
}

void nutCreateNUID(nHyperItem* hi) {
	sprintf(hi->NUID.String, "%hd%02hd%02hd%02hd%02hd%02hd%08X", NUL_HYPER_CREATED_TIME(hi), hi);
}
void nutMakeHyperData(nHyperItem* hi) {
	struct tm *time;
	time = nulGetFullTime();
	hi->CreatedBy = &MAIN.Author;
	hi->TimeCreated.Year = time->tm_year + 1900;
	hi->TimeCreated.Month = time->tm_mon + 1;
	hi->TimeCreated.Day = time->tm_mday;
	hi->TimeCreated.Hour = time->tm_hour;
	hi->TimeCreated.Minute = time->tm_min;
	hi->TimeCreated.Second = time->tm_sec;
	memcpy(&hi->TimeModified, &hi->TimeCreated, sizeof(nTimeInfo));
	nutCreateNUID(hi);
}
void* nutCallocHyper(int size, int num){
	nHyperItem *hi = calloc(size, num);
	nutMakeHyperData(hi);
	hi->DBInst = nul_AddDBInst(hi, hi->NUID.String, hi);
	return hi;
}

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

//inline void lstAppendItem(nListHandle* Handle, void* Item){
//
//	nListItem* li = Item;
//
//	li->pNext = li->pPrev = 0;
//
//	if (!Handle->pFirst)
//		Handle->pFirst = Item;
//
//	if (Handle->pLast)
//		((nListItem*)Handle->pLast)->pNext = li;
//
//	li->pPrev = Handle->pLast;
//	li->pNext = 0;
//	Handle->pLast = li;
//
//};
//inline void lstPushItem(nListHandle* Handle, void* Item){
//
//	nListItem* li = Item;
//
//	li->pNext = li->pPrev = 0;
//
//	if (!Handle->pLast)
//		Handle->pLast = Item;
//
//	li->pNext = Handle->pFirst;
//
//	if (Handle->pFirst)
//		((nListItem*)Handle->pFirst)->pPrev = Item;
//
//	Handle->pFirst = li;
//
//};
//inline void* lstPopItem(nListHandle* Handle){
//	void* popitem;
//	nListItem* next;
//	if (!Handle->pFirst)
//		return 0;
//
//	popitem = Handle->pFirst;
//
//	next = ((nListItem*)Handle->pFirst)->pNext;
//	if (!next){
//		Handle->pFirst = 0;
//		Handle->pLast = 0;
//	}
//	else{
//		Handle->pFirst = next;
//		if (next)
//			next->pPrev = 0;
//	};
//
//	return popitem;
//};
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
	lstRemoveItem(h, lip);
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
	lstRemoveItem(h, lip);
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
	lstRemoveItem(h, lip);
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
		lstRemoveItem(hlst, eli);
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


void memResetByteCount() {
	MAIN.ByteCount = 0;
	MAIN.TotalByteCount = 0;
}
int memGetByteCount() {
	return MAIN.ByteCount;
}
int memGetTotalByteCount() {
	return MAIN.TotalByteCount;
}
void memInitPool(nMemoryPool* mph, int NodeSize) {
	int Count = 4096;

	if (mph->NodeSize) return;

	mph->NodeSize = NodeSize;
	mph->CountPerPool = Count;

	return mph;
}
nMemoryPool* memInitGlobalPool(int NodeSize) {
	nMemoryPool* mph = calloc(1, sizeof(nMemoryPool));
	memInitPool(mph, NodeSize);

	u8bit Buckle = NodeSize;
	lstAppendItem(&MAIN.GlobalMemPool.Entries[Buckle], mph);

	return mph;
}
void memInitPoolSmall(nMemoryPool* mph, int NodeSize) {
	int Count = 16;

	if (mph->NodeSize) return;

	mph->NodeSize = NodeSize;
	mph->CountPerPool = Count;

	return mph;
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

	mpn = ((BYTE*)BeginMem);
	mpn->InPool = mp;

	for (i = 1; i < NodeCount; i++) {
		mpn = ((BYTE*)BeginMem) + RealNodeSize*i;
		mpn->InPool = mp;
		lstAppendItem(&mp->FreeMemoryNodes, mpn);
	}
	lstPushItem(&mph->Pools, mp);

	MAIN.TotalByteCount += TotalSize;

	return mp;
}
inline void* memAquireH(nMemoryPool* Handle) {
	nMemoryPoolPart* mp = Handle->Pools.pFirst;
	nMemoryPoolNode* mpn;
	
	if (!mp || !mp->FreeMemoryNodes.pFirst) {
		//if (!mp) {
			mp = memNewPoolPart(Handle);
		//	break;
		//}
		//mp = mp->Item.pNext;
	}

	if (!mp) return 0;

	mpn = mp->FreeMemoryNodes.pFirst;

	lstRemoveItem(&mp->FreeMemoryNodes, mpn);
	lstAppendItem(&mp->MemoryNodes, mpn);

	MAIN.ByteCount += Handle->NodeSize;
	
	return (((BYTE*)mpn) + sizeof(nMemoryPoolNode));
}
void* memAquireOnly(int Size) {
	nMemoryPool* mp;
	u8bit Buckle = Size;

	mp = MAIN.GlobalMemPool.Entries[Buckle].pFirst;

	while (mp && mp->NodeSize != Size) mp = mp->Item.pNext;

	if (!mp) mp = memInitGlobalPool(Size);

	return memAquireH(mp);
}
void* memAquire(int Size) {
	char* mem = memAquireOnly(Size);
	nMemoryPoolNode* mpn = mem - sizeof(nMemoryPoolNode);
	mpn->DBInst = nul_AddDBInst(mem, 0, mem);
	return mem;
}
void* memAquireHyper(int Size) {
	char* mem = memAquireOnly(Size);
	nMemoryPoolNode* mpn = mem - sizeof(nMemoryPoolNode);
	nutMakeHyperData(mem);
	((nHyperItem*)mem)->DBInst = mpn->DBInst;
	mpn->DBInst = nul_AddDBInst(mem, ((nHyperItem*)mem)->NUID.String, mem);
	return mem;
}
void* memAquireHyper1(int Size) {
	char* mem = memAquireOnly(Size);
	nMemoryPoolNode* mpn = mem - sizeof(nMemoryPoolNode);
	mpn->DBInst = nul_AddDBInst(mem, 0, mem);
	((nHyperLevel1*)mem)->DBInst = mpn->DBInst;
	return mem;
}
nDBInst* memGetDBInst(void* mem) {
	nMemoryPoolNode* mpn = ((BYTE*)mem) - sizeof(nMemoryPoolNode);
	return mpn->DBInst;
}
void memAssignDBInst(void* mem, nDBInst* DBInst) {
	nMemoryPoolNode* mpn = ((BYTE*)mem) - sizeof(nMemoryPoolNode);
	mpn->DBInst = DBInst;
}
void memFree(void* Data) {
	if (!Data)return 0;
	nMemoryPoolNode* mpn = ((BYTE*)Data) - sizeof(nMemoryPoolNode);
	nMemoryPoolPart* mp = mpn->InPool;
	nMemoryPool* mph = mp->PoolRoot;
	lstRemoveItem(&mp->MemoryNodes, mpn);
	lstAppendItem(&mp->FreeMemoryNodes, mpn);
	memset(Data, 0, mph->NodeSize);
	MAIN.ByteCount -= mph->NodeSize;
	if (!mp->MemoryNodes.pFirst) {
		lstRemoveItem(&mph->Pools, mp);
		FreeMem(mp);
		MAIN.TotalByteCount -= ((mph->NodeSize + sizeof(nMemoryPoolNode))*mph->CountPerPool+sizeof(nMemoryPoolPart));
	}
	//if (!mph->Pools.pFirst) {
	//	mph->CountPerPool = 0;
	//	mph->NodeSize = 0;
	//}
}
void memDestroyPool(nMemoryPool* Handle) {
	nMemoryPool* mp;
	while ((mp = lstPopItem(Handle))) {
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

	EnterCriticalSection(&smp->csMem);

	if (!smpn || (smpn->UsedByte + size) > NUL_MEMORY_POOL_128MB)
		smpn = memNewStaticPool(smp);

	ret = ((BYTE*)smpn) + smpn->UsedByte;

	smpn->UsedByte += size;

	LeaveCriticalSection(&smp->csMem);

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


//=======================================================================[ AVL Tree ]

#define AVL_HEIGHT(node)\
(node?node->Height:0)

nAVLNodeReal64* nAVRTreeRightRotate(nAVLNodeReal64* y) {
	nAVLNodeReal64* x = y->Smaller;
	nAVLNodeReal64* T2 = x->Greater;

	x->Greater = y;
	y->Smaller = T2;

	x->Height = max(AVL_HEIGHT(x->Smaller), AVL_HEIGHT(x->Greater)) + 1;
	y->Height = max(AVL_HEIGHT(y->Smaller), AVL_HEIGHT(y->Greater)) + 1;

	return x;
}
nAVLNodeReal64* nAVRTreeLeftRotate(nAVLNodeReal64* x) {
	nAVLNodeReal64* y = x->Greater;
	nAVLNodeReal64* T2 = y->Smaller;

	y->Smaller = x;
	x->Greater = T2;

	x->Height = max(AVL_HEIGHT(x->Smaller), AVL_HEIGHT(x->Greater)) + 1;
	y->Height = max(AVL_HEIGHT(y->Smaller), AVL_HEIGHT(y->Greater)) + 1;

	return y;
}

nAVLNodeReal64* nAVLTreeInsertNode(nAVLNodeReal64* Root, nAVLNodeReal64* Node) {
	if (!Root) return Node;

	if (Node->Value < Root->Value) {
		Root->Smaller = nAVLTreeInsertNode(Root->Smaller, Node);
		if (Root->Smaller == Root->Smaller) return Root;
	}elif(Node->Value > Root->Value) {
		Root->Greater = nAVLTreeInsertNode(Root->Greater, Node);
		if (Root->Greater == Root->Greater) return Root;
	}else return Root;

	int HS = AVL_HEIGHT(Root->Smaller);
	int HG = AVL_HEIGHT(Root->Greater);

	Root->Height = max(HS, HG) + 1;

	int BalanceStatus = HS - HG;

	if (BalanceStatus > 1 && Node->Value > Root->Smaller->Value) {
		Root->Smaller = nAVRTreeLeftRotate(Root->Smaller);
		return nAVRTreeRightRotate(Root);
	}

	if (BalanceStatus < -1 && Node->Value < Root->Greater->Value) {
		Root->Greater = nAVRTreeRightRotate(Root->Greater);
		return nAVRTreeLeftRotate(Root);
	}

	if (BalanceStatus > 1 && Node->Value <= Root->Smaller->Value) {
		return nAVRTreeRightRotate(Root);
	}

	if (BalanceStatus < -1 && Node->Value >= Root->Greater->Value) {
		return nAVRTreeLeftRotate(Root);
	}

	return Root;
}
nAVLNodeReal64* nAVLTreeLookupJustSmaller(nAVLNodeReal64* Root, real Value) {
	nAVLNodeReal64* anr = Root;
	while (anr) {
		if (anr->Value >= Value) {
			if (!anr->Smaller) return 0;
			if (anr->Smaller->Value <= Value) {
				anr = anr->Smaller;
				continue;
			}
		}else {
			if (!anr->Greater) return anr;
			else {
				if (anr->Greater->Value <= Value) {
					anr = anr->Greater;
					continue;
				}else {
					if (anr->Greater->Smaller) { anr = anr->Greater->Smaller; continue; }
					else return anr;
				}
			}
		}
	}
	return 0;
}
nAVLNodeReal64* nAVLTreeLookupJustGreater(nAVLNodeReal64* Root, real Value) {
	nAVLNodeReal64* anr = Root;
	while (anr) {
		if (anr->Value <= Value) {
			if (!anr->Greater) return 0;
			if (anr->Greater->Value >= Value) {
				anr = anr->Greater;
				continue;
			}
		}
		else {
			if (!anr->Smaller) return anr;
			else {
				if (anr->Smaller->Value >= Value) {
					anr = anr->Smaller;
					continue;
				}
				else {
					if (anr->Smaller->Greater) { anr = anr->Smaller->Greater; continue; }
					else return anr;
				}
			}
		}
	}
	return 0;
}
nAVLNodeReal64* nAVLTreeInsert(nAVLTreeReal64* Tree, real Value) {
	nAVLNodeReal64* anr = memAquireH(&Tree->MemoryPool,sizeof(nAVLNodeReal64));
	nAVLNodeReal64* Result;
	anr->Value = Value;
	anr->Height = 0;

	Result = nAVLTreeInsertNode(Tree->Root, anr);
		
	Tree->Root = Result;
	Tree->ItemCount++;

	return anr;
};

void nAVLTreePrintPreOrderRecursive(nAVLNodeReal64* Root) {
	if (!Root) return;
	printf("%lf ", Root->Value);
	nAVLTreePrintPreOrderRecursive(Root->Smaller);
	nAVLTreePrintPreOrderRecursive(Root->Greater);
}
void nAVLTreePrintPreOrder(nAVLTreeReal64* Tree) {
	nAVLTreePrintPreOrderRecursive(Tree->Root);
}
void nAVLTreeMakeListUpOrderRecursive(nAVLNodeReal64* Root, nListHandle* Result) {
	if (!Root) return;
	nAVLTreeMakeListUpOrderRecursive(Root->Smaller,Result);
	lstAppendPointer(Result, Root->Pointer);
	nAVLTreeMakeListUpOrderRecursive(Root->Greater,Result);
}
void nAVLTreeMakeListUpOrder(nAVLTreeReal64* Tree, nListHandle* Result) {
	nAVLTreePrintPreOrderRecursive(Tree->Root);
}

void nAVLTreeDestroyRecursive(nAVLNodeReal64* Root) {
	if (!Root) return;
	nAVLTreeDestroyRecursive(Root->Smaller);
	nAVLTreeDestroyRecursive(Root->Greater);
}
void nAVLTreeDestroy(nAVLTreeReal64* Tree) {
	nAVLTreeDestroyRecursive(Tree->Root);
	memDestroyPool(&Tree->MemoryPool);
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
	if (!pivot || !(*pivot) || !(*(*pivot)) || !content) return;

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
void strToWideChar(wchar_t* destBuf,char* srcBuf) {
	int len = strlen(srcBuf);
	MultiByteToWideChar(CP_ACP, 0, srcBuf, strlen(srcBuf), destBuf, 128);
	destBuf[len] = '\0';
}

int  strIsTheSame(char* src, char*dest) {
	return (src && !strcmp(src, dest));
}


void strSafeDestroy(nSafeString** ss) {
	if (!*ss)
		return;
	lstRemoveItem(&SSC.SafeStrings, *ss);
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



nStringEdit* strBeginEdit(char* FullStr) {
	char* p = FullStr;
	nStringEdit* se = CreateNew(nStringEdit);
	if (FullStr && FullStr[0]) {
		while ((*p)) {
			nStringLine* sl = CreateNew(nStringLine);
			p += (strGetStringTerminateBy(p, '\n', sl->Buf) + 1);
			lstAppendItem(&se->Lines, sl);
		}
	}
	if (!se->Lines.pFirst) {
		nStringLine* sl = CreateNew(nStringLine);
		lstAppendItem(&se->Lines, sl);
	}
	se->BeginLine = -1;
	se->BeginBefore = -1;
	se->EndLine = -1;
	se->EndBefore = -1;
	se->CusorLine = 0;
	se->CusorBefore = strlen(((nStringLine*)se->Lines.pFirst)->Buf);
	return se;
}
void strEndEdit(nStringEdit* se,char* Result) {
	char* p;
	nStringLine* sl, *NextSl;
	if (!se) return;
	sl = se->Lines.pFirst;
	while (sl) {
		NextSl = sl->Item.pNext;
		strcat(Result, sl->Buf);
		lstRemoveItem(&se->Lines,sl);
		FreeMem(sl);
		sl = NextSl;
	}
	FreeMem(se);
}
void strRemoveLine(nStringEdit* se, nStringLine* sl) {
	lstRemoveItem(&se->Lines, sl);
	FreeMem(sl);
}
void strRemoveLineI(nStringEdit* se, int LineIndex) {
	int i = 0;
	nStringLine* sl = se->Lines.pFirst, *NextSl;
	while (sl) {
		NextSl = sl->Item.pNext;
		if (i == LineIndex) {
			strRemoveLine(se, sl);
			break;
		}
		i++;
		sl = NextSl;
	}
}
void strSetCusor(nStringEdit* se, int LineIndex, int BeforeIndex) {
	int maxbefore;
	if (!se) return;

	maxbefore = strlen(strGetCursorLine(se)->Buf) + 1;

	BeforeIndex = BeforeIndex < 0 ? 0 : BeforeIndex>maxbefore ? maxbefore : BeforeIndex;

	se->CusorBefore = BeforeIndex;
	se->CusorLine = LineIndex;
	se->BeginLine = -1;
	se->BeginBefore = -1;
	se->EndLine = -1;
	se->EndBefore = -1;
}
void strMoveCusor(nStringEdit* se, int Left) {
	int maxbefore;
	int BeforeIndex;
	int width=1;
	nStringLine* sl;
	if (!se) return;

	sl = strGetCursorLine(se);

	maxbefore = strlen(sl->Buf) + 1;

	BeforeIndex = se->CusorBefore - (Left ? 1 : -1);
	if (Left) {
		if (BeforeIndex && sl->Buf[BeforeIndex]<0 && sl->Buf[BeforeIndex - 1] < 0)
			BeforeIndex--;
	}else{
		if (sl->Buf[BeforeIndex - 1] < 0)
			BeforeIndex++;
	}

	BeforeIndex = BeforeIndex < 0 ? 0 : BeforeIndex>maxbefore ? maxbefore : BeforeIndex;

	se->CusorBefore = BeforeIndex;
	se->BeginLine = -1;
	se->BeginBefore = -1;
	se->EndLine = -1;
	se->EndBefore = -1;
}

void strSelect(nStringEdit* se, int BeginLine, int BeginBefore, int EndLine, int EndBefore) {
	if (!se) return;
	se->BeginLine = BeginLine;
	se->BeginBefore = BeginBefore;
	se->EndLine = EndLine;
	se->EndBefore = EndBefore;
	se->CusorBefore = -1;
	se->CusorLine = -1;
}
void strSelectLineAll(nStringEdit* se) {
	if (!se) return;
	nStringLine* sl;
	int len;
	if (se->CusorLine == -1) sl = strGetBeginLine(se);
	else sl = strGetCursorLine(se);
	len = strlen(sl->Buf);

	se->EndBefore = len;
	se->BeginBefore = 0;

	se->CusorBefore = -1;
	se->CusorLine = -1;
}
void strDeselectAll(nStringEdit* se) {
	if (!se) return;
	nStringLine* sl;
	int len;
	if (se->CusorLine == -1) sl = strGetBeginLine(se);
	else sl = strGetCursorLine(se);
	len = strlen(sl->Buf);

	se->EndBefore = -1;
	se->BeginBefore = -1;

	se->CusorBefore = len;
	se->CusorLine = -1;
}
void strPanFoward(char* str, int Before, int Offset) {
	int len = strlen(str);
	int i = len + 1;
	for (i; i >= Before; i--) {
		str[i + Offset] = str[i];
	}
}
void strSquishBackward(char* str, int Before,int EndBefore) {
	int len = strlen(str);
	int i = Before;
	int Offset = Before - EndBefore;
	if (Before <= 0) return;
	for (i; i <=len; i++) {
		str[i - Offset] = str[i];
	}
}
void strClearSelection(nStringEdit* se) {
	//if (se->EndLine == -1) return;
	if (se->BeginLine != se->EndLine) {
		int i = 0;
		nStringLine* sl = se->Lines.pFirst, *NextSl;
		while (sl) {
			NextSl = sl->Item.pNext;
			if (i == se->BeginLine) {
				sl->Buf[i] = '\n';
				sl->Buf[i + 1] = '\0';
			}
			else if (i > se->BeginLine && i < se->EndLine) {
				strRemoveLine(se, sl);
			}
			else if (i == se->EndLine) {
				strSquishBackward(sl->Buf, se->EndBefore, 0);
				se->CusorLine = i;
				se->CusorBefore = 0;
				se->BeginLine = -1;
				se->BeginBefore = -1;
				se->EndLine = -1;
				se->EndBefore = -1;
			}
			if (i > se->EndLine) break;
			i++;
			sl = NextSl;
		}
	}else{
		int i = 0;
		nStringLine* sl = se->Lines.pFirst, *NextSl;
		while (sl) {
			NextSl = sl->Item.pNext;
			//if (i == se->EndLine) {
				strSquishBackward(sl->Buf, se->EndBefore, se->BeginBefore);
				se->CusorLine = i;
				se->CusorBefore = se->BeginBefore;
				se->BeginLine = -1;
				se->BeginBefore = -1;
				se->EndLine = -1;
				se->EndBefore = -1;
				break;
			//}
			i++;
			sl = NextSl;
		}
	}
}
nStringLine* strGetCursorLine(nStringEdit* se) {
	if (!se || se->CusorBefore <= -1) return  se->Lines.pFirst;
	int i = 0;
	nStringLine* sl = se->Lines.pFirst, *NextSl;
	while (sl) {
		NextSl = sl->Item.pNext;
		if (i == se->CusorLine) {
			return sl;
		}
		i++;
		sl = NextSl;
	}
	return se->Lines.pFirst;
}
nStringLine* strGetBeginLine(nStringEdit* se) {
	if (!se || se->BeginLine <= -1) return  se->Lines.pFirst;
	int i = 0;
	nStringLine* sl = se->Lines.pFirst, *NextSl;
	while (sl) {
		NextSl = sl->Item.pNext;
		if (i == se->BeginLine) {
			return sl;
		}
		i++;
		sl = NextSl;
	}
	return se->Lines.pFirst;
}
void strInsertChar(nStringEdit* se,char a) {
	nStringLine* sl;
	if (se->CusorBefore == -1) {
		strClearSelection(se);
	}
	sl = strGetCursorLine(se);
	strPanFoward(sl->Buf, se->CusorBefore, 1);
	sl->Buf[se->CusorBefore] = a;
	se->CusorBefore += 1;
}
void strBackspace(nStringEdit* se) {
	nStringLine* sl;
	int width = 1;
	if (se->CusorBefore == -1) {
		strClearSelection(se);
	}else {
		nStringLine* sl;
		sl = strGetCursorLine(se);
		if (se->CusorBefore > 1 && sl->Buf[se->CusorBefore - 2] < 0)width = 2;
		strSquishBackward(sl->Buf, se->CusorBefore, se->CusorBefore-width);
		se->CusorBefore -= width;
		if (se->CusorBefore <= -1) se->CusorBefore = 0;
	}
}


//======================================================[ translation ]

void transNewLanguage(const char* LanguageID) {
	nTranslationNode* tn = CreateNew(nTranslationNode);
	strSafeSet(&tn->LanguageName, LanguageID);

	lstAppendItem(&MAIN.Translation.Languages, tn);
	
	MAIN.Translation.CurrentLanguage = tn;
}
void transSetLanguage(const char* LanguageID) {
	nTranslationNode* tn;

	if (!LanguageID) {
		MAIN.Translation.CurrentLanguage = 0;
		return;
	}

	for (tn = MAIN.Translation.Languages.pFirst; tn; tn = tn->Item.pNext) {
		if (!strcmp(tn->LanguageName->Ptr, LanguageID)) {
			MAIN.Translation.CurrentLanguage = tn;
			break;
		}
	}
}
void transDumpMissMatchRecord(const char* filename) {
	nTranslationMatch* tm;
	nListHandle* lst;
	int i;

	FILE* f = fopen(filename, "w");
	if (!f) return;

	for (i = 0; i < 256; i++) {
		lst = &MAIN.Translation.MisMatches.Entries[i];
		for (tm = lst->pFirst; tm; tm = tm->Item.pNext) {
			fprintf(f, "transNewEntry(\"%s\",\"\");\n",tm->Target);
		}
	}

	fclose(f);
}
int IsThisTranslationMatch(nTranslationMatch* tm, char* p) {
	return (tm->Target && !strcmp(tm->Target, p));
}
void transNewEntry(const char* Target, const char* replacement) {
	nTranslationMatch* tm = memAquireOnly(sizeof(nTranslationMatch));
	tm->Target = Target;
	tm->Replacement = replacement;
	hsh256InsertItemCSTR(&MAIN.Translation.CurrentLanguage->Matches, tm, Target);
}
void transNewMissEntry(const char* Target) {
	if (!hsh256FindItemCSTR(&MAIN.Translation.MisMatches, IsThisTranslationMatch, Target)) {
		nTranslationMatch* tm = memAquireOnly(sizeof(nTranslationMatch));
		tm->Target = Target;
		hsh256InsertItemCSTR(&MAIN.Translation.MisMatches, tm, Target);
	}
}
char* transLate(char* Target) {
	if (!MAIN.Translation.CurrentLanguage || !MAIN.Translation.EnableTranslation) return Target;
	nTranslationMatch* tm = hsh256FindItemCSTR(&MAIN.Translation.CurrentLanguage->Matches, IsThisTranslationMatch, Target);
	if (!tm) {
		transNewMissEntry(Target);
		return Target;
	}
	return tm->Replacement;
}
void transState(void* UNUSED, int val) {
	if (val) MAIN.Translation.EnableTranslation = 1;
	else MAIN.Translation.EnableTranslation = 0;

	nulRedrawCurrentWindow();
}



//=========================================================[ Internet ]

void nulOpenInternetLink(char* link) {
	HKEY hkRoot, hSubKey;
	char ValueName[256];
	unsigned char DataValue[256];
	unsigned long cbValueName = 256;
	unsigned long cbDataValue = 256;
	char ShellChar[512];
	DWORD dwType;

	ShellExecute(0, "open", link, 0, 0, SW_SHOWNORMAL);

	return;
}


//===========================================================[ deprecated ]

void nulSendPanic(char* message){
	MessageBox(0, \
		message  
		, "SYSTEM_ERROR", 0);                                                       
	exit(0);
}

//===========================================================[ file util ]

char* txtReadFileAsString(char* FileName) {
	FILE* f = fopen(FileName, "r");
	int length;
	fseek(f, 0, SEEK_END);
	length = ftell(f);
	fseek(f, 0, SEEK_SET);

	char* buffer = CreateNewBuffer(char, length);

	fread(buffer, sizeof(char), length, f);

	fclose(f);

	return buffer;
}