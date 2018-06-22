#pragma once

#include <string.h>
//#include "lanpr_all.h"
#include "BLI_listbase.h"
#include "BLI_linklist.h"
#include "BLI_threads.h"

/*

   Ported from NUL4.0

   Author(s):WuYiming - xp8110@outlook.com

 */

#define _CRT_SECURE_NO_WARNINGS
#define BYTE unsigned char

typedef double real;
typedef unsigned long long u64bit;
typedef unsigned int u32bit;
typedef unsigned short u16bit;
typedef unsigned short ushort;
typedef unsigned char u8bit;
typedef char nShortBuf[16];

typedef float tnsMatrix44f[16];

typedef real tnsMatrix44d[16];
typedef real tnsVector2d[2];
typedef real tnsVector3d[3];
typedef real tnsVector4d[4];
typedef float tnsVector3f[3];
typedef float tnsVector4f[4];
typedef int tnsVector2i[2];

#define TNS_PI 3.1415926535897932384626433832795
#define deg(r) r / TNS_PI * 180.0
#define rad(d) d *TNS_PI / 180.0


#define inline __inline

#define NEED_STRUCTURE(a) \
	typedef struct _##a a;

#define STRUCTURE(a) \
	typedef struct _##a a; \
	struct _##a

#define DBL_TRIANGLE_LIM 1e-11
#define DBL_EDGE_LIM 1e-9


typedef struct _nListSingle nListSingle;
struct _nListSingle {
	void *pNext;
};

typedef struct _nListHandle nListHandle;
struct _nListHandle {
	void *pFirst;
	void *pLast;
};

typedef struct _nListWithPivot nListWithPivot;
struct _nListWithPivot {
	void *pFirst;
	void *pLast;
	void *Pivot;
};

typedef struct _nListItem nListItem;
struct _nListItem {
	void *pPrev;
	void *pNext;
};

typedef struct _nListItem2 nListItem2;
struct _nListItem2 {
	void *O1;
	void *O2;
	void *pPrev;
	void *pNext;
};

typedef struct _nListItem3 nListItem3;
struct _nListItem3 {
	void *O1;
	void *O2;
	void *O3;
	void *O4;
	void *pPrev;
	void *pNext;
};

typedef struct _nListItemPointer nListItemPointer;
struct _nListItemPointer {
	void *pPrev;
	void *pNext;
	void *p;
};

NEED_STRUCTURE(nSafeString);
STRUCTURE(nAuthorInfo)
{
	nListItem Item;
	nSafeString *Name;
	nSafeString *CopyrightString;
};
STRUCTURE(nTimeInfo)
{
	u16bit Year;//Also Used As Timer [ms] counter
	u8bit Month;
	u8bit Day;
	u8bit Hour;
	u8bit Minute;
	u8bit Second;
};

NEED_STRUCTURE(nPropContainer);

NEED_STRUCTURE(nDBInst);

typedef void (*nUserRemoveFunc)(void *, void *);//User,to be Destroyed
NEED_STRUCTURE(nProp);
typedef struct _nItemUserLinker nItemUserLinker;
struct _nItemUserLinker {
	nListItemPointer Pointer;
	nUserRemoveFunc Remove;
	nProp *Which;
	unsigned int FrameDistinguish;
};
typedef struct _nItemUsingLinker nItemUsingLinker;
struct _nItemUsingLinker {
	void *pNext;
	void *p;
};

STRUCTURE(nHyperLevel1)
{
	void *pPrev;
	void *pNext;

	void *SecondPrev;
	void *SecondNext;

	nDBInst *DBInst;

	nListHandle Users;
};

typedef struct _nElementListItem nElementListItem;
struct _nElementListItem {
	nListItem Item;
	void *Ext;
};

typedef struct _nListNonRecursiveRoot nListNonRecursiveRoot;
struct _nListNonRecursiveRoot {
	nListHandle NSItems;
};

typedef int (*nCompareFunc)(void *, void *);
typedef void (*nListDoFunc)(void *);
typedef void (*nListNonRecursiveDoFunc)(nListNonRecursiveRoot *, void *, void *);//item,custom
typedef void (*nListNonRecursiveCopyFunc)(nListNonRecursiveRoot *, void *, void *, void *);//old,new,custom
typedef void (*nListDoFuncArgp)(void *, void *);
typedef void (*nCopyListFunc)(void *, void *);
typedef void (*nListCustomDataRemover)(void *);
//typedef void(*ListMatcherFunc)(void*,void*);//gotten value,enumed curent lst item.

typedef struct _nListNonRecursiveItem nListNonRecursiveItem;
struct _nListNonRecursiveItem {
	nListItem Item;
	nListHandle handle;
	nListHandle            *ToHandle;//This Is Pointer!
	nListNonRecursiveDoFunc func;
	nListNonRecursiveCopyFunc CopyFunc;
	nListCustomDataRemover remover;
	void *CustomData;
	int bFreeList;
	int SizeEachNode;
};


typedef struct _nHash256 nHash256;
struct _nHash256 {
	nListHandle Entries[256];
};

typedef struct _nHash65536 nHash65536;
struct _nHash65536 {
	nListHandle Entries[65536];
	//nHash256 HashHandles[256];
};

typedef struct _nHash16M nHash16M;
struct _nHash16M {
	nListHandle Entries[16777216];
};

typedef struct _nSafeString nSafeString;
struct _nSafeString {
	nListItem Item;
	char *Ptr;
};

typedef struct _nSafeStringCollection nSafeStringCollection;
struct _nSafeStringCollection {
	nListHandle SafeStrings;
};

typedef struct _nStringSplitor nStringSplitor;
struct _nStringSplitor {
	int NumberParts;
	nListHandle parts;
};

typedef struct _nStringPart nStringPart;
struct _nStringPart {
	nListItem Item;
	char *Content;
	int IntValue;
	real FloatValue;
	char Type;
};

STRUCTURE(nStringLine)
{
	nListItem Item;
	char Buf[1024];
};

STRUCTURE(nStringEdit)
{
	nListHandle Lines;
	int CusorLine, CusorBefore;
	int BeginLine, BeginBefore;
	int EndLine, EndBefore;
};


#define NUL_MEMORY_POOL_1MB   1048576
#define NUL_MEMORY_POOL_128MB 134217728
#define NUL_MEMORY_POOL_256MB 268435456
#define NUL_MEMORY_POOL_512MB 536870912

STRUCTURE(nMemoryPool)
{
	nListItem Item;
	int NodeSize;
	int CountPerPool;
	nListHandle Pools;
};

STRUCTURE(nMemoryPoolPart)
{
	nListItem Item;
	nListHandle MemoryNodes;
	nListHandle FreeMemoryNodes;
	nMemoryPool *PoolRoot;
	//  <------Mem Begin Here.
};

NEED_STRUCTURE(nDBInst);

STRUCTURE(nMemoryPoolNode)
{
	nListItem Item;
	nMemoryPoolPart *InPool;
	nDBInst *DBInst;
	//  <------User Mem Begin Here
};

STRUCTURE(nStaticMemoryPoolNode)
{
	nListItem Item;
	int UsedByte;
	//  <----------- User Mem Start Here
};

STRUCTURE(nStaticMemoryPool)
{
	int EachSize;
	nListHandle Pools;
	SpinLock csMem;
};

#define CreateNew(Type) \
	MEM_callocN(sizeof(Type), "VOID")//nutCalloc(sizeof(Type),1)

#define CreateNew_Size(size) \
	nutCalloc(size, 1)

#define CreateNewBuffer(Type, Num) \
	MEM_callocN(sizeof(Type) *Num, "VOID BUFFER")//nutCalloc(sizeof(Type),Num);

#define FreeMem(ptr) \
	MEM_freeN(ptr)//nutFreeMem((&ptr))

#ifndef elif
#define elif \
	else if
#endif


void *nutCalloc(int size, int num);

void *nutCallocHyper(int size, int num);
void nutFreeMem(void **ptr);
int nutFloatCompare(real l, real r);

int nutSameAddress(void *l, void *r);


void lstEmptyDirect(nListHandle *h);

void lstPushSingle(void **Head, nListSingle *Item);
void *lstPopSingle(void **Head, nListSingle *Item);

void lstClearPrevNext(nListItem *li);
inline void lstAppendItem(nListHandle *Handle, void *Item) {

	nListItem *li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pFirst)
		Handle->pFirst = Item;

	if (Handle->pLast)
		((nListItem *)Handle->pLast)->pNext = li;

	li->pPrev = Handle->pLast;
	li->pNext = 0;
	Handle->pLast = li;

};
inline void lstPushItem(nListHandle *Handle, void *Item) {

	nListItem *li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pLast)
		Handle->pLast = Item;

	li->pNext = Handle->pFirst;

	if (Handle->pFirst)
		((nListItem *)Handle->pFirst)->pPrev = Item;

	Handle->pFirst = li;

};
inline void *lstPopItem(nListHandle *Handle) {
	void *popitem;
	nListItem *next;
	if (!Handle->pFirst)
		return 0;

	popitem = Handle->pFirst;

	next = ((nListItem *)Handle->pFirst)->pNext;
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

void  lstAppendItem2(nListHandle *Handle, void *Item);
void *lstPopItem2(nListHandle *Handle);
void  lstPushItem2(nListHandle *Handle, void *Item);
void  lstAppendItem3(nListHandle *Handle, void *Item);
void *lstPopItem3(nListHandle *Handle);
void  lstPushItem3(nListHandle *Handle, void *Item);

inline int lstRemoveItem(nListHandle *Handle, nListItem *li) {
	if (!li->pPrev && Handle->pFirst != li) return 0;

	if (!li->pPrev)
		Handle->pFirst = li->pNext;
	else
		((nListItem *)li->pPrev)->pNext = li->pNext;

	if (!li->pNext)
		Handle->pLast = li->pPrev;
	else
		((nListItem *)li->pNext)->pPrev = li->pPrev;

	li->pNext = li->pPrev = 0;
};
int lstRemoveItem2(nListHandle *Handle, nListItem2 *li);

int lstRemoveSegment(nListHandle *Handle, nListItem *Begin, nListItem *End);
void lstInsertItemBefore(nListHandle *Handle, nListItem *toIns, nListItem *pivot);
void lstInsertItemAfter(nListHandle *Handle, nListItem *toIns, nListItem *pivot);
void lstInsertSegmentBefore(nListHandle *Handle, nListItem *Begin, nListItem *End, nListItem *pivot);
void lstInsertSegmentAfter(nListHandle *Handle, nListItem *Begin, nListItem *End, nListItem *pivot);
int   lstHaveItemInList(nListHandle *Handle);
void *lstGetTop(nListHandle *Handle);

void lstPushSimpleItem(void **first, nItemUserLinker *iul);
void *lstPushItemUser(void **first, void *p);
void *lstPushItemUsing(void **first, void *p);

void *lstAppendPointerOnly(nListHandle *h, void *p);
void *lstAppendPointerSizedOnly(nListHandle *h, void *p, int size);
void *lstPushPointerOnly(nListHandle *h, void *p);
void *lstPushPointerSizedOnly(nListHandle *h, void *p, int size);

void *lstAppendPointer(nListHandle *h, void *p);
void *lstAppendPointerSized(nListHandle *h, void *p, int size);
void *lstPushPointer(nListHandle *h, void *p);
void *lstPushPointerSized(nListHandle *h, void *p, int size);

void *lstAppendPointerStatic(nListHandle *h, nStaticMemoryPool *smp, void *p);
void *lstAppendPointerStaticSized(nListHandle *h, nStaticMemoryPool *smp, void *p, int size);
void *lstPushPointerStatic(nListHandle *h, nStaticMemoryPool *smp, void *p);
void *lstPushPointerStaticSized(nListHandle *h, nStaticMemoryPool *smp, void *p, int size);

void *lstPopPointerOnly(nListHandle *h);
void lstRemovePointerItemOnly(nListHandle *h, nListItemPointer *lip);
void lstRemovePointerOnly(nListHandle *h, void *p);
void lstClearPointerOnly(nListHandle *h);
void lstGeneratePointerListOnly(nListHandle *from1, nListHandle *from2, nListHandle *to);

void *lstPopPointer(nListHandle *h);
void lstRemovePointerItem(nListHandle *h, nListItemPointer *lip);
void lstRemovePointer(nListHandle *h, void *p);
void lstClearPointer(nListHandle *h);
void lstGeneratePointerList(nListHandle *from1, nListHandle *from2, nListHandle *to);

void lstCopyHandle(nListHandle *target, nListHandle *src);

void *lstAppendPointerStaticPool(nStaticMemoryPool *mph, nListHandle *h, void *p);
void *lstPopPointerNoFree(nListHandle *h);
void lstRemovePointerItemNoFree(nListHandle *h, nListItemPointer *lip);

void lstMoveUp(nListHandle *h, nListItem *li);
void lstMoveDown(nListHandle *h, nListItem *li);

void  lstForAllItemsDo(nListDoFunc func, nListHandle *hList);
void lstForAllItemsDoLNRR(nListNonRecursiveDoFunc func, nListHandle *hList);
void lstForAllItemsDo_DirectFree(nListDoFunc func, nListHandle *hList);
void lstForAllItemsDo_arg_ptr(nListDoFuncArgp func, nListHandle *hList, void *arg);
void lstForAllItemsDo_NonRecursive_Root(nListHandle *FirstHandle, nListNonRecursiveDoFunc func, int bFreeItem, void *custom_data, nListCustomDataRemover remover);
void lstCopy_NonRecursive_Root(nListHandle *FromHandle, nListHandle *ToHandle, int SizeEachNode, nListNonRecursiveCopyFunc func, void *custom_data, nListCustomDataRemover remover);
void lstAddNonRecursiveListHandle(nListNonRecursiveRoot *root, nListHandle *newHandle, nListNonRecursiveDoFunc nrFunc, int bFreeList, void *custom_data, nListCustomDataRemover remover);
void lstAddNonRecursiveListCopier(nListNonRecursiveRoot *root, nListHandle *oldHandle, nListHandle *newHandle, int sizeEach, nListNonRecursiveCopyFunc nrCpyFunc, void *custom_data, nListCustomDataRemover remover);
void *lstFindItem(void *CmpData, nCompareFunc func, nListHandle *hList);
void lstCombineLists(nListHandle *dest, nListHandle *src);
void lstDestroyList(nListHandle *hlst);
void *lstReMatch(nListHandle *SearchHandle, nListHandle *CurrentHandle, void *ItemToFind);
typedef int (*MatcherFunc)(void *, void *);
void *lstReMatchEx(nListHandle *SearchHandle, nListHandle *CurrentHandle, void *ItemToFind, MatcherFunc func);

void lstAddElement(nListHandle *hlst, void *ext);
void lstDestroyElementList(nListHandle *hlst);


nListItem *hsh256FindItemCSTR(nHash256 *hash, nCompareFunc func, char *buckle);
unsigned char hsh256DoHashCSTR(char *buckle);


void memResetByteCount();
int memGetByteCount();
void memInitPool(nMemoryPool *mph, int NodeSize);
void memInitPoolSmall(nMemoryPool *mph, int NodeSize);
nMemoryPoolPart *memNewPoolPart(nMemoryPool *mph);

#define memAquireOnly(a) \
	MEM_callocN(a, "NONE")

#define memAquire  memAquireOnly

void memAssignDBInst(void *mem, nDBInst *DBInst);
nDBInst *memGetDBInst(void *mem);
void memFree(void *Data);
void memDestroyPool(nMemoryPool *Handle);

nStaticMemoryPoolNode *memNewStaticPool(nStaticMemoryPool *smp);
void *memStaticAquire(nStaticMemoryPool *smp, int size);
void *memStaticAquireThread(nStaticMemoryPool *smp, int size);
void *memStaticDestroy(nStaticMemoryPool *smp);

int  strGetStringTerminateBy(char *content, char terminator, char *Out);

int strHeadOfStringMatch(char *Str, char *SubStr);
int strSkipSegmet(char **pivot, char *content);
char *strgetLastSegmentSeperateBy(char *Content, char Seperator);
void strDiscardLastSegmentSeperateBy(char *Content, char Seperator);
void strDiscardSameBeginningSeperatedBy(char *s1, char *s2, char **Result1, char **Result2, char Seperator);
int strCountSegmentSeperateBy(char *Content, char Seperator);
void strMakeDifferentName(char *Target);

void strReplaceCharacter(char *Str, char Find, char Replace);
void strToUpperCase(char *Str);
void strToLowerCase(char *Str);

nStringSplitor *strSplitPath(char *path);
int strMakeInstructions(nStringSplitor **result, char *content);
nStringPart *strGetArgument(nStringSplitor *ss, char *content);
char *strGetArgumentString(nStringSplitor *ss, char *content);
int strArgumentMatch(nStringSplitor *ss, char *id, char *value);
int strDestroyStringSplitor(nStringSplitor **ss);

int strGetIntSimple(char *content);
real strGetFloatSimple(char *content);

void strConvInt_CString(int src, char *dest, int lenth);
void strConvFloat_CString(real src, char *dest, int lenth);

void strCopyFull(char *dest, char *src);
void strCopySized(char *dest, int LenthLim, char *src);
#define strAppend strcat_s
void strPrintFloatAfter(char *dest, int LenthLim, int bits, real data);
void strPrintIntAfter(char *dest, int LenthLim, int data);
void strToWideChar(wchar_t *destBuf, char *srcBuf);

int strIsTheSame(char *src, char *dest);

void strSafeDestroy(nSafeString **ss);
void strSafeSet(nSafeString **ss, char *Content);

void tMatObmatTo16d(float obmat[4][4], tnsMatrix44d out);

real tMatDistIdv2(real x1, real y1, real x2, real y2);
real tMatDist3dv(tnsVector3d l, tnsVector3d r);
real tMatDist2dv(tnsVector2d l, tnsVector2d r);

real tMatLength3d(tnsVector3d l); real tMatLength2d(tnsVector3d l);
void tMatNormalize3d(tnsVector3d result, tnsVector3d l);
void tMatNormalize3f(tnsVector3f result, tnsVector3f l);
void tMatNormalizeSelf3d(tnsVector3d result);
real tMatDot3d(tnsVector3d l, tnsVector3d r, int normalize);
real tMatDot3df(tnsVector3d l, tnsVector3f r, int normalize);
real tMatDot2d(tnsVector2d l, tnsVector2d r, int normalize);
real tMatVectorCross3d(tnsVector3d result, tnsVector3d l, tnsVector3d r);
real tMatAngleRad3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference);
void tMatApplyRotation33d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyRotation43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyNormalTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyNormalTransform43df(tnsVector3d result, tnsMatrix44d mat, tnsVector3f v);
void tMatApplyTransform44d(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v);
void tMatApplyTransform43df(tnsVector4d result, tnsMatrix44d mat, tnsVector3f v);
void tMatApplyTransform44dTrue(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v);


void tMatLoadIdentity44d(tnsMatrix44d m);
void tMatMakeOrthographicMatrix44d(tnsMatrix44d mProjection, real xMin, real xMax, real yMin, real yMax, real zMin, real zMax);
void tMatMakePerspectiveMatrix44d(tnsMatrix44d mProjection, real fFov_rad, real fAspect, real zMin, real zMax);
void tMatMakeTranslationMatrix44d(tnsMatrix44d mTrans, real x, real y, real z);
void tMatMakeRotationMatrix44d(tnsMatrix44d m, real angle_rad, real x, real y, real z);
void tMatMakeScaleMatrix44d(tnsMatrix44d m, real x, real y, real z);
void tMatMakeViewportMatrix44d(tnsMatrix44d m, real w, real h);
void tMatMultiply44d(tnsMatrix44d result, tnsMatrix44d l, tnsMatrix44d r);
void tMatInverse44d(tnsMatrix44d inverse, tnsMatrix44d mat);
void tMatMakeRotationXMatrix44d(tnsMatrix44d m, real angle_rad);
void tMatMakeRotationYMatrix44d(tnsMatrix44d m, real angle_rad);
void tMatMakeRotationZMatrix44d(tnsMatrix44d m, real angle_rad);
void tMatRemoveTranslation44d(tnsMatrix44d result, tnsMatrix44d mat);
void tMatClearTranslation44d(tnsMatrix44d mat);

real tMatAngleRad3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference);
real tMatLength3d(tnsVector3d l);
void tMatNormalize2d(tnsVector2d result, tnsVector2d l);
void tMatNormalize3d(tnsVector3d result, tnsVector3d l);
void tMatNormalizeSelf3d(tnsVector3d result);
real tMatDot3d(tnsVector3d l, tnsVector3d r, int normalize);
real tMatVectorCross3d(tnsVector3d result, tnsVector3d l, tnsVector3d r);
void tMatVectorCrossOnly3d(tnsVector3d result, tnsVector3d l, tnsVector3d r);

