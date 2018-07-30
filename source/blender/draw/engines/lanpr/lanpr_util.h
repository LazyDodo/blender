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

#define NEED_STRUCTURE(a) \
	typedef struct _##a a;

#define STRUCTURE(a) \
	typedef struct _##a a; \
	struct _##a

#define DBL_TRIANGLE_LIM 1e-11
#define DBL_EDGE_LIM 1e-9


typedef struct _nListItem nListItem;
struct _nListItem {
	void *pNext;
	void *pPrev;
};

typedef struct _nListItem2 nListItem2;
struct _nListItem2 {
	void *O1;
	void *O2;
	void *pNext;
	void *pPrev;
};

typedef struct _nListItemPointer nListItemPointer;
struct _nListItemPointer {
	void *pNext;
	void *pPrev;
	void *p;
};

typedef struct _nHash256 nHash256;
struct _nHash256 {
	ListBase Entries[256];
};

typedef struct _nHash65536 nHash65536;
struct _nHash65536 {
	ListBase Entries[65536];
	//nHash256 HashHandles[256];
};

typedef struct _nHash16M nHash16M;
struct _nHash16M {
	ListBase Entries[16777216];
};

typedef struct _nSafeString nSafeString;
struct _nSafeString {
	nListItem Item;
	char *Ptr;
};

typedef struct _nSafeStringCollection nSafeStringCollection;
struct _nSafeStringCollection {
	ListBase SafeStrings;
};

typedef struct _nStringSplitor nStringSplitor;
struct _nStringSplitor {
	int NumberParts;
	ListBase parts;
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
	ListBase Lines;
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
	ListBase Pools;
};

STRUCTURE(nMemoryPoolPart)
{
	nListItem Item;
	ListBase MemoryNodes;
	ListBase FreeMemoryNodes;
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
	ListBase Pools;
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


void list_handle_empty(ListBase *h);

void list_clear_prev_next(nListItem *li);

void list_insert_item_before(ListBase *Handle, nListItem *toIns, nListItem *pivot);
void list_insert_item_after(ListBase *Handle, nListItem *toIns, nListItem *pivot);
void list_insert_segment_before(ListBase *Handle, nListItem *Begin, nListItem *End, nListItem *pivot);
void lstInsertSegmentAfter(ListBase *Handle, nListItem *Begin, nListItem *End, nListItem *pivot);
int   lstHaveItemInList(ListBase *Handle);
void *lst_get_top(ListBase *Handle);


void *list_append_pointer_only(ListBase *h, void *p);
void *list_append_pointer_sized_only(ListBase *h, void *p, int size);
void *list_push_pointer_only(ListBase *h, void *p);
void *list_push_pointer_sized_only(ListBase *h, void *p, int size);

void *list_append_pointer(ListBase *h, void *p);
void *list_append_pointer_sized(ListBase *h, void *p, int size);
void *list_push_pointer(ListBase *h, void *p);
void *list_push_pointer_sized(ListBase *h, void *p, int size);

void *list_append_pointer_static(ListBase *h, nStaticMemoryPool *smp, void *p);
void *list_append_pointer_static_sized(ListBase *h, nStaticMemoryPool *smp, void *p, int size);
void *list_push_pointer_static(ListBase *h, nStaticMemoryPool *smp, void *p);
void *list_push_pointer_static_sized(ListBase *h, nStaticMemoryPool *smp, void *p, int size);

void *list_pop_pointer_only(ListBase *h);
void list_remove_pointer_item_only(ListBase *h, nListItemPointer *lip);
void list_remove_pointer_only(ListBase *h, void *p);
void list_clear_pointer_only(ListBase *h);
void list_generate_pointer_list_only(ListBase *from1, ListBase *from2, ListBase *to);

void *list_pop_pointer(ListBase *h);
void list_remove_pointer_item(ListBase *h, nListItemPointer *lip);
void list_remove_pointer(ListBase *h, void *p);
void list_clear_pointer(ListBase *h);
void list_generate_pointer_list(ListBase *from1, ListBase *from2, ListBase *to);

void list_copy_handle(ListBase *target, ListBase *src);

void *list_append_pointer_static_pool(nStaticMemoryPool *mph, ListBase *h, void *p);
void *list_pop_pointer_no_free(ListBase *h);
void list_remove_pointer_item_no_free(ListBase *h, nListItemPointer *lip);

void list_move_up(ListBase *h, nListItem *li);
void list_move_down(ListBase *h, nListItem *li);

void lstAddElement(ListBase *hlst, void *ext);
void lstDestroyElementList(ListBase *hlst);

void mem_init_pool(nMemoryPool *mph, int NodeSize);
void mem_init_pool_small(nMemoryPool *mph, int NodeSize);
nMemoryPoolPart *mem_new_pool_part(nMemoryPool *mph);

#define memAquireOnly(a) \
	MEM_callocN(a, "NONE")

#define memAquire  memAquireOnly

void mem_free(void *Data);
void mem_destroy_pool(nMemoryPool *Handle);

nStaticMemoryPoolNode *mem_new_static_pool(nStaticMemoryPool *smp);
void *mem_static_aquire(nStaticMemoryPool *smp, int size);
void *mem_static_aquire_thread(nStaticMemoryPool *smp, int size);
void *mem_static_destroy(nStaticMemoryPool *smp);

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
void tMatApplyTransform43dfND(tnsVector4d result, tnsMatrix44d mat, tnsVector3f v);
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
void tMatMakeViewportMatrix44d(tnsMatrix44d m, real w, real h, real Far, real Near);
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

