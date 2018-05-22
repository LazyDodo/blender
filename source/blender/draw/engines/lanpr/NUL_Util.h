#pragma once

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/

#define _CRT_SECURE_NO_WARNINGS


#include "GL\glew.h"
#include "GL\wglew.h"
#include "GL\GL.h"
#include <Windows.h>
#include <shellapi.h>
#include "ft2build.h"
#include "freetype/freetype.h"
#include "NUL_IconDefine.h"

#define inline __inline

#define NEED_STRUCTURE(a)\
typedef struct _##a a;

#define STRUCTURE(a)\
typedef struct _##a a;\
struct _##a

#define DBL_TRIANGLE_LIM 1e-11
#define DBL_EDGE_LIM 1e-9


#define NUL_HYPER_CREATED_TIME(hi)\
hi->TimeCreated.Year,hi->TimeCreated.Month,hi->TimeCreated.Day,hi->TimeCreated.Hour,hi->TimeCreated.Minute,hi->TimeCreated.Second


typedef double real;
typedef unsigned long long u64bit;
typedef unsigned int       u32bit;
typedef unsigned short     u16bit;
typedef unsigned short     ushort;
typedef unsigned char      u8bit;
typedef char nShortBuf[16];

typedef struct _nListSingle nListSingle;
struct _nListSingle {
	void* pNext;
};

typedef struct _nListHandle nListHandle;
struct _nListHandle {
	void* pFirst;
	void* pLast;
};

typedef struct _nListWithPivot nListWithPivot;
struct _nListWithPivot {
	void* pFirst;
	void* pLast;
	void* Pivot;
};

typedef struct _nListItem nListItem;
struct _nListItem {
	void* pPrev;
	void* pNext;
};

typedef struct _nListItem2 nListItem2;
struct _nListItem2 {
	void* O1;
	void* O2;
	void* pPrev;
	void* pNext;
};

typedef struct _nListItem3 nListItem3;
struct _nListItem3 {
	void* O1;
	void* O2;
	void* O3;
	void* O4;
	void* pPrev;
	void* pNext;
};


NEED_STRUCTURE(nSafeString);
STRUCTURE(nAuthorInfo) {
	nListItem    Item;
	nSafeString* Name;
	nSafeString* CopyrightString;
};
STRUCTURE(nTimeInfo) {
	u16bit Year;//Also Used As Timer [ms] counter
	u8bit  Month;
	u8bit  Day;
	u8bit  Hour;
	u8bit  Minute;
	u8bit  Second;
};

NEED_STRUCTURE(nPropContainer);

typedef struct _nUID nUID;
struct _nUID {
	char String[32];//a simplified uuid, example: 0E3F9BA4802FDDC2-20160601123546 [\0]
};

NEED_STRUCTURE(nDBInst);

typedef struct _nHyperItem nHyperItem;
struct _nHyperItem {
	void*            pPrev;
	void*            pNext;

	void*            SecondPrev;
	void*            SecondNext;

	nDBInst*         DBInst;

	nListHandle      Users;

	//================= Strictly satisfy entries above,
	//================= then the structure is HyperSafeLevel[1]

	char             SaveNewVersion;//not necessarily used

	int              Version;

	char             Linked;

	nUID             NUID;
	nUID             OldNUID;//unused currently

	nAuthorInfo*     CreatedBy;
	nTimeInfo        TimeCreated;
	nTimeInfo        TimeModified;
	nSafeString*     Description;

	//================= All satisfied,
	//================= then the structure is HyperSafeLevel[2]
};

typedef struct _nListItemPointer nListItemPointer;
struct _nListItemPointer {
	void* pPrev;
	void* pNext;
	void* p;
};

typedef void(*nUserRemoveFunc)(void*,void*);//User,to be Destroyed
NEED_STRUCTURE(nProp);
typedef struct _nItemUserLinker nItemUserLinker;
struct _nItemUserLinker {
	nListItemPointer Pointer;
	nUserRemoveFunc  Remove;
	nProp*           Which;
	unsigned int     FrameDistinguish;
};
typedef struct _nItemUsingLinker nItemUsingLinker;
struct _nItemUsingLinker {
	void*            pNext;
	void*            p;
};

STRUCTURE(nHyperLevel1) {
	void*            pPrev;
	void*            pNext;

	void*            SecondPrev;
	void*            SecondNext;

	nDBInst*         DBInst;

	nListHandle      Users;
};

typedef struct _nElementListItem nElementListItem;
struct _nElementListItem {
	nListItem Item;
	void*     Ext;
};

typedef struct _nListNonRecursiveRoot nListNonRecursiveRoot;
struct _nListNonRecursiveRoot {
	nListHandle NSItems;
};

typedef int(*nCompareFunc)(void*, void*);
typedef void(*nListDoFunc)(void*);
typedef void(*nListNonRecursiveDoFunc)(nListNonRecursiveRoot*, void*, void*);//item,custom
typedef void(*nListNonRecursiveCopyFunc)(nListNonRecursiveRoot*, void*, void*, void*);//old,new,custom
typedef void(*nListDoFuncArgp)(void*, void*);
typedef void(*nCopyListFunc)(void*, void*);
typedef void(*nListCustomDataRemover)(void*);
//typedef void(*ListMatcherFunc)(void*,void*);//gotten value,enumed curent lst item.

typedef struct _nListNonRecursiveItem nListNonRecursiveItem;
struct _nListNonRecursiveItem {
	nListItem               Item;
	nListHandle             handle;
	nListHandle            *ToHandle;//This Is Pointer!
	nListNonRecursiveDoFunc func;
	nListNonRecursiveCopyFunc CopyFunc;
	nListCustomDataRemover    remover;
	void*                  CustomData;
	int                    bFreeList;
	int                    SizeEachNode;
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
	char*     Ptr;
};

typedef struct _nSafeStringCollection nSafeStringCollection;
struct _nSafeStringCollection {
	nListHandle SafeStrings;
};

typedef struct _nStringSplitor nStringSplitor;
struct _nStringSplitor {
	int         NumberParts;
	nListHandle parts;
};

typedef struct _nStringPart nStringPart;
struct _nStringPart {
	nListItem  Item;
	char*      Content;
	int        IntValue;
	real      FloatValue;
	char       Type;
};

STRUCTURE(nStringLine) {
	nListItem Item;
	char      Buf[1024];
};

STRUCTURE(nStringEdit) {
	nListHandle Lines;
	int CusorLine, CusorBefore;
	int BeginLine, BeginBefore;
	int EndLine, EndBefore;
};


#define NUL_MEMORY_POOL_1MB   1048576
#define NUL_MEMORY_POOL_128MB 134217728
#define NUL_MEMORY_POOL_256MB 268435456
#define NUL_MEMORY_POOL_512MB 536870912

STRUCTURE(nMemoryPool) {
	nListItem   Item;
	int         NodeSize;
	int         CountPerPool;
	nListHandle Pools;
};

STRUCTURE(nMemoryPoolPart) {
	nListItem   Item;
	nListHandle MemoryNodes;
	nListHandle FreeMemoryNodes;
	nMemoryPool* PoolRoot;
	//  <------Mem Begin Here.
};

NEED_STRUCTURE(nDBInst);

STRUCTURE(nMemoryPoolNode) {
	nListItem        Item;
	nMemoryPoolPart* InPool;
	nDBInst*         DBInst;
	//  <------User Mem Begin Here
};

STRUCTURE(nStaticMemoryPoolNode) {
	nListItem Item;
	int       UsedByte;
	//  <----------- User Mem Start Here
};

STRUCTURE(nStaticMemoryPool) {
	int         EachSize;
	nListHandle Pools;
	CRITICAL_SECTION csMem;
};



STRUCTURE(nAVLNodeReal64) {
	nAVLNodeReal64* Parent;
	u64bit          Index;
	real            Value;
	//real            SmallestValue;
	//real            GreatestValue;
	nAVLNodeReal64* Smaller;
	nAVLNodeReal64* Greater;
	char            Height;
	void*           Pointer;
};	

STRUCTURE(nAVLTreeReal64) {
	nAVLNodeReal64*   Root;
	u64bit            ItemCount;
	nMemoryPool       MemoryPool;
};


STRUCTURE(nTimeRecorder) {
	SYSTEMTIME Time;
	time_t     At;
};

STRUCTURE(nTranslationNode) {
	nHyperItem   Item;
	nSafeString* LanguageName;
	nHash256     Matches;
};

STRUCTURE(nTranslation) {
	int               EnableTranslation;
	nListHandle       Languages;
	nTranslationNode* CurrentLanguage;
	nHash256          MisMatches;
};

STRUCTURE(nTranslationMatch) {
	nListItem Item;
	const char* Target;
	const char* Replacement;
};


char* txtReadFileAsString(char* FileName);


void nulSendPanic(char* message);

#define SEND_PANIC_ERROR(message)\
    nulSendPanic(message)

#define SEND_NORMAL_ERROR(message)                                                \
    MessageBox(0,\
    "An error message arrived:\n_______________________________\n"message         \
    "\n________________________________\nThis error is not caused by GUI system,\n\
Please check your own operation code",                                            \
"USER_CAUSED_ERROR",0);

#define SEND_USER_MESSAGE(message)                                                \
    MessageBox(0,message,"Message here:",0)

#define CreateNew(Type)\
    nutCalloc(sizeof(Type),1)

#define HyperNew(Type)\
    nutCallocHyper(sizeof(Type),1)

#define HyperNew_Size(Size)\
    nutCallocHyper(Size,1)


#define CreateNew_Size(size)\
    nutCalloc(size,1)

#define CreateNewBuffer(Type,Num)\
    nutCalloc(sizeof(Type),Num);

#define FreeMem(ptr)\
    nutFreeMem((&ptr))

#define elif\
	    else if


#define NUL_UNAVAILABLE_NAME "- Unknown -"

struct tm* nulGetFullTime();

void nulRecordTime(nTimeRecorder* tr);
real nulTimeElapsedSecondsf(nTimeRecorder* End, nTimeRecorder* Begin);
int nulTimeElapsedMilliseconds(nTimeRecorder* End, nTimeRecorder* Begin);

void nulSetAuthorInfo(char* Name, char* CopyrightString);

void nutCreateNUID(nHyperItem* hi);

void* nutCalloc(int size, int num);
void nutMakeHyperData(nHyperItem* hi);
void* nutCallocHyper(int size, int num);
void nutFreeMem(void** ptr);
int nutFloatCompare(real l, real r);

int nutSameAddress(void* l, void* r);


void lstEmptyDirect(nListHandle* h);

void lstPushSingle(void** Head, nListSingle* Item);
void* lstPopSingle(void** Head, nListSingle* Item);

void lstClearPrevNext(nListItem* li);
inline void lstAppendItem(nListHandle* Handle, void* Item) {

	nListItem* li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pFirst)
		Handle->pFirst = Item;

	if (Handle->pLast)
		((nListItem*)Handle->pLast)->pNext = li;

	li->pPrev = Handle->pLast;
	li->pNext = 0;
	Handle->pLast = li;

};
inline void lstPushItem(nListHandle* Handle, void* Item) {

	nListItem* li = Item;

	li->pNext = li->pPrev = 0;

	if (!Handle->pLast)
		Handle->pLast = Item;

	li->pNext = Handle->pFirst;

	if (Handle->pFirst)
		((nListItem*)Handle->pFirst)->pPrev = Item;

	Handle->pFirst = li;

};
inline void* lstPopItem(nListHandle* Handle) {
	void* popitem;
	nListItem* next;
	if (!Handle->pFirst)
		return 0;

	popitem = Handle->pFirst;

	next = ((nListItem*)Handle->pFirst)->pNext;
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

void  lstAppendItem2(nListHandle* Handle, void* Item);
void* lstPopItem2(nListHandle* Handle);
void  lstPushItem2(nListHandle* Handle, void* Item);
void  lstAppendItem3(nListHandle* Handle, void* Item);
void* lstPopItem3(nListHandle* Handle);
void  lstPushItem3(nListHandle* Handle, void* Item);

inline int lstRemoveItem(nListHandle* Handle, nListItem* li) {
	if (!li->pPrev && Handle->pFirst != li) return;

	if (!li->pPrev)
		Handle->pFirst = li->pNext;
	else
		((nListItem*)li->pPrev)->pNext = li->pNext;

	if (!li->pNext)
		Handle->pLast = li->pPrev;
	else
		((nListItem*)li->pNext)->pPrev = li->pPrev;

	li->pNext = li->pPrev = 0;
};
int lstRemoveItem2(nListHandle* Handle, nListItem2* li);

int lstRemoveSegment(nListHandle* Handle, nListItem* Begin, nListItem* End);
void lstInsertItemBefore(nListHandle* Handle, nListItem* toIns, nListItem* pivot);
void lstInsertItemAfter(nListHandle* Handle, nListItem* toIns, nListItem* pivot);
void lstInsertSegmentBefore(nListHandle* Handle, nListItem* Begin, nListItem* End, nListItem* pivot);
void lstInsertSegmentAfter(nListHandle* Handle, nListItem* Begin, nListItem* End, nListItem* pivot);
int   lstHaveItemInList(nListHandle* Handle);
/**/ void* lstGetTop(nListHandle* Handle);

void lstPushSimpleItem(void** first, nItemUserLinker* iul);
void* lstPushItemUser(void** first, void* p);
void* lstPushItemUsing(void** first, void* p);

void* lstAppendPointerOnly(nListHandle* h, void* p);
void* lstAppendPointerSizedOnly(nListHandle* h, void* p, int size);
void* lstPushPointerOnly(nListHandle* h, void* p);
void* lstPushPointerSizedOnly(nListHandle* h, void* p, int size);

void* lstAppendPointer(nListHandle* h, void* p);
void* lstAppendPointerSized(nListHandle* h, void* p, int size);
void* lstPushPointer(nListHandle* h, void* p);
void* lstPushPointerSized(nListHandle* h, void* p, int size);

void* lstAppendPointerStatic(nListHandle* h, nStaticMemoryPool* smp, void* p);
void* lstAppendPointerStaticSized(nListHandle* h, nStaticMemoryPool* smp, void* p, int size);
void* lstPushPointerStatic(nListHandle* h, nStaticMemoryPool* smp, void* p);
void* lstPushPointerStaticSized(nListHandle* h, nStaticMemoryPool* smp, void* p, int size);

void* lstPopPointerOnly(nListHandle* h);
void lstRemovePointerItemOnly(nListHandle* h, nListItemPointer* lip);
void lstRemovePointerOnly(nListHandle* h, void* p);
void lstClearPointerOnly(nListHandle* h);
void lstGeneratePointerListOnly(nListHandle* from1, nListHandle* from2, nListHandle* to);

void* lstPopPointer(nListHandle* h);
void lstRemovePointerItem(nListHandle* h, nListItemPointer* lip);
void lstRemovePointer(nListHandle* h, void* p);
void lstClearPointer(nListHandle* h);
void lstGeneratePointerList(nListHandle* from1, nListHandle* from2, nListHandle* to);

void lstCopyHandle(nListHandle* target, nListHandle* src);

void* lstAppendPointerStaticPool(nStaticMemoryPool* mph, nListHandle* h, void* p);
void* lstPopPointerNoFree(nListHandle* h);
void lstRemovePointerItemNoFree(nListHandle* h, nListItemPointer* lip);

void lstMoveUp(nListHandle* h, nListItem* li);
void lstMoveDown(nListHandle* h, nListItem* li);

void  lstForAllItemsDo(nListDoFunc func, nListHandle* hList);
void lstForAllItemsDoLNRR(nListNonRecursiveDoFunc func, nListHandle* hList);
void lstForAllItemsDo_DirectFree(nListDoFunc func, nListHandle* hList);
void lstForAllItemsDo_arg_ptr(nListDoFuncArgp func, nListHandle* hList, void* arg);
void lstForAllItemsDo_NonRecursive_Root(nListHandle* FirstHandle, nListNonRecursiveDoFunc func, int bFreeItem, void* custom_data, nListCustomDataRemover remover);
void lstCopy_NonRecursive_Root(nListHandle* FromHandle, nListHandle* ToHandle, int SizeEachNode, nListNonRecursiveCopyFunc func, void* custom_data, nListCustomDataRemover remover);
void lstAddNonRecursiveListHandle(nListNonRecursiveRoot* root, nListHandle* newHandle, nListNonRecursiveDoFunc nrFunc, int bFreeList, void* custom_data, nListCustomDataRemover remover);
void lstAddNonRecursiveListCopier(nListNonRecursiveRoot* root, nListHandle* oldHandle, nListHandle* newHandle, int sizeEach, nListNonRecursiveCopyFunc nrCpyFunc, void* custom_data, nListCustomDataRemover remover);
void* lstFindItem(void* CmpData, nCompareFunc func, nListHandle* hList);
void lstCombineLists(nListHandle* dest, nListHandle* src);
void lstDestroyList(nListHandle* hlst);
void* lstReMatch(nListHandle* SearchHandle, nListHandle* CurrentHandle, void* ItemToFind);
typedef int(*MatcherFunc)(void*, void*);
void* lstReMatchEx(nListHandle* SearchHandle, nListHandle* CurrentHandle, void* ItemToFind, MatcherFunc func);

void lstAddElement(nListHandle* hlst, void* ext);
void lstDestroyElementList(nListHandle* hlst);

void hsh256InsertItemCSTR(nHash256* hash, nListItem* li, char* buckle);
void hsh256InsertItem(nHash256* hash, nListItem* li, unsigned char buckle);
void hsh65536InsertItem(nHash65536* hash, nListItem* li, long buckle);
__inline nListHandle* hsh65536DoHashLongPtr(nHash65536* hash, long buckle) {
	return &hash->Entries[(unsigned short)((buckle >> 10))];
}
__inline nListHandle* hsh65536DoHashNUID(nHash65536* hash, char* NUID) {
	u64bit Hash;
	sscanf(NUID, "%I64d", &Hash);
	return hsh65536DoHashLongPtr(hash, (long)Hash);
}
__inline nListHandle* hsh16MDoHashLongPtr(nHash16M* hash, long long buckle) {
	return &hash->Entries[(buckle>>6)&0x00FFFFFF];
}
__inline nListHandle* hsh16MDoHashNUID(nHash16M* hash, char* NUID) {
	u64bit Hash;
	sscanf(NUID, "%I64d", &Hash);
	return hsh65536DoHashLongPtr(hash, (long)Hash);
}

nListItem* hsh256FindItemCSTR(nHash256* hash, nCompareFunc func, char* buckle);
unsigned char hsh256DoHashCSTR(char* buckle);


void memResetByteCount();
int memGetByteCount();
void memInitPool(nMemoryPool* mph, int NodeSize);
void memInitPoolSmall(nMemoryPool* mph, int NodeSize);
nMemoryPoolPart* memNewPoolPart(nMemoryPool* mph);
void* memAquireH(nMemoryPool* Handle);
void* memAquireOnly(int Size);
void* memAquire(int Size);
void* memAquireHyper(int Size);
void* memAquireHyper1(int Size);
void memAssignDBInst(void* mem, nDBInst* DBInst);
nDBInst* memGetDBInst(void* mem);
void memFree(void* Data);
void memDestroyPool(nMemoryPool* Handle);

nStaticMemoryPoolNode* memNewStaticPool(nStaticMemoryPool* smp);
void* memStaticAquire(nStaticMemoryPool*smp, int size);
void* memStaticAquireThread(nStaticMemoryPool*smp, int size);
void* memStaticDestroy(nStaticMemoryPool*smp);

nAVLNodeReal64* nAVLTreeLookupJustSmaller(nAVLNodeReal64* Root, real Value);
nAVLNodeReal64* nAVLTreeLookupJustGreater(nAVLNodeReal64* Root, real Value);
nAVLNodeReal64* nAVLTreeInsert(nAVLTreeReal64* Tree, real Value);
void nAVLTreePrintPreOrder(nAVLTreeReal64* Tree);
void nAVLTreeMakeListUpOrder(nAVLTreeReal64* Tree, nListHandle* Result);
void nAVLTreeDestroy(nAVLNodeReal64* Root);

int  strGetStringTerminateBy(char* content, char terminator, char* Out);

int strHeadOfStringMatch(char* Str, char* SubStr);
int strSkipSegmet(char** pivot, char* content);
char* strgetLastSegmentSeperateBy(char* Content, char Seperator);
void strDiscardLastSegmentSeperateBy(char* Content, char Seperator);
void strDiscardSameBeginningSeperatedBy(char* s1, char* s2, char** Result1, char** Result2, char Seperator);
int strCountSegmentSeperateBy(char* Content, char Seperator);
void strMakeDifferentName(char* Target);

void strReplaceCharacter(char* Str, char Find, char Replace);
void strToUpperCase(char* Str);
void strToLowerCase(char* Str);

nStringSplitor* strSplitPath(char* path);
int strMakeInstructions(nStringSplitor** result,char* content);
nStringPart* strGetArgument(nStringSplitor* ss, char* content);
char* strGetArgumentString(nStringSplitor* ss, char* content);
int strArgumentMatch(nStringSplitor* ss, char* id, char* value);
int strDestroyStringSplitor(nStringSplitor** ss);

int strGetIntSimple(char* content);
real strGetFloatSimple(char* content);

void strConvInt_CString(int src, char* dest, int lenth);
void strConvFloat_CString(real src, char* dest, int lenth);

void strCopyFull(char* dest, char* src);
void strCopySized(char* dest, int LenthLim, char* src);
#define strAppend strcat_s
void strPrintFloatAfter(char* dest, int LenthLim, int bits, real data);
void strPrintIntAfter(char* dest, int LenthLim, int data);
void strToWideChar(wchar_t* destBuf, char* srcBuf);

int strIsTheSame(char* src, char*dest);


void strSafeDestroy(nSafeString** ss);
void strSafeSet(nSafeString** ss, char* Content);




nStringEdit* strBeginEdit(char* FullStr);
void strEndEdit(nStringEdit* se, char* Result);
void strRemoveLine(nStringEdit* se, nStringLine* sl);
void strRemoveLineI(nStringEdit* se, int LineIndex);
void strSetCusor(nStringEdit* se, int LineIndex, int BeforeIndex);
void strMoveCusor(nStringEdit* se, int Left);
void strSelect(nStringEdit* se, int BeginLine, int BeginBefore, int EndLine, int EndBefore);
void strSelectLineAll(nStringEdit* se);
void strDeselectAll(nStringEdit* se);
void strPanFoward(char* str, int Before, int Offset);
void strSquishBackward(char* str, int Before, int EndBefore);
void strClearSelection(nStringEdit* se);
nStringLine* strGetCursorLine(nStringEdit* se);
nStringLine* strGetBeginLine(nStringEdit* se);
void strInsertChar(nStringEdit* se, char a);
void strBackspace(nStringEdit* se);

void transNewLanguage(const char* LanguageID);
void transSetLanguage(const char* LanguageID);
void transDumpMissMatchRecord(const char* filename);
void transNewEntry(const char* Target, const char* replacement);
char* transLate(char* Target);
void transState(void* UNUSED, int val);
void transInitTranslation_zh_cn();


void nulOpenInternetLink(char* link);


