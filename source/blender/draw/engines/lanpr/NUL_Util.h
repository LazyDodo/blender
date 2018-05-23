#pragma once

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/

#define _CRT_SECURE_NO_WARNINGS


typedef double real;
typedef unsigned long long u64bit;
typedef unsigned int       u32bit;
typedef unsigned short     u16bit;
typedef unsigned short     ushort;
typedef unsigned char      u8bit;
typedef char nShortBuf[16];

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


char* txtReadFileAsString(char* FileName);


void nulSendPanic(char* message);

#define SEND_PANIC_ERROR(message)\
    nulSendPanic(message)

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
