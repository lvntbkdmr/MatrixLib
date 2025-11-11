#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#ifdef _WIN32
#include <windows.h>
#else
// For non-Windows platforms, define INT64 if needed
#ifndef __int64
typedef long long __int64;
#endif
#endif

typedef int INT32;
typedef char CHAR;
typedef char* CHARP;
typedef const char* CCHARP;
typedef double REAL64;
typedef short INT16;
typedef unsigned short UINT16;
typedef unsigned int UINT32;
typedef unsigned char UINT8;
typedef unsigned char* UINT8P;
typedef float REAL32;
typedef float* REAL32P;
typedef unsigned short* UINT16P;
typedef unsigned int* UINT32P;
typedef unsigned char* UINT8P;
typedef short* INT16P;
typedef int* INT32P;
typedef double* REAL64P;
typedef signed char INT8;
typedef signed char* INT8P;
typedef void* VOIDP;
typedef __int64 INT64;

#endif