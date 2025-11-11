/**
 * @file CommonTypes.h
 * @brief Common type definitions for the Matrix Library
 * 
 * This header provides platform-independent type definitions used throughout
 * the matrix library. It ensures consistent data types across different
 * platforms (Windows, Linux, macOS) and provides type aliases for better
 * code readability and portability.
 * 
 * @note This library is designed for avionics software where dynamic memory
 * allocation is prohibited. All types are designed to work with stack-allocated
 * data structures.
 */

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

// Integer types
typedef int INT32;              ///< 32-bit signed integer
typedef short INT16;            ///< 16-bit signed integer
typedef signed char INT8;       ///< 8-bit signed integer
typedef __int64 INT64;          ///< 64-bit signed integer

// Unsigned integer types
typedef unsigned int UINT32;    ///< 32-bit unsigned integer
typedef unsigned short UINT16;  ///< 16-bit unsigned integer
typedef unsigned char UINT8;    ///< 8-bit unsigned integer

// Floating point types
typedef double REAL64;          ///< 64-bit floating point (double precision)
typedef float REAL32;           ///< 32-bit floating point (single precision)

// Character types
typedef char CHAR;              ///< Character type
typedef char* CHARP;            ///< Pointer to character (C-style string)
typedef const char* CCHARP;     ///< Pointer to constant character (const C-style string)

// Pointer types
typedef void* VOIDP;            ///< Generic void pointer
typedef float* REAL32P;         ///< Pointer to 32-bit float
typedef double* REAL64P;        ///< Pointer to 64-bit float
typedef unsigned short* UINT16P; ///< Pointer to 16-bit unsigned integer
typedef unsigned int* UINT32P;  ///< Pointer to 32-bit unsigned integer
typedef unsigned char* UINT8P;  ///< Pointer to 8-bit unsigned integer
typedef short* INT16P;          ///< Pointer to 16-bit signed integer
typedef int* INT32P;            ///< Pointer to 32-bit signed integer
typedef signed char* INT8P;     ///< Pointer to 8-bit signed integer

#endif