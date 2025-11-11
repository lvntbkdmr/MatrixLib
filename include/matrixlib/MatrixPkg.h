/**
 * @file MatrixPkg.h
 * @brief Package header that includes all matrix library components
 * 
 * This header serves as the main entry point for the matrix library.
 * It includes all necessary headers in the correct order to resolve
 * dependencies between classes.
 * 
 * @note Include this header to use the complete matrix library functionality.
 * 
 * @example
 * #include "MatrixPkg.h"
 * 
 * MatrixCls<6, 6, REAL64> A;
 * VectorCls<6, REAL64> b, x;
 * A.LeftDivide(b, x);
 */

#ifndef MATRIX_PKG_H
#define MATRIX_PKG_H

#include "matrixlib/CommonTypes.h"

// Forward declarations to resolve circular dependencies
template <UINT32 RowNum, UINT32 ColNum, class ElementType> class MatrixCls;
template <UINT32 RowNum, UINT32 ColNum, class ElementType> class SubMatPtrCls;

// Include class definitions in dependency order
#include "matrixlib/MatrixCls.h"
#include "matrixlib/SubMatPtrCls.h"
#include "matrixlib/VectorCls.h"

#endif