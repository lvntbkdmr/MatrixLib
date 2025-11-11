#ifndef MATRIX_PKG_H
#define MATRIX_PKG_H

#include "CommonTypes.h"

// Forward declarations
template <UINT32 RowNum, UINT32 ColNum, class ElementType> class MatrixCls;
template <UINT32 RowNum, UINT32 ColNum, class ElementType> class SubMatPtrCls;

// Include class definitions
#include "MatrixCls.h"
#include "SubMatPtrCls.h"
#include "VectorCls.h"

#endif