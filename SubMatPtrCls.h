/**
 * @file SubMatPtrCls.h
 * @brief Submatrix pointer class for non-owning matrix views
 * 
 * This header defines a class that provides a non-owning view (pointer-based)
 * into a portion of a larger matrix. This allows efficient access to submatrices
 * without copying data.
 * 
 * Key features:
 * - Non-owning view: points to elements in the original matrix
 * - Compile-time maximum dimensions (RowNum, ColNum)
 * - Runtime variable dimensions (SubRowNum, SubColNum)
 * - Vector-like access for column vectors (ColNum == 1)
 * 
 * @tparam RowNum Maximum number of rows (compile-time constant)
 * @tparam ColNum Maximum number of columns (compile-time constant)
 * @tparam ElementType Type of matrix elements (e.g., REAL64, REAL32)
 * 
 * @note This class is designed for avionics software where dynamic memory
 * allocation is prohibited. All pointers point to stack-allocated data.
 * 
 * @example
 * MatrixCls<6, 6, REAL64> mat;
 * mat.setSubMatrixSize(4, 4);
 * // Get a 3x3 submatrix view starting at (1, 1)
 * SubMatPtrCls<6, 6, REAL64> sub = mat(1, 1, 3, 3);
 * sub(0, 0) = 5.0;  // Modifies mat(1, 1)
 */

#ifndef SUB_MAT_PTR_CLS_H
#define SUB_MAT_PTR_CLS_H

#include "CommonTypes.h"
#include <cstring>

/**
 * @class SubMatPtrCls
 * @brief Non-owning view into a submatrix
 * 
 * This class stores pointers to elements in a parent matrix, allowing
 * efficient access to submatrices without copying data. Modifications
 * through this view directly affect the original matrix.
 * 
 * The class supports both matrix access (row, col) and vector-like access
 * (index) for column vectors (when ColNum == 1).
 */
template <UINT32 RowNum, UINT32 ColNum, class ElementType>
class SubMatPtrCls
{
public:
    /**
     * @brief Default constructor - creates an empty view
     * 
     * Initializes all pointers to nullptr and sets dimensions to 0.
     */
    SubMatPtrCls();
    
    /**
     * @brief Destructor
     * 
     * No cleanup needed since this class doesn't own the data.
     */
    ~SubMatPtrCls();

    ElementType* ElementsPtr[RowNum][ColNum]; ///< Array of pointers to matrix elements
    UINT32 SubRowNum;                         ///< Runtime number of rows (<= RowNum)
    UINT32 SubColNum;                         ///< Runtime number of columns (<= ColNum)

    /**
     * @brief Access matrix element by row and column (non-const)
     * @param row Row index (0-based)
     * @param col Column index (0-based)
     * @return Reference to the element at (row, col)
     * 
     * @note No bounds checking is performed. Ensure row < SubRowNum and col < SubColNum.
     */
    ElementType& operator()(UINT32 row, UINT32 col);
    
    /**
     * @brief Access matrix element by row and column (const)
     * @param row Row index (0-based)
     * @param col Column index (0-based)
     * @return Value of the element at (row, col)
     */
    ElementType operator()(UINT32 row, UINT32 col) const;
    
    /**
     * @brief Vector-like access for column vectors (non-const)
     * @param index Element index (0-based, row index when ColNum == 1)
     * @return Reference to the element at the specified index
     * 
     * This operator is provided for convenience when working with column vectors.
     * It is equivalent to operator()(index, 0).
     * 
     * @note Only valid when ColNum == 1 (column vector).
     */
    ElementType& operator()(UINT32 index);
    
    /**
     * @brief Vector-like access for column vectors (const)
     * @param index Element index (0-based, row index when ColNum == 1)
     * @return Value of the element at the specified index
     */
    ElementType operator()(UINT32 index) const;
    
    /**
     * @brief Get the number of rows in the submatrix
     * @return SubRowNum (runtime number of rows)
     */
    UINT32 getRowCount() const { return SubRowNum; }
    
    /**
     * @brief Get the number of columns in the submatrix
     * @return SubColNum (runtime number of columns)
     */
    UINT32 getColCount() const { return SubColNum; }
    
    /**
     * @brief Get the length (for vector compatibility)
     * @return SubRowNum (number of rows, useful for column vectors)
     * 
     * This method is provided for compatibility with vector interfaces.
     * For column vectors (ColNum == 1), this returns the number of elements.
     */
    UINT32 getLength() const { return SubRowNum; }
};

// ============================================================================
// Constructor Implementation
// ============================================================================

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
SubMatPtrCls<RowNum, ColNum, ElementType>::SubMatPtrCls() : SubRowNum(0), SubColNum(0)
{
    // Initialize all pointers to nullptr
    // This ensures the view is empty and safe to use
    memset(&ElementsPtr, 0, sizeof(ElementType *) * RowNum * ColNum);
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
SubMatPtrCls<RowNum, ColNum, ElementType>::~SubMatPtrCls()
{
    // No cleanup needed - this class doesn't own the data
    // The pointers point to elements in the parent matrix
}

// ============================================================================
// Element Access Implementations
// ============================================================================

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType& SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 row, UINT32 col)
{
    // Dereference the pointer to get the actual element
    // This directly accesses the element in the parent matrix
    return *(ElementsPtr[row][col]);
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 row, UINT32 col) const
{
    // Const version: returns a copy of the element value
    return *(ElementsPtr[row][col]);
}

// ============================================================================
// Vector-like Access Implementations (for column vectors)
// ============================================================================

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType& SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 index)
{
    // Vector-like access: treat as column vector (column 0)
    // This is convenient when working with column vectors (ColNum == 1)
    return (*this)(index, 0);
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 index) const
{
    // Const version of vector-like access
    return (*this)(index, 0);
}


#endif