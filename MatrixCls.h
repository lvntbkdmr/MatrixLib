#ifndef MATRIX_CLS_H
#define MATRIX_CLS_H

#include "CommonTypes.h"
#include "SubMatPtrCls.h"
#include "VectorCls.h"
#include <cstring>
#include <cstdlib>
#include <type_traits>
#include <cmath>



template <UINT32 RowNum = 3, UINT32 ColNum = RowNum, class ElementType = REAL64>
class MatrixCls
{
public:
    MatrixCls();
    MatrixCls(const ElementType& value);
    MatrixCls(VOIDP basicMatrixVoidAdress);
    MatrixCls(CHARP option);
    MatrixCls(const ElementType& val0,
              const ElementType& val1);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3,
              const ElementType& val4);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3,
              const ElementType& val4,
              const ElementType& val5);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3,
              const ElementType& val4,
              const ElementType& val5,
              const ElementType& val6);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3,
              const ElementType& val4,
              const ElementType& val5,
              const ElementType& val6,
              const ElementType& val7);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3,
              const ElementType& val4,
              const ElementType& val5,
              const ElementType& val6,
              const ElementType& val7,
              const ElementType& val8);
    MatrixCls(const ElementType& val0,
              const ElementType& val1,
              const ElementType& val2,
              const ElementType& val3,
              const ElementType& val4,
              const ElementType& val5,
              const ElementType& val6,
              const ElementType& val7,
              const ElementType& val8,
              const ElementType& val9);
    //operator overloads for =, +, -, /, *, !, ==, !=, ~, and ()

    SubMatPtrCls<RowNum, ColNum, ElementType> operator()(UINT32 rowStart, 
                                                         UINT32 colStart, 
                                                         UINT32 numRows, 
                                                         UINT32 numCols);

    MatrixCls<RowNum, ColNum, ElementType> operator()(UINT32 rowStart, 
                                                       UINT32 colStart, 
                                                       UINT32 numRows, 
                                                       UINT32 numCols) const;
    
    inline ElementType& operator()(UINT32 row, UINT32 col);

    ElementType operator()(UINT32 row, UINT32 col) const;

    ElementType& operator()(UINT32 index);

    ElementType* operator[](UINT32 rowIndex);

    UINT32 getRowCount() const { return SubRowNum; }
    UINT32 getColCount() const { return SubColNum; }

    void setSubMatrixSize(UINT32 numRows, UINT32 numCols);

    void LUDecompose(MatrixCls<RowNum, ColNum, ElementType>& L, 
                     MatrixCls<RowNum, ColNum, ElementType>& U) const;
    
    // Instance methods for matrix algorithms
    void ForwardSubstitution(const VectorCls<RowNum, ElementType>& b, VectorCls<RowNum, ElementType>& y);
    
    void BackwardSubstitution(const VectorCls<RowNum, ElementType>& y, VectorCls<RowNum, ElementType>& x);
    
    void LUDecomp(MatrixCls<RowNum, ColNum, ElementType>& L, MatrixCls<RowNum, ColNum, ElementType>& U);
    
    // Matrix property checks
    bool IsLowTri() const;
    bool IsUpTri() const;
    bool IsSymmetric() const;
    
    // Matrix operations
    void Transpose(MatrixCls<ColNum, RowNum, ElementType>& result) const;
    
    // Decomposition methods
    void CholeskyDecomp(MatrixCls<RowNum, ColNum, ElementType>& L);
    void QrDecomp(MatrixCls<RowNum, RowNum, ElementType>& Q, MatrixCls<RowNum, ColNum, ElementType>& R);
    
    bool LeftDivide(const VectorCls<RowNum, ElementType>& b, VectorCls<ColNum, ElementType>& x);
    
    bool LeftDivide(const MatrixCls<RowNum, ColNum, ElementType>& b, MatrixCls<RowNum, ColNum, ElementType>& x);
    
    ~MatrixCls();

    ElementType Elements[RowNum][ColNum];
    UINT32 SubRowNum;
    UINT32 SubColNum;

    void Zeros();
    void Ones();
    void Random();
    void Identity();
    void Eye();
};

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls() : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& value) : SubRowNum(RowNum), SubColNum(ColNum)
{
    for (UINT32 i = 0; i < RowNum; i++)
    {
        for (UINT32 j = 0; j < ColNum; j++)
        {
            Elements[i][j] = value;
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(VOIDP basicMatrixVoidAdress) : SubRowNum(RowNum), SubColNum(ColNum)
{
    // Copy from void pointer (assuming it points to RowNum*ColNum elements)
    ElementType* src = static_cast<ElementType*>(basicMatrixVoidAdress);
    for (UINT32 i = 0; i < RowNum; i++)
    {
        for (UINT32 j = 0; j < ColNum; j++)
        {
            Elements[i][j] = src[i * ColNum + j];
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(CHARP option) : SubRowNum(RowNum), SubColNum(ColNum)
{
    if (option && strcmp(option, "zeros") == 0)
    {
        Zeros();
    }
    else if (option && strcmp(option, "ones") == 0)
    {
        Ones();
    }
    else if (option && strcmp(option, "eye") == 0)
    {
        Eye();
    }
    else if (option && strcmp(option, "identity") == 0)
    {
        Identity();
    }
    else
    {
        Zeros();
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3, const ElementType& val4) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
    if (RowNum * ColNum >= 5) Elements[0][4] = val4;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3, const ElementType& val4, const ElementType& val5) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
    if (RowNum * ColNum >= 5) Elements[0][4] = val4;
    if (RowNum * ColNum >= 6) Elements[0][5] = val5;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3, const ElementType& val4, const ElementType& val5, const ElementType& val6) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
    if (RowNum * ColNum >= 5) Elements[0][4] = val4;
    if (RowNum * ColNum >= 6) Elements[0][5] = val5;
    if (RowNum * ColNum >= 7) Elements[0][6] = val6;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3, const ElementType& val4, const ElementType& val5, const ElementType& val6, const ElementType& val7) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
    if (RowNum * ColNum >= 5) Elements[0][4] = val4;
    if (RowNum * ColNum >= 6) Elements[0][5] = val5;
    if (RowNum * ColNum >= 7) Elements[0][6] = val6;
    if (RowNum * ColNum >= 8) Elements[0][7] = val7;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3, const ElementType& val4, const ElementType& val5, const ElementType& val6, const ElementType& val7, const ElementType& val8) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
    if (RowNum * ColNum >= 5) Elements[0][4] = val4;
    if (RowNum * ColNum >= 6) Elements[0][5] = val5;
    if (RowNum * ColNum >= 7) Elements[0][6] = val6;
    if (RowNum * ColNum >= 8) Elements[0][7] = val7;
    if (RowNum * ColNum >= 9) Elements[0][8] = val8;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& val0, const ElementType& val1, const ElementType& val2, const ElementType& val3, const ElementType& val4, const ElementType& val5, const ElementType& val6, const ElementType& val7, const ElementType& val8, const ElementType& val9) 
    : SubRowNum(RowNum), SubColNum(ColNum)
{
    Zeros();
    if (RowNum * ColNum >= 1) Elements[0][0] = val0;
    if (RowNum * ColNum >= 2) Elements[0][1] = val1;
    if (RowNum * ColNum >= 3) Elements[0][2] = val2;
    if (RowNum * ColNum >= 4) Elements[0][3] = val3;
    if (RowNum * ColNum >= 5) Elements[0][4] = val4;
    if (RowNum * ColNum >= 6) Elements[0][5] = val5;
    if (RowNum * ColNum >= 7) Elements[0][6] = val6;
    if (RowNum * ColNum >= 8) Elements[0][7] = val7;
    if (RowNum * ColNum >= 9) Elements[0][8] = val8;
    if (RowNum * ColNum >= 10) Elements[0][9] = val9;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
SubMatPtrCls<RowNum, ColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator()(UINT32 rowStart, 
    UINT32 colStart, 
    UINT32 numRows, 
    UINT32 numCols)
{
    SubMatPtrCls<RowNum, ColNum, ElementType> subMatPtr;

    for(UINT32 i = 0; i < numRows; i++)
    {
        for(UINT32 j = 0; j < numCols; j++)
        {
            subMatPtr.ElementsPtr[i][j] = &(Elements[rowStart + i][colStart + j]);
        }
    }

    subMatPtr.SubRowNum = numRows;
    subMatPtr.SubColNum = numCols;

    return subMatPtr;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator()(UINT32 rowStart, 
    UINT32 colStart, 
    UINT32 numRows, 
    UINT32 numCols) const
{
    MatrixCls<RowNum, ColNum, ElementType> returnMat;

    for(UINT32 i = 0; i < numRows; i++)
    {
        for(UINT32 j = 0; j < numCols; j++)
        {
            returnMat.Elements[i][j] = Elements[rowStart + i][colStart + j];
        }
    }

    returnMat.SubRowNum = numRows;
    returnMat.SubColNum = numCols;

    return returnMat;
}

// Implementation of LUDecompose using LUDecomp static method
template <UINT32 RowNum, UINT32 ColNum, class ElementType>
inline ElementType& MatrixCls<RowNum, ColNum, ElementType>::operator()(UINT32 row, UINT32 col)
{
    return Elements[row][col];
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType MatrixCls<RowNum, ColNum, ElementType>::operator()(UINT32 row, UINT32 col) const
{
    return Elements[row][col];
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType& MatrixCls<RowNum, ColNum, ElementType>::operator()(UINT32 index)
{
    // Treat matrix as row-major 1D array
    return Elements[index / ColNum][index % ColNum];
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType* MatrixCls<RowNum, ColNum, ElementType>::operator[](UINT32 rowIndex)
{
    return Elements[rowIndex];
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::~MatrixCls()
{
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::Zeros()
{
    for (UINT32 i = 0; i < RowNum; i++)
    {
        for (UINT32 j = 0; j < ColNum; j++)
        {
            Elements[i][j] = ElementType(0);
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::Ones()
{
    for (UINT32 i = 0; i < RowNum; i++)
    {
        for (UINT32 j = 0; j < ColNum; j++)
        {
            Elements[i][j] = ElementType(1);
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::Random()
{
    // Simple random number generation (using modulo for portability)
    // For avionics, you might want to use a deterministic PRNG
    for (UINT32 i = 0; i < RowNum; i++)
    {
        for (UINT32 j = 0; j < ColNum; j++)
        {
            // Placeholder - implement with your preferred RNG
            Elements[i][j] = ElementType(rand() % 100);
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::Identity()
{
    Zeros();
    UINT32 minDim = (RowNum < ColNum) ? RowNum : ColNum;
    for (UINT32 i = 0; i < minDim; i++)
    {
        Elements[i][i] = ElementType(1);
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::Eye()
{
    Identity();
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::setSubMatrixSize(UINT32 numRows, UINT32 numCols)
{
    SubRowNum = numRows;
    SubColNum = numCols;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::LUDecompose(
    MatrixCls<RowNum, ColNum, ElementType>& L, 
    MatrixCls<RowNum, ColNum, ElementType>& U) const
{
    // Create a mutable copy for LU decomposition
    MatrixCls<RowNum, ColNum, ElementType> A_copy = *this;
    A_copy.LUDecomp(L, U);
}

// Instance method implementations
template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::ForwardSubstitution(const VectorCls<RowNum, ElementType>& b, VectorCls<RowNum, ElementType>& y)
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    
    // Check if matrix is square
    if (n != m)
    {
        return; // Matrix must be square for forward substitution
    }
    
    // Bounds check: ensure runtime size doesn't exceed vector compile-time sizes
    if (n > b.getLength() || n > y.getLength())
    {
        return; // Dimension mismatch - runtime size exceeds vector capacity
    }
    
    for (UINT32 i = 0; i < n; i++)
    {
        y(i) = b(i);
        for (UINT32 j = 0; j < i; j++)
        {
            y(i) = y(i) - (*this)(i, j) * y(j);
        }

        y(i) = y(i) / (*this)(i, i);
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::BackwardSubstitution(const VectorCls<RowNum, ElementType>& y, VectorCls<RowNum, ElementType>& x)
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    
    // Check if matrix is square
    if (n != m)
    {
        return; // Matrix must be square for backward substitution
    }
    
    // Bounds check: ensure runtime size doesn't exceed vector compile-time sizes
    if (n > y.getLength() || n > x.getLength())
    {
        return; // Dimension mismatch - runtime size exceeds vector capacity
    }
    
    for (INT32 i = n - 1; i >= 0; i--)
    {
        x(i) = y(i);
        for (UINT32 j = i + 1; j < n; j++)
        {
            x(i) = x(i) - (*this)(i, j) * x(j);
        }
        x(i) = x(i) / (*this)(i, i);
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::LUDecomp(MatrixCls<RowNum, ColNum, ElementType>& L, MatrixCls<RowNum, ColNum, ElementType>& U)
{
    UINT32 n = this->getRowCount();
    
    if (n != this->getColCount() || n != L.getRowCount() || n != L.getColCount() || 
        n != U.getRowCount() || n != U.getColCount())
    {
        return; // Dimension mismatch
    }

    // Initialize L as zero matrix (lower triangle will be filled)
    // Initialize U as identity matrix (diagonal = 1, upper triangle will be filled)
    for (UINT32 i = 0; i < n; i++)
    {
        for (UINT32 j = 0; j < n; j++)
        {
            L(i, j) = ElementType(0); // Initialize L
            if (i == j)
            {
                U(i, j) = ElementType(1); // Identity diagonal for U
            }
            else
            {
                U(i, j) = ElementType(0); // Will be computed for upper triangle
            }
        }
    }

    // Perform LU decomposition using Crout's method
    // In Crout's method: L has actual values on diagonal, U has 1s on diagonal
    for (UINT32 k = 0; k < n; k++)
    {
        // Compute L(i, k) for i >= k
        for (UINT32 i = k; i < n; i++)
        {
            ElementType sum = (*this)(i, k);
            for (UINT32 p = 0; p < k; p++)
            {
                sum = sum - L(i, p) * U(p, k);
            }
            L(i, k) = sum;
        }

        // Check for singular matrix (zero diagonal in L)
        if (L(k, k) == ElementType(0))
        {
            return; // Singular matrix
        }

        // Compute U(k, j) for j > k
        for (UINT32 j = k + 1; j < n; j++)
        {
            ElementType sum = (*this)(k, j);
            for (UINT32 p = 0; p < k; p++)
            {
                sum = sum - L(k, p) * U(p, j);
            }
            U(k, j) = sum / L(k, k);
        }
    }
}

// Matrix property check implementations
template <UINT32 RowNum, UINT32 ColNum, class ElementType>
bool MatrixCls<RowNum, ColNum, ElementType>::IsLowTri() const
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    
    if (n != m) return false;
    
    const ElementType tolerance = ElementType(1e-10);
    for (UINT32 i = 0; i < n; i++)
    {
        for (UINT32 j = i + 1; j < n; j++)
        {
            ElementType val = (*this)(i, j);
            if (val < ElementType(0)) val = -val;
            if (val > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
bool MatrixCls<RowNum, ColNum, ElementType>::IsUpTri() const
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    
    if (n != m) return false;
    
    const ElementType tolerance = ElementType(1e-10);
    for (UINT32 i = 0; i < n; i++)
    {
        for (UINT32 j = 0; j < i; j++)
        {
            ElementType val = (*this)(i, j);
            if (val < ElementType(0)) val = -val;
            if (val > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
bool MatrixCls<RowNum, ColNum, ElementType>::IsSymmetric() const
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    
    if (n != m) return false;
    
    const ElementType tolerance = ElementType(1e-10);
    for (UINT32 i = 0; i < n; i++)
    {
        for (UINT32 j = i + 1; j < n; j++)
        {
            ElementType diff = (*this)(i, j) - (*this)(j, i);
            if (diff < ElementType(0)) diff = -diff;
            if (diff > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::Transpose(MatrixCls<ColNum, RowNum, ElementType>& result) const
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    
    result.setSubMatrixSize(m, n);
    
    for (UINT32 i = 0; i < n; i++)
    {
        for (UINT32 j = 0; j < m; j++)
        {
            result(j, i) = (*this)(i, j);
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::CholeskyDecomp(MatrixCls<RowNum, ColNum, ElementType>& L)
{
    UINT32 n = this->getRowCount();
    
    if (n != this->getColCount() || n != L.getRowCount() || n != L.getColCount())
    {
        return; // Dimension mismatch
    }
    
    L.setSubMatrixSize(n, n);
    L.Zeros();
    
    // Cholesky decomposition: A = L * L^T
    for (UINT32 i = 0; i < n; i++)
    {
        for (UINT32 j = 0; j <= i; j++)
        {
            ElementType sum = (*this)(i, j);
            for (UINT32 k = 0; k < j; k++)
            {
                sum -= L(i, k) * L(j, k);
            }
            
            if (i == j)
            {
                if (sum <= ElementType(0))
                {
                    return; // Matrix is not positive definite
                }
                L(i, j) = sqrt(sum);
            }
            else
            {
                if (L(j, j) == ElementType(0))
                {
                    return; // Division by zero
                }
                L(i, j) = sum / L(j, j);
            }
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
void MatrixCls<RowNum, ColNum, ElementType>::QrDecomp(MatrixCls<RowNum, RowNum, ElementType>& Q, MatrixCls<RowNum, ColNum, ElementType>& R)
{
    UINT32 m = this->getRowCount();
    UINT32 n = this->getColCount();
    
    Q.setSubMatrixSize(m, m);
    R.setSubMatrixSize(m, n);
    R.Zeros();
    
    // Initialize Q as identity
    Q.Identity();
    Q.setSubMatrixSize(m, m);
    
    // Create working copy
    MatrixCls<RowNum, ColNum, ElementType> A_orig = *this;
    A_orig.setSubMatrixSize(m, n);
    
    // Simplified QR using modified Gram-Schmidt
    UINT32 minDim = (m < n) ? m : n;
    
    for (UINT32 k = 0; k < minDim; k++)
    {
        // Compute R(k,k) and Q column k
        ElementType norm = ElementType(0);
        for (UINT32 i = 0; i < m; i++)
        {
            Q(i, k) = A_orig(i, k);
            norm += Q(i, k) * Q(i, k);
        }
        norm = sqrt(norm);
        
        if (norm == ElementType(0)) continue;
        
        R(k, k) = norm;
        for (UINT32 i = 0; i < m; i++)
        {
            Q(i, k) /= norm;
        }
        
        // Compute R(k,j) for j > k and update A
        for (UINT32 j = k + 1; j < n; j++)
        {
            ElementType dot = ElementType(0);
            for (UINT32 i = 0; i < m; i++)
            {
                dot += Q(i, k) * A_orig(i, j);
            }
            R(k, j) = dot;
            
            for (UINT32 i = 0; i < m; i++)
            {
                A_orig(i, j) -= Q(i, k) * dot;
            }
        }
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
bool MatrixCls<RowNum, ColNum, ElementType>::LeftDivide(const VectorCls<RowNum, ElementType>& b, VectorCls<ColNum, ElementType>& x)
{
    UINT32 n = this->getRowCount();
    UINT32 m = this->getColCount();
    bool isSquare = (n == m);
    
    if (isSquare)
    {
        
        // Lower Triangular test
        bool isLowTri = IsLowTri();
        if (isLowTri)
        {
            ForwardSubstitution(b, x);
            return true;
        }
        else
        {
            // Upper triangular check
            bool isUpTri = IsUpTri();
            if (isUpTri)
            {
                BackwardSubstitution(b, x);
                return true;
            }
            else
            {
                bool isSymmetric = IsSymmetric();
                if (isSymmetric)
                {
                    MatrixCls<RowNum, RowNum, ElementType> lowTri;
                    lowTri.setSubMatrixSize(n, n);
                    CholeskyDecomp(lowTri);
                    
                    VectorCls<RowNum, ElementType> c;
                    lowTri.ForwardSubstitution(b, c);
                    
                    MatrixCls<RowNum, RowNum, ElementType> upTri;
                    upTri.setSubMatrixSize(n, n);
                    lowTri.Transpose(upTri);
                    upTri.BackwardSubstitution(c, x);
                    return true;
                }
                else
                {
                    MatrixCls<RowNum, ColNum, ElementType> l;
                    MatrixCls<RowNum, ColNum, ElementType> u;
                    l.setSubMatrixSize(n, n);
                    u.setSubMatrixSize(n, n);
                    
                    LUDecomp(l, u);
                    
                    VectorCls<RowNum, ElementType> c;
                    l.ForwardSubstitution(b, c);
                    u.BackwardSubstitution(c, x);
                    return true;
                }
            }
        }
    }
    else
    {
        MatrixCls<RowNum, RowNum, ElementType> q;
        MatrixCls<RowNum, ColNum, ElementType> r;
        QrDecomp(q, r);
        
        // Compute Q^T * b
        VectorCls<RowNum, ElementType> qTb;
        UINT32 m = this->getRowCount();
        for (UINT32 i = 0; i < m; i++)
        {
            ElementType sum = ElementType(0);
            for (UINT32 j = 0; j < m; j++)
            {
                sum += q(j, i) * b(j); // Q^T[i,j] = Q[j,i]
            }
            qTb(i) = sum;
        }
        
        r.BackwardSubstitution(qTb, x);
        return true;
    }
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
bool MatrixCls<RowNum, ColNum, ElementType>::LeftDivide(const MatrixCls<RowNum, ColNum, ElementType>& b, 
                                                        MatrixCls<RowNum, ColNum, ElementType>& x)
{
    UINT32 m = this->getRowCount();  // Number of rows in A
    UINT32 n = this->getColCount();  // Number of columns in A
    UINT32 bCols = b.getColCount();  // Number of columns in B
    
    // Check dimensions: B must have same number of rows as A
    if (m != b.getRowCount())
    {
        return false; // Dimension mismatch
    }
    
    // According to dimension rules: x should be n x bCols
    // where n is the number of columns in A
    x.setSubMatrixSize(n, bCols);
    
    // Solve for each column using the vector LeftDivide
    VectorCls<RowNum, ElementType> b_col;
    VectorCls<ColNum, ElementType> x_col;  // x_col should have ColNum elements
    for (UINT32 col = 0; col < bCols; col++)
    {
        // Extract column from b
        for (UINT32 i = 0; i < m; i++)
        {
            b_col(i) = b(i, col);
        }
        
        // Use vector LeftDivide for each column
        // Now that LeftDivide accepts VectorCls<ColNum>, we can use x_col directly
        // The vector LeftDivide handles both square and non-square cases internally
        if (!LeftDivide(b_col, x_col))
        {
            return false;
        }
        
        // Store result in corresponding column of x
        for (UINT32 i = 0; i < n; i++)
        {
            x(i, col) = x_col(i);
        }
    }
    
    return true;
}

#endif