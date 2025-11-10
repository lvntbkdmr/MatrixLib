#pragma once

#include <cstddef>
#include <stdexcept>

// Forward declarations
template <typename ElementType>
class SubMatPtrCls;

template <size_t Length, typename ElementType>
class VectorCls;

template <size_t RowNum, size_t ColNum, typename ElementType>
class MatrixCls
{
public:
    // Constructors
    MatrixCls();
    MatrixCls(const ElementType& initialValue);
    MatrixCls(const ElementType data[RowNum][ColNum]);
    
    // Copy constructor
    MatrixCls(const MatrixCls& other);
    
    // Assignment operator
    MatrixCls& operator=(const MatrixCls& other);
    
    // Arithmetic operators
    MatrixCls operator+(const MatrixCls& other) const;
    MatrixCls operator-(const MatrixCls& other) const;
    MatrixCls operator*(const ElementType& scalar) const;
    MatrixCls operator/(const ElementType& scalar) const;
    
    // Compound assignment operators
    MatrixCls& operator+=(const MatrixCls& other);
    MatrixCls& operator-=(const MatrixCls& other);
    MatrixCls& operator*=(const ElementType& scalar);
    MatrixCls& operator/=(const ElementType& scalar);
    
    // Matrix multiplication
    template <size_t OtherColNum>
    MatrixCls<RowNum, OtherColNum, ElementType> operator*(const MatrixCls<ColNum, OtherColNum, ElementType>& other) const;
    
    // Element access
    ElementType& operator()(size_t row, size_t col);
    const ElementType& operator()(size_t row, size_t col) const;
    
    // Runtime-dimension submatrix view
    // Use this when you need variable-sized submatrices at runtime
    SubMatPtrCls<ElementType> submatrix(size_t rowStart, size_t colStart, size_t numRows, size_t numCols) const;
    
    // Convenience operator for runtime submatrix (rowStart, colStart, rowEnd, colEnd)
    // Returns a view of the submatrix from (rowStart, colStart) to (rowEnd, colEnd) inclusive
    SubMatPtrCls<ElementType> operator()(size_t rowStart, size_t colStart, size_t rowEnd, size_t colEnd) const;
    // Access to elements array
    ElementType Elements[RowNum][ColNum];
    
    // Utility functions
    size_t rows() const { return RowNum; }
    size_t cols() const { return ColNum; }
    
    // Transpose - calls generic transpose function from MatrixPkg
    MatrixCls<ColNum, RowNum, ElementType> transpose() const;
    
    // Forward Substitution - calls generic forwardSubstitution function from MatrixPkg
    // Solves Lx = b where this matrix is L (lower triangular)
    // Only valid for square matrices
    VectorCls<RowNum, ElementType> forwardSubstitution(const VectorCls<RowNum, ElementType>& b) const;
    
    // Backward Substitution - calls generic backwardSubstitution function from MatrixPkg
    // Solves Ux = b where this matrix is U (upper triangular)
    // Only valid for square matrices
    VectorCls<RowNum, ElementType> backwardSubstitution(const VectorCls<RowNum, ElementType>& b) const;
    
    // Left Divide - calls generic leftDivide function from MatrixPkg
    // Solves Ax = b for x (equivalent to A\b in MATLAB)
    // For square matrices: uses LU decomposition
    // For rectangular matrices: uses least squares solution
    
    // Overload that writes result directly to output parameter
    template <size_t VecSize>
    void leftDivide(const VectorCls<VecSize, ElementType>& b, VectorCls<VecSize, ElementType>& x) const;
    
    // Overload that returns VectorCls
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> leftDivide(const VectorCls<VecSize, ElementType>& b) const;
};

// Scalar multiplication (left side)
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType> operator*(const ElementType& scalar, const MatrixCls<RowNum, ColNum, ElementType>& matrix);

// Template implementations
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls()
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] = ElementType(0);
        }
    }
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& initialValue)
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] = initialValue;
        }
    }
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType data[RowNum][ColNum])
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] = data[i][j];
        }
    }
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const MatrixCls& other)
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] = other.Elements[i][j];
        }
    }
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>& MatrixCls<RowNum, ColNum, ElementType>::operator=(const MatrixCls& other)
{
    if (this != &other) {
        for (size_t i = 0; i < RowNum; ++i) {
            for (size_t j = 0; j < ColNum; ++j) {
                Elements[i][j] = other.Elements[i][j];
            }
        }
    }
    return *this;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator+(const MatrixCls& other) const
{
    MatrixCls result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            result.Elements[i][j] = Elements[i][j] + other.Elements[i][j];
        }
    }
    return result;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator-(const MatrixCls& other) const
{
    MatrixCls result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            result.Elements[i][j] = Elements[i][j] - other.Elements[i][j];
        }
    }
    return result;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator*(const ElementType& scalar) const
{
    MatrixCls result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            result.Elements[i][j] = Elements[i][j] * scalar;
        }
    }
    return result;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator/(const ElementType& scalar) const
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    MatrixCls result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            result.Elements[i][j] = Elements[i][j] / scalar;
        }
    }
    return result;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>& MatrixCls<RowNum, ColNum, ElementType>::operator+=(const MatrixCls& other)
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] += other.Elements[i][j];
        }
    }
    return *this;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>& MatrixCls<RowNum, ColNum, ElementType>::operator-=(const MatrixCls& other)
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] -= other.Elements[i][j];
        }
    }
    return *this;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>& MatrixCls<RowNum, ColNum, ElementType>::operator*=(const ElementType& scalar)
{
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] *= scalar;
        }
    }
    return *this;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType>& MatrixCls<RowNum, ColNum, ElementType>::operator/=(const ElementType& scalar)
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            Elements[i][j] /= scalar;
        }
    }
    return *this;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
template <size_t OtherColNum>
MatrixCls<RowNum, OtherColNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator*(const MatrixCls<ColNum, OtherColNum, ElementType>& other) const
{
    MatrixCls<RowNum, OtherColNum, ElementType> result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < OtherColNum; ++j) {
            ElementType sum = ElementType(0);
            for (size_t k = 0; k < ColNum; ++k) {
                sum += Elements[i][k] * other.Elements[k][j];
            }
            result.Elements[i][j] = sum;
        }
    }
    return result;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
ElementType& MatrixCls<RowNum, ColNum, ElementType>::operator()(size_t row, size_t col)
{
    if (row >= RowNum || col >= ColNum) {
        throw std::out_of_range("Matrix index out of range");
    }
    return Elements[row][col];
}

template <size_t RowNum, size_t ColNum, typename ElementType>
const ElementType& MatrixCls<RowNum, ColNum, ElementType>::operator()(size_t row, size_t col) const
{
    if (row >= RowNum || col >= ColNum) {
        throw std::out_of_range("Matrix index out of range");
    }
    return Elements[row][col];
}

// Member function implementations that require MatrixPkg.h
// These are defined after MatrixPkg.h is included (at end of file)

template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<RowNum, ColNum, ElementType> operator*(const ElementType& scalar, const MatrixCls<RowNum, ColNum, ElementType>& matrix)
{
    return matrix * scalar;
}

// Runtime-dimension submatrix view implementation
template <size_t RowNum, size_t ColNum, typename ElementType>
SubMatPtrCls<ElementType> MatrixCls<RowNum, ColNum, ElementType>::submatrix(size_t rowStart, size_t colStart, size_t numRows, size_t numCols) const
{
    // Bounds checking
    if (rowStart + numRows > RowNum || colStart + numCols > ColNum) {
        throw std::out_of_range("Submatrix indices out of range");
    }
    if (numRows == 0 || numCols == 0) {
        throw std::invalid_argument("Submatrix dimensions must be greater than zero");
    }
    
    // Get pointer to the start of the submatrix
    ElementType* subData = const_cast<ElementType*>(&Elements[rowStart][colStart]);
    
    // Stride is the number of columns in the original matrix
    size_t stride = ColNum;
    
    // Create and return the submatrix view
    return SubMatPtrCls<ElementType>(subData, stride, numRows, numCols);
}

template <size_t RowNum, size_t ColNum, typename ElementType>
SubMatPtrCls<ElementType> MatrixCls<RowNum, ColNum, ElementType>::operator()(size_t rowStart, size_t colStart, size_t rowEnd, size_t colEnd) const
{
    // Bounds checking
    if (rowStart >= RowNum || colStart >= ColNum || rowEnd >= RowNum || colEnd >= ColNum) {
        throw std::out_of_range("Submatrix indices out of range");
    }
    if (rowStart > rowEnd || colStart > colEnd) {
        throw std::invalid_argument("Invalid submatrix range: start indices must be <= end indices");
    }
    
    // Calculate dimensions (inclusive end)
    size_t numRows = rowEnd - rowStart + 1;
    size_t numCols = colEnd - colStart + 1;
    
    return submatrix(rowStart, colStart, numRows, numCols);
}

// Include SubMatPtrCls implementation
#include "SubMatPtrCls.h"

// Include MatrixPkg for transpose, forwardSubstitution, backwardSubstitution, and leftDivide functions
// (must be after SubMatPtrCls)
#include "MatrixPkg.h"

// Member function implementations that require MatrixPkg.h
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<ColNum, RowNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::transpose() const
{
    // Call generic transpose function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::transpose(*this);
}

template <size_t RowNum, size_t ColNum, typename ElementType>
VectorCls<RowNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::forwardSubstitution(const VectorCls<RowNum, ElementType>& b) const
{
    // Only valid for square matrices
    static_assert(RowNum == ColNum, "forwardSubstitution requires a square matrix");
    
    // Call generic forwardSubstitution function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::forwardSubstitution(*this, b);
}

template <size_t RowNum, size_t ColNum, typename ElementType>
VectorCls<RowNum, ElementType> MatrixCls<RowNum, ColNum, ElementType>::backwardSubstitution(const VectorCls<RowNum, ElementType>& b) const
{
    // Only valid for square matrices
    static_assert(RowNum == ColNum, "backwardSubstitution requires a square matrix");
    
    // Call generic backwardSubstitution function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::backwardSubstitution(*this, b);
}

template <size_t RowNum, size_t ColNum, typename ElementType>
template <size_t VecSize>
void MatrixCls<RowNum, ColNum, ElementType>::leftDivide(const VectorCls<VecSize, ElementType>& b, VectorCls<VecSize, ElementType>& x) const
{
    // For square matrices, VecSize must match RowNum
    if constexpr (RowNum == ColNum) {
        static_assert(RowNum == VecSize, "For square matrices, vector size must match matrix size");
        x = ::leftDivide(*this, b);
    } else {
        // For rectangular matrices, VecSize must match RowNum (number of rows)
        static_assert(RowNum == VecSize, "Vector size must match number of matrix rows");
        x = ::leftDivide(*this, b);
    }
}

template <size_t RowNum, size_t ColNum, typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> MatrixCls<RowNum, ColNum, ElementType>::leftDivide(const VectorCls<VecSize, ElementType>& b) const
{
    // For square matrices, VecSize must match RowNum
    if constexpr (RowNum == ColNum) {
        static_assert(RowNum == VecSize, "For square matrices, vector size must match matrix size");
        return ::leftDivide(*this, b);
    } else {
        // For rectangular matrices, VecSize must match RowNum (number of rows)
        static_assert(RowNum == VecSize, "Vector size must match number of matrix rows");
        return ::leftDivide(*this, b);
    }
}
