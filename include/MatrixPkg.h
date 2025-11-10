#pragma once

#include "MatrixCls.h"
#include "VectorCls.h"
#include "SubMatPtrCls.h"
#include "SubVecPtrCls.h"
#include <cstddef>
#include <type_traits>

// Unified Forward Substitution: solves Lx = b where L is lower triangular
// Works with both MatrixCls/VectorCls and SubMatPtrCls/SubVecPtrCls
// Uses function overloading - compiler automatically selects the correct overload

// Overload for MatrixCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> forwardSubstitution(
    const MatrixCls<Size, Size, ElementType>& L,
    const VectorCls<Size, ElementType>& b);

// Overload for SubMatPtrCls with runtime SubVecPtrCls
// Note: Size must be specified as template parameter for return type
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> forwardSubstitution(
    const SubMatPtrCls<ElementType>& L,
    const SubVecPtrCls<ElementType>& b);

// Overload for SubMatPtrCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> forwardSubstitution(
    const SubMatPtrCls<ElementType>& L,
    const VectorCls<Size, ElementType>& b);

// Unified Backward Substitution: solves Ux = b where U is upper triangular
// Works with both MatrixCls/VectorCls and SubMatPtrCls/SubVecPtrCls

// Overload for MatrixCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> backwardSubstitution(
    const MatrixCls<Size, Size, ElementType>& U,
    const VectorCls<Size, ElementType>& b);

// Overload for SubMatPtrCls with runtime SubVecPtrCls
// Note: Size must be specified as template parameter for return type
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> backwardSubstitution(
    const SubMatPtrCls<ElementType>& U,
    const SubVecPtrCls<ElementType>& b);

// Overload for SubMatPtrCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> backwardSubstitution(
    const SubMatPtrCls<ElementType>& U,
    const VectorCls<Size, ElementType>& b);

// LU Decomposition: decomposes A into L and U such that A = LU
// Returns true if decomposition is successful, false otherwise
template <size_t Size, typename ElementType>
bool luDecomposition(
    const MatrixCls<Size, Size, ElementType>& A,
    MatrixCls<Size, Size, ElementType>& L,
    MatrixCls<Size, Size, ElementType>& U);

// Left Divide: solves Ax = b for x (equivalent to A\b in MATLAB)
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> leftDivide(
    const MatrixCls<Size, Size, ElementType>& A,
    const VectorCls<Size, ElementType>& b);

// Left Divide for rectangular matrices (least squares solution)
template <size_t RowNum, size_t ColNum, typename ElementType>
VectorCls<ColNum, ElementType> leftDivide(
    const MatrixCls<RowNum, ColNum, ElementType>& A,
    const VectorCls<RowNum, ElementType>& b);


// Unified transpose function that works for both MatrixCls and SubMatPtrCls
// Uses function overloading - the compiler automatically selects the correct overload
// This is the standard C++ way to have "one function" work with multiple types

// Overload for MatrixCls - dimensions extracted from template parameters
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<ColNum, RowNum, ElementType> transpose(const MatrixCls<RowNum, ColNum, ElementType>& mat);

// Overload for SubMatPtrCls - dimensions must be specified as template parameters
// (since SubMatPtrCls has runtime dimensions but MatrixCls return type needs compile-time dimensions)
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<ColNum, RowNum, ElementType> transpose(const SubMatPtrCls<ElementType>& mat);

// Template implementations
#include <cmath>
#include <stdexcept>

// Unified Forward Substitution implementation
// Overload for MatrixCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> forwardSubstitution(
    const MatrixCls<Size, Size, ElementType>& L,
    const VectorCls<Size, ElementType>& b)
{
    VectorCls<Size, ElementType> x;
    
    for (size_t i = 0; i < Size; ++i) {
        ElementType sum = b[i];
        for (size_t j = 0; j < i; ++j) {
            sum -= L(i, j) * x[j];
        }
        
        if (L(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in forward substitution");
        }
        
        x[i] = sum / L(i, i);
    }
    
    return x;
}

// Overload for SubMatPtrCls with runtime SubVecPtrCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> forwardSubstitution(
    const SubMatPtrCls<ElementType>& L,
    const SubVecPtrCls<ElementType>& b)
{
    if (L.rows() != Size || L.cols() != Size) {
        throw std::runtime_error("Matrix dimensions must match vector size");
    }
    if (b.size() != Size) {
        throw std::runtime_error("Vector size must match specified Size template parameter");
    }
    
    VectorCls<Size, ElementType> x;
    
    for (size_t i = 0; i < Size; ++i) {
        ElementType sum = b[i];
        for (size_t j = 0; j < i; ++j) {
            sum -= L(i, j) * x[j];
        }
        
        if (L(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in forward substitution");
        }
        
        x[i] = sum / L(i, i);
    }
    
    return x;
}

// Overload for SubMatPtrCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> forwardSubstitution(
    const SubMatPtrCls<ElementType>& L,
    const VectorCls<Size, ElementType>& b)
{
    if (L.rows() != Size || L.cols() != Size) {
        throw std::runtime_error("Matrix dimensions must match vector size");
    }
    
    VectorCls<Size, ElementType> x;
    
    for (size_t i = 0; i < Size; ++i) {
        ElementType sum = b[i];
        for (size_t j = 0; j < i; ++j) {
            sum -= L(i, j) * x[j];
        }
        
        if (L(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in forward substitution");
        }
        
        x[i] = sum / L(i, i);
    }
    
    return x;
}

// Unified Backward Substitution implementation
// Overload for MatrixCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> backwardSubstitution(
    const MatrixCls<Size, Size, ElementType>& U,
    const VectorCls<Size, ElementType>& b)
{
    VectorCls<Size, ElementType> x;
    
    for (int i = static_cast<int>(Size) - 1; i >= 0; --i) {
        ElementType sum = b[i];
        for (size_t j = i + 1; j < Size; ++j) {
            sum -= U(i, j) * x[j];
        }
        
        if (U(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in backward substitution");
        }
        
        x[i] = sum / U(i, i);
    }
    
    return x;
}

// Overload for SubMatPtrCls with runtime SubVecPtrCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> backwardSubstitution(
    const SubMatPtrCls<ElementType>& U,
    const SubVecPtrCls<ElementType>& b)
{
    if (U.rows() != Size || U.cols() != Size) {
        throw std::runtime_error("Matrix dimensions must match vector size");
    }
    if (b.size() != Size) {
        throw std::runtime_error("Vector size must match specified Size template parameter");
    }
    
    VectorCls<Size, ElementType> x;
    
    for (int i = static_cast<int>(Size) - 1; i >= 0; --i) {
        ElementType sum = b[i];
        for (size_t j = i + 1; j < Size; ++j) {
            sum -= U(i, j) * x[j];
        }
        
        if (U(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in backward substitution");
        }
        
        x[i] = sum / U(i, i);
    }
    
    return x;
}

// Overload for SubMatPtrCls with VectorCls
template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> backwardSubstitution(
    const SubMatPtrCls<ElementType>& U,
    const VectorCls<Size, ElementType>& b)
{
    if (U.rows() != Size || U.cols() != Size) {
        throw std::runtime_error("Matrix dimensions must match vector size");
    }
    
    VectorCls<Size, ElementType> x;
    
    for (int i = static_cast<int>(Size) - 1; i >= 0; --i) {
        ElementType sum = b[i];
        for (size_t j = i + 1; j < Size; ++j) {
            sum -= U(i, j) * x[j];
        }
        
        if (U(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in backward substitution");
        }
        
        x[i] = sum / U(i, i);
    }
    
    return x;
}

template <size_t Size, typename ElementType>
bool luDecomposition(
    const MatrixCls<Size, Size, ElementType>& A,
    MatrixCls<Size, Size, ElementType>& L,
    MatrixCls<Size, Size, ElementType>& U)
{
    L = MatrixCls<Size, Size, ElementType>(ElementType(0));
    U = MatrixCls<Size, Size, ElementType>(ElementType(0));
    
    for (size_t i = 0; i < Size; ++i) {
        L(i, i) = ElementType(1);
    }
    
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = i; j < Size; ++j) {
            ElementType sum = ElementType(0);
            for (size_t k = 0; k < i; ++k) {
                sum += L(i, k) * U(k, j);
            }
            U(i, j) = A(i, j) - sum;
        }
        
        for (size_t j = i + 1; j < Size; ++j) {
            ElementType sum = ElementType(0);
            for (size_t k = 0; k < i; ++k) {
                sum += L(j, k) * U(k, i);
            }
            
            if (U(i, i) == ElementType(0)) {
                return false;
            }
            
            L(j, i) = (A(j, i) - sum) / U(i, i);
        }
    }
    
    return true;
}

template <size_t Size, typename ElementType>
VectorCls<Size, ElementType> leftDivide(
    const MatrixCls<Size, Size, ElementType>& A,
    const VectorCls<Size, ElementType>& b)
{
    MatrixCls<Size, Size, ElementType> L, U;
    
    if (!luDecomposition(A, L, U)) {
        throw std::runtime_error("Matrix is singular, cannot solve system");
    }
    
    VectorCls<Size, ElementType> y = forwardSubstitution(L, b);
    VectorCls<Size, ElementType> x = backwardSubstitution(U, y);
    
    return x;
}

template <size_t RowNum, size_t ColNum, typename ElementType>
VectorCls<ColNum, ElementType> leftDivide(
    const MatrixCls<RowNum, ColNum, ElementType>& A,
    const VectorCls<RowNum, ElementType>& b)
{
    auto At = A.transpose();
    MatrixCls<ColNum, ColNum, ElementType> AtA = At * A;
    
    VectorCls<ColNum, ElementType> Atb;
    for (size_t i = 0; i < ColNum; ++i) {
        ElementType sum = ElementType(0);
        for (size_t j = 0; j < RowNum; ++j) {
            sum += At(i, j) * b[j];
        }
        Atb[i] = sum;
    }
    
    return leftDivide(AtA, Atb);
}


// Unified transpose implementation - works for both MatrixCls and SubMatPtrCls
// This single template function handles both cases using the same logic
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<ColNum, RowNum, ElementType> transpose(const MatrixCls<RowNum, ColNum, ElementType>& mat)
{
    MatrixCls<ColNum, RowNum, ElementType> result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            result.Elements[j][i] = mat.Elements[i][j];
        }
    }
    return result;
}

// Overload for SubMatPtrCls - same function name, different parameter type
// The template parameters specify the expected dimensions (must match runtime dimensions)
template <size_t RowNum, size_t ColNum, typename ElementType>
MatrixCls<ColNum, RowNum, ElementType> transpose(const SubMatPtrCls<ElementType>& mat)
{
    if (mat.rows() != RowNum || mat.cols() != ColNum) {
        throw std::runtime_error("SubMatPtrCls dimensions do not match template parameters");
    }
    
    MatrixCls<ColNum, RowNum, ElementType> result;
    for (size_t i = 0; i < RowNum; ++i) {
        for (size_t j = 0; j < ColNum; ++j) {
            result.Elements[j][i] = mat(i, j);
        }
    }
    return result;
}
