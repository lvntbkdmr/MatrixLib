#pragma once

#include <cstddef>
#include <stdexcept>

// Forward declarations
template <size_t RowNum, size_t ColNum, typename ElementType>
class MatrixCls;

template <size_t Length, typename ElementType>
class VectorCls;

template <typename ElementType>
class SubVecPtrCls;

// Include VectorCls and SubVecPtrCls before class definition
// (needed for member function declarations that use these types)
#include "VectorCls.h"
#include "SubVecPtrCls.h"

// Sub matrix view class - provides a view into a portion of a matrix
// Now supports runtime dimensions
template <typename ElementType>
class SubMatPtrCls
{
public:
    // Constructor - takes pointer to the start of the submatrix, stride, and runtime dimensions
    SubMatPtrCls(ElementType* data, size_t stride, size_t rows, size_t cols);
    
    // Copy constructor
    SubMatPtrCls(const SubMatPtrCls& other);
    
    // Assignment operator
    SubMatPtrCls& operator=(const SubMatPtrCls& other);
    
    // Assignment from MatrixCls (template to support any size)
    // Forward declaration needed - implementation after MatrixCls is included
    template <size_t RowNum, size_t ColNum>
    SubMatPtrCls& operator=(const MatrixCls<RowNum, ColNum, ElementType>& other);
    
    // Compound assignment operators
    SubMatPtrCls& operator+=(const SubMatPtrCls& other);
    SubMatPtrCls& operator-=(const SubMatPtrCls& other);
    SubMatPtrCls& operator*=(const ElementType& scalar);
    SubMatPtrCls& operator/=(const ElementType& scalar);
    
    // Element access
    ElementType& operator()(size_t row, size_t col);
    const ElementType& operator()(size_t row, size_t col) const;
    
    // Utility functions
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    
    // Get pointer to element (for direct access)
    ElementType* ptr(size_t row, size_t col) const;
    
    // Transpose - returns a MatrixCls with transposed data
    // Note: Since MatrixCls requires compile-time dimensions, this is a template function
    template <size_t RowNum, size_t ColNum>
    MatrixCls<ColNum, RowNum, ElementType> transpose() const;
    
    // Forward Substitution - calls generic forwardSubstitution function from MatrixPkg
    // Solves Lx = b where this matrix is L (lower triangular)
    // Overload that writes result directly to a SubVecPtrCls (for runtime-sized vectors)
    void forwardSubstitution(const SubVecPtrCls<ElementType>& b, SubVecPtrCls<ElementType>& x) const;
    
    // Overload that returns a proxy object for assignment to SubVecPtrCls
    // Usage: x_sub = L_sub.forwardSubstitution(b_sub);
    class ForwardSubstitutionResult;
    ForwardSubstitutionResult forwardSubstitution(const SubVecPtrCls<ElementType>& b) const;
    
    // Overload that returns VectorCls (requires compile-time size)
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> forwardSubstitution(const SubVecPtrCls<ElementType>& b) const;
    
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> forwardSubstitution(const VectorCls<VecSize, ElementType>& b) const;
    
    // Backward Substitution - calls generic backwardSubstitution function from MatrixPkg
    // Solves Ux = b where this matrix is U (upper triangular)
    // Overload that writes result directly to a SubVecPtrCls (for runtime-sized vectors)
    void backwardSubstitution(const SubVecPtrCls<ElementType>& b, SubVecPtrCls<ElementType>& x) const;
    
    // Overload that returns a proxy object for assignment to SubVecPtrCls
    // Usage: x_sub = L_sub.backwardSubstitution(b_sub);
    class BackwardSubstitutionResult;
    BackwardSubstitutionResult backwardSubstitution(const SubVecPtrCls<ElementType>& b) const;
    
    // Overload that returns VectorCls (requires compile-time size)
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> backwardSubstitution(const SubVecPtrCls<ElementType>& b) const;
    
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> backwardSubstitution(const VectorCls<VecSize, ElementType>& b) const;
    
    // Left Divide - calls generic leftDivide function from MatrixPkg
    // Solves Ax = b for x (equivalent to A\b in MATLAB)
    // Note: For avionics safety (no dynamic allocation), template versions are provided
    // These require compile-time size specification
    
    // Overloads that write result directly to output parameter (avionics-safe)
    template <size_t VecSize>
    void leftDivide(const SubVecPtrCls<ElementType>& b, SubVecPtrCls<ElementType>& x) const;
    
    // Overload that writes to a VectorCls - size is inferred from runtime dimensions
    // The OutputSize template parameter is automatically deduced from the VectorCls type
    // The actual problem size is determined at runtime from the matrix/vector dimensions
    template <size_t OutputSize>
    void leftDivide(const SubVecPtrCls<ElementType>& b, VectorCls<OutputSize, ElementType>& x) const;
    
    // Overload that writes to a VectorCls with explicit size (for when you want to specify VecSize explicitly)
    template <size_t VecSize, size_t OutputSize>
    void leftDivide(const SubVecPtrCls<ElementType>& b, VectorCls<OutputSize, ElementType>& x) const;
    
    template <size_t VecSize>
    void leftDivide(const VectorCls<VecSize, ElementType>& b, SubVecPtrCls<ElementType>& x) const;
    
    template <size_t VecSize>
    void leftDivide(const VectorCls<VecSize, ElementType>& b, VectorCls<VecSize, ElementType>& x) const;
    
    // Overloads that return VectorCls (requires compile-time size)
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> leftDivide(const SubVecPtrCls<ElementType>& b) const;
    
    template <size_t VecSize>
    VectorCls<VecSize, ElementType> leftDivide(const VectorCls<VecSize, ElementType>& b) const;
    
    // Helper class for assignment syntax: x_sub = L_sub.forwardSubstitution(b_sub);
    class ForwardSubstitutionResult
    {
    public:
        ForwardSubstitutionResult(const SubMatPtrCls<ElementType>* mat, const SubVecPtrCls<ElementType>* vec)
            : mat_(mat), vec_(vec) {}
        
        // Assignment operator that computes and writes to SubVecPtrCls
        void operator=(SubVecPtrCls<ElementType>& x) const
        {
            mat_->forwardSubstitution(*vec_, x);
        }
        
    private:
        const SubMatPtrCls<ElementType>* mat_;
        const SubVecPtrCls<ElementType>* vec_;
    };
    
    // Helper class for backward substitution assignment syntax
    class BackwardSubstitutionResult
    {
    public:
        BackwardSubstitutionResult(const SubMatPtrCls<ElementType>* mat, const SubVecPtrCls<ElementType>* vec)
            : mat_(mat), vec_(vec) {}
        
        // Assignment operator that computes and writes to SubVecPtrCls
        void operator=(SubVecPtrCls<ElementType>& x) const
        {
            mat_->backwardSubstitution(*vec_, x);
        }
        
    private:
        const SubMatPtrCls<ElementType>* mat_;
        const SubVecPtrCls<ElementType>* vec_;
    };
    
    // Helper class for left divide assignment syntax
    class LeftDivideResult
    {
    public:
        LeftDivideResult(const SubMatPtrCls<ElementType>* mat, const SubVecPtrCls<ElementType>* vec)
            : mat_(mat), vec_(vec) {}
        
        // Assignment operator that computes and writes to SubVecPtrCls
        void operator=(SubVecPtrCls<ElementType>& x) const
        {
            mat_->leftDivide(*vec_, x);
        }
        
    private:
        const SubMatPtrCls<ElementType>* mat_;
        const SubVecPtrCls<ElementType>* vec_;
    };

private:
    ElementType* data_;
    size_t stride_;
    size_t rows_;
    size_t cols_;
};

// Include MatrixCls for implementations
#include "MatrixCls.h"

// Include MatrixPkg for generic functions (must be after MatrixCls)
#include "MatrixPkg.h"

// Template implementations
template <typename ElementType>
SubMatPtrCls<ElementType>::SubMatPtrCls(ElementType* data, size_t stride, size_t rows, size_t cols)
    : data_(data), stride_(stride), rows_(rows), cols_(cols)
{
    if (data_ == nullptr) {
        throw std::invalid_argument("Data pointer cannot be null");
    }
    if (rows == 0 || cols == 0) {
        throw std::invalid_argument("Submatrix dimensions must be greater than zero");
    }
}

template <typename ElementType>
SubMatPtrCls<ElementType>::SubMatPtrCls(const SubMatPtrCls& other)
    : data_(other.data_), stride_(other.stride_), rows_(other.rows_), cols_(other.cols_)
{
}

template <typename ElementType>
SubMatPtrCls<ElementType>& SubMatPtrCls<ElementType>::operator=(const SubMatPtrCls& other)
{
    if (this != &other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::runtime_error("Cannot assign submatrices with different dimensions");
        }
        // Copy element by element
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                data_[i * stride_ + j] = other.data_[i * other.stride_ + j];
            }
        }
    }
    return *this;
}

template <typename ElementType>
template <size_t RowNum, size_t ColNum>
SubMatPtrCls<ElementType>& SubMatPtrCls<ElementType>::operator=(const MatrixCls<RowNum, ColNum, ElementType>& other)
{
    if (rows_ != RowNum || cols_ != ColNum) {
        throw std::runtime_error("Cannot assign matrix with different dimensions to submatrix view");
    }
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            data_[i * stride_ + j] = other.Elements[i][j];
        }
    }
    return *this;
}

template <typename ElementType>
SubMatPtrCls<ElementType>& SubMatPtrCls<ElementType>::operator+=(const SubMatPtrCls& other)
{
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::runtime_error("Cannot add submatrices with different dimensions");
    }
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            data_[i * stride_ + j] += other.data_[i * other.stride_ + j];
        }
    }
    return *this;
}

template <typename ElementType>
SubMatPtrCls<ElementType>& SubMatPtrCls<ElementType>::operator-=(const SubMatPtrCls& other)
{
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::runtime_error("Cannot subtract submatrices with different dimensions");
    }
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            data_[i * stride_ + j] -= other.data_[i * other.stride_ + j];
        }
    }
    return *this;
}

template <typename ElementType>
SubMatPtrCls<ElementType>& SubMatPtrCls<ElementType>::operator*=(const ElementType& scalar)
{
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            data_[i * stride_ + j] *= scalar;
        }
    }
    return *this;
}

template <typename ElementType>
SubMatPtrCls<ElementType>& SubMatPtrCls<ElementType>::operator/=(const ElementType& scalar)
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            data_[i * stride_ + j] /= scalar;
        }
    }
    return *this;
}

template <typename ElementType>
ElementType& SubMatPtrCls<ElementType>::operator()(size_t row, size_t col)
{
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("SubMatPtr index out of range");
    }
    return data_[row * stride_ + col];
}

template <typename ElementType>
const ElementType& SubMatPtrCls<ElementType>::operator()(size_t row, size_t col) const
{
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("SubMatPtr index out of range");
    }
    return data_[row * stride_ + col];
}

template <typename ElementType>
ElementType* SubMatPtrCls<ElementType>::ptr(size_t row, size_t col) const
{
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("SubMatPtr index out of range");
    }
    return &data_[row * stride_ + col];
}

template <typename ElementType>
template <size_t RowNum, size_t ColNum>
MatrixCls<ColNum, RowNum, ElementType> SubMatPtrCls<ElementType>::transpose() const
{
    // Call generic transpose function from MatrixPkg
    // MatrixPkg.h is included via MatrixCls.h
    // Use :: to explicitly call the global function, not the member function
    return ::transpose<RowNum, ColNum>(*this);
}

// Overload that writes result directly to a SubVecPtrCls (for runtime-sized vectors)
template <typename ElementType>
void SubMatPtrCls<ElementType>::forwardSubstitution(const SubVecPtrCls<ElementType>& b, SubVecPtrCls<ElementType>& x) const
{
    if (rows_ != cols_) {
        throw std::runtime_error("Matrix must be square for forwardSubstitution");
    }
    if (rows_ != b.size() || rows_ != x.size()) {
        throw std::runtime_error("Matrix and vector sizes must match for forwardSubstitution");
    }
    
    // Direct implementation for runtime-sized vectors
    for (size_t i = 0; i < rows_; ++i) {
        ElementType sum = b[i];
        for (size_t j = 0; j < i; ++j) {
            sum -= (*this)(i, j) * x[j];
        }
        
        if ((*this)(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in forward substitution");
        }
        
        x[i] = sum / (*this)(i, i);
    }
}

// Overload that returns a proxy object for assignment syntax
template <typename ElementType>
typename SubMatPtrCls<ElementType>::ForwardSubstitutionResult SubMatPtrCls<ElementType>::forwardSubstitution(const SubVecPtrCls<ElementType>& b) const
{
    return ForwardSubstitutionResult(this, &b);
}

// Overload that returns VectorCls (requires compile-time size)
template <typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> SubMatPtrCls<ElementType>::forwardSubstitution(const SubVecPtrCls<ElementType>& b) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match vector size for forwardSubstitution");
    }
    if (b.size() != VecSize) {
        throw std::runtime_error("Vector size must match VecSize template parameter");
    }
    
    // Call generic forwardSubstitution function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::forwardSubstitution<VecSize>(*this, b);
}

template <typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> SubMatPtrCls<ElementType>::forwardSubstitution(const VectorCls<VecSize, ElementType>& b) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match vector size for forwardSubstitution");
    }
    
    // Call generic forwardSubstitution function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::forwardSubstitution(*this, b);
}

// Overload that writes result directly to a SubVecPtrCls (for runtime-sized vectors)
template <typename ElementType>
void SubMatPtrCls<ElementType>::backwardSubstitution(const SubVecPtrCls<ElementType>& b, SubVecPtrCls<ElementType>& x) const
{
    if (rows_ != cols_) {
        throw std::runtime_error("Matrix must be square for backwardSubstitution");
    }
    if (rows_ != b.size() || rows_ != x.size()) {
        throw std::runtime_error("Matrix and vector sizes must match for backwardSubstitution");
    }
    
    // Direct implementation for runtime-sized vectors
    for (int i = static_cast<int>(rows_) - 1; i >= 0; --i) {
        ElementType sum = b[i];
        for (size_t j = i + 1; j < rows_; ++j) {
            sum -= (*this)(i, j) * x[j];
        }
        
        if ((*this)(i, i) == ElementType(0)) {
            throw std::runtime_error("Singular matrix in backward substitution");
        }
        
        x[i] = sum / (*this)(i, i);
    }
}

// Overload that returns a proxy object for assignment syntax
template <typename ElementType>
typename SubMatPtrCls<ElementType>::BackwardSubstitutionResult SubMatPtrCls<ElementType>::backwardSubstitution(const SubVecPtrCls<ElementType>& b) const
{
    return BackwardSubstitutionResult(this, &b);
}

// Overload that returns VectorCls (requires compile-time size)
template <typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> SubMatPtrCls<ElementType>::backwardSubstitution(const SubVecPtrCls<ElementType>& b) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match vector size for backwardSubstitution");
    }
    if (b.size() != VecSize) {
        throw std::runtime_error("Vector size must match VecSize template parameter");
    }
    
    // Call generic backwardSubstitution function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::backwardSubstitution<VecSize>(*this, b);
}

template <typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> SubMatPtrCls<ElementType>::backwardSubstitution(const VectorCls<VecSize, ElementType>& b) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match vector size for backwardSubstitution");
    }
    
    // Call generic backwardSubstitution function from MatrixPkg
    // Use :: to explicitly call the global function, not the member function
    return ::backwardSubstitution(*this, b);
}

// Overloads that write result directly to output parameter (avionics-safe)
// These copy data to compile-time sized matrices, solve, and write back
template <typename ElementType>
template <size_t VecSize>
void SubMatPtrCls<ElementType>::leftDivide(const SubVecPtrCls<ElementType>& b, SubVecPtrCls<ElementType>& x) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match VecSize for leftDivide");
    }
    if (rows_ != b.size() || rows_ != x.size()) {
        throw std::runtime_error("Matrix and vector sizes must match for leftDivide");
    }
    
    // Create temporary compile-time sized matrices and vectors
    MatrixCls<VecSize, VecSize, ElementType> A;
    VectorCls<VecSize, ElementType> b_vec;
    
    // Copy submatrix data
    for (size_t i = 0; i < VecSize; ++i) {
        for (size_t j = 0; j < VecSize; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Copy subvector data
    for (size_t i = 0; i < VecSize; ++i) {
        b_vec[i] = b[i];
    }
    
    // Solve using generic leftDivide (avionics-safe, uses stack-allocated matrices)
    VectorCls<VecSize, ElementType> x_result = ::leftDivide(A, b_vec);
    
    // Write result back to SubVecPtrCls
    for (size_t i = 0; i < VecSize; ++i) {
        x[i] = x_result[i];
    }
}

// Overload that infers size from runtime dimensions (no template parameter needed)
template <typename ElementType>
template <size_t OutputSize>
void SubMatPtrCls<ElementType>::leftDivide(const SubVecPtrCls<ElementType>& b, VectorCls<OutputSize, ElementType>& x) const
{
    // Get runtime size from matrix dimensions
    size_t runtime_size = rows_;
    
    if (rows_ != cols_) {
        throw std::runtime_error("leftDivide requires a square matrix");
    }
    if (rows_ != b.size()) {
        throw std::runtime_error("Matrix and vector sizes must match for leftDivide");
    }
    if (runtime_size > OutputSize) {
        throw std::runtime_error("Runtime size exceeds output vector size");
    }
    
    // Dispatch to appropriate compile-time size based on runtime size
    // We use OutputSize as the maximum and check at runtime
    if (runtime_size == 0) {
        throw std::runtime_error("Matrix size must be greater than zero");
    }
    
    // Use OutputSize for stack allocation (avionics-safe)
    // We need to work with the actual runtime_size, but allocate OutputSize x OutputSize
    // Create a MatrixCls of the exact runtime size by using a switch/if-else chain
    // For now, we'll use OutputSize and only work with the first runtime_size x runtime_size portion
    
    // Initialize matrices to zero
    MatrixCls<OutputSize, OutputSize, ElementType> A;
    VectorCls<OutputSize, ElementType> b_vec;
    
    // Zero out the entire matrix first
    for (size_t i = 0; i < OutputSize; ++i) {
        for (size_t j = 0; j < OutputSize; ++j) {
            A(i, j) = ElementType(0);
        }
    }
    
    // Copy submatrix data (only the actual runtime_size x runtime_size portion)
    for (size_t i = 0; i < runtime_size; ++i) {
        for (size_t j = 0; j < runtime_size; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Set the unused portion to identity matrix (diagonal = 1, off-diagonal = 0)
    // This ensures the unused portion doesn't affect the solution
    for (size_t i = runtime_size; i < OutputSize; ++i) {
        for (size_t j = 0; j < OutputSize; ++j) {
            if (i == j) {
                A(i, j) = ElementType(1);  // Identity on diagonal
            } else {
                A(i, j) = ElementType(0);
            }
        }
    }
    
    // Copy subvector data (zero out the rest)
    for (size_t i = 0; i < OutputSize; ++i) {
        if (i < runtime_size) {
            b_vec[i] = b[i];
        } else {
            b_vec[i] = ElementType(0);
        }
    }
    
    // Solve using generic leftDivide (avionics-safe, uses stack-allocated matrices)
    // The identity in the unused portion ensures it doesn't affect the solution
    VectorCls<OutputSize, ElementType> x_result = ::leftDivide(A, b_vec);
    
    // Write result to first runtime_size elements of output vector
    for (size_t i = 0; i < runtime_size; ++i) {
        x[i] = x_result[i];
    }
}

// Overload with explicit VecSize (for when you want to specify size explicitly)
template <typename ElementType>
template <size_t VecSize, size_t OutputSize>
void SubMatPtrCls<ElementType>::leftDivide(const SubVecPtrCls<ElementType>& b, VectorCls<OutputSize, ElementType>& x) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match VecSize for leftDivide");
    }
    if (rows_ != b.size()) {
        throw std::runtime_error("Matrix and vector sizes must match for leftDivide");
    }
    if (OutputSize < VecSize) {
        throw std::runtime_error("Output vector size must be at least VecSize");
    }
    
    // Create temporary compile-time sized matrices and vectors
    MatrixCls<VecSize, VecSize, ElementType> A;
    VectorCls<VecSize, ElementType> b_vec;
    
    // Copy submatrix data
    for (size_t i = 0; i < VecSize; ++i) {
        for (size_t j = 0; j < VecSize; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Copy subvector data
    for (size_t i = 0; i < VecSize; ++i) {
        b_vec[i] = b[i];
    }
    
    // Solve using generic leftDivide (avionics-safe, uses stack-allocated matrices)
    VectorCls<VecSize, ElementType> x_result = ::leftDivide(A, b_vec);
    
    // Write result to first VecSize elements of output vector
    for (size_t i = 0; i < VecSize; ++i) {
        x[i] = x_result[i];
    }
}

template <typename ElementType>
template <size_t VecSize>
void SubMatPtrCls<ElementType>::leftDivide(const VectorCls<VecSize, ElementType>& b, SubVecPtrCls<ElementType>& x) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match VecSize for leftDivide");
    }
    if (rows_ != x.size()) {
        throw std::runtime_error("Matrix and output vector sizes must match for leftDivide");
    }
    
    // Create temporary compile-time sized matrix
    MatrixCls<VecSize, VecSize, ElementType> A;
    
    // Copy submatrix data
    for (size_t i = 0; i < VecSize; ++i) {
        for (size_t j = 0; j < VecSize; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Solve using generic leftDivide (avionics-safe, uses stack-allocated matrices)
    VectorCls<VecSize, ElementType> x_result = ::leftDivide(A, b);
    
    // Write result back to SubVecPtrCls
    for (size_t i = 0; i < VecSize; ++i) {
        x[i] = x_result[i];
    }
}

template <typename ElementType>
template <size_t VecSize>
void SubMatPtrCls<ElementType>::leftDivide(const VectorCls<VecSize, ElementType>& b, VectorCls<VecSize, ElementType>& x) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match VecSize for leftDivide");
    }
    
    // Create temporary compile-time sized matrix
    MatrixCls<VecSize, VecSize, ElementType> A;
    
    // Copy submatrix data
    for (size_t i = 0; i < VecSize; ++i) {
        for (size_t j = 0; j < VecSize; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Solve using generic leftDivide (avionics-safe, uses stack-allocated matrices)
    x = ::leftDivide(A, b);
}

// Overload that returns VectorCls (requires compile-time size)
// Note: For avionics safety (no dynamic allocation), this version copies data to
// compile-time sized matrices and uses the safe implementation
template <typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> SubMatPtrCls<ElementType>::leftDivide(const SubVecPtrCls<ElementType>& b) const
{
    if (rows_ != VecSize) {
        throw std::runtime_error("Matrix rows must match vector size for leftDivide");
    }
    if (rows_ != cols_) {
        throw std::runtime_error("leftDivide requires a square matrix for SubMatPtrCls");
    }
    if (b.size() != VecSize) {
        throw std::runtime_error("Vector size must match VecSize template parameter");
    }
    
    // Create a temporary MatrixCls and VectorCls, copy data, solve, return result
    MatrixCls<VecSize, VecSize, ElementType> A;
    VectorCls<VecSize, ElementType> b_vec;
    
    // Copy submatrix data
    for (size_t i = 0; i < VecSize; ++i) {
        for (size_t j = 0; j < VecSize; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Copy subvector data
    for (size_t i = 0; i < VecSize; ++i) {
        b_vec[i] = b[i];
    }
    
    // Solve using generic leftDivide
    VectorCls<VecSize, ElementType> x = ::leftDivide(A, b_vec);
    
    return x;
}

template <typename ElementType>
template <size_t VecSize>
VectorCls<VecSize, ElementType> SubMatPtrCls<ElementType>::leftDivide(const VectorCls<VecSize, ElementType>& b) const
{
    if (rows_ != VecSize || cols_ != VecSize) {
        throw std::runtime_error("Matrix dimensions must match vector size for leftDivide");
    }
    
    // Create a temporary MatrixCls, copy data, solve
    MatrixCls<VecSize, VecSize, ElementType> A;
    for (size_t i = 0; i < VecSize; ++i) {
        for (size_t j = 0; j < VecSize; ++j) {
            A(i, j) = (*this)(i, j);
        }
    }
    
    // Call generic leftDivide function
    // Use :: to explicitly call the global function, not the member function
    return ::leftDivide(A, b);
}
