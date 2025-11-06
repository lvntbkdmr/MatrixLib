#ifndef MATRIX_LIB_MATRIX_HPP
#define MATRIX_LIB_MATRIX_HPP

#include <cstdint>

// Type alias for consistency
typedef uint32_t UINT32;

namespace avionics {
namespace math {

/**
 * @brief Return codes for matrix operations
 * 
 * Following avionics standards, operations return status codes
 * instead of throwing exceptions.
 */
enum class MatrixStatus : std::int32_t {
    SUCCESS = 0,
    INVALID_DIMENSIONS = -1,
    OUT_OF_BOUNDS = -2,
    SINGULAR_MATRIX = -3,
    INCOMPATIBLE_SIZES = -4
};

/**
 * @brief Matrix class template for fixed-size matrices
 * 
 * @tparam RowNum Number of rows (compile-time constant, default: 3)
 * @tparam ColNum Number of columns (compile-time constant, default: RowNum)
 * @tparam ElementType Element type (default: double)
 * 
 * This implementation follows avionics software standards:
 * - No dynamic memory allocation
 * - Fixed-size at compile time
 * - No exceptions (uses return codes)
 * - MISRA C++ compliant subset
 * - Embedded system compatible (no I/O dependencies)
 */
template <int RowNum = 3, int ColNum = RowNum, class ElementType = double> 
class MatrixCls
{
public:
    using value_type = ElementType;
    using size_type = UINT32;
    
    static constexpr int ROW_COUNT = RowNum;
    static constexpr int COL_COUNT = ColNum;
    static constexpr int ELEMENT_COUNT = RowNum * ColNum;

    /**
     * @brief Default constructor - initializes to zero
     */
    MatrixCls() noexcept {
        // Zero-initialize array
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] = ElementType{0};
        }
    }

    /**
     * @brief Constructor with initializer array
     * 
     * @param init Initial values (row-major order)
     */
    MatrixCls(const ElementType (&init)[ELEMENT_COUNT]) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] = init[i];
        }
    }

    /**
     * @brief Copy constructor
     */
    MatrixCls(const MatrixCls& other) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] = other.data_[i];
        }
    }

    /**
     * @brief Copy assignment operator
     */
    MatrixCls& operator=(const MatrixCls& other) noexcept {
        if (this != &other) {
            for (int i = 0; i < ELEMENT_COUNT; ++i) {
                data_[i] = other.data_[i];
            }
        }
        return *this;
    }

    /**
     * @brief Access element at (row, col)
     * 
     * @param row Row index [0, RowNum)
     * @param col Column index [0, ColNum)
     * @return Reference to element
     */
    ElementType& operator()(int row, int col) noexcept {
        return data_[row * ColNum + col];
    }

    /**
     * @brief Access element at (row, col) - const version
     */
    const ElementType& operator()(int row, int col) const noexcept {
        return data_[row * ColNum + col];
    }

    /**
     * @brief Get element with bounds checking
     * 
     * @param row Row index
     * @param col Column index
     * @param value Output parameter for the element value
     * @return MatrixStatus::SUCCESS if valid, MatrixStatus::OUT_OF_BOUNDS otherwise
     */
    MatrixStatus get(int row, int col, ElementType& value) const noexcept {
        if (row < 0 || row >= RowNum || col < 0 || col >= ColNum) {
            return MatrixStatus::OUT_OF_BOUNDS;
        }
        value = data_[row * ColNum + col];
        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief Set element with bounds checking
     */
    MatrixStatus set(int row, int col, const ElementType& value) noexcept {
        if (row < 0 || row >= RowNum || col < 0 || col >= ColNum) {
            return MatrixStatus::OUT_OF_BOUNDS;
        }
        data_[row * ColNum + col] = value;
        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief Get number of rows
     */
    static constexpr int rows() noexcept {
        return RowNum;
    }

    /**
     * @brief Get number of columns
     */
    static constexpr int cols() noexcept {
        return ColNum;
    }

    /**
     * @brief Get total number of elements
     */
    static constexpr int size() noexcept {
        return ELEMENT_COUNT;
    }

    /**
     * @brief Fill matrix with a constant value
     */
    void fill(const ElementType& value) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] = value;
        }
    }

    /**
     * @brief Set matrix to identity (must be square)
     */
    MatrixStatus identity() noexcept {
        if (RowNum != ColNum) {
            return MatrixStatus::INVALID_DIMENSIONS;
        }
        fill(ElementType{0});
        for (int i = 0; i < RowNum; ++i) {
            data_[i * ColNum + i] = ElementType{1};
        }
        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief Set matrix to zero
     */
    void zero() noexcept {
        fill(ElementType{0});
    }

    /**
     * @brief Get raw data pointer (for compatibility)
     */
    ElementType* data() noexcept {
        return data_;
    }

    const ElementType* data() const noexcept {
        return data_;
    }

    // Arithmetic operators
    MatrixCls& operator+=(const MatrixCls& other) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    MatrixCls& operator-=(const MatrixCls& other) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    MatrixCls& operator*=(const ElementType& scalar) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    MatrixCls& operator/=(const ElementType& scalar) noexcept {
        for (int i = 0; i < ELEMENT_COUNT; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    /**
     * @brief LU Decomposition (Crout's algorithm)
     * 
     * Decomposes the matrix into L (lower triangular) and U (upper triangular) matrices
     * such that A = L * U
     * 
     * @param L Output parameter for lower triangular matrix
     * @param U Output parameter for upper triangular matrix
     * @return MatrixStatus::SUCCESS if successful, MatrixStatus::INVALID_DIMENSIONS if not square
     */
    MatrixStatus LuDecomp(MatrixCls<RowNum, ColNum, ElementType>& L,
                          MatrixCls<RowNum, ColNum, ElementType>& U) const noexcept {
        // Matrix must be square for LU decomposition
        if (RowNum != ColNum) {
            return MatrixStatus::INVALID_DIMENSIONS;
        }

        const UINT32 N = static_cast<UINT32>(RowNum);

        // Initialize L to zeros and U to identity
        L.zero();
        U.identity();

        // L(0,0) = A(0,0)
        L(0, 0) = (*this)(0, 0);

        // Check for division by zero
        if (L(0, 0) == ElementType{0}) {
            return MatrixStatus::SINGULAR_MATRIX;
        }

        // for j = 1 to N-1
        for (UINT32 j = 1U; j < N; ++j) {
            // L(j,0) = A(j,0)
            L(j, 0) = (*this)(j, 0);
            // U(0,j) = A(0,j) / L(0,0)
            U(0, j) = (*this)(0, j) / L(0, 0);
        }

        // for j = 1 to N-2
        for (UINT32 j = 1U; j < N - 1U; ++j) {
            // for i = j to N-1
            for (UINT32 i = j; i < N; ++i) {
                // L(i,j) = A(i,j)
                L(i, j) = (*this)(i, j);
                // for k = 0 to j-1
                for (UINT32 k = 0U; k < j; ++k) {
                    // L(i,j) = L(i,j) - L(i,k) * U(k,j)
                    L(i, j) = L(i, j) - L(i, k) * U(k, j);
                }
            }

            // Check for division by zero
            if (L(j, j) == ElementType{0}) {
                return MatrixStatus::SINGULAR_MATRIX;
            }

            // for k = j+1 to N-1
            for (UINT32 k = j + 1U; k < N; ++k) {
                // U(j,k) = A(j,k)
                U(j, k) = (*this)(j, k);
                // for i = 0 to j-1
                for (UINT32 i = 0U; i < j; ++i) {
                    // U(j,k) = U(j,k) - L(j,i) * U(i,k)
                    U(j, k) = U(j, k) - L(j, i) * U(i, k);
                }
                // U(j,k) = U(j,k) / L(j,j)
                U(j, k) = U(j, k) / L(j, j);
            }
        }

        // L(N-1,N-1) = A(N-1,N-1)
        L(N - 1U, N - 1U) = (*this)(N - 1U, N - 1U);
        // for k = 0 to N-2
        for (UINT32 k = 0U; k < N - 1U; ++k) {
            // L(N-1,N-1) = L(N-1,N-1) - L(N-1,k) * U(k,N-1)
            L(N - 1U, N - 1U) = L(N - 1U, N - 1U) - L(N - 1U, k) * U(k, N - 1U);
        }

        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief LU Decomposition with runtime size parameter
     * 
     * Decomposes the top-left (size x size) portion of the matrix into 
     * L (lower triangular) and U (upper triangular) matrices such that A = L * U.
     * 
     * This is useful when working with runtime-determined submatrix sizes.
     * Only the top-left (size x size) portion of the matrix is used for decomposition.
     * 
     * @param size Size of the submatrix to decompose (1 to min(RowNum, ColNum))
     * @param L Output parameter for lower triangular matrix (must be at least size x size)
     * @param U Output parameter for upper triangular matrix (must be at least size x size)
     * @return MatrixStatus::SUCCESS if successful, MatrixStatus::INVALID_DIMENSIONS if not square, 
     *         MatrixStatus::OUT_OF_BOUNDS if size is invalid, MatrixStatus::SINGULAR_MATRIX if singular
     */
    MatrixStatus LuDecomp(
        int size,
        MatrixCls<RowNum, ColNum, ElementType>& L,
        MatrixCls<RowNum, ColNum, ElementType>& U) const noexcept {
        // Validate size
        if (size < 1 || size > RowNum || size > ColNum) {
            return MatrixStatus::OUT_OF_BOUNDS;
        }
        
        // Matrix must be square for LU decomposition
        if (RowNum != ColNum) {
            return MatrixStatus::INVALID_DIMENSIONS;
        }

        const UINT32 N = static_cast<UINT32>(size);

        // Initialize L to zeros and U to identity
        L.zero();
        U.identity();

        // L(0,0) = A(0,0)
        L(0, 0) = (*this)(0, 0);

        // Check for division by zero
        if (L(0, 0) == ElementType{0}) {
            return MatrixStatus::SINGULAR_MATRIX;
        }

        // for j = 1 to N-1
        for (UINT32 j = 1U; j < N; ++j) {
            // L(j,0) = A(j,0)
            L(j, 0) = (*this)(j, 0);
            // U(0,j) = A(0,j) / L(0,0)
            U(0, j) = (*this)(0, j) / L(0, 0);
        }

        // for j = 1 to N-2
        for (UINT32 j = 1U; j < N - 1U; ++j) {
            // for i = j to N-1
            for (UINT32 i = j; i < N; ++i) {
                // L(i,j) = A(i,j)
                L(i, j) = (*this)(i, j);
                // for k = 0 to j-1
                for (UINT32 k = 0U; k < j; ++k) {
                    // L(i,j) = L(i,j) - L(i,k) * U(k,j)
                    L(i, j) = L(i, j) - L(i, k) * U(k, j);
                }
            }

            // Check for division by zero
            if (L(j, j) == ElementType{0}) {
                return MatrixStatus::SINGULAR_MATRIX;
            }

            // for k = j+1 to N-1
            for (UINT32 k = j + 1U; k < N; ++k) {
                // U(j,k) = A(j,k)
                U(j, k) = (*this)(j, k);
                // for i = 0 to j-1
                for (UINT32 i = 0U; i < j; ++i) {
                    // U(j,k) = U(j,k) - L(j,i) * U(i,k)
                    U(j, k) = U(j, k) - L(j, i) * U(i, k);
                }
                // U(j,k) = U(j,k) / L(j,j)
                U(j, k) = U(j, k) / L(j, j);
            }
        }

        // L(N-1,N-1) = A(N-1,N-1)
        L(N - 1U, N - 1U) = (*this)(N - 1U, N - 1U);
        // for k = 0 to N-2
        for (UINT32 k = 0U; k < N - 1U; ++k) {
            // L(N-1,N-1) = L(N-1,N-1) - L(N-1,k) * U(k,N-1)
            L(N - 1U, N - 1U) = L(N - 1U, N - 1U) - L(N - 1U, k) * U(k, N - 1U);
        }

        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief Extract top-left submatrix of specified size
     * 
     * Extracts the top-left (size x size) submatrix from this matrix.
     * Useful for runtime-determined submatrix operations.
     * 
     * @param size Size of the submatrix to extract (must be <= min(RowNum, ColNum))
     * @param result Output parameter for the extracted submatrix
     * @return MatrixStatus::SUCCESS if valid, MatrixStatus::OUT_OF_BOUNDS if size is invalid
     */
    template<int SubSize>
    MatrixStatus extractSubmatrix(int size, MatrixCls<SubSize, SubSize, ElementType>& result) const noexcept {
        if (size < 1 || size > SubSize || size > RowNum || size > ColNum) {
            return MatrixStatus::OUT_OF_BOUNDS;
        }
        
        result.zero();
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result(i, j) = (*this)(i, j);
            }
        }
        
        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief Extract top-left submatrix into a fixed-size matrix (up to 6x6)
     * 
     * This is a convenience function that extracts the top-left (size x size) submatrix
     * into a fixed-size result matrix. The result matrix must be large enough to hold
     * the extracted submatrix.
     * 
     * @param size Size of the submatrix to extract (1-6)
     * @param result Output parameter for the extracted submatrix (must be at least size x size)
     * @return MatrixStatus::SUCCESS if valid, MatrixStatus::OUT_OF_BOUNDS if size is invalid
     */
    template<int ResultSize>
    MatrixStatus extractTopLeft(int size, MatrixCls<ResultSize, ResultSize, ElementType>& result) const noexcept {
        if (size < 1 || size > ResultSize || size > RowNum || size > ColNum) {
            return MatrixStatus::OUT_OF_BOUNDS;
        }
        
        result.zero();
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result(i, j) = (*this)(i, j);
            }
        }
        
        return MatrixStatus::SUCCESS;
    }

    /**
     * @brief Set top-left submatrix from another matrix
     * 
     * Sets the top-left (size x size) portion of this matrix from the source matrix.
     * 
     * @param size Size of the submatrix to set (must be <= min(RowNum, ColNum, SourceRows, SourceCols))
     * @param source Source matrix to copy from
     * @return MatrixStatus::SUCCESS if valid, MatrixStatus::OUT_OF_BOUNDS if size is invalid
     */
    template<int SourceRows, int SourceCols>
    MatrixStatus setTopLeft(int size, const MatrixCls<SourceRows, SourceCols, ElementType>& source) noexcept {
        if (size < 1 || size > RowNum || size > ColNum || 
            size > SourceRows || size > SourceCols) {
            return MatrixStatus::OUT_OF_BOUNDS;
        }
        
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                (*this)(i, j) = source(i, j);
            }
        }
        
        return MatrixStatus::SUCCESS;
    }

private:
    ElementType data_[ELEMENT_COUNT];  // Fixed-size array (no dynamic allocation)
};

// Free function operators (function overloading)

/**
 * @brief Matrix addition
 */
template <int RowNum, int ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType> operator+(
    const MatrixCls<RowNum, ColNum, ElementType>& lhs,
    const MatrixCls<RowNum, ColNum, ElementType>& rhs) noexcept {
    MatrixCls<RowNum, ColNum, ElementType> result = lhs;
    result += rhs;
    return result;
}

/**
 * @brief Matrix subtraction
 */
template <int RowNum, int ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType> operator-(
    const MatrixCls<RowNum, ColNum, ElementType>& lhs,
    const MatrixCls<RowNum, ColNum, ElementType>& rhs) noexcept {
    MatrixCls<RowNum, ColNum, ElementType> result = lhs;
    result -= rhs;
    return result;
}

/**
 * @brief Scalar multiplication
 */
template <int RowNum, int ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType> operator*(
    const MatrixCls<RowNum, ColNum, ElementType>& matrix,
    const ElementType& scalar) noexcept {
    MatrixCls<RowNum, ColNum, ElementType> result = matrix;
    result *= scalar;
    return result;
}

template <int RowNum, int ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType> operator*(
    const ElementType& scalar,
    const MatrixCls<RowNum, ColNum, ElementType>& matrix) noexcept {
    return matrix * scalar;
}

/**
 * @brief Matrix multiplication
 * 
 * @tparam RowNum1 Rows of first matrix
 * @tparam Common Common dimension
 * @tparam ColNum2 Columns of second matrix
 * @tparam ElementType Element type
 */
template <int RowNum1, int Common, int ColNum2, class ElementType>
MatrixCls<RowNum1, ColNum2, ElementType> operator*(
    const MatrixCls<RowNum1, Common, ElementType>& lhs,
    const MatrixCls<Common, ColNum2, ElementType>& rhs) noexcept {
    MatrixCls<RowNum1, ColNum2, ElementType> result;
    
    for (int i = 0; i < RowNum1; ++i) {
        for (int j = 0; j < ColNum2; ++j) {
            ElementType sum = ElementType{0};
            for (int k = 0; k < Common; ++k) {
                sum += lhs(i, k) * rhs(k, j);
            }
            result(i, j) = sum;
        }
    }
    
    return result;
}

/**
 * @brief Matrix comparison operators
 */
template <int RowNum, int ColNum, class ElementType>
bool operator==(
    const MatrixCls<RowNum, ColNum, ElementType>& lhs,
    const MatrixCls<RowNum, ColNum, ElementType>& rhs) noexcept {
    for (int i = 0; i < MatrixCls<RowNum, ColNum, ElementType>::ELEMENT_COUNT; ++i) {
        if (lhs.data()[i] != rhs.data()[i]) {
            return false;
        }
    }
    return true;
}

template <int RowNum, int ColNum, class ElementType>
bool operator!=(
    const MatrixCls<RowNum, ColNum, ElementType>& lhs,
    const MatrixCls<RowNum, ColNum, ElementType>& rhs) noexcept {
    return !(lhs == rhs);
}

// Matrix operations (static functions)

/**
 * @brief Transpose matrix
 */
template <int RowNum, int ColNum, class ElementType>
MatrixCls<ColNum, RowNum, ElementType> transpose(
    const MatrixCls<RowNum, ColNum, ElementType>& matrix) noexcept {
    MatrixCls<ColNum, RowNum, ElementType> result;
    
    for (int i = 0; i < RowNum; ++i) {
        for (int j = 0; j < ColNum; ++j) {
            result(j, i) = matrix(i, j);
        }
    }
    
    return result;
}

/**
 * @brief Determinant for 2x2 matrix
 */
template <class ElementType>
ElementType determinant(const MatrixCls<2, 2, ElementType>& matrix) noexcept {
    return matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0);
}

/**
 * @brief Determinant for 3x3 matrix
 */
template <class ElementType>
ElementType determinant(const MatrixCls<3, 3, ElementType>& matrix) noexcept {
    return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)) -
           matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0)) +
           matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
}

/**
 * @brief Inverse for 2x2 matrix
 */
template <class ElementType>
MatrixStatus inverse(
    const MatrixCls<2, 2, ElementType>& matrix,
    MatrixCls<2, 2, ElementType>& result) noexcept {
    const ElementType det = determinant(matrix);
    
    // Check for singular matrix (avoid division by zero)
    if (det == ElementType{0}) {
        return MatrixStatus::SINGULAR_MATRIX;
    }
    
    const ElementType inv_det = ElementType{1} / det;
    
    result(0, 0) = matrix(1, 1) * inv_det;
    result(0, 1) = -matrix(0, 1) * inv_det;
    result(1, 0) = -matrix(1, 0) * inv_det;
    result(1, 1) = matrix(0, 0) * inv_det;
    
    return MatrixStatus::SUCCESS;
}

/**
 * @brief Forward substitution: Solve L*y = B for y
 * 
 * Solves a lower triangular system L*y = B using forward substitution.
 * 
 * @param L Lower triangular matrix (must be square)
 * @param B Right-hand side matrix
 * @param y Output parameter for solution
 * @return MatrixStatus::SUCCESS if successful
 */
template <int Size, int Cols, class ElementType>
MatrixStatus forwardSubstitution(
    const MatrixCls<Size, Size, ElementType>& L,
    const MatrixCls<Size, Cols, ElementType>& B,
    MatrixCls<Size, Cols, ElementType>& y) noexcept {
    const UINT32 N = static_cast<UINT32>(Size);
    
    // Initialize y with B
    y = B;
    
    // Forward substitution: solve L*y = B
    // For each column in B
    for (UINT32 col = 0U; col < static_cast<UINT32>(Cols); ++col) {
        // For each row
        for (UINT32 i = 0U; i < N; ++i) {
            ElementType sum = ElementType{0};
            // Sum L(i,k) * y(k,col) for k = 0 to i-1
            for (UINT32 k = 0U; k < i; ++k) {
                sum += L(i, k) * y(k, col);
            }
            // Check for division by zero
            if (L(i, i) == ElementType{0}) {
                return MatrixStatus::SINGULAR_MATRIX;
            }
            // y(i,col) = (B(i,col) - sum) / L(i,i)
            y(i, col) = (B(i, col) - sum) / L(i, i);
        }
    }
    
    return MatrixStatus::SUCCESS;
}

/**
 * @brief Forward substitution with runtime size parameter
 * 
 * Solves a lower triangular system L*y = B for y, but only operates on
 * the top-left (size x size) portion of the matrices.
 * 
 * This is useful when working with runtime-determined submatrix sizes.
 * 
 * @param L Lower triangular matrix (must be at least size x size)
 * @param B Right-hand side matrix (must be at least size x BCols)
 * @param size Size of the submatrix to use (1 to min(LRows, LCols, BRows))
 * @param y Output parameter for solution (must be at least size x BCols)
 * @return MatrixStatus::SUCCESS if successful, MatrixStatus::OUT_OF_BOUNDS if size is invalid
 */
template <int LRows, int LCols, int BRows, int BCols, class ElementType>
MatrixStatus forwardSubstitution(
    const MatrixCls<LRows, LCols, ElementType>& L,
    const MatrixCls<BRows, BCols, ElementType>& B,
    int size,
    MatrixCls<BRows, BCols, ElementType>& y) noexcept {
    // Validate size
    if (size < 1 || size > LRows || size > LCols || size > BRows) {
        return MatrixStatus::OUT_OF_BOUNDS;
    }
    
    // Initialize y with B (only top-left size×BCols portion)
    y.zero();
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < BCols; ++j) {
            y(i, j) = B(i, j);
        }
    }
    
    // Forward substitution: solve L*y = B (only top-left size×size portion)
    // For each column in B
    for (int col = 0; col < BCols; ++col) {
        // For each row (only top-left size×size portion)
        for (int i = 0; i < size; ++i) {
            ElementType sum = ElementType{0};
            // Sum L(i,k) * y(k,col) for k = 0 to i-1
            for (int k = 0; k < i; ++k) {
                sum += L(i, k) * y(k, col);
            }
            // Check for division by zero
            if (L(i, i) == ElementType{0}) {
                return MatrixStatus::SINGULAR_MATRIX;
            }
            // y(i,col) = (B(i,col) - sum) / L(i,i)
            y(i, col) = (B(i, col) - sum) / L(i, i);
        }
    }
    
    return MatrixStatus::SUCCESS;
}

/**
 * @brief Backward substitution: Solve U*x = y for x
 * 
 * Solves an upper triangular system U*x = y using backward substitution.
 * 
 * @param U Upper triangular matrix (must be square)
 * @param y Right-hand side matrix
 * @param x Output parameter for solution
 * @return MatrixStatus::SUCCESS if successful
 */
template <int Size, int Cols, class ElementType>
MatrixStatus backwardSubstitution(
    const MatrixCls<Size, Size, ElementType>& U,
    const MatrixCls<Size, Cols, ElementType>& y,
    MatrixCls<Size, Cols, ElementType>& x) noexcept {
    const UINT32 N = static_cast<UINT32>(Size);
    
    // Initialize x with y
    x = y;
    
    // Backward substitution: solve U*x = y
    // For each column in y
    for (UINT32 col = 0U; col < static_cast<UINT32>(Cols); ++col) {
        // For each row (from bottom to top)
        for (int i = static_cast<int>(N) - 1; i >= 0; --i) {
            ElementType sum = ElementType{0};
            // Sum U(i,k) * x(k,col) for k = i+1 to N-1
            for (UINT32 k = static_cast<UINT32>(i) + 1U; k < N; ++k) {
                sum += U(i, k) * x(k, col);
            }
            // Check for division by zero
            if (U(i, i) == ElementType{0}) {
                return MatrixStatus::SINGULAR_MATRIX;
            }
            // x(i,col) = (y(i,col) - sum) / U(i,i)
            x(i, col) = (y(i, col) - sum) / U(i, i);
        }
    }
    
    return MatrixStatus::SUCCESS;
}

/**
 * @brief Backward substitution with runtime size parameter
 * 
 * Solves an upper triangular system U*x = y for x, but only operates on
 * the top-left (size x size) portion of the matrices.
 * 
 * This is useful when working with runtime-determined submatrix sizes.
 * 
 * @param U Upper triangular matrix (must be at least size x size)
 * @param y Right-hand side matrix (must be at least size x YCols)
 * @param size Size of the submatrix to use (1 to min(URows, UCols, YRows))
 * @param x Output parameter for solution (must be at least size x YCols)
 * @return MatrixStatus::SUCCESS if successful, MatrixStatus::OUT_OF_BOUNDS if size is invalid
 */
template <int URows, int UCols, int YRows, int YCols, class ElementType>
MatrixStatus backwardSubstitution(
    const MatrixCls<URows, UCols, ElementType>& U,
    const MatrixCls<YRows, YCols, ElementType>& y,
    int size,
    MatrixCls<YRows, YCols, ElementType>& x) noexcept {
    // Validate size
    if (size < 1 || size > URows || size > UCols || size > YRows) {
        return MatrixStatus::OUT_OF_BOUNDS;
    }
    
    // Initialize x with y (only top-left size×YCols portion)
    x.zero();
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < YCols; ++j) {
            x(i, j) = y(i, j);
        }
    }
    
    // Backward substitution: solve U*x = y (only top-left size×size portion)
    // For each column in y
    for (int col = 0; col < YCols; ++col) {
        // For each row (from bottom to top, only top-left size×size portion)
        for (int i = size - 1; i >= 0; --i) {
            ElementType sum = ElementType{0};
            // Sum U(i,k) * x(k,col) for k = i+1 to size-1
            for (int k = i + 1; k < size; ++k) {
                sum += U(i, k) * x(k, col);
            }
            // Check for division by zero
            if (U(i, i) == ElementType{0}) {
                return MatrixStatus::SINGULAR_MATRIX;
            }
            // x(i,col) = (y(i,col) - sum) / U(i,i)
            x(i, col) = (y(i, col) - sum) / U(i, i);
        }
    }
    
    return MatrixStatus::SUCCESS;
}

/**
 * @brief Backslash operator: Solve A*x = B for x using LU decomposition
 * 
 * This is equivalent to MATLAB's A\B operation.
 * Solves the linear system A*x = B for x, where:
 * - A is a square matrix (Size x Size)
 * - B is a matrix (Size x Cols)
 * - x is the solution (Size x Cols)
 * 
 * Algorithm:
 * 1. Decompose A into L and U: A = L*U
 * 2. Solve L*y = B for y (forward substitution)
 * 3. Solve U*x = y for x (backward substitution)
 * 
 * @param A Coefficient matrix (must be square)
 * @param B Right-hand side matrix
 * @param x Output parameter for solution
 * @return MatrixStatus::SUCCESS if successful, MatrixStatus::SINGULAR_MATRIX if A is singular
 */
template <int Size, int Cols, class ElementType>
MatrixStatus backslash(
    const MatrixCls<Size, Size, ElementType>& A,
    const MatrixCls<Size, Cols, ElementType>& B,
    MatrixCls<Size, Cols, ElementType>& x) noexcept {
    // Step 1: LU decomposition
    MatrixCls<Size, Size, ElementType> L, U;
    MatrixStatus status = A.LuDecomp(L, U);
    if (status != MatrixStatus::SUCCESS) {
        return status;
    }
    
    // Step 2: Forward substitution: solve L*y = B
    MatrixCls<Size, Cols, ElementType> y;
    status = forwardSubstitution(L, B, y);
    if (status != MatrixStatus::SUCCESS) {
        return status;
    }
    
    // Step 3: Backward substitution: solve U*x = y
    status = backwardSubstitution(U, y, x);
    if (status != MatrixStatus::SUCCESS) {
        return status;
    }
    
    return MatrixStatus::SUCCESS;
}

/**
 * @brief Backslash operator overload: A \ B
 * 
 * Solves A*x = B for x using LU decomposition.
 * This is equivalent to MATLAB's A\B operation.
 * 
 * @param A Coefficient matrix (must be square)
 * @param B Right-hand side matrix
 * @return Solution matrix x such that A*x = B
 */
template <int Size, int Cols, class ElementType>
MatrixCls<Size, Cols, ElementType> operator/(
    const MatrixCls<Size, Cols, ElementType>& B,
    const MatrixCls<Size, Size, ElementType>& A) noexcept {
    MatrixCls<Size, Cols, ElementType> x;
    backslash(A, B, x);
    return x;
}

/**
 * @brief Backslash operation with runtime size parameter
 * 
 * Solves A*x = B for x using LU decomposition, but only operates on
 * the top-left (size x size) portion of the matrices.
 * 
 * This is useful when working with runtime-determined submatrix sizes.
 * The matrices A and B can be larger (e.g., 6x6), but only the top-left
 * (size x size) portion will be used for the computation.
 * 
 * Example usage:
 * @code
 * MatrixCls<6, 6, float> mat, mat2;
 * int x = 3;  // Runtime value
 * MatrixCls<6, 6, float> result;
 * 
 * // Extract top-left x×x submatrices
 * MatrixCls<6, 6, float> submat, submat2;
 * mat.extractTopLeft(x, submat);
 * mat2.extractTopLeft(x, submat2);
 * 
 * // Solve submat * x = submat2 (only using top-left x×x portion)
 * backslash(submat, submat2, x, result);
 * // Result is in result(0:x-1, 0:x-1)
 * @endcode
 * 
 * @param A Coefficient matrix (must be at least size x size)
 * @param B Right-hand side matrix (must be at least size x BCols)
 * @param size Size of the submatrix to use (1 to min(ARows, ACols, BRows))
 * @param x Output parameter for solution (must be at least size x BCols)
 * @return MatrixStatus::SUCCESS if successful, MatrixStatus::OUT_OF_BOUNDS if size is invalid
 */
template <int ARows, int ACols, int BRows, int BCols, class ElementType>
MatrixStatus backslash(
    const MatrixCls<ARows, ACols, ElementType>& A,
    const MatrixCls<BRows, BCols, ElementType>& B,
    int size,
    MatrixCls<ARows, BCols, ElementType>& x) noexcept {
    // Validate size
    if (size < 1 || size > ARows || size > ACols || size > BRows) {
        return MatrixStatus::OUT_OF_BOUNDS;
    }
    
    // Create temporary matrices for the size×size submatrix operations
    // We'll use a switch-like approach with template specialization, but for simplicity,
    // we'll work directly with the top-left portion
    
    // Extract top-left size×size from A into a working matrix
    MatrixCls<ARows, ACols, ElementType> A_work;
    A.extractTopLeft(size, A_work);
    
    // Extract top-left size×BCols from B
    MatrixCls<BRows, BCols, ElementType> B_work;
    B.extractTopLeft(size, B_work);
    
    // Step 1: LU decomposition of top-left size×size portion
    MatrixCls<ARows, ACols, ElementType> L, U;
    MatrixStatus status = A_work.LuDecomp(size, L, U);
    if (status != MatrixStatus::SUCCESS) {
        return status;
    }
    
    // Step 2: Forward substitution: solve L*y = B (only top-left size×size of L, size×BCols of B)
    MatrixCls<BRows, BCols, ElementType> y;
    status = forwardSubstitution(L, B_work, size, y);
    if (status != MatrixStatus::SUCCESS) {
        return status;
    }
    
    // Step 3: Backward substitution: solve U*x = y (only top-left size×size of U)
    status = backwardSubstitution(U, y, size, x);
    if (status != MatrixStatus::SUCCESS) {
        return status;
    }
    
    return MatrixStatus::SUCCESS;
}

// Type aliases for common matrix sizes
template <class ElementType = double>
using Matrix2x2 = MatrixCls<2, 2, ElementType>;

template <class ElementType = double>
using Matrix3x3 = MatrixCls<3, 3, ElementType>;

template <class ElementType = double>
using Matrix4x4 = MatrixCls<4, 4, ElementType>;

template <class ElementType = double>
using Vector2 = MatrixCls<2, 1, ElementType>;

template <class ElementType = double>
using Vector3 = MatrixCls<3, 1, ElementType>;

template <class ElementType = double>
using Vector4 = MatrixCls<4, 1, ElementType>;

/**
 * @brief Helper function to extract submatrix with runtime size selection
 * 
 * This function helps extract a top-left submatrix when the size is determined at runtime.
 * It uses a switch statement to select the appropriate fixed-size matrix type.
 * 
 * Example usage:
 * @code
 * MatrixCls<6, 6, float> mat;
 * int x = 3;  // Runtime value
 * MatrixCls<6, 6, float> submat;
 * 
 * switch(x) {
 *     case 1: { MatrixCls<1, 1, float> m; mat.extractTopLeft(x, m); submat.setTopLeft(x, m); break; }
 *     case 2: { MatrixCls<2, 2, float> m; mat.extractTopLeft(x, m); submat.setTopLeft(x, m); break; }
 *     // ... etc
 * }
 * @endcode
 * 
 * For a cleaner approach, use extractTopLeft with a fixed-size result matrix:
 * @code
 * MatrixCls<6, 6, float> mat;
 * int x = 3;
 * MatrixCls<6, 6, float> result;  // Use max size
 * mat.extractTopLeft(x, result);  // Extract top-left x×x
 * // Now work with result(0:x-1, 0:x-1) portion
 * @endcode
 */
template <int MaxSize, class ElementType>
struct SubmatrixHelper {
    /**
     * @brief Extract submatrix into a fixed-size result (recommended approach)
     * 
     * Extracts the top-left (size x size) submatrix into a fixed-size result matrix.
     * Only the top-left portion of the result matrix is populated.
     * 
     * @param source Source matrix
     * @param size Size of submatrix to extract (1 to MaxSize)
     * @param result Output matrix (must be at least MaxSize x MaxSize)
     * @return MatrixStatus::SUCCESS if successful
     */
    template<int SourceRows, int SourceCols>
    static MatrixStatus extract(
        const MatrixCls<SourceRows, SourceCols, ElementType>& source,
        int size,
        MatrixCls<MaxSize, MaxSize, ElementType>& result) noexcept {
        return source.extractTopLeft(size, result);
    }
};

} // namespace math
} // namespace avionics

#endif // MATRIX_LIB_MATRIX_HPP
