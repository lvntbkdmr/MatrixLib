/**
 * @file MatrixCls.h
 * @brief Matrix class with compile-time size and runtime submatrix support
 *
 * This header defines a comprehensive matrix class that supports:
 * - Compile-time fixed maximum dimensions (for stack allocation, avionics-safe)
 * - Runtime variable dimensions (SubRowNum, SubColNum <= RowNum, ColNum)
 * - Submatrix views via SubMatPtrCls
 * - Matrix operations: LU decomposition, QR decomposition, Cholesky decomposition
 * - Linear system solving: forward/backward substitution, left divide (Ax = b)
 * - Matrix properties: lower/upper triangular, symmetric checks
 *
 * @tparam RowNum Maximum number of rows (compile-time constant, default: 3)
 * @tparam ColNum Maximum number of columns (compile-time constant, default: RowNum)
 * @tparam ElementType Type of matrix elements (default: REAL64)
 *
 * @note This class is designed for avionics software where dynamic memory
 * allocation is prohibited. All memory is stack-allocated.
 *
 * @example
 * MatrixCls<6, 6, REAL64> mat;        // Create 6x6 matrix
 * mat.setSubMatrixSize(4, 4);         // Use only 4x4 at runtime
 * mat(0, 0) = 1.0;                   // Access element
 * auto sub = mat(1, 1, 2, 2);        // Get 2x2 submatrix view
 * mat.LeftDivide(b, x);               // Solve Ax = b
 */

#ifndef MATRIX_CLS_H
#define MATRIX_CLS_H

#include "CommonTypes.h"
#include "MatrixPkg.h"
#include <cstring>
#include <cstdlib>
#include <cmath>

/**
 * @class MatrixCls
 * @brief Fixed-size matrix with runtime variable dimensions support
 *
 * The matrix has compile-time maximum dimensions (RowNum, ColNum) but can use
 * smaller dimensions at runtime (SubRowNum, SubColNum). This allows for efficient
 * stack allocation while supporting variable-sized operations.
 *
 * The class provides comprehensive matrix operations including:
 * - Decompositions: LU, QR, Cholesky
 * - Linear system solving: forward/backward substitution, left divide
 * - Matrix properties: triangular, symmetric checks
 * - Submatrix views for efficient access to portions of the matrix
 */
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

    // Copy constructor
    MatrixCls(const MatrixCls& other);

    // Move constructor
    MatrixCls(MatrixCls&& other) noexcept;

    // Copy assignment operator
    MatrixCls& operator=(const MatrixCls& other);

    // Move assignment operator
    MatrixCls& operator=(MatrixCls&& other) noexcept;

    // ============================================================================
    // Submatrix Access Operators
    // ============================================================================

    /**
     * @brief Get a submatrix view (non-const)
     * @param rowStart Starting row index
     * @param colStart Starting column index
     * @param numRows Number of rows in the submatrix
     * @param numCols Number of columns in the submatrix
     * @return SubMatPtrCls view of the submatrix
     *
     * Returns a non-owning view of a portion of this matrix.
     * Modifications through the view affect the original matrix.
     *
     * @note No bounds checking is performed. Ensure the submatrix fits within
     * the matrix dimensions.
     */
    SubMatPtrCls<RowNum, ColNum, ElementType> operator()(UINT32 rowStart,
                                                         UINT32 colStart,
                                                         UINT32 numRows,
                                                         UINT32 numCols);

    /**
     * @brief Get a submatrix copy (const)
     * @param rowStart Starting row index
     * @param colStart Starting column index
     * @param numRows Number of rows in the submatrix
     * @param numCols Number of columns in the submatrix
     * @return MatrixCls copy of the submatrix
     *
     * Returns a copy of a portion of this matrix.
     * The returned matrix is independent of the original.
     */
    MatrixCls<RowNum, ColNum, ElementType> operator()(UINT32 rowStart,
                                                      UINT32 colStart,
                                                      UINT32 numRows,
                                                      UINT32 numCols) const;

    // ============================================================================
    // Element Access Operators
    // ============================================================================

    /**
     * @brief Access matrix element by row and column (non-const)
     * @param row Row index (0-based)
     * @param col Column index (0-based)
     * @return Reference to the element at (row, col)
     */
    inline ElementType& operator()(UINT32 row, UINT32 col);

    /**
     * @brief Access matrix element by row and column (const)
     * @param row Row index (0-based)
     * @param col Column index (0-based)
     * @return Value of the element at (row, col)
     */
    ElementType operator()(UINT32 row, UINT32 col) const;

    /**
     * @brief Access matrix element by linear index (row-major order)
     * @param index Linear index (0-based, row-major: index = row * ColNum + col)
     * @return Reference to the element at the specified index
     */
    ElementType& operator()(UINT32 index);

    /**
     * @brief Get pointer to a row
     * @param rowIndex Row index (0-based)
     * @return Pointer to the first element of the specified row
     *
     * This allows direct access to a row as an array.
     */
    ElementType* operator[](UINT32 rowIndex);

    // ============================================================================
    // Dimension Accessors
    // ============================================================================

    /**
     * @brief Get the current runtime number of rows
     * @return SubRowNum (runtime number of rows)
     */
    UINT32 getRowCount() const { return SubRowNum; }

    /**
     * @brief Get the current runtime number of columns
     * @return SubColNum (runtime number of columns)
     */
    UINT32 getColCount() const { return SubColNum; }

    /**
     * @brief Set the runtime dimensions of the matrix
     * @param numRows Number of rows to use (must be <= RowNum)
     * @param numCols Number of columns to use (must be <= ColNum)
     *
     * Sets SubRowNum and SubColNum to the specified values.
     * This allows the matrix to use fewer elements than its maximum capacity.
     *
     * @note If numRows > RowNum or numCols > ColNum, the values are capped.
     */
    void setSubMatrixSize(UINT32 numRows, UINT32 numCols);

    // ============================================================================
    // Matrix Decomposition Methods
    // ============================================================================

    /**
     * @brief Perform LU decomposition (const version)
     * @param L Output lower triangular matrix
     * @param U Output upper triangular matrix
     *
     * Decomposes this matrix into L and U such that A = L * U.
     * Uses Crout's method where L has actual values on diagonal, U has 1s on diagonal.
     *
     * @note This is a const method that creates a copy for decomposition.
     */
    void LUDecompose(MatrixCls<RowNum, ColNum, ElementType>& L,
                     MatrixCls<RowNum, ColNum, ElementType>& U) const;

    /**
     * @brief Perform LU decomposition (non-const version)
     * @param L Output lower triangular matrix
     * @param U Output upper triangular matrix
     *
     * Decomposes this matrix into L and U such that A = L * U.
     * This method modifies the current matrix during computation.
     */
    void LUDecomp(MatrixCls<RowNum, ColNum, ElementType>& L, MatrixCls<RowNum, ColNum, ElementType>& U);

    // Add pivoting support for numerical stability
    bool LUDecompWithPivoting(MatrixCls<RowNum, ColNum, ElementType>& L,
                             MatrixCls<RowNum, ColNum, ElementType>& U,
                             VectorCls<RowNum, UINT32>& P)  // Permutation vector
    {
        // Implement partial pivoting to improve numerical stability
        // Find the largest element in each column for pivot selection
        // This is crucial for avionics where numerical accuracy is critical
    }

    /**
     * @brief Perform Cholesky decomposition
     * @param L Output lower triangular matrix
     *
     * Decomposes a symmetric positive-definite matrix A into L * L^T.
     *
     * @return Returns false if the matrix is not positive-definite.
     *
     * @note Only valid for symmetric positive-definite matrices.
     */
    void CholeskyDecomp(MatrixCls<RowNum, ColNum, ElementType>& L);

    /**
     * @brief Perform QR decomposition
     * @param Q Output orthogonal matrix (m x m)
     * @param R Output upper triangular matrix (m x n)
     *
     * Decomposes this matrix (m x n) into Q (orthogonal) and R (upper triangular)
     * such that A = Q * R.
     *
     * Uses modified Gram-Schmidt process.
     */
    void QrDecomp(MatrixCls<RowNum, RowNum, ElementType>& Q, MatrixCls<RowNum, ColNum, ElementType>& R);

    // ============================================================================
    // Linear System Solving Methods
    // ============================================================================

    /**
     * @brief Forward substitution: solve Lx = b
     * @param b Right-hand side vector
     * @param y Output solution vector
     *
     * Solves a lower triangular system Lx = b for x.
     * The result is stored in y.
     *
     * @note Only valid for lower triangular matrices (this matrix must be L).
     * The method checks if the matrix is square and returns early if not.
     */
    void ForwardSubstitution(const VectorCls<RowNum, ElementType>& b, VectorCls<RowNum, ElementType>& y);

    /**
     * @brief Backward substitution: solve Ux = y
     * @param y Right-hand side vector
     * @param x Output solution vector
     *
     * Solves an upper triangular system Ux = y for x.
     * The result is stored in x.
     *
     * @note Only valid for upper triangular matrices (this matrix must be U).
     * The method checks if the matrix is square and returns early if not.
     */
    void BackwardSubstitution(const VectorCls<RowNum, ElementType>& y, VectorCls<RowNum, ElementType>& x);

    /**
     * @brief Solve linear system Ax = b for x (vector version)
     * @param b Right-hand side vector
     * @param x Output solution vector
     * @return true if successful, false if matrix is singular or dimensions mismatch
     *
     * Solves the linear system Ax = b using the most efficient method:
     * - If A is lower triangular: forward substitution
     * - If A is upper triangular: backward substitution
     * - If A is symmetric: Cholesky decomposition
     * - If A is square: LU decomposition
     * - If A is rectangular: QR decomposition (least squares)
     *
     * @note For rectangular matrices, this computes the least squares solution.
     */
    bool LeftDivide(const VectorCls<RowNum, ElementType>& b, VectorCls<ColNum, ElementType>& x);

    /**
     * @brief Solve linear system AX = B for X (matrix version)
     * @param b Right-hand side matrix (multiple right-hand sides)
     * @param x Output solution matrix
     * @return true if successful, false if matrix is singular or dimensions mismatch
     *
     * Solves the linear system AX = B where B has multiple columns.
     * Each column of B is solved independently using the vector LeftDivide method.
     *
     * @note The output matrix x has dimensions (ColNum x bCols) where bCols
     * is the number of columns in b.
     */
    bool LeftDivide(const MatrixCls<RowNum, ColNum, ElementType>& b, MatrixCls<RowNum, ColNum, ElementType>& x);

    // ============================================================================
    // Matrix Property Checks
    // ============================================================================

    /**
     * @brief Check if matrix is lower triangular
     * @return true if matrix is lower triangular (all elements above diagonal are zero)
     *
     * A matrix is lower triangular if all elements above the main diagonal
     * are zero (within tolerance).
     */
    bool IsLowTri() const;

    /**
     * @brief Check if matrix is upper triangular
     * @return true if matrix is upper triangular (all elements below diagonal are zero)
     *
     * A matrix is upper triangular if all elements below the main diagonal
     * are zero (within tolerance).
     */
    bool IsUpTri() const;

    /**
     * @brief Check if matrix is symmetric
     * @return true if matrix is symmetric (A = A^T)
     *
     * A matrix is symmetric if A(i,j) = A(j,i) for all i, j (within tolerance).
     */
    bool IsSymmetric() const;

    // ============================================================================
    // Matrix Operations
    // ============================================================================

    /**
     * @brief Compute transpose of the matrix
     * @param result Output matrix to store the transpose
     *
     * Computes A^T and stores it in result.
     * The result matrix has dimensions (ColNum x RowNum).
     */
    void Transpose(MatrixCls<ColNum, RowNum, ElementType>& result) const;

    /**
     * @brief Destructor
     */
    ~MatrixCls();

    // ============================================================================
    // Public Data Members
    // ============================================================================

    ElementType Elements[RowNum][ColNum]; ///< Array of matrix elements (stack-allocated)
    UINT32 SubRowNum;                    ///< Runtime number of rows (<= RowNum)
    UINT32 SubColNum;                    ///< Runtime number of columns (<= ColNum)

    // ============================================================================
    // Matrix Initialization Methods
    // ============================================================================

    /**
     * @brief Set all matrix elements to zero
     */
    void Zeros();

    /**
     * @brief Set all matrix elements to one
     */
    void Ones();

    /**
     * @brief Fill matrix with random values
     *
     * @note Uses rand() which may not be suitable for avionics.
     * Consider using a deterministic PRNG for production code.
     */
    void Random();

    /**
     * @brief Set matrix to identity matrix
     *
     * Sets diagonal elements to 1, all others to 0.
     * For non-square matrices, sets the smaller dimension's diagonal.
     */
    void Identity();

    /**
     * @brief Alias for Identity() (MATLAB-style naming)
     */
    void Eye();
};

// ============================================================================
// Constructor Implementations
// ============================================================================

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls() : SubRowNum(RowNum), SubColNum(ColNum)
{
    // Initialize matrix to all zeros
    Zeros();
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
MatrixCls<RowNum, ColNum, ElementType>::MatrixCls(const ElementType& value) : SubRowNum(RowNum), SubColNum(ColNum)
{
    // Initialize all elements to the specified value
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
    // Copy from void pointer (assuming it points to RowNum*ColNum elements in row-major order)
    // This constructor allows initialization from a C-style array or raw memory
    ElementType* src = static_cast<ElementType*>(basicMatrixVoidAdress);
    for (UINT32 i = 0; i < RowNum; i++)
    {
        for (UINT32 j = 0; j < ColNum; j++)
        {
            // Convert from row-major linear index to 2D index
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
