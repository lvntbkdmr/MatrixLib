/**
 * @file VectorCls.h
 * @brief Vector class with compile-time size and runtime subvector support
 *
 * This header defines a vector class that supports:
 * - Compile-time fixed maximum size (for stack allocation, avionics-safe)
 * - Runtime variable length (SubLength <= Length)
 * - Subvector views via SubMatPtrCls
 * - Matrix operations (bidiagonal difference matrix construction)
 *
 * @tparam Length Maximum vector length (compile-time constant)
 * @tparam ElementType Type of vector elements (e.g., REAL64, REAL32)
 *
 * @note This class is designed for avionics software where dynamic memory
 * allocation is prohibited. All memory is stack-allocated.
 *
 * @example
 * VectorCls<6, REAL64> vec;        // Create vector with max size 6
 * vec.setLength(3);                 // Use only first 3 elements at runtime
 * vec(0) = 1.0;                    // Access element 0
 * auto sub = vec(1, 2);            // Get subvector view of elements 1-2
 */

#ifndef VECTOR_CLS_H
#define VECTOR_CLS_H

#include "CommonTypes.h"
#include "MatrixPkg.h"

/**
 * @class VectorCls
 * @brief Fixed-size vector with runtime variable length support
 *
 * The vector has a compile-time maximum size (Length) but can use a smaller
 * number of elements at runtime (SubLength). This allows for efficient
 * stack allocation while supporting variable-sized operations.
 */
template <UINT32 Length, class ElementType = REAL64>
class VectorCls
{
public:
    /**
     * @brief Default constructor - initializes vector to zeros
     *
     * Creates a vector with all elements set to zero and SubLength = Length.
     */
    VectorCls();

    /**
     * @brief Constructor with initial value
     * @param value Value to initialize all elements to
     *
     * Creates a vector with all elements set to the specified value
     * and SubLength = Length.
     */
    VectorCls(const ElementType& value);

    /**
     * @brief Destructor
     */
    ~VectorCls();

    /**
     * @brief Get a subvector view
     * @param start Starting index of the subvector
     * @param length Number of elements in the subvector
     * @return SubMatPtrCls view of the subvector (as a column vector)
     *
     * Returns a non-owning view of a portion of this vector.
     * The view points to the original vector's elements, so modifications
     * through the view affect the original vector.
     *
     * @note If start + length exceeds SubLength or length > Length,
     * returns an empty view (SubRowNum = 0, SubColNum = 0).
     *
     * @example
     * VectorCls<6, REAL64> vec;
     * vec.setLength(5);
     * auto sub = vec(1, 3);  // Get view of elements 1, 2, 3
     */
    SubMatPtrCls<Length, 1, ElementType> operator()(UINT32 start, UINT32 length);

    /**
     * @brief Access element by index (non-const)
     * @param index Element index (0-based)
     * @return Reference to the element at the specified index
     *
     * @note If index >= SubLength, returns reference to first element
     * as a fallback. Consider using bounds checking in production code.
     */
    ElementType& operator()(UINT32 index);

    /**
     * @brief Access element by index (const)
     * @param index Element index (0-based)
     * @return Value of the element at the specified index
     *
     * @note If index >= SubLength, returns zero as a fallback.
     */
    ElementType operator()(UINT32 index) const;

    /**
     * @brief Set the runtime length of the vector
     * @param length New runtime length
     *
     * Sets SubLength to the specified value, but caps it at Length
     * (the compile-time maximum). This allows the vector to use
     * fewer elements than its maximum capacity at runtime.
     *
     * @note If length > Length, SubLength is set to Length.
     *
     * @example
     * VectorCls<6, REAL64> vec;  // Max size 6
     * vec.setLength(3);          // Use only 3 elements
     */
    void setLength(UINT32 length);

    /**
     * @brief Get the current runtime length
     * @return Current SubLength (number of active elements)
     */
    UINT32 getLength() const { return SubLength; }

    void BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const;

    UINT32 SubLength;           ///< Runtime length (number of active elements, <= Length)
    ElementType Elements[Length]; ///< Array of vector elements (stack-allocated)
};

// ============================================================================
// Constructor Implementations
// ============================================================================

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::VectorCls() : SubLength(Length)
{
    // Initialize all elements to zero
    for(UINT32 i = 0; i < SubLength; i++)
    {
        Elements[i] = ElementType(0);
    }
}

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::VectorCls(const ElementType& value) : SubLength(Length)
{
    // Initialize all elements to the specified value
    for(UINT32 i = 0; i < SubLength; i++)
    {
        Elements[i] = value;
    }
}

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::~VectorCls()
{
    // No dynamic memory to clean up (stack-allocated array)
}

// ============================================================================
// Subvector View Implementation
// ============================================================================

template <UINT32 Length, class ElementType>
SubMatPtrCls<Length, 1, ElementType> VectorCls<Length, ElementType>::operator()(UINT32 start, UINT32 length)
{
    SubMatPtrCls<Length, 1, ElementType> subVecPtr;

    // Bounds checking: ensure the subvector is within valid range
    if (start + length > SubLength || length > Length)
    {
        // Return empty view if bounds are invalid
        subVecPtr.SubRowNum = 0;
        subVecPtr.SubColNum = 0;
        return subVecPtr;
    }

    // Set up pointers to the subvector elements
    // The subvector is treated as a column vector (1 column)
    for(UINT32 i = 0; i < length; i++)
    {
        subVecPtr.ElementsPtr[i][0] = &(Elements[start + i]);
    }

    // Set the dimensions of the subvector view
    subVecPtr.SubRowNum = length;
    subVecPtr.SubColNum = 1;

    return subVecPtr;
}

// ============================================================================
// Element Access Implementations
// ============================================================================

template <UINT32 Length, class ElementType>
ElementType& VectorCls<Length, ElementType>::operator()(UINT32 index)
{
    #ifdef DEBUG
    assert(index < SubLength && "Vector index out of bounds");
    #elif defined(AVIONICS_SAFE)
    if (index >= SubLength) {
        // Log error or set error flag instead of silent fallback
        // For now, clamp to last valid element
        index = SubLength > 0 ? SubLength - 1 : 0;
    }
    #endif
    return Elements[index];
}

template <UINT32 Length, class ElementType>
ElementType VectorCls<Length, ElementType>::operator()(UINT32 index) const
{
    // Bounds check: ensure index is within runtime size (SubLength)
    if (index >= SubLength)
    {
        // Return zero as fallback for out-of-bounds access
        // Note: In production code, consider throwing an exception or using assert
        return ElementType(0);
    }
    return Elements[index];
}

// ============================================================================
// Length Management Implementation
// ============================================================================

template <UINT32 Length, class ElementType>
void VectorCls<Length, ElementType>::setLength(UINT32 length)
{
    // Bounds check: ensure length doesn't exceed compile-time capacity
    if (length > Length)
    {
        // Cap at maximum capacity to prevent out-of-bounds access
        SubLength = Length;
    }
    else
    {
        SubLength = length;
    }
}

// ============================================================================
// Matrix Construction Implementation
// ============================================================================

/** \fn template <UINT32 Length, class ElementType> void VectorCls<Length, ElementType>::BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const
 *  \brief Construct a bidiagonal difference matrix from this vector
 *  \param result Output matrix (will be resized to (SubLength-1) x SubLength)
 *
 *  Constructs a matrix where:
 *  - Main diagonal (i, i) contains Elements[i] for i = 0 to SubLength-2
 *  - First superdiagonal (i, i+1) contains -Elements[i+1] for i = 0 to SubLength-2
 *  - All other elements are zero
 *
 *  The resulting matrix has dimensions (SubLength-1) x SubLength.
 *
 *  \example
 *  VectorCls<4, REAL64> vec;
 *  vec.setLength(3);
 *  vec(0) = 2.0; vec(1) = 3.0; vec(2) = 4.0;
 *  MatrixCls<3, 4, REAL64> mat;
 *  vec.BidiagonalDifference(mat);
 *  // Result: [2.0, -3.0,  0.0,  0.0]
 *  //         [0.0,  3.0, -4.0,  0.0]
 */
template <UINT32 Length, class ElementType>
void VectorCls<Length, ElementType>::BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const
{
    // Set matrix size: (SubLength - 1) rows x SubLength columns
    // This creates a matrix with one fewer row than the vector length
    result.setSubMatrixSize(SubLength - 1, SubLength);

    // Initialize entire matrix to zero
    result.Zeros();

    // Fill main diagonal with vector elements (rows 0 to SubLength-2)
    // Diagonal element at (i, i) gets Elements[i]
    for (UINT32 i = 0; i < SubLength - 1; i++)
    {
        result(i, i) = Elements[i];
    }

    // Fill first superdiagonal with negative vector elements
    // Element at (i, i+1) gets -Elements[i+1]
    // This creates the "difference" part of the bidiagonal difference matrix
    for (UINT32 i = 0; i < SubLength - 1; i++)
    {
        result(i, i + 1) = ElementType(-1) * Elements[i + 1];
    }
}

#endif
