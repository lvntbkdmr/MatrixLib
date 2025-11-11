#ifndef VECTOR_CLS_H
#define VECTOR_CLS_H

#include "CommonTypes.h"
#include "SubMatPtrCls.h"
#include "MatrixCls.h"

template <UINT32 Length, class ElementType>
class VectorCls
{
public:
    VectorCls();
    VectorCls(const ElementType& value);
    ~VectorCls();

    SubMatPtrCls<Length, 1, ElementType> operator()(UINT32 start, UINT32 length);
    
    ElementType& operator()(UINT32 index);
    ElementType operator()(UINT32 index) const;

    void setLength(UINT32 length);
    
    UINT32 getLength() const { return SubLength; }
    
    // Construct bidiagonal difference matrix
    void BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const;

    UINT32 SubLength;

    ElementType Elements[Length];
};

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::VectorCls() : SubLength(Length)
{
    for(UINT32 i = 0; i < SubLength; i++)
    {
        Elements[i] = ElementType(0);
    }
}

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::VectorCls(const ElementType& value) : SubLength(Length)
{
    for(UINT32 i = 0; i < SubLength; i++)
    {
        Elements[i] = value;
    }
}

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::~VectorCls()
{
}

template <UINT32 Length, class ElementType>
SubMatPtrCls<Length, 1, ElementType> VectorCls<Length, ElementType>::operator()(UINT32 start, UINT32 length)
{
    SubMatPtrCls<Length, 1, ElementType> subVecPtr;
    
    // Bounds check
    if (start + length > SubLength || length > Length)
    {
        subVecPtr.SubRowNum = 0;
        subVecPtr.SubColNum = 0;
        return subVecPtr;
    }
    
    for(UINT32 i = 0; i < length; i++)
    {
        subVecPtr.ElementsPtr[i][0] = &(Elements[start + i]);
    }
    
    subVecPtr.SubRowNum = length;
    subVecPtr.SubColNum = 1;
    
    return subVecPtr;
}

template <UINT32 Length, class ElementType>
ElementType& VectorCls<Length, ElementType>::operator()(UINT32 index)
{
    // Bounds check: ensure index is within runtime size
    if (index >= SubLength)
    {
        // Return first element as fallback (or could throw/assert)
        return Elements[0];
    }
    return Elements[index];
}

template <UINT32 Length, class ElementType>
ElementType VectorCls<Length, ElementType>::operator()(UINT32 index) const
{
    // Bounds check: ensure index is within runtime size
    if (index >= SubLength)
    {
        // Return zero as fallback (or could throw/assert)
        return ElementType(0);
    }
    return Elements[index];
}

template <UINT32 Length, class ElementType>
void VectorCls<Length, ElementType>::setLength(UINT32 length)
{
    // Bounds check: ensure length doesn't exceed compile-time capacity
    if (length > Length)
    {
        SubLength = Length; // Cap at maximum capacity
    }
    else
    {
        SubLength = length;
    }
}

template <UINT32 Length, class ElementType>
void VectorCls<Length, ElementType>::BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const
{
    // Set matrix size: (SubLength - 1) rows x SubLength columns
    // Main diagonal uses rows 0 to SubLength-2, superdiagonal uses columns 0 to SubLength-1
    result.setSubMatrixSize(SubLength - 1, SubLength);
    
    // Initialize to zero
    result.Zeros();
    
    // Fill main diagonal with vector elements (up to SubLength - 1)
    for (UINT32 i = 0; i < SubLength - 1; i++)
    {
        result(i, i) = Elements[i];
    }
    
    // Fill first superdiagonal with -1 * Elements[i+1] at position (i, i+1)
    for (UINT32 i = 0; i < SubLength - 1; i++)
    {
        result(i, i + 1) = ElementType(-1) * Elements[i + 1];
    }
}

#endif