#ifndef VECTOR_CLS_H
#define VECTOR_CLS_H

#include "CommonTypes.h"
#include "SubMatPtrCls.h"

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
    
    UINT32 getLength() const { return Length; }

    ElementType Elements[Length];
};

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::VectorCls()
{
    for(UINT32 i = 0; i < Length; i++)
    {
        Elements[i] = ElementType(0);
    }
}

template <UINT32 Length, class ElementType>
VectorCls<Length, ElementType>::VectorCls(const ElementType& value)
{
    for(UINT32 i = 0; i < Length; i++)
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
    return Elements[index];
}

template <UINT32 Length, class ElementType>
ElementType VectorCls<Length, ElementType>::operator()(UINT32 index) const
{
    return Elements[index];
}

#endif