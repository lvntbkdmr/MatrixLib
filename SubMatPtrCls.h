#ifndef SUB_MAT_PTR_CLS_H
#define SUB_MAT_PTR_CLS_H

#include "CommonTypes.h"
#include <cstring>


template <UINT32 RowNum, UINT32 ColNum, class ElementType>
class SubMatPtrCls
{
public:
    SubMatPtrCls();
    ~SubMatPtrCls();

    ElementType* ElementsPtr[RowNum][ColNum];
    UINT32 SubRowNum;
    UINT32 SubColNum;

    ElementType& operator()(UINT32 row, UINT32 col);
    ElementType operator()(UINT32 row, UINT32 col) const;
    
    // Vector-like access: operator()(index) for column vectors (ColNum == 1)
    ElementType& operator()(UINT32 index);
    ElementType operator()(UINT32 index) const;
    
    UINT32 getRowCount() const { return SubRowNum; }
    UINT32 getColCount() const { return SubColNum; }
    
    // For vector compatibility: getLength() returns row count for column vectors
    UINT32 getLength() const { return SubRowNum; }
};

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
SubMatPtrCls<RowNum, ColNum, ElementType>::SubMatPtrCls() : SubRowNum(0), SubColNum(0)
{
    memset(&ElementsPtr, 0, sizeof(ElementType *) * RowNum * ColNum);
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
SubMatPtrCls<RowNum, ColNum, ElementType>::~SubMatPtrCls()
{
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType& SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 row, UINT32 col)
{
    return *(ElementsPtr[row][col]);
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 row, UINT32 col) const
{
    return *(ElementsPtr[row][col]);
}

// Vector-like access: operator()(index) for column vectors (ColNum == 1)
template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType& SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 index)
{
    return (*this)(index, 0);
}

template <UINT32 RowNum, UINT32 ColNum, class ElementType>
ElementType SubMatPtrCls<RowNum, ColNum, ElementType>::operator()(UINT32 index) const
{
    return (*this)(index, 0);
}


#endif