#pragma once

#include <cstddef>
#include <stdexcept>

// Forward declaration
template <size_t Length, typename ElementType>
class VectorCls;

// Forward declaration for SubMatPtrCls result types
template <typename ElementType>
class SubMatPtrCls;

// Sub vector view class - provides a view into a portion of a vector
// Now supports runtime dimensions
template <typename ElementType>
class SubVecPtrCls
{
public:
    // Constructor - takes pointer to the start of the subvector and runtime size
    SubVecPtrCls(ElementType* data, size_t size);
    
    // Copy constructor
    SubVecPtrCls(const SubVecPtrCls& other);
    
    // Assignment operator
    SubVecPtrCls& operator=(const SubVecPtrCls& other);
    
    // Assignment from VectorCls (template to support any size)
    // Checks at runtime that sizes match
    template <size_t Length>
    SubVecPtrCls& operator=(const VectorCls<Length, ElementType>& other);
    
    // Assignment from SubMatPtrCls forward/backward substitution results
    // These are defined after SubMatPtrCls is included (at end of file)
    SubVecPtrCls& operator=(const typename SubMatPtrCls<ElementType>::ForwardSubstitutionResult& result);
    SubVecPtrCls& operator=(const typename SubMatPtrCls<ElementType>::BackwardSubstitutionResult& result);
    
    // Arithmetic operators (return VectorCls since we can't create new views)
    // These are templates because VectorCls requires compile-time size
    template <size_t ResultSize>
    VectorCls<ResultSize, ElementType> operator+(const SubVecPtrCls& other) const;
    
    template <size_t ResultSize>
    VectorCls<ResultSize, ElementType> operator-(const SubVecPtrCls& other) const;
    
    template <size_t ResultSize>
    VectorCls<ResultSize, ElementType> operator*(const ElementType& scalar) const;
    
    template <size_t ResultSize>
    VectorCls<ResultSize, ElementType> operator/(const ElementType& scalar) const;
    
    // Compound assignment operators
    SubVecPtrCls& operator+=(const SubVecPtrCls& other);
    SubVecPtrCls& operator-=(const SubVecPtrCls& other);
    SubVecPtrCls& operator*=(const ElementType& scalar);
    SubVecPtrCls& operator/=(const ElementType& scalar);
    
    // Dot product
    ElementType dot(const SubVecPtrCls& other) const;
    
    // Element access
    ElementType& operator[](size_t index);
    const ElementType& operator[](size_t index) const;
    
    // Utility functions
    size_t size() const { return size_; }
    
private:
    ElementType* data_;
    size_t size_;
};

// Include VectorCls for implementations
#include "VectorCls.h"

// Template implementations
template <typename ElementType>
SubVecPtrCls<ElementType>::SubVecPtrCls(ElementType* data, size_t size)
    : data_(data), size_(size)
{
    if (data_ == nullptr) {
        throw std::invalid_argument("Data pointer cannot be null");
    }
    if (size == 0) {
        throw std::invalid_argument("Subvector size must be greater than zero");
    }
}

template <typename ElementType>
SubVecPtrCls<ElementType>::SubVecPtrCls(const SubVecPtrCls& other)
    : data_(other.data_), size_(other.size_)
{
}

template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator=(const SubVecPtrCls& other)
{
    if (this != &other) {
        if (size_ != other.size_) {
            throw std::runtime_error("Cannot assign SubVecPtrCls of different sizes");
        }
        for (size_t i = 0; i < size_; ++i) {
            data_[i] = other.data_[i];
        }
    }
    return *this;
}

template <typename ElementType>
template <size_t Length>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator=(const VectorCls<Length, ElementType>& other)
{
    if (size_ != Length) {
        throw std::runtime_error("Cannot assign VectorCls of different size to SubVecPtrCls");
    }
    for (size_t i = 0; i < size_; ++i) {
        data_[i] = other.Elements[i];
    }
    return *this;
}

template <typename ElementType>
template <size_t ResultSize>
VectorCls<ResultSize, ElementType> SubVecPtrCls<ElementType>::operator+(const SubVecPtrCls& other) const
{
    if (size_ != other.size_ || size_ != ResultSize) {
        throw std::runtime_error("Vector sizes must match for addition");
    }
    VectorCls<ResultSize, ElementType> result;
    for (size_t i = 0; i < size_; ++i) {
        result.Elements[i] = data_[i] + other.data_[i];
    }
    return result;
}

template <typename ElementType>
template <size_t ResultSize>
VectorCls<ResultSize, ElementType> SubVecPtrCls<ElementType>::operator-(const SubVecPtrCls& other) const
{
    if (size_ != other.size_ || size_ != ResultSize) {
        throw std::runtime_error("Vector sizes must match for subtraction");
    }
    VectorCls<ResultSize, ElementType> result;
    for (size_t i = 0; i < size_; ++i) {
        result.Elements[i] = data_[i] - other.data_[i];
    }
    return result;
}

template <typename ElementType>
template <size_t ResultSize>
VectorCls<ResultSize, ElementType> SubVecPtrCls<ElementType>::operator*(const ElementType& scalar) const
{
    if (size_ != ResultSize) {
        throw std::runtime_error("Vector size must match ResultSize for scalar multiplication");
    }
    VectorCls<ResultSize, ElementType> result;
    for (size_t i = 0; i < size_; ++i) {
        result.Elements[i] = data_[i] * scalar;
    }
    return result;
}

template <typename ElementType>
template <size_t ResultSize>
VectorCls<ResultSize, ElementType> SubVecPtrCls<ElementType>::operator/(const ElementType& scalar) const
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    if (size_ != ResultSize) {
        throw std::runtime_error("Vector size must match ResultSize for scalar division");
    }
    VectorCls<ResultSize, ElementType> result;
    for (size_t i = 0; i < size_; ++i) {
        result.Elements[i] = data_[i] / scalar;
    }
    return result;
}

template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator+=(const SubVecPtrCls& other)
{
    if (size_ != other.size_) {
        throw std::runtime_error("Vector sizes must match for compound addition");
    }
    for (size_t i = 0; i < size_; ++i) {
        data_[i] += other.data_[i];
    }
    return *this;
}

template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator-=(const SubVecPtrCls& other)
{
    if (size_ != other.size_) {
        throw std::runtime_error("Vector sizes must match for compound subtraction");
    }
    for (size_t i = 0; i < size_; ++i) {
        data_[i] -= other.data_[i];
    }
    return *this;
}

template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator*=(const ElementType& scalar)
{
    for (size_t i = 0; i < size_; ++i) {
        data_[i] *= scalar;
    }
    return *this;
}

template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator/=(const ElementType& scalar)
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    for (size_t i = 0; i < size_; ++i) {
        data_[i] /= scalar;
    }
    return *this;
}

template <typename ElementType>
ElementType SubVecPtrCls<ElementType>::dot(const SubVecPtrCls& other) const
{
    if (size_ != other.size_) {
        throw std::runtime_error("Vector sizes must match for dot product");
    }
    ElementType result = ElementType(0);
    for (size_t i = 0; i < size_; ++i) {
        result += data_[i] * other.data_[i];
    }
    return result;
}

template <typename ElementType>
ElementType& SubVecPtrCls<ElementType>::operator[](size_t index)
{
    if (index >= size_) {
        throw std::out_of_range("SubVecPtr index out of range");
    }
    return data_[index];
}

template <typename ElementType>
const ElementType& SubVecPtrCls<ElementType>::operator[](size_t index) const
{
    if (index >= size_) {
        throw std::out_of_range("SubVecPtr index out of range");
    }
    return data_[index];
}

// Include SubMatPtrCls to access result types
#include "SubMatPtrCls.h"

// Assignment operators for SubMatPtrCls substitution results
// These are defined here because they need the full definition of SubMatPtrCls
template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator=(const typename SubMatPtrCls<ElementType>::ForwardSubstitutionResult& result)
{
    result.operator=(*this);
    return *this;
}

template <typename ElementType>
SubVecPtrCls<ElementType>& SubVecPtrCls<ElementType>::operator=(const typename SubMatPtrCls<ElementType>::BackwardSubstitutionResult& result)
{
    result.operator=(*this);
    return *this;
}
