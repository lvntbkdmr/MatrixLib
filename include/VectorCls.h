#pragma once

#include <cstddef>
#include <cmath>
#include <stdexcept>

// Forward declaration
template <typename ElementType>
class SubVecPtrCls;

template <size_t Length, typename ElementType>
class VectorCls
{
public:
    // Constructors
    VectorCls();
    VectorCls(const ElementType& initialValue);
    VectorCls(const ElementType data[Length]);
    
    // Copy constructor
    VectorCls(const VectorCls& other);
    
    // Assignment operator
    VectorCls& operator=(const VectorCls& other);
    
    // Arithmetic operators
    VectorCls operator+(const VectorCls& other) const;
    VectorCls operator-(const VectorCls& other) const;
    VectorCls operator*(const ElementType& scalar) const;
    VectorCls operator/(const ElementType& scalar) const;
    
    // Compound assignment operators
    VectorCls& operator+=(const VectorCls& other);
    VectorCls& operator-=(const VectorCls& other);
    VectorCls& operator*=(const ElementType& scalar);
    VectorCls& operator/=(const ElementType& scalar);
    
    // Dot product
    ElementType dot(const VectorCls& other) const;
    
    // Element access
    ElementType& operator[](size_t index);
    const ElementType& operator[](size_t index) const;
    
    // Access to elements array
    ElementType Elements[Length];
    
    // Utility functions
    size_t size() const { return Length; }
    
    // Runtime-dimension subvector view
    // Use this when you need variable-sized subvectors at runtime
    SubVecPtrCls<ElementType> subvector(size_t start, size_t numElements) const;
    
    // Norm
    ElementType norm() const;
    ElementType normSquared() const;
};

// Scalar multiplication (left side)
template <size_t Length, typename ElementType>
VectorCls<Length, ElementType> operator*(const ElementType& scalar, const VectorCls<Length, ElementType>& vector);

// Template implementations
template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>::VectorCls()
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] = ElementType(0);
    }
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>::VectorCls(const ElementType& initialValue)
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] = initialValue;
    }
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>::VectorCls(const ElementType data[Length])
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] = data[i];
    }
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>::VectorCls(const VectorCls& other)
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] = other.Elements[i];
    }
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>& VectorCls<Length, ElementType>::operator=(const VectorCls& other)
{
    if (this != &other) {
        for (size_t i = 0; i < Length; ++i) {
            Elements[i] = other.Elements[i];
        }
    }
    return *this;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType> VectorCls<Length, ElementType>::operator+(const VectorCls& other) const
{
    VectorCls result;
    for (size_t i = 0; i < Length; ++i) {
        result.Elements[i] = Elements[i] + other.Elements[i];
    }
    return result;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType> VectorCls<Length, ElementType>::operator-(const VectorCls& other) const
{
    VectorCls result;
    for (size_t i = 0; i < Length; ++i) {
        result.Elements[i] = Elements[i] - other.Elements[i];
    }
    return result;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType> VectorCls<Length, ElementType>::operator*(const ElementType& scalar) const
{
    VectorCls result;
    for (size_t i = 0; i < Length; ++i) {
        result.Elements[i] = Elements[i] * scalar;
    }
    return result;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType> VectorCls<Length, ElementType>::operator/(const ElementType& scalar) const
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    VectorCls result;
    for (size_t i = 0; i < Length; ++i) {
        result.Elements[i] = Elements[i] / scalar;
    }
    return result;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>& VectorCls<Length, ElementType>::operator+=(const VectorCls& other)
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] += other.Elements[i];
    }
    return *this;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>& VectorCls<Length, ElementType>::operator-=(const VectorCls& other)
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] -= other.Elements[i];
    }
    return *this;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>& VectorCls<Length, ElementType>::operator*=(const ElementType& scalar)
{
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] *= scalar;
    }
    return *this;
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType>& VectorCls<Length, ElementType>::operator/=(const ElementType& scalar)
{
    if (scalar == ElementType(0)) {
        throw std::runtime_error("Division by zero");
    }
    for (size_t i = 0; i < Length; ++i) {
        Elements[i] /= scalar;
    }
    return *this;
}

template <size_t Length, typename ElementType>
ElementType VectorCls<Length, ElementType>::dot(const VectorCls& other) const
{
    ElementType result = ElementType(0);
    for (size_t i = 0; i < Length; ++i) {
        result += Elements[i] * other.Elements[i];
    }
    return result;
}

template <size_t Length, typename ElementType>
ElementType& VectorCls<Length, ElementType>::operator[](size_t index)
{
    if (index >= Length) {
        throw std::out_of_range("Vector index out of range");
    }
    return Elements[index];
}

template <size_t Length, typename ElementType>
const ElementType& VectorCls<Length, ElementType>::operator[](size_t index) const
{
    if (index >= Length) {
        throw std::out_of_range("Vector index out of range");
    }
    return Elements[index];
}

template <size_t Length, typename ElementType>
ElementType VectorCls<Length, ElementType>::norm() const
{
    ElementType sum = ElementType(0);
    for (size_t i = 0; i < Length; ++i) {
        sum += Elements[i] * Elements[i];
    }
    return std::sqrt(sum);
}

template <size_t Length, typename ElementType>
ElementType VectorCls<Length, ElementType>::normSquared() const
{
    ElementType sum = ElementType(0);
    for (size_t i = 0; i < Length; ++i) {
        sum += Elements[i] * Elements[i];
    }
    return sum;
}

// Include SubVecPtrCls for subvector implementation
#include "SubVecPtrCls.h"

template <size_t Length, typename ElementType>
SubVecPtrCls<ElementType> VectorCls<Length, ElementType>::subvector(size_t start, size_t numElements) const
{
    // Bounds checking
    if (start + numElements > Length) {
        throw std::out_of_range("Subvector indices out of range");
    }
    if (numElements == 0) {
        throw std::invalid_argument("Subvector size must be greater than zero");
    }
    
    // Get pointer to the start of the subvector
    ElementType* subData = const_cast<ElementType*>(&Elements[start]);
    
    // Create and return the subvector view
    return SubVecPtrCls<ElementType>(subData, numElements);
}

template <size_t Length, typename ElementType>
VectorCls<Length, ElementType> operator*(const ElementType& scalar, const VectorCls<Length, ElementType>& vector)
{
    return vector * scalar;
}
