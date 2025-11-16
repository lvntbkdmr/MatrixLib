/**
 * @file test_VectorCls.cpp
 * @brief Unit tests for VectorCls template class
 * 
 * This file contains comprehensive unit tests for the VectorCls class,
 * covering constructors, element access, subvector operations, and
 * matrix construction methods.
 */

#include <gtest/gtest.h>
#include "matrixlib/VectorCls.h"
#include "matrixlib/CommonTypes.h"

// ============================================================================
// Constructor Tests
// ============================================================================

/**
 * @brief Test default constructor
 * Verifies that default constructor initializes all elements to zero
 * and sets SubLength to Length
 */
TEST(VectorClsTest, DefaultConstructor) {
    VectorCls<5, REAL64> vec;
    
    EXPECT_EQ(vec.getLength(), 5);
    
    for (UINT32 i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ(vec(i), 0.0);
    }
}

/**
 * @brief Test constructor with initial value
 * Verifies that constructor initializes all elements to the specified value
 */
TEST(VectorClsTest, ConstructorWithValue) {
    VectorCls<4, REAL64> vec(3.14);
    
    EXPECT_EQ(vec.getLength(), 4);
    
    for (UINT32 i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(vec(i), 3.14);
    }
}

/**
 * @brief Test constructor with different element types
 */
TEST(VectorClsTest, ConstructorWithDifferentTypes) {
    VectorCls<3, REAL32> vecFloat(2.5f);
    VectorCls<3, INT32> vecInt(42);
    
    EXPECT_FLOAT_EQ(vecFloat(0), 2.5f);
    EXPECT_EQ(vecInt(0), 42);
}

// ============================================================================
// Length Management Tests
// ============================================================================

/**
 * @brief Test setLength within valid bounds
 */
TEST(VectorClsTest, SetLengthValid) {
    VectorCls<10, REAL64> vec;
    
    vec.setLength(5);
    EXPECT_EQ(vec.getLength(), 5);
    
    vec.setLength(10);
    EXPECT_EQ(vec.getLength(), 10);
    
    vec.setLength(1);
    EXPECT_EQ(vec.getLength(), 1);
}

/**
 * @brief Test setLength with value exceeding maximum
 * Should cap at compile-time maximum Length
 */
TEST(VectorClsTest, SetLengthExceedsMax) {
    VectorCls<8, REAL64> vec;
    
    vec.setLength(15);
    EXPECT_EQ(vec.getLength(), 8);  // Should be capped at 8
}

/**
 * @brief Test setLength to zero
 */
TEST(VectorClsTest, SetLengthZero) {
    VectorCls<5, REAL64> vec;
    
    vec.setLength(0);
    EXPECT_EQ(vec.getLength(), 0);
}

// ============================================================================
// Element Access Tests
// ============================================================================

/**
 * @brief Test element access (non-const)
 */
TEST(VectorClsTest, ElementAccessNonConst) {
    VectorCls<5, REAL64> vec;
    
    vec(0) = 1.5;
    vec(1) = 2.5;
    vec(2) = 3.5;
    
    EXPECT_DOUBLE_EQ(vec(0), 1.5);
    EXPECT_DOUBLE_EQ(vec(1), 2.5);
    EXPECT_DOUBLE_EQ(vec(2), 3.5);
}

/**
 * @brief Test element access (const)
 */
TEST(VectorClsTest, ElementAccessConst) {
    const VectorCls<3, REAL64> vec(7.0);
    
    EXPECT_DOUBLE_EQ(vec(0), 7.0);
    EXPECT_DOUBLE_EQ(vec(1), 7.0);
    EXPECT_DOUBLE_EQ(vec(2), 7.0);
}

/**
 * @brief Test out-of-bounds access (const version)
 * Const version should return zero for out-of-bounds access
 */
TEST(VectorClsTest, ElementAccessConstOutOfBounds) {
    VectorCls<5, REAL64> vec(3.0);
    vec.setLength(3);  // Only 3 elements are valid
    
    const VectorCls<5, REAL64>& constVec = vec;
    
    EXPECT_DOUBLE_EQ(constVec(0), 3.0);
    EXPECT_DOUBLE_EQ(constVec(2), 3.0);
    EXPECT_DOUBLE_EQ(constVec(3), 0.0);  // Out of bounds, should return 0
    EXPECT_DOUBLE_EQ(constVec(4), 0.0);  // Out of bounds, should return 0
}

/**
 * @brief Test element modification and retrieval
 */
TEST(VectorClsTest, ElementModification) {
    VectorCls<4, REAL64> vec;
    
    for (UINT32 i = 0; i < 4; i++) {
        vec(i) = i * 1.1;
    }
    
    for (UINT32 i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(vec(i), i * 1.1);
    }
}

// ============================================================================
// Subvector Tests
// ============================================================================

/**
 * @brief Test valid subvector creation
 */
TEST(VectorClsTest, SubvectorValid) {
    VectorCls<6, REAL64> vec;
    
    for (UINT32 i = 0; i < 6; i++) {
        vec(i) = i * 10.0;
    }
    
    auto subVec = vec(1, 3);  // Elements 1, 2, 3
    
    EXPECT_EQ(subVec.SubRowNum, 3);
    EXPECT_EQ(subVec.SubColNum, 1);
    
    // Verify the subvector points to the correct elements
    EXPECT_DOUBLE_EQ(*subVec.ElementsPtr[0][0], 10.0);
    EXPECT_DOUBLE_EQ(*subVec.ElementsPtr[1][0], 20.0);
    EXPECT_DOUBLE_EQ(*subVec.ElementsPtr[2][0], 30.0);
}

/**
 * @brief Test that modifying subvector affects original vector
 */
TEST(VectorClsTest, SubvectorModification) {
    VectorCls<5, REAL64> vec(1.0);
    
    auto subVec = vec(1, 2);
    *subVec.ElementsPtr[0][0] = 99.0;
    *subVec.ElementsPtr[1][0] = 88.0;
    
    // Original vector should be modified
    EXPECT_DOUBLE_EQ(vec(1), 99.0);
    EXPECT_DOUBLE_EQ(vec(2), 88.0);
    EXPECT_DOUBLE_EQ(vec(0), 1.0);  // Unaffected
    EXPECT_DOUBLE_EQ(vec(3), 1.0);  // Unaffected
}

/**
 * @brief Test subvector with start = 0
 */
TEST(VectorClsTest, SubvectorFromStart) {
    VectorCls<4, REAL64> vec;
    for (UINT32 i = 0; i < 4; i++) {
        vec(i) = i + 1.0;
    }
    
    auto subVec = vec(0, 2);  // First two elements
    
    EXPECT_EQ(subVec.SubRowNum, 2);
    EXPECT_DOUBLE_EQ(*subVec.ElementsPtr[0][0], 1.0);
    EXPECT_DOUBLE_EQ(*subVec.ElementsPtr[1][0], 2.0);
}

/**
 * @brief Test subvector exceeding bounds
 * Should return empty view (SubRowNum = 0, SubColNum = 0)
 */
TEST(VectorClsTest, SubvectorOutOfBounds) {
    VectorCls<5, REAL64> vec;
    vec.setLength(3);
    
    auto subVec = vec(2, 3);  // start=2, length=3 exceeds SubLength=3
    
    EXPECT_EQ(subVec.SubRowNum, 0);
    EXPECT_EQ(subVec.SubColNum, 0);
}

/**
 * @brief Test subvector with length exceeding max
 */
TEST(VectorClsTest, SubvectorLengthExceedsMax) {
    VectorCls<4, REAL64> vec;
    
    auto subVec = vec(0, 5);  // length=5 exceeds Length=4
    
    EXPECT_EQ(subVec.SubRowNum, 0);
    EXPECT_EQ(subVec.SubColNum, 0);
}

/**
 * @brief Test single element subvector
 */
TEST(VectorClsTest, SubvectorSingleElement) {
    VectorCls<5, REAL64> vec;
    vec(2) = 42.0;
    
    auto subVec = vec(2, 1);
    
    EXPECT_EQ(subVec.SubRowNum, 1);
    EXPECT_EQ(subVec.SubColNum, 1);
    EXPECT_DOUBLE_EQ(*subVec.ElementsPtr[0][0], 42.0);
}

// ============================================================================
// BidiagonalDifference Tests
// ============================================================================

/**
 * @brief Test BidiagonalDifference with simple vector
 */
TEST(VectorClsTest, BidiagonalDifferenceBasic) {
    VectorCls<4, REAL64> vec;
    vec.setLength(3);
    vec(0) = 2.0;
    vec(1) = 3.0;
    vec(2) = 4.0;
    
    MatrixCls<4, 4, REAL64> mat;
    vec.BidiagonalDifference(mat);
    
    // Result should be (SubLength-1) x SubLength = 2x3
    EXPECT_EQ(mat.SubRowNum, 2);
    EXPECT_EQ(mat.SubColNum, 3);
    
    // First row: [2.0, -3.0, 0.0]
    EXPECT_DOUBLE_EQ(mat(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(mat(0, 1), -3.0);
    EXPECT_DOUBLE_EQ(mat(0, 2), 0.0);
    
    // Second row: [0.0, 3.0, -4.0]
    EXPECT_DOUBLE_EQ(mat(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(mat(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(mat(1, 2), -4.0);
}

/**
 * @brief Test BidiagonalDifference with length 2 vector
 */
TEST(VectorClsTest, BidiagonalDifferenceLength2) {
    VectorCls<3, REAL64> vec;
    vec.setLength(2);
    vec(0) = 5.0;
    vec(1) = 7.0;
    
    MatrixCls<3, 3, REAL64> mat;
    vec.BidiagonalDifference(mat);
    
    // Result should be 1x2
    EXPECT_EQ(mat.SubRowNum, 1);
    EXPECT_EQ(mat.SubColNum, 2);
    
    // Only row: [5.0, -7.0]
    EXPECT_DOUBLE_EQ(mat(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(mat(0, 1), -7.0);
}

/**
 * @brief Test BidiagonalDifference with larger vector
 */
TEST(VectorClsTest, BidiagonalDifferenceLarger) {
    VectorCls<6, REAL64> vec;
    vec.setLength(5);
    for (UINT32 i = 0; i < 5; i++) {
        vec(i) = (i + 1) * 1.0;  // 1, 2, 3, 4, 5
    }
    
    MatrixCls<6, 6, REAL64> mat;
    vec.BidiagonalDifference(mat);
    
    // Result should be 4x5
    EXPECT_EQ(mat.SubRowNum, 4);
    EXPECT_EQ(mat.SubColNum, 5);
    
    // Check diagonal elements
    EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(mat(2, 2), 3.0);
    EXPECT_DOUBLE_EQ(mat(3, 3), 4.0);
    
    // Check superdiagonal elements (negative)
    EXPECT_DOUBLE_EQ(mat(0, 1), -2.0);
    EXPECT_DOUBLE_EQ(mat(1, 2), -3.0);
    EXPECT_DOUBLE_EQ(mat(2, 3), -4.0);
    EXPECT_DOUBLE_EQ(mat(3, 4), -5.0);
    
    // Check some zero elements
    EXPECT_DOUBLE_EQ(mat(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(mat(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(mat(2, 0), 0.0);
}

/**
 * @brief Test BidiagonalDifference with negative values
 */
TEST(VectorClsTest, BidiagonalDifferenceNegativeValues) {
    VectorCls<4, REAL64> vec;
    vec.setLength(3);
    vec(0) = -2.0;
    vec(1) = 3.0;
    vec(2) = -4.0;
    
    MatrixCls<4, 4, REAL64> mat;
    vec.BidiagonalDifference(mat);
    
    // First row: [-2.0, -3.0, 0.0]
    EXPECT_DOUBLE_EQ(mat(0, 0), -2.0);
    EXPECT_DOUBLE_EQ(mat(0, 1), -3.0);
    
    // Second row: [0.0, 3.0, 4.0]  (note: -(-4.0) = 4.0)
    EXPECT_DOUBLE_EQ(mat(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(mat(1, 2), 4.0);
}

// ============================================================================
// Edge Cases and Integration Tests
// ============================================================================

/**
 * @brief Test vector with length 1
 */
TEST(VectorClsTest, VectorLength1) {
    VectorCls<5, REAL64> vec;
    vec.setLength(1);
    vec(0) = 99.0;
    
    EXPECT_EQ(vec.getLength(), 1);
    EXPECT_DOUBLE_EQ(vec(0), 99.0);
}

/**
 * @brief Test multiple operations in sequence
 */
TEST(VectorClsTest, SequentialOperations) {
    VectorCls<10, REAL64> vec;
    
    // Initialize
    vec.setLength(5);
    for (UINT32 i = 0; i < 5; i++) {
        vec(i) = i * 2.0;
    }
    
    // Get subvector and modify
    auto subVec = vec(1, 3);
    *subVec.ElementsPtr[0][0] = 100.0;
    
    // Verify modification
    EXPECT_DOUBLE_EQ(vec(1), 100.0);
    
    // Change length
    vec.setLength(3);
    EXPECT_EQ(vec.getLength(), 3);
    
    // Out of bounds access on const should return 0
    const VectorCls<10, REAL64>& constVec = vec;
    EXPECT_DOUBLE_EQ(constVec(4), 0.0);
}

/**
 * @brief Test with REAL32 type
 */
TEST(VectorClsTest, Real32Type) {
    VectorCls<4, REAL32> vec(1.5f);
    
    vec(0) = 2.5f;
    vec(1) = 3.5f;
    
    EXPECT_FLOAT_EQ(vec(0), 2.5f);
    EXPECT_FLOAT_EQ(vec(1), 3.5f);
    EXPECT_FLOAT_EQ(vec(2), 1.5f);
}

/**
 * @brief Test with INT32 type
 */
TEST(VectorClsTest, Int32Type) {
    VectorCls<5, INT32> vec(10);
    
    vec(0) = 100;
    vec(1) = 200;
    
    EXPECT_EQ(vec(0), 100);
    EXPECT_EQ(vec(1), 200);
    EXPECT_EQ(vec(2), 10);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
