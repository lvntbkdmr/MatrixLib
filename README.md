# Matrix Library for Avionics Software

A C++ matrix library designed for avionics software applications, conforming to safety-critical software standards including DO-178C and MISRA C++ guidelines.

## Features

- **Template-based design** - Support for various numeric types (float, double, int32_t, etc.)
- **Function overloading** - Intuitive operator overloading for matrix operations
- **Fixed-size matrices** - Compile-time dimensions (no dynamic memory allocation)
- **Avionics compliant** - No exceptions, return codes for error handling
- **MISRA C++ compliant** - Uses restricted C++ subset suitable for safety-critical code
- **Deterministic behavior** - Predictable execution with bounded operations
- **Type safety** - Strong typing with template constraints

## Avionics Compliance Features

1. **No Dynamic Memory Allocation** - All matrices use fixed-size arrays
2. **No Exceptions** - Error handling via return codes (`MatrixStatus`)
3. **Compile-time Size Checking** - Matrix dimensions are template parameters
4. **Deterministic Execution** - All operations have predictable behavior
5. **MISRA C++ Friendly** - Uses restricted C++ subset
6. **Bounds Checking** - Optional bounds checking via `get()` and `set()` methods
7. **Comprehensive Documentation** - Doxygen-style comments throughout

## Requirements

- C++14 or later
- CMake 3.12 or later (for building examples and tests)

## Embedded System Compatibility

The library is designed for embedded systems and **does not require**:
- `printf()` or any standard I/O functions
- `<iostream>` or `<cstdio>` headers
- Standard library I/O streams

The library only requires:
- `<cstddef>` - For `std::size_t` (standard size type)
- `<cstdint>` - For `std::int32_t` (error code type)
- `<type_traits>` - For `std::is_arithmetic` (type checking)

All matrix operations are performed in-memory with no I/O dependencies. The examples and tests use `std::cout` for demonstration purposes only, but the library itself is completely I/O-free and suitable for bare-metal embedded systems.

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Usage Example

```cpp
#include "matrix.hpp"
using namespace avionics::math;

// Create a 3x3 identity matrix
Matrix<float, 3U, 3U> mat1;
mat1.identity();

// Create a matrix with initial values
const float init_data[9] = {1.0f, 2.0f, 3.0f,
                            4.0f, 5.0f, 6.0f,
                            7.0f, 8.0f, 9.0f};
Matrix<float, 3U, 3U> mat2(init_data);

// Matrix operations using function overloading
auto mat3 = mat1 + mat2;        // Addition
auto mat4 = mat2 * 2.0f;         // Scalar multiplication
auto mat5 = transpose(mat2);     // Transpose
auto mat6 = matA * matB;         // Matrix multiplication

// Access elements
float value = mat2(1, 2);  // Gets element at row 1, col 2

// Safe access with bounds checking
float safe_value;
if (mat2.get(1, 2, safe_value) == MatrixStatus::SUCCESS) {
    // Use safe_value
}
```

## Type Aliases

Common matrix sizes are available as type aliases:

```cpp
Matrix2x2<float>  // 2x2 matrix
Matrix3x3<float>  // 3x3 matrix
Matrix4x4<float>  // 4x4 matrix
Vector2<float>    // 2x1 vector
Vector3<float>    // 3x1 vector
Vector4<float>    // 4x1 vector
```

## Supported Operations

- **Arithmetic**: Addition, subtraction, scalar multiplication/division
- **Matrix Operations**: Matrix multiplication, transpose
- **Linear Algebra**: Determinant (2x2, 3x3), inverse (2x2)
- **Utility**: Identity matrix, zero matrix, fill with constant

## Error Handling

All operations that can fail return a `MatrixStatus` enum:

```cpp
enum class MatrixStatus {
    SUCCESS = 0,
    INVALID_DIMENSIONS = -1,
    OUT_OF_BOUNDS = -2,
    SINGULAR_MATRIX = -3,
    INCOMPATIBLE_SIZES = -4
};
```

## Running Examples

```bash
./matrix_example
```

## Running Tests

```bash
./matrix_tests
```

Or with CMake test framework:

```bash
cmake -DBUILD_TESTS=ON ..
make
ctest
```

## Avionics Standards Compliance

This library is designed to support compliance with:

- **DO-178C** - Software Considerations in Airborne Systems and Equipment Certification
- **MISRA C++** - Guidelines for the use of C++ in critical systems
- **ARINC 653** - Avionics Application Software Standard Interface

### Key Design Decisions for Compliance

1. **No Dynamic Allocation**: All memory is allocated at compile-time
2. **No Exceptions**: Error handling via return codes
3. **Bounded Operations**: All operations have compile-time known bounds
4. **Deterministic**: No undefined behavior or non-deterministic operations
5. **Type Safety**: Strong typing with template constraints
6. **Documentation**: Comprehensive comments for verification

## License

[Specify your license here]

## Contributing

[Contributing guidelines if applicable]

## References

- DO-178C: Software Considerations in Airborne Systems and Equipment Certification
- MISRA C++: Guidelines for the use of the C++ language in critical systems
- ARINC 653: Avionics Application Software Standard Interface

