# Forward Substitution with Runtime-Sized Submatrices

This example demonstrates how to use `forwardSubstitution` with runtime-sized submatrices extracted from compile-time sized matrices.

## Key Concepts

### 1. Compile-Time vs Runtime Dimensions

- **MatrixCls**: Has compile-time dimensions (e.g., `MatrixCls<6, 6, double>`)
- **SubMatPtrCls**: Has runtime dimensions (e.g., `SubMatPtrCls<double>`)
- You can extract a variable-sized submatrix from a fixed-size matrix at runtime

### 2. Creating Submatrices

There are two ways to create a submatrix view:

#### Method 1: `submatrix()` function
```cpp
MatrixCls<6, 6, double> L;
size_t size = 3;  // Runtime decision
SubMatPtrCls<double> sub = L.submatrix(0, 0, size, size);  // Top-left 3x3
```

#### Method 2: `operator()()` convenience operator
```cpp
SubMatPtrCls<double> sub = L(0, 0, 2, 2);  // From (0,0) to (2,2) inclusive
```

### 3. Using Forward Substitution

#### With VectorCls
```cpp
SubMatPtrCls<double> L_sub = matrix.submatrix(0, 0, 3, 3);
VectorCls<3, double> b;
VectorCls<3, double> x = L_sub.forwardSubstitution<3>(b);
```

#### With SubVecPtrCls
```cpp
VectorCls<6, double> b_full;
SubVecPtrCls<3, double> b_sub(&b_full.Elements[0]);  // View of first 3 elements
VectorCls<3, double> x = L_sub.forwardSubstitution<3>(b_sub);
```

#### Using Generic Function
```cpp
VectorCls<3, double> x = forwardSubstitution(L_sub, b);
```

## Important Notes

1. **Template Parameter**: The template parameter `<3>` in `forwardSubstitution<3>()` must match the vector size. This is a compile-time requirement.

2. **Runtime Dimension Check**: When using `SubMatPtrCls`, the function checks at runtime that the matrix dimensions match the vector size. If they don't match, it throws a `std::runtime_error`.

3. **Zero-Copy Views**: `SubMatPtrCls` and `SubVecPtrCls` are views - they don't copy data, they just point to the original matrix/vector.

4. **Lower Triangular Requirement**: The matrix must be lower triangular for forward substitution to work correctly.

## Compilation

```bash
clang++ -std=c++17 -I include examples/forward_substitution_example.cpp -o forward_substitution_example
```

## Running

```bash
./forward_substitution_example
```

The example demonstrates:
- Extracting fixed-size submatrices (3x3, 4x4)
- Extracting variable-size submatrices based on runtime decisions
- Using both `VectorCls` and `SubVecPtrCls` as right-hand side vectors
- Using both member functions and generic functions
- Extracting non-top-left submatrices

