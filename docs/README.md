# Matrix Library

Avionics-safe matrix and vector library with compile-time dimensions and runtime variable sizes.

## Features

- **Compile-time dimensions**: Maximum sizes specified at compile time for stack allocation
- **Runtime flexibility**: Use smaller dimensions at runtime than the compile-time maximum
- **No dynamic allocation**: Completely stack-allocated, avionics-safe
- **Comprehensive operations**: Matrix decompositions, linear system solving, and more
- **Submatrix views**: Efficient non-owning views into portions of matrices

## Project Structure

```
MatrixLib/
├── include/
│   └── matrixlib/          # Public header files
│       ├── CommonTypes.h   # Type definitions
│       ├── MatrixCls.h     # Matrix class
│       ├── VectorCls.h     # Vector class
│       ├── SubMatPtrCls.h  # Submatrix view class
│       └── MatrixPkg.h     # Main include file
├── src/                    # Source files (if any)
├── examples/               # Example code
│   └── app.cpp            # Test application
├── tests/                  # Test files
├── docs/                   # Documentation
│   ├── sphinx/            # Sphinx documentation
│   ├── doxygen_xml/       # Doxygen XML (generated)
│   └── html/              # Doxygen HTML (generated)
├── scripts/                # Build and utility scripts
├── build/                  # Build artifacts (generated)
├── CMakeLists.txt         # CMake configuration
└── README.md              # This file
```

## Quick Start

### Include the Library

```cpp
#include "matrixlib/MatrixPkg.h"

// Create a 6x6 matrix (max size), use 4x4 at runtime
MatrixCls<6, 6, REAL64> A;
A.setSubMatrixSize(4, 4);

// Create vectors
VectorCls<6, REAL64> b, x;
b.setLength(4);

// Solve Ax = b
A.LeftDivide(b, x);
```

### Building

```bash
mkdir build
cd build
cmake ..
make
```

### Running Examples

```bash
cd build
./MatrixPkg
```

## Documentation

### Generate Doxygen Documentation

```bash
doxygen Doxyfile
open docs/html/index.html
```

### Generate Sphinx Documentation

```bash
./scripts/generate_sphinx_docs.sh
open docs/sphinx/_build/html/index.html
```

## Requirements

- C++17 or later
- CMake 3.10 or later
- Doxygen (for documentation generation)
- Python 3 + Sphinx + Breathe (for Sphinx documentation)

## License

[Add your license here]

## Contributing

[Add contribution guidelines here]

