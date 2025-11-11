# Project Structure

This document describes the organization of the Matrix Library project.

## Directory Layout

```
MatrixLib/
├── include/                    # Public header files
│   └── matrixlib/             # Library namespace directory
│       ├── CommonTypes.h      # Type definitions (UINT32, REAL64, etc.)
│       ├── MatrixCls.h        # Matrix class implementation
│       ├── VectorCls.h        # Vector class implementation
│       ├── SubMatPtrCls.h     # Submatrix view class
│       └── MatrixPkg.h        # Main include file (includes all headers)
│
├── src/                       # Source files (currently empty - header-only library)
│
├── examples/                   # Example code and demos
│   └── app.cpp                # Test application demonstrating LeftDivide
│
├── tests/                     # Unit tests (to be added)
│
├── docs/                      # Documentation
│   ├── sphinx/                # Sphinx documentation source
│   │   ├── conf.py            # Sphinx configuration
│   │   ├── index.rst          # Main documentation page
│   │   ├── classes.rst        # Class documentation
│   │   ├── files.rst          # File documentation
│   │   └── examples.rst       # Usage examples
│   ├── doxygen_xml/           # Doxygen XML output (generated)
│   ├── html/                  # Doxygen HTML output (generated)
│   └── *.md                   # Documentation guides
│
├── scripts/                    # Build and utility scripts
│   ├── generate_docs.sh       # Generate Doxygen documentation
│   └── generate_sphinx_docs.sh # Generate Sphinx documentation
│
├── build/                     # Build artifacts (generated, gitignored)
│
├── venv/                      # Python virtual environment (generated, gitignored)
│
├── CMakeLists.txt             # CMake build configuration
├── Doxyfile                   # Doxygen configuration
├── requirements.txt           # Python dependencies for Sphinx
├── .gitignore                 # Git ignore rules
├── .clang-format              # Code formatting configuration
└── README.md                   # Project overview and quick start
```

## Key Directories

### `include/matrixlib/`
Contains all public header files. Users include the library via:
```cpp
#include "matrixlib/MatrixPkg.h"
```

This directory structure allows for:
- Clear namespace organization
- Easy installation to system include directories
- Prevention of header name conflicts

### `examples/`
Contains example code demonstrating library usage. The `app.cpp` file shows:
- Matrix and vector creation
- Runtime dimension management
- Linear system solving
- Submatrix operations

### `docs/`
All documentation-related files:
- **Sphinx source**: RST files for Sphinx documentation
- **Doxygen output**: Generated XML and HTML from Doxygen
- **Markdown guides**: README files and documentation guides

### `scripts/`
Utility scripts for:
- Generating documentation (Doxygen and Sphinx)
- Build automation
- Development tasks

## Build System

The project uses CMake for building:

```bash
mkdir build
cd build
cmake ..
make
```

The CMake configuration:
- Sets C++17 standard
- Creates an interface library target (`matrixlib`)
- Builds example executable (`matrixlib_example`)
- Configures include directories

## Usage

### Including the Library

```cpp
#include "matrixlib/MatrixPkg.h"

// Use the library
MatrixCls<6, 6, REAL64> A;
VectorCls<6, REAL64> b, x;
A.LeftDivide(b, x);
```

### Building Examples

```bash
cd build
make
./matrixlib_example
```

## File Organization Principles

1. **Separation of concerns**: Headers in `include/`, examples in `examples/`, docs in `docs/`
2. **Namespace organization**: Headers in `include/matrixlib/` to prevent conflicts
3. **Build artifacts**: All generated files in `build/`, `docs/*/`, `venv/`
4. **Documentation**: All docs in `docs/` with clear subdirectories
5. **Scripts**: All utility scripts in `scripts/` for easy discovery

## Future Additions

- `src/`: For any non-template implementation files (if needed)
- `tests/`: Unit tests using a testing framework
- `benchmarks/`: Performance benchmarks
- `tools/`: Development tools and utilities
