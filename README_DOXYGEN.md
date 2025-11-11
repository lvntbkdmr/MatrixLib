# Generating Doxygen Documentation

This guide explains how to generate HTML documentation from the source code comments using Doxygen.

## Prerequisites

### Installing Doxygen

**On macOS (using Homebrew):**
```bash
brew install doxygen
```

**On Ubuntu/Debian:**
```bash
sudo apt-get install doxygen
```

**On Windows:**
Download and install from: https://www.doxygen.nl/download.html

**Verify installation:**
```bash
doxygen --version
```

## Generating Documentation

### Quick Start

1. **Generate documentation:**
   ```bash
   doxygen Doxyfile
   ```

2. **View the documentation:**
   Open `docs/html/index.html` in your web browser.

### What Gets Generated

The documentation will be generated in the `docs/html/` directory and includes:

- **Main Page**: Overview of the Matrix Library
- **Classes**: Documentation for `MatrixCls`, `VectorCls`, `SubMatPtrCls`
- **Files**: Documentation for each header file
- **Class Hierarchy**: Visual representation of class relationships
- **File List**: List of all documented files
- **Search**: Full-text search functionality

### Configuration

The `Doxyfile` is pre-configured with the following settings:

- **Output Format**: HTML (with tree view)
- **Extract All**: Yes (documents all public members)
- **Source Browser**: Enabled (links to source code)
- **Search Engine**: Enabled
- **Class Diagrams**: Enabled (if Graphviz is installed)

### Advanced Options

**To enable class diagrams (requires Graphviz):**

1. Install Graphviz:
   ```bash
   # macOS
   brew install graphviz
   
   # Ubuntu/Debian
   sudo apt-get install graphviz
   ```

2. Edit `Doxyfile` and change:
   ```
   HAVE_DOT = YES
   ```

3. Regenerate documentation:
   ```bash
   doxygen Doxyfile
   ```

**To generate PDF documentation:**

1. Install LaTeX (required for PDF generation)

2. Edit `Doxyfile` and change:
   ```
   GENERATE_LATEX = YES
   ```

3. Generate documentation:
   ```bash
   doxygen Doxyfile
   cd docs/latex
   make
   ```

## Documentation Structure

The generated documentation includes:

### Main Components

- **MatrixCls**: Main matrix class with compile-time dimensions
- **VectorCls**: Vector class with runtime length support
- **SubMatPtrCls**: Non-owning submatrix view class
- **CommonTypes**: Type definitions for the library

### Key Features Documented

- All public methods with parameters and return values
- Usage examples for common operations
- Notes on avionics safety (no dynamic allocation)
- Template parameters and their meanings
- Algorithm descriptions (LU, QR, Cholesky decompositions)

## Troubleshooting

**Issue**: Doxygen not found
- **Solution**: Install Doxygen using the instructions above

**Issue**: Warnings about undocumented items
- **Solution**: These are normal. The code is fully documented, but Doxygen may warn about implementation details.

**Issue**: Class diagrams not showing
- **Solution**: Install Graphviz and set `HAVE_DOT = YES` in Doxyfile

**Issue**: Documentation looks incomplete
- **Solution**: Ensure `EXTRACT_ALL = YES` is set in Doxyfile (it is by default)

## Continuous Integration

You can integrate documentation generation into your build process:

```bash
# In your build script or CI/CD pipeline
doxygen Doxyfile
# Optionally: deploy docs/html/ to a web server
```

## Updating Documentation

After modifying source code comments:

1. Regenerate documentation:
   ```bash
   doxygen Doxyfile
   ```

2. Refresh your browser to view updates

The documentation is automatically generated from the comments in the source code, so keeping comments up-to-date ensures accurate documentation.

