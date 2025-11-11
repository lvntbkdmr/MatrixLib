# Sphinx + Breathe Documentation

This guide explains how to generate documentation using Sphinx with Breathe integration.

## Overview

Sphinx + Breathe provides a modern documentation system that:
- Uses Sphinx for beautiful HTML output (Read the Docs theme)
- Integrates Doxygen comments via Breathe
- Supports multiple output formats (HTML, PDF, LaTeX)
- Provides better navigation and search

## Prerequisites

1. **Python 3** (already installed ✓)
2. **Doxygen** (already installed ✓)
3. **Sphinx and Breathe** (will be installed automatically)

## Quick Start

### Option 1: Using the Script (Recommended)

```bash
./generate_sphinx_docs.sh
```

This script will:
1. Install required Python packages
2. Generate Doxygen XML output
3. Build Sphinx documentation
4. Open the documentation in your browser

### Option 2: Manual Steps

1. **Install Python dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```

2. **Generate Doxygen XML:**
   ```bash
   doxygen Doxyfile
   ```

3. **Build Sphinx documentation:**
   ```bash
   cd docs/sphinx
   make html
   ```

4. **View the documentation:**
   ```bash
   open _build/html/index.html  # macOS
   # or
   xdg-open _build/html/index.html  # Linux
   ```

## Project Structure

```
MatrixLib/
├── docs/
│   ├── doxygen_xml/          # Doxygen XML output (generated)
│   ├── sphinx/
│   │   ├── _build/           # Sphinx output (generated)
│   │   ├── _static/          # Static files
│   │   ├── _templates/       # Custom templates
│   │   ├── conf.py           # Sphinx configuration
│   │   ├── index.rst         # Main documentation page
│   │   ├── classes.rst       # Class documentation
│   │   ├── files.rst         # File documentation
│   │   └── examples.rst      # Usage examples
│   └── html/                 # Doxygen HTML (if generated)
├── Doxyfile                  # Doxygen configuration
├── requirements.txt          # Python dependencies
└── generate_sphinx_docs.sh  # Build script
```

## Configuration

### Sphinx Configuration (`docs/sphinx/conf.py`)

Key settings:
- **Theme**: Read the Docs theme (modern and clean)
- **Breathe project**: Points to `../doxygen_xml`
- **Extensions**: Breathe, autodoc, viewcode, intersphinx

### Doxygen Configuration (`Doxyfile`)

- **XML output**: Enabled (`GENERATE_XML = YES`)
- **XML directory**: `doxygen_xml/`
- This XML is what Breathe reads to generate Sphinx docs

## Building Documentation

### HTML (Default)

```bash
cd docs/sphinx
make html
```

Output: `docs/sphinx/_build/html/index.html`

### PDF (Requires LaTeX)

```bash
cd docs/sphinx
make latexpdf
```

Output: `docs/sphinx/_build/latex/MatrixLibrary.pdf`

### Other Formats

```bash
make help  # See all available formats
```

Available formats:
- `html` - HTML documentation
- `latex` - LaTeX source
- `latexpdf` - PDF via LaTeX
- `epub` - EPUB e-book
- `man` - Manual pages
- `text` - Plain text

## Customizing Documentation

### Adding Content

Edit the RST files in `docs/sphinx/`:
- `index.rst` - Main page
- `classes.rst` - Class documentation
- `files.rst` - File documentation
- `examples.rst` - Usage examples

### Changing Theme

Edit `docs/sphinx/conf.py`:

```python
html_theme = 'sphinx_rtd_theme'  # Read the Docs theme
# or
html_theme = 'classic'           # Classic theme
# or
html_theme = 'sphinxdoc'        # Sphinxdoc theme
```

### Adding Custom Pages

1. Create a new `.rst` file in `docs/sphinx/`
2. Add it to the `toctree` in `index.rst`:

```rst
.. toctree::
   :maxdepth: 2
   
   classes
   files
   examples
   your_new_page  # Add here
```

## Troubleshooting

### Issue: "No module named 'breathe'"

**Solution**: Install dependencies:
```bash
pip3 install -r requirements.txt
```

### Issue: "Doxygen XML not found"

**Solution**: Generate Doxygen XML first:
```bash
doxygen Doxyfile
```

### Issue: "No such file or directory: _build"

**Solution**: The `_build` directory is created automatically. If missing, run:
```bash
cd docs/sphinx
make html
```

### Issue: Breathe can't find classes

**Solution**: Ensure:
1. Doxygen XML is generated (`docs/doxygen_xml/` exists)
2. `breathe_projects` in `conf.py` points to the correct path
3. Class names match exactly (case-sensitive)

## Advantages of Sphinx + Breathe

Compared to Doxygen HTML:

✅ **Better themes**: Modern, responsive design (Read the Docs theme)
✅ **Better navigation**: Improved sidebar and search
✅ **Multiple formats**: HTML, PDF, EPUB, etc.
✅ **Custom content**: Easy to add tutorials, guides, examples
✅ **Python ecosystem**: Integrates with other Python tools
✅ **Version control**: Better integration with Git workflows

## Continuous Integration

You can integrate Sphinx documentation into CI/CD:

```yaml
# Example GitHub Actions workflow
- name: Generate Documentation
  run: |
    pip3 install -r requirements.txt
    doxygen Doxyfile
    cd docs/sphinx
    make html
- name: Deploy Documentation
  uses: peaceiris/actions-gh-pages@v3
  with:
    publish_dir: ./docs/sphinx/_build/html
```

## Updating Documentation

After modifying source code comments:

1. Regenerate documentation:
   ```bash
   ./generate_sphinx_docs.sh
   ```

2. Or manually:
   ```bash
   doxygen Doxyfile
   cd docs/sphinx
   make html
   ```

The documentation automatically reflects changes in your source code comments.

