#!/bin/bash
# Script to generate Sphinx documentation with Breathe

set -e  # Exit on error

echo "=========================================="
echo "Generating Sphinx Documentation"
echo "=========================================="
echo ""

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed."
    echo "Please install Python 3 first."
    exit 1
fi

# Check if pip is installed
if ! command -v pip3 &> /dev/null; then
    echo "Error: pip3 is not installed."
    echo "Please install pip3 first."
    exit 1
fi

# Install/upgrade required packages using virtual environment
echo "Step 1: Setting up Python virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Install/upgrade required packages
echo "Installing Python dependencies..."
pip install -q --upgrade -r requirements.txt

# Check if Doxygen is installed
if ! command -v doxygen &> /dev/null; then
    echo "Error: Doxygen is not installed."
    echo "Please install it using: brew install doxygen"
    exit 1
fi

# Generate Doxygen XML (required for Breathe)
echo ""
echo "Step 2: Generating Doxygen XML output..."
doxygen Doxyfile > /dev/null 2>&1

if [ ! -d "docs/doxygen_xml" ]; then
    echo "Error: Doxygen XML output not found."
    echo "Please check Doxygen configuration."
    exit 1
fi

echo "✓ Doxygen XML generated successfully"

# Change to Sphinx directory
cd docs/sphinx

# Build Sphinx documentation
echo ""
echo "Step 3: Building Sphinx documentation..."
# Use sphinx-build directly (works better with venv)
sphinx-build -b html . _build/html

if [ $? -eq 0 ]; then
    # Deactivate virtual environment
    deactivate
    
    echo ""
    echo "=========================================="
    echo "✓ Documentation generated successfully!"
    echo "=========================================="
    echo ""
    echo "To view the documentation, open:"
    echo "  docs/sphinx/_build/html/index.html"
    echo ""
    echo "Or run:"
    echo "  open docs/sphinx/_build/html/index.html    # macOS"
    echo "  xdg-open docs/sphinx/_build/html/index.html  # Linux"
    echo ""
else
    deactivate
    echo ""
    echo "✗ Error generating documentation"
    exit 1
fi

