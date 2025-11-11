#!/bin/bash
# Script to generate Doxygen documentation

echo "Generating Doxygen documentation..."

# Check if Doxygen is installed
if ! command -v doxygen &> /dev/null; then
    echo "Error: Doxygen is not installed."
    echo "Please install it using:"
    echo "  macOS: brew install doxygen"
    echo "  Ubuntu/Debian: sudo apt-get install doxygen"
    exit 1
fi

# Generate documentation
doxygen Doxyfile

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Documentation generated successfully!"
    echo ""
    echo "To view the documentation, open:"
    echo "  docs/html/index.html"
    echo ""
    echo "Or run:"
    echo "  open docs/html/index.html    # macOS"
    echo "  xdg-open docs/html/index.html  # Linux"
else
    echo ""
    echo "✗ Error generating documentation"
    exit 1
fi
