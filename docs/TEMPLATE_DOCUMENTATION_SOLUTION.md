# Solution: Template Class Member Function Documentation

## Problem Solved ✓

Doxygen was not extracting documentation comments for template class member functions. This affected both:
- Doxygen HTML output
- Sphinx + Breathe output (since Breathe reads Doxygen XML)

## Solution: Use `\fn` Command at Implementation Site

**The fix:** Document template member functions at their **implementation site** (where the template is defined) using the `\fn` command.

### Before (Didn't Work):
```cpp
// In class definition (header)
/**
 * @brief Construct a bidiagonal difference matrix
 * @param result Output matrix
 */
void BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const;
```

### After (Works!):
```cpp
// In class definition (header) - just declaration, no docs
void BidiagonalDifference(MatrixCls<Length, Length, ElementType>& result) const;

// At implementation site (template definition)
/** \fn template <UINT32 Length, class ElementType> 
 *      void VectorCls<Length, ElementType>::BidiagonalDifference(...) const
 *  \brief Construct a bidiagonal difference matrix from this vector
 *  \param result Output matrix (will be resized to (SubLength-1) x SubLength)
 * 
 *  Constructs a matrix where:
 *  - Main diagonal (i, i) contains Elements[i] for i = 0 to SubLength-2
 *  - First superdiagonal (i, i+1) contains -Elements[i+1] for i = 0 to SubLength-2
 *  - All other elements are zero
 * 
 *  The resulting matrix has dimensions (SubLength-1) x SubLength.
 * 
 *  \example
 *  VectorCls<4, REAL64> vec;
 *  vec.setLength(3);
 *  vec(0) = 2.0; vec(1) = 3.0; vec(2) = 4.0;
 *  MatrixCls<3, 4, REAL64> mat;
 *  vec.BidiagonalDifference(mat);
 *  // Result: [2.0, -3.0,  0.0,  0.0]
 *  //         [0.0,  3.0, -4.0,  0.0]
 */
template <UINT32 Length, class ElementType>
void VectorCls<Length, ElementType>::BidiagonalDifference(...) const
{
    // implementation
}
```

## Why This Works

1. **`\fn` command**: Explicitly tells Doxygen which function to document
2. **Full signature**: Provides the complete template signature so Doxygen can match it
3. **Implementation site**: Doxygen sometimes extracts better from template definitions
4. **LaTeX-style commands**: Using `\brief`, `\param`, `\example` instead of `@brief`, `@param`, `@example`

## Key Points

- **This is a Doxygen limitation**, not a Sphinx/Breathe issue
- **Both Doxygen HTML and Sphinx** now show the documentation correctly
- **The `\fn` approach** is the recommended workaround for template member functions
- **Documentation appears in both outputs** once Doxygen extracts it properly

## Verification

After applying this fix:

1. **Doxygen XML** contains the documentation:
   ```bash
   grep -A 50 "BidiagonalDifference" docs/doxygen_xml/classVectorCls.xml
   ```

2. **Doxygen HTML** shows the documentation:
   ```bash
   open docs/html/classVectorCls.html
   ```

3. **Sphinx + Breathe** shows the documentation:
   ```bash
   open docs/sphinx/_build/html/classes.html
   ```

## For Other Template Functions

Apply the same pattern to any other template member functions that aren't being documented:

1. Remove documentation from the class definition (header)
2. Add `\fn` documentation at the template implementation site
3. Use LaTeX-style commands (`\brief`, `\param`, etc.)
4. Regenerate documentation

## References

- [Doxygen `\fn` Command](https://www.doxygen.nl/manual/commands.html#cmdfn)
- [Doxygen Template Documentation](https://www.doxygen.nl/manual/commands.html#cmdtparam)

