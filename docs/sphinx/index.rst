Matrix Library Documentation
============================

Welcome to the Matrix Library documentation!

This library provides a comprehensive set of matrix and vector classes designed
for avionics software where dynamic memory allocation is prohibited. All data
structures are stack-allocated with compile-time maximum dimensions and runtime
variable sizes.

.. note::

   This library is designed for safety-critical avionics applications where
   dynamic memory allocation is not allowed. All memory is stack-allocated.

Features
--------

* **Compile-time dimensions**: Maximum sizes specified at compile time for stack allocation
* **Runtime flexibility**: Use smaller dimensions at runtime than the compile-time maximum
* **No dynamic allocation**: Completely stack-allocated, avionics-safe
* **Comprehensive operations**: Matrix decompositions, linear system solving, and more
* **Submatrix views**: Efficient non-owning views into portions of matrices

Quick Start
-----------

.. code-block:: cpp

   #include "MatrixPkg.h"
   
   // Create a 6x6 matrix (max size), use 4x4 at runtime
   MatrixCls<6, 6, REAL64> A;
   A.setSubMatrixSize(4, 4);
   
   // Create vectors
   VectorCls<6, REAL64> b, x;
   b.setLength(4);
   
   // Solve Ax = b
   A.LeftDivide(b, x);

Library Components
------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   classes
   files
   examples

Class Documentation
--------------------

The main classes are documented in detail in the :doc:`classes` section.

File Documentation
------------------

.. doxygenfile:: CommonTypes.h
   :project: MatrixLib

.. doxygenfile:: MatrixCls.h
   :project: MatrixLib

.. doxygenfile:: VectorCls.h
   :project: MatrixLib

.. doxygenfile:: SubMatPtrCls.h
   :project: MatrixLib

.. doxygenfile:: MatrixPkg.h
   :project: MatrixLib

Examples
--------

The library includes comprehensive examples demonstrating various operations:

* Matrix operations (transpose, decompositions)
* Linear system solving
* Submatrix and subvector views
* Runtime dimension management

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

