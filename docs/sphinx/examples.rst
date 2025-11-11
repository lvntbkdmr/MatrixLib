Examples
========

This section provides usage examples for the Matrix Library.

Basic Matrix Operations
------------------------

Creating and using matrices:

.. code-block:: cpp

   #include "MatrixPkg.h"
   
   // Create a 6x6 matrix
   MatrixCls<6, 6, REAL64> A;
   A.setSubMatrixSize(4, 4);  // Use 4x4 at runtime
   
   // Initialize matrix
   A.Zeros();
   A(0, 0) = 1.0;
   A(1, 1) = 2.0;
   
   // Get submatrix view
   auto sub = A(1, 1, 2, 2);

Linear System Solving
---------------------

Solving Ax = b:

.. code-block:: cpp

   MatrixCls<6, 6, REAL64> A;
   VectorCls<6, REAL64> b, x;
   
   A.setSubMatrixSize(4, 4);
   b.setLength(4);
   
   // Set up system
   // ... populate A and b ...
   
   // Solve
   bool success = A.LeftDivide(b, x);

Vector Operations
-----------------

Working with vectors:

.. code-block:: cpp

   VectorCls<6, REAL64> vec;
   vec.setLength(3);  // Use 3 elements at runtime
   
   vec(0) = 1.0;
   vec(1) = 2.0;
   vec(2) = 3.0;
   
   // Get subvector view
   auto sub = vec(1, 2);

