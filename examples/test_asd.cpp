#include "MatrixCls.h"
#include "VectorCls.h"
#include "SubMatPtrCls.h"
#include "SubVecPtrCls.h"
#include "MatrixPkg.h"

int main()
{
    MatrixCls<6, 6, double> L_full;
    VectorCls<6, double> b_full;
    VectorCls<6, double> x_full;

    size_t size = 3;  // Runtime decision

    SubMatPtrCls<double> L_sub = L_full.submatrix(0, 0, size, size);
    SubVecPtrCls<double> b_sub = b_full.subvector(0, size);
    SubVecPtrCls<double> x_sub = x_full.subvector(0, size);

    x_sub = L_sub.forwardSubstitution(b_sub);

    x_full = L_full.forwardSubstitution(b_full);

    x_sub = L_sub.backwardSubstitution(b_sub);

    x_full = L_full.backwardSubstitution(b_full);

    MatrixCls<6, 6, double> L2_full;

    L2_full = L_full.transpose();

    SubMatPtrCls<double> L2_sub = L2_full.submatrix(0, 0, size, size);

    L2_sub = L_sub.transpose<3, 3>();

    // Lx = b
    // Option 1: Write directly to SubVecPtrCls (requires compile-time size template parameter)
    L2_sub.leftDivide<3>(b_sub, x_sub);
    
    // Option 2: Write to full VectorCls - size is automatically inferred from runtime dimensions!
    // No template parameter needed - OutputSize (6) is deduced from x_full type
    L2_sub.leftDivide(b_sub, x_full);
    SubVecPtrCls<double> x_sub2 = x_full.subvector(0, size);

    x_full = L2_full.leftDivide(b_full);

    return 0;
}

