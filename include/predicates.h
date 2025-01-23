
#ifndef ORIENT2D_H
#define ORIENT2D_H

#ifdef __cplusplus
extern "C" {
#endif

typedef double REAL; // Assuming REAL is a double, update as needed

void exactinit();

// Function declaration
REAL orient2d(REAL *pa, REAL *pb, REAL *pc);

REAL orient2dexact(REAL *pa, REAL *pb, REAL *pc);
REAL orient2dslow(REAL *pa, REAL *pb, REAL *pc);

#ifdef __cplusplus
}
#endif

#endif // ORIENT2D_H
