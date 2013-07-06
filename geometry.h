#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>

// Structure holding x, y, z coordinates of an atom
struct xyz {
    double x, y, z;
};

// Functions for computing dihedral angle value
void Normalize (xyz* v);
void CrossProduct (const xyz* a, const xyz* b, xyz* aCrossB);
double DotProduct (const xyz* a, const xyz* b);
double CalcDihedral(const xyz* a1, const xyz* a2, const xyz* a3, const xyz* a4);

#endif

