/*  MDVibCor, a program for computing vibrational corrections for gas-phase
    electron diffraction using molecular dynamics simulation data

    Copyright (C) 2008-2013 Alexander Zakharov
    E-mail: Alexander.V.Zakharov@gmail.com

    This file is part of MDVibCor program.

    MDVibCor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MDVibCor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MDVibCor.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "geometry.h"

void Normalize (xyz* v) {
    double length;
    length = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
    if (length > 0.0) {
        v->x /= length; v->y /= length; v->z /= length;
    }
    else {
        v->x = v->y = v->z = 0;
    }
}

void CrossProduct (const xyz* a, const xyz* b, xyz* aCrossB) {
    aCrossB->x = a->y * b->z - a->z * b->y;
    aCrossB->y = a->z * b->x - a->x * b->z;
    aCrossB->z = a->x * b->y - a->y * b->x;
}

double DotProduct (const xyz* a, const xyz* b) {
    return (a->x*b->x+a->y*b->y+a->z*b->z);
}

double CalcDihedral(const xyz* a1, const xyz* a2, const xyz* a3, const xyz* a4) {
    xyz offset1, offset2, offset3;
    offset1.x = a2->x-a1->x; offset1.y = a2->y-a1->y; offset1.z = a2->z-a1->z;
    offset2.x = a3->x-a2->x; offset2.y = a3->y-a2->y; offset2.z = a3->z-a2->z;
    offset3.x = a4->x-a3->x; offset3.y = a4->y-a3->y; offset3.z = a4->z-a3->z;
    Normalize(&offset1); Normalize(&offset2); Normalize(&offset3);
    xyz normal1, normal2;
    CrossProduct(&offset1, &offset2, &normal1);
    CrossProduct(&offset2, &offset3, &normal2);
    double dotPJ = -DotProduct(&offset1, &offset2);
    double dotPK = -DotProduct(&offset2, &offset3);
    double sinPJ = sqrt(1.0-dotPJ*dotPJ);
    double sinPK = sqrt(1.0-dotPK*dotPK);
    double dot = DotProduct(&normal1, &normal2)/(sinPJ*sinPK);
    if (dot > 1.0) dot += 1.0-dot;
    else if (dot < -1.0) dot -= 1.0+dot;
    double dihedral = acos(dot);
    double sense = DotProduct(&normal2, &offset1);
    if (sense < 0.0) dihedral = -dihedral;
    return(dihedral);
}

