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

#include "atom.h"

Atom::Atom(string nname, double nx, double ny, double nz) {
    name = nname; eq.x = nx; eq.y = ny; eq.z = nz;
}

const xyz* Atom::GetC() {
    return &eq;
}

const xyz* Atom::GetAbinitCoord() {
    return &eq_ab;
}

string Atom::GetName() {
    return name;
}

const xyz* Atom::GetKthCoord(int step) {
    return &trajectory[step];
}

long int Atom::GetTrajSize() {
    return trajectory.size();
}

void Atom::SetAbinitCoord (double nx, double ny, double nz) {
    eq_ab.x = nx; eq_ab.y = ny; eq_ab.z = nz;
}

void Atom::AddStep(const xyz* c) {
    trajectory.push_back(*c);
}

void Atom::SetTrajCoord(int step, const xyz* c) {
    trajectory[step].x = c->x;
    trajectory[step].y = c->y;
    trajectory[step].z = c->z;
}

double Atom::GetEqDist(Atom* a) {
    const xyz* c = a->GetC();
    double xixj = eq.x-c->x;
    double yiyj = eq.y-c->y;
    double zizj = eq.z-c->z;
    return (sqrt(xixj*xixj+yiyj*yiyj+zizj*zizj));
}
double Atom::GetAbinitEqDist(Atom* a) {
    const xyz* c = a->GetAbinitCoord();
    double xixj = eq_ab.x-c->x;
    double yiyj = eq_ab.y-c->y;
    double zizj = eq_ab.z-c->z;
    return (sqrt(xixj*xixj+yiyj*yiyj+zizj*zizj));
}

double Atom::GetDist(int step, Atom* a) {
    const xyz* c = a->GetKthCoord(step);
    double xixj = trajectory[step].x-c->x;
    double yiyj = trajectory[step].y-c->y;
    double zizj = trajectory[step].z-c->z;
    return (sqrt(xixj*xixj+yiyj*yiyj+zizj*zizj));
}

