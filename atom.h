#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include "geometry.h"

using std::string; using std::vector;

// Class for description of an atom in a molecule
class Atom {
    string name;             // Atom name
    xyz eq;                  // Coordinates in the equilibrium geometry
    xyz eq_ab;               // Coordinates in the ab initio eq. geometry
    vector<xyz> trajectory;  // Coordinates for all steps of the trajectory
public:
    Atom(string nname, double nx, double ny, double nz);
    ~Atom() {}
    const xyz* GetC();
    const xyz* GetAbinitCoord();
    string GetName();
    const xyz* GetKthCoord(int step);
    long int GetTrajSize();
    void SetAbinitCoord (double nx, double ny, double nz);
    void AddStep(const xyz* c);
    void SetTrajCoord (int step, const xyz* c);
    double GetEqDist(Atom* a);          // Eqilibrium distance
    double GetAbinitEqDist(Atom* a);    // Ab initio equilibrium distance
    double GetDist(int step, Atom* a);  // Distance at kth step
};

#endif

