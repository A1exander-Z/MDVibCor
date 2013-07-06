#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include "atom.h"
#include "pair.h"
#include "geometry.h"

using namespace std;

const unsigned int CONF_NO = 10;  // Maximum number of conformers

// Class describing a group of different conformers defined by dihedral
class ConfGroup {
    int* dihedral;          // Dihedral angle defining a conformer
    vector<double> eq_dih;  // Equilibrium values of dihedral angles
public:
    ConfGroup() { dihedral = NULL; }
    ~ConfGroup() { if (dihedral) delete[] dihedral; }
    unsigned int GetEqDihSize() { return eq_dih.size(); }
    double GetEqDih(int i) { return eq_dih[i]; }
    const int* GetDih() { return dihedral; }
    void SetDih(int* new_dih ) { dihedral = new_dih; }
    void AddEqDih(double new_dih) { eq_dih.push_back(new_dih); }
};

// Class describing the whole simulation, may include several trajectories
class Simulation {
    int* renumber;                    // Renumbering table
    unsigned int atom_no;             // Number of atoms in a molecule
    unsigned long int tr_size;        // Trajectory size
    vector<Atom*> atom_list;          // List of atoms
    vector< vector<int> > rot_list;   // List of rotating atom groups
    vector< vector<int> > eq_list;    // List of groups of equivalent atoms
    vector<ConfGroup*> conf_list;     // List of conformer groups
    vector<int*> conformer;           // Number of conformer in given MD step
    vector< vector<int> > internals_list;  // List of internal coordinates
    long int conf_dist[360];          // Conformer distribution over angle
    vector< vector<double> > rot_equil;    // Equilibrium values of dihedrals
    vector<PairGroup*> pair_list;     // List of groups of equiv. atomic pairs
    double* dist_eq;    // Matrix of distances at the equilibrium geometry
    double* dist_eq_ab; // Matrix of distances at the ab initio eq. geometry
    double* dist_ra;    // Matrix of r_a distances
    double* dist_rg;    // Matrix of r_g distances
    double* u;          // Matrix of r.m.s. amplitudes
    double* a_M;        // Matrix of Morse constants
    double* kappa;      // Matrix of asymmetry constants
    double tolerance;   // Tolerance for symmetrically eqivalent distances
    bool abinit;        // Is ab initio geometry present?
    int chosen_conf[CONF_NO];
    double energy, eq_energy, T;
    double CalcEqDist(int i, int j);
    double CalcAbInitEqDist(int i, int j);
    double CalcTrDist(int i, int j, int k);
    bool EquivalentAtoms(int no_1, int no_2);
    bool EquivalentDistances(int no_1, int no_2);
    int CheckRotation(int step);
    int GetConformerNo(ConfGroup* conf_group, double* angle);
    void CheckConformer(int step, double* angle, int* conformation);
public:
    Simulation();
    ~Simulation();
    void SetTolerance(double new_tol);
    int ReadEqGeometry(const char* filename);
    void PrintCoordinates();
    int ReadAbInitioGeometry(const char* filename);
    void PrintAbInitioCoordinates();
    int ReadRenumbering(const char* filename);
    void PrintRenumbering();
    int ReadEquivalent(const char* filename);
    void PrintEquivalent();
    int ReadRotations(const char* filename);
    void PrintRotations();
    int ReadConformers(const char* filename);
        int ReadInternals(const char* filename);
    int ReadTrajectory(const char* filename, int skip, unsigned int* PI_steps);
    unsigned long int GetTrajectorySize();
    void ComformationalAnalysis();
    void CalcRaRg(unsigned int PI_steps);
    void CalcProbabilities(bool print_P, bool debug);
    bool CompareConf(int step);
    void PrintGroups();
    void PrintEDatED();
    void PrintKCED();
    void PrintUNEX();
};

#endif

