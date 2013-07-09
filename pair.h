#ifndef PAIR_H
#define PAIR_H

#include <string>
#include <vector>
#include <math.h>

using std::string; using std::vector;

const unsigned int P_NO = 100;  // Number of probabilty function points

// Class for description of group of symmetrically equivalent atomic pairs
// (internuclear distances)
class PairGroup {
    string name_i, name_j;
    vector<double> dist_ra_list;  // List of all r_a distances in group
    vector<double> dist_rg_list;  // List of all r_g distances in group
    vector<double> dist_eq_list;  // List of all r_e distances in group
    vector<double> dist_ab_list;  // List of all ab initio r_e dist. in group
    vector<double> u_list;        // List of all r.m.s. amplitudes in group
    vector<double> D_list;        // List of all corrections in group
    vector<double> a_list;        // List of all Morse constants in group
    vector<double> k_list;        // List of all asymmetry constants in group
    vector<int> i_eq_list;        // List of equivalent atom indices
    vector<int> j_eq_list;        // List of equivalent atom indices
    double dist_eq_av, dist_ab_av, dist_ra_av, dist_rg_av;       // Av. values
    double u_av, D_av, a_av, k_av;                               // (contd)
    double dist_eq_S, dist_ra_S, dist_rg_S, u_S, D_S, a_S, k_S;  // Std. dev.
    double P[P_NO];
public:
    PairGroup(string i, string j);
    ~PairGroup() {}
    string GetNameI();
    string GetNameJ();
    void AddRa(double v);
    void AddRg(double v);
    void AddEq(double v);
    void AddAb(double v);
    void AddU(double u);
    void AddD(double D);
    void AddA(double A);
    void AddK(double K);
    double GetRa(int i);
    double GetRg(int i);
    double GetEq(int i);
    double GetAb(int i);
    double GetU(int i);
    double GetD(int i);
    double GetA(int i);
    double GetK(int i);
    void AddEq(int i, int j);
    int GetEqI(int i);
    int GetEqJ(int i);
    unsigned int GetEqSize();
    double GetAEq();
    double GetAAbEq();
    double GetARg();
    double GetARa();
    double GetAU();
    double GetAD();
    double GetAA();
    double GetAK();
    double GetSEq();
    double GetSRg();
    double GetSRa();
    double GetSU();
    double GetSD();
    double GetSA();
    double GetSK();
    void SetAA(double a);
    void SetAK(double k);
    void Average();
    void Statistics();
    double GetProb(unsigned int i);
    void SetProb(unsigned int i, double new_P);
};

#endif

