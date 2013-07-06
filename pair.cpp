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

#include "pair.h"

PairGroup::PairGroup(string i, string j) {
    name_i = i; name_j = j;
    dist_eq_av = 0; dist_ab_av = 0; dist_ra_av = 0; dist_rg_av = 0;
    u_av = 0; D_av = 0; a_av = 0; k_av = 0;
    dist_eq_S = 0; dist_ra_S = 0; dist_rg_S = 0;
    u_S = 0; D_S = 0; a_S = 0; k_S = 0;
    for (unsigned int i = 0; i < P_NO; i++) P[i] = 0;
}

string PairGroup::GetNameI() { return name_i; }

string PairGroup::GetNameJ() { return name_j; }

void PairGroup::AddRa(double v) { dist_ra_list.push_back(v); }

void PairGroup::AddRg(double v) { dist_rg_list.push_back(v); }

void PairGroup::AddEq(double v) { dist_eq_list.push_back(v); }

void PairGroup::AddAb(double v) { dist_ab_list.push_back(v); }

void PairGroup::AddU(double u) { u_list.push_back(u); }

void PairGroup::AddD(double D) { D_list.push_back(D); }

void PairGroup::AddA(double A) { a_list.push_back(A); }

void PairGroup::AddK(double K) { k_list.push_back(K); }

double PairGroup::GetRa(int i) { return dist_ra_list[i]; }

double PairGroup::GetRg(int i) { return dist_rg_list[i]; }

double PairGroup::GetEq(int i) { return dist_eq_list[i]; }

double PairGroup::GetAb(int i) { return dist_ab_list[i]; }

double PairGroup::GetU(int i) { return u_list[i]; }

double PairGroup::GetD(int i) { return D_list[i]; }

double PairGroup::GetA(int i) { return a_list[i]; }

double PairGroup::GetK(int i) { return k_list[i]; }

void PairGroup::AddEq(int i, int j) {
    i_eq_list.push_back(i); j_eq_list.push_back(j);
}

int PairGroup::GetEqI(int i) { return i_eq_list[i]; }

int PairGroup::GetEqJ(int i) { return j_eq_list[i]; }

unsigned int PairGroup::GetEqSize() { return i_eq_list.size(); }

double PairGroup::GetAEq() { return dist_eq_av; }

double PairGroup::GetAAbEq() { return dist_ab_av; }

double PairGroup::GetARg() { return dist_rg_av; }

double PairGroup::GetARa() { return dist_ra_av; }

double PairGroup::GetAU() { return u_av; }

double PairGroup::GetAD() { return D_av; }

double PairGroup::GetAA() { return a_av; }

double PairGroup::GetAK() { return k_av; }

double PairGroup::GetSEq() { return dist_eq_S; }

double PairGroup::GetSRg() { return dist_rg_S; }

double PairGroup::GetSRa() { return dist_ra_S; }

double PairGroup::GetSU() { return u_S; }

double PairGroup::GetSD() { return D_S; }

double PairGroup::GetSA() { return a_S; }

double PairGroup::GetSK() { return k_S; }

void PairGroup::SetAA(double a) { a_av = a; }

void PairGroup::SetAK(double k) { k_av = k; }

void PairGroup::Average() {
    for (unsigned int j = 0; j < i_eq_list.size(); j++) {
        dist_ra_av += dist_ra_list[j];
        dist_rg_av += dist_rg_list[j];
        dist_eq_av += dist_eq_list[j];
        if (dist_ab_list.size() != 0) dist_ab_av += dist_ab_list[j];
        u_av += u_list[j];
        D_av += D_list[j];
        a_av += a_list[j];
        k_av += k_list[j];
    }
    dist_rg_av /= i_eq_list.size();
    dist_ra_av /= i_eq_list.size();
    dist_eq_av /= i_eq_list.size();
    if (dist_ab_list.size() != 0) dist_ab_av /= i_eq_list.size();
    u_av /= i_eq_list.size();
    D_av /= i_eq_list.size();
    a_av /= i_eq_list.size();
    k_av /= i_eq_list.size();
}

void PairGroup::Statistics() {
    for (unsigned int j = 0; j < i_eq_list.size(); j++) {
        dist_ra_S += pow(dist_ra_list[j]-dist_ra_av, 2);
        dist_rg_S += pow(dist_rg_list[j]-dist_rg_av, 2);
        dist_eq_S += pow(dist_eq_list[j]-dist_eq_av, 2);
        u_S += pow(u_list[j]-u_av, 2);
        D_S += pow(D_list[j]-D_av, 2);
        a_S += pow(a_list[j]-a_av, 2);
        k_S += pow(k_list[j]-k_av, 2);
    }
    if (i_eq_list.size() != 1) {
        dist_rg_S = sqrt(dist_rg_S/(i_eq_list.size()-1));
        dist_ra_S = sqrt(dist_ra_S/(i_eq_list.size()-1));
        dist_eq_S = sqrt(dist_eq_S/(i_eq_list.size()-1));
        u_S = sqrt(u_S/(i_eq_list.size()-1));
        D_S = sqrt(D_S/(i_eq_list.size()-1));
        a_S = sqrt(a_S/(i_eq_list.size()-1));
        k_S = sqrt(k_S/(i_eq_list.size()-1));
    }
    else {
        dist_rg_S = 0; dist_ra_S = 0; dist_eq_S = 0;
        u_S = 0; D_S = 0; a_S = 0; k_S = 0;
    }
}

double PairGroup::GetProb(unsigned int i) {
    if (i < P_NO) return(P[i]);
}

double PairGroup::SetProb(unsigned int i, double new_P) {
    if (i < P_NO) P[i] = new_P;
}

