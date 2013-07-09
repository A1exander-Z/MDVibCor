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

#include "simulation.h"

// MS Visual Studio does not define 'round' function, so we have to provide one
#ifdef _MSC_VER
    double round(double d) {
        return floor(d + 0.5);
    }
#endif

// Parameters of function for approximating probability distributions
class P_dist_par {
public:
    double P_step, min, eq;
    double* PP;     // Pointer to current PairGroup.P array
    double* VP;     // Pointer to V(r) array
    double* rP;     // Pointer to r-r_e array
};

// Function for approximating probability distributions using
// Nelder-Mead algorithm (from GSL)
double P_dist (const gsl_vector *v, void *params) {
    double A, a0, sigma, c1, c2;
    P_dist_par* p = (P_dist_par*)params;
    a0 = gsl_vector_get(v, 0); A = gsl_vector_get(v, 1);
    sigma = gsl_vector_get(v, 2); c1 = gsl_vector_get(v, 3);
    c2 = gsl_vector_get(v, 4);
    double sum = 0;
    for (int i = 0; i < P_NO; i++) {
        double a = i*p->P_step+p->min;
        double P = A*exp(-(pow(a-a0, 2)/(2*pow(sigma, 2))))*
                   (1+c1*(a-a0)+c2*pow(a-a0, 2));
        sum += pow(p->PP[i]-P, 2);
    }
    return(sum);
}

// Function for approximating potential by Morse potential
// using Nelder-Mead algorithm (from GSL)
double V_dist (const gsl_vector *v, void *params) {
    double a, D;
    P_dist_par* p = (P_dist_par*)params;
    D = gsl_vector_get(v, 0); a = gsl_vector_get(v, 1);
    double sum = 0;
    for (int i = 0; i < P_NO; i++) {
        double V = D*pow((1-exp(-a*p->rP[i])), 2);
        sum += pow(p->VP[i]-V, 2);
    }
    return(sum);
}

Simulation::Simulation() {
    renumber = NULL;
    dist_eq = dist_ra = dist_rg = dist_eq_ab = NULL;
    atom_no = tr_size = 0; energy = eq_energy = T = 0;
    abinit = false;
}

Simulation::~Simulation() {
    for (unsigned int i = 0; i < conformer.size(); i++) {
        delete[] conformer[i]; conformer[i] = NULL;
    }
    conformer.clear();
    for (unsigned int i = 0; i < conf_list.size(); i++) {
        delete conf_list[i]; conf_list[i] = NULL;
    }
    conf_list.clear();
    if (!renumber) {
        delete[] renumber; renumber = NULL;
    }
    if (!dist_eq) {
        delete[] dist_eq; dist_eq = NULL;
    }
    if (!dist_ra) {
        delete[] dist_ra; dist_ra = NULL;
    }
    if (!dist_rg) {
        delete[] dist_rg; dist_rg = NULL;
    }
    if (abinit && !dist_eq_ab) {
        delete[] dist_eq_ab; dist_eq_ab = NULL;
    }
    delete[] u; u = NULL; delete[] a_M; a_M = NULL;
    delete[] kappa; kappa = NULL;
    for (unsigned int i = 0; i < pair_list.size(); i++) {
        delete pair_list[i]; pair_list[i] = NULL;
    }
    pair_list.clear();
    for (unsigned int i = 0; i < atom_list.size(); i++) {
        delete atom_list[i]; atom_list[i] = NULL;
    }
    atom_list.clear();
}

void Simulation::SetTolerance(double new_tol) {
    tolerance = new_tol;
}

// Function to calculate equilibrium internuclear distance between i, j pair
double Simulation::CalcEqDist(int i, int j) {
    return(atom_list[i]->GetEqDist(atom_list[j]));
}

// Function to calculate ab initio equilibrium
// internuclear distance between i, j pair
double Simulation::CalcAbInitEqDist(int i, int j) {
    return(atom_list[i]->GetAbinitEqDist(atom_list[j]));
}

// Function to calculate internuclear distance between i, j pair at
// kth step of the MD trajectory
double Simulation::CalcTrDist(int i, int j, int k) {
    return(atom_list[i]->GetDist(k, atom_list[j]));
}

// Function to check equivalence of two atoms using the equivalence table
bool Simulation::EquivalentAtoms(int no_1, int no_2) {
    for (unsigned i = 0; i < eq_list.size(); i++)        // Cycle over rows
        for (unsigned j = 0; j < eq_list[i].size(); j++) // Search in ith row
            if (no_1 == eq_list[i][j])
                // No. 1 is present in the ith row of the table
                for (unsigned k = 0; k < eq_list[i].size(); k++)
                    if (no_2 == eq_list[i][k])
                        // No. 2 is present in the same row
                        return(true);
    return(false);
}

// Function to check the equivalence of two distances using equilibrium
// geometry (MD or ab initio)
bool Simulation::EquivalentDistances(int no_1, int no_2) {
    double diff;
    if (abinit) {  // Ab initio geometry present
        diff = fabs(dist_eq_ab[no_1]-dist_eq_ab[no_2]);
        return(diff <= dist_eq_ab[no_1]*(tolerance/100.0) ? true : false);
    }
    else {         // MD geometry
        diff = fabs(dist_eq[no_1]-dist_eq[no_2]);
        return(diff <= dist_eq[no_1]*(tolerance/100.0) ? true : false);
    }
}

// Function to get the number of conformer in its group of conformers
// (distinguished by a dihedral angle) using the angle value
int Simulation::GetConformerNo(ConfGroup* conf_group, double* angle) {
    int number;
    double min = 2*M_PI;
    // Search for conformer with minimal angle difference
    for (unsigned int i = 0; i < conf_group->GetEqDihSize(); i++) {
        // Sine and cosine of the given angle
        double sin_x = sin(*angle);
        double cos_x = cos(*angle);
        // Sine and cosine of the angle of ith conformer in the group
        double sin_y = sin(conf_group->GetEqDih(i)/180.0*M_PI);
        double cos_y = cos(conf_group->GetEqDih(i)/180.0*M_PI);
        double dev_2 = acos(cos_x*cos_y+sin_x*sin_y);   // Diff. of angles
        if (dev_2 < min) { min = dev_2; number = i; }
    }
    // Same as above, but for mirrored configuration
    // (by changing signs of the angles) FIXME (Is it _always_ necessary?)
    for (unsigned int i = 0; i < conf_group->GetEqDihSize(); i++) {
        double sin_x = sin(*angle);
        double cos_x = cos(*angle);
        double sin_y = sin(-conf_group->GetEqDih(i)/180.0*M_PI);
        double cos_y = cos(-conf_group->GetEqDih(i)/180.0*M_PI);
        double dev_2 = acos(cos_x*cos_y+sin_x*sin_y);
        if (dev_2 < min) { min = dev_2; number = i; }
    }
    return(number);
}


// Function to get the number of conformer at the given step of
// the MD trajectory
void Simulation::CheckConformer(int step, double* angle, int* conf) {
    for (unsigned int i = 0; i < conf_list.size(); i++) {
        const int* conformer_dihedral = conf_list[i]->GetDih();
        int a1i = conformer_dihedral[0]-1;      // First atom
        const xyz* a1 = atom_list[a1i]->GetKthCoord(step);
        int a2i = conformer_dihedral[1]-1;      // Second atom
        const xyz* a2 = atom_list[a2i]->GetKthCoord(step);
        int a3i = conformer_dihedral[2]-1;      // Third atom
        const xyz* a3 = atom_list[a3i]->GetKthCoord(step);
        int a4i = conformer_dihedral[3]-1;      // Fourth atom
        const xyz* a4 = atom_list[a4i]->GetKthCoord(step);
        *angle = CalcDihedral(a1, a2, a3, a4);
        conf[i] = GetConformerNo(conf_list[i], angle);
    }
}

// Function for treating internal rotations
// ('step' is the number of the MD trajectory step)
// FIXME This function needs refactoring...
int Simulation::CheckRotation(int step) {
    int rot_total = 0;
    for (unsigned int i = 0; i < rot_list.size(); i++) {
        double rot_r = 0;   // Parameters characterising rotation, CW and CCW
        double rot_l = 0;
        double no_rot = 0;
        // Tables of atoms in rotating group
        int* rot_from = new int[rot_list.size()-3];
        int* rot_to_r = new int[rot_list.size()-3];
        int* rot_to_l = new int[rot_list.size()-3];
        int a2i = rot_list[i][rot_list[i].size()-3]-1;  // Second atom
        const xyz* a2 = atom_list[a2i]->GetKthCoord(step);
        int a3i = rot_list[i][rot_list[i].size()-2]-1;  // Third atom
        const xyz* a3 = atom_list[a3i]->GetKthCoord(step);
        int a4i = rot_list[i][rot_list[i].size()-1]-1;  // Fourth atom
        const xyz* a4 = atom_list[a4i]->GetKthCoord(step);
        // Deviations from equilibrium angles
        for (unsigned int j = 0; j < rot_list[i].size()-3; j++) {
            int a1i = rot_list[i][j]-1;         // First atom
            const xyz* a1 = atom_list[a1i]->GetKthCoord(step);
            double angle = CalcDihedral(a1, a2, a3, a4);
            double sin_x = sin(angle);
            double cos_x = cos(angle);
            double sin_y = sin(rot_equil[i][j]);
            double cos_y = cos(rot_equil[i][j]);
            double dev_2 = acos(cos_x*cos_y+sin_x*sin_y);
            no_rot += dev_2*dev_2; rot_from[j] = a1i;
        }
        // Check shift to the right (1->2, 2->3, ..., n->1) 
        for (unsigned int j = 0; j < rot_list[i].size()-3; j++) {
            int a1i_n, j_n;
            (j < rot_list[i].size()-4) ? j_n = j+1 : j_n = 0;
            a1i_n = rot_list[i][j_n]-1;
            const xyz* a1_n = atom_list[a1i_n]->GetKthCoord(step);
            double angle_n = CalcDihedral(a1_n, a2, a3, a4);
            double sin_x = sin(angle_n);
            double cos_x = cos(angle_n);
            double sin_y = sin(rot_equil[i][j]);
            double cos_y = cos(rot_equil[i][j]);
            double dev_2 = acos(cos_x*cos_y+sin_x*sin_y);
            rot_r += dev_2*dev_2; rot_to_r[j] = a1i_n;
        }
        // Check shift to the left (2->1, 3->2, ..., 1->n)
        for (unsigned int j = 0; j < rot_list[i].size()-3; j++) {
            int a1i_p, j_p;
            (j > 0) ? j_p = j-1 : j_p = rot_list[i].size()-4;
            a1i_p = rot_list[i][j_p]-1;
            const xyz* a1_p = atom_list[a1i_p]->GetKthCoord(step);
            double angle_p = CalcDihedral(a1_p, a2, a3, a4);
            double sin_x = sin(angle_p);
            double cos_x = cos(angle_p);
            double sin_y = sin(rot_equil[i][j]);
            double cos_y = cos(rot_equil[i][j]);
            double dev_2 = acos(cos_x*cos_y+sin_x*sin_y);
            rot_l += dev_2*dev_2; rot_to_l[j] = a1i_p;
        }
        double rot = 0;
        int* rot_to;
        if (rot_r < rot_l) { rot = rot_r; rot_to = rot_to_r; }
        else { rot = rot_l; rot_to = rot_to_l; }
        // Actual exchange of coordinates
        if (rot < no_rot) {
            rot_total++;
            xyz* rotated = new xyz[rot_list[i].size()-3];
            // Storing coordinates
            for (unsigned int j = 0; j < rot_list[i].size()-3; j++)
                rotated[j] = *atom_list[rot_to[j]]->GetKthCoord(step);
            // Now the real exchange of coordinates
            for (unsigned int j = 0; j < rot_list[i].size()-3; j++)
                atom_list[rot_from[j]]->SetTrajCoord(step, &rotated[j]);
            delete[] rotated;
        }
        // Checking for tree-like structures (like rotating tert-butyl groups
        // with rotating methyl groups within them)
        if (rot < no_rot) {
            // Building indexes
            vector< vector<int> > tree_indexes;
            vector<int> bonded_index;
            for (unsigned int j = i+1; j < rot_list.size(); j++) {
                int next_level = rot_list[j][rot_list[j].size()-3];
                for (unsigned int k = 0; k < rot_list[i].size()-3; k++)
                    if (next_level == rot_list[i][k]) {
                        vector<int> tree_entry;
                        for (unsigned int n = 0; n < rot_list[j].size()-3;
                             n++) {
                            tree_entry.push_back(rot_list[j][n]);
                        }
                        tree_indexes.push_back(tree_entry);
                        bonded_index.push_back(next_level);
                    }
            }
            if (tree_indexes.size() != 0) {
                // Doing real exchange of coordinates
                xyz* coord_matrix = new xyz[tree_indexes.size()*
                                        tree_indexes[0].size()];
                for (unsigned int j = 0; j < tree_indexes.size(); j++) {
                    vector<int> tree_entry = tree_indexes[j];
                    int jm = j-(int)(j/(rot_list[i].size()-3))*
                             (rot_list[i].size()-3);
                    for (unsigned int k = j-jm; k < bonded_index.size(); k++)
                        if (bonded_index[k] == rot_from[jm]+1) {
                            for (unsigned int n = 0;
                                 n < tree_entry.size(); n++) {
                                coord_matrix[k*tree_entry.size()+n] =
                                *atom_list[tree_entry[n]-1]->GetKthCoord(step);
                            }
                            break;
                        }
                }
                for (unsigned int j = 0; j < tree_indexes.size(); j++) {
                    vector<int> tree_entry = tree_indexes[j];
                    int jm = j-(int)(j/(rot_list[i].size()-3))*
                             (rot_list[i].size()-3);
                    for (unsigned int k = j-jm; k < bonded_index.size(); k++)
                        if (bonded_index[k] == rot_to[jm]+1) {
                            for (unsigned int n = 0;
                                 n < tree_entry.size(); n++) {
                                atom_list[tree_entry[n]-1]->SetTrajCoord(step,
                                &coord_matrix[k*tree_entry.size()+n]);
                            }
                            break;
                        }
                }
                delete[] coord_matrix;
            }
        }
        delete[] rot_from; delete[] rot_to_l; delete[] rot_to_r;
    }
    return(rot_total);
}

// Are two conformers the same?
bool Simulation::CompareConf(int step) {
    for (unsigned int i = 0; i < conf_list.size(); i++)
        if (conformer[step][i] != chosen_conf[i]) return(false);
    return(true);
}

// Function for reading equilibrium geometry (in XYZ format) from
// the MD trajectory file
int Simulation::ReadEqGeometry(const char* filename) {
    ifstream in_tr(filename);
    if (in_tr.fail()) {
        cout << "Cannot open '" << filename << "' trajectory file." << endl;
        in_tr.close(); return(0);
    }
    string str_in;
    getline(in_tr, str_in);
    istringstream stream_in(str_in);
    // Reading the number of atoms
    if (!(stream_in >> atom_no)) {
        // Something is wrong, no valid number of atoms in the first line
        cout << "'" << filename << "' does not seem to be a trajectory file."
             << endl;
        in_tr.close(); return(0);
    }
    getline(in_tr, str_in); getline(in_tr, str_in);
    // Reading equilibrium geometry
    for (unsigned int i = 0; i < atom_no && !in_tr.eof(); i++) {
        istringstream stream_in(str_in);
        string name;
        double x, y, z;         
        stream_in >> name >> x >> y >> z;  // Reading coordinates
        if (name != "") atom_list.push_back(new Atom(name, x, y, z));
        else {
            cout << "Error in the atom list of the equilibrium geometry"
                 << endl;
            in_tr.close(); return(0);
        }
        getline(in_tr, str_in);
    }
    in_tr.close();
    if (atom_no != atom_list.size()) {
        // Something is wrong, number of atoms read is not equal to the given
        cout << "Error in the atom list of the equilibrium geometry" << endl;
        return(0);
    }
    // Allocating renumbering table
    renumber = new int[atom_no+1];
    // Filling renumbering table with default values (no renumbering)
    for (unsigned int i = 0; i <= atom_no; i++) renumber[i] = i;
    return(atom_no);
}

// Function for printing equilibrium (MD) coordinates of atoms
void Simulation::PrintCoordinates() {
    cout << endl << "Atomic coordinates:" << endl;
    for (unsigned int i = 0; i < atom_no; i++) {
        const xyz* C = atom_list[i]->GetC();
        cout << atom_list[i]->GetName() << "\t" << 
                C->x << "\t" << C->y << "\t" << C->z << endl;
    }
    cout << endl;
}

// Function for reading ab initio equilibrium geometry from an XYZ file
int Simulation::ReadAbInitioGeometry(const char* filename) {
    ifstream in_ab(filename);
    if (in_ab.fail()) {
        cout << "Cannot open '" << filename << "' ab initio geometry." << endl;
        in_ab.close(); return(0);
    }
    string str_in;                  // String holding input lines
    getline(in_ab, str_in);         // Reading number of atoms
    unsigned int atom_no_ab;
    istringstream stream_in(str_in); stream_in >> atom_no_ab;
    if (atom_no_ab == 0) {
        // Something is wrong, no valid number of atoms in the first line
        cout << "'" << filename << "' does not seem to be an XYZ file" << endl;
        in_ab.close(); return(0);
    }
    if (atom_no_ab != atom_no) {
        // Number of atoms is not equal to the number of atoms in the MD traj.
        cout << "Number of atoms is not equal to the number of atoms in the "
             << "MD trajectory" << endl;
        in_ab.close(); return(0);
    }
    getline(in_ab, str_in); getline(in_ab, str_in); // Skipping one line
    // Reading equilibrium geometry
    unsigned int i;
    for (i = 0; i < atom_no && !in_ab.eof(); i++) {
        istringstream stream_in(str_in);
        string name;
        double x, y, z;         
        stream_in >> name >> x >> y >> z;   // Reading coordinates
        // Adding atom coordinates to list
        if (name == atom_list[i]->GetName())
            atom_list[i]->SetAbinitCoord(x, y, z);
        else {
            cout << "One of the atom names in '" << filename
                 << "' does not match to name in the trajectory file" << endl;
            in_ab.close(); return(0);
        }
        getline(in_ab, str_in);
    }
    in_ab.close();
    if (i != atom_no) {
        // Something is wrong, number of atoms read is not equal to the given
        cout << "Error in the atom list of the ab initio geometry" << endl;
        return(0);
    }
    abinit = true; return(atom_no);
}

// Function for printing equilibrium (ab initio) coordinates of atoms
void Simulation::PrintAbInitioCoordinates() {
    cout << endl << "Atomic coordinates:" << endl;
    for (unsigned int i = 0; i < atom_no; i++) {
        const xyz* C = atom_list[i]->GetAbinitCoord();
        cout << atom_list[i]->GetName() << "\t" <<
                C->x << "\t" << C->y << "\t" << C->z << endl;
    }
    cout << endl;
}

// Function for reading renumbering table from file
int Simulation::ReadRenumbering(const char* filename) {
    ifstream in_renum(filename);
    if (in_renum.fail()) {
        cout << "Cannot open '" << filename << "' renumbering file." << endl;
        in_renum.close(); return(0);
    }
    string str_in;
    getline(in_renum, str_in);
    while (!in_renum.eof()) {
        istringstream stream_in(str_in);
        unsigned int in_no_value = 0;
        unsigned int out_no_value = 0;
        stream_in >> in_no_value >> out_no_value;  // Reading numbers
        // Numbers should be non-zero
        if (in_no_value != 0 && out_no_value != 0) {
            // Overflow check
            if (in_no_value > atom_no)
                cout << "One of numbers in the renumbering table is larger "
                     << "than the number of atoms in the molecule." << endl;
            else renumber[in_no_value] = out_no_value;
        }
        getline(in_renum, str_in);
    }
    in_renum.close(); return(0);
}

// Function for printing the renumbering table
void Simulation::PrintRenumbering() {
    cout << endl << "Renumbering table:" << endl;
    for (unsigned int i = 1; i <= atom_no; i++)
        cout << i << "\t" << renumber[i] << endl;
    cout << endl;
}

// Function for reading the lsit of equivalent atoms from file
int Simulation::ReadEquivalent(const char* filename) {
    ifstream in_equiv(filename);
    if (in_equiv.fail()) {
        cout << "Cannot open '" << filename << "' file with eqivalent atoms."
             << endl;
        in_equiv.close(); return(0);
    }
    string str_in;
    getline(in_equiv, str_in);
    while (!in_equiv.eof()) {
        istringstream stream_in(str_in);
        vector<int> equiv_atoms;  // Creating group of equivalent atoms
        int equiv_atom_no;
        while (stream_in >> equiv_atom_no)
            equiv_atoms.push_back(equiv_atom_no); // Adding atom
        if (equiv_atoms.size() != 0)  // Adding a group
             eq_list.push_back(equiv_atoms);
        getline(in_equiv, str_in);
    }
    in_equiv.close();
    return(eq_list.size());
}

// Function for printing the list of equivalent atoms
void Simulation::PrintEquivalent() {
    cout << endl << "List of groups of equivalent atoms:" << endl;
    for (unsigned int i = 0; i < eq_list.size(); i++)
        for (unsigned int j = 0; j < eq_list[i].size(); j++)
            cout << eq_list[i][j] << " "
                 << (j < eq_list[i].size()-1 ? " " : "\n");
    cout << endl;
}

// Function for reading the list of internal rotations from file
int Simulation::ReadRotations(const char* filename) {
    ifstream in_rot(filename);
    if (in_rot.fail()) {
        cout << "Cannot open '" << filename << "' file with rotations."
             << endl;
        in_rot.close(); return(0);
    }
    string str_in;
    getline(in_rot, str_in);
    while (!in_rot.eof()) {
        istringstream stream_in(str_in);
        vector<int> rot_atoms;     // Atoms in a group which may rotate
        int rot_atom_no = 0;
        int old_rot_atom_no = 0;
        stream_in >> rot_atom_no;  // Reading atom number
        while (rot_atom_no != old_rot_atom_no) {
            rot_atoms.push_back(rot_atom_no);
            old_rot_atom_no = rot_atom_no;
            stream_in >> rot_atom_no;
        }
        // There should be at least two atoms in rotating group and
        // the last 3rd for defining the dihedral angle
        if (rot_atoms.size() > 3) { 
            rot_list.push_back(rot_atoms);
            vector<double> equil_dih;
            int a2i = rot_atoms[rot_atoms.size()-3]-1;  // Second atom
            const xyz* a2 = atom_list[a2i]->GetC();
            int a3i = rot_atoms[rot_atoms.size()-2]-1;  // Third atom
            const xyz* a3 = atom_list[a3i]->GetC();
            int a4i = rot_atoms[rot_atoms.size()-1]-1;  // Fourth atom
            const xyz* a4 = atom_list[a4i]->GetC();
            for (unsigned int j = 0; j < rot_atoms.size()-3; j++) {
                int a1i = rot_atoms[j]-1;               // First atom
                const xyz* a1 = atom_list[a1i]->GetC();
                double angle = CalcDihedral(a1, a2, a3, a4);
                equil_dih.push_back(angle);  // Storing equilibrium dihedral
            }
            rot_equil.push_back(equil_dih);
        }
        getline(in_rot, str_in);
    }
    in_rot.close(); return(rot_list.size());
}

// Function for printing the list of rotations
void Simulation::PrintRotations() {
    cout << endl << "List of rotating atom groups:" << endl;
    for (unsigned int i = 0; i < rot_list.size(); i++) {
        cout << "[";
        for (unsigned int j = 0; j < rot_list[i].size()-3; j++)
            cout << rot_list[i][j]
                 << (j < rot_list[i].size()-3 ? " " : "]");
        for (unsigned int j = rot_list[i].size()-3; j < rot_list[i].size(); j++)
            cout << "-" << rot_list[i][j];
        cout << " (equilibrium values: ";
        for (unsigned int j = 0; j < rot_list[i].size()-3; j++)
            cout << rot_equil[i][j]*180.0/M_PI
                 << (j < rot_list[i].size()-3 ? " " : ")\n");
    }
    cout << endl;
}

// Function for reading the list of conformers from file
int Simulation::ReadConformers(const char* filename) {
    ifstream in_conformers(filename);
    if (in_conformers.fail()) {
        cout << "Cannot open '" << filename << "' file with conformers."
             << endl;
        in_conformers.close(); return(0);
    }
    string str_in;
    getline(in_conformers, str_in);
    ConfGroup* conf_group = new ConfGroup();
    bool conf_end = false;
    while (!in_conformers.eof()) {
        if (str_in.find("Dihedral") != string::npos) {
            str_in = str_in.substr(str_in.find("Dihedral")+8);
            istringstream stream_in(str_in);
            int* conf_dihedral = new int[4];
            int atom_ct = 0;
            int dih_atom = 0;
            int old_dih_atom = 0;
            stream_in >> dih_atom;  // Reading atom number
            while (dih_atom != old_dih_atom || atom_ct < 4) {
                conf_dihedral[atom_ct] = dih_atom;
                atom_ct++; old_dih_atom = dih_atom;
                stream_in >> dih_atom;
            }
            if (atom_ct == 4) conf_group->SetDih(conf_dihedral);
            else delete[] conf_dihedral;
        }
        if (str_in.find("Equilibrium") != string::npos) {
            str_in = str_in.substr(str_in.find("Equilibrium")+11);
            istringstream stream_in(str_in);
            double angle, old_angle = 1000;
            stream_in >> angle;  // Reading angle
            while (angle != old_angle) {
                conf_group->AddEqDih(angle);
                old_angle = angle;
                stream_in >> angle;
            }
        }
        getline(in_conformers, str_in);
        conf_end = (str_in == "" || in_conformers.eof() ? true : false);
        if (conf_end && conf_group->GetDih() && conf_group->GetEqDihSize() > 1) {
            // There should be 4 atoms for defining a dihedral and at least
            // two values of a dihedral angle
            // Adding group to list
            conf_list.push_back(conf_group);
            cout << "Conformers are determined by the value of "
                 << conf_group->GetDih()[0] << "-"
                 << conf_group->GetDih()[1] << "-"
                 << conf_group->GetDih()[2] << "-"
                 << conf_group->GetDih()[3] << " dihedral angle" << endl
                 << "The number of conformers is "
                 << conf_group->GetEqDihSize() << endl
                 << "Corresponding angle values are ";
            for (unsigned int j = 0; j < conf_group->GetEqDihSize(); j++) {
                cout << j+1 << ": " << conf_group->GetEqDih(j);
                if (j != conf_group->GetEqDihSize()-1) cout << ", ";
            }
            cout << endl;
            int a1i = conf_group->GetDih()[0]-1;    // First atom
            const xyz* a1 = atom_list[a1i]->GetC();
            int a2i = conf_group->GetDih()[1]-1;    // Second atom
            const xyz* a2 = atom_list[a2i]->GetC();
            int a3i = conf_group->GetDih()[2]-1;    // Third atom
            const xyz* a3 = atom_list[a3i]->GetC();
            int a4i = conf_group->GetDih()[3]-1;    // Fourth atom
            const xyz* a4 = atom_list[a4i]->GetC();
            double angle = CalcDihedral(a1, a2, a3, a4);
            chosen_conf[conf_list.size()-1] =
                GetConformerNo(conf_list[conf_list.size()-1], &angle);
            cout << "Equilibrium value of this dihedral angle is " <<
                    angle*180.0/M_PI << ", this is conformer " <<
                    chosen_conf[conf_list.size()-1]+1 << endl << endl;
            // Creating a new group
            conf_group = new ConfGroup();
        }
    }
    if (conf_group) delete conf_group;
    in_conformers.close();
    cout << "All parameters will be calculated for conformer ";
    for (unsigned int i = 0; i < conf_list.size(); i++) {
        cout << chosen_conf[i]+1;
        if (i != conf_list.size()-1) cout << "-";
    }
    cout << endl << endl; return conf_list.size();
}

// Function for reading internal coordinates
int Simulation::ReadInternals(const char* filename) {
    ifstream in_internals(filename);
    string str_in;
    if (in_internals.fail()) {
        cout << "Cannot open '" << filename << "' file with internal "
             "coordinates." << endl;
        in_internals.close(); return(0);
    }
    getline(in_internals, str_in);
    while (!in_internals.eof()) {
        istringstream stream_in(str_in);
        vector<int> internal;  // Atoms in a group which may rotate
        int internal_atom_no = 0;
        int old_internal_atom_no = 0;
        stream_in >> internal_atom_no;  // Reading atom number
        while (internal_atom_no != old_internal_atom_no) {
            internal.push_back(internal_atom_no);
            old_internal_atom_no = internal_atom_no;
            stream_in >> internal_atom_no;
        }
        if (internal.size() == 4) internals_list.push_back(internal);
        else if (internal.size() == 3)
            cout << "Currently only dihedral/torsional angles are suppoted"
                 << endl;
        getline(in_internals, str_in);
    }
    in_internals.close(); return internals_list.size();
}

// Function for reading trajectory from file
int Simulation::ReadTrajectory(const char* filename, int skip,
                               unsigned int* PI_steps) {
    energy = 0;
    int skip_ct = 1;
    cout << "Skipping first " << skip << " steps from the trajectory" << endl;
    ifstream in_tr(filename);
    if (in_tr.fail()) {
        cout << "Cannot open '" << filename << "' trajectory file." << endl;
        in_tr.close(); return(0);
    }
    string str_in;
    // Skipping equilibrium geometry at the beginning
    for (unsigned int i = 0; i < atom_list.size()+3; i++)
        getline(in_tr, str_in);
    unsigned int tr_atom_no;  // Number of atoms for the current step
    *PI_steps = 1;
    while (!in_tr.eof()) {
        istringstream stream_in_1(str_in);
        // Reading the number of atoms at the current step
        stream_in_1 >> tr_atom_no;
        if (atom_no != tr_atom_no) {
            if ((unsigned int)(round(tr_atom_no/atom_no)) ==
                tr_atom_no/atom_no) {
                if (*PI_steps == 1) {
                    *PI_steps = tr_atom_no/atom_no;
                    cout << "This is a PIMD trajectory with "
                         << *PI_steps << " beads per step" << endl;
                }
                else
                    if (*PI_steps != tr_atom_no/atom_no) {
                        cout << "PI error!" << endl; return(1);
                    }
            }
            else {
                cout << "Number of atoms in current step is not equal to the "
                     << "number of atoms at the equilibrium geometry" << endl;
                return(1);
            }
        }
        // Reading step number
        getline(in_tr, str_in);
        istringstream stream_in_2(str_in);
        string S;
        int step_no;
        double t, en = 0;
        stream_in_2 >> S >> S >> step_no >> S >> S >> S >> t
                    >> S >> S >> S >> en;
        if (skip_ct > skip) energy += en;
        // Reading coordinates
        for (unsigned int j = 0; j < *PI_steps; j++) {
            for (unsigned int i = 0; i < tr_atom_no/(*PI_steps); i++) {
                getline(in_tr, str_in);
                istringstream stream_in(str_in);
                string name;
                xyz coords;
                stream_in >> name >> coords.x >> coords.y >> coords.z;
                // Storing coordinates
                if (name == atom_list[i]->GetName()) {
                    if (skip_ct > skip)
                        atom_list[i]->AddStep(&coords);
                }
                else {
                    cout << name << "\t" << atom_list[i]->GetName()
                         << " Error" << endl;
                    return(4);
                }
            }
            // Checking for rotations
            if (atom_list[0]->GetTrajSize() != 0)
                CheckRotation(atom_list[0]->GetTrajSize()-1);
        }
        getline(in_tr, str_in); skip_ct++;
    }
    in_tr.close();
    int tr_size_this = atom_list[0]->GetTrajSize()-tr_size;
    tr_size = atom_list[0]->GetTrajSize();
    energy /= tr_size_this; return(tr_size_this);
}

unsigned long int Simulation::GetTrajectorySize() {
    return(tr_size);
}

void Simulation::ComformationalAnalysis() {
    // Conformational analysis
    if (conf_list.size() != 0) {
        for (unsigned int k = 0; k < tr_size; k++) {
            // Checking conformers
            double angle;
            conformer.push_back(new int[CONF_NO]);
            CheckConformer(k, &angle, conformer[k]);
            int angle_index = (int)(round(angle*180.0/M_PI));
            if (angle_index < 0) angle_index = -angle_index;
        }
    }
    else for (unsigned long int k = 0; k < tr_size; k++) {
        conformer.push_back(new int[1]); conformer[k][0] = 0;
    }
}

// Function for calculating average distances, amplitudes, etc.
void Simulation::CalcRaRg(unsigned int PI_steps) {
    // Creating matrices of internuclear distances, amplitudes, etc.
    dist_eq = new double[atom_no*atom_no];         // Equilibrium distances
    if (abinit)
        dist_eq_ab = new double[atom_no*atom_no];  // Ab initio distances
    dist_ra = new double[atom_no*atom_no];         // r_a distances
    dist_rg = new double[atom_no*atom_no];         // r-g distances
    u = new double[atom_no*atom_no];               // r.m.s. amplitudes
    a_M = new double[atom_no*atom_no];             // Morse constants
    kappa = new double[atom_no*atom_no];           // asymmetry constants
    // Initialising internuclear distances, amplitudes, etc.
    for (unsigned int i = 1; i <= atom_no; i++)
        for (unsigned int j = i+1; j <= atom_no; j++) {
            int index = (i-1)*atom_no+j-1;
            // Zeroing r_a, r_g, amplitudes, Morse and asymmetry constants
            dist_ra[index] = 0; dist_rg[index] = 0; u[index] = 0;
            a_M[index] = 0; kappa[index] = 0;
            // Calculating equlibrium internuclear distances
            dist_eq[index] = CalcEqDist(i-1, j-1);
            if (abinit) dist_eq_ab[index] = CalcAbInitEqDist(i-1, j-1);
        }
    cout << "Calculating average internuclear distances and amplitudes..."
         << endl << endl;
    unsigned int current_conf_ct = 0;
    for (unsigned int i = 1; i <= atom_no; i++)
        for (unsigned int j = i+1; j <= atom_no; j++) {
            int index = (i-1)*atom_no+j-1;
            // Calculating the sums
            for (unsigned int k = 0; k < PI_steps; k++) {
                long int sum_ct = 0;
                double r_ak = 0;
                double r_gk = 0;
                for (unsigned int l = 0; l < tr_size/PI_steps; l++) {
                    unsigned int idx = l*PI_steps+k;
                    if (CompareConf(idx)) {
                        r_ak += 1.0/CalcTrDist(i-1, j-1, idx);
                        r_gk += CalcTrDist(i-1, j-1, idx);
                        sum_ct++;
                    }
                }
                // Averaging
                dist_ra[index] += 1.0/(r_ak/sum_ct);
                dist_rg[index] += r_gk/sum_ct;
                current_conf_ct += sum_ct;
            }
            dist_ra[index] /= PI_steps; dist_rg[index] /= PI_steps;
        }
    current_conf_ct /= (int)(atom_no*(atom_no-1)/2.0);
    if (conf_list.size() != 0) {
        cout << "Conformer ";
        for (unsigned int i = 0; i < conf_list.size(); i++) {
            cout << chosen_conf[i]+1; if (i != conf_list.size()-1) cout << "-";
        }
        cout << " is present in " << current_conf_ct << " MD steps of "
             << tr_size << endl << endl;
    }
    // Calculating amplitudes
    for (unsigned int i = 1; i <= atom_no; i++)
        for (unsigned int j = i+1; j <= atom_no; j++) {
            int index = (i-1)*atom_no+j-1;
            for (unsigned int k = 0; k < PI_steps; k++) {
                long int sum_ct = 0;
                double uk = 0;
                double kk = 0;
                for (unsigned long int l = 0; l < tr_size/PI_steps; l++) {
                    int idx = l*PI_steps+k;
                    if (CompareConf(idx)) {
                        double x = CalcTrDist(i-1, j-1, idx);
                        double dist = x-dist_rg[index];
                        uk += pow(dist, 2);
                        kk += pow(dist, 3);
                        sum_ct++;
                    }
                }
                u[index] += sqrt(uk/sum_ct);
                kappa[index] += kk/sum_ct/6.0;
            }
            u[index] /= PI_steps; kappa[index] /= PI_steps;
            a_M[index] = kappa[index]/pow(u[index], 4)*6.0;
        }
    // Merging symmetrically equivalent internuclear distances
    for (unsigned int i = 1; i <= atom_no; i++)
        for (unsigned int j = i+1; j <= atom_no; j++) {
            // Check if the equivalent pair is in the list
            bool equivalent_pair = false;
            for (unsigned int k = 0; k < pair_list.size(); k++) {
                bool eq_ii = EquivalentAtoms(i, pair_list[k]->GetEqI(0));
                bool eq_jj = EquivalentAtoms(j, pair_list[k]->GetEqJ(0));
                bool eq_ij = EquivalentAtoms(i, pair_list[k]->GetEqJ(0));
                bool eq_ji = EquivalentAtoms(j, pair_list[k]->GetEqI(0));
                if ((eq_ii && eq_jj) || (eq_ij && eq_ji)) {
                    int index_1 = (i-1)*atom_no+j-1;
                    int index_2 = (pair_list[k]->GetEqI(0)-1)*atom_no+
                                  pair_list[k]->GetEqJ(0)-1;
                    if (EquivalentDistances(index_1, index_2))
                        equivalent_pair = true;
                }
                // If a group of equiv. atomic pairs is present, add data to it
                if (equivalent_pair) {
                    pair_list[k]->AddEq(i, j);
                    int index = (i-1)*atom_no+j-1;
                    pair_list[k]->AddRa(dist_ra[index]);
                    pair_list[k]->AddRg(dist_rg[index]);
                    pair_list[k]->AddEq(dist_eq[index]);
                    pair_list[k]->AddU(u[index]);
                    pair_list[k]->AddD(dist_ra[index]-dist_eq[index]);
                    pair_list[k]->AddA(a_M[index]);
                    pair_list[k]->AddK(kappa[index]);
                    if (abinit)
                        pair_list[k]->AddAb(dist_eq_ab[index]);
                    break;
                }
            }
            // Create new group of equivalent atomic pairs
            if (!equivalent_pair) {
                PairGroup* pair = new PairGroup(atom_list[i-1]->GetName(),
                                                atom_list[j-1]->GetName());
                pair->AddEq(i, j);
                int index = (i-1)*atom_no+j-1;
                pair->AddRa(dist_ra[index]);
                pair->AddRg(dist_rg[index]);
                pair->AddEq(dist_eq[index]);
                pair->AddU(u[index]);
                pair->AddD(dist_ra[index]-dist_eq[index]);
                pair->AddA(a_M[index]);
                pair->AddK(kappa[index]);
                if (abinit) pair->AddAb(dist_eq_ab[index]);
                pair_list.push_back(pair);
            }
        }
    // Statistics for symmetrically equivalent distances
    for (unsigned int i = 0; i < pair_list.size(); i++)
        pair_list[i]->Average();
    // Standard deviations
    for (unsigned int i = 0; i < pair_list.size(); i++)
        pair_list[i]->Statistics();
}

void Simulation::CalcProbabilities(bool print_P, bool debug) {
    cout << "Calculating probability distributions..." << endl << endl;
    for (unsigned int i = 0; i < pair_list.size(); i++) {
        if (print_P || debug) cout << "Group " << i+1 << endl << endl;
        if (print_P && debug) cout << "Looking for the boundaries..."  << endl;
        PairGroup* pair = pair_list[i];
        double min = 0;
        double max = 0;
        // This is not very efficient -
        // an additional pass to find the boundaries...
        for (unsigned int j = 0; j < pair->GetEqSize(); j++) {
            int i_index = pair->GetEqI(j)-1;
            int j_index = pair->GetEqJ(j)-1;
            for (unsigned long int k = 0; k < tr_size; k++) {
                if (CompareConf(k)) {
                    double delta = CalcTrDist(i_index, j_index, k)-
                                   CalcEqDist(i_index, j_index);
                    if (delta < min) min = delta;
                    if (delta > max) max = delta;
                }
            }
        }
        max -= max/5.0; min -= min/5.0;
        double P_step = (max-min)/(float)(P_NO-1);
        // Calculating probabilities
        if (print_P && debug) cout << "Calculating probabilities..."  << endl;
        for (unsigned int j = 0; j < pair->GetEqSize(); j++) {
            int i_index = pair->GetEqI(j)-1;
            int j_index = pair->GetEqJ(j)-1;
            for (unsigned long int k = 0; k < tr_size; k++) {
                if (CompareConf(k)) {
                    double delta = CalcTrDist(i_index, j_index, k)-
                                   CalcEqDist(i_index, j_index);
                    int index = (int)(round((delta-min)/P_step));
                    if (index >= 0 && index < P_NO)
                        pair->SetProb(index, pair->GetProb(index)+1);
                }
            }
        }
        double sum = 0;
        for (int j = 0; j < P_NO; j++) sum += pair->GetProb(j);
        for (int j = 0; j < P_NO; j++) pair->SetProb(j, pair->GetProb(j)/sum);
        // Inversion of probability distribution
        if (print_P && debug) cout << "Calculating potential..."  << endl;
        double V[P_NO];
        double r[P_NO];
        int index = (int)(round((-min)/P_step));
        for (unsigned int j = 0; j < P_NO; j++) {
            if (pair->GetProb(j) == 0) {
                if ((j != 0) & (j != P_NO-1))
                    pair->SetProb(j, (pair->GetProb(j-1)+
                                      pair->GetProb(j+1))/2.0);
                else {
                    if (j == 0) pair->SetProb(j, pair->GetProb(j+1)/2.0);
                    else pair->SetProb(j, pair->GetProb(j-1)/2.0);
                }
            }
            V[j] = -(log(pair->GetProb(j))-log(pair->GetProb(index)));
            r[j] = j*P_step+min;
        }
        // Morse potential regression (using GSL)
        if (print_P && debug)
            cout << "Morse potential regression..." << endl << endl;
        P_dist_par par;
        par.P_step = P_step; par.min = min; par.eq = pair->GetAEq();
        par.VP = V; par.rP = r;  
        const gsl_multimin_fminimizer_type *T =
              gsl_multimin_fminimizer_nmsimplex;
        gsl_multimin_fminimizer *s = NULL;
        gsl_vector *ss, *a;
        gsl_multimin_function min_func;
        size_t iter = 0;
        int status;
        double size;
        a = gsl_vector_alloc (2);
        gsl_vector_set (a, 0, 1.0); gsl_vector_set (a, 1, 1.0);
        ss = gsl_vector_alloc (2);
        gsl_vector_set (ss, 0, 1.0); gsl_vector_set (ss, 1, 1.0);
        min_func.n = 2; min_func.f = &V_dist;
        min_func.params = (void *)&par;
        s = gsl_multimin_fminimizer_alloc (T, 2);
        gsl_multimin_fminimizer_set (s, &min_func, a, ss);
        do {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);
            if (status) break;    
            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-2);
        }
        while (status == GSL_CONTINUE && iter < 1000);
        double D = gsl_vector_get(s->x, 0);
        double a_ = gsl_vector_get(s->x, 1);
        gsl_vector_free (a); gsl_vector_free(ss);
        gsl_multimin_fminimizer_free (s);
        pair->SetAA(a_); pair->SetAK(a_*pow(pair->GetAU(), 4)/6.0);
        // Printing
        if (print_P) {
            cout << setprecision(5) << fixed;
            for (int j = 0; j < P_NO; j++) {
                double a = j*P_step+min;
                cout << a+pair->GetAEq() << "\t";
                cout << pair->GetProb(j) << endl;
            }
            cout << endl;
            for (int j = 0; j < P_NO; j++) {
                double a = j*P_step+min;
                cout << a+pair->GetAEq() << "\t";
                cout << V[j] << "\t";
                cout << D*pow((1-exp(-a_*a)), 2) << endl;
            }
            cout << endl;
            cout << "Morse constant is " << pair->GetAA() << endl;
            cout << "Asymmetry constant is ";
            cout << setw(10) << setprecision(4) << scientific << pair->GetAK();
            cout << endl << endl;
        }
    }
}

// Prints all groups of symmetrically equivalent atoms
void Simulation::PrintGroups() {
    cout << "List of groups of symmetrically equivalent atomic pairs" << endl;
    cout << "i-j   r_a       r_g       r_e       u         D         "
            "a         kappa" << endl;
    cout << "pair  distance  distance  distance  ampl.     r_a-r_e   "
            "Morse c.  asym. c." << endl << endl;
    for (unsigned int i = 0; i < pair_list.size(); i++) {
        cout << "Group " << i+1 << endl;
        PairGroup* pair = pair_list[i];
        for (unsigned int j = 0; j < pair->GetEqSize(); j++) {
            cout << renumber[pair->GetEqI(j)] << "-";
            cout << renumber[pair->GetEqJ(j)];
            cout << setw(10) << setprecision(5) << fixed << pair->GetRa(j);
            cout << setw(10) << setprecision(5) << fixed << pair->GetRg(j);
            cout << setw(10) << setprecision(5) << fixed << pair->GetEq(j);
            cout << setw(10) << setprecision(5) << fixed << pair->GetU(j);
            cout << setw(10) << setprecision(5) << fixed << pair->GetD(j);
            cout << setw(10) << setprecision(5) << fixed << pair->GetA(j);
            cout << setw(12) << setprecision(3) << scientific << pair->GetK(j);
            cout << endl;
        }
        cout << "---" << endl;
        cout << "Av.";
        cout << setw(10) << setprecision(5) << fixed << pair->GetARa();
        cout << setw(10) << setprecision(5) << fixed << pair->GetARg();
        cout << setw(10) << setprecision(5) << fixed << pair->GetAEq();
        cout << setw(10) << setprecision(5) << fixed << pair->GetAU();
        cout << setw(10) << setprecision(5) << fixed << pair->GetAD();
        cout << setw(10) << setprecision(5) << fixed << pair->GetAA();
        cout << setw(12) << setprecision(3) << scientific << pair->GetAK();
        cout << endl << "S  ";
        cout << setw(10) << setprecision(5) << fixed << pair->GetSRa();
        cout << setw(10) << setprecision(5) << fixed << pair->GetSRg();
        cout << setw(10) << setprecision(5) << fixed << pair->GetSEq();
        cout << setw(10) << setprecision(5) << fixed << pair->GetSU();
        cout << setw(10) << setprecision(5) << fixed << pair->GetSD();
        cout << setw(10) << setprecision(5) << fixed << pair->GetSA();
        cout << setw(12) << setprecision(3) << scientific << pair->GetSK();
        cout << endl << endl;
    }
}

// Output of the data in the ED@ED format
void Simulation::PrintEDatED() {
    unsigned int ct = 0;
    bool* printed = new bool[pair_list.size()];
    for (unsigned int i = 0; i < pair_list.size(); i++) {
        printed[i] = false;
    }
    cout << "List of the data in the ED@ED format:" << endl << endl;
    while (ct < pair_list.size()) {
        double shortest = 1E9;
        int shortest_index = 0;
        for (unsigned int i = 0; i < pair_list.size(); i++) {
            if (pair_list[i]->GetAAbEq() < shortest && !printed[i]) {
                shortest = pair_list[i]->GetAAbEq();
                shortest_index = i;
            }
        }
        printed[shortest_index] = true; ct++;
        PairGroup* pair = pair_list[shortest_index];
        int atom_1 = renumber[pair->GetEqI(0)];
        int atom_2 = renumber[pair->GetEqJ(0)];
        if (atom_1 < atom_2)
            cout << setw(6) << left << atom_1 << setw(6) << left << atom_2;
        else cout << setw(6) << left << atom_2 << setw(6) << left << atom_1;
        cout << setw(6) << left << pair->GetEqSize();
        cout << setw(12) << setprecision(4) << left << fixed << pair->GetAU();
        cout << setw(12) << setprecision(4) << left << fixed
             << pair->GetARa()-pair->GetAEq();
        cout << setw(12) << setprecision(4) << left << fixed
             << pair->GetAA() << endl;
    }
    cout << endl; delete[] printed;
}

// Output of the data in the KCED format
void Simulation::PrintKCED() {
    unsigned int ct = 0;
    bool* printed = new bool[pair_list.size()];
    for (unsigned int i = 0; i < pair_list.size(); i++) {
        printed[i] = false;
    }
    cout << "List of the data in the KCED format:" << endl << endl;
    while (ct < pair_list.size()) {
        double shortest = 1E9;
        int shortest_index = 0;
        for (unsigned int i = 0; i < pair_list.size(); i++) {
            if (pair_list[i]->GetAAbEq() < shortest && !printed[i]) {
                shortest = pair_list[i]->GetAAbEq();
                shortest_index = i;
            }
        }
        printed[shortest_index] = true; ct++;
        PairGroup* pair = pair_list[shortest_index];
        cout << "     ";
        cout << setw(2) << right << pair->GetNameI();
        cout << setw(2) << right << pair->GetNameJ();
        int atom_1 = renumber[pair->GetEqI(0)];
        int atom_2 = renumber[pair->GetEqJ(0)];
        cout << setw(3) << right << atom_1;
        cout << setw(3) << right << atom_2;
        cout << setw(3) << right << pair->GetEqSize();
        cout << " 1.0";
        cout << setw(14) << setprecision(5) << fixed << pair->GetAAbEq();
        cout << setw(8) << setprecision(5) << fixed << pair->GetAU();
        cout << setw(10) << setprecision(1) << scientific << pair->GetAK();
        cout << setw(10) << setprecision(5) << fixed << -pair->GetAD();
        cout << "    0  1  0  0" << endl;
    }
    cout << endl; delete[] printed;
}

// Output of the data in the UNEX format
void Simulation::PrintUNEX() {
    unsigned int ct = 0;
    bool* printed = new bool[pair_list.size()];
    for (unsigned int i = 0; i < pair_list.size(); i++) {
        printed[i] = false;
    }
    cout << "List of the data in the UNEX format:" << endl << endl;
    while (ct < pair_list.size()) {
        double shortest = 1E9;
        int shortest_index = 0;
        for (unsigned int i = 0; i < pair_list.size(); i++) {
            if (pair_list[i]->GetAAbEq() < shortest && !printed[i]) {
                shortest = pair_list[i]->GetAAbEq();
                shortest_index = i;
            }
        }
        printed[shortest_index] = true; ct++;
        PairGroup* pair = pair_list[shortest_index];
        for (unsigned int j = 0; j < pair->GetEqSize(); j++) {
            int atom_1 = renumber[pair->GetEqI(j)];
            int atom_2 = renumber[pair->GetEqJ(j)];
            cout << setw(1) << left << pair->GetNameI();
            cout << setw(4) << left << atom_1;
            cout << setw(1) << left << pair->GetNameJ();
            cout << setw(4) << left << atom_2;
            cout << setw(10) << right << setprecision(5) << fixed
                 << pair->GetAAbEq();
            cout << setw(10) << right << setprecision(5) << fixed
                 << pair->GetAU();
            cout << setw(10) << right << setprecision(5) << fixed
                 << -pair->GetAD() << endl;
        }
    }
    cout << endl; delete[] printed;
}

