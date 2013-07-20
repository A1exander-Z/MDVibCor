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

// GSL error handler function
void new_gsl_error_handler(const char* reason, const char* file, int line,
                           int gsl_errno) {
    cout << "GSL error:" << endl << reason << endl << endl;
}

// Helper function informing about required parameter
void ReqParam(const char c) {
    cout << "Error: -" << c << " option requires a parameter" << endl;
}

// Function printing help/usage message
void PrintUsage() {
    cout << "Usage: mdvibcor [options] <trajectory_1> [trajectory_2] ...\n"
    " trajectory - one or several files containing trajectories "
    "from\n              MD simulations by CP2K program\n"
    "Options:\n"
    " -r renumbering - file with a list of atoms to be renumbered\n"
    " -e equivalents - file with a list of symmetrically equivalent atoms\n"
    " -R rotations - file with a list of subunits which can rotate\n"
    " -c conformers - file describing different conformers\n"
    " -a ab_initio - XYZ file containing ab initio equilibrium geometry\n"
    " -i internals - file describing internal coordinates for which\n"
    "    distributions will be computed\n"
    " -s skip - number of steps in the beginning of MD simulation"
    " to be skipped\n    (default: 2000)\n"
    " -t tolerance - maximum allowed difference between "
    "symmetrically equivalent\n    distances, in % (default: 0.1%)\n"
    " -P - print probability distributions\n"
    " -d - print additional output for debugging\n"
    "See the program documentation for descriptions of file formats\n";
}

// Main program starts here
int main(int argc, const char* argv[]) {
    clock_t st_time = clock();
    cout << "MDVibCor version 0.8.0, copyright (C) 2008-2013 "
         << "Alexander Zakharov" << endl
         << "Multithreaded version, running on a CPU with "
         << boost::thread::hardware_concurrency() << " cores";
    cout << endl << endl;
    double tolerance = 0.1; // Default maximum allowed difference between
                            // symmetrically equivalent distances (in %)
    int skip = 2000;        // Default number of initial MD steps to skip
    bool debug = false;     // Should I print additional output for debugging?
    bool print_P = false;   // Should I print probability distributions?
    // Starting command line parsing
    bool bad_command = false;
    int traj_arg = 0;
    int renumber_arg = 0;
    int equiv_arg = 0;
    int rot_arg = 0;
    int conf_arg = 0;
    int abinit_arg = 0;
    int internals_arg = 0;
    int argv_no = 1;
    while (argv_no < argc) {
        if (argv[argv_no][0] == '-' && argv_no != argc-1 &&
            strlen(argv[argv_no]) == 2) {  // This is a command line option
            switch (argv[argv_no][1]) {
                case 'r':  // Renumbering file
                    renumber_arg = argv_no+1; argv_no++;
                    if (argv[renumber_arg][0] == '-') {
                        ReqParam('r'); argv_no--; bad_command = true;
                    }
                    break;
                case 'e':  // Equivalent atoms file
                    equiv_arg = argv_no+1; argv_no++;
                    if (argv[equiv_arg][0] == '-') {
                        ReqParam('e'); argv_no--; bad_command = true;
                    }
                    break;
                case 'R':  // Rotating atoms file
                    rot_arg = argv_no+1; argv_no++;
                    if (argv[rot_arg][0] == '-') {
                        ReqParam('R'); argv_no--; bad_command = true;
                    }
                    break;
                case 'c':  // Conformers file
                    conf_arg = argv_no+1; argv_no++;
                    if (argv[conf_arg][0] == '-') {
                        ReqParam('c'); argv_no--; bad_command = true;
                    }
                    break;
                case 'a':  // Ab initio geometry file
                    abinit_arg = argv_no+1; argv_no++;
                    if (argv[conf_arg][0] == '-') {
                        ReqParam('a'); argv_no--; bad_command = true;
                    }
                    break;
                case 's': {  // Skipping MD steps
                    string skip_arg = argv[argv_no+1]; argv_no++;
                    if (skip_arg[0] == '-') {
                        ReqParam('s'); argv_no--; bad_command = true;
                    }
                    else {
                        istringstream skip_arg_str(skip_arg);
                        skip_arg_str >> skip;
                    }
                    break;
                }
                case 't': {  // Tolerance
                    string tol_arg = argv[argv_no+1]; argv_no++;
                    if (tol_arg[0] == '-') {
                        ReqParam('t'); argv_no--; bad_command = true;
                    }
                    else {
                        istringstream tol_arg_str(tol_arg);
                        tol_arg_str >> tolerance;
                    }
                    break;
                }
                case 'i': {  // Internal coordinates
                    internals_arg = argv_no+1; argv_no++;
                    if (argv[internals_arg][0] == '-') {
                        ReqParam('i'); argv_no--; bad_command = true;
                    }
                    break;
                }
                case 'd':  // Debug
                    debug = true; break;
                case 'P':  // Compute probabilities
                    print_P = true; break;
            }
            traj_arg = argv_no+1;
        }
        if ((argv[argv_no][0] == '-' && argv_no == argc-1) ||
           (argv[argv_no][0] == '-' && strlen(argv[argv_no]) != 2))
            bad_command = true;
        argv_no++;
    }
    if (traj_arg >= argc) bad_command = true;
    if (traj_arg == 0 && argc == 1) bad_command = true;
    if (traj_arg == 0 && argc == 2) traj_arg++;
    if (bad_command) {
        PrintUsage(); return(0);
    }
    // Starting to read and process parameters, trajectory files, etc.
    Simulation sim;
    sim.SetTolerance(tolerance);
    // Use two-digit exponent format under Windows
    #ifdef _WIN32
        int output_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    #endif
    // Opening a trajectory file (containing equilibrium geometry as the first
    // entry) and reading equilibrium geometry
    cout << "Reading equilibrium geometry from '" << argv[traj_arg] << "'..."
         << endl;
    int atom_no;
    atom_no = sim.ReadEqGeometry(argv[traj_arg]);
    if (atom_no == 0) {
        cout << "Canont read equilibrium geometry from a trajectory, "
             << "unable to proceed." << endl;
        return(0);
    }
    if (debug) sim.PrintCoordinates();
    cout << atom_no << " atomic coordinates have been read" << endl << endl;
    // Reading ab initio equilibrium geometry (if present)
    if (abinit_arg != 0) {
        cout << "Reading ab initio equilibrium geometry from '"
             << argv[abinit_arg] << "'..." << endl << endl;
        atom_no = sim.ReadAbInitioGeometry(argv[abinit_arg]);
        if (atom_no != 0 && debug) sim.PrintAbInitioCoordinates();
    }
    // Reading atom renumbering table (if any)
    if (renumber_arg != 0) {
        cout << "Reading atom renumbering table from '"
             << argv[renumber_arg] << "'..." << endl << endl;
        sim.ReadRenumbering(argv[renumber_arg]);
        if (debug) sim.PrintRenumbering();
    }
    // Reading list of groups of equivalent atoms (if any)
    if (equiv_arg != 0) {
        sim.ReadEquivalent(argv[equiv_arg]);
        cout << "Reading list of groups of equivalent atoms from '"
             << argv[equiv_arg] << "'..." << endl << endl;
        if (debug) sim.PrintEquivalent();
    }
    // Reading list of groups of atoms which may exchange places
    // (rotating subunits), if any
    if (rot_arg != 0) {
        cout << "Reading list of rotating atoms groups from '"
             << argv[rot_arg] << "'..." << endl;
        int result = sim.ReadRotations(argv[rot_arg]);
        cout << "Data about " << result
             << " rotating groups of atoms have been read" << endl << endl;
        if (debug) sim.PrintRotations();
    }
    // Reading list of conformers. Conformer of the equilibrium geometry
    // will be used to calculate all parameters
    if (conf_arg != 0) {
        cout << "Reading list of conformers from '"
             << argv[conf_arg] << "'..." << endl << endl;
        sim.ReadConformers(argv[conf_arg]);
    }
    // Reading list of internal coordinates
    if (internals_arg != 0) {
        cout << "Reading list of internal coordinates from '"
             << argv[internals_arg] << "'..." << endl << endl;
        sim.ReadInternals(argv[internals_arg]);
    }
    cout << "Elapsed time: " << (clock()-st_time)/CLOCKS_PER_SEC
         << " s" << endl;
    // Reading trajectory (trajectories)
    unsigned int PI_steps = 0;
    while (traj_arg < argc) {
        cout << "Reading trajectory from '" << argv[traj_arg]
             << "'..." << endl;
        int result = sim.ReadTrajectory(argv[traj_arg], skip, &PI_steps);
        if (result == 0) {
            cout << "Canont read trajectory, unable to proceed." << endl;
            return(0);
        }
        cout << "Read " << result/PI_steps
             << " trajectory steps" << endl << endl;
        traj_arg++;
    }
    cout << "Elapsed time: " << (clock()-st_time)/CLOCKS_PER_SEC
         << " s" << endl;
    int tr_size = sim.GetTrajectorySize();
    if (PI_steps == 1)
        cout << "Total number of MD steps is " << tr_size << endl << endl;
    else
        cout << "Total number of MD steps is " << tr_size/PI_steps
             << " (" << PI_steps << " beads per step)" << endl << endl;
    // Doing conformational analysis (this takes time)
    sim.ComformationalAnalysis();
    cout << "Elapsed time: " << (clock()-st_time)/CLOCKS_PER_SEC
         << " s" << endl;
    // Computing thermal-average distances, etc.
    sim.Statistics(PI_steps);
    cout << "Elapsed time: " << (clock()-st_time)/CLOCKS_PER_SEC
         << " s" << endl;
    // Printing the results
    sim.PrintGroups();
    gsl_set_error_handler (&new_gsl_error_handler);
    if (PI_steps == 1) {
        sim.CalcProbabilities(print_P, debug);
    }
    sim.PrintKCED(); sim.PrintEDatED(); sim.PrintUNEX();
    cout << "Elapsed time: " << (clock()-st_time)/CLOCKS_PER_SEC
         << " s" << endl;
}

