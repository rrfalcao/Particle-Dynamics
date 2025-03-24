/**
 * @file mdmP.cpp
 * @brief Molecular Dynamics Simulation Code
 * 
 * This code simulates molecular dynamics using Lennard-Jones interactions between particles.
 * It includes functions for initializing particle positions and velocities, computing forces,
 * updating positions and velocities, applying boundary conditions, and computing system properties.
 * 
 * @author Ridge
 * @date 2025-03-13
 */

#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include <limits>
#include <fstream> 
#include <boost/program_options.hpp>

#include "init_testing.h"
#include "compute_parallel.h"
namespace po = boost::program_options;
using namespace std;



/**
 * @brief Computes forces on particles using the Lennard-Jones potential.
 * 
 * @param x X-coordinates array.
 * @param y Y-coordinates array.
 * @param z Z-coordinates array.
 * @param fx Force components in X direction.
 * @param fy Force components in Y direction.
 * @param fz Force components in Z direction.
 * @param type Particle type array.
 * @param N Total number of particles.
 * @param min_sep Minimum separation distance (used for unit testing).
 * @param test Flag for unit testing.
 */

 





/**
 * @brief Updates particle positions using velocity integration.
 *
 * @param x X-coordinates of particles.
 * @param y Y-coordinates of particles.
 * @param z Z-coordinates of particles.
 * @param vx X-velocity components.
 * @param vy Y-velocity components.
 * @param vz Z-velocity components.
 * @param N Number of particles.
 * @param dt Time step size.
 */

    
/**
 * @brief Updates particle velocities using computed forces and mass properties.
 * @param vx X-velocity array.
 * @param vy Y-velocity array.
 * @param vz Z-velocity array.
 * @param fx X-force array.
 * @param fy Y-force array.
 * @param fz Z-force array.
 * @param m Mass of particles.
 * @param N Number of particles.
 * @param dt Time step.
 */

/**
 * @brief Applies boundary conditions by reflecting particles off walls.
 * @param x Array of x-coordinates.
 * @param y Array of y-coordinates.
 * @param z Array of z-coordinates.
 * @param vx Array of x-velocity components.
 * @param vy Array of y-velocity components.
 * @param vz Array of z-velocity components.
 * @param N Number of particles.
 * @param Lx Box length in x-direction.
 * @param Ly Box length in y-direction.
 * @param Lz Box length in z-direction.
 */


/**
 * @brief Computes the system temperature based on kinetic energy.
 * @param vx X-velocity array.
 * @param vy Y-velocity array.
 * @param vz Z-velocity array.
 * @param m Mass of particles.
 * @param N Number of particles.
 * @return Computed temperature.
 */

/**
 * @brief Computes the system temperature based on kinetic energy.
 * @param vx X-velocity array.
 * @param vy Y-velocity array.
 * @param vz Z-velocity array.
 * @param m Mass of particles.
 * @param N Number of particles.
 * @return Computed temperature.
 */

/**
 * @brief Scales velocities to achieve a target temperature.
 * @param vx X-velocity array.
 * @param vy Y-velocity array.
 * @param vz Z-velocity array.
 * @param m Mass of particles.
 * @param N Number of particles.
 * @param T_target Target temperature.
 */


/**
 * @brief Performs unit tests for verifying particle behavior.
 *
 * This function checks the correctness of simulation outcomes based on 
 * predefined test conditions. It validates parameters such as velocity,
 * particle collisions, and minimum separation distances.
 *
 * @param test_flag The name of the test case being executed.
 * @param time The current simulation time.
 * @param min_separation The minimum separation observed between particles.
 * @param x Array of x positions of particles.
 * @param y Array of y positions of particles.
 * @param vx Array of x velocities of particles.
 * @param vy Array of y velocities of particles.
 * @param N Number of particles.
 * @param Lx Length of the simulation box in x direction.
 * @param Ly Length of the simulation box in y direction.
 * @param end Boolean flag indicating if this is the final test run.
 */

/**
 * @brief Main function for running the molecular dynamics simulation.
 */

int main(int argc, char** argv) {

    double start_time, end_time;
    
    start_time = omp_get_wtime();  // Start timer

    double Lx = 20.0, Ly = 20.0, Lz = 20.0;
    double dt = 0.001, T_tot = 50.0;
    int N = 8;
    double percent_type1 = 10.0;
    string initial_condition;
    bool end = false;
    double temperature = 0.0;
    

    ////////////////////////// Command Line Configuration ///////////////////////////////////
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "Print available options")
        ("Lx", po::value<double>(&Lx)->default_value(20.0), "x length (Angstroms)")
        ("Ly", po::value<double>(&Ly)->default_value(20.0), "y length (Angstroms)")
        ("Lz", po::value<double>(&Lz)->default_value(20.0), "z length (Angstroms)")
        ("dt", po::value<double>(&dt)->default_value(0.001), "Time step")
        ("T", po::value<double>(&T_tot)->default_value(50.0), "Final time")
        ("N", po::value<int>(&N)->default_value(8), "Number of particles")
        ("percent-type1", po::value<double>(&percent_type1)->default_value(10.0), "Percentage of type 1 particles")
        ("ic-one", "Initial condition: one stationary particle")
        ("ic-one-vel", "Initial condition: one moving particle")
        ("ic-two", "Initial condition: two bouncing particles")
        ("ic-two-pass1", "Initial condition: two passing particles")
        ("ic-two-pass2", "Initial condition: two passing particles close")
        ("ic-two-pass3", "Initial condition: two passing particles close, heavy")
        ("ic-random", "Initial condition: N random particles")
        ("temp", po::value<double>(), "Temperature (optional)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return 0;
    }

    
    if (vm.count("ic-one")) { initial_condition = "ic-one"; N = 1; }
    if (vm.count("ic-one-vel")) { initial_condition = "ic-one-vel"; N = 1; }
    if (vm.count("ic-two")) { initial_condition = "ic-two"; N = 2; }
    if (vm.count("ic-two-pass1")) { initial_condition = "ic-two-pass1"; N = 2; }
    if (vm.count("ic-two-pass2")) { initial_condition = "ic-two-pass2"; N = 2; }
    if (vm.count("ic-two-pass3")) { initial_condition = "ic-two-pass3"; N = 2; }
    if (vm.count("ic-random")) {
        initial_condition = "ic-random";
        if (!vm.count("N") || !vm.count("percent-type1")) {
            std::cerr << "Error: --N and --percent-type1 required for --ic-random.\n";
            return 1;
        }
    }

    if (initial_condition.empty()) {
        std::cerr << "Error: You must specify one --ic-* option.\n";
        return 1;
    }

    if (vm.count("temp")) {
        temperature = vm["temp"].as<double>();
    }
    //////////////////////////////////////////////////////////////////////


    // **Allocate Memory for Particles**
    
    double* x = new double[N]; 
    double* y = new double[N]; 
    double* z = new double[N];
    double* vx = new double[N]; 
    double* vy = new double[N]; 
    double* vz = new double[N];
    double* fx = new double[N]; 
    double* fy = new double[N]; 
    double* fz = new double[N];
    double* m = new double[N]; 
    int* type = new int[N];
 
    double m0 = 1.0;
    double m1 = 10.0;
    double min_sep = Lx;
    char test = 'y';
    if (initial_condition == "ic-random") {
        test='n';  }
    if (test == 'y') {
        cout<<initial_condition<<endl;    }

    init_particle(x, y, z, vx, vy, vz, type, m, N, Lx, Ly, Lz, m0, m1,percent_type1, initial_condition);
    std::ofstream kinetic_file("kinetic_energy.txt", std::ofstream::trunc);
    
    /////////////////////// Numerical Loop /////////////////////////
    int steps = T_tot / dt;
    double time=0.0;
    int writestep=static_cast<int>(0.1 / dt); // Cast double division to int
    for (int t = 0; t < steps; t++) {
        
        if (t%writestep==0) {  // Write data every 0.1 time units
            double K = compute_KE(vx, vy, vz, m, N);
            kinetic_file << time << " " << K << "\n";
        }
        if (test=='y') {
            compute_forces_test(x, y, z, fx, fy, fz, type, N, min_sep,test);
            
        }
        else{
            compute_forces(x, y, z, fx, fy, fz, type, N);
        }
        update_velocities(vx, vy, vz, fx, fy, fz, m, N, dt);
        // Temperature Change - only if temp is set
        if (temperature > 0.0) {
            scale_velocities(vx, vy, vz, m, N, temperature);
        }
        update_positions(x, y, z, vx, vy, vz, N, dt);
        
        if (test=='y') {
            unit_tests(initial_condition, time ,min_sep, x, y, vx, vy, N, Lx, Ly,end);
        }

        apply_boundary_conditions(x, y, z, vx, vy, vz, N, Lx, Ly, Lz);

        time += dt;
    }

    end=true;
    if (test=='y') {
        unit_tests(initial_condition, time, min_sep, x, y, vx, vy, N, Lx, Ly,end);
    }
    /////// Cleanup Section ///////////
    delete[] x; delete[] y; delete[] z;
    delete[] vx; delete[] vy; delete[] vz;
    delete[] fx; delete[] fy; delete[] fz;
    delete[] type; delete[] m;
    kinetic_file.close();
    end_time = omp_get_wtime();  // End timer
    
    double elapsed_time = end_time - start_time;
    cout << "Execution time: " << elapsed_time << " seconds\n";
    ///////////////////////////////////
    
}

