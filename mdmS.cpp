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
#include <chrono>
#include <iostream>
// #include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include <limits>
#include <fstream> 
#include <random>
#include "compute_serial.h"
#include "init_testing.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

/**
 * @brief Main function for running the molecular dynamics simulation.
 */

int main(int argc, char** argv) {
    auto start_time = std::chrono::high_resolution_clock::now();
    double Lx = 20.0, Ly = 20.0, Lz = 20.0;
    double dt = 0.001, T_tot = 50.0;
    int N = 8;
    double percent_type1 = 10.0;
    string initial_condition;
    bool end = false;
    double temperature = 0.0;
    double min_sep = Lx;

    ////////////////////////// Command Line Configuration ///////////////////////////////////
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "Print available options")
        ("Lx", po::value<double>(&Lx)->default_value(20.0), "x length (Angstroms)")
        ("Ly", po::value<double>(&Ly)->default_value(20.0), "y length (Angstroms)")
        ("Lz", po::value<double>(&Lz)->default_value(20.0), "z length (Angstroms)")
        ("dt", po::value<double>(&dt)->default_value(0.001), "Time step")
        ("T", po::value<double>(&T_tot)->default_value(40.0), "Final time")
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
    int* type = new int[N];
    double* m = new double[N]; 

    double m0 = 1.0;
    double m1 = 10.0;
    char test = 'y';
    if (initial_condition == "ic-random") {
        test='n';  }
    if (test == 'y') {
        cout<<initial_condition<<endl;    }

    init_particle(x, y, z, vx, vy, vz, type, m, N, Lx, Ly, Lz, m0, m1, percent_type1,initial_condition);    
    std::ofstream kinetic_file("kinetic_energy.txt", std::ofstream::trunc);
    std::ofstream particle_file("particles.txt", std::ofstream::trunc);

        
    // Create text files in overwrite mode
    
    if (temperature > 0.0) {
        scale_velocities(vx, vy, vz, m, N, temperature);
    }
    /////////////////////// Numerical Loop /////////////////////////
    int steps = T_tot / dt;

    double time=0.0;
    int writestep=static_cast<int>(0.1 / dt); // Cast double division to int

    for (int t = 0; t < steps; t++) {
        
        if (t%writestep==0) {  // Write data every 0.1 time units
            double K = compute_KE(vx, vy, vz, m, N);
            kinetic_file << time << " " << K << "\n";
            if(test=='y'){
                for (int i = 0; i < N; i++) {
                    particle_file << time << " " << i << " " << type[i] << " "
                                  << x[i] << " " << y[i] << " " << z[i] << " "
                                  << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
                }
            }
        }
        compute_forces(x, y, z, fx, fy, fz, type, N, min_sep,test);
        update_velocities(vx, vy, vz, fx, fy, fz, m, N, dt);
        // Temperature Change - only if temp is set
        if (temperature > 0.0) {
            scale_velocities(vx, vy, vz, m, N, temperature);
        }
        update_positions(x, y, z, vx, vy, vz, N, dt);
        
        
        if (test=='y') {
            unit_tests(initial_condition, time, min_sep, x, y, vx, vy, N, Lx, Ly,end);
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
    delete[] type;
    auto end_time = std::chrono::high_resolution_clock::now();
    double time_new = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    cout<<"Time taken: "<<time_new/1000<<" s"<<endl;
    ///////////////////////////////////
    
}

