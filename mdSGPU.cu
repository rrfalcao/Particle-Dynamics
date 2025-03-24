/**
 * @file mdSGPU.cpp
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
#include <boost/program_options.hpp>
#include <cuda_runtime.h>
#include <random>
#include "init_testing.h"
namespace po = boost::program_options;
using namespace std;
#define BLOCK_SIZE 32
/**
 * @brief Checks if a new particle is far enough from existing particles before placing it.
 * 
 * @param x X-coordinate of the new particle.
 * @param y Y-coordinate of the new particle.
 * @param z Z-coordinate of the new particle.
 * @param x_arr Array of existing particle X-coordinates.
 * @param y_arr Array of existing particle Y-coordinates.
 * @param z_arr Array of existing particle Z-coordinates.
 * @param init Number of already placed particles.
 * @param R Minimum allowed separation distance (default 0.5).
 * @return True if the new particle is far enough, False otherwise.
 */


/**
 * @brief Initializes particle positions and velocities based on a chosen initial condition.
 * 
 * @param x X-coordinates array.
 * @param y Y-coordinates array.
 * @param z Z-coordinates array.
 * @param vx X-velocities array.
 * @param vy Y-velocities array.
 * @param vz Z-velocities array.
 * @param type Array to store particle types.
 * @param N Total number of particles.
 * @param Lx Simulation box length in X direction.
 * @param Ly Simulation box length in Y direction.
 * @param Lz Simulation box length in Z direction.
 * @param m0 Mass of type 0 particles.
 * @param m1 Mass of type 1 particles.
 * @param percent_type1 Percentage of type 1 particles.
 * @param initial_condition Selected initial condition.
 */

 
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

// Lennard-Jones constants


/**
 * @brief CUDA Kernel to compute Lennard-Jones forces in parallel.
 */

/**
 * @brief Host function to call the GPU kernel.
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
 * @param m0 Mass of type 0 particles.
 * @param m1 Mass of type 1 particles.
 * @param N0 Number of type 0 particles.
 * @param N1 Number of type 1 particles.
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
 * @param m0 Mass of type 0 particles.
 * @param m1 Mass of type 1 particles.
 * @param N0 Number of type 0 particles.
 * @param N1 Number of type 1 particles.
 * @return Computed temperature.
 */

/**
 * @brief Computes the system temperature based on kinetic energy.
 * @param vx X-velocity array.
 * @param vy Y-velocity array.
 * @param vz Z-velocity array.
 * @param m0 Mass of type 0 particles.
 * @param m1 Mass of type 1 particles.
 * @param N0 Number of type 0 particles.
 * @param N1 Number of type 1 particles.
 * @return Computed temperature.
 */


/**
 * @brief Scales velocities to achieve a target temperature.
 * @param vx X-velocity array.
 * @param vy Y-velocity array.
 * @param vz Z-velocity array.
 * @param m0 Mass of type 0 particles.
 * @param m1 Mass of type 1 particles.
 * @param N0 Number of type 0 particles.
 * @param N1 Number of type 1 particles.
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
    auto start_time = std::chrono::high_resolution_clock::now();
    double Lx = 20.0, Ly = 20.0, Lz = 20.0;
    double dt = 0.001, T_tot = 50.0;
    int N = 8;
    double percent_type1 = 10.0;
    string initial_condition;
    bool end=true;
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
    double min_sep = Lx;
    
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


    //////////////////////// Memory allocation /////////////////////////
    
    double* h_x = new double[N];
    double* h_y = new double[N];
    double* h_z = new double[N];
    double* h_vx = new double[N];
    double* h_vy = new double[N];
    double* h_vz = new double[N];
    double* h_m  = new double[N];
    int* h_type = new int[N];
    double m0 = 1.0;
    double m1 = 10.0;

    double *d_x, *d_y, *d_z, *d_vx, *d_vy, *d_vz, *d_fx, *d_fy, *d_fz, *d_m;
    int* d_type;

    cudaMalloc(&d_x, N * sizeof(double));
    cudaMalloc(&d_y, N * sizeof(double));
    cudaMalloc(&d_z, N * sizeof(double));
    cudaMalloc(&d_vx, N * sizeof(double));
    cudaMalloc(&d_vy, N * sizeof(double));
    cudaMalloc(&d_vz, N * sizeof(double));
    cudaMalloc(&d_fx, N * sizeof(double));
    cudaMalloc(&d_fy, N * sizeof(double));
    cudaMalloc(&d_fz, N * sizeof(double));
    cudaMalloc(&d_m,  N * sizeof(double));
    cudaMalloc(&d_type, N * sizeof(int));

    ///////////////////////////////////////////////////
    bool test = (initial_condition != "ic-random");

    if (test) {
        cout<<initial_condition<<endl;    }

    
    

    ////////// Initislise on CPU and copy to GPU /////////////////////////
    init_particle(h_x, h_y, h_z, h_vx, h_vy, h_vz, h_type, h_m, N, Lx, Ly, Lz, m0, m1, percent_type1, initial_condition);
    cudaMemcpy(d_x, h_x, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, h_y, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_z, h_z, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vx, h_vx, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vy, h_vy, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vz, h_vz, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_m, h_m, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_type, h_type, N * sizeof(int), cudaMemcpyHostToDevice);
    /////////////////////////////////////////////////////////////////////
    
    std::ofstream kinetic_file("kinetic_energy.txt", std::ofstream::trunc);

    
    
    if (temperature > 0.0) {
        scale_velocities(d_vx, d_vy, d_vz, d_m, N, temperature);

    }
    /////////////////////// Numerical Loop /////////////////////////
    int steps = T_tot / dt;

    double time=0.0;
    int writestep=static_cast<int>(0.1 / dt); // Cast double division to int

    for (int t = 0; t < steps; t++) {
        
        if (t%writestep==0) {  // Write data every 0.1 time units
            double K = compute_KE(d_vx, d_vy, d_vz, d_m, N);
            kinetic_file << time << " " << K << "\n";
        }
        compute_forces(d_x, d_y, d_z, d_fx, d_fy, d_fz, d_type, N);

        update_velocities(d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, d_m, N, dt);
        // Temperature Change - only if temp is set
        if (temperature > 0.0) {
            scale_velocities(d_vx, d_vy, d_vz, d_m, N, temperature);
        }
        update_positions(d_x, d_y, d_z, d_vx, d_vy, d_vz, N, dt);

        apply_boundary_conditions(d_x, d_y, d_z, d_vx, d_vy, d_vz, N, Lx, Ly, Lz);
        
        if (test) {
            cudaMemcpy(h_x, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_y, d_y, N * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_vx, d_vx, N * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_vy, d_vy, N * sizeof(double), cudaMemcpyDeviceToHost);
            double r;
            for (int i = 0; i < N; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    double dx = h_x[i] - h_x[j];
                    double dy = h_y[i] - h_y[j];
                    double dz = h_z[i] - h_z[j];
                    r = sqrt(dx * dx + dy * dy + dz * dz);
                    if (r < min_sep) {
                        min_sep = r;
                    }
                }
            }
            
            unit_tests(initial_condition, time, min_sep, h_x, h_y, h_vx, h_vy, N, Lx, Ly,end);
        }

        


        time += dt;
    }
    end=false;
    if (test) {
        unit_tests(initial_condition, time, min_sep, h_x, h_y, h_vx, h_vy, N, Lx, Ly,end);
    }
    
    

    /////// Cleanup Section ///////////
    cudaFree(d_x); cudaFree(d_y); cudaFree(d_z);
    cudaFree(d_vx); cudaFree(d_vy); cudaFree(d_vz);
    cudaFree(d_fx); cudaFree(d_fy); cudaFree(d_fz);
    cudaFree(d_m); cudaFree(d_type);
    delete[] h_x; delete[] h_y; delete[] h_z;
    delete[] h_vx; delete[] h_vy; delete[] h_vz;
    delete[] h_m; delete[] h_type;
    kinetic_file.close();
    auto end_time = std::chrono::high_resolution_clock::now();
    double time_new = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    cout<<"Time taken: "<<time_new/1000<<" s"<<endl;
    ///////////////////////////////////
    
}

