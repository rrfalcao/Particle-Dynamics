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
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

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
bool far_enough(double x, double y, double z, double* x_arr, double* y_arr, double* z_arr, int init,double R=0.5){ 
        for (int i = 0; i < init; i++) {
            double dx = x - x_arr[i];
            double dy = y - y_arr[i];
            double dz = z - z_arr[i];
            double dist2 = dx * dx + dy * dy + dz * dz;
            if (dist2 < R * R) {
                
                return false;  // Too close, reject
            }
        }
        
        return true;  // Far enough, accept
    } 

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

void init_particle(double* x, double* y, double* z, double* vx, double* vy, double* vz, int* type,int N, double Lx, double Ly, double Lz,double m0, double m1, double& percent_type1,std::string initial_condition) {
    srand(time(0));

    if (initial_condition == "ic-one") {
        // **Test Case 1: One stationary particle**
        x[0] = 10.0; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 0.0; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; 
        percent_type1=0.0;
    } 
    else if (initial_condition == "ic-one-vel") {
        // **Test Case 2: One moving particle**
        x[0] = 10.0; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 5.0; vy[0] = 2.0; vz[0] = 1.0;
        type[0] = 0; 
        percent_type1=0.0;
    } 
    else if (initial_condition == "ic-two") {
        // **Test Case 3: Two bouncing particles**
        x[0] = 8.5; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 0.0; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; 

        x[1] = 11.5; y[1] = 10.0; z[1] = 10.0;
        vx[1] = 0.0; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; 
        percent_type1=0.0;
    }
    else if (initial_condition == "ic-two-pass1") {
        // **Test Case 4: Two passing particles**
        x[0] = 8.5; y[0] = 11.5; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; 

        x[1] = 11.5; y[1] = 8.5; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; 
        percent_type1=0.0;
    }
    else if (initial_condition == "ic-two-pass2") {
        // **Test Case 5: Two passing particles close**
        x[0] = 8.5; y[0] = 11.3; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; 

        x[1] = 11.5; y[1] = 8.7; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; 
        percent_type1=0.0;
    }
    else if (initial_condition == "ic-two-pass3") {
        // **Test Case 6: Two passing heavy particles**
        x[0] = 8.5; y[0] = 11.3; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 1; 

        x[1] = 11.5; y[1] = 8.7; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 1; 
        percent_type1=100.0;
    }
    else if (initial_condition == "ic-random") {
        // **Random initialization for N particles**
        int inited=0;
        for (int i=0; i<N; i++){
        

            if (i < N * percent_type1 / 100.0) {
                type[i] = 1;
                
            } else {
                type[i] = 0;
                
            }

            x[i]=Lx*(double)rand()/RAND_MAX;
            y[i]=Ly*(double)rand()/RAND_MAX;
            z[i]=Lz*(double)rand()/RAND_MAX;
            vx[i]=(double)rand()/RAND_MAX - 0.5;
            vy[i]=(double)rand()/RAND_MAX - 0.5;
            vz[i]=(double)rand()/RAND_MAX - 0.5;
            
            if (far_enough(x[i], y[i], z[i], x, y, z, inited)){
                inited++;
            }
            else{
                i--;
            }

    }
    
}}
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

void compute_forces(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep,char test) {

    // Reset forces
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);

    // Lennard-Jones parameters
    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double r = 0.0;
    double eps24 =0.0;
    double sig6 = 0.0;
    int t1 = 0;
    int t2 = 0;
    double f = 0.0;

    double r6 = 0.0;
    double r2 = 0.0;
    // Loop over all pairs of particles

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            // Compute distance components
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dz = z[i] - z[j];

            // Compute squared distance
            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0) {
                
                if (test=='y'){
                    r= sqrt(r2);
                
                    if (r<min_sep){
                        min_sep=r; //For Unit Testing
                    }
                }

                t1 = type[i];
                t2 = type[j];
                
                sig6 = sigma6[t1][t2];
                eps24 = epsilon24[t1][t2];
                r6=r2*r2*r2;
                f = eps24 * sig6*(2*sig6 -r6) / (r6*r6*r2);
                 
                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;    
                fz[j] -= f * dz;
    }}}}

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
void update_positions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double dt) {
    cblas_daxpy(N, dt, vx, 1, x, 1);  // x = dt * vx + x
    cblas_daxpy(N, dt, vy, 1, y, 1);  // y = dt * vy + y
    cblas_daxpy(N, dt, vz, 1, z, 1);  // z = dt * vz + z
}
    
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
void update_velocities(double* vx, double* vy, double* vz, double* fx, double* fy, double* fz, double m0,double m1, int N0,int N1, double dt) {
    

    double factor0 = dt / m0;   // For type 0 particles (mass = 1)
    double factor1 = dt / m1;  // For type 1 particles (mass = 10)

    // Apply daxpy separately for type 1 and type 0 particles
    if (N1 > 0) {
        cblas_daxpy(N1, factor1, fx, 1, vx, 1);
        cblas_daxpy(N1, factor1, fy, 1, vy, 1);
        cblas_daxpy(N1, factor1, fz, 1, vz, 1);
    }
    
    if (N0 > 0) {
        cblas_daxpy(N0, factor0, fx + N1, 1, vx + N1, 1);
        cblas_daxpy(N0, factor0, fy + N1, 1, vy + N1, 1);
        cblas_daxpy(N0, factor0, fz + N1, 1, vz + N1, 1);
    }}

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

void apply_boundary_conditions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double Lx, double Ly, double Lz) {
    for (int i = 0; i < N; i++) {
        if (x[i] < 0) { x[i] = -x[i]; vx[i] = -vx[i]; }
        if (x[i] > Lx) { x[i] = 2 * Lx - x[i]; vx[i] = -vx[i]; }
        if (y[i] < 0) { y[i] = -y[i]; vy[i] = -vy[i]; }
        if (y[i] > Ly) { y[i] = 2 * Ly - y[i]; vy[i] = -vy[i]; }
        if (z[i] < 0) { z[i] = -z[i]; vz[i] = -vz[i]; }
        if (z[i] > Lz) { z[i] = 2 * Lz - z[i]; vz[i] = -vz[i]; }
    }
}
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
double compute_temperature(double* vx, double* vy, double* vz, double m0,double m1,int N0, int N1) {
    double kinetic_energy = 0.0;
    double boltz=0.8314459920816467;
    for (int i = 0; i < N0; i++) {  // Type 0 (mass m0)
        kinetic_energy += 0.5 * m0 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    for (int i = N0; i < (N0 + N1); i++) {  // Type 1 (mass m1)
        kinetic_energy += 0.5 * m1 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    return (2.0 / (3.0 * boltz * (N0 + N1))) * kinetic_energy;  
}
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
double compute_KE(double* vx, double* vy, double* vz, double m0,double m1, int N0, int N1) {
    double kinetic_energy = 0.0;
    
    for (int i = 0; i < N0; i++) {  // Type 0 (mass m0)
        kinetic_energy += 0.5 * m0 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    for (int i = N0; i < (N0 + N1); i++) {  // Type 1 (mass m1)
        kinetic_energy += 0.5 * m1 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    return kinetic_energy;  
}

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
void scale_velocities(double* vx, double* vy, double* vz, double m0,double m1,int N0, int N1, double T_target) {
    double T_current = compute_temperature(vx, vy, vz, m0,m1,N0, N1);

    if (T_current == 0.0) return;  // Avoid division by zero

        double scale_factor = sqrt(T_target / T_current);

    for (int i = 0; i < (N0+N1); i++) {
        vx[i] *= scale_factor;
        vy[i] *= scale_factor;
        vz[i] *= scale_factor;
    }
}

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
void unit_tests(std::string test_flag, double time, double min_separation, double* x, double* y, double* vx, double* vy, int N, double Lx, double Ly,bool end) {

    
    static bool first_collision1_detected = false;
    static bool first_collision2_detected = false;
    const double tol = 1e-4;

    if (test_flag == "ic-one") {
        // Test 1: Particle must remain stationary
        
        if (abs(vx[0]) > tol || abs(vy[0]) > tol) {
            cout << "FAIL: Particle has moved! Velocity: ("<< vx[0] << ", " << vy[0] << ")\n";
            return;
        }
        
        else if(end){cout << "PASS: Particle remained stationary.\n";}
        
    } 
    
    else if (test_flag == "ic-one-vel") {
        // Test 2: First wall collision
        if (!first_collision1_detected) {
            if (x[0] <= 0.0 || x[0] >= Lx || y[0] <= 0.0 || y[0] >= Ly) {
               
                first_collision1_detected = true;
                cout << "PASS: First wall collision at ("  << x[0] << ", " << y[0]  << ") at t = " << time << "\n";
            }
        }
    } 
    
    else if (test_flag == "ic-two") {
        // Test 3: Minimum separation check
        if(end){
            if (abs(min_separation-1.00022)<tol){
                cout << "PASS: Minimum separation: " << min_separation 
                        << " (Expected: ~1.0022)\n";
            } else {
                cout << "FAIL: Minimum separation: " << min_separation 
                        << " (Expected: ~1.0022)\n";
            }
    } }
    
    else if (test_flag == "ic-two-pass1") {
        // Test 4 & 5: Particle 1 and 2 wall collision + min separation
        if (!first_collision1_detected) {
            if (x[0] <= 0.0 || x[0] >= Lx || y[0] <= 0.0 || y[0] >= Ly) {
                first_collision1_detected = true;
                if (abs(x[0]-20.0)<(tol*100) && abs(y[0]-8.99)<(tol*100)){
                    cout << "PASS: Particle 1 hit wall at (" << x[0] << ", " << y[0] << ")\n";
                } else {
                    cout << "FAIL: Particle 1 hit wall at (" << x[0] << ", " <<y[0]  << ")\n";
                }
            }
        }
        if (!first_collision2_detected) {
            if (x[1] <= 0.0 || x[1] >= Lx || y[1] <= 0.0 || y[1] >= Ly) {
                first_collision2_detected = true;
                if (abs(x[1]-0.0)<(tol*100) && abs(y[1]-11.0)<(tol*1000)){
                    
                    cout << "PASS: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                } else {
        
                    cout << "FAIL: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                }
            }
        }
        if(end){
        if (abs(min_separation-2.89619)<tol){
            cout << "PASS: Minimum separation: " << min_separation 
                        << " (Expected: ~2.89619)\n";
            } else {
            cout << "FAIL: Minimum separation: " << min_separation 
                        << " (Expected: ~2.89619)\n";
            }
            } 
    }
    else if (test_flag == "ic-two-pass2") {
        // Test 5: Particle 1 and 2 wall collision + min separation
        if (!first_collision1_detected) {
            ;
            if (x[0] <= 0.0 || x[0] >= Lx || y[0] <= 0.0 || y[0] >= Ly) {
                first_collision1_detected = true;
                
                if (abs(x[0]-20.0)<(tol*100) && abs(y[0]-15.33)<(tol*100)){
                    cout << "PASS: Particle 1 hit wall at (" << x[0] << ", " << y[0] << ")\n";
                } else {
                    cout << "FAIL: Particle 1 hit wall at (" << x[0] << ", " <<y[0]  << ")\n";
                }
            }
        }
        if (!first_collision2_detected) {
            if (x[1] <= 0.0 || x[1] >= Lx || y[1] <= 0.0 || y[1] >= Ly) {
                first_collision2_detected = true;
                if (abs(x[1]-0.0)<(tol*100) && abs(y[1]-4.67)<(tol*100)){
                    cout << "PASS: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                } else {
                    cout << "FAIL: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                }
            }
        }
        if(end){
        if (abs(min_separation-1.02368)<tol){
                cout << "PASS: Minimum separation: " << min_separation 
                            << " (Expected: ~1.02368)\n";
            } else {
                cout << "FAIL: Minimum separation: " << min_separation 
                            << " (Expected: ~1.02368)\n";
            }
            } 
    }
    
    else if (test_flag == "ic-two-pass3") {
        // Test 6: Minimum separation only
        if(end){
            if (abs(min_separation-3.10166)<tol){
                cout << "PASS: Minimum separation: " << min_separation 
                        << " (Expected: ~3.10166)\n";
            } else {
                cout << "FAIL: Minimum separation: " << min_separation 
                        << " (Expected: ~3.10166)\n";
            }
        }
}
}


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
    int* type = new int[N];
    double m0 = 1.0;
    double m1 = 10.0;
    char test = 'y';
    if (initial_condition == "ic-random") {
        test='n';  }
    if (test == 'y') {
        cout<<initial_condition<<endl;    }

    init_particle(x, y, z, vx, vy, vz, type, N, Lx, Ly, Lz, m0, m1, percent_type1,initial_condition);
    int N1 = static_cast<int>(N * percent_type1 / 100.0); // Number of type 1 particles
    int N0 = N - N1; // Number of type 0 particles
    
    std::ofstream kinetic_file("kinetic_energy.txt", std::ofstream::trunc);


    // Create text files in overwrite mode
    
    if (temperature > 0.0) {
        scale_velocities(vx, vy, vz, m0, m1, N0, N1, temperature);
    }
    /////////////////////// Numerical Loop /////////////////////////
    int steps = T_tot / dt;

    double time=0.0;
    int writestep=static_cast<int>(0.1 / dt); // Cast double division to int

    for (int t = 0; t < steps; t++) {
        
        if (t%writestep==0) {  // Write data every 0.1 time units
            double K = compute_KE(vx, vy, vz, m0, m1, N0, N1);
            kinetic_file << time << " " << K << "\n";
        }
        compute_forces(x, y, z, fx, fy, fz, type, N, min_sep,test);
        update_velocities(vx, vy, vz, fx, fy, fz, m0, m1, N0, N1, dt);
        // Temperature Change - only if temp is set
        if (temperature > 0.0) {
            scale_velocities(vx, vy, vz, m0, m1, N0, N1, temperature);
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

