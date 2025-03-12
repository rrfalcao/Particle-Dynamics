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
#define F77NAME(x) x##_

extern "C" {

    }

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

void init_particle(double* x, double* y, double* z, double* vx, double* vy, double* vz, int* type, double* m,int N, double Lx, double Ly, double Lz,double m0, double m1, double percent_type1,std::string initial_condition) {
    srand(time(0));

    if (initial_condition == "ic-one") {
        // **Test Case 1: One stationary particle**
        x[0] = 10.0; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 0.0; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;
    } 
    else if (initial_condition == "ic-one-vel") {
        // **Test Case 2: One moving particle**
        x[0] = 10.0; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 5.0; vy[0] = 2.0; vz[0] = 1.0;
        type[0] = 0; m[0] = m0;
    } 
    else if (initial_condition == "ic-two") {
        // **Test Case 3: Two bouncing particles**
        x[0] = 8.5; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 0.0; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;

        x[1] = 11.5; y[1] = 10.0; z[1] = 10.0;
        vx[1] = 0.0; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; m[1] = m0;
    }
    else if (initial_condition == "ic-two-pass1") {
        // **Test Case 4: Two passing particles**
        x[0] = 8.5; y[0] = 11.5; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;

        x[1] = 11.5; y[1] = 8.5; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; m[1] = m0;
    }
    else if (initial_condition == "ic-two-pass2") {
        // **Test Case 5: Two passing particles close**
        x[0] = 8.5; y[0] = 11.3; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;

        x[1] = 11.5; y[1] = 8.7; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; m[1] = m0;
    }
    else if (initial_condition == "ic-two-pass3") {
        // **Test Case 6: Two passing heavy particles**
        x[0] = 8.5; y[0] = 11.3; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 1; m[0] = m1;

        x[1] = 11.5; y[1] = 8.7; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 1; m[1] = m1;
    }
    else if (initial_condition == "ic-random") {
        // **Random initialization for N particles**
        int inited=0;
        for (int i=0; i<N; i++){
        

            if (i < N * percent_type1 / 100.0) {
                type[i] = 1;
                m[i] = m1;
            } else {
                type[i] = 0;
                m[i] = m0;
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

void compute_forces(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep) {

    // Reset forces
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);

    // Lennard-Jones parameters
    const double epsilon[2][2] = {{3.0, 15.0}, {15.0, 60.0}};
    const double sigma[2][2] = {{1.0, 2.0}, {2.0, 3.0}};

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double r = 0.0;
    double eps =0.0;
    double sig = 0.0;
    int t1 = 0;
    int t2 = 0;
    double f = 0.0;
    
    double sig_r6 = 0.0;
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
                r= sqrt(r2);

                if (r<min_sep){
                    min_sep=r; //For Unit Testing
                }

                t1 = type[i];
                t2 = type[j];
                eps = epsilon[t1][t2];
                sig = sigma[t1][t2];
                sig_r6 = pow(sig / r, 6);
                f = 24 * eps * (2 * sig_r6 * sig_r6 - sig_r6) / (r2);

                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;    
                fz[j] -= f * dz;
    }}}}

void update_positions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double dt) {
    for (int i = 0; i < N; i++) {
        x[i] += dt * vx[i];
        y[i] += dt * vy[i];
        z[i] += dt * vz[i];
    }
}
    

void update_velocities(double* vx, double* vy, double* vz, double* fx, double* fy, double* fz, double* m, int N, double dt) {
    for (int i = 0; i < N; i++) {
        vx[i] += dt * fx[i] / m[i];
        vy[i] += dt * fy[i] / m[i];
        vz[i] += dt * fz[i] / m[i];
    }
}


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

double compute_temperature(double* vx, double* vy, double* vz, double* m, int N) {
    double kinetic_energy = 0.0;
    double boltz=0.8314459920816467;
    for (int i = 0; i < N; i++) {
        kinetic_energy += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    return (2.0 / (3.0 *boltz* N)) * kinetic_energy;  
}
double compute_KE(double* vx, double* vy, double* vz, double* m, int N) {
    double kinetic_energy = 0.0;
    
    for (int i = 0; i < N; i++) {
        kinetic_energy += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    return kinetic_energy;  
}
void scale_velocities(double* vx, double* vy, double* vz, double* m, int N, double T_target) {
    double T_current = compute_temperature(vx, vy, vz, m, N);

    if (T_current == 0.0) return;  // Avoid division by zero

        double scale_factor = sqrt(T_target / T_current);

    for (int i = 0; i < N; i++) {
        vx[i] *= scale_factor;
        vy[i] *= scale_factor;
        vz[i] *= scale_factor;
    }
}

int main(int argc, char** argv) {
    // **Command-line Configuration**
    double Lx = 20.0, Ly = 20.0, Lz = 20.0;
    double dt = 0.001, T_tot = 50.0;
    int N = 8;
    double percent_type1 = 10.0;
    string initial_condition;
    double temperature = 0.0;
    double min_sep = Lx;
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "Print available options")
        ("Lx", po::value<double>(&Lx)->default_value(20.0), "x length (Angstroms)")
        ("Ly", po::value<double>(&Ly)->default_value(20.0), "y length (Angstroms)")
        ("Lz", po::value<double>(&Lz)->default_value(20.0), "z length (Angstroms)")
        ("dt", po::value<double>(&dt)->default_value(0.001), "Time step")
        ("T", po::value<double>(&T_tot)->default_value(10.0), "Final time")
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

    // **Determine Initial Condition**
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


    

    
    init_particle(x, y, z, vx, vy, vz, type, m, N, Lx, Ly, Lz, m0, m1,percent_type1, initial_condition);


    // Temperature Change - only if temp is set
    if (temperature > 0.0) {
        scale_velocities(vx, vy, vz, m, N, temperature);
    }

    // Create text files in overwrite mode
    std::ofstream particle_file("particles.txt", std::ofstream::trunc);
    std::ofstream kinetic_file("kinetic_energy.txt", std::ofstream::trunc);


    /////////////////////// Numerical Loop /////////////////////////
    int steps = T_tot / dt;
    double time=0.0;
    int writestep=static_cast<int>(0.1 / dt); // Cast double division to int
    for (int t = 0; t < steps; t++) {
        
        if (t%writestep==0) {  // Write data every 0.1 time units
            double K = compute_KE(vx, vy, vz, m, N);
            kinetic_file << time << " " << K << "\n";

            for (int i = 0; i < N; i++) {
                particle_file << time << " " << i << " " << type[i] << " "
                              << x[i] << " " << y[i] << " " << z[i] << " "
                              << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
            }
        }
        compute_forces(x, y, z, fx, fy, fz, type, N, min_sep);
        update_velocities(vx, vy, vz, fx, fy, fz, m, N, dt);
        update_positions(x, y, z, vx, vy, vz, N, dt);
        apply_boundary_conditions(x, y, z, vx, vy, vz, N, Lx, Ly, Lz);

        time += dt;
    }
    cout<<min_sep<<endl;


    /////// Cleanup Section ///////////
    delete[] x; delete[] y; delete[] z;
    delete[] vx; delete[] vy; delete[] vz;
    delete[] fx; delete[] fy; delete[] fz;
    delete[] m; delete[] type;
    particle_file.close();
    kinetic_file.close();
    ///////////////////////////////////
    
}