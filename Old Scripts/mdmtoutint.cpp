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

void compute_forces(double* x, double* y, double* z, 
    double* fx, double* fy, double* fz, 
    int N0, int N1, double& min_sep, char test) {

    // Reset forces
    fill(fx, fx + (N0 + N1), 0.0);
    fill(fy, fy + (N0 + N1), 0.0);
    fill(fz, fz + (N0 + N1), 0.0);

    // Lennard-Jones parameters
    

    double dx, dy, dz, r2, r6, f;

    // Pointer offsets:
    double* x0 = x + N1;  // Type 0 particles start at index N1
    double* y0 = y + N1;
    double* z0 = z + N1;
    double* fx0 = fx + N1;
    double* fy0 = fy + N1;
    double* fz0 = fz + N1;
    double r=0.0;
    // **Step 1: Type 1 - Type 1 interactions**
    for (int i = 0; i < N1; i++) {
        for (int j = i + 1; j < N1; j++) {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        dz = z[i] - z[j];

        r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > 0) {
            
            if (test=='y'){
                r= sqrt(r2);
            
                if (r<min_sep){
                    min_sep=r; //For Unit Testing
                }
            }
            
            
            r6 = r2 * r2 * r2;
            f = 1049760 * (1458 - r6) / (r6 * r6 * r2);

            fx[i] += f * dx;
            fy[i] += f * dy;
            fz[i] += f * dz;
            fx[j] -= f * dx;
            fy[j] -= f * dy;
            fz[j] -= f * dz;
        }
        }
    }

// **Step 2: Type 1 - Type 0 interactions**
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N0; j++) {
            dx = x[i] - x0[j];
            dy = y[i] - y0[j];
            dz = z[i] - z0[j];

            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0) {
                if (test=='y'){
                    r= sqrt(r2);
                
                    if (r<min_sep){
                        min_sep=r; //For Unit Testing
                    }
                }
                
                
                r6 = r2 * r2 * r2;
                f = 23040 * (128 - r6) / (r6 * r6 * r2);

                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx0[j] -= f * dx;
                fy0[j] -= f * dy;
                fz0[j] -= f * dz;
            }
        }
    }

    // **Step 3: Type 0 - Type 0 interactions**
    for (int i = 0; i < N0; i++) {
        for (int j = i + 1; j < N0; j++) {
            dx = x0[i] - x0[j];
            dy = y0[i] - y0[j];
            dz = z0[i] - z0[j];

            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0) {
                if (test=='y'){
                    r= sqrt(r2);
                
                    if (r<min_sep){
                        min_sep=r; //For Unit Testing
                    }
                }
                

                r6 = r2 * r2 * r2;
                
                f = 72 * (2 - r6) / (r6 * r6 * r2);
                fx0[i] += f * dx;
                fy0[i] += f * dy;
                fz0[i] += f * dz;
                fx0[j] -= f * dx;
                fy0[j] -= f * dy;
                fz0[j] -= f * dz;
            }
        }
    }
}


void update_positions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double dt) {
    cblas_daxpy(N, dt, vx, 1, x, 1);  // x = dt * vx + x
    cblas_daxpy(N, dt, vy, 1, y, 1);  // y = dt * vy + y
    cblas_daxpy(N, dt, vz, 1, z, 1);  // z = dt * vz + z
}
    

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




int main(int argc, char** argv) {
    // **Command-line Configuration**
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
    
    // Temperature Change - only if temp is set
    if (temperature > 0.0) {
        scale_velocities(vx, vy, vz, m0, m1, N0, N1, temperature);
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
            double K = compute_KE(vx, vy, vz, m0, m1, N0, N1);
            kinetic_file << time << " " << K << "\n";

            for (int i = 0; i < N; i++) {
                particle_file << time << " " << i << " " << type[i] << " "
                              << x[i] << " " << y[i] << " " << z[i] << " "
                              << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
            }
        }
        compute_forces(x, y, z, fx, fy, fz, N0,N1, min_sep,test);
        update_velocities(vx, vy, vz, fx, fy, fz, m0, m1, N0, N1, dt);

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
    particle_file.close();
    kinetic_file.close();
    ///////////////////////////////////
    
}

// void compute_forcesopt(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
//     fill(fx, fx + N, 0.0);
//     fill(fy, fy + N, 0.0);
//     fill(fz, fz + N, 0.0);

//     const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
//     const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

//     for (int i = 0; i < N - 1; i++) {
//         int len = N - i - 1;  // Remaining pairs
//         double dx[len], dy[len], dz[len];

        
//         cblas_dcopy(len, x + i + 1, 1, dx, 1);
//         cblas_daxpy(len, -1.0, x + i, 0, dx, 1);

//         cblas_dcopy(len, y + i + 1, 1, dy, 1);
//         cblas_daxpy(len, -1.0, y + i, 0, dy, 1);

//         cblas_dcopy(len, z + i + 1, 1, dz, 1);
//         cblas_daxpy(len, -1.0, z + i, 0, dz, 1);

        
//         double r2[len];
//         for (int j = 0; j < len; j++) {
//             r2[j] = dx[j] * dx[j] + dy[j] * dy[j] + dz[j] * dz[j];
            
//             if (r2[j] > 0) {
//                 if (test=='y'){
//                             double r= sqrt(r2[j]);       
//                             if (r<min_sep){
//                             min_sep=r; 
//                             }
//                 }

//                 int t1 = type[i], t2 = type[i + j + 1];
//                 double eps24 = epsilon24[t1][t2];
//                 double sig6 = sigma6[t1][t2];
                
                
//                 double r6 = r2[j] * r2[j] * r2[j];
//                 double f = -eps24 * sig6 *((2 * sig6 - r6) / (r6*r6*r2[j]));
                
//                 fx[i] += f * dx[j];
//                 fy[i] += f * dy[j];
//                 fz[i] += f * dz[j];

//                 fx[i + j + 1] -= f * dx[j];
//                 fy[i + j + 1] -= f * dy[j];
//                 fz[i + j + 1] -= f * dz[j];
//             }
//         }
//     }
// }