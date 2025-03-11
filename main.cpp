#include <iostream>
// #include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include <fstream> 
using namespace std;
#define F77NAME(x) x##_


std::ofstream particle_file("particles.txt");
std::ofstream kinetic_file("kinetic_energy.txt");

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

void init_particle(double* x, double* y, double* z, double * vx, double * vy, double* vz, int N, double Lx, double Ly, double Lz){
    srand(time(0));
    int inited=0;
    for (int i=0; i<N; i++){
    
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
        
    
    
}}

void compute_forces(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N) {

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
int main(){


    int N=2;
    const double dt = 0.001;  // Time step
    const double Lx = 20.0, Ly = 20.0, Lz = 20.0;  // Box dimensions
    double T_tot = 10.0;  // Total simulation time
    int steps = T_tot / dt;  // Number of time steps
    double p1=0.1;
    double* x=new double[N];  
    double* y=new double[N];
    double* z=new double[N];
    double* vx=new double[N];
    double* vy=new double[N];
    double* vz=new double[N];
    double* fx=new double[N];
    double* fy=new double[N];
    double* fz=new double[N];
    double* m=new double[N];
    double* r=new double[N];
    int m0 = 1;
    int m1 = 10;

    int* type=new int[N];

    //Assign types and masses
    for (int i = 0; i < N; i++) {
        if (i < N * p1) {
            type[i] = 0;
            m[i] = m0;
        } else {
            type[i] = 1;
            m[i] = m1;
        }
    }

    init_particle(x, y, z, vx, vy, vz, N, Lx, Ly, Lz);
    
    cout << "Initial temperature: " << compute_temperature(vx, vy, vz, m, N) << endl;
    for (int t = 0; t < steps; t++) {
        compute_forces(x, y, z, fx, fy, fz, type, N);

        double time = t * dt;  // Current time
    if (fmod(time, 0.1) < dt) {  // Every 0.1 time units, write data
        double K = compute_KE(vx, vy, vz, m, N);
        
        // Write kinetic energy data
        kinetic_file << time << " " << K << "\n";

        // Write particle data
        for (int i = 0; i < N; i++) {
            particle_file << time << " " << i << " "
                          << x[i] << " " << y[i] << " " << z[i] << " "
                          << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
        }
    }
    }

}