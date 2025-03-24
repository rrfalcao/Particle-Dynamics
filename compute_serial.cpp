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
using namespace std;
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
void update_positions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double dt) {
    
        for(int i=0; i<N; i++){
            x[i] += dt * vx[i];
            y[i] += dt * vy[i];
            z[i] += dt * vz[i];
    }}
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
    
        return (2.0 / (3.0 * boltz * (N))) * kinetic_energy;  
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