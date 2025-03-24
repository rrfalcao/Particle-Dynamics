#include <chrono>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include <limits>
#include <fstream> 
#include <random>
using namespace std;
void compute_forces_test(double* x, double* y, double* z,double* fx, double* fy, double* fz,int* type, int N, double& min_sep, char test) {

    // Reset global forces
    #pragma omp parallel for simd
    for (int i = 0; i < N; i++) {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
    }

    // Lennard-Jones parameters
    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    int total_pairs = N * (N - 1) / 2;
    int num_threads = omp_get_max_threads();

    // Allocate flat buffers: one block of N for each thread
    double* fx_priv = new double[num_threads * N]();
    double* fy_priv = new double[num_threads * N]();
    double* fz_priv = new double[num_threads * N]();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        #pragma omp for schedule(static)
        for (int k = 0; k < total_pairs; ++k) {
            int i = static_cast<int>(floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2));
            int j = k - (i * (2 * N - i - 1) / 2) + i + 1;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0.0) {
                if (test == 'y') {
                    double r = sqrt(r2);
                    #pragma omp critical
                    {
                        if (r < min_sep) min_sep = r;
                    }
                }

                int t1 = type[i];
                int t2 = type[j];
                double sig6 = sigma6[t1][t2];
                double eps24 = epsilon24[t1][t2];
                double r6 = r2 * r2 * r2;
                double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

                fx_priv[tid * N + i] += f * dx;
                fy_priv[tid * N + i] += f * dy;
                fz_priv[tid * N + i] += f * dz;

                fx_priv[tid * N + j] -= f * dx;
                fy_priv[tid * N + j] -= f * dy;
                fz_priv[tid * N + j] -= f * dz;
            }
        }
    }

    // Final reduction
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int t = 0; t < num_threads; t++) {
            fx[i] += fx_priv[t * N + i];
            fy[i] += fy_priv[t * N + i];
            fz[i] += fz_priv[t * N + i];
        }
    }

    // Clean up
    delete[] fx_priv;
    delete[] fy_priv;
    delete[] fz_priv;
}
void compute_forces(double* x, double* y, double* z,double* fx, double* fy, double* fz,int* type, int N) {

    // Reset global forces
    #pragma omp parallel for simd
    for (int i = 0; i < N; i++) {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
    }

    // Lennard-Jones parameters
    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    int total_pairs = N * (N - 1) / 2;
    int num_threads = omp_get_max_threads();

    // Force buffers, each row is a thread
    double* fx_priv = new double[num_threads * N]();
    double* fy_priv = new double[num_threads * N]();
    double* fz_priv = new double[num_threads * N]();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        #pragma omp for schedule(static)
        for (int k = 0; k < total_pairs; ++k) {
            int i = (floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2));
            int j = k - (i * (2 * N - i - 1) / 2) + i + 1;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0.0) {
        
                int t1 = type[i];
                int t2 = type[j];
                double sig6 = sigma6[t1][t2];
                double eps24 = epsilon24[t1][t2];
                double r6 = r2 * r2 * r2;
                double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

                fx_priv[tid * N + i] += f * dx;
                fy_priv[tid * N + i] += f * dy;
                fz_priv[tid * N + i] += f * dz;

                fx_priv[tid * N + j] -= f * dx;
                fy_priv[tid * N + j] -= f * dy;
                fz_priv[tid * N + j] -= f * dz;
            }
        }
    }

    // Final reduction
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int t = 0; t < num_threads; t++) {
            fx[i] += fx_priv[t * N + i];
            fy[i] += fy_priv[t * N + i];
            fz[i] += fz_priv[t * N + i];
        }
    }

    // Clean up
    delete[] fx_priv;
    delete[] fy_priv;
    delete[] fz_priv;
}
void update_positions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double dt) {
    
    #pragma omp parallel for simd
    for(int i=0; i<N; i++){
        
        x[i] += dt * vx[i];
        y[i] += dt * vy[i];
        z[i] += dt * vz[i]; 
        
}}
void update_velocities(double* vx, double* vy, double* vz, double* fx, double* fy, double* fz, double* m, int N, double dt) {
    #pragma omp parallel for simd
    for (int i = 0; i < N; i++) {
        vx[i] += dt * fx[i] / m[i];
        vy[i] += dt * fy[i] / m[i];
        vz[i] += dt * fz[i] / m[i];
    }
}
void apply_boundary_conditions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double Lx, double Ly, double Lz) {
    #pragma omp parallel for
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
    #pragma omp parallel for simd reduction(+:kinetic_energy)
    for (int i = 0; i < N; i++) {
        kinetic_energy += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    return (2.0 / (3.0 *boltz* N)) * kinetic_energy;  
}
double compute_KE(double* vx, double* vy, double* vz, double* m, int N) {
    double kinetic_energy = 0.0;
    #pragma omp parallel for simd reduction(+:kinetic_energy)
    for (int i = 0; i < N; i++) {
        kinetic_energy += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    return kinetic_energy;  
}
void scale_velocities(double* vx, double* vy, double* vz, double* m, int N, double T_target) {
    double T_current = compute_temperature(vx, vy, vz, m, N);

    if (T_current == 0.0) return;  // Avoid division by zero

        double scale_factor = sqrt(T_target / T_current);
    #pragma omp parallel for simd
    for (int i = 0; i < N; i++) {
        vx[i] *= scale_factor;
        vy[i] *= scale_factor;
        vz[i] *= scale_factor;
    }
}