#include <chrono>
#include <random>
#include <iostream>
#include <cmath>
using namespace std;


void compute_forces_v0(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);

    const double epsilon[2][2] = {{3.0, 15.0}, {15.0, 60.0}};
    const double sigma[2][2] = {{1.0, 2.0}, {2.0, 3.0}};

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > 0.0) {
                double r = sqrt(r2);
                if (test == 'y' && r < min_sep) min_sep = r;

                int t1 = type[i], t2 = type[j];
                double eps = epsilon[t1][t2];
                double sig = sigma[t1][t2];
                double sig_r6 = pow(sig / r, 6);
                double f = 24.0 * eps * (2 * sig_r6 * sig_r6 - sig_r6) / r2;

                fx[i] += f * dx; fy[i] += f * dy; fz[i] += f * dz;
                fx[j] -= f * dx; fy[j] -= f * dy; fz[j] -= f * dz;
            }
        }
    }
}
void compute_forces_v1(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);

    const double epsilon[2][2] = {{3.0, 15.0}, {15.0, 60.0}};
    const double sigma[2][2] = {{1.0, 2.0}, {2.0, 3.0}};

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > 0.0) {
                double r = sqrt(r2);
                if (test == 'y' && r < min_sep) min_sep = r;

                int t1 = type[i], t2 = type[j];
                double eps = epsilon[t1][t2];
                double sig = sigma[t1][t2];
                double inv_r = 1.0 / r;
                double sig_r = sig * inv_r;
                double sig_r2 = sig_r * sig_r;
                double sig_r6 = sig_r2 * sig_r2 * sig_r2;
                double f = 24.0 * eps * (2 * sig_r6 * sig_r6 - sig_r6) / r2;

                fx[i] += f * dx; fy[i] += f * dy; fz[i] += f * dz;
                fx[j] -= f * dx; fy[j] -= f * dy; fz[j] -= f * dz;
            }
        }
    }
}
void compute_forces_v2(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);

    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > 0.0) {
                double r = sqrt(r2);
                if (test == 'y' && r < min_sep) min_sep = r;

                int t1 = type[i], t2 = type[j];
                double eps24 = epsilon24[t1][t2];
                double sig6 = sigma6[t1][t2];
                double r6 = r2 * r2 * r2;
                double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

                fx[i] += f * dx; fy[i] += f * dy; fz[i] += f * dz;
                fx[j] -= f * dx; fy[j] -= f * dy; fz[j] -= f * dz;
            }
        }
    }
}
void compute_forces_v3(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);
//IGNORE THIS FUNCTION
    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > 0.0) {
                if (test == 'y') {
                    double r = sqrt(r2);
                    if (r < min_sep) min_sep = r;
                }

                int t1 = type[i], t2 = type[j];
                double eps24 = epsilon24[t1][t2];
                double sig6 = sigma6[t1][t2];
                double r6 = r2 * r2 * r2;
                double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

                fx[i] += f * dx; fy[i] += f * dy; fz[i] += f * dz;
                fx[j] -= f * dx; fy[j] -= f * dy; fz[j] -= f * dz;
            }
        }
    }
}
void compute_forces_v4(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
    fill(fx, fx + N, 0.0);
    fill(fy, fy + N, 0.0);
    fill(fz, fz + N, 0.0);

    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 > 0.0) {
                if (test == 'y') {
                    double r = sqrt(r2);
                    if (r < min_sep) min_sep = r;
                }

                int t1 = type[i];
                int t2 = type[j];
                double sig6 = sigma6[t1][t2];
                double eps24 = epsilon24[t1][t2];
                double r6 = r2 * r2 * r2;
                double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;
                fz[j] -= f * dz;
            }
        }
    }
}
void compute_forces_v5(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
    

    const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
    const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};

    for (int i = 0; i < N; ++i) {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 > 0.0) {
                if (test == 'y') {
                    double r = sqrt(r2);
                    if (r < min_sep) min_sep = r;
                }

                int t1 = type[i];
                int t2 = type[j];
                double sig6 = sigma6[t1][t2];
                double eps24 = epsilon24[t1][t2];
                double r6 = r2 * r2 * r2;
                double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;
                fz[j] -= f * dz;
            }
        }
    }
}
void compute_forces_v6(double* x, double* y, double* z,
                       double* fx, double* fy, double* fz,
                       int* type, int N, double& min_sep, char test) {
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

    for (int k = 0; k < N * (N - 1) / 2; k++) {
        int i = floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2);
        int j = k - (i * (2 * N - i - 1) / 2) + i + 1;
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
    }}
}


void benchmark_compute_forces(int N, int version) {
    // Allocate memory
    double* x = new double[N];
    double* y = new double[N];
    double* z = new double[N];
    double* fx = new double[N];
    double* fy = new double[N];
    double* fz = new double[N];
    int* type = new int[N];
    double min_sep = 1e9;
    char test = 'n';

    // Random init
    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<double> dist(0.0, 20.0);
    for (int i = 0; i < N; ++i) {
        x[i] = dist(gen);
        y[i] = dist(gen);
        z[i] = dist(gen);
        type[i] = i % 2;
    }

    const int num_runs = 100;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_runs; ++i) {
        switch(version) {
            case 0: compute_forces_v0(x, y, z, fx, fy, fz, type, N, min_sep, test); break;
            case 1: compute_forces_v1(x, y, z, fx, fy, fz, type, N, min_sep, test); break;
            case 2: compute_forces_v2(x, y, z, fx, fy, fz, type, N, min_sep, test); break;
            case 3: compute_forces_v3(x, y, z, fx, fy, fz, type, N, min_sep, test); break;
            case 4: compute_forces_v4(x, y, z, fx, fy, fz, type, N, min_sep, test); break;
            case 5: compute_forces_v5(x, y, z, fx, fy, fz, type, N, min_sep, test); break;
            case 6: compute_forces_v6(x, y, z, fx, fy, fz, type, N, min_sep, test); break;}
        }
    auto end = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "Version " << version << ": "
              << (duration_ms / num_runs) << " ms per call (avg over " << num_runs << " runs)\n";

    delete[] x; delete[] y; delete[] z;
    delete[] fx; delete[] fy; delete[] fz;
    delete[] type;
}
int main() {
    int N = 1000;
    for (int v = 0; v <= 6; ++v) {
        benchmark_compute_forces(N, v);
    }
    return 0;
}
