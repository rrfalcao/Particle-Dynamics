#include <chrono>  // For high-resolution timing
#include <iostream> // For printing
#include <cmath>
#include <omp.h>
#include <vector>
using namespace std;


    void compute_forces_opt(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep,char test) {

        // Reset forces
        // #pragma omp parallel for simd
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
        }
    
        // Lennard-Jones parameters
        const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
        const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};
    
       
        
        
        // Loop over all pairs of particles
    
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; i++) {
            #pragma omp simd
            for (int j = i + 1; j < N; j++) {
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];
    
                double r2 = dx * dx + dy * dy + dz * dz;
                // if (test=='y'){
                //   double r= sqrt(r2);
                
                //     if (r<min_sep){
                //         min_sep=r; //For Unit Testing
                //     }
                // }
                double r6 = r2 * r2 * r2;
                double f = epsilon24[type[i]][type[j]] * sigma6[type[i]][type[j]] * (2 * sigma6[type[i]][type[j]] - r6) / (r6 * r6 * r2);
    
                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;
                fz[j] -= f * dz;
            }
        }
    
    }
        

    void compute_forces_sing(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
        // Reset forces efficiently using SIMD
        #pragma omp parallel for simd
        for (int i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
        }
    
        // Lennard-Jones parameters
        const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
        const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};
    
        int total_pairs = N * (N - 1) / 2;  // Total interactions
    
        // **Use Reduction to Eliminate the Critical Section**
        #pragma omp parallel for schedule(dynamic) reduction(+:fx[:N], fy[:N], fz[:N])
        for (int k = 0; k < total_pairs; k++) {
            int i = floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2);
            int j = k - (i * (2 * N - i - 1) / 2) + i + 1;
    
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
    
            double r2 = dx * dx + dy * dy + dz * dz;
            double r6 = r2 * r2 * r2;
    
            // Avoid recalculating values by storing them
            double epsilon = epsilon24[type[i]][type[j]];
            double sigma = sigma6[type[i]][type[j]];
            double sigma_r6 = sigma * sigma - r6;
            double f = epsilon * sigma_r6 / (r6 * r6 * r2);
    
            fx[i] += f * dx;
            fy[i] += f * dy;
            fz[i] += f * dz;
            fx[j] -= f * dx;
            fy[j] -= f * dy;
            fz[j] -= f * dz;
        }
    }
    void compute_forces_new(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
        // Reset forces in parallel
        #pragma omp parallel for simd
        for (int i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
        }
    
        // Lennard-Jones parameters
        const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
        const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};
    
        // Parallel reduction across all threads
        #pragma omp parallel
        {
            // Allocate thread-local force accumulators (avoid memory allocation overhead)
            std::vector<double> fx_private(N, 0.0);
            std::vector<double> fy_private(N, 0.0);
            std::vector<double> fz_private(N, 0.0);
    
            // Distribute iterations among threads
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < N; i++) {
                #pragma omp simd
                for (int j = i + 1; j < N; j++) {
                    double dx = x[i] - x[j];
                    double dy = y[i] - y[j];
                    double dz = z[i] - z[j];
                    double r2 = dx * dx + dy * dy + dz * dz;
                    double r6 = r2 * r2 * r2;
                    double f = epsilon24[type[i]][type[j]] * sigma6[type[i]][type[j]] * (2 * sigma6[type[i]][type[j]] - r6) / (r6 * r6 * r2);
    
                    // Accumulate forces in private buffers
                    fx_private[i] += f * dx;
                    fy_private[i] += f * dy;
                    fz_private[i] += f * dz;
                    fx_private[j] -= f * dx;
                    fy_private[j] -= f * dy;
                    fz_private[j] -= f * dz;
                }
            }
    
            // Reduce thread-private arrays into global force arrays
            #pragma omp for simd reduction(+:fx[:N], fy[:N], fz[:N])
            for (int i = 0; i < N; i++) {
                fx[i] += fx_private[i];
                fy[i] += fy_private[i];
                fz[i] += fz_private[i];
            }
        }
    }
    void compute_forces_opt_sing(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep, char test) {
        // Reset forces efficiently using SIMD
        #pragma omp parallel for simd
        for (int i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
        }
    
        // Lennard-Jones parameters
        const double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
        const double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};
    
        int num_threads = omp_get_max_threads();  // Get available threads
        int total_pairs = N * (N - 1) / 2;  // Total interactions
        int chunk_size = total_pairs / num_threads; // Even work division
    
        // ** Use Thread-Private Buffers for Forces **
        #pragma omp parallel
        {
            
            double* fx_private = new double[N]();  // Local thread buffers
            double* fy_private = new double[N]();
            double* fz_private = new double[N]();
    
            // ** Compute Forces with Static Scheduling **
            #pragma omp for schedule(static, chunk_size)
            for (int k = 0; k < total_pairs; k++) {
                int i = floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2);
                int j = k - (i * (2 * N - i - 1) / 2) + i + 1;
    
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];
    
                double r2 = dx * dx + dy * dy + dz * dz;
                double r6 = r2 * r2 * r2;
                double f = epsilon24[type[i]][type[j]] * sigma6[type[i]][type[j]] * (2 * sigma6[type[i]][type[j]] - r6) / (r6 * r6 * r2);
    
                // ** Use SIMD to accelerate force updates **
    
                    fx_private[i] += f * dx;    
                    fy_private[i] += f * dy;
                    fz_private[i] += f * dz;
                    fx_private[j] -= f * dx;
                    fy_private[j] -= f * dy;
                    fz_private[j] -= f * dz;
                
            }
    
            // ** Reduction: Sum Thread-Private Buffers Into Global Force Arrays **
            #pragma omp critical
            {
                for (int i = 0; i < N; i++) {
                    fx[i] += fx_private[i];
                    fy[i] += fy_private[i];
                    fz[i] += fz_private[i];
                }
            }
    
            // Cleanup Thread-Local Memory
            delete[] fx_private;
            delete[] fy_private;
            delete[] fz_private;
        }
    }
void compute_forces2(double* x, double* y, double* z,
                    double* fx, double* fy, double* fz,
                    int* type, int N, double& min_sep, char test) {

    // Reset forces
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

    // Thread-private force arrays
    vector<vector<double>> fx_private(num_threads, vector<double>(N, 0.0));
    vector<vector<double>> fy_private(num_threads, vector<double>(N, 0.0));
    vector<vector<double>> fz_private(num_threads, vector<double>(N, 0.0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        

        #pragma omp for schedule(static)
        for (int k = 0; k < total_pairs; ++k) {
            // Map flat index to pair (i, j)
            int i = (floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2));
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

                fx_private[tid][i] += f * dx;
                fy_private[tid][i] += f * dy;
                fz_private[tid][i] += f * dz;
                fx_private[tid][j] -= f * dx;
                fy_private[tid][j] -= f * dy;
                fz_private[tid][j] -= f * dz;
            }
        }
    }

    // Final reduction into global force arrays
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int t = 0; t < num_threads; t++) {
            fx[i] += fx_private[t][i];
            fy[i] += fy_private[t][i];
            fz[i] += fz_private[t][i];
        }
    }
}
void compute_forces3(double* x, double* y, double* z,
                    double* fx, double* fy, double* fz,
                    int* type, int N, double& min_sep, char test) {

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
    double* fx_private = new double[num_threads * N]();
    double* fy_private = new double[num_threads * N]();
    double* fz_private = new double[num_threads * N]();

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

                fx_private[tid * N + i] += f * dx;
                fy_private[tid * N + i] += f * dy;
                fz_private[tid * N + i] += f * dz;

                fx_private[tid * N + j] -= f * dx;
                fy_private[tid * N + j] -= f * dy;
                fz_private[tid * N + j] -= f * dz;
            }
        }
    }

    // Final reduction
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int t = 0; t < num_threads; t++) {
            fx[i] += fx_private[t * N + i];
            fy[i] += fy_private[t * N + i];
            fz[i] += fz_private[t * N + i];
        }
    }

    // Clean up
    delete[] fx_private;
    delete[] fy_private;
    delete[] fz_private;
}

int main() {
    // **Allocate Memory for Particles (example values)**
    int N = 10000; // Change this based on your simulation size
    double* x = new double[N]; 
    double* y = new double[N]; 
    double* z = new double[N];
    double* fx = new double[N]; 
    double* fy = new double[N]; 
    double* fz = new double[N];
    int* type = new int[N];
    double min_sep = 0.0;
    char test = 'n';

    // **Initialize particles (dummy values for testing)**
    for (int i = 0; i < N; i++) {
        srand(time(0));
        x[i] = static_cast<double>(rand()) / RAND_MAX * 50.0;
        y[i] = static_cast<double>(rand()) / RAND_MAX * 50.0;
        z[i] = static_cast<double>(rand()) / RAND_MAX * 50.0;
        type[i] = rand() % 2; // Randomly assign types 0 or 1
    }
    double start_time1, end_time1;
    double start_time2, end_time2;
    double start_time3, end_time3;
    
    
    // **Timing two loop compute_forces**
    start_time1 = omp_get_wtime();  // Start timer
    compute_forces3(x, y, z, fx, fy, fz, type, N, min_sep,test);
    end_time1 = omp_get_wtime();  // End timer

    double elapsed_time1 = end_time1 - start_time1;
    
    cout << "Execution time one loop: " << elapsed_time1 << " seconds\n";

    // **Timing one loop compute_forces **
    start_time2 = omp_get_wtime();  // Start timer
    compute_forces_opt_sing(x, y, z, fx, fy, fz, type, N, min_sep,test);
    end_time2 = omp_get_wtime();  // End timer

    double elapsed_time2 = end_time2 - start_time2;
    
    cout << "Execution time one loop improved: " << elapsed_time2 << " seconds\n";

    start_time3 = omp_get_wtime();  // Start timer
    compute_forces2(x, y, z, fx, fy, fz, type, N, min_sep,test);
    end_time3 = omp_get_wtime();  // End timer
    double elapsed_time3 = end_time3 - start_time3;
    cout << "Execution time new two loop improved: " << elapsed_time3 << " seconds\n";


    // **Cleanup Memory**
    delete[] x; delete[] y; delete[] z;
    delete[] fx; delete[] fy; delete[] fz;
    delete[] type;

    return 0;
}
