#include <chrono>  // For high-resolution timing
#include <iostream> // For printing
#include <cmath>
using namespace std;
void compute_forces(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep) {

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
                
                
                r= sqrt(r2);
                
                if (r<min_sep){
                   min_sep=r; //For Unit Testing    
                }

                t1 = type[i];
                t2 = type[j];
                eps24 = epsilon24[t1][t2];
                sig6 = sigma6[t1][t2];
                r6=r2*r2*r2;
                f = eps24 * sig6*(2*sig6 -r6) / (r6*r6*r2);

                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;    
                fz[j] -= f * dz;
    }}}}

    void compute_forces_opt(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep,char test) {

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
        

void compute_forces_optsing(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep,char test) {

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
    
        for (int k = 0; k < N * (N - 1) / 2; k++) {
            int i = floor((2 * N - 1 - sqrt((2 * N - 1) * (2 * N - 1) - 8 * k)) / 2);
            int j = k - (i * (2 * N - i - 1) / 2) + i + 1;
    
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
    
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0) {
                if (test=='y'){
                    r= sqrt(r2);
                
                    if (r<min_sep){
                        min_sep=r; //For Unit Testing
                    }
                }
                double r6 = r2 * r2 * r2;
                double f = epsilon24[type[i]][type[j]] * sigma6[type[i]][type[j]] * (2 * sigma6[type[i]][type[j]] - r6) / (r6 * r6 * r2);
    
                fx[i] += f * dx;
                fy[i] += f * dy;
                fz[i] += f * dz;
                fx[j] -= f * dx;
                fy[j] -= f * dy;
                fz[j] -= f * dz;
            }}}
       
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
        x[i] = static_cast<double>(rand()) / RAND_MAX * 10.0;
        y[i] = static_cast<double>(rand()) / RAND_MAX * 10.0;
        z[i] = static_cast<double>(rand()) / RAND_MAX * 10.0;
        type[i] = rand() % 2; // Randomly assign types 0 or 1
    }

    // **Timing OLD compute_forces**
    auto start_old = std::chrono::high_resolution_clock::now();
    compute_forces(x, y, z, fx, fy, fz, type, N, min_sep);
    auto end_old = std::chrono::high_resolution_clock::now();
    double time_old = std::chrono::duration<double, std::milli>(end_old - start_old).count();
    
    std::cout << "Old compute_forces execution time: " << time_old << " ms\n";

    // **Timing NEW compute_forces (optimized)**
    auto start_new = std::chrono::high_resolution_clock::now();
    compute_forces_opt(x, y, z, fx, fy, fz, type, N, min_sep, test);
    auto end_new = std::chrono::high_resolution_clock::now();
    double time_new = std::chrono::duration<double, std::milli>(end_new - start_new).count();

    std::cout << "New compute_forces execution time: " << time_new << " ms\n";

    // **Timing NEW compute_forces (optimized)**
    auto start_newsing = std::chrono::high_resolution_clock::now();
    compute_forces_optsing(x, y, z, fx, fy, fz, type, N, min_sep, test);
    auto end_newsing = std::chrono::high_resolution_clock::now();
    double time_newsing = std::chrono::duration<double, std::milli>(end_newsing - start_newsing).count();
    std::cout << "New single loop compute_forces execution time: " << time_new << " ms\n";
    // **Compare speedup**
    if (time_new < time_old) {
        std::cout << "Optimized compute_forces is " << (time_old / time_new) << "x faster!\n";
    } else {
        std::cout << "Warning: Optimized compute_forces is slower than the original.\n";
    }

    // **Cleanup Memory**
    delete[] x; delete[] y; delete[] z;
    delete[] fx; delete[] fy; delete[] fz;
    delete[] type;

    return 0;
}
