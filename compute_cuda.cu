#include "compute_cuda.h"

__constant__ double epsilon24[2][2] = {{72.0, 360.0}, {360.0, 1440.0}};
__constant__ double sigma6[2][2] = {{1.0, 64.0}, {64.0, 729.0}};
__global__ void compute_forces_gpu(double* x, double* y, double* z,double* fx, double* fy, double* fz,int* type, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;

    for (int j = i + 1; j < N; j++) {  // Avoid redundant calculations
        
        double dx = x[i] - x[j];
        double dy = y[i] - y[j];
        double dz = z[i] - z[j];

        double r2 = dx * dx + dy * dy + dz * dz;
        double r6 = r2 * r2 * r2;
        double sig6 = sigma6[type[i]][type[j]];
        double eps24 = epsilon24[type[i]][type[j]];

        double f = eps24 * sig6 * (2 * sig6 - r6) / (r6 * r6 * r2);

        fx_i += f * dx;
        fy_i += f * dy;
        fz_i += f * dz;
        
        fx[j] -= f * dx;
        fy[j] -= f * dy;
        fz[j] -= f * dz;
        
    }

    // Store computed forces back in global memory
    fx[i] += fx_i;
    fy[i] += fy_i;
    fz[i] += fz_i;
    
}
void compute_forces(double* x, double* y, double* z,double* fx, double* fy, double* fz,int* type, int N) {
    int num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // Reset force arrays
    cudaMemset(fx, 0, N * sizeof(double));
    cudaMemset(fy, 0, N * sizeof(double));
    cudaMemset(fz, 0, N * sizeof(double));

    compute_forces_gpu<<<num_blocks, BLOCK_SIZE>>>(x, y, z, fx, fy, fz, type, N);
    
}
__global__ void update_positions_gpu(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    x[i] += dt * vx[i];
    y[i] += dt * vy[i];
    z[i] += dt * vz[i];
}
void update_positions(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double dt) {
    int num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    update_positions_gpu<<<num_blocks, BLOCK_SIZE>>>(x, y, z, vx, vy, vz, N, dt);
}
__global__ void update_velocities_gpu(double* vx, double* vy, double* vz,double* fx, double* fy, double* fz,double* m, int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    vx[i] += dt * fx[i] / m[i];
    vy[i] += dt * fy[i] / m[i];
    vz[i] += dt * fz[i] / m[i];
}
void update_velocities(double* vx, double* vy, double* vz,double* fx, double* fy, double* fz,double* m, int N, double dt) {
    int num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    update_velocities_gpu<<<num_blocks, BLOCK_SIZE>>>(vx, vy, vz, fx, fy, fz, m, N, dt);
    
}
__global__ void apply_boundary_conditions_gpu(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double Lx, double Ly, double Lz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    if (x[i] < 0) { x[i] = -x[i]; vx[i] = -vx[i]; }
    if (x[i] > Lx) { x[i] = 2 * Lx - x[i]; vx[i] = -vx[i]; }

    if (y[i] < 0) { y[i] = -y[i]; vy[i] = -vy[i]; }
    if (y[i] > Ly) { y[i] = 2 * Ly - y[i]; vy[i] = -vy[i]; }

    if (z[i] < 0) { z[i] = -z[i]; vz[i] = -vz[i]; }
    if (z[i] > Lz) { z[i] = 2 * Lz - z[i]; vz[i] = -vz[i]; }
}
void apply_boundary_conditions(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double Lx, double Ly, double Lz) {
    int num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    apply_boundary_conditions_gpu<<<num_blocks, BLOCK_SIZE>>>(x, y, z, vx, vy, vz, N, Lx, Ly, Lz);
    
}
__global__ void compute_temperature_gpu(const double* vx, const double* vy, const double* vz,const double* m, double* partial_sum, int N) {
    __shared__ double temp_sum[BLOCK_SIZE];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double ke = 0.0;
    if (i < N) {
        double v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        ke = 0.5 * m[i] * v2;
    }
    temp_sum[tid] = ke;

    __syncthreads();

    // Reduce within the block
    for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
        if (tid < s)
        temp_sum[tid] += temp_sum[tid + s];
        __syncthreads();
}

// Write result from each block to global memory
    if (tid == 0)
        partial_sum[blockIdx.x] = temp_sum[0];
}
double compute_KE(double* vx, double* vy, double* vz, double* m, int N) {
    int num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    double* d_partial_sum;
    double* h_partial_sum = new double[num_blocks];
    cudaMallocManaged(&d_partial_sum, num_blocks * sizeof(double));

    compute_temperature_gpu<<<num_blocks, BLOCK_SIZE>>>(vx, vy, vz, m, d_partial_sum, N);
    cudaDeviceSynchronize();
    cudaMemcpy(h_partial_sum, d_partial_sum, num_blocks * sizeof(double), cudaMemcpyDeviceToHost);

    double KE_Tot = 0.0;
    for (int i = 0; i < num_blocks; ++i) {
        KE_Tot += h_partial_sum[i];
    }

    cudaFree(d_partial_sum);
    delete[] h_partial_sum;
    return KE_Tot;
}
__global__ void scale_velocities_gpu(double* vx, double* vy, double* vz, int N, double scale_factor) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        vx[i] *= scale_factor;
        vy[i] *= scale_factor;
        vz[i] *= scale_factor;
    }
}
void scale_velocities(double* vx, double* vy, double* vz, double* m, int N, double T_target) {
    const double boltz = 0.8314459920816467;
    int num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    double* d_partial_sum;
    double* h_partial_sum = new double[num_blocks];
    cudaMalloc(&d_partial_sum, num_blocks * sizeof(double));

    compute_temperature_gpu<<<num_blocks, BLOCK_SIZE>>>(vx, vy, vz, m, d_partial_sum, N);
    cudaDeviceSynchronize();
    cudaMemcpy(h_partial_sum, d_partial_sum, num_blocks * sizeof(double), cudaMemcpyDeviceToHost);

    double KE_Tot = 0.0;
    for (int i = 0; i < num_blocks; ++i) {
        KE_Tot += h_partial_sum[i];
    }

    cudaFree(d_partial_sum);
    delete[] h_partial_sum;

    double T_current = (2.0 / (3.0 * boltz * N)) * KE_Tot;
    if (T_current == 0.0) return;

    double scale_factor = sqrt(T_target / T_current);
    scale_velocities_gpu<<<num_blocks, BLOCK_SIZE>>>(vx, vy, vz, N, scale_factor);
}