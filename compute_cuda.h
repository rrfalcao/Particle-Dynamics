#define COMPUTE_CUDA


__global__ void compute_forces_gpu(double* x, double* y, double* z,double* fx, double* fy, double* fz,int* type, int N);
void compute_forces(double* x, double* y, double* z,double* fx, double* fy, double* fz,int* type, int N)
__global__ void update_positions_gpu(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double dt)
void update_positions(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double dt)
__global__ void update_velocities_gpu(double* vx, double* vy, double* vz,double* fx, double* fy, double* fz,double* m, int N, double dt)
void update_velocities(double* vx, double* vy, double* vz,double* fx, double* fy, double* fz,double* m, int N, double dt)
__global__ void apply_boundary_conditions_gpu(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double Lx, double Ly, double Lz)
void apply_boundary_conditions(double* x, double* y, double* z,double* vx, double* vy, double* vz,int N, double Lx, double Ly, double Lz)
__global__ void compute_temperature_gpu(const double* vx, const double* vy, const double* vz,const double* m, double* partial_sum, int N)
double compute_KE(double* vx, double* vy, double* vz, double* m, int N)
__global__ void scale_velocities_gpu(double* vx, double* vy, double* vz, int N, double scale_factor)
void scale_velocities(double* vx, double* vy, double* vz, double* m, int N, double T_target) 