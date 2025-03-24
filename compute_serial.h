#define COMPUTE_SERIAL

void update_positions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double dt);
void compute_forces(double* x, double* y, double* z, double* fx, double* fy, double* fz, int* type, int N, double& min_sep,char test);
void update_velocities(double* vx, double* vy, double* vz, double* fx, double* fy, double* fz, double* m, int N, double dt);
void apply_boundary_conditions(double* x, double* y, double* z, double* vx, double* vy, double* vz, int N, double Lx, double Ly, double Lz);
double compute_temperature(double* vx, double* vy, double* vz, double* m, int N);
double compute_KE(double* vx, double* vy, double* vz, double* m, int N);
void scale_velocities(double* vx, double* vy, double* vz, double* m, int N, double T_target);
