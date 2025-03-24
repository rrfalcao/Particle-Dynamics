#define INIT_TESTING
#include <string>
bool far_enough(double x, double y, double z, double* x_arr, double* y_arr, double* z_arr, int init,double R=0.5);
void init_particle(double* x, double* y, double* z, double* vx, double* vy, double* vz, int* type, double* m,int N, double Lx, double Ly, double Lz,double m0, double m1, double& percent_type1,std::string initial_condition); 
void unit_tests(std::string test_flag, double time, double min_separation, double* x, double* y, double* vx, double* vy, int N, double Lx, double Ly,bool end) ;
