#include <iostream>
// #include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cblas.h>

using namespace std;
#define F77NAME(x) x##_

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


int main(){


    int N=3;
    const double dt = 0.001;  // Time step
    const double Lx = 20.0, Ly = 20.0, Lz = 20.0;  // Box dimensions
    const double epsilon[2][2] = {{3.0, 15.0}, {15.0, 60.0}}; // LJ epsilon values
    const double sigma[2][2] = {{1.0, 2.0}, {2.0, 3.0}};      // LJ sigma values
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

    init_particle(x, y, z, vx, vy, vz, N, Lx, Ly, Lz);
}