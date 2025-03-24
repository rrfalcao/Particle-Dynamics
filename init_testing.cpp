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
bool far_enough(double x, double y, double z, double* x_arr, double* y_arr, double* z_arr, int init,double R){ 
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

void init_particle(double* x, double* y, double* z, double* vx, double* vy, double* vz, int* type, double* m,int N, double Lx, double Ly, double Lz,double m0, double m1, double& percent_type1,std::string initial_condition) {
    

    if (initial_condition == "ic-one") {
        // **Test Case 1: One stationary particle**
        x[0] = 10.0; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 0.0; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;
        percent_type1=0.0;
    } 
    else if (initial_condition == "ic-one-vel") {
        // **Test Case 2: One moving particle**
        x[0] = 10.0; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 5.0; vy[0] = 2.0; vz[0] = 1.0;
        type[0] = 0; m[0] = m0;
        percent_type1=0.0;
    } 
    else if (initial_condition == "ic-two") {
        // **Test Case 3: Two bouncing particles**
        x[0] = 8.5; y[0] = 10.0; z[0] = 10.0;
        vx[0] = 0.0; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;

        x[1] = 11.5; y[1] = 10.0; z[1] = 10.0;
        vx[1] = 0.0; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; m[1] = m0;
        percent_type1=0.0;
    }
    else if (initial_condition == "ic-two-pass1") {
        // **Test Case 4: Two passing particles**
        x[0] = 8.5; y[0] = 11.5; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;

        x[1] = 11.5; y[1] = 8.5; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; m[1] = m0;
        percent_type1=0.0;
    }
    else if (initial_condition == "ic-two-pass2") {
        // **Test Case 5: Two passing particles close**
        x[0] = 8.5; y[0] = 11.3; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 0; m[0] = m0;

        x[1] = 11.5; y[1] = 8.7; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 0; m[1] = m0;
        percent_type1=0.0;
    }
    else if (initial_condition == "ic-two-pass3") {
        // **Test Case 6: Two passing heavy particles**
        x[0] = 8.5; y[0] = 11.3; z[0] = 10.0;
        vx[0] = 0.5; vy[0] = 0.0; vz[0] = 0.0;
        type[0] = 1; m[0] = m1;

        x[1] = 11.5; y[1] = 8.7; z[1] = 10.0;
        vx[1] = -0.5; vy[1] = 0.0; vz[1] = 0.0;
        type[1] = 1; m[1] = m1;
        percent_type1=100.0;
    }
    else if (initial_condition == "ic-random") {
        // **Random initialization for N particles**
        int inited=0;
        srand(time(0));
        mt19937 gen(time(0));
        for (int i=0; i<N; i++){
        

            if (i < N * percent_type1 / 100.0) {
                type[i] = 1;
                m[i] = m1;
            } else {
                type[i] = 0;
                m[i] = m0;
            }

            x[i]=uniform_real_distribution<double>(0.0, Lx)(gen);
            y[i]=uniform_real_distribution<double>(0.0, Ly)(gen);
            z[i]=uniform_real_distribution<double>(0.0, Lz)(gen);
            vx[i]=(double)rand()/RAND_MAX - 0.5;
            vy[i]=(double)rand()/RAND_MAX - 0.5;
            vz[i]=(double)rand()/RAND_MAX - 0.5;
            
            if (far_enough(x[i], y[i], z[i], x, y, z, inited,0.5)){
                inited++;
            }
            else{
                i--;
            }

    }
    
}}

void unit_tests(string test_flag, double time, double min_separation, double* x, double* y, double* vx, double* vy, int N, double Lx, double Ly,bool end) {

    
    static bool first_collision1_detected = false;
    static bool first_collision2_detected = false;
    const double tol = 1e-4;

    if (test_flag == "ic-one") {
        // Test 1: Particle must remain stationary
        
        if (abs(vx[0]) > tol || abs(vy[0]) > tol) {
            cout << "FAIL: Particle has moved! Velocity: ("<< vx[0] << ", " << vy[0] << ")\n";
            return;
        }
        
        else if(end){cout << "PASS: Particle remained stationary.\n";}
        
    } 
    
    else if (test_flag == "ic-one-vel") {
        // Test 2: First wall collision
        if (!first_collision1_detected) {
            if (x[0] <= 0.0 || x[0] >= Lx || y[0] <= 0.0 || y[0] >= Ly) {
               
                first_collision1_detected = true;
                cout << "PASS: First wall collision at ("  << x[0] << ", " << y[0]  << ") at t = " << time << "\n";
            }
        }
    } 
    
    else if (test_flag == "ic-two") {
        // Test 3: Minimum separation check
        if(end){
            if (abs(min_separation-1.00022)<tol){
                cout << "PASS: Minimum separation: " << min_separation 
                        << " (Expected: ~1.0022)\n";
            } else {
                cout << "FAIL: Minimum separation: " << min_separation 
                        << " (Expected: ~1.0022)\n";
            }
    } }
    
    else if (test_flag == "ic-two-pass1") {
        // Test 4 & 5: Particle 1 and 2 wall collision + min separation
        if (!first_collision1_detected) {
            if (x[0] <= 0.0 || x[0] >= Lx || y[0] <= 0.0 || y[0] >= Ly) {
                first_collision1_detected = true;
                if (abs(x[0]-20.0)<(tol*100) && abs(y[0]-8.99)<(tol*100)){
                    cout << "PASS: Particle 1 hit wall at (" << x[0] << ", " << y[0] << ")\n";
                } else {
                    cout << "FAIL: Particle 1 hit wall at (" << x[0] << ", " <<y[0]  << ")\n";
                }
            }
        }
        if (!first_collision2_detected) {
            if (x[1] <= 0.0 || x[1] >= Lx || y[1] <= 0.0 || y[1] >= Ly) {
                first_collision2_detected = true;
                if (abs(x[1]-0.0)<(tol*100) && abs(y[1]-11.0)<(tol*1000)){
                    
                    cout << "PASS: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                } else {
        
                    cout << "FAIL: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                }
            }
        }
        if(end){
        if (abs(min_separation-2.89619)<tol){
            cout << "PASS: Minimum separation: " << min_separation 
                        << " (Expected: ~2.89619)\n";
            } else {
            cout << "FAIL: Minimum separation: " << min_separation 
                        << " (Expected: ~2.89619)\n";
            }
            } 
    }
    else if (test_flag == "ic-two-pass2") {
        // Test 5: Particle 1 and 2 wall collision + min separation
        if (!first_collision1_detected) {
            ;
            if (x[0] <= 0.0 || x[0] >= Lx || y[0] <= 0.0 || y[0] >= Ly) {
                first_collision1_detected = true;
                
                if (abs(x[0]-20.0)<(tol*100) && abs(y[0]-15.33)<(tol*100)){
                    cout << "PASS: Particle 1 hit wall at (" << x[0] << ", " << y[0] << ")\n";
                } else {
                    cout << "FAIL: Particle 1 hit wall at (" << x[0] << ", " <<y[0]  << ")\n";
                }
            }
        }
        if (!first_collision2_detected) {
            if (x[1] <= 0.0 || x[1] >= Lx || y[1] <= 0.0 || y[1] >= Ly) {
                first_collision2_detected = true;
                if (abs(x[1]-0.0)<(tol*100) && abs(y[1]-4.67)<(tol*100)){
                    cout << "PASS: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                } else {
                    cout << "FAIL: Particle 2 hit wall at (" << x[1] << ", " << y[1] << ")\n";
                }
            }
        }
        if(end){
        if (abs(min_separation-1.02368)<tol){
                cout << "PASS: Minimum separation: " << min_separation 
                            << " (Expected: ~1.02368)\n";
            } else {
                cout << "FAIL: Minimum separation: " << min_separation 
                            << " (Expected: ~1.02368)\n";
            }
            } 
    }
    
    else if (test_flag == "ic-two-pass3") {
        // Test 6: Minimum separation only
        if(end){
            if (abs(min_separation-3.10166)<tol){
                cout << "PASS: Minimum separation: " << min_separation 
                        << " (Expected: ~3.10166)\n";
            } else {
                cout << "FAIL: Minimum separation: " << min_separation 
                        << " (Expected: ~3.10166)\n";
            }
        }
}
}
