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



const double dt = 0.001;  // Time step
const double Lx = 20.0, Ly = 20.0, Lz = 20.0;  // Box dimensions
const double epsilon[2][2] = {{3.0, 15.0}, {15.0, 60.0}}; // LJ epsilon values
const double sigma[2][2] = {{1.0, 2.0}, {2.0, 3.0}};      // LJ sigma values