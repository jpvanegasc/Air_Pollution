/**
 * Constants & Main libraries used
 */
#include<iostream>
#include<cmath>
#include<fstream>
#include<string>

#include<cuda_runtime.h>

// Geometry
#define Lx 256
#define Ly 64
#define Lz 64
#define D 3
#define Q 19

// Obstacle
#define Lx3 (Lx/3)
#define Ly3 (Ly/3)
#define Lz2 (Lz/2Z2

#define T_Lx3 (2*Lx/3)
#define T_Ly3 (2*Ly/3)

// 3D to 1D
#define size (Lx*Ly*Lz*Q)
#define x_mult (Ly*Lz*Q)
#define y_mult (Lz*Q)
#define z_mult Q

/* Parallelization constants */
int pitch;
// Threads per Block & Blocks per Grid
#define Tx 8
#define Ty 8
#define Tz 8

const int Bx = (Lx+Tx-1)/Tx;
const int By = (Ly+Ty-1)/Ty;
const int Bz = (Lz+Tz-1)/Tz;
// Not sure if this is useful, given that the implementation of the f array turns it into 1D
const dim3 TpB(Tx, Ty, Tz);
const dim3 BpG(Bx, By, Bz);

// LB constants
#define tau 0.55
#define Utau (1.0/tau)
#define UmUtau (1.0-Utau)
