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
const int LX3 = Lx/3, LY3 = Ly/3, LZ2 = Lz/2;
const int T_LX3 = 2*Lx/3, T_LY3 = 2*Ly/3;

#define Lx3 LX3
#define Ly3 LY3
#define Lz2 LZ2

#define T_Lx3 T_LX3
#define T_Ly3 T_LY3

// 3D to 1D
const unsigned int SIZE = Lx*Ly*Lz*Q;
const unsigned int X_MULT = Ly*Lz*Q;
const unsigned int Y_MULT = Lz*Q;

#define size SIZE
#define x_mult X_MULT
#define y_mult Y_MULT
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
const double TAU = 0.55, UTAU = 1.0/TAU, UMUTAU = 1.0-UTAU;

#define tau TAU
#define Utau UTAU
#define UmUtau UMUTAU
