/**
 * Constants & Main libraries used
 */
#include<iostream>
#include<cmath>
#include<fstream>
#include<string>

#include<cuda_runtime.h>

/* Size */
#define Lx 16
#define Ly 16
#define Lz 16
#define Q 19
#define D 3

#define size Lx*Ly*Lz*Q
#define x_mult Ly*Lz*Q
#define y_mult Lz*Q
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

/* LB Constants */
#define W0 (1.0/3.0)
#define C (0.5)
#define TresC2 (0.75)
#define AUX0 (0.5)

#define tau (0.5)
#define Utau (2.0)
#define UmUtau (-1.0)