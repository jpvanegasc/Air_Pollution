/**
 * LB_D3Q19.h
 */
#ifndef __LB_CPP_LB_D3Q19_H
#define __LB_CPP_LB_D3Q19_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

#include "constants.h"

// Geometry
#define D 3
#define Q 19

// #define C_S ... // speed of sound c_s = 1/std::sqrt(3)
// #define C_S2 0.33333

// 3D to 1D
#define size (Lx*Ly*Lz*Q)
#define x_mult (Ly*Lz*Q)
#define y_mult (Lz*Q)
#define z_mult Q

/**
 * Transform from 3D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy, iz) ((ix*x_mult) + (iy*y_mult) + (iz*z_mult))

// #undef f_eq

// #if LB_TYPE == FLUIDS
// #define f_eq() ()
// #endif

class LatticeBoltzmann3D{
    private:
        double w[Q]; int V[D][Q];
        double *f = NULL;

        // #if EVOLUTION_ALGORITHM == TWO_STEP
        double *f_new = NULL;
        // #endif
    public:
        LatticeBoltzmann3D(void);
        ~LatticeBoltzmann3D();

        double rho(int position);
        double Jx(int position);
        double Jy(int position);
        double Jz(int position);

        double Jx_new(int ix, int iy, int iz);
        double Jy_new(int ix, int iy, int iz);
        double Jz_new(int ix, int iy, int iz);
};

#endif