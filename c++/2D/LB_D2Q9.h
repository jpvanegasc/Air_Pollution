/**
 * LB_D2Q9.h
 */
#ifndef __LB_CPP_LB_D2Q9_H
#define __LB_CPP_LB_D2Q9_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

#include "constants.h"

// Geometry
#define D 2
#define Q 9

// #define C_S ... // speed of sound c_s = 1/std::sqrt(3)
// #define C_S2 0.33333

// 2D to 1D
#define size (Lx*Ly*Q)
#define x_mult (Ly*Q)
#define y_mult Q

/**
 * Transform from 2D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy) ((ix*x_mult) + (iy*y_mult))

// #undef f_eq

// #if LB_TYPE == FLUIDS
// #define f_eq() ()
// #endif

class LatticeBoltzmann2D{
    private:
        double w[Q]; int V[D][Q];
        double *f = NULL;

        // #if EVOLUTION_ALGORITHM == TWO_STEP
        double *f_new = NULL;
        // #endif
    public:
        LatticeBoltzmann2D(void);
        ~LatticeBoltzmann2D();

        double rho(int position);
        double Jx(int position);
        double Jy(int position);

        double Jx_new(int ix, int iy);
        double Jy_new(int ix, int iy);
};

#endif