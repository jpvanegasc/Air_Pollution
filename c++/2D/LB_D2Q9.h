/**
 * 2D Lattice Boltzmann module.
 * 
 * This module requires a constants.h file, where geometry, relaxation time (tau), LB type and 
 * evolution algorithm are specified. Everything you would need to modify should be there.
 * 
 * constants.h requirements:
 * - Lx, Ly: Integers defining number of nodes in x and y, respectively
 * - tau, o_tau, o_m_o_tou: Doubles defining relaxation time tau, (1/tau) and (1 - 1/tau), respectiveley
 * - LB_TYPE: Either 1 (fluids) or 2 (diffusion). Defines if LB evaluates fluids or diffusion.
 * - EVOLUTION_ALGORITHM: Either 2 (two steps) or 1 (one step). Defines if LB evolution is done in a single
 *      step, with collision and streaming done in a single pass over the f array, or with collision
 *      and streaming having separate passes, respectiveley. Note: Two step evolution requires more 
 *      allocated memory than single step.
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

#define C_S 0.5773502691896258 // speed of sound c_s = 1/sqrt(3)
#define C_S2 0.3333333333333333 // (c_s)^2
#define O_C_S 1.7320508075688772 // (1/c_s)
#define O_C_S2 3.0 // (1/c_s)^2


// 2D to 1D
#define SIZE (Lx*Ly*Q)
#define X_MULT (Ly*Q)
#define Y_MULT Q

/**
 * Transform from 2D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy) ((ix*X_MULT) + (iy*Y_MULT))


// 1D to 2D
/**
 * Transform from 1D notation to 2D notation
 * @return x coordinate
 */
#define get_ix(index) (index/(X_MULT))
/**
 * Transform from 1D notation to 2D notation
 * @return y coordinate
 */
#define get_iy(index) ((index/Y_MULT)%Ly)


// Equilibrium function
#undef f_eq

#if LB_TYPE == 1
/**
 * Equilibrium function for fluids
 * @param rho: density at position.
 * @param U_Vi: velocity field U dot velocity vector V_i. (U . V_i).
 * @param U2: velocity field U norm squared.
 */
#define f_eq(rho, U_Vi, U2, i) (rho*w[i]*(1.0 + 3.0*U_Vi + 4.5*U_Vi*U_Vi - 1.5*U2))

#elif LB_TYPE == 2
/**
 * Equilibrium function for diffusion
 * @param rho: density at position.
 * @param U_Vi: velocity field U dot velocity vector V_i. (U . V_i).
 * @param U2: velocity field U norm squared. (Not used for diffusion, you can pass any value here, so there's no need to actually calculate it).
 */
#define f_eq(rho, U_Vi, U2, i) (rho*w[i]*(1.0 + U_Vi*O_C_S2))

#endif // LB_TYPE

class LatticeBoltzmann2D{
    private:
        double w[Q]; int V[D][Q];
        double *f = NULL;

        #if EVOLUTION_ALGORITHM == 2
        double *f_new = NULL;
        #endif // EVOLUTION_ALGORITHM

        int opposite_of[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    public:
        LatticeBoltzmann2D(void);
        ~LatticeBoltzmann2D();

        void initialize(void);

        void collide(void);
        void stream(void);
        void evolve(void);

        double rho(int position);
        double Jx(int position);
        double Jy(int position);

        void save(std::string filename, double mult=1);

        double Jx_new(int ix, int iy);
        double Jy_new(int ix, int iy);
};

#endif // __LB_CPP_LB_D2Q9_H