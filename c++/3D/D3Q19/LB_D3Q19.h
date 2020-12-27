/**
 * 3D Lattice Boltzmann module.
 * 
 * This module requires a constants.h file, where geometry, relaxation time (tau), LB type and 
 * evolution algorithm are specified. Everything you would need to modify should be there.
 * 
 * constants.h requirements:
 * - Lx, Ly, Lz: Integers defining number of nodes in x, y and z, respectively
 * - tau, o_tau, o_m_o_tou: Doubles defining relaxation time tau, (1/tau) and (1 - 1/tau), respectiveley
 * - EVOLUTION_ALGORITHM: Either 2 (two steps) or 1 (one step). Defines if LB evolution is done in a single
 *      step, with collision and streaming done in a single pass over the f array, or with collision
 *      and streaming having separate passes, respectiveley. Note: Two step evolution requires more 
 *      allocated memory than single step.
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

#define C_S 0.5773502691896258 // speed of sound c_s = 1/sqrt(3)
#define C_S2 0.3333333333333333 // (c_s)^2
#define O_C_S 1.7320508075688772 // (1/c_s)
#define O_C_S2 3.0 // (1/c_s)^2


// Equilibrium function
#undef f_eq

/**
 * Equilibrium function for fluids
 * @param rho: density at position.
 * @param U_Vi: velocity field U dot velocity vector V_i. (U . V_i).
 * @param U2: velocity field U norm squared.
 */
#define f_eq(rho, U_Vi, U2, i) (rho*w[i]*(1.0 + 3.0*U_Vi + 4.5*U_Vi*U_Vi - 1.5*U2))


class LatticeBoltzmann3D{
    private:
        double w[Q]; int V[D][Q];
        double *f = NULL;

        #if EVOLUTION_ALGORITHM == 2
        double *f_new = NULL;
        #endif // EVOLUTION_ALGORITHM

        int opposite_of[Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};

    public:
        LatticeBoltzmann3D(void);
        ~LatticeBoltzmann3D();

        void initialize(void);

        void collide(void);
        void stream(void);
        void evolve(void);

        double rho(int position);
        double Jx(int position);
        double Jy(int position);
        double Jz(int position);

        void save(std::string filename, double mult=1);
        void save_2D(std::string filename, int position, bool x=false, bool y=false, bool z=true, double mult=1);

};


// 3D to 1D
#define SIZE (Lx*Ly*Lz*Q)
#define X_MULT (Ly*Lz*Q)
#define Y_MULT (Lz*Q)
#define Z_MULT Q

/**
 * Transform from 3D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy, iz) ((ix*X_MULT) + (iy*Y_MULT) + (iz*Z_MULT))


// 1D to 3D
/**
 * Transform from 1D notation to 3D notation
 * @return y coordinate
 */
#define get_iz(index) ((index/Z_MULT)%Lz)
/**
 * Transform from 1D notation to 3D notation
 * @return y coordinate
 */
#define get_iy(index) ( ((index - get_iz(index)*Z_MULT)/Y_MULT)%Ly )
/**
 * Transform from 1D notation to 3D notation
 * @return x coordinate
 */
#define get_ix(index) ( (pos - get_iy(pos)*Lz*Z_MULT - get_iz(pos)*Z_MULT)/X_MULT )


#endif // __LB_CPP_LB_D3Q19_H