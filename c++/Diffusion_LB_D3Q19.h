#ifndef __AIR_POLLUTION_CPP_DIFFUSION_LB_D3Q19_H
#define __AIR_POLLUTION_CPP_DIFFUSION_LB_D3Q19_H

#include "constants.h"
#include "LB_D3Q19.h"

#define c_s 0.5773502692 // = 1/sqrt(3)
#define c_s2 (c_s*c_s)

#undef f_eq

// eq function for diffusion
#define f_eq(rho0, U_Vi, i) (rho0*w[i]*(1.0 + U_Vi/c_s2))

class Diffusion : public LatticeBoltzmann{
    public:
        void collide(void);
        void propagate(void);
        void initialize(double rho0, double Ux0, double Uy0, double Uz0);
        void save(std::string filename, double v);
        void save_2D(std::string filename, int z_pos, double v);
        void print(double v);
};

#endif
