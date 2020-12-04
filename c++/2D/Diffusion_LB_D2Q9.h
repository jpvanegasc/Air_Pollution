#ifndef __LB_CPP_DIFFUSION_LB_D2Q9_H
#define __LB_CPP_DIFFUSION_LB_D2Q9_H

#include "LB_D2Q9.h"

#define c_s 0.5773502692 // Speed of sound = 1/sqrt(3)
#define c_s2 (c_s*c_s)

#undef f_eq

// eq function for diffusion
#define f_eq(rho0, U_Vi, i) (rho0*w[i]*(1.0 + U_Vi/c_s2))

class Diffusion : public LatticeBoltzmann{
    public:
        void collide(void);
        void propagate(void);
        void impose_fields(double v);
        void initialize(double rho0, double Ux0, double Uy0, double Uz0);
        double detector(int x_pos);
        double sigmax2(void);
        void save(std::string filename, double v);
        void print(double v);
};

#endif
