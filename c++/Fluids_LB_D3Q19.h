#ifndef __AIR_POLLUTION_CPP_FLUIDS_LB_D3Q19_H
#define __AIR_POLLUTION_CPP_FLUIDS_LB_D3Q19_H

#include "constants.h"
#include "LB_D3Q19.h"

#undef f_eq

// eq function for fluids
#define f_eq(rho0, U_Vi, U_2, i) (rho0*w[i]*(1.0 + 3.0*U_Vi + 4.5*U_Vi*U_Vi - 1.5*U_2))

class Fluids : public LatticeBoltzmann{
    private:
        int opposite_of[19] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
    public:
        void collide(void);
        void propagate(void);
        void initialize(double rho0, double Ux0, double Uy0, double Uz0);
        void impose_fields(double v);
        void save(std::string filename, double v);
        void save_2D(std::string filename, int z_pos, double v);
        void print(double v);
};

#endif
