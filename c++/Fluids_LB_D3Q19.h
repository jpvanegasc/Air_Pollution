#include"constants.h"

#ifndef __AIR_POLLUTION_CPP_FLUIDS_LB_D3Q19_H
#define __AIR_POLLUTION_CPP_FLUIDS_LB_D3Q19_H

class Fluids{
    private:
        double w[Q]; int V[D][Q];
        double *f = NULL,   *f_new = NULL;
        int opposite_of[19] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
    public:
        Fluids(void);
        ~Fluids();
        double rho(int position);
        double Jx(int position);
        double Jy(int position);
        double Jz(int position);
        double Jx_new(int ix, int iy, int iz);
        double Jy_new(int ix, int iy, int iz);
        double Jz_new(int ix, int iy, int iz);
        double f_neq(void);
        void collide(void);
        void propagate(void);
        void initialize(double rho0, double Ux0, double Uy0, double Uz0);
        void impose_fields(double v);
        void save(std::string filename, double v);
        void save_2D(std::string filename, int z_pos, double v);
        void print(double v);
};

#endif
