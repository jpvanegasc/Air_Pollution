#include "Diffusion_LB_D2Q9.h"


void Diffusion::collide(void){
    double rho0, Ux0, Uy0, Uz0; int pos;
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++){
                pos = get_1D(ix, iy);

                rho0 = rho(pos);
                Ux0 = Jx(pos)/rho0; Uy0 = Jy(pos)/rho0;
                double U2 = Ux0*Ux0 + Uy0*Uy0;

                for(int i=0; i<Q; i++){
                    double UdotVi = Ux0*V[0][i] + Uy0*V[1][i];
                    f_new[pos + i] = UmUtau*f[pos + i] + Utau*f_eq(rho0, UdotVi, i);
                }
            }
}

void Diffusion::propagate(void){
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++){
                int pos_new = get_1D(ix, iy);
                    for(int i=0; i<Q; i++){
                        int x_pos = (Lx + ix + V[0][i])%Lx, y_pos = (Ly + iy + V[1][i])%Ly;
                        int pos = get_1D(x_pos, y_pos);
                        f[pos + i] = f_new[pos_new + i];
                    }
        }
}

void Diffusion::initialize(double rho0, double Ux0, double Uy0, double Uz0){
    double mx = Lx/2.0, my = Ly/2.0, sx = Lx/15.0, sy = Ly/15.0;
    double first_part = 1.0/std::sqrt(2*M_PI*sx*sy);

    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++){
            double rho_ = first_part*exp(-0.5*( (((ix-mx)/sx)*((ix-mx)/sx)) + (((iy-my)/sy)*((iy-my)/sy)) ));
        for(int iy=0; iy<Ly; iy++){
                int pos = get_1D(ix, iy);
                double U2 = Ux0*Ux0 + Uy0*Uy0;

                for(int i=0; i<Q; i++){
                    double UdotVi = Ux0*V[0][i] + Uy0*V[1][i];
                    f[pos + i] = f_eq(rho_, UdotVi, i);
                }
            }
    }
}

void Diffusion::save(std::string filename, double v){
    if(v == 0.0) std::cout << "v = 0" << std::endl;
    std::ofstream File(filename); double rho0, Ux0, Uy0, Uz0;
    for(int ix=0; ix<Lx; ix+=4){
        for(int iy=0; iy<Ly; iy+=4){
                int pos = get_1D(ix, iy);

                rho0 = rho(pos);
                //Ux0 = Jx(pos)/rho0; Uy0 = Jy(pos)/rho0; Uz0 = Jz(pos)/rho0;
                // File << ix << '\t' << iy << '\t' << iz << '\t' << 4*(Ux0)/v << '\t'
                // << 4*Uy0/v << '\t' << 4*Uz0/v << '\n';
                File << ix << '\t' << iy << '\t' << '\t' << rho0 << '\n';
        }
        File << '\n';
    }
    File << std::endl;
    File.close();
}

void Diffusion::print(double v){
    if(v == 0.0) std::cout << "v = 0" << std::endl;
    double rho0, Ux0, Uy0, Uz0;
    for(int ix=0; ix<Lx; ix+=4)
        for(int iy=0; iy<Ly; iy+=4){
                int pos = get_1D(ix, iy);

                rho0 = rho(pos);
                Ux0 = Jx(pos)/rho0; Uy0 = Jy(pos)/rho0;
                std::cout << ix << '\t' << iy << '\t' << "\t\t" << 4*(Ux0)/v << '\t'
                << 4*Uy0/v << '\t' << '\n';
            }
    std::cout << std::endl;
}
