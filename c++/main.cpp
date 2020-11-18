#include"Fluids_LB_D3Q19.h"
#include"Diffusion_LB_D3Q19.h"

std::string filename(int t);

int main(void){
    // Fluids Bogota;
    Diffusion Bogota;
    int t_max = 500;
    double rho0 = 1.0, v = 0.1;

    Bogota.initialize(rho0, 0, 0, 0);
    Bogota.save_2D("diffusion_initial.txt", Lz/2, v);

    for(int t=0; t<t_max; t++){
        Bogota.collide();
        //Bogota.impose_fields(v);
        Bogota.propagate();
    }

    Bogota.save_2D("diffusion_final.txt", Lz/2, v);

    return 0;
}

std::string filename(int t){
    std::string name; std::stringstream t_s; t_s << t;
    name = t_s.str() + ".txt";
    return name;
}