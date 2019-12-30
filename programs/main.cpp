#include<iostream>
#include<cmath>
#include<string>
#include<sstream>
#include"Fluids_D3Q19.hpp"

std::string filename(int t);

int main(void){
    LatticeBoltzmann Bogota;
    int t_max = 100;
    double rho0 = 1.0, v = 10.0;

    Bogota.initialize(rho0, 0, 0, 0);

    for(int t=0; t<t_max; t++){
        Bogota.collide();
        Bogota.impose_fields(v);
        Bogota.propagate();
    }
    Bogota.save(filename(0), 1);
    return 0;
}

std::string filename(int t){
    std::string name; std::stringstream t_s; t_s << t;
    name = t_s.str() + ".txt";
    return name;
}