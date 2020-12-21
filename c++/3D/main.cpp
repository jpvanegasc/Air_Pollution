#include<iostream>

#include "LB_D3Q19.h"

#define TMAX 100

int main(void){
    LatticeBoltzmann3D Fluids;

    Fluids.initialize();

    for(int t=0; t<TMAX; t++){
        Fluids.collide();
        Fluids.stream();
    }

    return 0;
}