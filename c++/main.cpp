#include<iostream>

#include "2D/D2Q9/LB_D2Q9.h"
#define TMAX 100

int main(void){
    LatticeBoltzmann2D Fluids;

    Fluids.initialize();

    for(int t=0; t<TMAX; t++){
        Fluids.collide();
        Fluids.stream();
    }

    return 0;
}