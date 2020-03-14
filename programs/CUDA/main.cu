#include"LB_D3Q19.h"

int main(void){
    LatticeBoltzmann Bogota;
    int t_max = 10;

    Bogota.initialize(1.0, 1.0, 0, 0);

    for(int t=0; t<t_max; t++){
        Bogota.collide();
        Bogota.impose_fields();
        Bogota.propagate();
    }

    Bogota.save("test.txt");

    return 0;
}
