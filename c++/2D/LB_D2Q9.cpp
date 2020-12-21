#include "LB_D2Q9.h"



/**
 * Macroscopic quantities and utilities
 */

/* Initialize weights and basis vectors, allocate memory for arrays. */
LatticeBoltzmann2D::LatticeBoltzmann2D(void){
    // weights
    w[0] = 4.0/9.0;
    w[1] = w[2] = w[3] = w[4] = 1.0/9.0;
    w[5] = w[6] = w[7] = w[8] = 1.0/36.0;

    // basis vectors
    V[0][0]=0; //x
    V[1][0]=0; //y

    V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0;
    V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

    V[0][5]=1;   V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1;
    V[1][5]=1;   V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1;

    // f and f_new
    f = new double[SIZE];

    #if EVOLUTION_ALGORITHM == TWO_STEP
    f_new = new double[SIZE];
    #endif // EVOLUTION_ALGORITHM
}

/* Free memory */
LatticeBoltzmann2D::~LatticeBoltzmann2D(){
    delete[] f; delete[] f_new;
}

/* System density */
double LatticeBoltzmann2D::rho(int position){
    double r = 0;

    r += f[position + 0];

    r += f[position + 1];
    r += f[position + 2];
    r += f[position + 3];
    r += f[position + 4];

    r += f[position + 5];
    r += f[position + 6];
    r += f[position + 7];
    r += f[position + 8];

    return r;
}

/* Velocity field in the x axis, times the density. (i.e., U_x * rho) */
double LatticeBoltzmann2D::Jx(int position){
    double J_x = 0;

    J_x += f[position + 1];

    J_x -= f[position + 3];

    J_x += f[position + 5];
    J_x -= f[position + 6];
    J_x -= f[position + 7];
    J_x += f[position + 8];

    return J_x;
}

/* Velocity field in the y axis, times the density. (i.e., U_y * rho) */
double LatticeBoltzmann2D::Jy(int position){
    double J_y = 0;

    J_y += f[position + 2];

    J_y -= f[position + 4];
    J_y += f[position + 5];
    J_y += f[position + 6];
    J_y -= f[position + 7];
    J_y -= f[position + 8];

    return J_y;
}

#if EVOLUTION_ALGORITHM == TWO_STEP
// Using f_new
double LatticeBoltzmann2D::Jx_new(int ix, int iy){
    double J_x = 0; int pos = get_1D(ix, iy);
    for(int i=0; i<Q; i++)
        J_x += f_new[pos + i] * V[0][i];
    return J_x;
}

// Using f_new
double LatticeBoltzmann2D::Jy_new(int ix, int iy){
    double J_y = 0; int pos = get_1D(ix, iy);
    for(int j=0; j<Q; j++)
        J_y += f_new[pos+j] * V[1][j];
    return J_y;
}

#endif // EVOLUTION_ALGORITHM
