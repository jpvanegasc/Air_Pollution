#include "LB_D3Q19.h"


void LatticeBoltzmann3D::collide(void){
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++){
                unsigned int pos = get_1D(ix, iy, iz);

                double rho0 = rho(pos);

                double Ux = Jx(pos)/rho0;
                double Uy = Jy(pos)/rho0;
                double Uz = Jz(pos)/rho0;

                double U2 = Ux*Ux + Uy*Uy + Uz*Uz;

                for(int i=0; i<Q; i++){
                    double UdotVi = Ux*V[0][i] + Uy*V[1][i] + Uz*V[2][i];

                    f_new[pos + i] = o_m_o_tau*f[pos + i] + o_tau*f_eq(rho0, UdotVi, U2, i);
                }
            }
}

void LatticeBoltzmann3D::stream(void){
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++){
                unsigned int pos = get_1D(ix, iy, iz);

                for(int i=0; i<Q; i++){
                    unsigned int x = ix + V[0][i];
                    unsigned int y = iy + V[1][i];
                    unsigned int z = iz + V[2][i];

                    if ( // Walls by halfway bounce back
                        (x > Lx-1) || (y > Ly-1) || (z > Lz-1)
                    ){
                        f_new[pos + opposite_of[i]] = f[pos+i];
                    }
                    else{ // Fluid site
                        unsigned int streamed_pos = get_1D(x, y, z);
                        f[streamed_pos + i] = f_new[pos + i];
                    }
                }
            }
}

// Initialize population using the mei et al. scheme
void LatticeBoltzmann3D::initialize(void){
    #define V0 0.0
    #define rho0 1.0

    // Load initial density
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++){
                unsigned int pos = get_1D(ix, iy, iz);

                for(int i=0; i<Q; i++) f[pos + i] = f_eq(rho0, V0, V0, i);
            }
    #undef rho0

    // Collide & propagate just the density
    #define STEPS 100
    for(int t=0; t<STEPS; t++){
        for(unsigned int pos=0; pos<SIZE; pos+=Q){
            double rho0 = rho(pos);

            for(int i=0; i<Q; i++) f_new[pos + i] = o_m_o_tau*f[pos + i] + o_tau*f_eq(rho0, V0, V0, i);
        }
        stream();
    }
    #undef V0
    #undef STEPS
}



/**
 * Macroscopic quantities and utilities
 */

/* Initialize weights and basis vectors, allocate memory for arrays. */
LatticeBoltzmann3D::LatticeBoltzmann3D(void){
    // weights
    w[0] = 1.0/3.0;
    w[1] = w[2] = w[3] = w[4] = w[5] = w[6] = 1.0/18.0;
    w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = 1.0/36.0;
    w[13] = w[14] = w[15] = w[16] = w[17] = w[18] = 1.0/36.0;

    // basis vectors
    V[0][0]=0; //x
    V[1][0]=0; //y
    V[2][0]=0; //z

    V[0][1]=1;   V[0][2]=-1;  V[0][3]=0;   V[0][4]=0;   V[0][5]=0;   V[0][6]=0;
    V[1][1]=0;   V[1][2]=0;   V[1][3]=1;   V[1][4]=-1;  V[1][5]=0;   V[1][6]=0;
    V[2][1]=0;   V[2][2]=0;   V[2][3]=0;   V[2][4]=0;   V[2][5]=1;   V[2][6]=-1;

    V[0][7]=1;   V[0][8]=-1;  V[0][9]=1;   V[0][10]=-1;  V[0][11]=0;   V[0][12]=0;
    V[1][7]=1;   V[1][8]=-1;  V[1][9]=0;   V[1][10]=0;   V[1][11]=1;   V[1][12]=-1;
    V[2][7]=0;   V[2][8]=0;   V[2][9]=1;   V[2][10]=-1;  V[2][11]=1;   V[2][12]=-1;

    V[0][13]=1;   V[0][14]=-1;  V[0][15]=1;   V[0][16]=-1;  V[0][17]=0;   V[0][18]=0;
    V[1][13]=-1;  V[1][14]=1;   V[1][15]=0;   V[1][16]=0;   V[1][17]=1;   V[1][18]=-1;
    V[2][13]=0;   V[2][14]=0;   V[2][15]=-1;  V[2][16]=1;   V[2][17]=-1;  V[2][18]=1;

    // f and f_new
    f = new double[SIZE];

    #if EVOLUTION_ALGORITHM == TWO_STEP
    f_new = new double[SIZE];
    #endif // EVOLUTION_ALGORITHM
}

/* Free memory */
LatticeBoltzmann3D::~LatticeBoltzmann3D(){
    delete[] f; delete[] f_new;
}

/* System density */
double LatticeBoltzmann3D::rho(int position){
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
    r += f[position + 9];
    r += f[position + 10];
    r += f[position + 11];
    r += f[position + 12];

    r += f[position + 13];
    r += f[position + 14];
    r += f[position + 15];
    r += f[position + 16];
    r += f[position + 17];
    r += f[position + 18];

    return r;
}

/* Velocity field in the x axis, times the density. (i.e., U_x * rho) */
double LatticeBoltzmann3D::Jx(int position){
    double J_x = 0;

    J_x += f[position + 1];
    J_x -= f[position + 2];

    J_x += f[position + 7];
    J_x -= f[position + 8];
    J_x += f[position + 9];
    J_x -= f[position + 10];

    J_x += f[position + 13];
    J_x -= f[position + 14];
    J_x += f[position + 15];
    J_x -= f[position + 16];

    return J_x;
}

/* Velocity field in the y axis, times the density. (i.e., U_y * rho) */
double LatticeBoltzmann3D::Jy(int position){
    double J_y = 0;

    J_y += f[position + 3];
    J_y -= f[position + 4];

    J_y += f[position + 7];
    J_y -= f[position + 8];

    J_y += f[position + 11];
    J_y -= f[position + 12];
    J_y -= f[position + 13];
    J_y += f[position + 14];

    J_y += f[position + 17];
    J_y -= f[position + 18];

    return J_y;
}

/* Velocity field in the z axis, times the density. (i.e., U_z * rho) */
double LatticeBoltzmann3D::Jz(int position){
    double J_z = 0;

    J_z += f[position + 5];
    J_z -= f[position + 6];

    J_z += f[position + 9];
    J_z -= f[position + 10];
    J_z += f[position + 11];
    J_z -= f[position + 12];

    J_z -= f[position + 15];
    J_z += f[position + 16];
    J_z -= f[position + 17];
    J_z += f[position + 18];

    return J_z;
}

void LatticeBoltzmann3D::save(std::string filename, double mult){
    std::ofstream file(filename);

    for(int ix=0; ix<Lx; ix+=mult){
        for(int iy=0; iy<Ly; iy+=mult){
            for(int iz=0; iz<Lz; iz+=mult){
                unsigned int pos = get_1D(ix, iy, iz);

                double rho0 = rho(pos);
                double Ux = Jx(pos)/rho0;
                double Uy = Jy(pos)/rho0;
                double Uz = Jz(pos)/rho0;

                file << ix << ',' << iy << ',' << iz << ',' << mult*Ux << ',' << mult*Uy << ',' << mult*Uz << '\n';
            }
            file << '\n';
        }
        file << '\n';
    }

    file.close();
}

// Saves a 2D view from a fixed x position
void LatticeBoltzmann3D::save_2D(std::string filename, int pos, bool x, bool y, bool z, double mult){
    std::ofstream file(filename);

    if (x){
        for(int iy=0; iy<Ly; iy+=mult){
            for(int iz=0; iz<Lz; iz+=mult){
                unsigned int pos = get_1D(pos, iy, iz);

                double rho0 = rho(pos);
                double Uy = Jy(pos)/rho0;
                double Uz = Jz(pos)/rho0;

                file << iy << ',' << iz << ',' << mult*Uy << ',' << mult*Uz << ',' << rho0 << '\n';
            }
            file << '\n';
        }
    }
    else if (y){
        for(int ix=0; ix<Lx; ix+=mult){
            for(int iz=0; iz<Lz; iz+=mult){
                unsigned int pos = get_1D(ix, pos, iz);

                double rho0 = rho(pos);
                double Ux = Jx(pos)/rho0;
                double Uz = Jz(pos)/rho0;

                file << ix << ',' << iz << ',' << mult*Ux << ',' << mult*Uz << ',' << rho0 << '\n';
            }
            file << '\n';
        }
    }
    else if (z){
        for(int ix=0; ix<Lx; ix+=mult){
            for(int iy=0; iy<Ly; iy+=mult){
                unsigned int pos = get_1D(ix, iy, pos);

                double rho0 = rho(pos);
                double Ux = Jx(pos)/rho0;
                double Uy = Jy(pos)/rho0;

                file << ix << ',' << iy << ',' << mult*Ux << ',' << mult*Uy << ',' << rho0 << '\n';
            }
            file << '\n';
        }
    }

    file.close();
}

#if EVOLUTION_ALGORITHM == TWO_STEP
// Using f_new
double LatticeBoltzmann3D::Jx_new(int ix, int iy, int iz){
    double J_x = 0; int pos = get_1D(ix, iy, iz);
    for(int i=0; i<Q; i++)
        J_x += f_new[pos + i] * V[0][i];
    return J_x;
}

// Using f_new
double LatticeBoltzmann3D::Jy_new(int ix, int iy, int iz){
    double J_y = 0; int pos = get_1D(ix, iy, iz);
    for(int j=0; j<Q; j++)
        J_y += f_new[pos+j] * V[1][j];
    return J_y;
}

// Using f_new
double LatticeBoltzmann3D::Jz_new(int ix, int iy, int iz){
    double J_z = 0; int pos = get_1D(ix, iy, iz);
    for(int k=0; k<Q; k++)
        J_z += f_new[pos+k] * V[2][k];
    return J_z;
}
#endif // EVOLUTION_ALGORITHM
