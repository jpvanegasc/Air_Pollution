#include "LB_D3Q7.h"


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

                    f_new[pos + i] = o_m_o_tau*f[pos + i] + o_tau*f_eq(rho0, UdotVi, i);
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

                for(int i=0; i<Q; i++) f[pos + i] = f_eq(rho0, V0, i);
            }
    #undef rho0

    // Collide & propagate just the density
    #define STEPS 100
    for(int t=0; t<STEPS; t++){
        for(unsigned int pos=0; pos<SIZE; pos+=Q){
            double rho0 = rho(pos);

            for(int i=0; i<Q; i++) f_new[pos + i] = o_m_o_tau*f[pos + i] + o_tau*f_eq(rho0, V0, i);
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
    w[0] = 1.0/4.0;
    w[1] = w[2] = w[3] = w[4] = w[5] = w[6] = 1.0/8.0;

    // basis vectors
    V[0][0]=0; //x
    V[1][0]=0; //y
    V[2][0]=0; //z

    V[0][1]=1;   V[0][2]=-1;  V[0][3]=0;   V[0][4]=0;   V[0][5]=0;   V[0][6]=0;
    V[1][1]=0;   V[1][2]=0;   V[1][3]=1;   V[1][4]=-1;  V[1][5]=0;   V[1][6]=0;
    V[2][1]=0;   V[2][2]=0;   V[2][3]=0;   V[2][4]=0;   V[2][5]=1;   V[2][6]=-1;

    // f and f_new
    f = new double[SIZE];

    #if EVOLUTION_ALGORITHM == 2
    f_new = new double[SIZE];
    #endif // EVOLUTION_ALGORITHM
}

/* Free memory */
LatticeBoltzmann3D::~LatticeBoltzmann3D(){
    delete[] f;
    #if EVOLUTION_ALGORITHM == 2
    delete[] f_new;
    #endif // EVOLUTION_ALGORITHM
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

    return r;
}

/* Momentum field in the x axis. (i.e., U_x * rho) */
double LatticeBoltzmann3D::Jx(int position){
    double J_x = 0;

    J_x += f[position + 1];
    J_x -= f[position + 2];

    return J_x;
}

/* Momentum field in the y axis. (i.e., U_y * rho) */
double LatticeBoltzmann3D::Jy(int position){
    double J_y = 0;

    J_y += f[position + 3];
    J_y -= f[position + 4];

    return J_y;
}

/* Momentum field in the z axis. (i.e., U_z * rho) */
double LatticeBoltzmann3D::Jz(int position){
    double J_z = 0;

    J_z += f[position + 5];
    J_z -= f[position + 6];

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
void LatticeBoltzmann3D::save_2D(std::string filename, int position, bool x, bool y, bool z, double mult){
    std::ofstream file(filename);

    if (x){
        for(int iy=0; iy<Ly; iy+=(int)mult){
            for(int iz=0; iz<Lz; iz+=(int)mult){
                unsigned int pos = get_1D(position, iy, iz);

                double rho0 = rho(pos);
                double Uy = Jy(pos)/rho0;
                double Uz = Jz(pos)/rho0;

                file << iy << ',' << iz << ',' << mult*Uy << ',' << mult*Uz << ',' << rho0 << '\n';
            }
            file << '\n';
        }
    }
    else if (y){
        for(int ix=0; ix<Lx; ix+=(int)mult){
            for(int iz=0; iz<Lz; iz+=(int)mult){
                unsigned int pos = get_1D(ix, position, iz);

                double rho0 = rho(pos);
                double Ux = Jx(pos)/rho0;
                double Uz = Jz(pos)/rho0;

                file << ix << ',' << iz << ',' << mult*Ux << ',' << mult*Uz << ',' << rho0 << '\n';
            }
            file << '\n';
        }
    }
    else if (z){
        for(int ix=0; ix<Lx; ix+=(int)mult){
            for(int iy=0; iy<Ly; iy+=(int)mult){
                unsigned int pos = get_1D(ix, iy, position);

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
