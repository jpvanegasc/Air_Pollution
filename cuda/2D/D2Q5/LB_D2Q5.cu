#include"LB_D2Q5.h"

__global__ void d_collide(float *d_f, float *d_f_new){
    unsigned int pos = get_1D(threadIdx.x, threadIdx.y);

    float rho0 = d_rho(d_f, pos); 
    float Ux0 = d_Jx(d_f, pos)/rho0;
    float Uy0 = d_Jy(d_f, pos)/rho0;

    float U2 = Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0;

    for(int i=0; i<Q; i++){
        float UdotVi = Ux0*d_Vx[i] + Uy0*d_Vy[i] + Uz0*d_Vz[i];

        d_f_new[pos + i] = o_m_o_tau*d_f[pos + i] + o_tau*f_eq(rho0, UdotVi, i);
    }
}

__global__ void d_stream(float *d_f, float *d_f_new){
    unsigned int ix = threadIdx.x, iy = threadIdx.y;
    unsigned int pos = get_1D(ix, iy);

    for(int i=0; i<Q; i++){
        unsigned int x = ix + V[0][i];
        unsigned int y = iy + V[1][i];

        if( // Walls by halfway bounce back
            (x > Lx-1) || (y > Ly-1)
        ){
            f_new[pos + d_opposite_of[i]] = f[pos+i];
        }
        else{ // Fluid site
            unsigned int streamed_pos = get_1D(x, y);
            f[streamed_pos + i] = f_new[pos + i];
        }
    }
}

// Initialize population using the mei et al. scheme
void LatticeBoltzmann2D::initialize(void){
    #define V0 0.0
    #define rho0 1.0

    // Load initial density
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++){
                unsigned int pos = get_1D(ix, iy);

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
    #undef STEPS
    #undef V0

    cudaMemcpy(d_f, f, SIZE*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_f_new, f_new, SIZE*sizeof(float), cudaMemcpyHostToDevice);
}

/**
 * Macroscopic quantities and utilities
 */

/* System density */
__device__ float d_rho(float *d_f, unsigned int position){
    float r = 0;

    r += d_f[position + 0];

    r += d_f[position + 1];
    r += d_f[position + 2];
    r += d_f[position + 3];
    r += d_f[position + 4];

    return r;
}

/* Momentum field in the x axis. (i.e., U_x * rho) */
__device__ float d_Jx(float *d_f, unsigned int position){
    float J_x = 0;

    J_x += d_f[position + 1];

    J_x -= d_f[position + 3];

    return J_x;
}

/* Momentum field in the y axis. (i.e., U_y * rho) */
__device__ float d_Jy(float *d_f, unsigned int position){
    float J_y = 0;

    J_y += d_f[position + 2];

    J_y -= d_f[position + 4];

    return J_y;
}

/* Initialize weights and basis vectors, allocate memory for arrays. */
LatticeBoltzmann2D::LatticeBoltzmann2D(void){
    // weights
    w[0] = 1.0/3.0;
    w[1] = w[2] = w[3] = w[4] = 1.0/6.0;

    // basis vectors
    V[0][0]=0; //x
    V[1][0]=0; //y

    V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0;
    V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

    // f and f_new
    f = new float[SIZE];
    cudaMalloc((void **) &d_f, SIZE*sizeof(float));

    #if EVOLUTION_ALGORITHM == 2
    f_new = new float[SIZE];
    cudaMalloc((void **) &d_f_new, SIZE*sizeof(float));
    #endif // EVOLUTION_ALGORITHM

    // LB constants on device
    cudaMemcpyToSymbol(d_w, w, Q*sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_Vx, V[0], Q*sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_Vy, V[1], Q*sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_opposite_of, opposite_of, Q*sizeof(int), 0, cudaMemcpyHostToDevice);
}

/* Free memory */
LatticeBoltzmann2D::~LatticeBoltzmann2D(){
    delete[] f; cudaFree(d_f);
    #if EVOLUTION_ALGORITHM == 2
    delete[] f_new; cudaFree(d_f_new);
    #endif // EVOLUTION_ALGORITHM
}

/* System density */
float LatticeBoltzmann2D::rho(int position){
    float r = 0;

    r += f[position + 0];

    r += f[position + 1];
    r += f[position + 2];
    r += f[position + 3];
    r += f[position + 4];

    return r;
}

/* Momentum field in the x axis. (i.e., U_x * rho) */
float LatticeBoltzmann2D::Jx(int position){
    float J_x = 0;

    J_x += f[position + 1];

    J_x -= f[position + 3];

    return J_x;
}

/* Momentum field in the y axis. (i.e., U_y * rho) */
float LatticeBoltzmann2D::Jy(int position){
    float J_y = 0;

    J_y += f[position + 2];

    J_y -= f[position + 4];

    return J_y;
}

void LatticeBoltzmann2D::save(std::string filename, double mult){
    cudaMemcpy(f, d_f, SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    std::ofstream file(filename);

    for(int ix=0; ix<Lx; ix+=(int)mult){
        for(int iy=0; iy<Ly; iy+=(int)mult){
            int pos = get_1D(ix, iy);

            float rho0 = rho(pos);
            float Ux0 = Jx(pos)/rho0;
            float Uy0 = Jy(pos)/rho0;

            file << ix << ',' << iy << ',' << mult*Ux0 << ',' << mult*Uy0 << ',' <<  rho0 << '\n';
        }
        file << '\n';
    }
    file << std::endl;
    file.close();
}
