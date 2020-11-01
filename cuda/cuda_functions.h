/**
 * CUDA constants and kernel implementations
 */
#include"constants.h"

/*==================
 CUDA constants & macros
 ==================*/

__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];
__constant__ int d_Vz[Q];

/**
 * Transform from 3D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy, iz) ((ix*x_mult) + (iy*y_mult) + (iz*z_mult))

// Equilibrium function for fluids
#define d_f_eq(rho0, U_Vi, U_2, i) (rho0*d_w[i]*(1.0 + 3.0*U_Vi + 4.5*U_Vi*U_Vi - 1.5*U_2))

/*==================
 Macroscopic fields 
 ===================*/

__device__ float d_rho(float *f0, unsigned int pos){
    float rho = 0;
    for(int i=0; i<Q; i++) rho += f0[pos+i];
    return rho;
}

__device__ float d_Jx(float *f0, int pos){
    float Jx = 0;
    for(int i=0; i<Q; i++) Jx += f0[pos+i]*d_Vx[i];
    return Jx;
}

__device__ float d_Jy(float *f0, int pos){
    float Jy = 0;
    for(int i=0; i<Q; i++) Jy += f0[pos+i]*d_Vy[i];
    return Jy;
}

__device__ float d_Jz(float *f0, int pos){
    float Jz = 0;
    for(int i=0; i<Q; i++) Jz += f0[pos+i]*d_Vz[i];
    return Jz;
}

/*================= 
 Evolution Kernels 
 =================*/

__global__ void d_collide(float *d_f, float *d_f_new){
    unsigned int pos = get1D(threadIdx.x, threadIdx.y, threadIdx.z);

    float rho0 = d_rho(d_f, pos); 
    float Ux0 = d_Jx(d_f, pos)/rho0, Uy0 = d_Jy(d_f, pos)/rho0, Uz0 = d_Jz(d_f, pos)/rho0;

    float U2 = Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0;

    for(int i=0; i<Q; i++){
        float UdotVi = Ux0*d_Vx[i] + Uy0*d_Vy[i] + Uz0*d_Vz[i];
        d_f_new[pos+i] = UmUtau*d_f[pos+i] + Utau*d_f_eq(rho0, UdotVi, U2, i);
    }
}

__global__ void d_propagate(float *d_f, float *d_f_new){
    unsigned int ix = threadIdx.x, iy = threadIdx.y, iz = threadIdx.x;
    unsigned int pos_new = get1D(ix, iy, iz);

    for(int i=0; i<Q; i++){
        unsigned int x_pos = (Lx + ix + d_Vx[i])%Lx, y_pos = (Ly + iy + d_Vy[i])%Ly, z_pos = (Lz + iz + d_Vz[i])%Lz;
        unsigned int pos = get1D(x_pos, y_pos, z_pos);
        d_f[pos + i] = d_f_new[pos_new + i];
    }
}

__global__ void d_impose_fields(float *d_f, float *d_f_new, float v){
    unsigned int ix = threadIdx.x, iy = threadIdx.y, iz = threadIdx.x;

    if(ix==0){
        unsigned int pos = get_1D(ix, iy, iz);
        float rho0 = d_rho(d_f, pos);
        for(int i=0; i<Q; i++){
            float UdotVi = v*d_Vx[i];
            float v2 = v*v;
            d_f_new[pos + i] = d_f_eq(rho0, UdotVi, v2, i);
        }
    }
    else if((ix-Lx/2)*(ix-Lx/2) + (iy-Ly/2)*(iy-Ly/2) + (iz-Lz/2)*(iz-Lz/2) <= Ly*Ly/9.0){
        unsigned int pos = get_1D(ix, iy, iz);
        float rho0 = d_rho(d_f, pos);
        for(int i=0; i<Q; i++) d_f_new[pos + i] = d_f_eq(rho0, 0.0, 0.0, i);
    }
}