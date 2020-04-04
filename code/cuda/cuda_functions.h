/**
 * CUDA constants and kernel implementations
 */
#include"constants.h"

// CUDA constants

__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];
__constant__ int d_Vz[Q];

/* 3D to 1D */
__device__ unsigned int d_get1D(unsigned int ix, unsigned int iy, unsigned int iz){
    return ix*x_mult + iy*y_mult + iz*z_mult;
}

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
/* Eq equation for fluids */
__device__ float d_f_eq(double rho0, double Ux0, double Uy0, double Uz0, int i){
    double UdotVi = Ux0*d_Vx[i] + Uy0*d_Vy[i] + Uz0*d_Vz[i];
    double U2 = Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0;
    return rho0*d_w[i]*(1 + 3*UdotVi + 4.5*UdotVi*UdotVi - 1.5*U2);
}

/*================= 
 Evolution Kernels 
 =================*/

__global__ void d_collide(float *d_f, float *d_f_new){
    unsigned int pos = d_get1D(threadIdx.x, threadIdx.y, threadIdx.z);

    float rho0, Jx0, Jy0, Jz0;
    rho0 = d_rho(d_f, pos); Jx0 = d_Jx(d_f, pos); Jy0 = d_Jy(d_f, pos); Jz0 = d_Jz(d_f, pos);

    for(int i=0; i<Q; i++) 
        d_f_new[pos+i] = UmUtau*d_f[pos+i] + Utau*d_f_eq(rho0,Jx0,Jy0,Jz0,i);
}
/* I have a feeling this is unstable. This function accesses memory in a weird way */
__global__ void d_propagate(float *d_f, float *d_f_new){
    unsigned int ix = threadIdx.x, iy = threadIdx.y, iz = threadIdx.x;

    for(int i=0; i<Q; i++){
        // this should be unsigned int
        int x_pos = (ix + d_Vx[i]), y_pos = (iy + d_Vy[i]), z_pos = (iz + d_Vz[i]);
        if(x_pos<0 || x_pos>=Lx) continue;
        if(y_pos<0 || y_pos>=Ly) continue;
        if(z_pos<0 || z_pos>=Lz) continue;
        
        unsigned int pos = d_get1D(ix, iy, iz);
        unsigned int new_pos = d_get1D(x_pos, y_pos, z_pos);
        d_f[new_pos + i] = d_f_new[pos + i];
    }
}

__global__ void d_impose_fields(float *d_f, float *d_f_new){
    unsigned int ix = threadIdx.x, iy = threadIdx.y, iz = threadIdx.x;
    unsigned int pos = d_get1D(ix, iy, iz);
    unsigned int pos_new = d_get1D(threadIdx.x, threadIdx.y, threadIdx.z);

    float rho0, Jx0, Jy0, Jz0;
    rho0 = d_rho(d_f, pos); Jx0 = d_Jx(d_f, pos); Jy0 = d_Jy(d_f, pos); Jz0 = d_Jz(d_f, pos);
}