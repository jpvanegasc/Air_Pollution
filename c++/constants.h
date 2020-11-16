#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>


#ifndef AIR_POLLUTION_CPP_CONSTANTS_H
#define AIR_POLLUTION_CPP_CONSTANTS_H

// Geometry
#define Lx 256
#define Ly 64
#define Lz 64
#define D 3
#define Q 19

// Obstacle
#define Lx3 (Lx/3)
#define Ly3 (Ly/3)
#define Lz2 (Lz/2)

#define T_Lx3 (2*Lx/3)
#define T_Ly3 (2*Ly/3)

// 3D to 1D
#define size (Lx*Ly*Lz*Q)
#define x_mult (Ly*Lz*Q)
#define y_mult (Lz*Q)
#define z_mult Q

/**
 * Transform from 3D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy, iz) ((ix*x_mult) + (iy*y_mult) + (iz*z_mult))

// LB constants
#define tau 0.55
#define Utau (1.0/tau)
#define UmUtau (1.0-Utau)

#endif
