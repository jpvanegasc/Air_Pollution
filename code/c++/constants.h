#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

#define cube 64
#define Lx 256
#define Ly cube
#define Lz cube
#define D 3
#define Q 19

#define size (Lx*Ly*Lz*Q)
#define x_mult (Ly*Lz*Q)
#define y_mult (Lz*Q)
#define z_mult Q

const double tau = 0.55, Utau = 1.0/tau, UmUtau = 1-Utau;