#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

const int cube = 64;
const int Lx = cube, Ly = cube, Lz = cube;
const int D = 3, Q = 19;

const int size = Lx*Ly*Lz*Q;
const int x_mult = Ly*Lz*Q;
const int y_mult = Lz*Q;
const int z_mult = Q;

const double tau = 0.55, Utau = 1.0/tau, UmUtau = 1-Utau;