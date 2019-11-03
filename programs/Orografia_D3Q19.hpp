#include<iostream>
#include<cmath>
#include<fstream>
#include<string>

const int cube = 8;
const int Lx = cube, Ly = cube, Lz = cube;
const int D = 3, Q = 19;

const double tau = 0.55, Utau = 1.0/tau, UmUtau = 1-Utau;

class LatticeBoltzmann{
    private:
        double w[Q]; int V[D][Q];
        double *f = NULL,   *f_new = NULL;
    public:
        LatticeBoltzmann(void);
        ~LatticeBoltzmann(void);
        // density
        double rho(int ix, int iy, int iz);
        // U_x * rho
        double Jx(int ix, int iy, int iz);
        // U_y * rho
        double Jy(int ix, int iy, int iz);
        // U_z * rho
        double Jz(int ix, int iy, int iz);
        // eq function for fluids
        double f_eq(double rho0, double Ux0, double Uy0, double Uz0, int i);
        void collide(void);
        void propagate(void);
        void initialize(double rho0, double Ux0, double Uy0, double Uz0);
        void impose_fields(double v);
        void print(std::string filename, double v);
        void test(void){ std::cout << "LB working ok" << std::endl;}
};

LatticeBoltzmann::LatticeBoltzmann(void){
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
    f = new double[Lx*Ly*Lz*Q];
    f_new = new double[Lx*Ly*Lz*Q];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f; delete[] f_new;
}

double LatticeBoltzmann::rho(int ix, int iy, int iz){
    double rho = 0; int pos = ix+iy+iz;
    for(int i=0; i<Q; i++){
        rho += f[pos+i];
    }
    return rho;
}

double LatticeBoltzmann::Jx(int ix, int iy, int iz){
    double J_x = 0; int pos = ix+iy+iz;
    for(int i=0; i<Q; i++){
        J_x += f[pos+i] * V[0][i];
    }
    return J_x;
}

double LatticeBoltzmann::Jy(int ix, int iy, int iz){
    double J_y = 0; int pos = ix+iy+iz;
    for(int j=0; j<Q; j++){
        J_y += f[pos+j] * V[1][j];
    }
    return J_y;
}

double LatticeBoltzmann::Jz(int ix, int iy, int iz){
    double J_z = 0; int pos = ix+iy+iz;
    for(int k=0; k<Q; k++){
        J_z += f[pos+k] * V[2][k];
    }
    return J_z;
}

double LatticeBoltzmann::f_eq(double rho0, double Ux0, double Uy0, double Uz0, int i){
    double UdotVi = Ux0*V[0][i] + Uy0*V[1][i] + Uz0*V[2][i];
    double U2 = Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0;
    return rho0*w[i]*(1 + 3*UdotVi + 4.5*UdotVi*UdotVi - 1.5*U2);
}

void LatticeBoltzmann::collide(void){
    double rho0, Ux0, Uy0, Uz0;
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++){
                rho0 = rho(ix, iy, iz);
                Ux0 = Jx(ix, iy, iz)/rho0; Uy0 = Jy(ix, iy, iz)/rho0; Uz0 = Jz(ix, iy, iz)/rho0;
                for(int i=0; i<Q; i++)
                    f_new[ix+iy+iz+i] = UmUtau*f[ix+iy+iz+i] + Utau*f_eq(rho0, Ux0, Uy0, Uz0, i);
            }
}

void LatticeBoltzmann::propagate(void){
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++)
                for(int i=0; i<Q; i++){
                    int x_pos = (ix + V[0][i] + Lx)%Lx;
                    int y_pos = (iy + V[1][i] + Ly)%Ly;
                    int z_pos = (iz + V[2][i] + Lz)%Lz;
                    f[z_pos] = f_new[ix+iy+iz+i];
                }
}

void LatticeBoltzmann::initialize(double rho0, double Ux0, double Uy0, double Uz0){
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lx; iz++)
                for(int i=0; i<Q; i++){
                    f[ix+iy+iz+i] = f_eq(rho0, Ux0, Uy0, Uz0, i);
                }
}

void LatticeBoltzmann::impose_fields(double v){
    double rho0;
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lx; iz++){
                rho0 = rho(ix, iy, iz);
                
                if(ix == 0)
                    for(int i=0; i<Q; i++) f_new[ix+iy+iz+i] = f_eq(rho0, v, 0, 0, i);
            }
}

void LatticeBoltzmann::print(std::string filename, double v){
    std::ofstream File(filename); double rho0, Ux0, Uy0, Uz0;
    for(int ix=0; ix<Lx; ix+=4)
        for(int iy=0; iy<Ly; iy+=4)
            for(int iz=0; iz<Lx; iz+=4){
                rho0 = rho(ix, iy, iz);
                Ux0 = Jx(ix, iy, iz)/rho0; Uy0 = Jy(ix, iy, iz)/rho0; Uz0 = Jz(ix, iy, iz)/rho0;
                File << ix << '\t' << iy << '\t' << iz << '\t' << 4*(Ux0)/v << '\t'
                << 4*Uy0/v << '\t' << 4*Uz0/v << '\n';
            }
    File << std::endl;
    File.close();
}