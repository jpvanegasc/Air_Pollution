#include<iostream>
#include<cmath>

const int cube = 64;
const int Lx = cube, Ly = cube, Lz = cube;
const int Q = 19;

const double tau = 0.55, Utau = 1.0/tau, UmUtau = 1-Utau;

class LatticeBoltzmann{
    private:
        double w[Q]; int V[3][Q];
        double *f = NULL,   *f_new = NULL;
    public:
        LatticeBoltzmann(void);
        ~LatticeBoltzmann(void);
        // macro
        double rho(void);
        double Jx(void);
        double Jy(void);
        double f_eq(void);
        void collide(void);
        void propagate(void);
        void initialize(void);
        void impose_fields(void);
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
    V[0][13]=0;   V[2][14]=0;   V[2][15]=-1;  V[2][16]=1;   V[2][17]=-1;  V[2][18]=1;
    // f and f_new
    f = new double[Lx*Ly*Lz*Q];
    f_new = new double[Lx*Ly*Lz*Q];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f; delete[] f_new;
}

double LatticeBoltzmann::rho(void){
    return 0;
}

double LatticeBoltzmann::Jx(void){
    return 0;
}

double LatticeBoltzmann::Jy(void){
    return 0;
}

double LatticeBoltzmann::f_eq(void){
    return 0;
}
