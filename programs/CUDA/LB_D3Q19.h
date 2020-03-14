/**
 * Main module for the CUDA implementation. This module runs on host, and controls the 
 * device.
 */
#include"cuda_functions.h"

class LatticeBoltzmann{
private:
    double w[Q]; int V[D][Q];
    double *f = NULL, *f_new = NULL;
    float *d_f = NULL, *d_f_new = NULL;
public:
    LatticeBoltzmann(void);
    ~LatticeBoltzmann(void);
    int get_1D(int ix, int iy, int iz);
    double rho(int ix, int iy, int iz);
    double Jx(int ix, int iy, int iz);
    double Jy(int ix, int iy, int iz);
    double Jz(int ix, int iy, int iz);
    double Jx_new(int ix, int iy, int iz);
    double Jy_new(int ix, int iy, int iz);
    double Jz_new(int ix, int iy, int iz);
    double f_eq(double rho0, double Ux0, double Uy0, double Uz0, int i);
    void collide(void);
    void propagate(void);
    void initialize(void);
    void impose_fields(void);
    void save(std::string filename);
};
/* Initialize weights and basis vectors, allocate memory for arrays */
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

    // allocate memory on host & device
    f = new double[size]; f_new = new double[size];
    // .
}
/* Free momory on host & device */
LatticeBoltzmann::~LatticeBoltzmann(){
    delete[] f; delete[] f_new;
    cudaFree(d_f); cudaFree(d_f_new);
}
/**
 * Transform from 3D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
 int LatticeBoltzmann::get_1D(int ix, int iy, int iz){
    return ix*x_mult + iy*y_mult + iz*z_mult;
}
// Density
double LatticeBoltzmann::rho(int ix, int iy, int iz){
    double rho = 0; int pos = get_1D(ix, iy, iz);
    for(int i=0; i<Q; i++)
        rho += f[pos + i];
    return rho;
}
// U_x * rho
double LatticeBoltzmann::Jx(int ix, int iy, int iz){
    double J_x = 0; int pos = get_1D(ix, iy, iz);
    for(int i=0; i<Q; i++)
        J_x += f[pos + i] * V[0][i];
    return J_x;
}
// U_y * rho
double LatticeBoltzmann::Jy(int ix, int iy, int iz){
    double J_y = 0; int pos = get_1D(ix, iy, iz);
    for(int j=0; j<Q; j++)
        J_y += f[pos+j] * V[1][j];
    return J_y;
}
// U_z * rho
double LatticeBoltzmann::Jz(int ix, int iy, int iz){
    double J_z = 0; int pos = get_1D(ix, iy, iz);
    for(int k=0; k<Q; k++)
        J_z += f[pos+k] * V[2][k];
    return J_z;
}
// Using f_new
double LatticeBoltzmann::Jx_new(int ix, int iy, int iz){
    double J_x = 0; int pos = get_1D(ix, iy, iz);
    for(int i=0; i<Q; i++)
        J_x += f_new[pos + i] * V[0][i];
    return J_x;
}
// Using f_new
double LatticeBoltzmann::Jy_new(int ix, int iy, int iz){
    double J_y = 0; int pos = get_1D(ix, iy, iz);
    for(int j=0; j<Q; j++)
        J_y += f_new[pos+j] * V[1][j];
    return J_y;
}
// Using f_new
double LatticeBoltzmann::Jz_new(int ix, int iy, int iz){
    double J_z = 0; int pos = get_1D(ix, iy, iz);
    for(int k=0; k<Q; k++)
        J_z += f_new[pos+k] * V[2][k];
    return J_z;
}
// Eq function for fluids
double LatticeBoltzmann::f_eq(double rho0, double Ux0, double Uy0, double Uz0, int i){
    double UdotVi = Ux0*V[0][i] + Uy0*V[1][i] + Uz0*V[2][i];
    double U2 = Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0;
    return rho0*w[i]*(1 + 3*UdotVi + 4.5*UdotVi*UdotVi - 1.5*U2);
}
/* Host */
void LatticeBoltzmann::collide(void){
    d_collide<<<BpG, TpB>>>(d_f, d_f_new);
}
/* Host */
void LatticeBoltzmann::propagate(void){
    d_propagate<<<BpG, TpB>>>(d_f, d_f_new);
}
/* Host */
void LatticeBoltzmann::initialize(void){
    // Load in host
    int nada = 0;
    // Send to device
    // .
}
/* Host */
void LatticeBoltzmann::impose_fields(void){
    int nada = 0;
    // run impose fields kernel
}

void LatticeBoltzmann::save(std::string filename){
    // load data from host
    //.
    // save on file
    int nada = 0;
}
