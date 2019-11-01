#include<iostream>
#include<fstream>
#include<cmath>
#include"Random64.h"
using namespace std;

const int Lx=256, Ly=64;

const int Q=9;

const double tau=0.55;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
	private:
		double w[Q], f[Lx][Ly][Q], fnew[Lx][Ly][Q]; 
		int V[2][Q]; // V[0][i]=V_ix   V[1][i]=V_iy
	public:
		LatticeBoltzmann(void){
			//Cargar los pesos
			w[0]=4/9.0; 
			w[1]=w[2]=w[3]=w[4]=1/9.0;
			w[5]=w[6]=w[7]=w[8]=1/36.0;
			//Cargar los vectores
			V[0][0]=0;
			V[1][0]=0;

			V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
			V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

			V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;
			V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
		}
		/**
		 * 
		 */
		double rho(int ix, int iy, bool use_new){
				int i; double suma;
			for(suma=0,i=0;i<Q;i++)
				if(use_new) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
			return suma;
		}
		/**
		 * Corriente en dirección x
		 */
		double Jx(int ix, int iy, bool use_new){
			int i; double suma;
			for(suma=0,i=0;i<Q;i++)
				if(use_new) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
			return suma;
		}
		/**
		 * Corriente en dirección y
		 */
		double Jy(int ix, int iy, bool use_new){
			int i; double suma;
			for(suma=0,i=0;i<Q;i++)
				if(use_new) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
			return suma;
		}
		/**
		 * Súper importante
		 */
		double f_eq(double rho0, double Ux0, double Uy0, int i){
			double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
			return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
		}
		/**
		 * 
		 */
		void colisione(void){
			int ix,iy,i; double rho0,Ux0,Uy0;
			//Para cada celda
			for(ix=0;ix<Lx;ix++)
				for(iy=0;iy<Ly;iy++){
					//Calcular las cantidades macroscópicas
					rho0=rho(ix,iy,false);	Ux0=Jx(ix,iy,false)/rho0;	Uy0=Jy(ix,iy,false)/rho0;
					for(i=0;i<Q;i++)
			fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*f_eq(rho0,Ux0,Uy0,i);}
		}
		/**
		 * 
		 */
		void adveccione(void){
			for(int ix=0;ix<Lx;ix++)
				for(int iy=0;iy<Ly;iy++)
					for(int i=0;i<Q;i++)
						f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
			//imponer_fronteras();
		}
		/**
		 * 
		 */
		void inicie(double rho0,double Ux0,double Uy0){
			for(int ix=0;ix<Lx;ix++)
				for(int iy=0;iy<Ly;iy++)
					for(int i=0;i<Q;i++)
			f[ix][iy][i]=f_eq(rho0,Ux0,Uy0,i);
		}
		/**
		 * Campos de velocidades a los que someto el "Fluido". Un túnel de viento sería un campo constante hacía un lado solamente
		 * 
		 * @param t dependencia temporal de los campos que quiero
		 */
		void imponer_campos(int t, double velocidad_ventilador){
			int i,ix,iy; double rho0; int ixc=Lx/8,iyc=Ly/2; int R=Ly/5, R2=R*R;

			for(ix=0;ix<Lx;ix++)
				for(iy=0;iy<Ly;iy++){
					rho0=rho(ix,iy,false);
					//ventilador
					if(ix==0) 
						for(i=0;i<Q;i++) fnew[ix][iy][i]=f_eq(rho0,velocidad_ventilador,0,i);
					//Obstáculo cilíndrico
					else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2)
						for(i=0;i<Q;i++) fnew[ix][iy][i]=f_eq(rho0,0,0,i);
					//Un puntito extra
					//else if(ix==ixc && iy==iyc+R+1)
						for(i=0;i<Q;i++) fnew[ix][iy][i]=f_eq(rho0,0,0,i);
					//Condiciones de frontera
					//else if(ix==Lx-1)
					//	for(i=0;i<Q;i++) fnew[ix][iy][i]=f_eq(rho0,velocidad_ventilador,0,i);
				}
		}
		/**
		 * Genera un archivo con los datos para graficar
		 * 
		 * @param nombre_archivo
		 * @param velocidad_ventilador utilizo la velocidad para determinar el largo de los vectores
		 */
		void imprimase( const char *nombre_archivo, double velocidad_ventilador, bool all);
		void imponer_fronteras(void);
};



int main(){
	LatticeBoltzmann Aire;
	int t, tmax=5000;
	double rho_inicial=1.0, V_ventilador=0.1;

	Aire.inicie(rho_inicial,V_ventilador,0);

	for(t=0; t<tmax; t++){
		Aire.colisione();
		Aire.imponer_campos(t, V_ventilador);
		Aire.adveccione();
	}
	Aire.imprimase("aire.dat", V_ventilador, true);
	
	return 0;
}


void LatticeBoltzmann::imprimase( const char * nombre_archivo, double velocidad_ventilador, bool all){
			ofstream MiArchivo(nombre_archivo); double rho0, Ux0, Uy0;
			if(all==true)
				for(int ix=0; ix<Lx; ix++)
					for(int iy=0; iy<Ly; iy++){
						rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix, iy, true)/rho0;
						MiArchivo<<ix<<'\t'<<iy<<'\t'
							<<(Ux0-velocidad_ventilador)/velocidad_ventilador<<'\t'
							<<Uy0/velocidad_ventilador<<'\n';
					}
			else
				for(int ix=0; ix<Lx; ix+=4)
					for(int iy=0; iy<Ly; iy+=4){
						rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix, iy, true)/rho0;
						MiArchivo<<ix<<'\t'<<iy<<'\t'
							<<4*(Ux0-velocidad_ventilador)/velocidad_ventilador<<'\t'
							<<4*Uy0/velocidad_ventilador<<'\n';
					}
			MiArchivo << endl;
			MiArchivo.close();
}

void LatticeBoltzmann::imponer_fronteras(void){
	for(int i=0; i<Lx; i++)
		for(int j=0; j<Q; j++){ //fronteras x
			//f[i][0][j]=0;
			f[i][Ly-1][j]=0; 
			//fnew[i][0][j]=0;
			fnew[i][Ly-1][j]=0;
		}
	for(int i=0; i<Ly; i++)
		for(int j=0; j<Q; j++){ //fronteras y
			//f[0][i][j]=0; 
			f[Lx-1][i][j]=0; 
			//fnew[0][i][j]=0; 
			fnew[Lx-1][i][j]=0;
		}
}
