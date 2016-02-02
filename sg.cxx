#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);


void sonderling(const int Nx, const double dx, const double dt, const double omega, cmplx* const f0,cmplx* const f1,const cmplx I,const double xmin);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
        const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = 2.0*xmax/Nx;
	const double dt = dx/1000;
        double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
        const double omega = 0.2;
        const double alpha = sqrt(omega);
        
        const complex <double> I (0,1);

        stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
        cmplx* psi1 = new cmplx[Nx];
        cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
                    
                    sonderling(Nx,dx,dt,omega,psi0,psi1,I,xmin);
                    
                    
                   
                    h=psi0;
                    psi0=psi1;
                    psi1=h;
                    
                   // for(int l=0; l<Nx;l++) psi0[l]=psi1[l];
                    
                    
                    t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
  
  
        delete[] psi0;
        delete[] psi1;

	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}


//--------------------------------------------------------------

void sonderling(const int Nx, const double dx, const double dt, const double omega, cmplx* const f0,cmplx* const f1,const cmplx I, const double xmin){
    
    cmplx* a1 = new cmplx[Nx];
    cmplx* a2 = new cmplx[Nx];
    cmplx* d = new cmplx[Nx];
    cmplx* a1c = new cmplx[Nx];
    cmplx* a2c = new cmplx[Nx];
    cmplx* dc = new cmplx[Nx];
    cmplx* f00 = new cmplx[Nx];
    double V;
    
    
    for(int i=0; i< Nx; i++) a1[i]= -I* dt/(4.0*pow(dx,2));
    for(int i=0; i< Nx; i++) a2[i]= -I* dt/(4.0*pow(dx,2));
    for(int i=0; i< Nx; i++){    
        V = 0.5*pow(omega,2)*pow(xmin+i*dx,2);
        d[i]= 1.0+I*(dt/(2.0*pow(dx,2))+ (dt/2.0)*V);
        }
        
    
    for(int i=0; i< Nx; i++) a1c[i]= I* dt/(4.0*pow(dx,2));
    for(int i=0; i< Nx; i++) a2c[i]= I* dt/(4.0*pow(dx,2));
    for(int i=0; i< Nx; i++){    
        V = 0.5*pow(omega,2)*pow(xmin+i*dx,2);
        dc[i]= 1.0-I*(dt/(2.0*pow(dx,2))+ (dt/2.0)*V);
        }
        
        
        
        
    f00[0] = dc[0]*f0[0]+a1c[0]*f0[1];
    f00[Nx-1] = a2c[Nx-1]*f0[Nx-2]+dc[Nx-1]*f0[Nx-1];
    for(int i=1;i<Nx-1;i++) f00[i] = a2c[i]*f0[i-1]+dc[i]*f0[i]+a1c[i]*f0[i+1];
    
    
    
    
    for(int i=1; i<Nx;i++){
        d[i] -=  a2[i]/d[i-1]*a1[i-1];
        f00[i]-= a2[i]/d[i-1]*f00[i-1];
    }
    
    
    
    f1[Nx-1]=f00[Nx-1]/d[Nx-1];
    for(int i=Nx-2;i>=0;i--) f1[i]=(f00[i]-a1[i]*f1[i+1])/d[i];
        
        
        
    delete[] a1;
    delete[] a2;
    delete[] d;
    delete[] a1c;
    delete[] a2c;
    delete[] dc;
    delete[] f00;

}
