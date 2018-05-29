#include<iostream>
#include<cmath>
#include <vector>


const double g = -0.5 ; //gamma
const double u = 5*M_PI/4 ; // miu
 
double Fr (double t, const std::vector<double> & y, const std::vector<double> & x,int id);          // Funcion para R
double Fz (double t, const std::vector<double> & y, const std::vector<double> & x ,int id);         // Funcion para z
double ddR (double R, double Z);                                                                    // Segunda derivada de R = R"
double ddZ (double R, double Z);                                                                    // Segunda derivada de Z = Z"
double dphi (double R, double r);                                                                   // Primera dericada de phi = phi'
double eulern(double n0,double h, double v0);                                                       // Euler calcula el 1er paso
double rk4r(double ta, double tb, double dt, std::vector<double> & y,std::vector<double> & x);      // Runge-Kutta calcula 1er paso de R
double rk4z(double ta, double tb, double dt, std::vector<double> & x,std::vector<double> & y);      // Runge-Kutta calcula 1er paso de Z
double Integracion (double h, double (*f)(double, double), double nr, double nz, double phi0);      // Integra phi'
void stormer1(double ta,double tb, double h, std::vector<double> & y,std::vector<double> & x,std::vector<double> & p); // Metodo Stormer


int main(void)
{
  std::cout.precision(6);

  //condiciones iniciales
  
  double dt= 0.0001 ;                             // Cambio en el tiempo
  const double ta = 0.0;                          // Limite inf intervalo
  const double tb = 0.3;                          // Limite sup intervalo
  int n = (tb - ta)/dt ;                          // Numero de pasos (iteraciones)
  double R0 = 0.257453;                           // R inicial
  double Z0 = 0.314687;                           // Z inicial
  double r0 = std::hypot (R0,Z0) ;                // r inicial
  double c = (2*g/R0) + ( R0/std::pow(r0,3) );    // parte del calculo de Q0
  double Q0 = 1 - std::pow( (c) , 2 );            // Q inicial
  double dR0 = (std::sqrt(Q0)) * (std::cos(u)) ;  // R' inicial
  double dZ0 = (std::sqrt(Q0)) * (std::sin(u)) ;  // Z' inicial
  double phi0 = 0.0;                              // phi inicial
  

  // Vectores con posiciones y velocidades iniciales
  
  std::vector<double> r = {R0,dR0};
  std::vector<double> z = {Z0,dZ0};

  
  //Pasos de tiempo con Euler
  
  double nr, nz, np;                           
  nr=eulern(R0,dt,dR0);                       // Primer paso R1
  nz=eulern(Z0,dt,dZ0);                       // Primer paso Z1
  np= Integracion (dt, dphi, nr, nz, phi0);   // Primer paso phi1
  std::vector<double> y = {R0,nr};            // Vector con n-1 y n para R
  std::vector<double> x = {Z0,nz};            // Vector con n-1 y n para Z
  std::vector<double> p = {phi0,np};          // Vector con n-1 y n para phi  

  
  //Pasos de tiempo con Runge kutta
  
  double s,q;
  s=rk4r(ta,tb,dt,r,z);               // Primer paso R1
  q=rk4z(ta,tb,dt,r,z);               // Primer paso Z1
  std::vector<double> ry = {R0,s};    // Vector con n-1 y n para R
  std::vector<double> rx = {Z0,q};    // Vector con n-1 y n para Z


  
  //Metodo de Stormer-Verlet
  
  stormer1(ta,tb,dt,y,x,p);         //Posiciones en cada tiempo
  
  return 0;
}



 
// Funcion para R
double Fr (double t, const std::vector<double> & y, const std::vector<double> & x,int id) 
{
  if (0 == id)
    {
      return y[1];         // R en cada paso
    }
  else if (1 == id) 
    {
      double h;
      h = ddR (y[0],x[0]); // función que calcula Z"
      return h;            // Z en cada paso
    }
  else
    { 
      std::cerr << "It does not exists -> " << id << std::endl;
      exit(1);
    }
}

// Funcion para Z
double Fz (double t, const std::vector<double> & y, const std::vector<double> & x ,int id) 
{
  if (0 == id)
    {
      return x[1];          // R en cada paso
    }
  else if (1 == id) 
    {
      double h;
      h = ddZ (y[0],x[0]);  // función que calcula Z"
      return h;             // Z en cada paso
    }
  else
    {
      std::cerr << "It does not exists -> " << id << std::endl;
      exit(1);
    }
}

//Funcion primera derivada de phi = phi'
double dphi (double R, double r)
{
  return ( 2*g/R + (R/std::pow(r,3)) )* (1/R) ;
}


//Funcion segunda derivada de R = R"
double ddR (double R, double Z)
{
    double r = std::hypot (R,Z) ;                                                         // r 
    double d = ((2*g/R) + (R/std::pow(r,3)));                                             // parte de R"
    double e = ((2*g/std::pow(R,2))+(3*std::pow(R,2)/std::pow(r,5))-(1/std::pow(r,3))) ;  // parte de R"    
    double f = d*e;                                                                       // R"
    
    return f;
}


//Funcion segunda derivada de Z = Z"
double ddZ (double R, double Z)
{
    double r = std::hypot (R,Z) ;                   // r 
    double f = ( (2*g/R) + (R/std::pow(r,3)) );     // parte de Z"
    double g = (3*R*Z/std::pow(r,5));               // parte de Z"
    double h = f*g;                                 // Z"
    
    return h; // Z"
}


//Metodo de Euler para calcular primeros pasos
double eulern(double n0,double h, double v0)
{
  return n0+(h*v0);
}



// Runge-Kutta calcula 1er paso de Z
double  rk4r(double ta, double tb, double dt, std::vector<double> & y,std::vector<double> & x)  
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  double t = ta + dt;
  // k1
  for(int ii = 0; ii < y.size(); ++ii)
    {
      k1[ii] = dt*Fr(t, y,x, ii);
    }
  // k2 aux
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      aux[ii] = y[ii] + k1[ii]/2;
    }
  //k2
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      k2[ii] = dt*Fr(t + dt/2, aux,x, ii);
    }
  // k3 aux
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      aux[ii] = y[ii] + k2[ii]/2;
    }
  //k3
  for(int ii = 0; ii < y.size(); ++ii) 	
    {
      k3[ii] = dt*Fr(t + dt/2, aux,x, ii);
    }
  // k4 aux
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      aux[ii] = y[ii] + k3[ii];
    }
  //k4
  for(int ii = 0; ii < y.size(); ++ii) 	
    {
      k4[ii] = dt*Fr(t + dt, aux,x, ii);
    }
  // write new y
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }
	    
  return y[0] ;
}

// Runge-Kutta calcula 1er paso de Z
double  rk4z(double ta, double tb, double dt, std::vector<double> & x,std::vector<double> & y) 
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  double t = ta + dt;
  // k1
  for(int ii = 0; ii < y.size(); ++ii)
    {
      k1[ii] = dt*Fz(t, y,x, ii);
    }
  // k2 aux
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      aux[ii] = y[ii] + k1[ii]/2;
    }
  //k2
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      k2[ii] = dt*Fz(t + dt/2, aux,x, ii);
    }
  // k3 aux
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      aux[ii] = y[ii] + k2[ii]/2;
    }
  //k3
  for(int ii = 0; ii < y.size(); ++ii) 	
    {
      k3[ii] = dt*Fz(t + dt/2, aux,x, ii);
    }
  // k4 aux
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      aux[ii] = y[ii] + k3[ii];
    }
  //k4
  for(int ii = 0; ii < y.size(); ++ii) 	
    {
      k4[ii] = dt*Fz(t + dt, aux,x, ii);
    }
  // write new y
  for(int ii = 0; ii < y.size(); ++ii) 
    {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }
	    
  return y[0] ;
}


//Metodo de Stormer
void stormer1(double ta,double tb, double h, std::vector<double> & y,std::vector<double> & x,std::vector<double> & p)
{
  const int N =(tb-ta)/h;                                // Numero de pasos (iteraciones)
  std::vector<double>yaux(y.size());                     // Vector Auxiliar de R
  std::vector<double>xaux(x.size());                     // Vector Auxiliar de Z
  
  std::cout<<"T"<<"\t"<<"X"<<"\t"<<"\t"<<"Y"<<"\t"<<"\t"<<"Z"<<std::endl;                           //Titulo
  std::cout<<"0"<<"\t"<<y[0]*std::cos(p[0])<<"\t"<<y[0]*std::sin(p[0])<<"\t"<<x[0]<<std::endl;      //Imprime Condiciones iniciales
  std::cout<< h <<"\t"<<y[1]*std::cos(p[1])<<"\t"<<y[1]*std::sin(p[1])<<"\t"<<x[1]<<std::endl;      //Imprime Primer paso
  
  for(int nt =2; nt<N;++nt)                             // Imprime siguientes pasos en coor Cartesianas
    {
      double t=ta+h*nt;
      double r=std::hypot(y[1],x[1]);
      double k1=0.0,k2=0.0,k3=0.0,k4=0.0;
  
      k1=dphi(y[1],r);
      k2=dphi(y[1]+(h/2),r+(k1/2));
      k3=dphi(y[1]+(h/2),r+(k2/2));
      k4=dphi(y[1]+h,r+k3);
      p[1] =p[1]+((h/6)*(k1+(2*k2)+(2*k3)+k4));     // Imprime Phi (Runge Kutta)

      std::copy(y.begin(),y.end(),yaux.begin());
      std::copy(x.begin(),x.end(),xaux.begin());
      
      y[1]=(h*h*ddR(y[1],x[1]))+(2*y[1])-y[0];      //Stormer para R
      x[1]=(h*h*ddZ(yaux[1],x[1])+ (2*x[1])-x[0]);  //Stormer para Z
      x[0]=xaux[1];
      y[0]=yaux[1];
       
      
      std::cout<<t<<"\t"<<y[1]*std::cos(p[1])<<"\t"<<y[1]*std::sin(p[1])<<"\t"<<x[1]<<std::endl; 
    }
}


// Integracion de phi - Metodo Runge Kutta
double Integracion (double h, double (*f)(double, double), double nr, double nz, double phi0)
{
  double r=std::hypot(nr,nz);
  double k1=0.0,k2=0.0,k3=0.0,k4=0.0;
  
  k1=dphi(nr,r);
  k2=dphi(nr+(h/2),r+(k1/2));
  k3=dphi(nr+(h/2),r+(k2/2));
  k4=dphi(nr+h,r+k3);
  
  double phin=phi0+((h/6)*(k1+(2*k2)+(2*k3)+k4));  // Imprime Phi (Runge Kutta)
  return phin;
}
