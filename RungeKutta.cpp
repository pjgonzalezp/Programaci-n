

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

const double m = 90.0;
const double g = 9.81;
const double Cd = 0.225;


void rk4(double ta, double tb, double dt, std::vector<double> & y);
double f(double t, const std::vector<double> & y, int id);

int main (void)

{
  const double ta = 0.0;
  const double tb = 20.4;
  const double dt = 0.001;
  
  
  std::vector<double> y = {1000.0, 0}; // {x0, v0} 

  rk4(ta,tb, dt, y);

  return 0;
}


void rk4(double ta, double tb, double dt, std::vector<double> & y)
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  const int N = (tb-ta)/dt;
  for (int nt = 0; nt < N; ++nt) {
    double t = ta + dt*nt;
    // k1
    for(int ii = 0; ii < y.size(); ++ii) {
      k1[ii] = dt*f(t, y, ii);
    }
    // k2 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k1[ii]/2;
    }
    //k2
    for(int ii = 0; ii < y.size(); ++ii) {
      k2[ii] = dt*f(t + dt/2, aux, ii);
    }
    // k3 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k2[ii]/2;
    }
    //k3
    for(int ii = 0; ii < y.size(); ++ii) {
      k3[ii] = dt*f(t + dt/2, aux, ii);
    }
    // k4 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k3[ii];
    }
    //k4
    for(int ii = 0; ii < y.size(); ++ii) {
      k4[ii] = dt*f(t + dt, aux, ii);
    }
    // write new y
    for(int ii = 0; ii < y.size(); ++ii) {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    } if (y[0] < 0.0) break;
    
    std::cout << t << "\t" << y[0] << "\t" << y[1] << std::endl;
  }
}


double f(double t, const std::vector<double> & y, int id)
{
  if (0 == id) {
    return y[1];
  }
  else if (1 == id) {
    return -g + (Cd*pow(y[1],2))/m;
    
  }
  else {
    std::cerr << "ERROR!!!!!!!! Id does not exists -> " << id << std::endl;
    exit(1);
  }
}

