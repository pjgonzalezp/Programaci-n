#include <iostream>
#include <vector>
#include <cmath>

void fill_vector(std::vector<double> & vec);
void fill_matrix(std::vector<double> & Mat);
void print(const std::vector<double> & a);
void printMat(const std::vector<double> & a, int N);
void multiply (const std::vector<double> & A,const std::vector<double> & x, std::vector<double> & b, int C);
// declare la funcion multiply
int main(void)
{
  const int N = 100;
  std::vector<double> A, x, b;
  A.resize(N*N); // modela la matriz con un arreglo unidimensional
  x.resize(N);
  b.resize(N);
  fill_vector(x);
  fill_matrix(A);
  //print(x);
  //printMat(A, N);
  multiply(A, x, b, N);
  std::cout << b[0] << std::endl; // prints : 6.41436
  std::cout << b[N/2] << std::endl; // prints : 0.000385495
  std::cout << b[N-1] << std::endl; // prints : 9.96416e-05
  
return 0;
}
void fill_vector(std::vector<double> & vec)
{
  for(int ii = 0; ii < vec.size(); ++ii) {
    vec[ii] = 2*ii;
  }
}
void fill_matrix(std::vector<double> & Mat)
{
  for(int ii = 0; ii < Mat.size(); ++ii) {
    Mat[ii] = 1.0/((ii+2.0)*(ii+1.0));
  }
}



void print(const std::vector<double> & a)
{
  // print the array
  for(int i = 0; i < a.size(); ++i){
    std::cout << a[i] << std::endl;
  }
}

void printMat (const std::vector<double> & a, int N) //para imprimir la matriz original
{
  std::cout<< "matriz A : \n";
  for (int ii=0; ii<N;++ii){
    for(int jj=0; jj<N; ++jj){

      std::cout << a[ii*N+jj]<<"  ";
    }
    std::cout<<"\n";
  } 
}

void multiply (const std::vector<double> & A,const std::vector<double> & x, std::vector<double> & b, int C)
{
  for(int ii=0; ii<b.size(); ++ii)
    {
      for(int jj = 0; jj < C; ++jj)
	{
	  b[ii] += (A[ii*C+jj])*(x[jj]);
	}
    }
  
}
