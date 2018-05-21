#include <iostream>
#include <vector>
#include <cmath>

void fill_vector (std::vector<double> & a);
void print(const std::vector<double> & a);
double vector_p_norm (const std::vector<double> &a, int p);

// Escriba aca las declaraciones de las funciones
int main(void)
{
  const int N = 57;
  std::vector<double> data;
  data.resize(N);
  fill_vector(data);
  //print(data);

  std::cout << vector_p_norm(data, 1 ) << std::endl; // prints: 53.9967
  std::cout << vector_p_norm(data, 2 ) << std::endl; // prints: 7.22653
  std::cout << vector_p_norm(data, 3 ) << std::endl; // prints: 3.69935
  std::cout << vector_p_norm(data, 10) << std::endl; // prints: 1.45413
  std::cout << vector_p_norm(data, 90) << std::endl; // prints: 1.02753
 
 
  return 0;
}
// Escriba aca las implementaciones de las funciones

void fill_vector (std::vector<double> & a)
{
  for(int i=0; i<a.size(); ++i)
    {
      // a[i]=2*i;
      a[i]=(2.0*i)/((2.0*i)+1);
    }
}

void print(const std::vector<double> & a)
{
  // print the array
  for(int i = 0; i < a.size(); ++i){
    std::cout << a[i] << std::endl;
  }
}

double vector_p_norm (const std::vector<double> &a, int p)
{
  double  norma=0.0, b=0.0, sum=0.0;
  for (int i=0; i<a.size(); ++i)
    {
      b=std::pow(a[i],p);
      sum += b;
    }
  norma=std::pow(sum, 1.0/p);

  return norma;
}
