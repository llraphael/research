#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "func.h"

using namespace std;



// Multiplication between a complex vector and a complex matrix by given order
vector<complex<double> > vectorMultiplyMatrix(ComplexMatrix &matrix, vector<complex<double> > &vec, string order)
{
  //Check if they can multiply, assume matrix is a square matrix
  if(vec.size() != matrix.size) {
    cout << " Size does not match, multiplication failed!" << endl;
    exit(0);
  }

  if(order == "left") {
    
    vector<complex<double> > res(matrix.size, complex<double>(0, 0));
    for(int i=0;i<res.size();++i)
      for(int j=0;j<matrix.size;++j)
	res[i] += vec[j] * matrix.getEntryVal(j, i);
    return res;
  }
  
  if(order == "right") {

    vector<complex<double> > res(matrix.size, complex<double>(0, 0));
    for(int i=0;i<res.size();++i)
      for(int j=0;j<matrix.size;++j)
	res[i] += vec[j] * matrix.getEntryVal(i, j);
    return res;
  }
}
// Multiplication between a real vector and a complex matrix by given order
vector<complex<double> > vectorMultiplyMatrix(ComplexMatrix &matrix, vector<double> &vec, string order)
{
  // Transform the double vector into complex number vector.
  vector<complex<double> > newvec;
  for(int i=0;i<vec.size();++i)
    newvec.push_back(complex<double>(vec[i], 0));

  return vectorMultiplyMatrix(matrix, newvec, order);
  
}
    

// Generate needed numbers with gaussian distribution
vector<double> generateGaussianNum(double mean, double sigma, int num, int seed)
{
      vector<double> gauNum(num);

      double N0=2*pow(sigma,2);
  
      int range=10000;
      double u1,u2,v1,v2,s;
      double x,y;

      //srand(time(NULL));
      srand(seed);

      for(int index=0;index<num;++index)
	{
      
	  do{
	
	    u1=double(rand())/RAND_MAX;
	    // u1/=RAND_MAX;
	    u2=double(rand())/RAND_MAX;
	    //	u2/=RAND_MAX;
	    v1=2*u1-1;
	    v2=2*u2-1;
	    s=v1*v1+v2*v2;
	  }while(s>=1);

	  x=sqrt((-2)*log(s)/s)*v1;
	  y=sqrt((-2)*log(s)/s)*v2;

	  x=sqrt(N0/2)*x;
	  y=sqrt(N0/2)*y;
    
	  gauNum[index] = x;
	}
      return gauNum;
}

// Real AWGN channel
void AWGNChannel(vector<double> &signal, double sigma, int seed)
{
  // Get the number of transmitted signals 
  int length = signal.size();
  
  // Generate Gaussian noise
  vector<double> guassianNoise = generateGaussianNum(0, sigma, length, seed);

  // Add noise to the signal
  for(int i=0;i<length;++i)
    signal[i] += guassianNoise[i];

}

// Complex AWGN channel, sigma^2 is the variance of the noise in each dimension
void complexAWGNChannel(vector<complex<double> > &signal, double sigma, int seed)
{
  int length = signal.size();
  
  vector<double> realNoise = generateGaussianNum(0, sigma, length, seed);
  vector<double> imaginaryNoise = generateGaussianNum(0, sigma, length, seed+200);

  for(int i=0;i<length;++i)
    {
      signal[i].real(signal[i].real() + realNoise[i]);
      signal[i].imag(signal[i].imag() + imaginaryNoise[i]);
    }
  
}
