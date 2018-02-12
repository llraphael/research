/* This is the declaration of a collection of functions in
   channel.cpp and function.cpp */

#ifndef FUNC_H
#define FUNC_H

#include <vector>
#include <complex>
#include <stdlib.h>
#include "node.h"



// Complex valued matrix, assumed square.
class ComplexMatrix {

 public:

  std::vector< std::vector<std::complex<double> > > matrix;
  std::vector< std::vector<double> > realMatrix;
  int size;

 public:

  ComplexMatrix() {}

  ComplexMatrix(int matrixsize) 
    {
      matrix = std::vector<std::vector<std::complex<double> > >(matrixsize, std::vector<std::complex<double> >(matrixsize));
      size = matrixsize;
    }
  void writeMatrix(std::vector<std::vector<std::complex<double> > > data)
  {
    int size = data.size();
    
    for(int i=0;i<size;++i)
      for(int j=0;j<size;++j)
	matrix[i][j] = data[i][j];
  }
  void writeMatrix(std::vector<std::complex<double> > values) 
    {
      // Check given values
      if(values.size() != size * size) {
	std::cout << "The size of given values do not match the size of channel!" << std::endl;
	exit(0);
      }
      
      // Start to read values
      int readi = 0;
      for(int i=0;i<size;++i)
	for(int j=0;j<size;++j)
	  matrix[i][j] = values[readi++];
    };

  void writeMatrix(std::complex<double> values[], int length)
    {
      if(length != size * size) {
	std::cout << "The size of given values do not match the size of channel!" << std::endl;
	exit(0);
      }

      // Start to read values
      int readi = 0;
      for(int i=0;i<size;++i)
	for(int j=0;i<size;++j)
	  matrix[i][j] = values[readi++];
    };

  std::complex<double> getEntryVal(int i, int j)
    {
      if(i >= size || j >= size) {
	std::cout << " Given location is out of the matrix!" << std::endl;
	exit(0);
      }

      return matrix[i][j];
    };

  void setEntryVal(int i, int j, std::complex<double> val)
    {
      if(i >= size || j >= size) {
	std::cout << "Given location is out of the matrix!" << std::endl;
	exit(0);
      }
      
      matrix[i][j] = val;
    };

  ComplexMatrix Hermitian()
  {
    ComplexMatrix hmatrix(size);
    for(int rowi=0;rowi<size;++rowi)
      for(int coli=0;coli<size;++coli)
	hmatrix.setEntryVal(rowi, coli, std::conj(matrix[coli][rowi]));
    return hmatrix;
  }

  std::complex<double> determinant() {
    
    std::complex<double> a1 = getEntryVal(0, 0);
    std::complex<double> a2 = getEntryVal(0, 1);
    std::complex<double> a3 = getEntryVal(1, 0);
    std::complex<double> a4 = getEntryVal(1, 1);

    std::complex<double> a5 = a1 * a2 - a3 * a4;
    return a5;
  }

  ComplexMatrix adjacencyMatrix() {
    
    ComplexMatrix adMatrix(2);
    adMatrix.setEntryVal(0, 0, getEntryVal(1, 1));
    adMatrix.setEntryVal(1, 1, getEntryVal(0, 0));
    adMatrix.setEntryVal(0, 1, -getEntryVal(1,0));
    adMatrix.setEntryVal(1, 0, -getEntryVal(0,1));

    return adMatrix;
  }

  ComplexMatrix inverse() {
    
    std::complex<double> divider = determinant();
    
    // Get adjacency matrix
    ComplexMatrix adMatrix = adjacencyMatrix();

    ComplexMatrix inverseMa(2);
    
    for(int i=0;i<2;++i) {
      for(int j=0; j<2; ++j) {
	std::complex<double> ele = getEntryVal(i,j);
	ele = ele / divider;
	inverseMa.setEntryVal(i, j, ele);
      }
    }

    return inverseMa;
    
  }

  
};

std::vector<std::complex<double> > vectorMultiplyMatrix(ComplexMatrix &matrix, std::vector<std::complex<double> > &vec, std::string order);
std::vector<std::complex<double> > vectorMultiplyMatrix(ComplexMatrix &matrix, std::vector<double> &vec, std::string order);


//std::vector<complex> vectorMultiplyMatrix(const MIMOMatrix &matrix, const std::vector<complex> &vec, std::string order);

std::vector<double> generateGaussianNum(double mean, double sigma, int num, int seed);  // Generate a sequence of numbers based on given Gaussian distribution parameters.

void AWGNChannel(std::vector<double> &signal, double sigma, int seed=1000);  // Normal AWGN channel

void complexAWGNChannel(std::vector<std::complex<double> > &signal, double sigma, int seed=1000); // Complex AWGN channel.

/* set of communications function */
std::vector<int> sourceGenerator(int length, double p1, double seed);  // Generate source bits according to p1.

void encoder(std::vector<Coded_info> &codedNodeInfo, std::vector<int> &sysBits, std::vector<int> &codedBits, int startPos = 0);

double getNormalizationFactorPAM(std::vector<double> &symbolset, std::vector<double> &symbolpro);      // Return the normalization factor for PAM modulation.

double getNormalizationFactorQAM(std::vector<double> &symbolset, std::vector<double> &symbolpro);      // Return the normalization factor for QAM modulation.

std::vector<std::complex<double> > modulatorQAM(std::vector<double> &symbols, double normalizefactor);   // QAM modulation.

std::map<double, double> getSymbolSet(std::vector<double> &weightset, const double p1, double zeroProcess=0);   // Employed in the RCM scheme. Determine symbol values based on chosed weight set.

int linearEncoder(std::vector<Coded_info> &codedNodeInfo, std::vector<int> &sysBits, std::vector<double> &symbols, int startPos = 0, double zeroProcess=0);  //Generate coded bits for LDGM codes.

std::map<int, std::vector<std::vector<int> > > getMappingTable(std::vector<double> &weightSet, int zeroProcess=0); // Get mapping table

double gaussianFunc(double signal, double mean, double sigma);  // Gaussian pdf 

int decoderInitializerBPSK(std::vector<double> &channelOutput, std::vector<double> &codeword, std::vector<double> &fade_a, std::vector<double> &fade_b, double sigma, double E_sender = 1); //Decoder initializer for independent sources over Rayleigh fading channel

int sameDecision(std::vector<int> &predecision, std::vector<Sys_info> &sysinfo);  // Determine if the decoder has made same decision for all the bits as the last iteration.



/* Set of MIMO matrix manipulations */
void MIMO22Channel(std::vector<double> &x1, std::vector<double> &x2, ComplexMatrix matrix, double sigma); // For PAM signals.
void MIMO22Channel(std::vector<std::complex<double> > &x1, std::vector<std::complex<double> > &x2, ComplexMatrix matrix, double sigma);    // 2 streams of signals go through the MIMO channel
void MIMO22FastFadingChannel(std::vector<double> &x1, std::vector<double> &x2, std::vector<ComplexMatrix> &HMatrix, double sigma);
void MIMO22FastFadingChannel(std::vector<std::complex<double> > &x1, std::vector<std::complex<double> > &x2, std::vector<ComplexMatrix> &HMatrix, double sigma); 
void MIMOFastFadingChannel(std::vector<std::vector<std::complex<double> > > &x, std::vector<ComplexMatrix> &channel, double sigma);

void SVDPreCoding(std::vector<double> &x1, std::vector<double> &x2, double V_matrix[][2]);  // Precoding for the trasmitted vector X, to VX, where H = U * Z * V^T

void SVDProcessY(std::vector<double> &y1, std::vector<double> &y2, double UT_matrix[][2]);  // Necessary steps to process received signal when using SVD Decoding method.

std::vector<std::vector<double> >  MIMO22MMSEDetector(std::vector<std::complex<double> > y,  const std::vector<double> &symbolset, double norfactor, ComplexMatrix &channel, double sigma, std::string modulationMethod);


/* Set of pdf vector operations */
double complexVecMagnitude(std::vector<std::complex<double> > &vec);   // Compute the magnitude of a complex vector.
int vecNorm(std::vector<double> &origin, int sizeDesire = 0);  //Pdf vector normalization and clip.

int conv(std::vector<double> &vect1, std::vector<double> &vect2, std::vector<double> &res); // Get the pdf of z based on pdfs of x and y, where z = x + y.

double deconvolve(std::vector<double> &sumpdf, std::vector<double> &vpdf, double weight, std::vector<double> &respdf);  // A reasonable method of reverse the convolution process: Get pdf of x based on pdfs of z and y, where z = x + y.

double recurCalPro(std::vector<double> &sumpdf, std::vector<double> &vpdf, double weight, double n);  // Use a recursive way to calculate the value of each position in pdf vector, employed by deconvolve.

int weightPdf(double weight, std::vector<double> &pdf);  // Compute the pdf of constant * x.

double innerProduct(std::vector<double> &v1, std::vector<double> &v2);

std::vector<double> vecProjection(std::vector<double> &a, std::vector<double> &b);  // Compute the projection of vector a onto vector b.

std::vector<std::vector<double> > GramSchmit(std::vector<std::vector<double> > &ma); // Gram-Schmit process.

std::vector<std::vector<std::vector<double> > >  SDMatrixDecomp(ComplexMatrix &ma);

std::vector<std::vector<double> > LSDecoder(std::vector<std::complex<double> > &y, std::vector<double> &symbolset, ComplexMatrix HMatrix, double norfactor, double sigma);

void findSphereCandidate(std::vector<std::vector<double> > &candidates, double d, std::vector<double> &Y, std::vector<double> &symbolset, std::vector<std::vector<double> > &QMatrix, std::vector<std::vector<double> > &RMatrix, double index, std::vector<double> curCandidate, int candiNum);

std::vector<std::vector<double> >  sphereSoftEstimator(std::vector<std::complex<double> > y, std::vector<double> &symbolset, std::vector<std::vector<double> > &candidates, double norfactor, ComplexMatrix &channel, double sigma);

#endif
 
