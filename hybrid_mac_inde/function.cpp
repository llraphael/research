#include<iostream>
#include<fstream>
#include<cmath>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "func.h"

using namespace std;

// Binary digits source generator
vector<int> sourceGenerator(int length, double p1, double seed)
{
  vector<int> sourcebits(length);

   srand(seed);
    
    for(int index=0;index<length;++index)
      {
	int random = rand()%1000000;
	if(random<p1*1000000)
	  sourcebits[index] = 1;
	else
	  sourcebits[index] = 0;
      }
    return sourcebits;
}

//ldgm encoder
void encoder(vector<Coded_info> &codedNodeInfo, vector<int> &sysBits, vector<int> &codedBits,int startPos)
{
  int codedNo = codedBits.size();    
  
  for(int i=0;i<codedNo;++i)
    {  
      // Current coded node degree
      int codeddegree = codedNodeInfo[i].nodeData.degree;
      int sysbitSum=0;
      for(int index=0;index!=codeddegree;++index)
	{
	  int no=codedNodeInfo[i+startPos].nodeData.neighbourNum[index];
	  if(sysBits[no]==1)
	    sysbitSum+=1;
	  else
	    sysbitSum+=0;
	}
      codedBits[i]=sysbitSum%2;
      if(codedBits[i]==1)
	;
      else
	codedBits[i]=0;
      }
}

//RCM Encoder
int linearEncoder(vector<Coded_info> &codedNodeInfo, vector<int> &sysBits, vector<double> &symbols,int startPos, double zeroProcess)
{
  int codedNo = symbols.size();    

  if(codedNodeInfo.size() - codedNo != startPos) {
    cout<<" The length of coded bits is not right!"<<endl;
    return -1;
  }
  
  #pragma omp parallel for
  for(int i=0;i<codedNo;++i)
    {  
      double sysbitSum=0;
      int codeddegree = codedNodeInfo[i+startPos].nodeData.degree;
      for(int index=0;index!=codeddegree;++index)
	{
	  int no=codedNodeInfo[i+startPos].nodeData.neighbourNum[index];
	  if(zeroProcess == -1)
	    sysbitSum += codedNodeInfo[i+startPos].weightset[index]*(2*sysBits[no]-1);
	  else
	    sysbitSum += codedNodeInfo[i+startPos].weightset[index]*sysBits[no];

	}
      symbols[i]=sysbitSum;
    }
}

double getNormalizationFactorPAM(vector<double> &symbolset, vector<double> &symbolpro)
{
  // Normalize the energy
  int symbolnum = symbolset.size();

  double averEnergy = 0;
  for(int i=0;i<symbolnum;++i)
    averEnergy += symbolpro[i] * pow(symbolset[i], 2);

  double normalizedfactor = sqrt( 1 / averEnergy);
  return normalizedfactor;
}

double getNormalizationFactorQAM(vector<double> &symbolset, vector<double> &symbolpro)
{
  
  // Get the constellation and normalize the energy.
  int symbolnum = symbolset.size();

  double averener = 0;
  for(int xmag = 0; xmag < symbolnum; ++xmag)
    {
      for(int ymag = 0; ymag < symbolnum; ++ymag)
	{
	  double signalpro = symbolpro[xmag] * symbolpro[ymag];
	  averener += signalpro * (pow(symbolset[xmag], 2) + pow(symbolset[ymag], 2));
	}
    }
  
  double normalizedfactor = 1 / averener;                 // The value of d^2
  normalizedfactor = sqrt(normalizedfactor);

  return normalizedfactor;
}

vector<complex<double> > modulatorQAM(vector<double> &symbols, double norfactor)
{  
  if(symbols.size() % 2 != 0) {
    cout << " The number of symbols modulated should be even!" << endl;
    exit(0);
  }
  vector<complex<double> > channelsignals(symbols.size() / 2);
  // Group every two symbols into one wireless signals
  for(int i=0;i<symbols.size();i=i+2)
    channelsignals[i/2] = complex<double>(symbols[i]*norfactor, symbols[i+1]*norfactor);

  return channelsignals;
  
}

map<double, double> getSymbolSet(vector<double> &weightset, const double p1, double zeroProcess)
{

  int weightsize = weightset.size();

  vector<double> sysbits(weightsize, 0);
  vector<double> symbolset;
  vector<double> symbolpro;
  
  int i = 0;
  while(i<pow(2,weightsize)) {
    
    double symbol = 0;
    for(int sysindex=0;sysindex<weightsize;++sysindex) {
      if(zeroProcess == -1)
	symbol += (2*sysbits[sysindex]-1) * weightset[sysindex];
      else
	symbol += sysbits[sysindex] * weightset[sysindex];
    }

    symbolset.push_back(symbol);
    double pro = 1;
    for(int sysindex=0;sysindex<weightsize;++sysindex)
      {
	if(sysbits[sysindex] == 0)
	  pro *= (1-p1);
	else
	  pro *= p1;
      }
    symbolpro.push_back(pro);

    int j = 0;
    ++sysbits[j];
    while(sysbits[j] == 2 && j<weightsize ) {
      if(j == weightsize-1)
	sysbits[j] = 0;
      else {
	sysbits[j] = 0;
	++sysbits[++j];
      }
    }

    ++i;
  }

  vector<double> symbol, procal;
  for(i=0;i<symbolset.size();++i)
    {
      int j = 0;
      for(j;j<symbol.size();++j)
	{
	  if(symbolset[i] == symbol[j]) {
	    procal[j] += symbolpro[i];
	    break;
	  }
	  else 
	    continue;
	}

      if(j == symbol.size()) {                   // There is no symbol with the same value in the set
	symbol.push_back(symbolset[i]);
	procal.push_back(symbolpro[i]);
      }
    }
  
  vecNorm(procal);                               // Normalize the distribution of symbols.

  map<double, double> symboldata;
  for(i=0;i<symbol.size();++i)
    symboldata.insert( pair<double, double>(symbol[i], procal[i]) );

  return symboldata;


}

map<int, vector<vector<int> > > getMappingTable(vector<double> &weightSet, int zeroProcess) {
  int weightsize = weightSet.size();

  vector<int> sysbits(weightsize, 0);
  map<int, vector<vector<int> > > table;
 
  
  int i = 0;
  while(i<pow(2,weightsize)) {
    
    int symbol = 0;
    for(int sysindex=0;sysindex<weightsize;++sysindex) {
      if(zeroProcess == -1)
	symbol += (2*sysbits[sysindex]-1) * weightSet[sysindex];
      else
	symbol += sysbits[sysindex] * weightSet[sysindex];
    }
    
    if(table.count(symbol))
      table[symbol].push_back(sysbits);
    else {
      table[symbol] = vector<vector<int> >();
      table[symbol].push_back(sysbits);
    }

    int j = 0;
    ++sysbits[j];
    while(sysbits[j] == 2 && j<weightsize ) {
      if(j == weightsize-1)
	sysbits[j] = 0;
      else {
	sysbits[j] = 0;
	++sysbits[++j];
      }
    }

    ++i;
  }

  return table;
}

double gaussianFunc(double signal, double mean, double sigma)
{

  double res =  1 / sqrt(2*pi*pow(sigma,2)) *	    \
    exp(-pow(signal-mean, 2) / (2*pow(sigma,2)));

  return res;

}

//Add one to current number based on the given base system), leftmost position is lowest. if overflow, become all zero. 
void addOneToNumber(vector<int> &num, int base)
{
  int len = num.size();
  int carry = 1;

  for(int i=0;i<len;++i)
    {
      if(carry == 0)
	break;
      else if(carry == 1 && num[i] == base - 1)
	num[i] = 0;
      else {
	++num[i];
	carry = 0;
      }
    }
}

// Decoder initialization for independent sources over Rayleigh fading channel
int decoderInitializerBPSK(vector<double> &channelOutput, vector<double> &codeword, vector<double> &fade_a, vector<double> &fade_b, double sigma, double E_sender)
{
  //examine the input
  if(channelOutput.size()!=codeword.size()) {
    cout<<"input codeword not right!"<<endl;
    return -1;
  }
  else if (fade_a.size()!=fade_b.size()) {
    cout<<"fading coefficient size does not match!"<<endl;
    return -1;
  }


  int symbolnum = channelOutput.size();

  #pragma omp parallel for
  for(int j=0;j<symbolnum;++j)
    {
      double v0;
      double y=codeword[j];

      double h1 = fade_a[j];
      double h2 = fade_b[j];

      v0=log((exp(-pow(y-(-h1-h2)*sqrt(E_sender),2)/(2*pow(sigma,2)))	\
	      +exp(-pow(y-(-h1+h2)*sqrt(E_sender),2)/(2*pow(sigma,2))))/ \
	     (exp(-pow(y-(h1+h2)*sqrt(E_sender),2)/(2*pow(sigma,2)))	\
	      +exp(-pow(y-(h1-h2)*sqrt(E_sender),2)/(2*pow(sigma,2)))));       
	
      channelOutput[j]=v0;
    }
}

int sameDecision(vector<int> &predecision, vector<Sys_info> &sysinfo)
{

  int sysno = predecision.size();

  int flag = 1;

  for(int i=0;i<sysno;++i)
    {
      if(predecision[i] != sysinfo[i].decision) {
	flag = 0;
	break;
      }
    }
  
  #pragma omp parallel for
  for(int i=0;i<sysno;++i)
    predecision[i] = sysinfo[i].decision;

  return flag;

}

  
// Based on the CSI and symbol candidate set, compute soft estimations of the actual symbols from each size-2 received
// vectors. It returns the pdf for each actual symbols. 
vector<vector<double> >  MIMO22MMSEDetector(vector<complex<double> > y, const vector<double> &symbolset, double norfactor, ComplexMatrix &channel, double sigma, string modulationMethod)
{
  // The number of symbols to be estimated for each received vector based on the modulation method.
  int symbolnum;
  if(modulationMethod == "PAM")
    symbolnum = 2;
  if(modulationMethod == "QAM")
    symbolnum = 4;

  // Get the number of possible values for each symbol
  // and construct a traverse table
  int possiblenum = symbolset.size();
  vector<int> traversetable(symbolnum, 0);

  // Compute the pdf for each possible x based on received vector y = H * s + n.
  vector<double> originalSymbols(symbolnum, 0);
  map<vector<double>, double> prodata;
  double possibleOrigVecNum = pow(possiblenum, symbolnum);

  // If QAM modulation is adopted.
  if(modulationMethod == "QAM") {
    for(int vecindex=0;vecindex<possibleOrigVecNum;++vecindex)
      {
	for(int i=0;i<symbolnum;++i)
	  {
	    int valindicator = traversetable[i];
	    originalSymbols[i] = symbolset[valindicator];
	  }
      
	// Construct channel symbols from original symbols.
	vector<complex<double> > transmittedSymbols(2);
	transmittedSymbols[0] = complex<double>(originalSymbols[0]*norfactor, originalSymbols[1]*norfactor);
	transmittedSymbols[1] = complex<double>(originalSymbols[2]*norfactor, originalSymbols[3]*norfactor);
      
	// Compute the pdf: p(y = map(x)) = exp[ -1/(2*sigma^2) * ||y - H*s||^2 ] / (2*pi*sigma^2)^2
	vector<complex<double> > mHs = vectorMultiplyMatrix(channel, transmittedSymbols, "right");
	mHs[0] = -mHs[0];
	mHs[1] = -mHs[1];

	vector<complex<double> > yhs(2);
	yhs[0] = y[0] + mHs[0];
	yhs[1] = y[1] + mHs[1];
	double yhsSquare = norm(yhs[0]) + norm(yhs[1]);
	double pro = exp( (-1) * yhsSquare  / (2 * pow(sigma,2)) ) / pow(2*pi*pow(sigma,2), 2);
	  
	// Store the original symbol vector and corresponding pdf to the prodata
	prodata[originalSymbols] = pro;
	addOneToNumber(traversetable, possiblenum);
      }
  }

  // If PAM modulation is adopted.
  if(modulationMethod == "PAM") {
    
    vector<vector<double> > realMatrix(2, vector<double>(2));
    realMatrix[0][0] = channel.getEntryVal(0,0).real();
    realMatrix[0][1] = channel.getEntryVal(0,1).real();
    realMatrix[1][0] = channel.getEntryVal(1,0).real();
    realMatrix[1][1] = channel.getEntryVal(1,1).real();
    
    for(int vecindex=0;vecindex<possibleOrigVecNum;++vecindex)
      {
	// One possible vector of original symbols transmitted.
	for(int i=0;i<symbolnum;++i)
	  {
	    int valindicator = traversetable[i];
	    originalSymbols[i] = symbolset[valindicator];
	  }
	// Construct channel symbols from original symbols.
	vector<double> transmittedSymbols(2);
	transmittedSymbols[0] = originalSymbols[0] * norfactor;
	transmittedSymbols[1] = originalSymbols[1] * norfactor;

	// Compute the pdf:

	// Calculate || Y - H s || ^ 2
	vector<double> mHs(2);
	mHs[0] = (realMatrix[0][0] * transmittedSymbols[0] + realMatrix[0][1] * transmittedSymbols[1]);
	mHs[1] = (realMatrix[1][0] * transmittedSymbols[0] + realMatrix[1][1] * transmittedSymbols[1]);

	vector<double> yhs(2);
	yhs[0] = y[0].real() - mHs[0];
	yhs[1] = y[1].real() - mHs[1];
	double yhsSquare = pow(yhs[0], 2) + pow(yhs[1], 2);
	double pro = exp( (-1) * yhsSquare / (2 * pow(sigma, 2)) ) / pow(2 * pi * pow(sigma,2), 2);
	
	// Store the probability of this possible vector.
	prodata[originalSymbols] = pro;
	addOneToNumber(traversetable, possiblenum);
      }
  }
  // Organize the data to get the possibility for each value of each symbol in the original symbol vectors.
  // Each position of the original symbol vector
  
  //Determine the size of the pdf vector for each symbol in the vector
  int pdfsize;
  if(possiblenum == 2)  // Digital bit
    pdfsize = 3;
  else
    pdfsize = 2 * symbolset[symbolset.size()-1] + 1;   // Pdf size depends on the largest value in the set.
  int zeropoint = (pdfsize - 1) / 2;

  vector<double> pdfinitial(pdfsize, 0);
  vector<vector<double> > pdfdata(symbolnum, pdfinitial);
  traversetable = vector<int>(symbolnum, 0);

  // Iterate through all the possible sent vectors
  for(int vecindex=0;vecindex<possibleOrigVecNum;++vecindex)
    {
      // Construct the original vector for current case
      vector<double> curvec(symbolnum, 0);
      for(int symi=0;symi<symbolnum;++symi)
	{
	  int valindicator = traversetable[symi];
	  curvec[symi] = symbolset[valindicator];
	}
      double pro = prodata[curvec];

      // Process each each position in the symbol vector
      for(int symi=0;symi<symbolnum;++symi)
	{
	  int candivalue = curvec[symi];
	  int posOfCandi = zeropoint + candivalue;
	  pdfdata[symi][posOfCandi] += pro;
	}
      
      addOneToNumber(traversetable, possiblenum);
    }

  // Normalization
  for(int symi=0;symi<symbolnum;++symi)
    vecNorm(pdfdata[symi]);
  
  return pdfdata;
}


// Transmit bits through the 2 * 2 MIMO channel
void MIMO22Channel(vector<complex<double> > &x1, vector<complex<double> > &x2, ComplexMatrix channel, double sigma)
{
  // Check if two antennas transmit the same number of symbols.
  if(x1.size() != x2.size()) {
    cout << "There should be the same number of symbols on two antennas" << endl;
    exit(0);
  }

  //Y1 = h11*x1 + h12*x2 + n1, Y2 = h21*x1 + h22*x2 + n2
  vector<complex<double> > y1(x1.size());
  vector<complex<double> > y2(x1.size());

  //Fading first
  for(int i=0;i<x1.size();++i)
    {
      vector<complex<double> > tVector;
      tVector.push_back(x1[i]);
      tVector.push_back(x2[i]);

      // y = H * s
      vector<complex<double> > rVector = vectorMultiplyMatrix(channel, tVector, "right");
      y1[i] = rVector[0];
      y2[i] = rVector[1];
    }

  //Add AWGN noise
  int seed = 1000;
  complexAWGNChannel(y1, sigma, seed);
  complexAWGNChannel(y2, sigma, seed+10000);

  x1 = y1;
  x2 = y2;
}


void MIMO22Channel(vector<double> &x1, vector<double> &x2, ComplexMatrix channel, double sigma)
{
   // Check if two antennas transmit the same number of symbols.
  if(x1.size() != x2.size()) {
    cout << "There should be the same number of symbols on two antennas" << endl;
    exit(0);
  }

  //Y1 = h11*x1 + h12*x2 + n1, Y2 = h21*x1 + h22*x2 + n2
  vector<double> y1(x1.size());
  vector<double> y2(x1.size());
  
  //Fading first
  for(int i=0;i<x1.size();++i)
    {
      vector<double> tVector;
      tVector.push_back(x1[i]);
      tVector.push_back(x2[i]);

      // y = H * s
      vector<complex<double> > rVector = vectorMultiplyMatrix(channel, tVector, "right");
      y1[i] = rVector[0].real();
      y2[i] = rVector[1].real();
    }
  

  //Add AWGN noise
  int seed = 1000;
  AWGNChannel(y1, sigma, seed);
  AWGNChannel(y2, sigma, seed+10000);

  x1 = y1;
  x2 = y2;
}


void MIMO22FastFadingChannel(vector<double> &x1, vector<double> &x2, vector<ComplexMatrix> &channel, double sigma)
{
   // Check if two antennas transmit the same number of symbols.
  if(x1.size() != x2.size()) {
    cout << "There should be the same number of symbols on two antennas" << endl;
    exit(0);
  }

  //Y1 = h11*x1 + h12*x2 + n1, Y2 = h21*x1 + h22*x2 + n2
  vector<double> y1(x1.size());
  vector<double> y2(x1.size());
  
  //Fading first
  for(int i=0;i<x1.size();++i)
    {
      vector<double> tVector;
      tVector.push_back(x1[i]);
      tVector.push_back(x2[i]);

      // y = H * s
      vector<complex<double> > rVector = vectorMultiplyMatrix(channel[i], tVector, "right");
      y1[i] = rVector[0].real();
      y2[i] = rVector[1].real();
    }
  

  //Add AWGN noise
  int seed = 1000;
  AWGNChannel(y1, sigma, seed);
  AWGNChannel(y2, sigma, seed+10000);

  x1 = y1;
  x2 = y2;
}

void MIMO22FastFadingChannel(vector<complex<double> > &x1, vector<complex<double> > &x2, vector<ComplexMatrix> &channel, double sigma)
{
   // Check if two antennas transmit the same number of symbols.
  if(x1.size() != x2.size()) {
    cout << "There should be the same number of symbols on two antennas" << endl;
    exit(0);
  }

  //Y1 = h11*x1 + h12*x2 + n1, Y2 = h21*x1 + h22*x2 + n2
  vector<complex<double> > y1(x1.size());
  vector<complex<double> > y2(x1.size());
  
  //Fading first
  for(int i=0;i<x1.size();++i)
    {
      vector<complex<double> > tVector;
      tVector.push_back(x1[i]);
      tVector.push_back(x2[i]);

      // y = H * s
      vector<complex<double> > rVector = vectorMultiplyMatrix(channel[i], tVector, "right");
      y1[i] = rVector[0];
      y2[i] = rVector[1];
    }
  

  //Add AWGN noise
  int seed = 1000;
  complexAWGNChannel(y1, sigma, seed);
  complexAWGNChannel(y2, sigma, seed+10000);

  x1 = y1;
  x2 = y2;
}

// This is general MIMO channel compatible with any channel with equal number of transmitter and receiver.
void MIMOFastFadingChannel(vector<vector<complex<double> > > &x, vector<ComplexMatrix> & channel, double sigma)
{
  // Get the MIMO channel size
  int transNum = x.size();
  int transLen = x[0].size();

  // Received signals
  vector<vector<complex<double> > > y(transNum, vector<complex<double> >(transLen));

  // Fading first
  for(int i=0;i<transLen;++i) {
    vector<complex<double> > tVector;
    for(int transI=0;transI<transNum;++transI)
      tVector.push_back(x[transI][i]);

    // y = H * s
    vector<complex<double> > rVector = vectorMultiplyMatrix(channel[i], tVector, "right");
    
    for(int recI=0;recI<y.size();++recI)
      y[recI][i] = rVector[recI];
  }

  // Add Awgn noise
  int seed = 1000;
  for(int recI=0;recI<y.size();++recI) {
    complexAWGNChannel(y[recI], sigma, seed);
    seed += 10000;
  }

  for(int recI=0;recI<y.size();++recI)
    x[recI] = y[recI];
}
  

//Precoding for the transmitted vector X to VX, where H = U * Z * V^T
void SVDPreCoding(vector<complex<double> > &x1, vector<complex<double> > &x2, ComplexMatrix V)
{
  vector<complex<double> > codedx1(x1.size());
  vector<complex<double> > codedx2(x1.size());

  for(int i=0;i<x1.size();++i)
    {
      vector<complex<double> > sVector;
      sVector.push_back(x1[i]);
      sVector.push_back(x2[i]);
      
      // \hat{x} = V * X
      vector<complex<double> > cVector = vectorMultiplyMatrix(V, sVector, "right");
      codedx1[i] = cVector[0];
      codedx2[i] = cVector[1];
    }

  x1 = codedx1;
  x2 = codedx2;
}

//Process received signal Y in order to get soft information that can be feed to the decoder.
void SVDProcessY(vector<complex<double> > &y1, vector<complex<double> > &y2, ComplexMatrix U)
{
  vector<complex<double> > codedy1(y1.size());
  vector<complex<double> > codedy2(y2.size());
  
  ComplexMatrix UHermitian = U.Hermitian();
  for(int i=0;i<y1.size();++i)
    {
      vector<complex<double> > rVector;
      rVector.push_back(y1[i]);
      rVector.push_back(y2[i]);

      // \hat{Y} = U^H * Y
      vector<complex<double> > cVector = vectorMultiplyMatrix(UHermitian, rVector, "right");
      codedy1[i] = cVector[0];
      codedy2[i] = cVector[1];
    }
  
  y1 = codedy1;
  y2 = codedy2;

}



/* The set of defining pdf vector functions */

// Compute the magnitude of a complex vector
double complexVecMagnitude(vector<complex<double> > &vec)
{
  vector<double> vecEleMag(vec.size());

  for(int i=0;i<vec.size();++i)
    vecEleMag[i] = norm(vec[i]);

  double sum = 0;
  for(int i=0;i<vec.size();++i)
    sum += vecEleMag[i];

  return sqrt(sum);
}

// Compute <v1, v2>
double innerProduct(vector<double> &v1, vector<double> &v2)
{
  if(v1.size() != v2.size()) {
    cout << "Sizes of two vectors do not match. Cannot compute inner product!" << endl;
  }
  
  double sum = 0;
  for(int i=0;i<v1.size();++i)
    {
      double temp = v1[i] * v2[i];
      sum += temp;
    }

  return sum;
}

  

// Normalize pdf vector to make the sum of vector as 1.
int vecNorm(vector<double> &origin, int sizeDesire)
{
  if(origin.size()%2 == 0 || (sizeDesire && sizeDesire%2 == 0)) {
    cout<<"the size of pdf vector should be odd!"<<endl;
    return -1;
  }

  if(sizeDesire == 0)
    sizeDesire = origin.size();
  
  vector<double> vfinal(sizeDesire,0);
  double sum=0;

  int n=origin.size();
  int middleIndex=(n-1)/2;
  int beginIndex=middleIndex-(sizeDesire-1)/2;
  int endIndex=middleIndex+(sizeDesire-1)/2;

  for(int i=0;i<beginIndex+1;++i)
    vfinal[0]+=origin[i];
  for(int i=endIndex;i<n;++i)
    vfinal[sizeDesire-1]+=origin[i];
  
  int index=1;
  for(int i=beginIndex+1;i<endIndex;++i)
    {
      vfinal[index]=origin[i];
      ++index;
    }

  for(int i=0;i!=sizeDesire;++i)
    sum+=vfinal[i];
  for(int i=0;i!=sizeDesire;++i)
    vfinal[i]/=sum;

  origin = vfinal;
};  


/* The part of operation for the pdf */
// This fuction takes two vector and takes operation of convolution
int conv(vector<double> &vect1, vector<double> &vect2, vector<double> &res)
{

  if(vect1.size()%2 == 0 || vect2.size()%2 == 0) {
    cout<<"the size of vectors should be odd! Fail."<<endl;
    return -1;
  }
  
  int size1 = vect1.size();
  int size2 = vect2.size();
  int ressize = size1+size2-1;
  vector<double> result(ressize, 0);
 
  for(int l=0;l<ressize;++l)
    {
      for(int n=0;n<size2;++n)
	{
	  if(vect2[n] != 0) {
	    if(l-n>=0 && l-n<size1)
	      {
		double n1=vect2[n];
		double n2=vect1[l-n];
		if(n1 && n2)
		  result[l]+=n1*n2;
	      }
	  }
	}
    }
  res = result;
}

double recurCalPro(vector<double> &sumpdf, vector<double> &vpdf, double weight, double n)
{

  double pro = 0;

  if(n > sumpdf.size()-1 || n < 0)
    return pro; 


  if(vpdf[1]>vpdf[2]) {

    if( (n-weight < 0) || (n-weight > sumpdf.size()-1) )
      pro = sumpdf[n] / vpdf[1];
    else
      pro = ( sumpdf[n] - recurCalPro(sumpdf, vpdf, weight, n-weight) * vpdf[2]) / vpdf[1] ;
  }
  else {

    if( (n+weight > sumpdf.size()-1) || (n+weight<0) )
      pro = 0;
    else
      pro = ( sumpdf[n+weight] - recurCalPro(sumpdf, vpdf, weight, n+weight) * vpdf[1]) / vpdf[2];
  } 

  if(pro < 0)
    pro = 0;

  return pro;
    
}

double deconvolve(vector<double> &sumpdf, vector<double> &vpdf, double weight, vector<double> &respdf)
{
  respdf.resize(sumpdf.size());
  
  int w = weight;
   // if p(0) > p(1), update Pc/v from left to the right. Conversely, go the opposite
  if(vpdf[1] > vpdf[2]) {
    if(w > 0) {
      for(int i = 0; i < sumpdf.size(); i++) {
				if(i - w < 0)
	  			respdf[i] = sumpdf[i] / vpdf[1];
				else
	  			respdf[i] = (sumpdf[i] - respdf[i - w] * vpdf[2]) / vpdf[1];
      } 
    } else {
      for(int i = sumpdf.size() - 1; i >= 0; i--) {
				if(i - w >= sumpdf.size())
	  			respdf[i] = sumpdf[i]  / vpdf[1];
				else
	  			respdf[i] = (sumpdf[i] - respdf[i - w] * vpdf[2]) / vpdf[1];
      }
    }
  } else {
    if(w > 0) {
      for(int i = sumpdf.size() - 1; i >= 0; i--) {
				if(i + w >= sumpdf.size())
	  			respdf[i] = 0;
				else
	  			respdf[i] = (sumpdf[i + w] - respdf[i + w] * vpdf[1]) / vpdf[2];
      }
    } else {
      for(int i = 0; i < sumpdf.size(); i++) {
				if(i + w < 0)
	  			respdf[i] = 0;
				else
	  			respdf[i] = (sumpdf[i + w] - respdf[i + w] * vpdf[1]) / vpdf[2];
      }
    }
  }
  double l = pow(10, -9);
  for(int i=0; i<respdf.size(); i++) {
    if(respdf[i] < l)
      respdf[i] = 0;
  
  }
  
  vecNorm(respdf);
}

double deconvolve2(vector<double> &sumpdf, vector<double> &vpdf, double weight, vector<double> &respdf) {
 
 	vector<double> tmpRes = vector<double>(sumpdf.size());
  
  int w = weight;
	
	int zeroPoint = 2;

	double threshold = pow(10, -9);
	// Depends on w is positive or negative, the update direction is quite different. If we look at the 
	// equation carefully: P(Y = a) = P(X = a) * P(x1 = 0) + P(X = a - w1) * P(x1 = 1) + P(X = a - 2w1) * P(x1 = 2).
  
	// P(r = a) = (P(Y = a) - P(r = a - w1) * P(x1 = 1) - 
	//            P(r = a - 2 * w1) * P(x1 = 2)) / P(x1 = 0)
	if(vpdf[2] > vpdf[4]) {
		if(w > 0) {

			// Actually the value of i is the value of a.
			for(int i = 0; i < sumpdf.size(); i++) {
				double nominator = sumpdf[i];

				if(i - w >= 0)
					nominator -= tmpRes[i - w] * vpdf[3];
				if(i - 2 * w >= 0)
					nominator -= tmpRes[i - 2 * w] * vpdf[4];
				
				tmpRes[i] = nominator / vpdf[2];

				if(abs(tmpRes[i]) < threshold || tmpRes[i] < 0)
					tmpRes[i] = 0;

			} 
  	} else {
				for(int i = sumpdf.size() - 1; i >= 0; i--) {
					double nominator = sumpdf[i];		 		
					
					if(i - w < sumpdf.size())
						nominator -= tmpRes[i - w] * vpdf[3];
					if(i - 2 * w < sumpdf.size())
						nominator -= tmpRes[i - 2 * w] * vpdf[4];
				
					tmpRes[i] = nominator / vpdf[2];
					if(abs(tmpRes[i]) < threshold || tmpRes[i] < 0)
						tmpRes[i] = 0;
  	   	}
  		}
	} 
	// P(r = a) = (P(Y = a + 2 * w1) - P(r = a + 2 * w1) * P(x1 = 0) - 
	//            P(r = a + w1) * P(x1 = 1)) / P(x1 = 2)
	else if(vpdf[4] >= vpdf[2]) {
		if(w > 0) {
			for(int i=sumpdf.size()-1; i>=0; i--) {
				
				if(i + 2 * w >= sumpdf.size())  continue;

				double nominator = sumpdf[i + 2 * w];
				nominator -= tmpRes[i + 2 * w] * vpdf[2];
				nominator -= tmpRes[i + w] * vpdf[3];
				
				tmpRes[i] = nominator / vpdf[4];
				if(abs(tmpRes[i]) < threshold || tmpRes[i] < 0)
					tmpRes[i] = 0;
					
			}
		} else {
			for(int i=0; i<sumpdf.size(); i++) {
				
				if( i + 2 * w < 0)  continue;

				double nominator = sumpdf[i + 2 * w];
				nominator -= tmpRes[i + 2 * w] * vpdf[2];
				nominator -= tmpRes[i + w] * vpdf[3];

				tmpRes[i] = nominator / vpdf[4];
				if(abs(tmpRes[i]) < threshold || tmpRes[i] < 0)
					tmpRes[i] = 0;
			}
		}
	} 

	respdf = tmpRes;
	for(int i=0; i<respdf.size(); i++) {
		if(std::isnan(abs(respdf[i])) || std::isinf(respdf[i]))
			cout << "nan message!" << endl;
	}

  vecNorm(respdf);

	double sum = 0;
 	for(int i=0; i<respdf.size(); i++) {
		sum += respdf[i];
	}

	if(std::isnan(sum) || sum < 0.5)
		cout << "problem!" << endl;
}

int weightPdf(double weight, vector<double> &pdf)
{
  if(pdf.size()%2 == 0) {
    cout<<"the size of pdf vector should be odd!"<<endl;
    return -1;
  }

  int newsize = (pdf.size()-1) * abs(weight)+1;
  vector<double> newpdf(newsize);

  int midpoint_old = (pdf.size()-1)/2;        // 0 position in the old pdf
  int midpoint_new = (newsize-1)/2;           // 0 position in the weighted pdf

  for(int index=0;index<pdf.size();++index)	{
      double realvalue = index - midpoint_old; //the value of random variables
      double weightvalue = weight*realvalue;   //weighted value of random variables

      int newpos = weightvalue + midpoint_new;
        
			// In the case that weight is not an integer, only calculate the probability for meaningful
			// value, which are all integers.
      if( int(newpos) != newpos )            
				continue;                             
      else                                    
				newpdf[newpos] = pdf[index];
	}

  pdf = newpdf;
}

// Project vector a onto vector b
vector<double> vecProjection(vector<double> &a, vector<double> &b)
{
  double scalar = innerProduct(a, b) / innerProduct(b, b);
  vector<double> c(b);

  for(int i=0;i<c.size();++i)
    c[i] *= scalar;

  return c;
}
  

// Now it is for any square matrix
vector<vector<double> > GramSchmit(vector<vector<double> > &matrix)
{
  int maSize = matrix.size();
  
  // Start orthogonization
  
  vector<vector<double> > baseVectors;
  for(int i=0;i<maSize;++i)
    {
      // Get column vector ai
      vector<double> ai(matrix[i].size());
      for(int rowi=0;rowi<ai.size();++rowi)
	ai[rowi] = matrix[rowi][i];

      // Get the summation of projection of ai onto 
      // all the previous base vectors.
      vector<double> sumProj(ai.size(), 0);
      for(int basei=0;basei<baseVectors.size();++basei)
	{
	  vector<double> proj = vecProjection(ai, baseVectors[basei]);
	  
	  for(int veci=0;veci<sumProj.size();++veci)
	    sumProj[veci] += proj[veci];
	}

      // Get vector ui = ai - sum(projection)
      vector<double> ui(ai);
      for(int veci=0;veci<ai.size();++veci)
	ui[veci] = ai[veci] - sumProj[veci];

      // Get base vector ei
      double sum = 0;
      for(int veci=0;veci<ui.size();++veci)
	sum += pow(ui[veci], 2);
      sum = sqrt(sum);
      for(int veci=0;veci<ui.size();++veci)
				ui[veci] /= sum;

      baseVectors.push_back(ui);
    }

  return baseVectors;

}

// Decoupling of 2X2 MIMO channel
vector<vector<double> > MIMODecoupling(ComplexMatrix &ma)
{
  vector<vector<double> > decouplingMa(4, vector<double>(4));
  vector<vector<double> > real(2, vector<double>(2)), imaginery(2, vector<double>(2));
  
  for(int i=0;i<2;++i)
    {
      for(int j=0;j<2;++j)
	{
	  real[i][j] = ma.matrix[i][j].real();
	  imaginery[i][j] = ma.matrix[i][j].imag();
	}
    }

  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      decouplingMa[i][j] = real[i][j];

  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      decouplingMa[i][j+2] = -imaginery[i][j];

  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      decouplingMa[i+2][j] = imaginery[i][j];

  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      decouplingMa[i+2][j+2] = real[i][j];

  return decouplingMa;
}
  
vector<vector<vector<double> > >  SDMatrixDecomp(ComplexMatrix &ma)
{
  
  // Decouple the complex channel matrix into double size real matrix
  vector<vector<double> > SDMatrix = MIMODecoupling(ma);

  // Get base vectors of the decoupling matrix.
  vector<vector<double> > baseVector = GramSchmit(SDMatrix);

  // Decompose the matrix into upper triangular matrix and an othornomal matrix
  vector<vector<double> > QMatrix(baseVector[0].size(), vector<double>(baseVector.size()));
  for(int i=0;i<baseVector.size();++i)
    for(int j=0;j<baseVector[i].size();++j)
      QMatrix[j][i] = baseVector[i][j];

  vector<vector<double> > RMatrix(baseVector.size(), vector<double>(baseVector.size()));
  for(int j=0;j<RMatrix[0].size();++j)
    {
      // Get column vector ai
      vector<double> ai(SDMatrix[j].size());
      for(int rowi=0;rowi<ai.size();++rowi)
	ai[rowi] = SDMatrix[rowi][j];
      
      for(int i=0;i<SDMatrix.size();++i)
	{
	  if(j >= i)
	    RMatrix[i][j] = innerProduct(baseVector[i], ai);
	  else
	    RMatrix[i][j] = 0;
	}
    }

  vector<vector<vector<double> > > maDecomp;
  maDecomp.push_back(QMatrix);
  maDecomp.push_back(RMatrix);

  return maDecomp;

}


vector<vector<double> > LSDecoder(vector<complex<double> > &y, vector<double> &symbolset, ComplexMatrix HMatrix, double norfactor, double sigma)
{

  // Compute H' = alpha * H, where alpha is the normalization factor
  ComplexMatrix proH(HMatrix.size);
  for(int i=0;i<HMatrix.size;++i)
    for(int j=0;j<HMatrix.matrix[i].size();++j)
      proH.setEntryVal(i,j,complex<double>(norfactor,0) * HMatrix.getEntryVal(i,j)); 

  // Decompose MIMO channel matrix
  vector<vector<vector<double> > > realH = SDMatrixDecomp(proH);
  vector<vector<double> > QMatrix = realH[0];
  vector<vector<double> > RMatrix = realH[1];
  
  // Get the decoupled received signals
  vector<double> realY(y.size()*2);
  for(int i=0;i<y.size();++i)
    realY[i] = y[i].real();
  for(int i=0;i<y.size();++i)
    realY[i+y.size()] = y[i].imag();

  // Get y' = Q^T
  vector<double> proY(realY.size(), 0);
  for(int index = 0; index < proY.size(); ++index)
    for(int j = 0; j < QMatrix.size(); ++j)
      proY[index] += QMatrix[j][index] * realY[j];  
  
   // Candidate list
  vector<vector<double> > candidates;
  vector<double> oneCandidate;

  // Determine the radius
  double radius = 4;
  // According to the symbol set size, determine the necessary candidate number
  int needCandiNum;
  if(symbolset.size() == 2)
    needCandiNum = 16;
  else
    needCandiNum = 5000;
  int dimension = 4;



  // Depth-first search
  int index = 3;
  findSphereCandidate(candidates, radius, proY, symbolset, QMatrix, RMatrix, index, oneCandidate, needCandiNum);

  // If candidate number is too small, adapt the radius and choose candidates again.
  while(candidates.size() < needCandiNum) {

    radius += 1;

    // Clear candidates last time
    candidates.clear();
    
    findSphereCandidate(candidates, radius, proY, symbolset, QMatrix, RMatrix, index, oneCandidate, needCandiNum);
  }

  // Original candidates data also contains radius info(radius squared), which would be used in soft estimation
  // Another important thing is that the order of the candidate now is [Re(x), Im(x)], which is not the same
  // as transmitted symbols ([X1+iX2, X3+iX4]). Therefore, the order should be adjusted.
  for(int i=0;i<candidates.size();++i)
    swap(candidates[i][1], candidates[i][2]);
    

  vector<vector<double> > pdfMessage = sphereSoftEstimator(y, symbolset, candidates, norfactor, HMatrix, sigma);

  return pdfMessage;
}

void findSphereCandidate(vector<vector<double> > &candidates, double d, vector<double> &Y, vector<double> &symbolset, vector<vector<double> > &QMatrix, vector<vector<double> > &RMatrix, const double symIndex, vector<double> curCandidate, int candiNum)
{
  int dimension = 4;

  // The boundary for each S is as follows:
  //  (-d(i) + y(i/di)) / R[i][i] <= S(i) <= ((d(i) + y(i/di)) / R[i][i])
  // Therefore, before calculating boundary for the current level S(i),
  // y(i/di) and d for current depth should be obtained.


  // Calculate y(i/di) = y(i) - sum{j=i+1}{di} R[i][j]*S[j]
  double y_i = Y[symIndex];
  vector<double> S(curCandidate);
  for(int j=dimension-1;j >= symIndex+1; --j)
    {
      double upperS = S[S.size()-1];
      y_i -= RMatrix[symIndex][j] * upperS;
      S.erase(S.end()-1);
    }

  
  // Calculate upper bound and lower bound for s_i
  double upperbound = (d + y_i) / RMatrix[symIndex][symIndex];
  double lowerbound = (-d + y_i) / RMatrix[symIndex][symIndex];

  upperbound = floor(upperbound);
  lowerbound = ceil(lowerbound);

  // Check the index of upper and lower bound in the symbol set
  // Notice that symbol set is ordered.
  int upperindex = 0;
  for(int i=symbolset.size()-1;i>=0;--i)
    {
      if(symbolset[i] > upperbound)
	continue;
      else {
	upperindex = i;
	break;
      }
    }
  int lowerindex = 0;
  for(int i=0;i<symbolset.size();++i)
    {
      if(symbolset[i] < lowerbound)
	continue;
      else {
	lowerindex = i;
	break;
      }
    }

  // If there is no value for Sm to choose, just return to S_(m+1)
  if(lowerindex > upperindex) 
    return;
  
  // Choose a value for current s, which is lowerbound <= s <= upperbound
  for(int index = lowerindex; index <= upperindex; ++index) 
    {
      double s = symbolset[index];
     
      // If the process comes to s1, one candidate has been found
      if(symIndex == 0) {

	curCandidate.insert(curCandidate.begin(), s);

	// Calculate the radius of current candidate
	double curR = 0;
	for(int rowi = 0; rowi < dimension; ++rowi)
	  {
	    double curY = Y[rowi];
	    for(int coli = rowi; coli < dimension; ++coli)
	      curY -= RMatrix[rowi][coli] * curCandidate[coli];
	    curR += pow(curY, 2);
	  }
      
	// Push the radius squared of current vector into the end of candidate
	curCandidate.push_back(curR);

	
	// Find the corresponding index in the candidates list, and update the list
	int insertIndex = 0;
	while(insertIndex < candidates.size()) {
	    
	  if(curCandidate[dimension] > candidates[insertIndex][dimension])  
	    ++insertIndex;
	  else
	    break;
	}
	candidates.insert(candidates.begin()+insertIndex, curCandidate);
	if(candidates.size() > candiNum)
	  candidates.resize(candiNum);
	
	// Prepare for the next candidate on current level
	curCandidate.erase(curCandidate.begin());
	curCandidate.erase(curCandidate.end()-1);
	continue;
      }
      

      // If the process comes to other s, continue the process.
      curCandidate.insert(curCandidate.begin(), s);

      // Compute new d for the next candidate
      // d(i-1)^2 = d(i)^2 - (y(i) - sum{j=i}{di} R{i][j] * S[j])^2
      double newd = pow(d, 2);
      double extraComp = Y[symIndex];
      for(int j=symIndex;j<dimension;++j)
	extraComp -= RMatrix[symIndex][j] * curCandidate[j-symIndex];

      newd = sqrt(newd - pow(extraComp, 2));

      // Compute for the next candidate
      findSphereCandidate(candidates, newd, Y, symbolset, QMatrix, RMatrix, symIndex-1, curCandidate, candiNum);

      // Since all candidates based on chosen S(i) at this level have been found, current S(i) 
      // would be eliminated and be ready to choose next available S(i) at this level
      // while keep symbols of the upper level.
      curCandidate.erase(curCandidate.begin());
    }
}
  
  
vector<vector<double> >  sphereSoftEstimator(vector<complex<double> > y, vector<double> &symbolset, vector<vector<double> > &candidates, double norfactor, ComplexMatrix &channel, double sigma)
{
  // The number of symbols to be estimated for each received vector based on the modulation method.
  int symbolnum = 4;
  
  // Get the number of possible values for each symbol
  // and construct a traverse table
  int possiblenum = symbolset.size();
  
  // Compute the pdf for each possible x based on received vector y = H * s + n.
  vector<double> originalSymbols(symbolnum, 0);
  map<vector<double>, double> prodata;
  
  double checksum = 0;

  // Calculate the probability of each candidate vector
  for(int vecindex=0;vecindex<candidates.size();++vecindex) {
	
		double yhsSquare = candidates[vecindex][symbolnum];
		double pro = exp( (-1) * yhsSquare  / (2 * pow(sigma,2)) ) / pow(2*pi*pow(sigma,2), 2);
	  
		// Store the original symbol vector and corresponding pdf to the prodata
		candidates[vecindex].resize(symbolnum);
		originalSymbols = candidates[vecindex];
		prodata[originalSymbols] = pro;
	}

  //Determine the size of the pdf vector for each symbol in the vector
  int pdfsize;
  if(possiblenum == 2)  // Digital bit
    pdfsize = 3;
  else
    pdfsize = 2 * symbolset[symbolset.size()-1] + 1;   // Pdf size depends on the largest value in the set.
  int zeropoint = (pdfsize - 1) / 2;

  vector<double> pdfinitial(pdfsize, 0);
  vector<vector<double> > pdfdata(symbolnum, pdfinitial);


  // Iterate through all the possible sent vectors
  for(int vecindex=0;vecindex<candidates.size();++vecindex)
    {
      // Construct the original vector for current case
      vector<double> curvec(candidates[vecindex]);
      
      double pro = prodata[curvec];

      // Process each each position in the symbol vector
      for(int symi=0;symi<symbolnum;++symi)
	{
	  int candivalue = curvec[symi];
	  int posOfCandi = zeropoint + candivalue;
	  pdfdata[symi][posOfCandi] += pro;
	}
    }

  // Normalization
  for(int symi=0;symi<symbolnum;++symi)
    vecNorm(pdfdata[symi]);

  return pdfdata;
}

// Synthetic RP decoder with individual sys nodes.
vector<vector<double> > syntheticDecoderRP(Coded_info &symNodeInfo1, Coded_info &symNodeInfo2, int symIndex, vector<Sys_info>& sysNodeInfo1, vector<Sys_info>& sysNodeInfo2, vector<double>& channelPdf, double norFactor) {

	// Get LLR messages from systematic node group and transform messages
	// into synthetic message: P0, P1, P2: the probability of the sum
	// of two corresponding systematic bits being i. And the information
	// is represented in the form of pdf vector.
	
	int degree1 = symNodeInfo1.nodeData.degree;
	vector<double> inmessage1(degree1);

	for(int index=0; index<degree1; ++index) {
  	int neiNum = symNodeInfo1.nodeData.neighbourNum[index];
    int messageIndex = symNodeInfo1.nodeData.inmessageIndex[index];
    inmessage1[index] = sysNodeInfo1[neiNum].getMessage(messageIndex, 1);
	}

	int degree2 = symNodeInfo2.nodeData.degree;
	vector<double> inmessage2(degree2);

	for(int index=0; index<degree2; ++index) {
  	int neiNum = symNodeInfo2.nodeData.neighbourNum[index];
    int messageIndex = symNodeInfo2.nodeData.inmessageIndex[index];
    inmessage2[index] = sysNodeInfo2[neiNum].getMessage(messageIndex, 1);
	}

	// Make sure that the degrees of 2 RP symbol nodes are the same.
	// Since RCM matrix are the same. They also have the same weights.
	if(degree1 != degree2) {
		cerr << "The degree of two RP symbol nodes are different. Systhetic decoding is not able to perform!" << endl;
		exit(EXIT_FAILURE) ;
	}
	
	int degree = degree1;

	// Out going pdfs.
	vector<vector<double> > outPdf;

	// There are total degree pdfs from sys node and the length of
	// each one is 5.
	vector<vector<double> >sysPdf(degree, vector<double>(5));
	for(int index=0; index<degree; index++) {
		double mess1 = inmessage1[index];
		double mess2 = inmessage2[index];

		double p0_1 = exp(mess1) / (1 + exp(mess1));
		double p1_1 = 1 - p0_1;
		double p0_2 = exp(mess2) / (1 + exp(mess2));
		double p1_2 = 1 - p0_2;

		sysPdf[index][2] = p0_1 * p0_2;
		sysPdf[index][3] = p0_1 * p1_2 + p1_1 * p0_2;
		sysPdf[index][4] = p1_1 * p1_2;

	}

	// Get normalized weight set.
	vector<double> weightSet(symNodeInfo1.weightset);
	vector<double> norWeight(degree);
	for(int i=0; i<degree; i++)
		norWeight[i] = norFactor * weightSet[i]; 
	
	// Get the pdf of the value after synthetic values multiplied by
	// the weight.
	vector<vector<double> > sysAlias(sysPdf);
	for(int i=0; i<degree; i++)
		weightPdf(weightSet[i], sysAlias[i]);

	// Get the pdf of linear combination.
	vector<double> sumPdf(1, 1);
	for(int i=0; i<degree; i++)
		conv(sumPdf, sysAlias[i], sumPdf);

	vecNorm(sumPdf, channelPdf.size());

	for(int i=0; i<degree; i++) {
		
		vector<double> resPdf(1, 1);
		deconvolve2(sumPdf, sysPdf[i], weightSet[i], resPdf);

		// Compute p0, p1, p2, which is the probability of systhetic node
		// beging x. Then we can know the probability of 00,
		// 01, 10, 11 and we can calculate LLR message to each
		// connected systematic node.
		double zeroPoint = (resPdf.size() - 1) / 2;

		double p0 = 0, p1 = 0, p2 = 0;
		for(int symI = 0; symI < channelPdf.size(); symI++) {

				p0 += resPdf[symI] * channelPdf[symI];
				
				int comPos = symI - weightSet[i];
				if(comPos >= 0 && comPos < resPdf.size())
	  			p1 += resPdf[comPos] * channelPdf[symI];

				comPos = symI - 2 * weightSet[i];
				if(comPos >= 0 && comPos < resPdf.size())
					p2 += resPdf[comPos] * channelPdf[symI];
	
		}

		double p0Individual = p0 + p1 / 2;
		double p1Individual = p2 + p1 / 2;
		double message = 0;
		double threshold = 50;

		if(p0Individual == 0 || p1Individual == 0)
			message = 0;
    else {
			message = log(p0Individual / p1Individual);
			if(message > threshold)
	  		message = threshold;
			if(message < -threshold)
	  		message = -threshold;
			if(std::isnan(message))
	  		std::cout<<"invalid message!"<<std::endl;
      }

			symNodeInfo1.setMessage(message, i);
			symNodeInfo2.setMessage(message, i);
	
			// Out going pdf one the current edge.
			vector<double> oneOutMessage {0, 0, p0, p1, p2};
			vecNorm(oneOutMessage);
			outPdf.emplace_back(oneOutMessage);
	}
	return outPdf;
}

// Synthetic RP decoder with synthetic sys nodes.
vector<vector<double> > syntheticDecoderRP(Coded_info &symNodeInfo1, Coded_info &symNodeInfo2, int symIndex, vector<vector<vector<double> > > &pdfMessageFromSys, vector<double>& channelPdf, double norFactor) {

	int degree = symNodeInfo1.nodeData.degree;
	vector<vector<double> > sysPdf(degree, vector<double>(5, 0));

	for(int index=0; index<degree; ++index) {
  	int neiNum = symNodeInfo1.nodeData.neighbourNum[index];
    int messageIndex = symNodeInfo1.nodeData.inmessageIndex[index];
    sysPdf[index] = pdfMessageFromSys[neiNum][messageIndex];
	}
	
	// Out going pdfs.
	vector<vector<double> > outPdf;

	// Get normalized weight set.
	vector<double> weightSet(symNodeInfo1.weightset);
	vector<double> norWeight(degree);
	for(int i=0; i<degree; i++)
		norWeight[i] = norFactor * weightSet[i]; 
	
	// Get the pdf of the value after synthetic values multiplied by
	// the weight.
	vector<vector<double> > sysAlias(sysPdf);
	for(int i=0; i<degree; i++)
		weightPdf(weightSet[i], sysAlias[i]);

	// Get the pdf of linear combination.
	vector<double> sumPdf(1, 1);
	for(int i=0; i<degree; i++)
		conv(sumPdf, sysAlias[i], sumPdf);

	vecNorm(sumPdf, channelPdf.size());

	for(int i=0; i<degree; i++) {
		
		vector<double> resPdf(1, 1);
		deconvolve2(sumPdf, sysPdf[i], weightSet[i], resPdf);

		// Compute p0, p1, p2, which is the probability of systhetic node
		// beging x. Then we can know the probability of 00,
		// 01, 10, 11 and we can calculate LLR message to each
		// connected systematic node.
		double zeroPoint = (resPdf.size() - 1) / 2;

		double p0 = 0, p1 = 0, p2 = 0;
		for(int symI = 0; symI < channelPdf.size(); symI++) {

				p0 += resPdf[symI] * channelPdf[symI];
				
				int comPos = symI - weightSet[i];
				if(comPos >= 0 && comPos < resPdf.size())
	  			p1 += resPdf[comPos] * channelPdf[symI];

				comPos = symI - 2 * weightSet[i];
				if(comPos >= 0 && comPos < resPdf.size())
					p2 += resPdf[comPos] * channelPdf[symI];
	
		}

		// Out going pdf one the current edge.
		vector<double> oneOutMessage {0, 0, p0, p1, p2};
		vecNorm(oneOutMessage);
		outPdf.emplace_back(oneOutMessage);
	}
	return outPdf;
}

// For hybrid. Take pdf message from RP nodes.
vector<vector<double> > syntheticDecoderSys(Sys_info &sysNodeInfo1, Sys_info &sysNodeInfo2, vector<vector<vector<double> > >& pdfMessageFromRP, vector<Coded_info> &neiNode1, vector<Coded_info>& neiNode2, vector<double>& channelInfo, vector<double>& channelLLR, vector<double>& sideInfo, double y, double sigma, double p, int index, ofstream &proTrack, int sys1, int sys2, int iteration) {

	int degree = sysNodeInfo1.nodeData.degree;
	int paraDegree = sysNodeInfo1.nodeData.para_degree;

	vector<vector<double> > inPdfMessage(degree);
	vector<double> inLLR1(paraDegree);
	vector<double> inLLR2(paraDegree);

	for(int i=0; i<degree; i++) {
		int neiNum = sysNodeInfo1.nodeData.neighbourNum[i];
		int messageIndex = sysNodeInfo1.nodeData.inmessageIndex[i];
		inPdfMessage[i] = pdfMessageFromRP[neiNum][messageIndex];
	}

	for(int i=0; i<paraDegree; i++) {
		int neiNum1 = sysNodeInfo1.nodeData.para_neighbourNum[i];
		int messageIndex1 = sysNodeInfo1.nodeData.para_inmessageIndex[i];
		inLLR1[i] = neiNode1[neiNum1].getMessage(messageIndex1);

		int neiNum2 = sysNodeInfo2.nodeData.para_neighbourNum[i];
		int messageIndex2 = sysNodeInfo2.nodeData.para_inmessageIndex[i];
		inLLR2[i] = neiNode2[neiNum2].getMessage(messageIndex2);
	}

	double p0FromRP = 1;
	double p1FromRP = 1;
	double p2FromRP = 1;

	for(int i=0; i<degree; i++) {
		p0FromRP *= inPdfMessage[i][2];
		p1FromRP *= inPdfMessage[i][3];
		p2FromRP *= inPdfMessage[i][4];
	}

	// Normalization probability mass from RP.
	double pTotalFromRP = p0FromRP + p1FromRP + p2FromRP;
	p0FromRP /= pTotalFromRP;
	p1FromRP /= pTotalFromRP;
	p2FromRP /= pTotalFromRP;

	// Combine messages from coded nodes
	double LLRFromLDGM1 = 0, LLRFromLDGM2 = 0;
	for(int i=0; i<paraDegree; i++) {
		LLRFromLDGM1 += inLLR1[i]; 
		LLRFromLDGM2 += inLLR2[i];
	}
	
	double p1FromCoded1 = 1 / (1 + exp(LLRFromLDGM1));
	double p0FromCoded1 = 1 - p1FromCoded1;
	double p1FromCoded2 = 1 / (1 + exp(LLRFromLDGM2));
	double p0FromCoded2 = 1 - p1FromCoded2;
	
	double p0FromCoded = p0FromCoded1 * p0FromCoded2;
	double p1FromCoded = p0FromCoded1 * p1FromCoded2 + p1FromCoded1 * p0FromCoded2;
	double p2FromCoded = p1FromCoded1 * p1FromCoded2;

	// Combine everything and normalize them.
	double p0Final = p0FromCoded * p0FromRP * channelInfo[2] * (1 - p) * 0.5;
	double p1Final = p1FromCoded * p1FromRP * channelInfo[3] * p;
	double p2Final = p2FromCoded * p2FromRP * channelInfo[4] * (1 - p) * 0.5;
	double pTotalFinal = p0Final + p1Final + p2Final;
	p0Final /= pTotalFinal;
	p1Final /= pTotalFinal;
	p2Final /= pTotalFinal;
	
	// Compute the message to the synthetic RP nodes
	vector<vector<double> > outPdf(degree, vector<double>(5, 0));
	for(int i=0; i<degree; i++) {
		outPdf[i][2] = p0Final / inPdfMessage[i][2];
		outPdf[i][3] = p1Final / inPdfMessage[i][3];
		outPdf[i][4] = p2Final / inPdfMessage[i][4];
		vecNorm(outPdf[i]);
	}

	// Method1: Compute the probability on systematic nodes first.
	// Compute the message to the LDGM nodes.
	double p1Prior = 1 / (1 + exp(sideInfo[index]));
	double p0Prior = 1 - p1Prior;
	
	double pLDGM1_0 = p0FromCoded1 * (p0Prior * (1 - p) + p1Prior * p) * (channelInfo[2] + channelInfo[3] * p1Prior) * (p0FromRP + p1FromRP * p1Prior);
	double pLDGM1_1 = p1FromCoded1 * (p1Prior * (1 - p) + p0Prior * p) * (channelInfo[4] + channelInfo[3] * p0Prior) * (p2FromRP + p1FromRP * p0Prior);

	
	/* // Update LDGM2 without using channel and RP information.
	double pLDGM2_0 = p0FromCoded2 * (pLDGM1_0 * (1 - p) + pLDGM1_1 * p);
	double pLDGM2_1 = p1FromCoded2 * (pLDGM1_1 * (1 - p) + pLDGM1_0 * p);
*/

	// Use message from RP and Channel with LDGM1 pro
	double pLDGM2_0 = p0FromCoded2 * (pLDGM1_0 * (1 - p) + pLDGM1_1 * p) * (channelInfo[2] + channelInfo[3] * pLDGM1_1) * (p0FromRP + p1FromRP * pLDGM1_1);
	double pLDGM2_1 = p1FromCoded2 * (pLDGM1_0 * p + pLDGM1_1 * (1 - p)) * (channelInfo[4] + channelInfo[3] * pLDGM1_0) * (p2FromRP + p1FromRP * pLDGM1_0);

	
	/* // Use Message from RP and Channel Unbiased
	double pLDGM2_0 = p0FromCoded2 * (p0FromCoded1 * (1 - p) + p1FromCoded1 * p) * (channelInfo[2] + 0.5 * channelInfo[3]) * (p0FromRP + 0.5 * p1FromRP);
	double pLDGM2_1 = p1FromCoded2 * (p1FromCoded1 * (1 - p) + p0FromCoded1 * p) * (channelInfo[4] + 0.5 * channelInfo[3]) * (p2FromRP + 0.5 * p1FromRP);
	*/

	sideInfo[index] = log(pLDGM2_0 / pLDGM2_1);

	double LLRFinal1 = log(pLDGM1_0 / pLDGM1_1);
	double LLRFinal2 = log(pLDGM2_0 / pLDGM2_1); 

	double threshold = 30;
	for(int i=0; i<paraDegree; i++) {
		double outMessage = LLRFinal1 - inLLR1[i];

		if(outMessage < -threshold)  outMessage = -threshold;
		else if(outMessage > threshold)  outMessage = threshold;

		sysNodeInfo1.setMessage(outMessage, i, 2);
				
		if(std::isnan(outMessage))
		cout << "Outgoing messages from sys nodes are nan!" << endl;
	
		outMessage = LLRFinal2 - inLLR2[i];

		if(outMessage < -threshold)  outMessage = -threshold;
		else if(outMessage > threshold)  outMessage = threshold;

		sysNodeInfo2.setMessage(outMessage, i, 2);
				
		if(std::isnan(outMessage))
		cout << "Outgoing messages from sys nodes are nan!" << endl;
	}
	/*
//	// Initial LLR for sys node 1
//	double LLRForSys = log((p0RPandChannel + 0.5 * p1RPandChannel) / 
//												  (p2RPandChannel + 0.5 * p1RPandChannel));
//	*/
//	
//	double LLRFromRP = log((p0FromRP + 0.5 * p1FromRP) / (p2FromRP + 0.5 * p1FromRP));
//
//	/*	
//	if(sys1 == 0)
//		channelLLR[index] = 0;
//	else
//		channelLLR[index] = 0;
//	*/
//	double LLRFinal1 = channelLLR[index] + LLRFromRP + LLRFromLDGM1 + sideInfo[index];
//	
//	/*
//	// If P1 from RP is greater than the other two, this information
//	// would not be used.
//	if(p1FromRP >= p0FromRP && p1FromRP >= p2FromRP)
//		LLRFinal1 -= LLRFromRP;
//	if(channelInfo[3] >= channelInfo[2] && channelInfo[3] >= channelInfo[4])
//		LLRFinal1 -= channelLLR[index];
//	*/
//
//	/*
//	double LLRFinal2 = channelLLR[index] + LLRFromRP + LLRFromLDGM2;
//	if(p1FromRP >= p0FromRP && p1FromRP >= p2FromRP)
//		LLRFinal2 -= LLRFromRP;
//	if(channelInfo[3] >= channelInfo[2] && channelInfo[3] >= channelInfo[4])
//		LLRFinal2 -= channelLLR[index];
//*/
//		
//	
//	// Incorporate correlation into LLR.
//	//LLRFinal1 += log(((1 - p) * exp(LLRFinal2) + p) / (1 - p + p * exp(LLRFinal2)));
//	//LLRFinal2 += log(((1 - p) * exp(LLRFinal1) + p) / (1 - p + p * exp(LLRFinal1)));
//
//	double threshold = 30;
//	for(int i=0; i<paraDegree; i++) {
//		double outMessage = LLRFinal1 - inLLR1[i];
//
//		if(outMessage < -threshold)  outMessage = -threshold;
//		else if(outMessage > threshold)  outMessage = threshold;
//
//		sysNodeInfo1.setMessage(outMessage, i, 2);
//				
//		if(std::isnan(outMessage))
//		cout << "Outgoing messages from sys nodes are nan!" << endl;
//	
//		/*		
//		outMessage = LLRFinal2 - inLLR2[i];
//
//		if(outMessage < -threshold)  outMessage = -threshold;
//		else if(outMessage > threshold)  outMessage = threshold;
//
//		sysNodeInfo2.setMessage(outMessage, i, 2);
//				
//		if(std::isnan(outMessage))
//		cout << "Outgoing messages from sys nodes are nan!" << endl;
//		*/
//	}
//
//		// Decoding process for sys node 2.
//		
//		// Calculate channel state message to sys node 2
//		double Esender = 1;
//		double messageToState = LLRFinal1 - channelLLR[index];
//		double nominator1 = exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(messageToState);
//		double nominator2 = (exp(-pow(y+2/sqrt(2)*sqrt(Esender), 2) / (2 * pow(sigma, 2))));
//		double denominator1 = exp(-pow(y-2/sqrt(2)*sqrt(Esender), 2) / (2 * pow(sigma, 2))) / exp(messageToState);
//		double denominator2 = (exp(-pow(y, 2) / (2 * pow(sigma, 2))));
//
//		if(std::isinf(nominator1) && std::isinf(denominator1)) {
//			channelLLR[index] = log(exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * pow(sigma, 2))));
//		} else {
//				channelLLR[index] = log((nominator1 + nominator2) / (denominator1 + denominator2));
//		}
//
//		// Update the correlation message for sys node 2.
//		sideInfo[index] = log(((1 - p) * exp(LLRFinal1) + p) / (1 - p + p * exp(LLRFinal1)));
//
//		double LLRFinal2 = channelLLR[index] + LLRFromLDGM2 + sideInfo[index];
//
//		for(int i=0; i<paraDegree; i++) {
//			double outMessage = LLRFinal2 - inLLR2[i];
//
//			if(outMessage < -threshold)  outMessage = -threshold;
//			else if(outMessage > threshold)  outMessage = threshold;
//
//			sysNodeInfo2.setMessage(outMessage, i, 2);
//				
//			if(std::isnan(outMessage))
//			cout << "Outgoing messages from sys nodes are nan!" << endl;
//		}
//
//		// Compute the state message to the sys node 1.
//		messageToState = LLRFinal2 - channelLLR[index];
//		nominator1 = exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(messageToState);
//		nominator2 = (exp(-pow(y+2/sqrt(2)*sqrt(Esender), 2) / (2 * pow(sigma, 2))));
//		denominator1 = exp(-pow(y-2/sqrt(2)*sqrt(Esender), 2) / (2 * pow(sigma, 2))) / exp(messageToState);
//		denominator2 = (exp(-pow(y, 2) / (2 * pow(sigma, 2))));
//
//		if(std::isinf(nominator1) && std::isinf(denominator1)) {
//			channelLLR[index] = log(exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * pow(sigma, 2))));
//		} else {
//				channelLLR[index] = log((nominator1 + nominator2) / (denominator1 + denominator2));
//		}
//
//		// Update the correlation message for sys node 1.
//		sideInfo[index] = log(((1 - p) * exp(LLRFinal2) + p) / (1 - p + p * exp(LLRFinal2)));
//	

	// Make decision.
	sysNodeInfo1.llrEstimation = LLRFinal1;
	sysNodeInfo2.llrEstimation = LLRFinal2;

	double p00 = 0.5 * (1- p) * p0FromRP * p0FromCoded1 * p0FromCoded2 * channelInfo[2];
	double p01 = 0.5 * p * p1FromRP * channelInfo[3] * p0FromCoded1 * p1FromCoded2;
	double p10 = 0.5 * p * p1FromRP * channelInfo[3] * p1FromCoded1 * p0FromCoded2;
	double p11 = 0.5 * (1 - p) * p2FromRP * p1FromCoded1 * p1FromCoded2 * channelInfo[4];

	if(p00 >= p01 && p00 >= p10 && p00 >= p11) {
		sysNodeInfo1.decision = 0;
		sysNodeInfo2.decision = 0;
	} else if(p01 >= p00 && p01 >= p10 && p01 >= p11) {
		sysNodeInfo1.decision = 0;
		sysNodeInfo2.decision = 1;
	} else if(p10 >= p00 && p10 >= p01 && p10 >= p11) {
		sysNodeInfo1.decision = 1;
		sysNodeInfo2.decision = 0;
	} else if(p11 >= p00 && p11 >= p01 && p11 >= p10) {
		sysNodeInfo1.decision = 1;
		sysNodeInfo2.decision = 1;
	} else {
		cout << "Wrong decision!" << endl;
	}


/*	
	if(LLRFinal1 >= 0)
		sysNodeInfo1.decision = 0;
	else
		sysNodeInfo1.decision = 1;
	
	if(LLRFinal2 >= 0)
		sysNodeInfo2.decision = 0;
	else
		sysNodeInfo2.decision = 1;
	*/

	// Record some information.
	if(iteration > 2 && (sysNodeInfo1.decision != sys1 || sysNodeInfo2.decision != sys2)) {
		proTrack << "Iteration: " << iteration << endl;
		proTrack << "Index: " << index 
						 << " (sys1, sys2):" << "(" << sys1 << "," << sys2 << ")" << endl;
		proTrack << "Decision: " << "(" << sysNodeInfo1.decision << "," 
						 << sysNodeInfo2.decision << ")" << endl;
		proTrack << "Final (p00, p01, p10, p11): (" 
						 << p00 << "," << p01 << "," << p10 << "," << p11 << ")" << endl;
		proTrack << "(p0FromRP, p1FromRP, p2FromRP): (" 
						 << p0FromRP << "," << p1FromRP << "," << p2FromRP << ")" << endl;
		proTrack << "From LDGM1 (p0, p1): (" << p0FromCoded1 << "," 
						 << p1FromCoded1 << ")" << endl;
		proTrack << "From LDGM2 (p0, p1): (" << p0FromCoded2 << "," 
					   << p1FromCoded2 << ")" << endl;
		proTrack << "Channel information: " 
						 << "(" << channelInfo[2] << "," << channelInfo[3] << ","
						 << channelInfo[4] << ")" << endl;
		proTrack << "LLREstimation: (" << LLRFinal1 << "," << LLRFinal2 << ")" << endl;
		proTrack << endl;
	}
	
	return outPdf;
}

// For pure RCM.
vector<vector<double> > syntheticDecoderSys(Sys_info &sysNodeInfo1, Sys_info &sysNodeInfo2, vector<vector<vector<double> > >& pdfMessageFromRP, vector<double>& channelInfo, double p, int index, ofstream &proTrack, int sys1, int sys2, int iteration) {

	int degree = sysNodeInfo1.nodeData.degree;

	vector<vector<double> > inPdfMessage(degree);

	// Get incoming messages.
	for(int i=0; i<degree; i++) {
		int neiNum = sysNodeInfo1.nodeData.neighbourNum[i];
		int messageIndex = sysNodeInfo1.nodeData.inmessageIndex[i];
		inPdfMessage[i] = pdfMessageFromRP[neiNum][messageIndex];
	}

	double p0FromRP = 1;
	double p1FromRP = 1;
	double p2FromRP = 1;

	for(int i=0; i<degree; i++) {
		p0FromRP *= inPdfMessage[i][2];
		p1FromRP *= inPdfMessage[i][3];
		p2FromRP *= inPdfMessage[i][4];
	}

	// Probability from RP
	double pTotalFromRP = p0FromRP + p1FromRP + p2FromRP;
	p0FromRP /= pTotalFromRP;
	p1FromRP /= pTotalFromRP;
	p2FromRP /= pTotalFromRP;

	// Probability in total and normalization.
	double p0Final = p0FromRP * channelInfo[2] * 0.5 * (1 - p);
	double p1Final = p1FromRP * channelInfo[3] * p;
	double p2Final = p2FromRP * channelInfo[4] * 0.5 * (1 - p);
	
	double pTotalFinal = p0Final + p1Final + p2Final;
	p0Final /= pTotalFinal;
	p1Final /= pTotalFinal;
	p2Final /= pTotalFinal;
	

	// Compute the message to the synthetic RP nodes
	vector<vector<double> > outPdf(degree, vector<double>(5, 0));
	for(int i=0; i<degree; i++) {
		outPdf[i][2] = p0Final / inPdfMessage[i][2];
		outPdf[i][3] = p1Final / inPdfMessage[i][3];
		outPdf[i][4] = p2Final / inPdfMessage[i][4];
		vecNorm(outPdf[i]);
	}

	/*
	if(p0Final > p1Final && p0Final > p2Final) {
		sysNodeInfo1.decision = 0;
		sysNodeInfo2.decision = 0;
	} else if(p2Final > p0Final && p2Final > p1Final) {
		sysNodeInfo1.decision = 1;
		sysNodeInfo2.decision = 1;
	} else {
		int random = rand() % 100;
		if(random < 50) {
			sysNodeInfo1.decision = 0;
			sysNodeInfo2.decision = 1;
		} else {
			sysNodeInfo1.decision = 1;
			sysNodeInfo2.decision = 0;
		}
	}*/

	// ----- To use LLR to make decision
	// Start to calculate the message to coded nodes and make decision.
	double channelLLR = log((channelInfo[2] + 0.5 * channelInfo[3]) /
												  (channelInfo[4] + 0.5 * channelInfo[3]));
	double LLRFromRP = 0;
	for(int i=0; i<degree; i++) {
		LLRFromRP += log((inPdfMessage[i][2] + 0.5 * inPdfMessage[i][3]) /
										 (inPdfMessage[i][4] + 0.5 * inPdfMessage[i][3]));
	}	
	
	double LLRFinal1 = channelLLR + LLRFromRP;
	double LLRFinal2 = channelLLR + LLRFromRP;

	// Make decision.
	sysNodeInfo1.llrEstimation = LLRFinal1;
	sysNodeInfo2.llrEstimation = LLRFinal2;

	if(LLRFinal1 > 0)
		sysNodeInfo1.decision = 0;
	else
		sysNodeInfo1.decision = 1;
	
	if(LLRFinal2 > 0)
		sysNodeInfo2.decision = 0;
	else
		sysNodeInfo2.decision = 1;
	
	// Record some information.
	if(iteration > 2 && (sysNodeInfo1.decision != sys1 || sysNodeInfo2.decision != sys2)) {
		proTrack << "Iteration: " << iteration << endl;
		proTrack << "Index: " << index 
						 << " (sys1, sys2):" << "(" << sys1 << "," << sys2 << ")" << endl;
		proTrack << "Decision: " << "(" << sysNodeInfo1.decision << "," 
						 << sysNodeInfo2.decision << ")" << endl;
		proTrack << "Final (p0, p1, p2): (" 
						 << p0FromRP * channelInfo[2] << ","
						 << p1FromRP * channelInfo[3] << ","
						 << p2FromRP * channelInfo[4] << ")" << endl;
		proTrack << "(p0FromRP, p1FromRP, p2FromRP): (" 
						 << p0FromRP << "," << p1FromRP << "," << p2FromRP << ")" << endl;
		proTrack << "Channel information: " 
						 << "(" << channelInfo[2] << "," << channelInfo[3] << ","
						 << channelInfo[4] << ")" << endl;
		proTrack << endl;

	}
	
	return outPdf;
}
   
