#include<iostream>
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

  for(int index=0;index<pdf.size();++index)
    {
      double realvalue = index - midpoint_old; //the value of random variables
      double weightvalue = weight*realvalue;   //weighted value of random variables

      int newpos = weightvalue + midpoint_new;
            
      if( int(newpos) != newpos )             // In the case that weight is not an interger, only 
	continue;                             // calculate the probability for meaningful value, which
      else                                    // are all intergers.
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
  for(int vecindex=0;vecindex<candidates.size();++vecindex)
      {
	
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
  
  
  
    

