#include<iostream>
#include<vector>
#include<map>
#include<stdlib.h>
#include<cmath>
#include<fstream>
#include<algorithm>

#include"node.h"
#include"sparsematrix.h"
#include"function.h"

using namespace std;

void printPdf(const vector<double> &pdf)
{
  for(int i=0;i<pdf.size();++i)
    cout<<pdf[i]<<" ";
  cout<<endl;
};

void printMatrix(Matrix &ma);

void horizontalStack(vector<Matrix> &ma, int stacksequen[], int &seed, Matrix &result);

Matrix consEleMatrix(int row, int column, int value)
{
  Matrix elematrix(row, column);

  for(int i=0;i<row;++i)
    {
      int j= 2 * i;
      ElementData ele(i, j, value);
      elematrix.readData(ele);

      ++j;
      ele.setValue(i, j, -1 * value);
      elematrix.readData(ele);
    }
  
  return elematrix;
}


//Generate a N/2 x N  matrix with symbol degree 8
Matrix linearMatrixConsD8(int N, vector<int> &weight, int seed)
{

  int weightpair = weight.size();

  int elematrixcolumn = N / weightpair;
  int elematrixrow = elematrixcolumn / 2;
  
  int index = 0;

  vector<Matrix> ccm_matrix(weightpair);

  //Construct elementary matrix for each weight value
  for(int i=0;i<weightpair;++i)
    ccm_matrix[i] = consEleMatrix(elematrixrow, elematrixcolumn, weight[i]);


 	 
  //For weightset size of 8
  int horizontalseq[][4] = { {2, 3, 1, 0},
			     {1, 0, 2, 3},
			     {3, 2, 0, 1},
			     {0, 1, 3, 2} };
  
  /*
  //For weightset size of 16
   int horizontalseq[][8] = { {7, 6, 5, 4, 3, 2, 1, 0},
			     {3, 2, 1, 0, 7, 6, 5, 4},
			     {7, 6, 5, 4, 3, 2, 1, 0},
			     {3, 2, 1, 0, 7, 6, 5, 4} };
  */
  //Four row elementary matrices group forms a whole matrix
  Matrix matrixhori[4]; 

  for(int i=0;i<4;++i)
    horizontalStack(ccm_matrix, horizontalseq[i], seed, matrixhori[i]);
  
  Matrix finalmatrix = matrixhori[0];
  for(int i=1;i<4;++i)
    finalmatrix.matrixStackDown(matrixhori[i]);

  return finalmatrix;

}

void horizontalStack(vector<Matrix> &ma, int stacksequen[], int &seed, Matrix &result)
{
  int manum = ma.size();

  for(int i=0;i<manum;++i)
    ma[i].columnPermutation(seed++);
  
  result = ma[stacksequen[0]];
  for(int i=1;i<manum;++i)
    result.matrixStackRight(ma[stacksequen[i]]);
  
}

int main() 
{
  // Source info
  int sysNumber = 10000;
  int symbolNumber = 8000;
  int codedNumber = 2000;
  int sysDegree = 3;
  
  double p1 = 0.5;    //sparsity of the source.
  double p0 = 1 - p1;

  double entropy = -p0 * log2(p0) - p1 * log2(p1);
  double coderate = double(sysNumber) / (symbolNumber + codedNumber);
  double inforate = entropy * coderate * 2;
  double throughput = coderate * 2;
  /********/
  double capacity = 10 * log10( pow(2, inforate) - 1);     // Channel capacity in terms of Es/N0.

  cout<<"Code rate is: "<<coderate<<endl;
  cout<<"entropy is: "<<entropy<<endl;
  cout<<"Throughput is:"<<throughput<<endl;
  cout<<"channel capacity is: "<<capacity<<endl;

  // The indicator for the message passing method:
  // 1: Accurate message passing; 2: Simplified decoding.
  int messageMethod = 2;

  //weightset
  //int wval[] = {1, 1, 1, 1, 2, 2, 2, 2};
  int wval[] = {2, 3, 7, 10};
  
  int weightsize = 2 * sizeof(wval) / sizeof(int);
  vector<double> weightset(weightsize, 0);
  
  int j = 0;
  for(int i=0;i<weightsize;i=i+2)
    {
      weightset[i] = wval[j];
      weightset[i+1] = -wval[j];

      ++j;
    }

  vector<int> weightval(wval, wval+weightsize/2);
 
  map<int, double> symbolsetdata = getSymbolSet(weightset, p1, messageMethod == 1 ? 0 : 0);
  map<int, vector<vector<int> > > mappingTable = getMappingTable(weightset);
  vector<double> symbolset;
  vector<double> symbolpro;
  
   map<int, double>::iterator it = symbolsetdata.begin();

   for(it; it!=symbolsetdata.end(); ++it) {
      symbolset.push_back( it->first ); 
      symbolpro.push_back( it->second );
    }
  
   double symbol_entropy = 0;
   for(int i=0;i<symbolpro.size();++i)
     symbol_entropy += (-symbolpro[i] * log2(symbolpro[i]));

   cout << "Symbol entropy:" << symbol_entropy << endl;

   // RCM encoder
   // Determine how many G0 are needed to build G
   int G0RowNum = sysNumber / (weightsize / 2) / 2 * 4;
   int G0Num = 1;
   while(G0RowNum * G0Num < symbolNumber)
     ++G0Num;
  
   double gmatrixSeed = 200;
  Matrix gmatrix = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
  for(int i=0;i<G0Num;++i) {
    Matrix gmatrix1 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
    gmatrix.matrixStackDown(gmatrix1);
  }
  vector<Sys_info> sysNodeInfo(sysNumber);   // Establish the data for systematic bits
  vector<Coded_info> symNodeInfo(symbolNumber);         // Establish the data for RP symbols

  gmatrix.organizebyCol();
  for(int index=0;index<sysNumber;++index)
    sysNodeInfo[index].readLinearMatrix(gmatrix, index, symbolNumber, 1);
  for(int index=0;index<symbolNumber;++index)
    symNodeInfo[index].readLinearMatrix(gmatrix, index);

 
  // LDGM encoder
  vector<multimap<int, int> > madata = paritycm(sysNumber, codedNumber, sysDegree, sysNumber, 1500);
  multimap<int, int> onePositionbyRow(madata[0]);
  multimap<int, int> onePositionbyColumn(madata[1]);

  vector<Coded_info> codedNodeInfo(codedNumber);

  for(int index=0;index<sysNumber;++index)
    sysNodeInfo[index].readNeighbourNum(onePositionbyRow, index, 2);
  for(int index=0;index<codedNumber;++index)
    codedNodeInfo[index].readNeighbourNum(onePositionbyColumn, index);
  
  // Get neighbor index.
  for(int index=0;index<symbolNumber;++index)
    symNodeInfo[index].getsourceNum(sysNodeInfo, index, 1);
  for(int index=0;index<codedNumber;++index)
    codedNodeInfo[index].getsourceNum(sysNodeInfo, index, 2);
  for(int index=0;index<sysNumber;++index)
    sysNodeInfo[index].getsourceNum(symNodeInfo, codedNodeInfo, index);

  // Calculate the average value of transimitted symbols
  double Es = 1;
  double Eb = Es / inforate;

  cout<<"average energy is:"<<Es<<endl;

  // Start to generate sys bits and transmitted symbols
  vector<int> systematicbits(sysNumber,0);
  vector<double> rpsymbols(symbolNumber, 0);
  vector<int> codedbits(codedNumber, 0);

  // Our assumption is that the generated symbols are symmetric, i.e
  // the number of possible symbols is odd and symmetric to zero.
  // The length of the vector should be determined by the largest value
  // of the generated symbols.
  int sympdfsize = 0;
  // The following is under the contidition that the generated symbols are balanced
  sympdfsize = symbolset[symbolset.size()-1] * 2 + 1;   

  vector<double> pdfinitializer(sympdfsize,0);
  vector< vector<double> > channelPdf(symbolNumber, pdfinitializer);    // PDF message from channel(symbols)
  vector<double> sysChannelOutput(sysNumber, 0);                       // LLR message from channel(sysbits)
  vector<double> codedChannelOutput(codedNumber, 0);                   // LLR message from channel(codedbits)

	double gap = 1.5;
  double snr = 8;
  double sigma = sqrt( Es / (2 * pow(10, snr/10)) ); 
  double noiseVar = pow(sigma, 2);

  const int blocknum = 5000;
  int currentblock = 0;

  double bersum = 0;
  
  int seed = 2000;

  ofstream error;
  //error.open("res_hybrid_w4142_6000_7_01_Rc1_2p3db.txt", ios::out);
  error.open("test.txt", ios::out);

  ofstream erroreachblock;
  //erroreachblock.open("err_m4_2500_5_001_t10_2p5db.txt", ios::out);
  erroreachblock.open("test1.txt", ios::out);

  vector<double> error_record;
  double error_min = sysNumber;
  double error_max = 0;
  double error_ave;

  double norfactor = getNormalizationFactorQAM(symbolset, symbolpro);
  vector<double> normalizedSymbolset(symbolset);
  for(int i=0;i<symbolset.size();++i)
      normalizedSymbolset[i] *= norfactor;
 
  // Allocate space for outgoing mean and variance if simplified method
  // is used.
  if(messageMethod == 2) {
    for(int i=0; i<sysNumber; i++)  sysNodeInfo[i].allocateForMeanAndVar();
    for(int i=0; i<symbolNumber; i++)  symNodeInfo[i].allocateForMeanAndVar();
  }

  while(currentblock < blocknum) {

    ++currentblock;

    error<<"Block "<<currentblock<<":";
    
    srand(seed++);

    // Generate souce bits
    systematicbits = sourceGenerator(sysNumber, p1, seed);
    
    vector<int> systrans(systematicbits);
    // RCM encoding
    linearEncoder(symNodeInfo, systematicbits, rpsymbols, 0, messageMethod == 1 ? 0 : 0);       
    encoder(codedNodeInfo, systematicbits, codedbits);           // Digital bits Encoder
    
    vector<complex<double> > channelSignals;
    // Modulate RP symbols using QAM
		vector<complex<double> > rpSignals = modulatorQAM(rpsymbols, norfactor);   // QAM and get the normalized factor.
		channelSignals.insert(channelSignals.end(), rpSignals.begin(), rpSignals.end());
    
    // Modulate coded bits using 4-QAM.
    for(int i=0;i<codedNumber;i=i+2)                                 
      channelSignals.push_back(complex<double>( (2*codedbits[i]-1) / sqrt(2), (2*codedbits[i+1]-1) / sqrt(2) ));

    // Send symbols through complex AWGN channel
    complexAWGNChannel(channelSignals, sigma);                 

    // Calculate the pdf of corrupted RP symbols (Initialization for symbol nodes)    
    for(int index=0;index<symbolNumber;++index) {
			double y;
			if(index % 2 == 0)
	  		y = channelSignals[index / 2].real();
			else
	  		y = channelSignals[index / 2].imag();

			vector<double> symdis(sympdfsize);
      
			for(int symindex=0;symindex<symbolset.size();++symindex) {
	    	int zepoint = (sympdfsize - 1) / 2;
	    	double symva = symbolset[symindex];
	    	int sympos = symva + zepoint;
	    	symdis[sympos] =  gaussianFunc(y, normalizedSymbolset[symindex], sigma);
	  	}
			vecNorm(symdis);
			channelPdf[index] = symdis;
			symNodeInfo[index].clearMessage();
		}

    // Initialization for coded nodes
    for(int index=0;index<codedNumber;++index) {
			double y;
			if(index % 2 == 0)
	  		y = channelSignals[(symbolNumber + index)/2].real();
			else
	  		y = channelSignals[(symbolNumber + index)/2].imag();

			codedChannelOutput[index] = (-sqrt(2) * y )/(pow(sigma,2));
			codedNodeInfo[index].clearMessage();
		}
	
    // Initialization for sys nodes
    for(int index=0;index<sysNumber;++index) {
				double initialmessage = log(p0/p1);
				sysChannelOutput[index] = initialmessage; 
				sysNodeInfo[index].clearMessage();
      }

    int iterationtime = 200;
    int currenttime =0;  
    double watchllr[sysNumber];

    vector<int> predecision(sysNumber, 0);
    int sametime = 0;
    while(currenttime != iterationtime+1) {

      ++currenttime;
     
      // Standard message calculation.
      if(messageMethod == 1) {
        for(int i=0;i<symbolNumber;++i)
				  symNodeInfo[i].computePdf(sysNodeInfo, channelPdf[i], norfactor, noiseVar, 1, i); 
        for(int i=0;i<codedNumber;++i)
				  codedNodeInfo[i].computeMessage(sysNodeInfo, codedChannelOutput[i], 2);
        for(int i=0;i<sysNumber;++i)
				  sysNodeInfo[i].computeMessage(symNodeInfo, codedNodeInfo, sysChannelOutput[i], 0);
      }
      // Simplified message passing.
      else {
        for(int i=0; i<symbolNumber; i++) {
          double rpObservation = channelSignals[i/2].real();
          if(i % 2 == 1)  rpObservation = channelSignals[i/2].imag();
          symNodeInfo[i].computeInMeanAndVar(sysNodeInfo, rpObservation, symbolsetdata, norfactor, noiseVar, 1, "llr");
        }
        for(int i=0; i<codedNumber; i++)
          codedNodeInfo[i].computeMessage(sysNodeInfo, codedChannelOutput[i], 2);
        for(int i=0; i<sysNumber; i++) 
          sysNodeInfo[i].computeMessage(symNodeInfo, codedNodeInfo, sysChannelOutput[i], 0);
      }

      double errorbit = 0;
      for(int i=0;i<sysNumber;++i) {
	  		if(sysNodeInfo[i].decision != systrans[i])
	    	++errorbit;
			}   

      cout<<"error number is:"<<errorbit<<"("<<currenttime<<")"<<endl;

      error<<errorbit<<" ";      

      if(sameDecision(predecision, sysNodeInfo))
				++sametime;
      else
				sametime=0;

      if(sametime == 5) {
				cout<<" Make same decision for 5 consecutive iterations !"<<endl;
				break;
      }

    }
    /*
    // record combined pdf for specific bit for specific node
    ofstream combinedPdfFile;
    //combinedPdfFile.open("Node100_W2_P05_37000_10000_25db_10pass.txt", ios::out);
    
    //vector<double> testCombinedPdf;
    int zeroPoint = (testCombinedPdf.size() - 1) / 2;
    int tempDe = testCombinedPdf[testCombinedPdf.size()-1];
    for(int i=0;i<testCombinedPdf.size()-1;i++)
      combinedPdfFile << i-zeroPoint << " " << testCombinedPdf[i]/tempDe << endl;
    combinedPdfFile.close();
    */
    double errorbit = 0;

    for(int i=0;i<sysNumber;++i) {
			if(sysNodeInfo[i].decision != systrans[i])
	  		++errorbit;
		} 
  
    double blockber = errorbit / sysNumber;
    bersum += blockber;

    if(errorbit != 0) {
      error_record.push_back(errorbit);
      erroreachblock << errorbit << endl;
    }      
    if(error_max < errorbit)
      error_max = errorbit;
    if(error_min > errorbit && errorbit != 0)
      error_min = errorbit;

    error<<"("<<error_max<<", "<<error_min<<") ";

    error<<"sum("<<bersum<<","<<error_record.size()<<")"<<endl;
    error<<endl;
    cout<<endl;

    cout<<" BER for SNR = "<<snr<<" db is "<<blockber<<"  sum:"<<bersum<<"("<<currentblock<<")  "<<endl;
    cout<<endl;

  }
  double sum = 0;
  for(int i=0;i<error_record.size();++i)
      sum += error_record[i];
     
  error_ave = sum / error_record.size();

  error<<"The number of blocks that have errors: "<<error_record.size()<<endl;
  error<<"The average error number is: "<<error_ave<<endl;
  error.close();
  erroreachblock.close();
} 

