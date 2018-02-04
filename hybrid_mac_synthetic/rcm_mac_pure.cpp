#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "node.h"
#include "sparsematrix.h"
#include "func.h"
#include "pcm.cpp"

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

	/* 
  //For weightset size of 4
  int horizontalseq[][2] = { {0, 1},
			     {1, 0}, 
			     {0, 1},
			     {1, 0} };
  */

	/*	 
  //For weightset size of 8
  int horizontalseq[][4] = { {2, 3, 1, 0},
			     {1, 0, 2, 3},
			     {3, 2, 0, 1},
			     {0, 1, 3, 2} };
  */
  
  //For weightset size of 16 
  /*
   int horizontalseq[][8] = { {7, 6, 5, 4, 3, 2, 1, 0},
			     {3, 2, 1, 0, 7, 6, 5, 4},
			     {7, 6, 5, 4, 3, 2, 1, 0},
			     {3, 2, 1, 0, 7, 6, 5, 4} };
  */
	int horizontalseq[][8] = { {7, 3, 6, 2, 5, 1, 4, 0},
			     {0, 4, 1, 5, 2, 6, 3, 7},
			     {4, 3, 5, 2, 6, 0, 7, 1},
			     {2, 6, 3, 7, 0, 5, 1, 4} };
	
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

double getSymPdf(vector<double> symbolSet, vector<double> estimate, double candiVal, double norFactor, double y, double sigma);
int countError(vector<Sys_info> &sysNodeInfo, vector<int> &systrans);
void test();


int main() 
{

	//test();
	//return 1;
   // Source info
  int sysNumber = 10000;
  int symbolNumber = 20000;
	
  double p1 = 0.5;    //sparsity of the source.
  double p0 = 1 - p1;
  double p = 0.01; // correlation between source 1 and source 2.

	double H_U1_givenU2 = -p*log2(p) - (1 - p) * log2(1 - p);
	double H_U2 = 1;
	double entropy = H_U1_givenU2 + H_U2;
	
	double codeRate = double(sysNumber) / (sysNumber + symbolNumber);
	double R = (H_U1_givenU2 + H_U2) * codeRate * 2;
	double throughput = codeRate * 2;
  
	double capacity = 10 * log10((pow(2, R) - 1) / (4 * codeRate));   // in terms of Eso/N0
	//double capacity = 10 * log10((pow(2, 2 * R) - 1) / 2);     // Channel capacity in terms of Es/N0.

  cout<<"Code rate is: "<<codeRate<<endl;
  cout<<"entropy is: "<<entropy<<endl;
  cout<<"Throughput is:"<<throughput<<endl;
  cout<<"channel capacity is: "<<capacity<<endl;


  //weightset
  int wval[] = {1, 1, 1, 1, 2, 2, 2, 2};
  //int wval[] = {2, 3, 7, 10};
  //int wval[] = {1, 2};
  int weightSize = 2 * sizeof(wval) / sizeof(int);
  vector<double> weightset(weightSize, 0);
  
  int j = 0;
  for(int i=0;i<weightSize;i=i+2)
    {
      weightset[i] = wval[j];
      weightset[i+1] = -wval[j];

      ++j;
    }

  vector<int> weightval(wval, wval+weightSize/2);
 
  // We need to know the value range of RP symbol and its distribuiton
	// so that we can calculate the normalization factor and for non-synthetic
	// decoding method, we also use the symbolset to do the decoding.
	map<double, double> symbolsetdata = getSymbolSet(weightset, p1);
  vector<double> symbolset;
  vector<double> symbolpro;
  
	map<double, double>::iterator it = symbolsetdata.begin();

	for(it; it!=symbolsetdata.end(); ++it) {
      symbolset.push_back( it->first ); 
      symbolpro.push_back( it->second );
  }
  
  double symbol_entropy = 0;
  for(int i=0;i<symbolpro.size();++i)
    symbol_entropy += (-symbolpro[i] * log2(symbolpro[i]));

  cout << "Symbol entropy:" << symbol_entropy << endl;

	// In order to do the decoding sythetically, we have to know all the 
	// possbile values of synthetic nodes, which is the sum of two RP symbol
	// values.
	set<int> symSumSet;
	for(int i=0; i<symbolset.size(); i++) {
		for(int j=0; j<symbolset.size(); j++) {
			symSumSet.insert(symbolset[i] + symbolset[j]);
		}
	}
   
  // Generator RCM matrix for two RCM encoders
  // Determine how many G0 are needed to build G
  int G0RowNum = sysNumber / (weightSize / 2) / 2 * 4;
  int G0Num = 1;
  while(G0RowNum * G0Num < symbolNumber)
     ++G0Num;
  
  double gmatrixSeed = 200;
  Matrix gmatrix1 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
  for(int i=0;i<G0Num;++i)
    {
      Matrix gmatrix0 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
      gmatrix1.matrixStackDown(gmatrix0);
    }

	Matrix gmatrix2 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
	for(int i=0;i<G0Num;++i) {
		Matrix gmatrix0 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
		gmatrix2.matrixStackDown(gmatrix0);
    }


  	// Systematic nodes.
  vector<Sys_info> sysNodeInfo1(sysNumber);
  vector<Sys_info> sysNodeInfo2(sysNumber);
  
	// RCM nodes.
	vector<Coded_info> symNodeInfo1(symbolNumber);
  vector<Coded_info> symNodeInfo2(symbolNumber);

  gmatrix1.organizebyCol();
  for(int index=0;index<sysNumber;++index) {
    sysNodeInfo1[index].readLinearMatrix(gmatrix1, index, symbolNumber, 1);
    sysNodeInfo2[index].readLinearMatrix(gmatrix1, index, symbolNumber, 1);
	}

  for(int index=0;index<symbolNumber;++index) {
    symNodeInfo1[index].readLinearMatrix(gmatrix1, index);
    symNodeInfo2[index].readLinearMatrix(gmatrix1, index);
  }
 
	// Each Sys node get the storage index of its connections.
  for(int index=0;index<sysNumber;++index) {
    sysNodeInfo1[index].getsourceNum(symNodeInfo1, index);
    sysNodeInfo2[index].getsourceNum(symNodeInfo2, index);
  } 

	// Each RP symbol node get the storage index of its connections.  
  for(int index=0;index<symbolNumber;++index) {              
    symNodeInfo1[index].getsourceNum(sysNodeInfo1, index, 1);
    symNodeInfo2[index].getsourceNum(sysNodeInfo2, index, 1);
  }


	// Calculate the average value of transimitted symbols
	double Es = 2;
	double Esender = Es / 2;
	double Es1 = Esender, Es2 = Esender;

	double Eso = Es / (4 * codeRate);

	cout << "energy per source bit is:" << Eso << endl;

  // Start to generate sys bits and transmitted symbols
  vector<int> systematicBits1(sysNumber, 0);
  vector<int> systematicBits2(sysNumber, 0);
  vector<double> rpsymbols1(symbolNumber, 0);
  vector<double> rpsymbols2(symbolNumber, 0);
  
  
  // Our assumption is that the generated symbols are symmetric, i.e
  // the number of possible symbols is odd and symmetric to zero.
  // The length of the vector should be determined by the largest value
  // of the generated symbols.
  int sympdfsize = symbolset[symbolset.size()-1] * 2 + 1;
  
  int symSumPdfSize = symbolset[symbolset.size() - 1] * 2 * 2 + 1; 
  
	// PDF message from channel(synthetic symbols)
	vector< vector<double> > channelPdf(symbolNumber, vector<double>(symSumPdfSize, 0));        
  
	// Channel messages for systematic nodes.
	vector<vector<double> > sysChannelOutput(sysNumber, vector<double>(5, 0));  

	// Channel information.
	double gap = 1;
	double snr = capacity + gap;
	double sigma = sqrt( Eso / (2 * pow(10, snr/10)) ); 
	double noiseVar = pow(sigma, 2);

	cout << "Gap to the limit is: " << gap << endl;
	cout << "SNR is: " << snr << endl;

	// Calculate accumulated ber over blocks.
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

  double norFactor = getNormalizationFactorQAM(symbolset, symbolpro);

	// Because of synthetic decoding, messages passed between synthetic
	// RP nodes and synthetic sys nodes are pdfs. In order to facilitate
	// the implementation, these messages are stored outside of nodes.
	// Since sys node degree to RP symbol is not all the same, so the pdf
	// message should be initialized individually.
	vector<vector<vector<double> > > pdfMessageFromRP(symbolNumber, vector<vector<double> >(weightSize, vector<double>(5, 0)));
	vector<vector<vector<double> > > pdfMessageFromSys;
	for(int i=0; i<sysNumber; i++) {
		int curDegree = sysNodeInfo1[i].nodeData.degree;
		pdfMessageFromSys.emplace_back(vector<vector<double> >(curDegree, vector<double>(5, 0)));
	}

	// If the random generated number is smaller than the threshold, sysbit2 would be different from sysbit 1
	// while on the other hand, sysbit2 would be the same as sysbit 1.
  int sysDiffThreshold = p * 1000000;
  
	const int blockNum = 5000;
  int currentBlock = 0;

	while(currentBlock < blockNum) {

    ++currentBlock;

    error<<"Block "<<currentBlock<<":";
    
    srand(seed++);

    // Generate souce bits; Use identification bit
		int count = 0;
    systematicBits1 = sourceGenerator(sysNumber, p1, seed++);
    systematicBits1[0] = 1;

    for(int i=0; i<systematicBits1.size(); i++) {
      int random = rand() % 1000000;
      if(random < sysDiffThreshold) {
				systematicBits2[i] = 1 ^ systematicBits1[i];
				count++;
			}
      else
				systematicBits2[i] = systematicBits1[i];
    }
    systematicBits2[0] = 0;

		cout << "Number of different bits is: " << count << endl;

    vector<int> systrans1(systematicBits1);
    vector<int> systrans2(systematicBits2);

    // RP symbol encoder
    linearEncoder(symNodeInfo1, systematicBits1, rpsymbols1);
    linearEncoder(symNodeInfo2, systematicBits2, rpsymbols2);

    // Channel symbol
    vector<complex<double> > channelSignals;
    
    // Generate channel symbol for systematic bits. Notice that we for sender i,
		// we have to satisfy the requirement that the average energy is Es_i, which
		// is 1 here, for conveniece.
    for(int i=0; i<sysNumber; i=i+2) {
    	double tempAX = (2 * systematicBits1[i] - 1) / sqrt(2) * sqrt(Es1);
			double tempAY = (2 * systematicBits1[i+1] - 1) / sqrt(2) * sqrt(Es1);
			complex<double> digitSignalA = complex<double>(tempAX, tempAY);

			double tempBX = (2 * systematicBits2[i] - 1) / sqrt(2) * sqrt(Es2);
			double tempBY = (2 * systematicBits2[i+1] - 1) / sqrt(2) * sqrt(Es2);
			complex<double> digitSignalB = complex<double>(tempBX, tempBY);
      channelSignals.push_back(digitSignalA + digitSignalB);
		} 

    // Generate channel symbol for RP symbol, QAM. Since normalized RP symbol point energy is 1 and now each assigned point energy
    // is Es/2, normalized factor should be adjusted.
    vector<complex<double> > tempChannelSignals1 = modulatorQAM(rpsymbols1, norFactor);
    vector<complex<double> > tempChannelSignals2 = modulatorQAM(rpsymbols2, norFactor);
    for(int i=0; i<tempChannelSignals1.size(); i++)
      channelSignals.push_back(tempChannelSignals1[i] + tempChannelSignals2[i]);
   

		// Send symbols through complex AWGN channel
    complexAWGNChannel(channelSignals, sigma);                 

		// Initialization for sys nodes
		int wrongChannelCount = 0;
    for(int index=0;index<sysNumber;++index) {

      double y;
      if(index % 2 == 0)
				y = channelSignals[index/2].real();
      else
				y = channelSignals[index/2].imag();

      // 0,0 --> -2; 0, 1 or 1, 0 --> 0; 1, 1 --> 2.
			double pSumM2 =  gaussianFunc(y, (-2) / sqrt(2) * sqrt(Esender), sigma) * (1 - p);            
      double pSum0 = gaussianFunc(y, 0, sigma) * p;    
      double pSum2 = gaussianFunc(y, 2 / sqrt(2) * sqrt(Esender), sigma) * (1 - p);     

			sysChannelOutput[index] = vector<double> {0, 0, pSumM2, pSum0, pSum2};
			vecNorm(sysChannelOutput[index]);

			double maxPro = max(pSumM2, max(pSum0, pSum2));
			if(maxPro == pSumM2 && (systrans1[index] + systrans2[index] != 0))
				wrongChannelCount++;
			if(maxPro == pSum0 && (systrans1[index] + systrans2[index] != 1))
				wrongChannelCount++;
			if(maxPro == pSum2 && (systrans1[index] + systrans2[index] != 2))
				wrongChannelCount++;

			for(int d=0; d<sysNodeInfo1[index].nodeData.degree;d++)
				pdfMessageFromSys[index][d] = sysChannelOutput[index];

			sysNodeInfo1[index].clearMessage();
      sysNodeInfo2[index].clearMessage();
    }
    
		cout << "The number of wrong channel estimation is: " << wrongChannelCount << endl;

		// Calculate the pdf of corrupted synthetic RP symbols (Initialization for synthetic symbol nodes)    
		for(int index=0;index<symbolNumber;++index) {
			double y;
			if(index % 2 == 0)
	  		y = channelSignals[(sysNumber + index) / 2].real();
			else
	  		y = channelSignals[(sysNumber + index) / 2].imag();

			vector<double> symdis(symSumPdfSize);
      
			for(set<int>::iterator it=symSumSet.begin(); it!=symSumSet.end(); it++) {
	    	int zepoint = (symSumPdfSize - 1) / 2;
	    	int symPos = *it + zepoint;
	    	symdis[symPos] =  gaussianFunc(y, (*it) * norFactor, sigma); 
	  	}
			vecNorm(symdis);
			channelPdf[index] = symdis;
			symNodeInfo1[index].clearMessage();
			symNodeInfo2[index].clearMessage();
    }

   	// side information
    vector<double> sideInfo(sysNumber, 0);
    double csiThreshold = 0;
 
    int iterationtime = 200;
    int currentTime =0;  
    double watchllr[sysNumber];

    vector<int> predecision1(sysNumber, 0);
    vector<int> predecision2(sysNumber, 0);
    int sametime = 0;
		ofstream proTrack("pure_rcm_proTrack.txt", ios::out);
    while(currentTime != iterationtime+1) {

      ++currentTime;
     
			for(int i=0; i<symbolNumber; i++)
				pdfMessageFromRP[i] = syntheticDecoderRP(symNodeInfo1[i], symNodeInfo2[i], i, pdfMessageFromSys, channelPdf[i], norFactor);
			
			for(int i=0; i<sysNumber; i++) {
				pdfMessageFromSys[i] = syntheticDecoderSys(sysNodeInfo1[i], sysNodeInfo2[i], pdfMessageFromRP, sysChannelOutput[i], p, i, proTrack, systrans1[i], systrans2[i], currentTime);
			}
			
			// Count error number.
      double errorbit = 0;
      
      if(sysNodeInfo1[0].decision == 1) {
				errorbit += countError(sysNodeInfo1, systrans1);
				errorbit += countError(sysNodeInfo2, systrans2);
      } else if(sysNodeInfo1[0].decision == 0) {
				errorbit += countError(sysNodeInfo1, systrans2);
				errorbit += countError(sysNodeInfo2, systrans1);
      } 
	   
      cout<<"error number is:"<<errorbit<<"("<<currentTime<<")"<<endl;

      error<<errorbit<<" ";      

      double diffPos = 0;
      for(int i=0; i<sysNumber; i++) {
				if(sysNodeInfo1[i].decision != sysNodeInfo2[i].decision)
	  			diffPos++;
      }
      cout << "Number of different decisions: " << diffPos << endl;

      if(sameDecision(predecision1, sysNodeInfo1) && sameDecision(predecision2, sysNodeInfo2))
				++sametime;
      else
				sametime=0;

      if(sametime == 5) {
				cout<<" Make same decision for 5 consecutive iterations !"<<endl;
				break;
      }

    }
  
    double errorbit = 0;
    
    if(sysNodeInfo1[0].decision == 1) {
      errorbit += countError(sysNodeInfo1, systrans1);
      errorbit += countError(sysNodeInfo2, systrans2);
    } else if(sysNodeInfo1[0].decision == 0) {
      errorbit += countError(sysNodeInfo1, systrans2);
      errorbit += countError(sysNodeInfo2, systrans1);
    } 
       

    double blockber = errorbit / (2 * sysNumber);
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

    cout<<" BER for SNR = "<<snr<<" db is "<<blockber<<"  sum:"<<bersum<<"("<<currentBlock<<")  "<<endl;
    cout<<endl;

  }
  double sum = 0;
  for(int i=0;i<error_record.size();++i)
      sum += error_record[i];
     
  error_ave = sum / error_record.size();

  error << "The number of blocks that have errors: " << error_record.size()<<endl;
  error << "The average error number is: " << error_ave << endl;
  error.close();
  erroreachblock.close();
} 

vector<int> interleaver(vector<double> &source, double percentage, vector<Sys_info> &sysNode, vector<Coded_info> &codedNodeP, int codedNoP, int group)
{

  int interleavingNum=int(codedNoP*percentage);

  vector<int> newPos(codedNoP);

  for(int index=0; index<interleavingNum; index++)
    newPos[index] = index;


  //start to generate random position for para coded bits.
  srand(500);
  for(int i=interleavingNum;i!=1;--i)
    {
      int random = rand() % i;

      int intervalue=newPos[random];
      newPos[random]=newPos[i-1];
      newPos[i-1]=intervalue;
    }



  //start to interleave the para coded bits.
  vector<double> interSource(interleavingNum,0);
  for(int index=0;index<interleavingNum;++index)
    {

      int newPosition = newPos[index];
      int oldPosition = index;

      interSource[newPosition] = source[oldPosition];
    }

  for(int index=0; index<interleavingNum; ++index)
    source[index]=interSource[index];

  //adjust data in the correspoding node
  int degree=codedNodeP[0].nodeData.degree;
  vector<Coded_info> interCoded(interleavingNum, Coded_info(degree));

  for(int index=0;index<interleavingNum;++index)
    {

      int oldNum=index;
      int newNum=newPos[index];

      //first, find the sysNode Num that connects to the original coded Node
      for(int i=0;i<degree;++i)
	{
	  int sysNodeNum=codedNodeP[oldNum].nodeData.neighbourNum[i];
	  int coded_index=codedNodeP[oldNum].nodeData.inmessageIndex[i];
	  
	  if(group==1)
	    sysNode[sysNodeNum].nodeData.neighbourNum[coded_index]=newNum;
	  if(group==2)
	    sysNode[sysNodeNum].nodeData.para_neighbourNum[coded_index]=newNum;
	}

      interCoded[newNum]=codedNodeP[oldNum];
    }

  for(int index=0;index<interleavingNum;++index)
    codedNodeP[index]=interCoded[index];
 
  return newPos;
}

void deinterleaver(double percentage, int *newPos, vector<Sys_info> &sysNode, vector<Coded_info> &codedNodeP, int codedNoP, int group)
{


  int interleavingNum=int(codedNoP * percentage);

  /*
  vector<int> origin(interleavingNum,0);

  // deinterleave the source
  for(int index=0;index<interleavingNum;++index)
    {
      int originPos=newPos[index];
      origin[originPos]=source[index];
    }
  */

  int degree=codedNodeP[0].nodeData.degree;
  vector<Coded_info> interCoded(interleavingNum, Coded_info(degree));

  //deinterleave the data in the sys node
  for(int index=0;index<interleavingNum;++index)
    {
      int originNum=index;
      int newNum=newPos[index];

      for(int neighbour=0;neighbour<degree;++neighbour)
	{
	  int sysNum=codedNodeP[newNum].nodeData.neighbourNum[neighbour];
	  int codedIndex=codedNodeP[newNum].nodeData.inmessageIndex[neighbour];

	  if(group==1)
	    sysNode[sysNum].nodeData.neighbourNum[codedIndex]=originNum;
	  if(group==2)
	    sysNode[sysNum].nodeData.para_neighbourNum[codedIndex]=originNum;
	}

      interCoded[originNum]=codedNodeP[newNum];
    }


  for(int index=0;index<interleavingNum;++index)
    codedNodeP[index]=interCoded[index];
}

int countError(vector<Sys_info> &sysNodeInfo, vector<int> &systrans) {

  int errorNum = 0;
  for(int i=0; i<systrans.size(); i++) {
    if(sysNodeInfo[i].decision != systrans[i])
      errorNum++;
  }
  
  return errorNum;
}

void test() {

	// Test of deconvolve 2 function.
	vector<double> pdf1 = {0, 0, 0.1, 0.5, 0.4};
	vector<double> pdf2 = {0, 0, 0.3, 0.5, 0.2};
	vector<double> pdf3 = {0, 0, 0.4, 0.5, 0.1};

	double weight = -2;
	vector<double> weightedPdf1(pdf1);
	weightPdf(weight, weightedPdf1);

	vector<double> sumPdf(1, 1);
	
	conv(weightedPdf1, pdf2, sumPdf);
	conv(sumPdf, pdf3, sumPdf);
	vecNorm(sumPdf);

	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl;
	
	deconvolve2(sumPdf, pdf1, weight, sumPdf);
	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl;

	
	conv(pdf2, pdf3, sumPdf);
	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl;

	/*
	conv(sumPdf, pdf1, sumPdf);
	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl;
*/

	/*
	// Test of deconvolve function.
	pdf1 = {0, 0.6, 0.4};
	pdf2 = {0, 0.2, 0.8};
	pdf3 = {0, 0.3, 0.7};

	conv(pdf1, pdf2, sumPdf);
	conv(sumPdf, pdf3, sumPdf);
	vecNorm(sumPdf);

	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl;

	deconvolve(sumPdf, pdf1, 1, sumPdf);
	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl; */

}
