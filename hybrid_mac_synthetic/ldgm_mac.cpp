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
   int horizontalseq[][8] = { {7, 6, 5, 4, 3, 2, 1, 0},
			     {3, 2, 1, 0, 7, 6, 5, 4},
			     {7, 6, 5, 4, 3, 2, 1, 0},
			     {3, 2, 1, 0, 7, 6, 5, 4} };
  
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

int main() 
{
  // Source info
	int sysNumber = 10000;
	int codedNumber = 25000;
  int sysDegree = 10;
  
  double p1 = 0.5;    //sparsity of the source.
  double p0 = 1 - p1;
  double p = 0.01; // correlation between source 1 and source 2.

  double H_U1_givenU2 = -p*log2(p) - (1 - p) * log2(1 - p);
	double H_U2 = 1;
	double entropy = H_U1_givenU2 + H_U2;
	
  double codeRate = double(sysNumber) / (sysNumber + codedNumber);
  double R = (H_U1_givenU2 + H_U2) * codeRate;
 
 	double capacity = 10 * log10((pow(2, 2 * R) - 1) / (2 * R) * entropy / 2);  // in terms of Eso/N0.
	//double capacity = 10 * log10((pow(2, 2 * R) - 1) / 2);     // Channel capacity in terms of Es/N0.

  cout<<"Code rate is: "<<codeRate<<endl;
  cout<<"joint entropy is: "<<entropy<<endl;
  cout<<"Information rate is:"<<R<<endl;
  cout<<"channel capacity is: "<<capacity<<endl;


	// Generate LDGM generator matrices.
	vector<multimap<int, int> > maData1 = paritycm(sysNumber, codedNumber, sysDegree, sysNumber, 2000);
	multimap<int, int> onePositionByRow1(maData1[0]);
	multimap<int, int> onePositionByColumn1(maData1[1]);

	// Systematic nodes.
  vector<Sys_info> sysNodeInfo1(sysNumber);
  vector<Sys_info> sysNodeInfo2(sysNumber);
  
	// Coded nodes.
	vector<Coded_info> codedNodeInfo1(codedNumber);
	vector<Coded_info> codedNodeInfo2(codedNumber);

	for(int index=0;index<sysNumber;++index) {
		sysNodeInfo1[index].readNeighbourNum(onePositionByRow1, index, 1);
		sysNodeInfo2[index].readNeighbourNum(onePositionByRow1, index, 1);
	}
	
	for(int index=0; index<codedNumber; ++index) {
		codedNodeInfo1[index].readNeighbourNum(onePositionByColumn1, index);
		codedNodeInfo2[index].readNeighbourNum(onePositionByColumn1, index);
	}

	// Each Sys node get the storage index of its connections.
  for(int index=0;index<sysNumber;++index) {
		sysNodeInfo1[index].getsourceNum(codedNodeInfo1, index);
		sysNodeInfo2[index].getsourceNum(codedNodeInfo2, index);
	}
	
	// Each coded node get the storage index of its connections.
	for(int index=0; index<codedNumber; index++) {
		codedNodeInfo1[index].getsourceNum(sysNodeInfo1, index, 1);
		codedNodeInfo2[index].getsourceNum(sysNodeInfo2, index, 1);
	}
 
	// Calculate the average value of transimitted symbols
  double Es = 2;
  double Esender = Es / 2;
 	double Es1 = Esender, Es2 = Esender; 

  double Eso = Esender / codeRate;

  cout << "energy per source bit is:" << Eso << endl;

 
  // Start to generate sys bits and transmitted symbols
  vector<int> systematicBits1(sysNumber, 0);
  vector<int> systematicBits2(sysNumber, 0);
	vector<int> codedBits1(codedNumber, 0);
	vector<int> codedBits2(codedNumber, 0);
 		
	// Channel messages for systematic nodes.
	vector<double> sysChannelOutput(sysNumber, 0);  

	// Channel messages for coded nodes.
	vector<double> codedChannelOutput(codedNumber, 0);
	
	// Channel information.
	double gap = -1;
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

    // LDGM encoder
		encoder(codedNodeInfo1, systematicBits1, codedBits1);
		encoder(codedNodeInfo2, systematicBits2, codedBits2);

    // Channel symbol
    vector<double> channelSignals;
   
    // Generate channel symbol for systematic bits. Notice that we for sender i,
		// we have to satisfy the requirement that the average energy is Es_i, which
		// is 1 here, for conveniece.
    for(int i=0; i<sysNumber; i++) {
			double digitSignalA = (2 * systematicBits1[i] - 1) * Esender;
			double digitSignalB = (2 * systematicBits2[i] - 1) * Esender;
      channelSignals.push_back(digitSignalA + digitSignalB);
		} 

		// Generate channel symbols for coded bits. The process is the same as that in systematic bits.
		for(int i=0; i<codedNumber; i++) {
      double digitSignalA = (2 * codedBits1[i] - 1) * Esender;
			double digitSignalB = (2 * codedBits2[i] - 1) * Esender;
			channelSignals.push_back(digitSignalA + digitSignalB);
		}
		
    // Send symbols through complex AWGN channel
    AWGNChannel(channelSignals, sigma);                 

		// Initialization for sys nodes
    for(int index=0;index<sysNumber;++index) {

      double y = channelSignals[index];
      double v0 = log((exp(-pow(y + 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))) \
								+	exp(-pow(y, 2) / (2 * pow(sigma, 2)))) /	\
								(exp(-pow(y - 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))) \
								+ exp(-pow(y, 2) / (2 * pow(sigma, 2)))));
			sysChannelOutput[index] = v0;

			sysNodeInfo1[index].clearMessage();
      sysNodeInfo2[index].clearMessage();
    }
   
	 	
   	// Initialization for coded nodes.
		for(int index=0; index<codedNumber; ++index) {
			double y = channelSignals[sysNumber + index];
			double v0 = log((exp(-pow(y + 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))) \
								+	exp(-pow(y, 2) / (2 * pow(sigma, 2)))) /	\
								(exp(-pow(y - 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))) \
								+ exp(-pow(y, 2) / (2 * pow(sigma, 2)))));
			codedChannelOutput[index] = v0;
			
			codedNodeInfo1[index].clearMessage();
			codedNodeInfo2[index].clearMessage();
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
    while(currentTime != iterationtime+1) {

      ++currentTime;
      
						
			for(int i=0; i<codedNumber; i++)
				codedNodeInfo1[i].computeMessage(sysNodeInfo1, codedChannelOutput[i], 1);
		
			for(int i=0; i<sysNumber; i++)
				sysNodeInfo1[i].computeMessage(codedNodeInfo1, sysChannelOutput[i] , sideInfo[i]);
	
			for(int i=0; i<codedNumber; i++)
				codedNodeInfo1[i].computeMessage(sysNodeInfo1, codedChannelOutput[i], 1);

			// Update channel info
			for(int i=0; i<sysNumber; i++) {
				double y = channelSignals[i];
				
				double estimate = sysNodeInfo1[i].message_to_state;
				double nominator1 = exp(-pow(y,2) / (2*pow(sigma, 2))) / exp(estimate);
				double nominator2 = (exp(-pow(y+2*sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				double denominator1 = exp(-pow(y-2*sqrt(Esender), 2) / (2 * pow(sigma, 2))) / exp(estimate);
				double denominator2 = (exp(-pow(y, 2) / (2 * pow(sigma, 2))));

				if(isinf(nominator1) && isinf(denominator1)) {
					sysChannelOutput[i] = log(exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(-pow(y - 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				} else {
					sysChannelOutput[i] = log((nominator1 + nominator2) / (denominator1 + denominator2));
				}

			
				if(currentTime == 1) {
					for(int k=0; k<sysDegree; k++)
						sysNodeInfo2[i].setMessage(sysChannelOutput[i], k, 1);
				}
			}
			
			for(int i=0; i<codedNumber; i++) {
				double y = channelSignals[sysNumber + i];
				double estimate = codedNodeInfo1[i].message_to_state;
				double nominator1 = exp(-pow(y,2) / (2*pow(sigma, 2))) / exp(estimate);
				double nominator2 = (exp(-pow(y+2*sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				double denominator1 = exp(-pow(y-2*sqrt(Esender), 2) / (2 * pow(sigma, 2))) / exp(estimate);
				double denominator2 = (exp(-pow(y, 2) / (2 * pow(sigma, 2))));

				if(isinf(nominator1) && isinf(denominator1)) {
					codedChannelOutput[i] = log(exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(-pow(y - 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				} else {
					codedChannelOutput[i] = log((nominator1 + nominator2) / (denominator1 + denominator2));
				}

			}

			// Update side information.
			for(int i=0; i<sysNumber; ++i) {
				double decision = sysNodeInfo1[i].llrEstimation;

				double x = exp(decision);

				if(isinf(x))
					sideInfo[i] = log((1 - p) / p);
				else {
					sideInfo[i] = log(((1- p) * x + p) / (1 - p + p * x));
				}
			}

			for(int i=0; i<codedNumber; i++)
				codedNodeInfo2[i].computeMessage(sysNodeInfo2, codedChannelOutput[i], 1);

			for(int i=0; i<sysNumber; i++)
				sysNodeInfo2[i].computeMessage(codedNodeInfo2, sysChannelOutput[i] , sideInfo[i]);

			for(int i=0; i<codedNumber; i++)
				codedNodeInfo2[i].computeMessage(sysNodeInfo2, codedChannelOutput[i], 1);

			// Update channel info
			for(int i=0; i<sysNumber; i++) {
				double y = channelSignals[i];
				
				double estimate = sysNodeInfo2[i].message_to_state;
				double nominator1 = exp(-pow(y,2) / (2*pow(sigma, 2))) / exp(estimate);
				double nominator2 = (exp(-pow(y+2*sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				double denominator1 = exp(-pow(y-2*sqrt(Esender), 2) / (2 * pow(sigma, 2))) / exp(estimate);
				double denominator2 = (exp(-pow(y, 2) / (2 * pow(sigma, 2))));

				if(isinf(nominator1) && isinf(denominator1)) {
					sysChannelOutput[i] = log(exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(-pow(y - 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				} else {
					sysChannelOutput[i] = log((nominator1 + nominator2) / (denominator1 + denominator2));
				}

			}

			for(int i=0; i<codedNumber; i++) {
				double y = channelSignals[sysNumber + i];
				double estimate = codedNodeInfo2[i].message_to_state;
				double nominator1 = exp(-pow(y,2) / (2*pow(sigma, 2))) / exp(estimate);
				double nominator2 = (exp(-pow(y+2*sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				double denominator1 = exp(-pow(y-2*sqrt(Esender), 2) / (2 * pow(sigma, 2))) / exp(estimate);
				double denominator2 = (exp(-pow(y, 2) / (2 * pow(sigma, 2))));

				if(isinf(nominator1) && isinf(denominator1)) {
					codedChannelOutput[i] = log(exp(-pow(y, 2) / (2 * pow(sigma, 2))) / exp(-pow(y - 2 * sqrt(Esender), 2) / (2 * pow(sigma, 2))));
				} else {
					codedChannelOutput[i] = log((nominator1 + nominator2) / (denominator1 + denominator2));
				}

			}
			
			// Update side information.
			for(int i=0; i<sysNumber; ++i) {
				double decision = sysNodeInfo2[i].llrEstimation;

				double x = exp(decision);

				if(isinf(x))
					sideInfo[i] = log((1 - p) / p);
				else {
					sideInfo[i] = log(((1 - p) * x + p) / (1 - p + p * x));
				}
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
			int diffPosInSys = 0;
      for(int i=0; i<sysNumber; i++) {
				if(sysNodeInfo1[i].decision != sysNodeInfo2[i].decision) {
	  			diffPos++;
					if(systrans1[i] != systrans2[i])
						diffPosInSys++;
				}
      }
      cout << "Number of different decisions: " << diffPos << " " << diffPosInSys << endl;

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

vector<int> interleaver(vector<double> &source, double percentage, vector<Sys_info> &sysNumberde, vector<Coded_info> &codedNodeP, int codedNoP, int group)
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

      //first, find the sysNumberde Num that connects to the original coded Node
      for(int i=0;i<degree;++i)
	{
	  int sysNumberdeNum=codedNodeP[oldNum].nodeData.neighbourNum[i];
	  int coded_index=codedNodeP[oldNum].nodeData.inmessageIndex[i];
	  
	  if(group==1)
	    sysNumberde[sysNumberdeNum].nodeData.neighbourNum[coded_index]=newNum;
	  if(group==2)
	    sysNumberde[sysNumberdeNum].nodeData.para_neighbourNum[coded_index]=newNum;
	}

      interCoded[newNum]=codedNodeP[oldNum];
    }

  for(int index=0;index<interleavingNum;++index)
    codedNodeP[index]=interCoded[index];
 
  return newPos;
}

void deinterleaver(double percentage, int *newPos, vector<Sys_info> &sysNumberde, vector<Coded_info> &codedNodeP, int codedNoP, int group)
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
	    sysNumberde[sysNum].nodeData.neighbourNum[codedIndex]=originNum;
	  if(group==2)
	    sysNumberde[sysNum].nodeData.para_neighbourNum[codedIndex]=originNum;
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
