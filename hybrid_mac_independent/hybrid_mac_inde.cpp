#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <algorithm>

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

  /*
  //For weightset size of 8
  int horizontalseq[][4] = { {2, 3, 1, 0},
                             {1, 0, 2, 3},
                             {3, 2, 0, 1},
                             {0, 1, 3, 2} };
  */
  
  //For weightset size of 16
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
int countError(vector<Sys_info> &sysNodeInfo, vector<int> &sysTrans);

int main() 
{
  // Source info
  int sysNumber = 10000;
  int symbolNumber = 10000;
  int codedNumber = 10000;
  int sysDegree = 30;
  
  
  double p1 = 0.5;    //sparsity of the source.
  double p0 = 1 - p1;
  double p = 0.01; // correlation between source 1 and source 2.

  double H_U1_givenU2 = -p*log2(p) - (1 - p) * log2(1 - p);
  double H_U2 = 1;
  double entropy = H_U1_givenU2 + H_U2;

  double codeRate = double(sysNumber) / (sysNumber + symbolNumber + codedNumber);
  double R = (H_U1_givenU2 + H_U2) * codeRate * 2;
  double throughput = codeRate * 2;
  double capacity = 10 * log10((pow(2, R) - 1) / (4 * codeRate));   // in terms of Eso/N0

  cout<<"Code rate is: "<<codeRate<<endl;
  cout<<"entropy is: "<<entropy<<endl;
  cout<<"Throughput is:"<<throughput<<endl;
  cout<<"channel capacity is: "<<capacity<<endl;


  //weightset
  int wval[] = {1, 1, 1, 1, 2, 2, 2, 2};
  //int wval[] = {2, 3, 4, 8};
  
  int weightsize = 2 * sizeof(wval) / sizeof(int);
  vector<double> weightset(weightsize, 0);
  
  int j = 0;
  for(int i=0;i<weightsize;i=i+2) {
      weightset[i] = wval[j];
      weightset[i+1] = -wval[j];

      ++j;
    }

  vector<int> weightval(wval, wval+weightsize/2);
 

  // Get the general information about RP symbol
  map<int, double> symbolsetdata = getSymbolSet(weightset, p1);
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

  // Two sources use the same RCM encoder while different LDGM encoders.
  // Determine how many G0 are needed to build G
   int G0RowNum = sysNumber / (weightsize / 2) / 2 * 4;
   int G0Num = 1;
   while(G0RowNum * G0Num < symbolNumber)
     ++G0Num;
  
   double gmatrixSeed = 200;
  Matrix gmatrix1 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
  for(int i=0;i<G0Num;++i) {
    Matrix gmatrix0 = linearMatrixConsD8(sysNumber, weightval, gmatrixSeed++);
    gmatrix1.matrixStackDown(gmatrix0);
  }

  // Generate two different LDGM generator matrices.
  vector<multimap<int, int> > maData1 = paritycm(sysNumber, codedNumber, sysDegree, sysNumber, 2650);
  multimap<int, int> onePositionByRow1(maData1[0]);
  multimap<int, int> onePositionByColumn1(maData1[1]);

  vector<multimap<int, int> > maData2 = paritycm(sysNumber, codedNumber, sysDegree, sysNumber, 3350);
  multimap<int, int> onePositionByRow2(maData2[0]);
  multimap<int, int> onePositionByColumn2(maData2[1]);

  // Scourc bit nodes
  vector<Sys_info> sysNodeInfo1(sysNumber);
  vector<Sys_info> sysNodeInfo2(sysNumber);
  
  // RP symbol nodes
  vector<Coded_info> symNodeInfo1(symbolNumber);
  vector<Coded_info> symNodeInfo2(symbolNumber);

  // Coded bit nodes
  vector<Coded_info> codedNodeInfo1(codedNumber);
  vector<Coded_info> codedNodeInfo2(codedNumber);

  // Read RCM and LDGM generator matrix.
  gmatrix1.organizebyCol();
  for(int index=0;index<sysNumber;++index) {
    sysNodeInfo1[index].readLinearMatrix(gmatrix1, index, symbolNumber, 1);
    sysNodeInfo2[index].readLinearMatrix(gmatrix1, index, symbolNumber, 1);
    sysNodeInfo1[index].readNeighbourNum(onePositionByRow1, index, 2);
    sysNodeInfo2[index].readNeighbourNum(onePositionByRow2, index, 2);
  }
  
  for(int index=0;index<symbolNumber;++index) {
    symNodeInfo1[index].readLinearMatrix(gmatrix1, index);
    symNodeInfo2[index].readLinearMatrix(gmatrix1, index);
  }
 
  for(int index=0;index<codedNumber;++index) {
    codedNodeInfo1[index].readNeighbourNum(onePositionByColumn1, index);
    codedNodeInfo2[index].readNeighbourNum(onePositionByColumn2, index);
  }

  // Build graph.
  for(int index=0; index<sysNumber; ++index) {
    sysNodeInfo1[index].getsourceNum(symNodeInfo1, codedNodeInfo1, index);
    sysNodeInfo2[index].getsourceNum(symNodeInfo2, codedNodeInfo2, index);
  }
  for(int index=0;index<symbolNumber;++index) {
    symNodeInfo1[index].getsourceNum(sysNodeInfo1, index, 1);
    symNodeInfo2[index].getsourceNum(sysNodeInfo2, index, 1);
  }
  for(int index=0;index<codedNumber;++index) {
    codedNodeInfo1[index].getsourceNum(sysNodeInfo1, index, 2);
    codedNodeInfo2[index].getsourceNum(sysNodeInfo2, index, 2);
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
  vector<double> rpSymbols1(symbolNumber, 0);
  vector<double> rpSymbols2(symbolNumber, 0);
  vector<int> codedBits1(codedNumber, 0);
  vector<int> codedBits2(codedNumber, 0);

  // Our assumption is that the generated symbols are symmetric, i.e
  // the number of possible symbols is odd and symmetric to zero.
  // The length of the vector should be determined by the largest value
  // of the generated symbols.
  int sympdfsize = 0;
  // The following is under the contidition that the generated symbols are balanced
  sympdfsize = symbolset[symbolset.size()-1] * 2 + 1;   

  vector< vector<double> > channelPdf(symbolNumber, vector<double>(sympdfsize, 0));    // PDF message from channel(symbols)
  vector<double> sysChannelOutput(sysNumber, 0);                       // LLR message from channel(sysbits)
  vector<double> codedChannelOutput(codedNumber, 0);                   // LLR message from channel(codedBits)


  // Channel condition.
  double gap = 5;
  double snr = capacity + gap;
  double sigma = sqrt( Eso / (2 * pow(10, snr/10)) ); 
  double noiseVar = pow(sigma, 2);

  cout << "Gap to the limit is: " << gap << endl;
  cout << "SNR is: " << snr << endl;

  const int blocknum = 5000;
  int currentblock = 0;

  double bersum = 0;
  
  int seed = 2000;

  // Error file.
  ofstream error;
  stringstream stream1;
  stream1 << fixed << setprecision(2) << p;
  string pString = stream1.str();

  stringstream stream2;
  stream2 << fixed << setprecision(2) << gap;
  string gapString = stream2.str();
  string errorFileName = "mac_hybrid_inde_w4142_" + std::to_string(static_cast<long long>(symbolNumber)) + "_" + std::to_string(static_cast<long long>(codedNumber)) + "_" + std::to_string(static_cast<long long>(sysDegree)) + "_p" + pString + "_" + gapString + "db.txt";
	
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

  int geneThreshold = p * 1000000;
  while(currentblock < blocknum) {

    ++currentblock;

    error<<"Block "<<currentblock<<":";
    
    srand(seed++);

    // Generate souce bits; Use identification bit
    systematicBits1 = sourceGenerator(sysNumber, p1, seed++);
    systematicBits1[0] = 1;

    for(int i=0; i<systematicBits1.size(); i++) {
      int random = rand() % 1000000;
      if(random < geneThreshold)
        systematicBits2[i] = 1 ^ systematicBits1[i];
      else
        systematicBits2[i] = systematicBits1[i];
    }
    systematicBits2[0] = 0;


    vector<int> sysTrans1(systematicBits1);
    vector<int> sysTrans2(systematicBits2);

    // RCM Encoder
    linearEncoder(symNodeInfo1, systematicBits1, rpSymbols1);
    linearEncoder(symNodeInfo2, systematicBits2, rpSymbols2);

    // LDGM encoder
    encoder(codedNodeInfo1, systematicBits1, codedBits1);     
    encoder(codedNodeInfo2, systematicBits2, codedBits2);
    
    // Generate channel symbol for systematic bits.
    vector<complex<double> > channelSignals;
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
    vector<complex<double> > tempChannelSignals1 = modulatorQAM(rpSymbols1, norfactor);
    vector<complex<double> > tempChannelSignals2 = modulatorQAM(rpSymbols2, norfactor);
    for(int i=0; i<tempChannelSignals1.size(); i++)
      channelSignals.push_back(tempChannelSignals1[i] + tempChannelSignals2[i]);
    
    for(int i=0; i<codedNumber; i=i+2) {
      double tempAX = (2 * codedBits1[i] - 1) / sqrt(2) * sqrt(Es1);
      double tempAY = (2 * codedBits1[i+1] - 1) / sqrt(2) * sqrt(Es1);
      complex<double> digitSignalA = complex<double>(tempAX, tempAY);

      double tempBX = (2 * codedBits2[i] - 1) / sqrt(2) * sqrt(Es2);
      double tempBY = (2 * codedBits2[i+1] - 1) / sqrt(2) * sqrt(Es2);
      complex<double> digitSignalB = complex<double>(tempBX, tempBY);
      channelSignals.push_back(digitSignalA + digitSignalB);
    }
    
    // Send symbols through complex AWGN channel
    complexAWGNChannel(channelSignals, sigma);                 

    // Calculate the pdf of corrupted RP symbols (Initialization for symbol nodes)    
    vector<double> x2Prior(sympdfsize, 0);
    for(int index=0;index<symbolNumber;++index) {
      double y;
      if(index % 2 == 0)
        y = channelSignals[(sysNumber + index) / 2].real();
      else
        y = channelSignals[(sysNumber + index) / 2].imag();

      vector<double> symdis(sympdfsize);
      for(int symindex=0; symindex<symbolset.size(); ++symindex) {
        int zepoint = (sympdfsize - 1) / 2;
        double symVal = symbolset[symindex];
        int symPos = symVal + zepoint;
        symdis[symPos] =  getSymPdf(symbolset, x2Prior, symVal, norfactor, y, sigma); 
      }
      vecNorm(symdis);
      channelPdf[index] = symdis;
      symNodeInfo1[index].clearMessage();
      symNodeInfo2[index].clearMessage();
    }

    // Initialization for coded nodes
    for(int index=0;index<codedNumber;++index) {
      double y;
      if(index % 2 == 0)
        y = channelSignals[(sysNumber + symbolNumber + index)/2].real();
      else
        y = channelSignals[(sysNumber + symbolNumber + index)/2].imag();

      codedChannelOutput[index] = log((exp(-pow(y + 2 /sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)) + exp(-pow(y, 2) / (2 * noiseVar))) / (exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2)/ (2 * noiseVar)) + exp(-pow(y, 2) / (2 * noiseVar))));                
        codedNodeInfo1[index].clearMessage();
        codedNodeInfo2[index].clearMessage();
    }
        
    // Initialization for sys nodes
    for(int index=0;index<sysNumber;++index) {
      double y;
      if(index % 2 == 0)
        y = channelSignals[index/2].real();
      else
        y = channelSignals[index/2].imag();

      sysChannelOutput[index] = log((exp(-pow(y + 2 /sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)) + exp(-pow(y, 2) / (2 * noiseVar))) / (exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2)/(2 * noiseVar)) + exp(-pow(y, 2) / (2 * noiseVar))));    

      sysNodeInfo1[index].clearMessage();
      sysNodeInfo2[index].clearMessage();
    }
    
    // side information
    vector<double> sideInfo(sysNumber, 0);
    double csiThreshold = 20;
 
    int iterationtime = 200;
    int currenttime =0;  
    double watchllr[sysNumber];

    vector<int> predecision1(sysNumber, 0);
    vector<int> predecision2(sysNumber, 0);
    int sametime = 0;
    while(currenttime != iterationtime+1) {

      ++currenttime;
      
      // Decoder 1
      for(int i=0;i<symbolNumber;++i)
        symNodeInfo1[i].computePdf(sysNodeInfo1, channelPdf[i], norfactor, noiseVar, 1, i);
     
      for(int i=0;i<codedNumber;++i)
        codedNodeInfo1[i].computeMessage(sysNodeInfo1, codedChannelOutput[i], 2);

      for(int i=0;i<sysNumber;++i)
        sysNodeInfo1[i].computeMessage(symNodeInfo1, codedNodeInfo1, sysChannelOutput[i], sideInfo[i]);
   
      // Channel state node
      for(int index=0; index<symbolNumber; index++) {
        double y;
        if(index % 2 == 0)
          y = channelSignals[(sysNumber + index) / 2].real();
        else
          y = channelSignals[(sysNumber + index) / 2].imag();

        vector<double> estimate = symNodeInfo1[index].pdf_to_state;
        vector<double> symdis(sympdfsize);
        for(int symindex=0;symindex<symbolset.size(); symindex++) {
          int zepoint = (sympdfsize - 1) / 2;
          double symVal = symbolset[symindex];
          int symPos = symVal + zepoint;
          symdis[symPos] =  getSymPdf(symbolset, estimate, symVal, norfactor, y, sigma); 
        }
        vecNorm(symdis);
        channelPdf[index] = symdis;
      }
       
      for(int index=0;index<codedNumber; index++) {
        double y;
        if(index % 2 == 0)
          y = channelSignals[(sysNumber + symbolNumber + index)/2].real();
        else
          y = channelSignals[(sysNumber + symbolNumber + index)/2].imag();

        double estimate = codedNodeInfo1[index].message_to_state;
        double nominator_ele1 = exp(-pow(y, 2) / (2 * noiseVar)) / exp(estimate);
        double nominator_ele2 = exp(-pow(y + 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar));
        double denominator_ele1 = exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)) / exp(estimate);
        double denominator_ele2 = (exp(-pow(y,2)/(2 * noiseVar)));

        if( std::isinf(nominator_ele1) && std::isinf (denominator_ele1) ) 
          codedChannelOutput[index] = log( exp(-pow(y,2)/(2 * noiseVar)) / exp(-pow(y-2 /sqrt(2) * sqrt(Esender), 2)/(2 * noiseVar)));
        else
          codedChannelOutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

        if( std::isnan(codedChannelOutput[index]) || std::isinf(codedChannelOutput[index]) ) 
            cout<<"sys channel message not a number!"<<endl;
      }

      for(int index=0; index<sysNumber; index++) {
        double y;
        if(index % 2 == 0)
          y = channelSignals[index / 2].real();
        else
          y = channelSignals[index / 2].imag();

        double estimate = sysNodeInfo1[index].message_to_state;
        double nominator_ele1 = exp(-pow(y, 2) / (2 * noiseVar)) / exp(estimate);
        double nominator_ele2 = exp(-pow(y + 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar));
        double denominator_ele1 = exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)) / exp(estimate);
        double denominator_ele2 = (exp(-pow(y,2)/(2 * noiseVar)));

        if( std::isinf(nominator_ele1) && std::isinf(denominator_ele1) ) 
          sysChannelOutput[index] = log(exp(-pow(y, 2)/(2 * noiseVar)) / exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)));
        else
          sysChannelOutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

        if( std::isnan(sysChannelOutput[index]) || std::isinf(sysChannelOutput[index]) ) 
            cout<<"sys channel message not a number!"<<endl;
      }
      
      // Correlation link
      for(int index=0; index<sysNumber; index++) {
        double x = exp(sysNodeInfo1[index].llrEstimation);
        
        if( std::isinf(x) )
          x = log( (1-p)/p );
        else
          x = log((1 + (1 - p) / p * x) / (x + (1 - p) / p));

        if(abs(x) >= csiThreshold)
          sideInfo[index] = x > 0 ? csiThreshold : -csiThreshold;
        
        if( std::isnan(sideInfo[index]) || std::isinf(sideInfo[index]) )
            cout<<"side info2 wrong!"<<endl;    
      }
    

      // Decoder 2
      for(int i=0;i<symbolNumber;++i)
        symNodeInfo2[i].computePdf(sysNodeInfo2, channelPdf[i], norfactor, noiseVar, 1, i);
     
      for(int i=0;i<codedNumber;++i)
        codedNodeInfo2[i].computeMessage(sysNodeInfo2, codedChannelOutput[i], 2);

      for(int i=0;i<sysNumber;++i)
        sysNodeInfo2[i].computeMessage(symNodeInfo2, codedNodeInfo2, sysChannelOutput[i], sideInfo[i]);
      
      // Channel nodes
      for(int index=0; index<symbolNumber; index++) {
        double y;
        if(index % 2 == 0)
          y = channelSignals[(sysNumber + index) / 2].real();
        else
          y = channelSignals[(sysNumber + index) / 2].imag();

        vector<double> estimate = symNodeInfo2[index].pdf_to_state;
        vector<double> symdis(sympdfsize);
        for(int symindex=0;symindex<symbolset.size(); symindex++) {
          int zepoint = (sympdfsize - 1) / 2;
          double symVal = symbolset[symindex];
          int symPos = symVal + zepoint;
          symdis[symPos] =  getSymPdf(symbolset, estimate, symVal, norfactor, y, sigma); 
        }
        vecNorm(symdis);
        channelPdf[index] = symdis;
      }
       
      for(int index=0;index<codedNumber; index++) {
        double y;
        if(index % 2 == 0)
          y = channelSignals[(sysNumber + symbolNumber + index)/2].real();
        else
          y = channelSignals[(sysNumber + symbolNumber + index)/2].imag();

        double estimate = codedNodeInfo2[index].message_to_state;
        double nominator_ele1 = exp(-pow(y, 2) / (2 * noiseVar)) / exp(estimate);
        double nominator_ele2 = exp(-pow(y + 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar));
        double denominator_ele1 = exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)) / exp(estimate);
        double denominator_ele2 = (exp(-pow(y,2)/(2 * noiseVar)));

        if( std::isinf(nominator_ele1) && std::isinf (denominator_ele1) ) 
          codedChannelOutput[index] = log( exp(-pow(y,2)/(2 * noiseVar)) / exp(-pow(y-2 /sqrt(2) * sqrt(Esender), 2)/(2 * noiseVar)));
        else
          codedChannelOutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

        if( std::isnan(codedChannelOutput[index]) || std::isinf(codedChannelOutput[index]) ) 
            cout<<"sys channel message not a number!"<<endl;
      }

      for(int index=0; index<sysNumber; index++) {
        double y;
        if(index % 2 == 0)
          y = channelSignals[index / 2].real();
        else
          y = channelSignals[index / 2].imag();

        double estimate = sysNodeInfo2[index].message_to_state;
        double nominator_ele1 = exp(-pow(y, 2) / (2 * noiseVar)) / exp(estimate);
        double nominator_ele2 = exp(-pow(y + 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar));
        double denominator_ele1 = exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)) / exp(estimate);
        double denominator_ele2 = (exp(-pow(y,2)/(2 * noiseVar)));

        if( std::isinf(nominator_ele1) && std::isinf(denominator_ele1) ) 
          sysChannelOutput[index] = log(exp(-pow(y, 2)/(2 * noiseVar)) / exp(-pow(y - 2 / sqrt(2) * sqrt(Esender), 2) / (2 * noiseVar)));
        else
          sysChannelOutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

        if( std::isnan(sysChannelOutput[index]) || std::isinf(sysChannelOutput[index]) ) 
            cout<<"sys channel message not a number!"<<endl;
      }

      // Correlation link
      for(int index=0; index<sysNumber; index++) {
        double x = exp(sysNodeInfo2[index].llrEstimation);
        
        if( std::isinf(x) )
          x = log( (1-p)/p );
        else
          x = log((1 + (1 - p) / p * x) / (x + (1 - p) / p));

        if(abs(x) >= csiThreshold)
          sideInfo[index] = x > 0 ? csiThreshold : -csiThreshold;
        
        if( std::isnan(sideInfo[index]) || std::isinf(sideInfo[index]) )
            cout<<"side info2 wrong!"<<endl;    
      }

      // Count error number.
      double errorbit = 0;
      
      if(sysNodeInfo1[0].decision == 0) {
        errorbit += countError(sysNodeInfo1, sysTrans1);
        errorbit += countError(sysNodeInfo2, sysTrans2);
      } else if(sysNodeInfo1[0].decision == 1) {
        errorbit += countError(sysNodeInfo1, sysTrans2);
        errorbit += countError(sysNodeInfo2, sysTrans1);
      } 
           
      cout<<"error number is:"<<errorbit<<"("<<currenttime<<")"<<endl;

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
    
    if(sysNodeInfo1[0].decision == 0) {
      errorbit += countError(sysNodeInfo1, sysTrans1);
      errorbit += countError(sysNodeInfo2, sysTrans2);
    } else if(sysNodeInfo1[0].decision == 1) {
      errorbit += countError(sysNodeInfo1, sysTrans2);
      errorbit += countError(sysNodeInfo2, sysTrans1);
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

    cout<<" BER for SNR = "<<snr<<" db is "<<blockber<<"  sum:"<<bersum<<"("<<currentblock<<")  "<<endl;
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

double getSymPdf(vector<double> symbolSet, vector<double> estimate, double candiVal, double norFactor, double y, double sigma) {

  // y' = y - norfactor * x1.
  double yPrime = y - norFactor * candiVal;
  double resPdf = 0;

  // P{y | x1 = k} = sum_x2 P(x2) * P{x2 + n = y'}
  // P(x2) is the probability of x2 being the value.
  double zeroPos = (estimate.size() - 1) / 2;
  for(int i=0; i<symbolSet.size(); i++) {
    double x2Pro = estimate[symbolSet[i] + zeroPos];
    
    // If there is no estimate, uniform distribution is assumed.
    if(x2Pro == 0)
      x2Pro = double(1) / symbolSet.size();
    double x2 = symbolSet[i] * norFactor;
    resPdf += x2Pro * gaussianFunc(yPrime, x2, sigma);
  }

  return resPdf; 
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
          int sysNumberdeNum=codedNodeP[oldNum].nodeData.neighbourNum[i];
          int coded_index=codedNodeP[oldNum].nodeData.inmessageIndex[i];
          
          if(group==1)
            sysNode[sysNumberdeNum].nodeData.neighbourNum[coded_index]=newNum;
          if(group==2)
            sysNode[sysNumberdeNum].nodeData.para_neighbourNum[coded_index]=newNum;
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

int countError(vector<Sys_info> &sysNodeInfo, vector<int> &sysTrans) {

  int errorNum = 0;
  for(int i=0; i<sysTrans.size(); i++) {
    if(sysNodeInfo[i].decision != sysTrans[i])
      errorNum++;
  }
  
  return errorNum;
}
