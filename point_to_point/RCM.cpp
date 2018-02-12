#include<iostream>
#include<vector>
#include<map>
#include<stdlib.h>
#include<cmath>
#include<fstream>
#include<algorithm>
#include"node.h"
#include"sparsematrix.h"
#include"func.h"
#include"pcm.cpp"

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
			     {2, 3, 0, 1},
			     {0, 1, 2, 3} };
  
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
  int sys_no = 10000;
  int symbol_no = 20000;
  
  double p1 = 0.5;    //sparsity of the source.
  double p0 = 1 - p1;

  double entropy = -p0*log2(p0) - p1*log2(p1);
  double coderate = double(sys_no) / (symbol_no);
  double inforate = entropy * coderate * 2;
  double throughput = coderate * 2;
  double capacity = 10 * log10( pow(2, inforate) - 1);     // Channel capacity in terms of Es/N0.

  cout<<"Code rate is: "<<coderate<<endl;
  cout<<"entropy is: "<<entropy<<endl;
  cout<<"Throughput is:"<<throughput<<endl;
  cout<<"channel capacity is: "<<capacity<<endl;


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
 
  map<double, double> symbolsetdata = getSymbolSet(weightset, p1);
  vector<double> symbolset;
  vector<double> symbolpro;
  
  // Get symbol set and its distribution
   map<double, double>::iterator it = symbolsetdata.begin();

   for(it; it!=symbolsetdata.end(); ++it)
    {
      symbolset.push_back( it->first ); 
      symbolpro.push_back( it->second );
    }
  
   
   double symbol_entropy = 0;
   for(int i=0;i<symbolpro.size();++i)
     symbol_entropy += (-symbolpro[i] * log2(symbolpro[i]));

   cout << "Symbol entropy:" << symbol_entropy << endl;

   // Generator matrix for symbols
   // Determine how many G0 are needed to build G
   int G0RowNum = sys_no / (weightsize / 2) / 2 * 4;
   int G0Num = 1;
   while(G0RowNum * G0Num < symbol_no)
     ++G0Num;
  Matrix gmatrix = linearMatrixConsD8(sys_no, weightval, 100);
  double gmatrixSeed = 1000;
  for(int i=0;i<G0Num;++i)
    {
      Matrix gmatrix1 = linearMatrixConsD8(sys_no, weightval, gmatrixSeed++);
      gmatrix.matrixStackDown(gmatrix1);
    }
  vector<Sys_info> sysnodeinfo(sys_no);   // Establish the data for systematic bits
  vector<Coded_info> symnodeinfo(symbol_no);         // Establish the data for RP symbols

  gmatrix.organizebyCol();
  for(int index=0;index<sys_no;++index)
    sysnodeinfo[index].readLinearMatrix(gmatrix, index, symbol_no, 1);
  for(int index=0;index<symbol_no;++index)
    symnodeinfo[index].readLinearMatrix(gmatrix, index);

 
  
  
  // Get each neighbor's index in storage
  for(int index=0;index<symbol_no;++index)                 
    symnodeinfo[index].getsourceNum(sysnodeinfo, index, 1);
  for(int index=0;index<sys_no;++index)
    sysnodeinfo[index].getsourceNum(symnodeinfo, index);

 
  // Calculate the average value of transimitted symbols
  double E_ave = 1;
  double Eb = E_ave / inforate;
  cout<<"average energy is:"<<E_ave<<endl;

 
  // Start to generate sys bits and transmitted symbols
  vector<int> systematicbits(sys_no,0);
  vector<double> rpsymbols(symbol_no, 0);
  

  // Our assumption is that the generated symbols are symmetric, i.e
  // the number of possible symbols is odd and symmetric to zero.
  // The length of the vector should be determined by the largest value
  // of the generated symbols.
  int sympdfsize = 0;
  // The following is under the contidition that the generated symbols are balanced
  sympdfsize = symbolset[symbolset.size()-1] * 2 + 1;   

  vector<double> pdfinitializer(sympdfsize,0);
  vector< vector<double> > channelpdf(symbol_no, pdfinitializer);    // PDF message from channel(symbols)
  vector<double> sys_channeloutput(sys_no, 0);                       // LLR message from channel(sysbits)
  

  double snr = 9;
  double sigma = sqrt( E_ave / (2 * pow(10, snr/10)) ); 
 
  const int blocknum = 1000;
  int currentblock = 0;

  double bersum = 0;
  
  int seed = 2000;

  ofstream error;
  error.open("res_ana_23710_Rc05_9db_2000.txt", ios::out);
  //error.open("test.txt", ios::out);

  ofstream erroreachblock;
  //erroreachblock.open("err_m4_2500_5_001_t10_2p5db.txt", ios::out);
  erroreachblock.open("test1.txt", ios::out);

  vector<double> error_record;
  double error_min = sys_no;
  double error_max = 0;
  double error_ave;

  double norfactor = getNormalizationFactorQAM(symbolset, symbolpro);
  vector<double> normalized_symbolset(symbolset);
    for(int i=0;i<symbolset.size();++i)
      normalized_symbolset[i] *= norfactor;

  while(currentblock < blocknum) {

    ++currentblock;

    error<<"Block "<<currentblock<<":";
    
    srand(seed++);

    // Generate souce bits
    systematicbits = sourceGenerator(sys_no, p1, seed);
    
    // Get RP symbols
    linearEncoder(symnodeinfo, systematicbits, rpsymbols);       // Symbol Encoder 
    
    // Modulate RP symbols
    vector<ComplexNumber> channelsignals = modulatorQAM(rpsymbols, norfactor);   // QAM and get the normalized factor.
    
    // Complex AWGN channel
    complexAWGNChannel(channelsignals, sigma);                 

  
    // Calculate the pdf of corrupted RP symbols (Initialization for symbol nodes)    
    for(int index=0;index<symbol_no;++index)
      {
	double y;
	if(index % 2 == 0)
	  y = channelsignals[index/2].real;
	else
	  y = channelsignals[index/2].imaginary;

	vector<double> symdis(sympdfsize);
      
	for(int symindex=0;symindex<symbolset.size();++symindex)
	  {
	    int zepoint = (sympdfsize - 1) / 2;
	    double symva = symbolset[symindex];
	    int sympos = symva + zepoint;
	    symdis[sympos] =  gaussianFunc(y, normalized_symbolset[symindex], sigma);
	  }
	vecNorm(symdis);
	channelpdf[index] = symdis;
	symnodeinfo[index].clearMessage();
      }
	
    // Initialization for sys nodes
    for(int index=0;index<sys_no;++index)
      {
	double initialmessage = log(p0/p1);
	sys_channeloutput[index] = initialmessage; 
	sysnodeinfo[index].clearMessage();
      }

 
    // Start the message passing
    int iterationtime = 200;
    int currenttime =0;  
    double watchllr[sys_no];

    vector<int> predecision(sys_no, 0);
    int sametime = 0;

    while(currenttime != iterationtime+1) {

      ++currenttime;
    
      for(int i=0;i<symbol_no;++i)
	symnodeinfo[i].computePdf(sysnodeinfo, channelpdf[i], 1);

      for(int i=0;i<sys_no;++i)
	sysnodeinfo[i].computeMessage(symnodeinfo, sys_channeloutput[i]);

      double errorbit = 0;
      for(int i=0;i<sys_no;++i)
	{
	  if(sysnodeinfo[i].decision != systematicbits[i])
	    ++errorbit;
	}   

      cout<<"error number is:"<<errorbit<<"("<<currenttime<<")"<<endl;

      error<<errorbit<<" ";      

      if(sameDecision(predecision, sysnodeinfo))
	++sametime;
      else
	sametime=0;

      if(sametime == 5) {
	cout<<" Make same decision for 5 consecutive iterations !"<<endl;
	break;
      }

    }
    
    double errorbit = 0;

    for(int i=0;i<sys_no;++i)
      {
	if(sysnodeinfo[i].decision != systematicbits[i])
	  ++errorbit;
      } 
  
    double blockber = errorbit / sys_no;
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




