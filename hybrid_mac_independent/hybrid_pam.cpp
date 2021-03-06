#include<iostream>
#include<vector>
#include<map>
#include<stdlib.h>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<omp.h>
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

double getSymPdf(vector<double> symbolSet, vector<double> estimate, double candiVal, double norFactor, double y, double sigma);
int countError(vector<Sys_info> &sysnodeinfo, vector<int> &systrans);

int main() 
{
  // Source info
  int sys_no = 10000;
  int symbol_no = 10000;
  int coded_no = 10000;
  int para_sysdegree = 10;
  
  
  double p1 = 0.5;    //sparsity of the source.
  double p0 = 1 - p1;
  double p = 0.01; // correlation between source 1 and source 2.

  double entropy = -p0*log2(p0) - p1*log2(p1);
  double coderate = double(sys_no) / (symbol_no + coded_no);
  double inforate = entropy * coderate * 2;
  double R = inforate + inforate;
  double throughput = coderate * 2;
  /********/
  double capacity = 10 * log10( pow(2, inforate) - 1);     // Channel capacity in terms of Es/N0.

  cout<<"Code rate is: "<<coderate<<endl;
  cout<<"entropy is: "<<entropy<<endl;
  cout<<"Throughput is:"<<throughput<<endl;
  cout<<"channel capacity is: "<<capacity<<endl;


  //weightset
  //int wval[] = {1, 1, 1, 1, 2, 2, 2, 2};
  int wval[] = {2, 3, 4, 8};
  
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
 

  // Get the general information about RP symbol
  map<double, double> symbolsetdata = getSymbolSet(weightset, p1);
  map<int, vector<vector<int> > > mappingTable = getMappingTable(weightset);
  vector<double> symbolset;
  vector<double> symbolpro;
  
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

   // Generator RCM matrix for two RCM encoders
   // Determine how many G0 are needed to build G
   int G0RowNum = sys_no / (weightsize / 2) / 2 * 4;
   int G0Num = 1;
   while(G0RowNum * G0Num < symbol_no)
     ++G0Num;
  
   double gmatrixSeed = 200;
  Matrix gmatrix1 = linearMatrixConsD8(sys_no, weightval, gmatrixSeed++);
  for(int i=0;i<G0Num;++i)
    {
      Matrix gmatrix0 = linearMatrixConsD8(sys_no, weightval, gmatrixSeed++);
      gmatrix1.matrixStackDown(gmatrix0);
    }

  gmatrixSeed = 200;
  Matrix gmatrix2 = linearMatrixConsD8(sys_no, weightval, gmatrixSeed++);
  for(int i=0; i<G0Num; i++) {
    Matrix gmatrix0 = linearMatrixConsD8(sys_no, weightval, gmatrixSeed++);
    gmatrix2.matrixStackDown(gmatrix0);
  }


  vector<Sys_info> sysnodeinfo1(sys_no);   // Establish the data for systematic bits
  vector<Sys_info> sysnodeinfo2(sys_no);
  vector<Coded_info> symnodeinfo1(symbol_no);         // Establish the data for RP symbols
  vector<Coded_info> symnodeinfo2(symbol_no);

  gmatrix1.organizebyCol();
  gmatrix2.organizebyCol();
  for(int index=0;index<sys_no;++index) {
    sysnodeinfo1[index].readLinearMatrix(gmatrix1, index, symbol_no, 1);
    sysnodeinfo2[index].readLinearMatrix(gmatrix2, index, symbol_no, 1);
  }
  for(int index=0;index<symbol_no;++index) {
    symnodeinfo1[index].readLinearMatrix(gmatrix1, index);
    symnodeinfo2[index].readLinearMatrix(gmatrix2, index);
  }
 
  // Generator matrix for digit bits
  vector<multimap<int, int> > madata = paritycm(sys_no, coded_no, para_sysdegree, sys_no, 1500);
  multimap<int, int> onePositionbyRow(madata[0]);
  multimap<int, int> onePositionbyColumn(madata[1]);


  vector<Coded_info> codednodeinfo1(coded_no);
  vector<Coded_info> codednodeinfo2(coded_no);

  // This is for regular degree data reading    
  for(int index=0;index<sys_no;++index) {
    sysnodeinfo1[index].readNeighbourNum(onePositionbyRow, index, 2);
    sysnodeinfo2[index].readNeighbourNum(onePositionbyRow, index, 2);
  }
  for(int index=0;index<coded_no;++index) {
    codednodeinfo1[index].readNeighbourNum(onePositionbyColumn, index);
    codednodeinfo2[index].readNeighbourNum(onePositionbyColumn, index);
  }

  for(int index=0;index<symbol_no;++index) {                // Get its index in its neighbour's data.
    symnodeinfo1[index].getsourceNum(sysnodeinfo1, index, 1);
    symnodeinfo2[index].getsourceNum(sysnodeinfo2, index, 1);
  }
  for(int index=0;index<coded_no;++index) {
    codednodeinfo1[index].getsourceNum(sysnodeinfo1, index, 2);
    codednodeinfo2[index].getsourceNum(sysnodeinfo2, index, 2);
  }
  for(int index=0;index<sys_no;++index) {
    sysnodeinfo1[index].getsourceNum(symnodeinfo1, codednodeinfo1, index);
    sysnodeinfo2[index].getsourceNum(symnodeinfo2, codednodeinfo2, index);
  } 

 

  // Calculate the average value of transimitted symbols
  double Es = 2;
  double Es1 = Es / 2;
  double Es2 = Es / 2;

  double Eso = Es / coderate;

  cout << "energy per source bit is:" << Eso << endl;

 
  // Start to generate sys bits and transmitted symbols
  vector<int> systematicbits1(sys_no, 0);
  vector<int> systematicbits2(sys_no, 0);
  vector<double> rpsymbols1(symbol_no, 0);
  vector<double> rpsymbols2(symbol_no, 0);
  vector<int> codedbits1(coded_no, 0);
  vector<int> codedbits2(coded_no, 0);

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
  vector<double> coded_channeloutput(coded_no, 0);                   // LLR message from channel(codedbits)

  double snr = 10;
  double sigma = sqrt( Eso / (2 * pow(10, snr/10)) ); 
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
  double error_min = sys_no;
  double error_max = 0;
  double error_ave;

  double norfactor = getNormalizationFactorPAM(symbolset, symbolpro);

  int geneThreshold = p * 1000000;
  while(currentblock < blocknum) {

    ++currentblock;

    error<<"Block "<<currentblock<<":";
    
    srand(seed++);

    // Generate souce bits; Use identification bit
    systematicbits1 = sourceGenerator(sys_no, p1, seed++);
    systematicbits1[0] = 1;

    for(int i=0; i<systematicbits1.size(); i++) {
      int random = rand() % 1000000;
      if(random < geneThreshold)
	systematicbits2[i] = 1 ^ systematicbits1[i];
      else
	systematicbits2[i] = systematicbits1[i];
    }
    systematicbits2[0] = 0;


    vector<int> systrans1(systematicbits1);
    vector<int> systrans2(systematicbits2);

    // RP symbol Encoder
    linearEncoder(symnodeinfo1, systematicbits1, rpsymbols1);
    linearEncoder(symnodeinfo2, systematicbits2, rpsymbols2);

    // Digital bits Encoder
    encoder(codednodeinfo1, systematicbits1, codedbits1);     
    encoder(codednodeinfo2, systematicbits2, codedbits2);
    
    // Generate channel symbol for systematic bits.
    vector<double> channelsignals;
    for(int i=0; i<sys_no; i++) {
      channelsignals.push_back(2 * systematicbits1[i] - 1 + 2 * systematicbits2[i] - 1);
    }

    // Use PAM for RP symbols.
    for(int i=0; i<symbol_no; i++)
      channelsignals.push_back( (rpsymbols1[i] + rpsymbols2[i]) * norfactor);
    
    // Generate channel symbol for digital bits. 
    for(int i=0;i<coded_no;i++) {                                 
      channelsignals.push_back(2 * codedbits1[i] - 1 + 2 * codedbits2[i] - 1);
    }
    
    // Send symbols through complex AWGN channel
    AWGNChannel(channelsignals, sigma);                 

    // Calculate the pdf of corrupted RP symbols (Initialization for symbol nodes)    
    vector<double> allZeroVec(sympdfsize, 0);
    for(int index=0;index<symbol_no;++index) {
      double y = channelsignals[sys_no + index];

      vector<double> symdis(sympdfsize);
      
      for(int symindex=0;symindex<symbolset.size();++symindex)
	{
	  int zepoint = (sympdfsize - 1) / 2;
	  double symVal = symbolset[symindex];
	  int symPos = symVal + zepoint;
	  symdis[symPos] =  getSymPdf(symbolset, allZeroVec, symVal, norfactor, y, sigma); 
	}
      vecNorm(symdis);
      channelpdf[index] = symdis;
      symnodeinfo1[index].clearMessage();
      symnodeinfo2[index].clearMessage();
    }

    // Initialization for coded nodes
    for(int index=0;index<coded_no;++index) {
      
      double y = channelsignals[sys_no + symbol_no + index];

      coded_channeloutput[index] = log((exp(-pow(y + 2 * sqrt(Es1), 2) / (2 * pow(sigma, 2)))	
					+ exp(-pow(y, 2) / (2 * pow(sigma, 2)))) /		
				       (exp(-pow(y - 2 * sqrt(Es1),2)/(2*pow(sigma,2)))	
					+ exp(-pow(y, 2) / (2 * pow(sigma, 2)))));                
      
      codednodeinfo1[index].clearMessage();
      codednodeinfo2[index].clearMessage();
    }
	
    // Initialization for sys nodes
    for(int index=0;index<sys_no;++index) {

      double y = channelsignals[index];

      sys_channeloutput[index] = log((exp(-pow(y + 2 * sqrt(Es1), 2) / (2 * pow(sigma, 2)))	
				      + exp(-pow(y, 2) / (2 * pow(sigma, 2)))) /		
				     (exp(-pow(y - 2 * sqrt(Es1),2)/(2*pow(sigma,2)))	
				      + exp(-pow(y, 2) / (2 * pow(sigma, 2)))));    

      sysnodeinfo1[index].clearMessage();
      sysnodeinfo2[index].clearMessage();
    }
    
    // side information
    vector<double> sideInfo(sys_no, 0);
    double csiThreshold = 0;
 
    int iterationtime = 200;
    int currenttime =0;  
    double watchllr[sys_no];

    vector<int> predecision1(sys_no, 0);
    vector<int> predecision2(sys_no, 0);
    int sametime = 0;
    while(currenttime != iterationtime+1) {

      ++currenttime;
      
      // Decoder 1
      for(int i=0;i<symbol_no;++i)
	symnodeinfo1[i].computePdf(sysnodeinfo1, channelpdf[i], norfactor, noiseVar, 1, i);
     
      for(int i=0;i<coded_no;++i)
	codednodeinfo1[i].computeMessage(sysnodeinfo1, coded_channeloutput[i], 2);

      for(int i=0;i<sys_no;++i)
	sysnodeinfo1[i].computeMessage(symnodeinfo1, codednodeinfo1, sys_channeloutput[i], sideInfo[i]);
   
      // Channel state node
       for(int index=0; index<symbol_no; index++) {
	
	 double y = channelsignals[sys_no + index];
	 
	 vector<double> estimate = symnodeinfo1[index].pdf_to_state;
	 
	 vector<double> symdis(sympdfsize);
      
	 for(int symindex=0;symindex<symbolset.size(); symindex++)
	   {
	     int zepoint = (sympdfsize - 1) / 2;
	     double symVal = symbolset[symindex];
	     int symPos = symVal + zepoint;
	     symdis[symPos] =  getSymPdf(symbolset, estimate, symVal, norfactor, y, sigma); 
	   }
	 vecNorm(symdis);
	 channelpdf[index] = symdis;
	
       }
       
      for(int index=0;index<coded_no; index++) {
	
	double y = channelsignals[sys_no + symbol_no + index];
	
	double estimate=codednodeinfo1[index].message_to_state;
	double nominator_ele1 = exp(-pow(y,2)/(2*pow(sigma,2)))/exp(estimate);
	double nominator_ele2 = (exp(-pow(y+2*sqrt(Es1),2)/(2*pow(sigma,2))));
	double denominator_ele1 = exp(-pow(y-2*sqrt(Es1),2)/(2*pow(sigma,2)))/exp(estimate);
	double denominator_ele2 = (exp(-pow(y,2)/(2*pow(sigma,2))));

	if( isinf(nominator_ele1) && isinf (denominator_ele1) ) 
	  coded_channeloutput[index] = log( exp(-pow(y,2)/(2*pow(sigma,2))) / exp(-pow(y-2*sqrt(Es1),2)/(2*pow(sigma,2))));
	else
	  coded_channeloutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

	if( isnan(coded_channeloutput[index]) || isinf(coded_channeloutput[index]) ) 
	    cout<<"sys channel message not a number!"<<endl;

	
      }

      for(int index=0; index<sys_no; index++) {

	double y = channelsignals[index];

	double estimate = sysnodeinfo1[index].message_to_state;
	double nominator_ele1 = exp(-pow(y,2)/(2*pow(sigma,2)))/exp(estimate);
	double nominator_ele2 = (exp(-pow(y+2*sqrt(Es1),2)/(2*pow(sigma,2))));
	double denominator_ele1 = exp(-pow(y-2*sqrt(Es1),2)/(2*pow(sigma,2)))/exp(estimate);
	double denominator_ele2 = (exp(-pow(y,2)/(2*pow(sigma,2))));

	if( isinf(nominator_ele1) && isinf (denominator_ele1) ) 
	  sys_channeloutput[index] = log( exp(-pow(y,2)/(2*pow(sigma,2))) / exp(-pow(y-2*sqrt(Es1),2)/(2*pow(sigma,2))));
	else
	  sys_channeloutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

	if( isnan(coded_channeloutput[index]) || isinf(coded_channeloutput[index]) ) 
	    cout<<"sys channel message not a number!"<<endl;
	
      }

      for(int index=0; index<sys_no; index++) {
	
	double x = sysnodeinfo1[index].llrEstimation;
        
	x = exp(x);

	if( isinf(x) )
	  x = log( (1-p)/p );
	else
	  x = log((1+(1-p)/p*x)/(x+(1-p)/p));

	if( abs(x) >= csiThreshold )
	  sideInfo[index] = x;
	
	if( isnan(sideInfo[index]) || isinf(sideInfo[index]) )
	    cout<<"side info2 wrong!"<<endl;	
      }
    

      // Decoder 2
      for(int i=0;i<symbol_no;++i)
	symnodeinfo2[i].computePdf(sysnodeinfo2, channelpdf[i], norfactor, noiseVar, 1, i);
     
      for(int i=0;i<coded_no;++i)
	codednodeinfo2[i].computeMessage(sysnodeinfo2, coded_channeloutput[i], 2);

      for(int i=0;i<sys_no;++i)
	sysnodeinfo2[i].computeMessage(symnodeinfo2, codednodeinfo2, sys_channeloutput[i], sideInfo[i]);
      
      // Channel nodes
      for(int index=0; index<symbol_no; index++) {
	
	double y = channelsignals[sys_no + index];

	 vector<double> estimate=symnodeinfo2[index].pdf_to_state;
	 
	 vector<double> symdis(sympdfsize);
      
	 for(int symindex=0;symindex<symbolset.size(); symindex++)
	   {
	     int zepoint = (sympdfsize - 1) / 2;
	     double symVal = symbolset[symindex];
	     int symPos = symVal + zepoint;
	     symdis[symPos] =  getSymPdf(symbolset, estimate, symVal, norfactor, y, sigma); 
	   }
	 vecNorm(symdis);
	 channelpdf[index] = symdis;
	
       }
       
      for(int index=0;index<coded_no;++index) {
	double y = channelsignals[sys_no + symbol_no + index];

	double estimate=codednodeinfo2[index].message_to_state;
	double nominator_ele1 = exp(-pow(y,2)/(2*pow(sigma,2)))/exp(estimate);
	double nominator_ele2 = (exp(-pow(y+2*sqrt(Es2),2)/(2*pow(sigma,2))));
	double denominator_ele1 = exp(-pow(y-2*sqrt(Es2),2)/(2*pow(sigma,2)))/exp(estimate);
	double denominator_ele2 = (exp(-pow(y,2)/(2*pow(sigma,2))));

	if( isinf(nominator_ele1) && isinf (denominator_ele1) ) 
	  coded_channeloutput[index] = log( exp(-pow(y,2)/(2*pow(sigma,2))) / exp(-pow(y-2*sqrt(Es2),2)/(2*pow(sigma,2))));
	else
	  coded_channeloutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

	if( isnan(coded_channeloutput[index]) || isinf(coded_channeloutput[index]) ) 
	    cout<<"sys channel message not a number!"<<endl;
      }
      
      for(int index=0; index<sys_no; index++) {
	
	double y = channelsignals[index];
	
	double estimate = sysnodeinfo2[index].message_to_state;
	double nominator_ele1 = exp(-pow(y,2)/(2*pow(sigma,2)))/exp(estimate);
	double nominator_ele2 = (exp(-pow(y+2*sqrt(Es1),2)/(2*pow(sigma,2))));
	double denominator_ele1 = exp(-pow(y-2*sqrt(Es1),2)/(2*pow(sigma,2)))/exp(estimate);
	double denominator_ele2 = (exp(-pow(y,2)/(2*pow(sigma,2))));

	if( isinf(nominator_ele1) && isinf (denominator_ele1) ) 
	  sys_channeloutput[index] = log( exp(-pow(y,2)/(2*pow(sigma,2))) / exp(-pow(y-2*sqrt(Es1),2)/(2*pow(sigma,2))));
	else
	  sys_channeloutput[index] = log((nominator_ele1+nominator_ele2) / (denominator_ele1+denominator_ele2));

	if( isnan(sys_channeloutput[index]) || isinf(sys_channeloutput[index]) ) 
	    cout<<"sys channel message not a number!"<<endl;

      }

       for(int index=0; index<sys_no; index++) {
	
	double x = sysnodeinfo2[index].llrEstimation;
        
	x = exp(x);

	if( isinf(x) )
	  x = log( (1-p)/p );
	else
	  x = log((1+(1-p)/p*x)/(x+(1-p)/p));

	if( abs(x) >= csiThreshold )
	  sideInfo[index] = x;
	
	if( isnan(sideInfo[index]) || isinf(sideInfo[index]) )
	    cout<<"side info 1 wrong!"<<endl;	
      }

      // Count error number.
      double errorbit = 0;
      
      if(sysnodeinfo1[0].decision == 0) {
	errorbit += countError(sysnodeinfo1, systrans1);
	errorbit += countError(sysnodeinfo2, systrans2);
      } else if(sysnodeinfo1[0].decision == 1) {
	errorbit += countError(sysnodeinfo1, systrans2);
	errorbit += countError(sysnodeinfo2, systrans1);
      } 
	   
      cout<<"error number is:"<<errorbit<<"("<<currenttime<<")"<<endl;

      error<<errorbit<<" ";      

      double diffPos = 0;
      for(int i=0; i<sys_no; i++) {
	if(sysnodeinfo1[i].decision != sysnodeinfo2[i].decision)
	  diffPos++;
      }
      cout << "Number of different decisions: " << diffPos << endl;

      if(sameDecision(predecision1, sysnodeinfo1) && sameDecision(predecision2, sysnodeinfo2))
	++sametime;
      else
	sametime=0;

      if(sametime == 5) {
	cout<<" Make same decision for 5 consecutive iterations !"<<endl;
	break;
      }

    }
  
    double errorbit = 0;
    
    if(sysnodeinfo1[0].decision == 0) {
      errorbit += countError(sysnodeinfo1, systrans1);
      errorbit += countError(sysnodeinfo2, systrans2);
    } else if(sysnodeinfo1[0].decision == 1) {
      errorbit += countError(sysnodeinfo1, systrans2);
      errorbit += countError(sysnodeinfo2, systrans1);
    } 
       

    double blockber = errorbit / (2 * sys_no);
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

  // P{y | x1 = k} = sum_x2 P(x2) * P{x2 + n = y}
  // P(x2) belongs to prior knowldge, currently it is assumed that all the symbols
  // in the symbol set have the same probability.
  double zeroPos = (estimate.size() - 1) / 2;
  for(int i=0; i<symbolSet.size(); i++) {
    double x2Pro = estimate[symbolSet[i] + zeroPos];
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
	  int sys_nodeNum=codedNodeP[oldNum].nodeData.neighbourNum[i];
	  int coded_index=codedNodeP[oldNum].nodeData.inmessageIndex[i];
	  
	  if(group==1)
	    sysNode[sys_nodeNum].nodeData.neighbourNum[coded_index]=newNum;
	  if(group==2)
	    sysNode[sys_nodeNum].nodeData.para_neighbourNum[coded_index]=newNum;
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

int countError(vector<Sys_info> &sysnodeinfo, vector<int> &systrans) {

  int errorNum = 0;
  for(int i=0; i<systrans.size(); i++) {
    if(sysnodeinfo[i].decision != systrans[i])
      errorNum++;
  }
  
  return errorNum;
}
