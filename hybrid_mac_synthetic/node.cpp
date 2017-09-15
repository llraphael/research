#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include "node.h"
#include "func.h"


const double threshold = 30;

Node_info::Node_info(int x):
  degree(x),
  para_degree(0)
{
  neighbourNum.resize(x);
  inmessageIndex.resize(x);
  outMessage.resize(x);

}

Node_info::Node_info(int a, int b):
  degree(a),
  para_degree(b)
{

  neighbourNum.resize(a);
  inmessageIndex.resize(a);
  outMessage.resize(a);

  para_neighbourNum.resize(b);
  para_inmessageIndex.resize(b);
  para_outMessage.resize(b);

}

Coded_info :: Coded_info () :
  nodeData(),
  message_to_state(0),
  llrEstimation(0),
  decision(0)
{}

Coded_info :: Coded_info (int x) : 
nodeData(x), 
message_to_state(0),
llrEstimation(0),
decision(0) 
{};

int Coded_info::readNeighbourNum(std::multimap<int,int> &connection, int index)
{
  int degree = connection.count(index);

  nodeData.degree = degree;
  
  nodeData.neighbourNum.resize(nodeData.degree);
  nodeData.inmessageIndex.resize(nodeData.degree);
  nodeData.outMessage.resize(nodeData.degree);


  std::multimap<int,int>::iterator noReading;
  std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> pointer;
  

  //make sure which column or row it need to read.
  pointer=connection.equal_range(index);
  
  // If the matrix is for adding connection, start from new pos
  int i = 0;

  for(noReading=pointer.first;noReading!=pointer.second;++noReading)
    {
	 
      nodeData.neighbourNum[i]=(*noReading).second;
      ++i;
    }
}

int Coded_info::readLinearMatrix(Matrix &ma, int index)
{
  std::vector<ElementData> colpos(ma.matrix_data[index].colPos);

  int degree = colpos.size();
  nodeData.degree = degree;

  nodeData.neighbourNum.resize(degree);
  nodeData.inmessageIndex.resize(degree);
  nodeData.outMessage.resize(degree);

  weightset.resize(degree);

  for(int i=0;i<degree;++i)
    {
      nodeData.neighbourNum[i] = colpos[i].column_pos;
      weightset[i] = colpos[i].value;
    }
}

//neidata is the first node
void Coded_info::getsourceNum(std::vector<Sys_info> &neidata,int nodeNum, int identity)
{


  int degree=nodeData.degree;
  
  //suppose first group of coded bits want to find source index
  if(identity==1)
    {
      for(int neighbourIndex=0;neighbourIndex!=degree;++neighbourIndex)
	{
	  int neiNum=nodeData.neighbourNum[neighbourIndex];
	  int neiDegree = neidata[neiNum].nodeData.degree;
	  for(int i=0;i!=neiDegree;++i)
	    {
	      if(neidata[neiNum].nodeData.neighbourNum[i]==nodeNum)
		nodeData.inmessageIndex[neighbourIndex]=i;
	    }
	}
    }

  //suppose parallel group of coded bits want to find source index
  if(identity==2)
    {
      for(int neighbourIndex=0;neighbourIndex!=degree;++neighbourIndex)
	{
	  int neiNum=nodeData.neighbourNum[neighbourIndex];
	  int neiDegree = neidata[neiNum].nodeData.para_degree;
	  for(int i=0;i!=neiDegree;++i)
	    {
	      if(neidata[neiNum].nodeData.para_neighbourNum[i]==nodeNum)
		nodeData.inmessageIndex[neighbourIndex]=i;
	    }
	}
    }
}

void Coded_info :: assignWeight(std::vector<double> weight, int seed)
{
  if(weight.size() != nodeData.degree) {
    std::cout<<"weight num does not match the degree!"<<std::endl;
    exit(EXIT_FAILURE);
  }

  weightset.resize(nodeData.degree);

  std::vector<int> ranarr(nodeData.degree);
  for(int i=0;i<ranarr.size();++i)
    ranarr[i] = i;

  srand(seed);

  for(int i=0;i<ranarr.size();++i)
    {
      int random = rand()%(ranarr.size()-i);

      int temp = ranarr[ranarr.size()-1];
      ranarr[ranarr.size()-1] = ranarr[random];
      ranarr[random] = temp;

    }

  for(int i=0;i<weight.size();++i)
    weightset[i] = weight[ ranarr[i] ];

}


double Coded_info::getMessage(int messageIndex)
{
  return nodeData.outMessage[messageIndex];
}

void Coded_info::setMessage(double message, int messageIndex)
{
  nodeData.outMessage[messageIndex]=message;
}

void Coded_info::clearMessage()
{
  for(int i=0;i<nodeData.degree;++i)
    setMessage(0,i);
}

// parameter identity is used for getting the right message in neighbourgh sys nodes 
// in parallel case.
void Coded_info::computeMessage(std::vector<Sys_info> &neiNode, double channelInfo, int identity)
{
  int degree=nodeData.degree;
  std:: vector<double> inmessage(degree,0);
  
  for(int index=0;index<degree;++index)
    {
      int neiNum=nodeData.neighbourNum[index];
      int messageIndex=nodeData.inmessageIndex[index];
      double message=neiNode[neiNum].getMessage(messageIndex, identity);
 
      inmessage[index]=message;
    }
   
  double zero_message_num = 0;

  double x=tanh(channelInfo/2);
  for(int index=0;index<degree;++index)
    {
      if(inmessage[index]!=0) 
				x *= tanh(inmessage[index]/2);
      else
				++zero_message_num;
      
    }

  if(zero_message_num==1) {
    
    for(int index=0;index<degree;++index) {
			if(inmessage[index]==0) {
	  		double u = x;

	  		if(u==1) 
	    		u = threshold;
	  		else if(u==-1)
	    		u = -threshold;
	  		else 
	    		u = 2 * atanh(u);

	  		nodeData.outMessage[index] = u;
			}
			else 
				nodeData.outMessage[index] = 0;
    }   
	} 
  
  else if(zero_message_num>1) {
    for(int index=0;index<degree;++index)
      nodeData.outMessage[index] = 0;
  }
  
  else {

  //compute the outmessage(U) one by one for the code
  for(int outIndex=0;outIndex<degree;++outIndex)
    {
      double u = tanh(inmessage[outIndex]/2);
      
      u = x / u;
      
      if( u<-1 || u>1)
				std::cout<<"check message not senseful"<<std::endl;      

      if(u==1) 
				u = threshold;
      else if(u==-1)
				u = -threshold;
      else 
        u = 2 * atanh(u);

      if(isinf(u) || isnan(u))
				std::cout<<"warning! wrong messages!"<<u<<std::endl;

      nodeData.outMessage[outIndex]=u;
    }

  }
   
  //compute the message to the state node
  if(channelInfo != 0)
		message_to_state = x / tanh(channelInfo/2);
 	else if(zero_message_num >= 1)
		message_to_state = 0;
	else {
		message_to_state = 1;
		for(int i=0; i<degree; i++)
			message_to_state *= tanh(inmessage[i] / 2);
	}

  if(message_to_state==1) 
    message_to_state=threshold;
 
  else if(message_to_state==-1)
    message_to_state=-threshold;
  else 
    message_to_state=2*atanh(message_to_state);
  

  if(message_to_state >threshold)
    message_to_state=threshold;
  else if(message_to_state<-threshold)
    message_to_state=-threshold;
  else ;

      
  if(isinf(message_to_state) || isnan(message_to_state) )
    std::cout<<"warning! check channel message wrong!"<<std::endl;   
}

const bool Coded_info::compStruc(const Coded_info& rhs)
{ 

  bool returnvalue=0;
  if(nodeData.neighbourNum != rhs.nodeData.neighbourNum) 
    return returnvalue;
  else if (nodeData.inmessageIndex != rhs.nodeData.inmessageIndex)
    return returnvalue;
  else
    return !returnvalue;

}

double Coded_info :: get_ratio_estimation()
{
  return llrEstimation;
}

std::vector<double> testCombinedPdf;
void Coded_info :: computePdf(std::vector<Sys_info> &neiNode, std::vector<double> &channelpdf, double norFactor, double noriseVar, int identity, int number)
{

  int degree=nodeData.degree;
  std::vector<double> inmessage(degree);
  
  for(int index=0;index<degree;++index)
    {
      int neiNum = nodeData.neighbourNum[index];
      int messageIndex = nodeData.inmessageIndex[index];
      inmessage[index] = neiNode[neiNum].getMessage(messageIndex, identity);
 
    }

  std::vector<double> a(3);
  std::vector< std::vector<double> > syspdf(degree, a);     //construct pdf of sysbits from llr message.

  for(int index=0;index<degree;++index)
    {
      double mess = inmessage[index];
      syspdf[index][1] = exp(mess)/(1+exp(mess));       //probability of zero
      syspdf[index][2] = 1 - syspdf[index][1];          //probability of one
    }

  std::vector<double> norWeight(degree);
  for(int i=0; i<degree; i++)
    norWeight[i] = norFactor * weightset[i];
  
  // Calculate the variance of combined value
  double var = 0;
  for(int i=0; i<degree; i++) 
    var += pow(norWeight[i], 2) * (pow(syspdf[i][2], 2) * syspdf[i][1] + pow(1-syspdf[i][2], 2) * syspdf[i][2]); 
				   
  std::vector< std::vector<double> > sysalias(syspdf);
  for(int index=0;index<degree;++index)
    weightPdf(weightset[index], sysalias[index]);        // Get weighted pdf of sysbits
  
  
  std::vector<double> sumpdf(1,1);                     
  for(int index=0;index<degree;++index)                // convolute all the pdf
    conv(sumpdf, sysalias[index], sumpdf);

  vecNorm(sumpdf, channelpdf.size());

  pdf_to_state = sumpdf;
  for(int index=0;index<degree;++index)
    {

      std::vector<double> respdf(1, 1);
 
      deconvolve(sumpdf, syspdf[index], weightset[index], respdf);

      double p0 = 0, p1 = 0;
      double zeropoint = (respdf.size()-1)/2;
      double message = 0;
            
      // notice that respdf and channelpdf has the same length, so they
      // have the same reference point in the vector
      double rpThreshold = 0;//pow(10, -100);
      for(int symI = 0; symI < channelpdf.size(); symI++) {
				if(channelpdf[symI] < rpThreshold)
	  			continue;

				p0 += respdf[symI] * channelpdf[symI];
				int comPos = symI - weightset[index];
				if(comPos >= 0 && comPos < respdf.size())
	  		p1 += respdf[comPos] * channelpdf[symI];
	
      }
      
      if(p0 == 0 || p1 == 0)
				message = 0;
      else {
				message = log(p0 / p1);
				if(message > threshold)
	  			message = threshold;
				if(message < -threshold)
	  			message = -threshold;
				if(isnan(message))
	  			std::cout<<"invalid message!"<<std::endl;
				setMessage(message, index);
      }
    }
}

void Coded_info :: computePdfApp(std::vector<Sys_info> &neiNode, std::vector<double> &channelpdf, std::map<int, std::vector<std::vector<int> > > &mappingTable, std::vector<double> &weights, int identity) {
   
  int degree=nodeData.degree;
  std::vector<double> inmessage(degree,0);
   
  for(int index=0;index<degree;++index)
    {
      int neiNum=nodeData.neighbourNum[index];
      int messageIndex=nodeData.inmessageIndex[index];
      double message=neiNode[neiNum].getMessage(messageIndex, identity);
 
      inmessage[index]=message;
    }

  std::vector<double> inZeroPro(degree, 0);
  std::vector<double> inOnePro(degree, 0);
  for(int index=0;index<degree;++index)
    {
      double mess = inmessage[index];
      inZeroPro[index] = exp(mess)/(1+exp(mess));       //probability of zero
      inOnePro[index] = 1 - inZeroPro[index];          //probability of one
    }
  

  std::vector<double> outZeroPro(degree, 0);
  std::vector<double> outOnePro(degree, 0);

  // Real position of each bits due to the shuffling of the weight.
  std::vector<int> bitRealPos(degree, 0);
  std::vector<double> realWeights(weightset);
  for(int i=0; i<degree; i++) {
    int oneW = weights[i];

    // Search for the index of oneW in realWeights
    for(int j=0; j<degree; j++) {
      if(oneW == realWeights[j]) {
	bitRealPos[i] = j;
	realWeights[j] = 0;
	break;
      }
    }

  }
  // For each value of RP symbol, the corresponding bit values can be known from the table, so that the 
  // probability of that bit being that value can be computed from the pro of RP symbol and other bits.
  // And RP values with pro smaller than pow(10, -8) are ignored
  double rpThreshold = pow(10, -2);
  int zeroPoint = (channelpdf.size() - 1) / 2;
  for(int i = 0; i < channelpdf.size(); i++) {
    if(channelpdf[i] < rpThreshold)
      continue;

    int rpVal = i - zeroPoint;
    std::vector<std::vector<int> > bitComs = mappingTable[rpVal];
    for(int comI = 0; comI < bitComs.size(); comI++) {
      std::vector<int> temp = bitComs[comI];
      std::vector<int> oneCom(temp.size());
      for(int bitI = 0; bitI < degree; bitI++)
	oneCom[bitRealPos[bitI]] = temp[bitI];

      // compute the product of all the incoming pro and RP pro based on current RP value
      double proAll = channelpdf[i];
      for(int bitI = 0; bitI < degree; bitI++) {
	if(oneCom[bitI] == 0)
	  proAll *= inZeroPro[bitI];
	else
	  proAll *= inOnePro[bitI];
      }

      // start to compute the pro for each bit being that value
      for(int bitI = 0; bitI < degree; bitI++) {
	double bitPro = (oneCom[bitI] == 0 ? proAll / inZeroPro[bitI] : proAll / inOnePro[bitI]);
	if(oneCom[bitI] == 0)
	  outZeroPro[bitI] += bitPro;
	else
	  outOnePro[bitI] += bitPro;
      }
    }
  }

  for(int bitI = 0; bitI < degree; bitI++) {
    double pro0 = outZeroPro[bitI];
    double pro1 = outOnePro[bitI];
    double message = 0;
    if(pro0 == 0 && pro1 == 0)
      message = 0;
    else if(pro0 == 0)
      message = -threshold;
    else if(pro1 == 0)
      message = threshold;
    else {
      message = log(pro0 / pro1);
      if(message > threshold)
	message = threshold;
      if(message < -threshold)
	message = -threshold;
    }
    setMessage(message, bitI);
  }
  
} 

void Coded_info :: allocateForMeanAndVar() {
  nodeData.outMean.resize(nodeData.degree);
  nodeData.outVar.resize(nodeData.degree);
}

void Coded_info :: computeInMeanAndVar(std::vector<Sys_info> &neiNode, double receivedSym, std::map<int, double> &symPro, double norFactor, double noiseVar, int identity, std::string format) {
  
  int degree=nodeData.degree;
  
  std::vector<double> inMean(degree, 0);
  std::vector<double> inVar(degree, 0);
  
  if(format == "llr") {
    std::vector<double> inMessage(degree, 0);
    for(int index=0;index<degree;++index)
      {
	int neiNum = nodeData.neighbourNum[index];
	int messageIndex = nodeData.inmessageIndex[index];

	inMessage[index] = neiNode[neiNum].getMessage(messageIndex, identity);
      }
    
    for(int i=0; i<degree; i++) {
      double p1 = 1 / (1 + exp(inMessage[i]));
      double p0 = 1 - p1;
      inMean[i] = 0 * p0 + 1 * p1;
      //inVar[i] =  0.0001 + 1 + 2 * (p1 - p0) * inMean[i] + pow(inMean[i], 2); //Gau
      inVar[i] = pow(0 - inMean[i], 2) * p0 + pow(1 - inMean[i], 2) * p1;
    }
  }
  if(format == "Gau") {
    for(int index=0;index<degree;++index)
      {
	int neiNum = nodeData.neighbourNum[index];
	int messageIndex = nodeData.inmessageIndex[index];
	inMean[index] = neiNode[neiNum].nodeData.outMean[messageIndex];
	inVar[index] = neiNode[neiNum].nodeData.outVar[messageIndex];
      }
  }
  
  // weight * norFactor
  std::vector<double> tempW(degree, 0);
  for(int i=0; i<degree; i++)
    tempW[i] = weightset[i] * norFactor;

  /*
  // Get mean sum and mean variance
  double meanSum = -1 * receivedSym;
  for(int i=0; i<degree; i++)
    meanSum += tempW[i] * inMean[i];
  
  double varSum = noiseVar;
  for(int i=0; i<degree; i++)
    varSum += pow(tempW[i], 2) * inVar[i];

  // compute each outgoing mean and variance
  for(int i=0; i<degree; i++) {
    nodeData.outMean[i] = -1 / tempW[i] * (meanSum - tempW[i] * inMean[i]);
    nodeData.outVar[i] = 1 / pow(tempW[i], 2) * (varSum - pow(tempW[i], 2) * inVar[i]);
  }
  // compute each outgoing llr message
  for(int i=0; i<degree; i++) {
    double p0 = gaussianFunc(-1, nodeData.outMean[i], sqrt(nodeData.outVar[i]));
    double p1 = gaussianFunc(1, nodeData.outMean[i], sqrt(nodeData.outVar[i]));
   double message = log(p0 / p1);
    if(message > threshold)
      message = threshold;
    if(message < -threshold)
      message = -threshold;

    setMessage(message, i);
  }
  */

  
  // New method
  double symMinVal = symPro.begin() -> first;
  double symMaxVal = symPro.rbegin() -> first;
  double meanSum = 0;
  for(int i=0; i<degree; i++) 
    meanSum += weightset[i] * inMean[i];
  double varSum = 0;
  for(int i=0; i<degree; i++)
    varSum += pow(weightset[i], 2) * inVar[i];

  for(int i=0; i<degree; i++) {
    double p0 = 0, p1 = 0;
    double otherMean = meanSum - weightset[i] * inMean[i];
    double otherSig = sqrt(varSum - pow(weightset[i], 2) * inVar[i]);
    
    std::map<int, double>::iterator it = symPro.begin();
    for(it; it != symPro.end(); it++) {
      //if(it -> second < pow(10, -10))
      //continue;
      //else {
	double symVal = it -> first;
	// pro being 1
	double comVal = symVal - weightset[i] * 1;
	if(comVal >= symMinVal && comVal <= symMaxVal)
	  p1 += (it->second) * gaussianFunc(comVal, otherMean, otherSig);
	// pro being 0
	comVal = symVal - weightset[i] * (0);
	if(comVal >= symMinVal && comVal <= symMaxVal)
	  p0 += (it->second) * gaussianFunc(comVal, otherMean, otherSig);
	//}
    }

    double message = log(p0 / p1);
    if(p0 == 0 || p1 == 0)
      message = 0;

    if(message > threshold)
      message = threshold;
    if(message < -threshold)
      message = -threshold;
    
    if(isnan(message))
      std::cout << "problem" << std::endl;

    setMessage(message, i);
  }
  

}


Sys_info :: Sys_info () :
  nodeData(),
  message_to_state(0),
  llrEstimation(0),
  decision(0)
{};

Sys_info :: Sys_info (int a) : 
  nodeData(a),
  message_to_state(0),
  llrEstimation(0),
  decision(0)
{};

Sys_info :: Sys_info (int a, int b) : 
  nodeData(a, b),
  message_to_state(0),
  llrEstimation(0),
  decision(0)
{};

int Sys_info::readNeighbourNum(std::multimap<int,int> &connection, int index, int group, int addConnection)
{
  if(group == 1) {
    int degree = connection.count(index);
    int oridegree = nodeData.degree;

    // If the node already has connections, the reading would be 
    // adding connections to the node
    if(addConnection) 
      nodeData.degree += degree;
    else
      nodeData.degree = degree;

    resize(nodeData.degree, 1);

    std::multimap<int,int>::iterator noReading;
    std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> pointer;

    //make sure which column or row it need to read.
    pointer=connection.equal_range(index);
     
    int i=0;
    if(addConnection)
      i = oridegree;
    else
      i = 0;

    for(noReading=pointer.first;noReading!=pointer.second;++noReading)
      {
	nodeData.neighbourNum[i]=(*noReading).second;
	++i;
      }
  }
  if(group == 2) 
    readParaNeighbourNum(connection, index, addConnection);
}


// Get the new neighbour information in rateless scheme.
int Sys_info::readNeighbourNum(int current_coded_num, std::multimap<int,int> &connection, int index)
{
  int current_size = nodeData.degree;

  int increase_size = connection.count(index);
  resize(increase_size+current_size, 1);

  std::multimap<int,int>::iterator noReading;
  std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> pointer;

  //make sure which column or row it need to read.
  pointer=connection.equal_range(index);
     
  int i=0;
  for(noReading=pointer.first;noReading!=pointer.second;++noReading)
    {
      
      nodeData.neighbourNum[current_size+i]=(*noReading).second+current_coded_num;
      ++i;
    }


}


int Sys_info::readParaNeighbourNum(std::multimap<int,int> &connection, int index, int addConnection)
{
  int para_degree=connection.count(index);
  int oridegree = nodeData.para_degree;

  if(addConnection) 
    nodeData.para_degree += para_degree;
  else
    nodeData.para_degree = para_degree;

  
  resize(nodeData.para_degree, 2);

  std::multimap<int,int>::iterator noReading;
  std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> pointer;
  
  //make sure which column or row it need to read.
  pointer=connection.equal_range(index);
     
  int i=0;
  if(addConnection)
    i = oridegree;
  else
    i = 0;

  for(noReading=pointer.first;noReading!=pointer.second;++noReading)
    {
      nodeData.para_neighbourNum[i]=(*noReading).second;
      ++i;
    }
       
}


int Sys_info::readLinearMatrix(Matrix &ma, int index, int codedNo, int group)
{

  
  std::multimap<int,int>::iterator noReading;
  std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> pointer;

  resize(0, 1);

  //make sure which column or row it need to read.
  pointer=ma.matrixbycolumn.equal_range(index);
     
  for(noReading=pointer.first;noReading!=pointer.second;++noReading)
    {
      if( (*noReading).second < codedNo)
	nodeData.neighbourNum.push_back( (*noReading).second );
    }

  int degree = nodeData.neighbourNum.size();
  resize(degree, 1);

  /*
  if(group == 1)
    readNeighbourNum(ma.matrixbycolumn, index);
  if(group == 2)
    readParaNeighbourNum(ma.matrixbycolumn, index);
  */
}

 
	
void Sys_info::getsourceNum(std::vector<Coded_info> &neidata, int nodeNum)
{
  int degree=0;
  
  // get the degree of the neighbour coded node
   
   degree=nodeData.degree;

  // read the source index from first group
   for(int neighbourIndex=0;neighbourIndex!=degree;++neighbourIndex)
     {
       int neiNum=nodeData.neighbourNum[neighbourIndex];
       int neiDegree = neidata[neiNum].nodeData.degree;
       for(int i=0;i!=neiDegree;++i)
	 {
	   if(neidata[neiNum].nodeData.neighbourNum[i]==nodeNum)
	     nodeData.inmessageIndex[neighbourIndex]=i;
	 }
     }
}

// On sys node side, the value of the third parameter determines  
// source data should be write into which group
void Sys_info::getsourceNum(std::vector<Coded_info> &neidata1,std::vector<Coded_info> &neidata2, int nodeNum)
{

  
  int degree=0;
  
  degree=nodeData.degree;

  // read the source index from first group
  for(int neighbourIndex=0;neighbourIndex!=degree;++neighbourIndex)
    {
      int neiNum=nodeData.neighbourNum[neighbourIndex];
      int neiDegree = neidata1[neiNum].nodeData.degree;
      for(int i=0;i!=neiDegree;++i)
	{
	  if(neidata1[neiNum].nodeData.neighbourNum[i]==nodeNum)
	    nodeData.inmessageIndex[neighbourIndex]=i;
	}
    }
  degree=nodeData.para_degree;

   for(int neighbourIndex=0;neighbourIndex!=degree;++neighbourIndex)
     {
       int neiNum=nodeData.para_neighbourNum[neighbourIndex];
       int neiDegree = neidata2[neiNum].nodeData.degree;
       for(int i=0;i!=neiDegree;++i)
	 {
	   if(neidata2[neiNum].nodeData.neighbourNum[i]==nodeNum)
	     nodeData.para_inmessageIndex[neighbourIndex]=i;
	 }
     }

}
/*
int Sys_info::connect(int sysindex, int codedindex, Coded_info &newneighbour, int group)
{

  increaseSize(1, group);

  int neighbournum = nodeData.neighbourNum.size();
  int para_neighbournum = nodeData.para_neighbourNum.size();

  // Add the new neighbour number
  if(group == 1)
    nodeData.neighbourNum[neighbournum-1] = codedindex;
  if(group == 2)
    nodeData.para_neighbourNum[para_neighbournum-1] = codedindex;

  // Get the index of itself(sys node).
  int codeddegree = newneighbour.degree;
  for(int i=0;i<codeddegree;++i)
    {
      if(newneighbour.nodeData.neighbourNum[i] == sysindex) {

	if(group == 1)
	  nodeData.inmessageIndex[neighbournum-1] = i;
	if(group == 2)
	  nodeData.inmessageIndex[para_neighbournum-1] = i;

	break;
      }
    }




}
*/
double Sys_info::getMessage(int index, int group)
{
  if(group==1)
    return nodeData.outMessage[index];
  if(group==2)
    return nodeData.para_outMessage[index];
}

void Sys_info::setMessage(double message, int index, int group)
{
  if(group==1)
    nodeData.outMessage[index]=message;
  if(group==2)
    nodeData.para_outMessage[index]=message;
}

void Sys_info::clearMessage()
{
  for(int i=0;i<nodeData.degree;++i)
    setMessage(0, i, 1);
  for(int i=0;i<nodeData.para_degree;++i)
    setMessage(0, i, 2);
}


void Sys_info::computeMessage(std::vector<Coded_info> &neiNode, double channelInfo, double sideInfo, int group)
{

  int degree=nodeData.degree;
  if(group == 2)
    degree = nodeData.para_degree;
  
  std::vector<double> inmessage(degree,0);
  
  if(group == 1) {
    for(int index=0;index<degree;++index)
      {
				int neiNum=nodeData.neighbourNum[index];
				int messageIndex=nodeData.inmessageIndex[index];
				double message=neiNode[neiNum].getMessage(messageIndex);
				inmessage[index] = message;
      }
  	}
  	if(group == 2) {
    	for(int index=0;index<degree;++index) {
				int neiNum=nodeData.para_neighbourNum[index];
				int messageIndex=nodeData.para_inmessageIndex[index];
				double message=neiNode[neiNum].getMessage(messageIndex);
				inmessage[index] = message;
      }
  	}  
  
		double v = channelInfo + sideInfo;
  	for(int index=0;index<degree;++index)
    	v += inmessage[index];
  
		//compute the outmessage(U) one by one for the code
  	for(int outIndex=0;outIndex<degree;++outIndex) {
      double message = v - inmessage[outIndex];

      if(message > threshold)
				message = threshold;
      if(message < -threshold)
				message = -threshold;
 
      setMessage(message, outIndex, group);
    }

  //make decision for the node
  llrEstimation=v;
	message_to_state = v - channelInfo;

	if(message_to_state > threshold)
		message_to_state = threshold;
	if(message_to_state < -threshold)
		message_to_state = -threshold;

  if(v>0)
    decision=0;
  else
    decision=1;
}

void Sys_info::computeMessage(std::vector<Coded_info> &neiNode1, std::vector<Coded_info> &neiNode2, double channelInfo, double sideInfo)
{

  int degree = nodeData.degree;
  int para_degree = nodeData.para_degree;

  std::vector<double> inmessage1(degree,0);
  std::vector<double> inmessage2(para_degree, 0);

  for(int index=0;index<degree;++index)
    {
      int neiNum=nodeData.neighbourNum[index];
      int messageIndex=nodeData.inmessageIndex[index];
      double message=neiNode1[neiNum].getMessage(messageIndex);
      inmessage1[index] =  message ;
    }

  for(int index=0;index<para_degree;++index)
    {
      int neiNum=nodeData.para_neighbourNum[index];
      int messageIndex=nodeData.para_inmessageIndex[index];
      double message=neiNode2[neiNum].getMessage(messageIndex);
      inmessage2[index] = message;
    }

  double v = channelInfo + sideInfo;
  for(int inIndex=0;inIndex!=degree;++inIndex)
    v+=inmessage1[inIndex];
  for(int inIndex=0;inIndex!=para_degree;++inIndex)
    v+=inmessage2[inIndex];

  
  //compute the outmessage(U) for the first group of coded bits
  for(int outIndex=0;outIndex<degree;++outIndex)
    {
      double outMessage = v - inmessage1[outIndex];
  
      if(outMessage<-threshold)
				outMessage=-threshold;
      else if(outMessage>threshold)
				outMessage=threshold;
      else ;
      
      setMessage(outMessage, outIndex, 1);
    }

  //compute the outmessage for the parallel group of coded bits
  for(int outIndex=0;outIndex<para_degree;++outIndex)
    {
      double outMessage = v - inmessage2[outIndex];
    
      if(outMessage<-threshold)
				outMessage=-threshold;
      else if(outMessage>threshold)
				outMessage=threshold;
      else ;
      
      setMessage(outMessage, outIndex, 2);
    }


  llrEstimation = v;
  message_to_state = v - channelInfo;

  if( isnan(llrEstimation) || isnan(message_to_state) )
    std::cout<<"sys message wrong!"<<std::endl;
  
  if(message_to_state>threshold)
    message_to_state=threshold;
  else if(message_to_state<-threshold)
    message_to_state=-threshold;
  else ;

  if(llrEstimation>threshold) 
    llrEstimation=threshold;
  else if(llrEstimation<-threshold)
    llrEstimation=-threshold;
  
  
  if(v>0)
    decision=0;
  else
    decision=1;


}
void Sys_info :: allocateForMeanAndVar() {
  nodeData.outMean.resize(nodeData.degree);
  nodeData.outVar.resize(nodeData.degree);
}
//double gaussianFunc(double signal, double mean, double sigma);

void Sys_info :: computeInMeanAndVar(std::vector<Coded_info> &neiNode1, std::vector<Coded_info> &neiNode2, double p0) {
   
  int degree = nodeData.degree;
  int para_degree = nodeData.para_degree;

  std::vector<double> inMean(degree, 0);
  std::vector<double> inVar(degree, 0);
  std::vector<double> inmessage2(para_degree, 0);

  for(int index=0;index<degree;++index)
    {
      int neiNum=nodeData.neighbourNum[index];
      int messageIndex=nodeData.inmessageIndex[index];
      inMean[index] = neiNode1[neiNum].nodeData.outMean[messageIndex];
      inVar[index] = neiNode1[neiNum].nodeData.outVar[messageIndex];
    }

  for(int index=0;index<para_degree;++index)
    {
      int neiNum=nodeData.para_neighbourNum[index];
      int messageIndex=nodeData.para_inmessageIndex[index];
      double message=neiNode2[neiNum].getMessage(messageIndex);
      inmessage2[index] = message;
    }

  // Compute outgoing mean and variance
  double llrSum = log(p0 / (1 - p0));
  for(int i=0; i<para_degree; i++)
    llrSum += inmessage2[i];
  double tempP1 = 1 / (1 + exp(llrSum));
  double tempP0 = 1 - tempP1;
  double rightMean = tempP0 * (-1) + tempP1 * 1;
  double rightVar = 0.0001 + 1 + 2 * (tempP1 - tempP0) * rightMean + pow(rightMean, 2);

  double meanDivVarSum = rightMean / rightVar;
  double varInvSum = 1 / rightVar;
  for(int i=0; i<degree; i++) {
    meanDivVarSum += inMean[i] / inVar[i];
    varInvSum += 1 / inVar[i];
  }
  for(int i=0; i<degree; i++) {
    nodeData.outMean[i] = (meanDivVarSum - inMean[i] / inVar[i]) / (varInvSum - 1 / inVar[i]);
    nodeData.outVar[i] = 1 / (varInvSum - 1 / inVar[i]);
  }

  // Compute llr message to ldgm codes
  meanDivVarSum -= rightMean / rightVar;
  varInvSum -= 1 / rightVar;
  double tempMean = meanDivVarSum / varInvSum;
  double tempSig = sqrt(1 / varInvSum);
  tempP0 = gaussianFunc(-1, tempMean, tempSig);
  tempP1 = gaussianFunc(1, tempMean, tempSig);
  llrSum += log(tempP0 / tempP1);
  for(int i=0; i<para_degree; i++) {
    double message = llrSum - inmessage2[i];
    if(message > threshold)
      message = threshold;
    if(message < -threshold)
      message = -threshold;
    setMessage(message, i, 2);
  }

  //Make decision
  llrEstimation = llrSum;
  if(llrEstimation > 0)
    decision = 0;
  else
    decision = 1;
    
}


double Sys_info::get_ratio_estimation()
{
  return llrEstimation;
}

void Sys_info::resize(int size, int group)
{

  if(group==1)
    {
      nodeData.neighbourNum.resize(size);
      nodeData.inmessageIndex.resize(size);
      nodeData.outMessage.resize(size);
    
      nodeData.degree = size;
    }
  else
    {
      nodeData.para_neighbourNum.resize(size);
      nodeData.para_inmessageIndex.resize(size);
      nodeData.para_outMessage.resize(size);

      nodeData.para_degree = size;
    }
}

const bool Sys_info::compStruc(const Sys_info& rhs)
{ 

  bool returnvalue=0;
  if(nodeData.neighbourNum != rhs.nodeData.neighbourNum) 
    return returnvalue;
  else if (nodeData.inmessageIndex != rhs.nodeData.inmessageIndex)
    return returnvalue;
  else if (nodeData.para_neighbourNum != rhs.nodeData.para_neighbourNum)
    return returnvalue;
  else if (nodeData.para_inmessageIndex != rhs.nodeData.para_inmessageIndex)
    return returnvalue;
  else
    return !returnvalue;

}
