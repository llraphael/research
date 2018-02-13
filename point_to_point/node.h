#ifndef NODE_H
#define NODE_H

#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include "stdlib.h"
#include"sparsematrix.h"


extern std::vector<double> testCombinedPdf;

class Node_info;
class Coded_info;
class Sys_info;

const double pi = 3.1415926;

class Node_info {

 public:

  Node_info(): degree(0), para_degree(0) { };
  Node_info(int x);
  Node_info(int a1, int a2);

 public:

     int degree;
     std::vector<int> neighbourNum; 
     std::vector<int> inmessageIndex;
     std::vector<double> outMessage;

     int para_degree;
     std::vector<int> para_neighbourNum;
     std::vector<int> para_inmessageIndex;
     std::vector<double> para_outMessage;

 public:

     std::vector<double> outMean;
     std::vector<double> outVar;

};



class Coded_info {

public:
     
  Coded_info();
  Coded_info(int x);

public:

 

  //get neighbours' number. 
  int readNeighbourNum(std::multimap<int,int> &matrixbycolumn, int index);
  //For each node, get the storage index of itself in its neighbour's data
  int readLinearMatrix(Matrix &ma, int index);
  void getsourceNum(std::vector<Sys_info>& , int num, int identity);
  void assignWeight(std::vector<double> weight, int seed);  

 public:

  double getMessage(int index);
  void setMessage(double message, int index);

  void computeMessage(std::vector<Sys_info> &neiNode, double initialInfo, int identity);
  void computePdf(std::vector<Sys_info> &neiNode, std::vector<double> &channelpdf,  double norFactor, double noiseVar, int identity, int number);
  void computePdfApp(std::vector<Sys_info> &neiNode, std::vector<double> &channelpdf, std::map<int, std::vector<std::vector<int> > > &mappingTable, std::vector<double> &weights, int identity);
  void computeInMeanAndVar(std::vector<Sys_info> &neiNode, double RPsymbol, std::map<int, double> &table, double norFactor, double noiseVar, int identity, std::string format);
  double get_ratio_estimation();
  void allocateForMeanAndVar();

 public:
    
  void clearMessage();
  const bool compStruc (const Coded_info& rhs);

 public:
     
  Node_info nodeData;
  double message_to_state;
  std::vector<double> pdf_to_state;
  double llrEstimation;
  double decision;

 public:

  std::vector<double> weightset;
};


class Sys_info {

 public:
  
  Sys_info();
  Sys_info(int);	
  Sys_info(int,int);

 public:

  int readNeighbourNum(std::multimap<int, int> &matrixbyrow, int index, int group, int addConnection = 0);
  int readNeighbourNum(int, std::multimap<int,int>&, int);  //for the case of adding column

  int readLinearMatrix(Matrix &ma, int index, int codedNo, int group);
 
  void getsourceNum(std::vector<Coded_info>&, int index);
  void getsourceNum(std::vector<Coded_info>&, std::vector<Coded_info>&, int index);

 private:

  int readParaNeighbourNum(std::multimap<int,int> &, int index, int addConnection);


 public:

  double getMessage(int index, int group);
  void setMessage(double message, int index, int group);
  void clearMessage();
  double get_ratio_estimation();

  const bool compStruc(const Sys_info& rhs);

 public:

  void computeMessage(std::vector<Coded_info> &neiNode, double channelInfo,double sideInfo = 0, int group = 1);  //for single ldgm
  void computeMessage(std::vector<Coded_info> &neiNode1, std::vector<Coded_info> &neiNode2, double channelInfo, double sideInfo = 0); //for parallel ldgm
  void computeInMeanAndVar(std::vector<Coded_info> &neiNode1, std::vector<Coded_info> &neiNode2, double p0);

  void resize(int size, int group);
  void allocateForMeanAndVar();
  


 public:

  Node_info nodeData;
  double message_to_state;
  double llrEstimation;
  int decision;

};

#endif
