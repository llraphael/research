#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include<iostream>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<map>
#include<cmath>

class ElementData {

 public:
  
   ElementData(): row_pos(0), column_pos(0), value(0) {}; 
   ElementData(int row, int column, int initialvalue);

 public:

  int row_pos;
  int column_pos;
  int value;

 public:

  void setValue(int row, int column, int setValue);

};

class Matrix {

 public:

  struct RowData {

    std::vector<ElementData> colPos;
  };

 public: 

  Matrix(): row_num(0), column_num(0){};
  Matrix(int row, int column);

 public:

  int readData(const ElementData data);
  int searchPos(const int row, const int column);
  int posValue(const int row, const int column);
  int columnPosbyIndex(const int row, const int index);

 public:
  
  int matrixStackRight(const Matrix &right);
  int matrixStackDown(const Matrix &down);

 public:
  
  int organizebyCol();
  Matrix columnPermutation(double seed);

 public:

  std::vector<RowData> matrix_data;
  std::multimap<int, int> matrixbycolumn;
  int row_num;
  int column_num;

};

#endif
