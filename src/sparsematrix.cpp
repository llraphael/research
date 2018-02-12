#include<iostream>
#include<vector>
#include"sparsematrix.h"


using namespace std;

ElementData::ElementData(int row, int column, int initialvalue):
  row_pos(row),
  column_pos(column),
  value(initialvalue)
{}

void ElementData::setValue(int row, int column, int setValue)
{

  row_pos = row;
  column_pos = column;
  value = setValue;

}



Matrix::Matrix(int row, int column):
  matrix_data(row),
  row_num(row),
  column_num(column)
{
}

int Matrix::readData(const ElementData data)
{

  if(data.row_pos > row_num || data.column_pos > column_num) {
    cout<<"input position is beyond matrix range! Reading fails!"<<endl;
    return -1;
  }

  int row_pos = data.row_pos;
  int column_pos = data.column_pos;

  int size = matrix_data[row_pos].colPos.size();

  if(size == 0)
    matrix_data[row_pos].colPos.push_back(data);
  else {
    for(int i=0;i<size;++i)
      {
	if(column_pos > matrix_data[row_pos].colPos[i].column_pos) {
	  if(i == size-1)    //End of the column position
	    matrix_data[row_pos].colPos.push_back(data);
	  else 
	    continue;
	}
	else if(column_pos == matrix_data[row_pos].colPos[i].column_pos) {
	  matrix_data[row_pos].colPos[i] = data;
	  break;
	}
	else {
	  matrix_data[row_pos].colPos.insert(matrix_data[row_pos].colPos.begin()+i, data);
	  break;
	}
      }
  }

}


int Matrix::searchPos(const int row, const int column)
{

  if(row > row_num || column> column_num) {
    cout<<"input position is beyond matrix range! Reading fails!"<<endl;
    return -1;
  }

  int size = matrix_data[row].colPos.size();
  int index = -1;
  for(int i=0;i<size;++i)
    {
      if(matrix_data[row].colPos[i].column_pos == column) {
	index = i;
	break;
      }
    }

  return index;
}

int Matrix::posValue(int row, int column)
{

  int index = searchPos(row, column);

  if(index>=0) 
    return matrix_data[row].colPos[index].value;
  else 
    return 0;

}

int Matrix::columnPosbyIndex(const int row, const int index)
{

  if(row > row_num || index > matrix_data[row].colPos.size() ) {
    cout<<"input position is beyond matrix range! Reading fails!"<<endl;
    return -1;
  }

  return matrix_data[row].colPos[index].column_pos;
}





int Matrix::matrixStackRight(const Matrix &right)
{
  if(row_num != right.row_num) {
    cout<<"Two matrices do not have the same row number. Cannot stack horizontally!"<<endl;
    return -1;
  }

  int origin_column_num = column_num;
  column_num = column_num + right.column_num;          // column num of the new matrix

  for(int index=0;index<row_num;++index)
    {
      int ele_num = right.matrix_data[index].colPos.size();
      for(int ele_index=0;ele_index<ele_num;++ele_index)
	{
	  ElementData ele(right.matrix_data[index].colPos[ele_index]);
	  int ele_col = ele.column_pos;
	  int ele_row = ele.row_pos;
	  int value = ele.value;
	  ele.setValue(ele_row, ele_col+origin_column_num, value);

	  readData(ele);
	}
    }

}

int Matrix::matrixStackDown(const Matrix &down)
{

  if(column_num != down.column_num) {
    cout<<"Two matrices do not have the same column number. Cannot stack vertically"<<endl;
    return -1;
  }

  int origin_row_num = row_num;
  row_num += down.row_num;

  for(int i=0;i<down.row_num;++i)
    {
      RowData eachrow(down.matrix_data[i]);
      for(int j=0;j<eachrow.colPos.size();++j)
	{
	  int origrow = eachrow.colPos[j].row_pos;
	  eachrow.colPos[j].row_pos += origin_row_num;
	}
      matrix_data.push_back(eachrow);
    }
}

int Matrix::organizebyCol()
{
  for(int i=0;i<row_num;++i)
    {
      int size = matrix_data[i].colPos.size();
      for(int j=0;j<size;++j)
	{
	  int column = matrix_data[i].colPos[j].column_pos;
	  int row = matrix_data[i].colPos[j].row_pos;
	  matrixbycolumn.insert( std::pair<int, int>(column, row) );
	}
    }
}


Matrix Matrix::columnPermutation(double seed)
{

  vector<int> oldpos(column_num);
  for(int i=0;i<column_num;++i)
    oldpos[i] = i;

  vector<int> newpos(column_num, 0);              // New number for each column.
 
  srand(seed);

  for(int i=0;i<column_num;++i)
    {
      int ran = rand ()% (column_num-i);
      newpos[i] = oldpos[ran];

      swap(oldpos[ran], oldpos[column_num-i-1]);
    }

  for(int rowindex=0;rowindex<row_num;++rowindex)
    {
      int rowelenum = matrix_data[rowindex].colPos.size();
      for(int colindex=0;colindex<rowelenum;++colindex)
	{
	  int orirow = matrix_data[rowindex].colPos[colindex].column_pos;

	  matrix_data[rowindex].colPos[colindex].column_pos = newpos[orirow];
	}
    }
  return *this;

}
