/* Generate random matrix with desired degrees, regular or irregular*/
#include <iostream>
#include <map>
#include "stdlib.h"
#include <time.h>

using namespace std;

void scrambleArray(vector<int> &arr, int seed) {
  int len = arr.size();
  
  int randMax = len;
  int exchangePos = len - 1;

  srand(seed);
  while(exchangePos != 0) {
    int randi = rand ()% randMax;
    swap(arr[randi], arr[exchangePos]);
    --randMax;
    --exchangePos;
  }
}

vector<multimap<int, int> > paritycm(int sysNo,int codedNo,int sysDegree, int numWithSysDegree, int seed) {
  
  multimap<int,int> onePositionbyRow;
  multimap<int,int> onePositionbyColumn;

  //sparse matrix parameters
  const int rowNumber=sysNo;
  const int columnNumber=codedNo;

  // Determine constraints for the one number of each row
  vector<int> rowConstraint(numWithSysDegree, sysDegree);
  vector<int> temp(sysNo-numWithSysDegree, sysDegree+1);
  rowConstraint.insert(rowConstraint.end(), temp.begin(), temp.end());
  int scrambleSeed = seed + 1000;
  scrambleArray(rowConstraint, scrambleSeed++);
  
  
  // Determine constraints for the one number of each column
  int totalConnectionNum = numWithSysDegree * sysDegree + (sysNo - numWithSysDegree) * (sysDegree + 1); 
  int codedDegree =  totalConnectionNum / codedNo;
  int numWithCodedDegree = codedNo - totalConnectionNum % codedNo;
  vector<int> columnConstraint(numWithCodedDegree, codedDegree);
  temp = vector<int> (codedNo-numWithCodedDegree, codedDegree+1);
  columnConstraint.insert(columnConstraint.end(), temp.begin(), temp.end());
  scrambleArray(columnConstraint, scrambleSeed++);

  vector<int> columnNumberRec(columnNumber,0);//record of one's number in each column
  vector<int> columnFull;
  vector<int> availColumn(columnNumber,0);
  for(int i=0;i<columnNumber;++i) {
      availColumn[i] = i;
  }
 
   
  //initialize random seeds
  // first generate a random matrix having the right number ones
  if(seed)
    srand(seed);
  else
    srand((unsigned int)time (NULL));
  
  for(int row=0; row<rowNumber; ++row) {
      vector<int> colStoreEachRow;    //store the column number each row
      int colIndex;

      // Get the number of available one positions for the current row
      int colRemain=availColumn.size();

      if(colRemain < rowConstraint[row]) {
				// find the column with smallest number of ones
				int minCol;
	  		minCol=0;
	  		for(int colIndex=1; colIndex<columnNumber; ++colIndex) {
	      	if(onePositionbyColumn.count(colIndex) < onePositionbyColumn.count(minCol))
						minCol=colIndex;
	    	}
	  
	  		// Remove one column from the full list because of the adjustment of one.
	  		int colProvider=columnFull[columnFull.size()-1];
				columnFull.erase(columnFull.end()-1);

	  		//Since a one is removed from the chosen full column.
	  		//that column can be added into the available column list
	  		availColumn.push_back(colProvider);    
	  		--columnNumberRec[colProvider];

				// Move the one from the full column to the column with the smallest number of ones
	 			// Find that for which row the one can be moved from the full column to the column 
	  		// with the smallest number of ones.
	  		pair<multimap<int,int>::iterator, multimap<int,int>::iterator> equal1, equal2;
	  		equal1 = onePositionbyColumn.equal_range(colProvider);
	  		equal2 = onePositionbyColumn.equal_range(minCol);
	  		multimap<int, int>::iterator itDel;
	  		for(multimap<int,int>::iterator itPro = equal1.first;itPro != equal1.second; ++itPro) {
	      	int successFlag=1;
	      
	      	for(multimap<int,int>::iterator itMin = equal2.first; itMin != equal2.second; ++itMin) {
		   			// The minimum column already has one on that row.
	          if((*itPro).second == (*itMin).second) {
            	successFlag = 0;
		       		break;
		     		}
		 			}
	      	// Current row is availabe, so one can be removed from this row in the full column.
	      	if(successFlag == 1) {
						itDel = itPro;
						break;
	      	}
	    	}
	  		onePositionbyColumn.erase(itDel);            //delete the position in col record.

        int rowDeleted= itDel -> second;
	  		equal1 = onePositionbyRow.equal_range(rowDeleted);
	  		itDel = equal1.first;
	  		while(itDel->second != colProvider)  ++itDel;
	    	
	  		onePositionbyRow.erase(itDel);              //delete the position in row record.

	  		// Move the one.
	  		onePositionbyColumn.insert(pair<int,int>(minCol, rowDeleted));
	  		onePositionbyRow.insert(pair<int,int>(rowDeleted,minCol));
	  		++columnNumberRec[minCol];

	  		// Minimum column becomes full
	  		if(columnNumberRec[minCol] == columnConstraint[minCol]) {
	    		columnFull.push_back(minCol);
	    
	    		for(int i=0;i<availColumn.size();++i) {
						if(availColumn[i] == minCol)
		  			availColumn.erase(availColumn.begin()+i);
	      	}
	  		}

	  		cout<<"in row:"<<row<<" adjust one's position!"<<endl;

	  		--row;
	  		continue;
	  
		}

		// Generate all the needed column positions at once for the current row.    
		for(int index=0;index < rowConstraint[row];++index) {
			colIndex = rand() % colRemain;   //random a column postion in available columns
			colStoreEachRow.push_back(availColumn[colIndex]); 
			swap(availColumn[colIndex], availColumn[colRemain-1]);
			--colRemain;
		}
      
    for(int index=0;index < rowConstraint[row];++index) {
			int candidColumn = colStoreEachRow[index];
      
			// Check the column status if this column is chosen
			int onesNumber = columnNumberRec[candidColumn] + 1;
			if(onesNumber <= columnConstraint[candidColumn]) {
				columnNumberRec[candidColumn] = onesNumber;
        if(onesNumber == columnConstraint[candidColumn]) {
					for(int i=0;i<availColumn.size();++i) {
		    		if(availColumn[i] == candidColumn)
		      	availColumn.erase(availColumn.begin()+i);
		  		}
					columnFull.push_back(candidColumn);
	      }
	   	} 
	}//for column iteration
         
    // if the column number is right,store it as a position
    for(int i=0;i<colStoreEachRow.size();++i) {
	  	int column = colStoreEachRow[i];
	  	onePositionbyRow.insert(pair<int, int>(row, column));
	  	onePositionbyColumn.insert(pair<int, int>(column, row));
		}
	} //for row iteration
  
  
  /*check if the number of one in each row and column is correct */
  for(int index=0;index!=rowNumber;++index) {
		if(onePositionbyRow.count(index) != rowConstraint[index]) {
			cout << "One's number is not right at row " << index << endl;
	  	exit(0);
		}
	}
	for(int index=0;index!=columnNumber;++index) {
		if(onePositionbyColumn.count(index) != columnConstraint[index]) {
			cout <<" One's number is not right at column " << index << endl;
	  	exit(0);
		}
	}
  
  cout<<"generate a right random parity check matrix!"<<endl;
  vector<multimap<int, int> > res;
  res.push_back(onePositionbyRow);
  res.push_back(onePositionbyColumn);
  return res;
  /*  
  //-----test part for small matrix----
  int pmatrix[rowNumber][columnNumber];
  for(int i=0;i!=rowNumber;++i)
    for(int j=0;j!=columnNumber;++j)
      pmatrix[i][j]=0;
   
  multimap<int,int>::iterator onesExamer=onepositionbyRow.begin();
  for(;onesExamer!=onepositionbyRow.end();++onesExamer)
    {
      row=onesExamer->first;
      column=onesExamer->second;
      pmatrix[row][column]=1;
    }
   
  for(int i=0;i!=rowNumber;++i)
    {
    for(int j=0;j!=columnNumber;++j)
      {
       cout<<pmatrix[i][j]<<" ";
      }
       cout<<endl;
    }
  return 1;
  */
}
