#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>


using namespace std;



void awgn(vector<double>::iterator outputBegin, vector<double>::iterator outputEnd, double sigma)
{
	
	vector<double>::iterator begin=outputBegin;
	vector<double>::iterator end=outputEnd;
	vector<double>::iterator index;

        double N0=2*pow(sigma,2);

	int range=10000;
	double u1,u2,v1,v2,s;
	double x,y;

	

	

	//srand(time(NULL));
	srand(5000);

	for(index=begin;index!=end;++index)
	{

	do{
              
	  u1=double(rand())/RAND_MAX;
	  // u1/=RAND_MAX;
	  u2=double(rand())/RAND_MAX;
	  //	u2/=RAND_MAX;
		v1=2*u1-1;
		v2=2*u2-1;
		s=v1*v1+v2*v2;
	}while(s>=1);

	x=sqrt((-2)*log(s)/s)*v1;
	y=sqrt((-2)*log(s)/s)*v2;

	x=sqrt(N0/2)*x;
	y=sqrt(N0/2)*y;
    
	*index=*index+x;
	}

}

