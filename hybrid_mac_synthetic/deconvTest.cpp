#include <iostream>
#include "func.h"

using namespace std;

int main() {

	vector<double> pdf1 = {0, 0, 0.1, 0.5, 0.4};
	vector<double> pdf2 = {0, 0, 0.3, 0.5, 0.2};
	vector<double> pdf3 = {0, 0, 0.4, 0.5, 0.1};

	vector<double> sumPdf(1, 1);
	
	conv(pdf1, pdf2, sumPdf);
	conv(sumPdf, pdf3, sumPdf);
	vecNorm(sumPdf);

	cout << "[";
	for(int i=0; i<sumPdf.size(); i++)
		cout << sumPdf[i] << ", ";
	cout << "]" << endl;

}
