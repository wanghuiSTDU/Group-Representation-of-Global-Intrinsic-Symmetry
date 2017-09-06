/*=================================================================
* Generating the consistent pairs with intrinsic distances

*
* This is a MEX-file for MATLAB.  
* Hui Wang, Feb. 28, 2017, wanghui19841109@163.com 
*
=================================================================*/

#include <mex.h>
#include <math.h>
#include <vector>

using namespace std;
typedef vector<int> Vector_int;
typedef vector<Vector_int> Vector_int_int;

bool isConsistent(Vector_int P1, Vector_int P2, double *distance, int dimOfDis,double deta)
{
   bool flag;
   int i,j,n = P1.size();
   
   double min = 1e+5;
   double radio;
   for(i = 0; i < n; i++)
       for(j = 0; j < n; j++)
       {
           int in1i = P1.at(i);
           int in1j = P1.at(j);
           double d1 = distance[in1i * dimOfDis + in1j];
           
           int in2i = P2.at(i);
           int in2j = P2.at(j);
           double d2 = distance[in2i * dimOfDis + in2j];
           
           double radio1 = d1 / d2;
           double radio2 = d2 / d1;
           if(radio1 < radio2)
               radio = radio1;
           else
               radio = radio2;
           
           if(radio < min)
               min = radio;
       }
   
   if(min < deta)
       return false;
   else
       return true;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	int i,j,k,n1,n2,numOfPoints,dimOfDis;
	double *p1,*p2,*distance,*consistentPairs1,*consistentPairs2,*deta;
    
	n1 = mxGetM(prhs[0]);
	numOfPoints = mxGetN(prhs[0]);
	p1 = mxGetPr(prhs[0]);

    n2 = mxGetM(prhs[1]);
	p2 = mxGetPr(prhs[1]);
    
	distance = mxGetPr(prhs[2]);
    dimOfDis = mxGetM(prhs[2]);
    
    deta = mxGetPr(prhs[3]);
    
    //Process
    Vector_int_int Pairs1,Pairs2;
    int numberOfConsisPairs = 0;
	for(i = 0; i < n1; i++)
        for(j = 0; j < n2; j++)
        {  
            Vector_int pairs1,pairs2;
            pairs1.resize(numOfPoints);
            pairs2.resize(numOfPoints);
            for(k = 0; k < numOfPoints; k++)
            {
               pairs1[k] = p1[k * n1 + i] - 1;
               pairs2[k] = p2[k * n2 + j] - 1;
            }
         
           if(isConsistent(pairs1, pairs2, distance, dimOfDis, deta[0]) == true)
           {
               Pairs1.push_back(pairs1);
               Pairs2.push_back(pairs2);
               numberOfConsisPairs++;
           }
           
          pairs1.clear();
          pairs2.clear();
        }
    
    //Output data
    plhs[0] = mxCreateDoubleMatrix(numberOfConsisPairs, numOfPoints,mxREAL);
    consistentPairs1 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(numberOfConsisPairs, numOfPoints,mxREAL);
    consistentPairs2 = mxGetPr(plhs[1]);
    
    for(i = 0; i < numberOfConsisPairs; i++)
        for(k = 0; k < numOfPoints; k++)
        {
            consistentPairs1[k * numberOfConsisPairs + i] = Pairs1.at(i).at(k) + 1;
            consistentPairs2[k * numberOfConsisPairs + i] = Pairs2.at(i).at(k) + 1;
        }
    
   //Clear data
   for(i = 0; i < numberOfConsisPairs; i++)
   {
        Pairs1.at(i).clear();
        Pairs2.at(i).clear();
   }
   Pairs1.clear();
   Pairs2.clear();
}

