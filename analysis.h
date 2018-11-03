#ifndef ANALYSIS
#define ANALYSIS

#include <TString.h>

using namespace std;

string checkFilename(TString filename);
float CDF(TH1F* hWave,TF1* fTrigFit,float thr);
float CDFinvert(TH1F* hWave,float thr);
float integral(TH1F* hWave,float t1,float t2,float BL);
float* getBL(TH1F* hWave, float* BL, float t1, float t2);
double correction_function(double x);

#endif