#ifndef ANALYSIS
#define ANALYSIS

#include <TString.h>
#include <TPolyMarker.h>

using namespace std;

string checkFilename(TString filename);
float CDF(TH1F* hWave,float thr);
float CFD2(TH1F* hWave,float thr);
float CDFinvert(TH1F* hWave,float thr);
float CFDinvert2(TH1F* hWave,float thr);
float Integrate_50ns(TH1F* hWave, float BL);
float integral(TH1F* hWave,float t1,float t2,float BL);
float* getBL(TH1F* hWave, float* BL, float t1, float t2);
float* BL_fit(TH1F* hWave, float* BL_chi2, float t1, float t2);
float PE(TH1F* hWave, float calib_factor, float BL, float t1, float t2);
float max_inRange(TH1F* hWave,float t1, float t2);
double amp2pe(double y, float calib_factor, float BL_upper, float BL_lower, float BL_Chi2_upper, float BL_Chi2_lower);
double correction_function(double x);
void peakfinder(TH1F *hWave, float t1, float t2, int nPeaks, int sigma, double thr, double *Xarray, double *Yarray, TPolyMarker *pfMarker, bool pfON);

#endif