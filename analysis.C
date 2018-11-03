//root
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>

//C, C++
#include <math.h>
#include <vector>

//specific
#include "analysis.h"

using namespace std;

extern float SP;

/******** FUNCTIONS ********/
string checkFilename(TString filename){
  std::string first ("AB"), second ("CD"), str(filename);
  std::string delimeter("_");

  size_t position = str.rfind(delimeter);
  string ending = str.substr(position+1,2);
  if (ending == second){
    return "WC-Version:1.7";
  }
  else{
    return "WC-Version:1.14";
  }
}

float CDF(TH1F* hWave, TF1* fTrigFit,float thr){
  float peak=hWave->GetMaximum();
  int timePos=1;
  float val = 0;
  while(abs(val)<thr*peak){
    timePos+=1;
    val = hWave->GetBinContent(timePos);
  }
  
  double x1 = SP*(timePos-1);
  double x2 = SP*(timePos);
  double y1 = hWave->GetBinContent(timePos-1);
  double y2 = hWave->GetBinContent(timePos);
  double k = (x2-x1)/(y2-y1);
  return  x1+k*(thr*peak-y1);
}

float CDFinvert(TH1F* hWave,float thr){
  float peak=hWave->GetMaximum();
  int timePos=hWave->GetMaximumBin();
  float val = peak;
  while(val>thr*peak){
    val = hWave->GetBinContent(timePos);
    timePos-=1;
  }
  
  double x1 = SP*(timePos);
  double x2 = SP*(timePos+1);
  double y1 = hWave->GetBinContent(timePos);
  double y2 = hWave->GetBinContent(timePos+1);
  double k = (x2-x1)/(y2-y1);
  return  x1+k*(thr*peak-y1);
  
}

float integral(TH1F* hWave,float t1,float t2,float BL){
  float BW = hWave->GetXaxis()->GetBinWidth(1);
  int bin1 = hWave->FindBin(t1);
  int bin2 = hWave->FindBin(t2);
  float c1 = hWave->GetBinContent(bin1);
  float c2 = hWave->GetBinContent(bin2);
  return hWave->Integral(bin1,bin2,"width")-BL*(t2-t1)-c1*(t1-hWave->GetXaxis()->GetBinLowEdge(bin1))-c2*(hWave->GetXaxis()->GetBinUpEdge(bin2)-t2);
}

float* getBL(TH1F* hWave, float* BL, float t1, float t2){
  /*
  Function to calculate the baseline of the given TH1F-Object.
  Input: TH1F-Object to calculate baseline for; bool isNegative; float-array BL for the output
  Output: baseline and rms of baseline written to 1st and 2nd component of BL-array
  The baseline is calculated as the mean of the values in range of (t1,t2) of the TH1F-Object. The rms
  value is also calculated from the same values.
  
  The float-array BL that is used for the output must be declared before the function call using 'float BL[2];'.
  This is to insure that the output is stored on the heap and not deleted when the memory on the stack is freed up.
   
  Dependencies: function uses C++ vector-class and needs the TMath-header
  */
  
  vector<float> amp;
  for (int i = int(t1/SP); i < int(t2/SP); i++){
    
    amp.push_back(hWave->GetBinContent(i+1));
  }
  BL[0] = TMath::Mean(amp.begin(), amp.end());
  BL[1] = TMath::RMS(amp.begin(), amp.end());
  return BL;
}

double correction_function(double x){
  return (-142.761 + 0.976471* x ) - (-6498.75 + 2.76626*x - 
    0.000141752*TMath::Power(x,2) + 
    4.03526*TMath::Power(10,-9)*TMath::Power(x,3) - 
    6.92814*TMath::Power(10,-14)*TMath::Power(x,4) + 
    7.54846*TMath::Power(10,-19)*TMath::Power(x,5) - 
    5.33119*TMath::Power(10,-24)*TMath::Power(x,6) + 
    2.42945*TMath::Power(10,-29)*TMath::Power(x,7) - 
    6.88509*TMath::Power(10,-35)*TMath::Power(x,8) +
    1.10249*TMath::Power(10,-40)*TMath::Power(x,9) -
    7.61427*TMath::Power(10,-47)*TMath::Power(x,10));
}