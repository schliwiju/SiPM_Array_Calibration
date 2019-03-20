//root
#include <TLine.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <THStack.h>
#include <THistPainter.h>
#include <TText.h>
#include <TSpectrum.h> // peakfinder
#include <TPolyMarker.h> // peakfinder
#include <TError.h> // root verbosity level
//#include <TStyle.h>

//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
//#include <stdlib.h>
//#include <string>
//#include <iomanip>

//specific
#include "geometry.h"
#include "analysis.h"
#include "read.h"

float SP = 0.3125;
float pe = 47.46;//mV*ns
//vector<float> pe_SiPM = {32.14, 39.33, 34.20, 30.79, 34.09, 29.99, 30.69, 29.95}; //a,b,c,d,e,f,g,h  -  Gain-Baseline from fit
vector<float> pe_SiPM = {42.01, 34.67, 34.28, 33.84, 37.55, 34.68, 33.81, 38.84}; //sorted by Wavecatcher-Channel
vector<float> SiPM_shift = {2.679, 2.532, 3.594, 3.855, 3.354, 3.886, 3.865, 4.754};
// vector<float> calib_amp_AB = {10.0024,9.24254,9.08902,10.0149,9.95047,9.55901,10.1483,10.4179,10.0141,9.92513,10.4975,10.422,10.1208,10.1884,10.1682,1};
vector<float> calib_amp_AB = {6.748,6.16313,6.07082,6.68036,6.65783,6.37541,6.7711,6.85418,6.68469,6.58283,6.98329,6.97906,6.76493,6.75924,6.78279,1};
vector<float> calib_amp_AB_max = {9.91652,8.86927,8.88147,9.57771,9.58071,9.14965,9.53239,10.1344,9.62728,9.62879,10.0288,10.3354,9.75948,9.53048,9.68774,1};

/*
__ CALIBRATION FACTORS of indiv. arrays _______________________________
*/
vector<float> calib_amp_8SiPM_17 = {1,3.90,3.35,3.25,3.25,3.60,3.25,3.25,3.55};
vector<float> calib_amp_8SiPM_1 = {1,6.748,6.16313,6.14954,6.07082,6.68036,6.65783,6.37541,6.7711};
vector<float> calib_amp_8SiPM_1_v2 = {1,6.522099,5.992113,5.895799,5.845471,6.190067,6.401603,6.188478,6.546069};

vector<float> calib_amp_8SiPM_2 = {1,6.85418,6.68469,6.58283,6.98329,6.97906,6.76493,6.75924,6.78279};
vector<float> calib_amp_8SiPM_2_v2 = {1,6.813743,6.468777,6.242779,6.607681,6.813801,6.603786,6.476896,6.634614};

vector<float> calib_amp_40SiPM_1_sw1 = {1,5.051830,4.945279,4.925690,4.749162,4.949097,4.723241,4.707982,4.868597};
vector<float> calib_amp_40SiPM_1_sw2 = {1,4.819562,4.740458,4.655669,4.563672,4.697865,4.539159,4.613854,4.748918};
vector<float> calib_amp_40SiPM_1_sw4 = {1,4.682508,4.566132,4.469709,4.511808,4.556756,4.569269,4.513887,4.730196};

vector<float> calib_amp_40SiPM_2_sw1 = {1,4.504408,4.321436,4.349675,4.474370,4.773048,4.893937,4.735639,4.620717};
vector<float> calib_amp_40SiPM_2_sw2 = {1,4.427750,4.289465,4.360091,4.501094,4.655731,4.804050,4.808499,4.657945};
vector<float> calib_amp_40SiPM_2_sw4 = {1,4.396243,4.217127, 4.344094,4.416440,4.678121,4.678319,4.633572,4.705655};

// new procedure:
// 20 ns amplitude acceptance window
// constant baseline
// continuous gain fit
vector<float> g_8SiPM_1_constBL = {1,6.643101,5.937994,5.910624,5.840343,6.142987,6.347942,6.203394,6.606222};
vector<float> g_8SiPM_2_constBL = {1,6.838216,6.477332,6.248899,6.621680,6.813820,6.551395,6.484221,6.685122};
vector<float> g_8SiPM_17_constBL = {1,3.496053,3.189863,3.542560,3.110981,3.210311,3.576601,3.231899,3.257059};

vector<float> g_8SiPM_17_constBL_HV58 = {1,2.472151,2.216387,2.595906,2.191874,2.224226,2.598138,2.246901,2.250155};
vector<float> g_8SiPM_17_constBL_HV59 = {1,2.902316,2.627994,3.016038,2.613367,2.629431,3.043153,2.671950,2.697508};
vector<float> g_8SiPM_17_constBL_HV59_alt = {1,2.888785,2.615228,2.966872,2.555194,2.590216,3.037215,2.637975,2.638252};
vector<float> g_8SiPM_17_dynBL_HV59 = {1,2.937356,2.624025,3.017423,2.613488,2.639697,3.069505,2.650986,2.712761};
vector<float> g_8SiPM_17_constBL_Charge_HV59 = {1,31.8030,28.2290,33.2173,27.3521,27.1600,31.3153,27.4802,27.7432};
vector<float> g_8SiPM_17_constBL_HV60 = {1,3.477412,3.172458,3.543589,3.105156,3.209005,3.558541,3.210445,3.255698};
vector<float> g_8SiPM_17_constBL_HV61 = {1,3.828382,3.523702,3.865647,3.438725,3.538159,3.862610,3.525762,3.600609};

vector<float> g_40SiPM_1_sw4_constBL = {1,4.761221,4.719708,4.592255,4.580244,4.561942,4.508876,4.519476,4.744483};
vector<float> g_40SiPM_1_sw2_Int590_constBL = {1,4.800180,4.712084,4.578888,4.521083,4.783934,4.477283,4.640081,4.785927};
vector<float> g_40SiPM_1_sw1_constBL = {1,5.057210,4.952975,4.914601,4.745170,5.006344,4.695190,4.671079,4.845346};

vector<float> g_40SiPM_2_sw4_constBL = {1,4.367106,4.278874,4.324155,4.502402,4.770954,4.676074,4.684658,4.699855};
vector<float> g_40SiPM_2_sw2_constBL = {1,4.384197,4.294646,4.372352,4.479390,4.659250,4.775097,4.811665,4.681038};
vector<float> g_40SiPM_2_sw1_constBL = {1,4.484410,4.304996,4.311017,4.451629,4.779645,4.873369,4.759235,4.637783};

vector<float> calib_amp_dummy(16,1);

vector<float> calib_amp = g_8SiPM_17_dynBL_HV59;

/*
__ BASELINE SHIFT VALUES of indiv. arrays _______________________________
*/
vector<float> BL_shift_8SiPM_1_run51 = {0,-1.347619,-0.503938,-0.563570,-0.421087,-0.527824,-0.626668,-0.498493,-0.923212};
vector<float> BL_shift_8SiPM_2_run52 = {0,-1.338460,-0.708778,-0.553848,-0.777716,-0.929497,-0.588786,-0.525951,-0.583667};
vector<float> BL_shift_8SiPM_17_run32 = {0,-0.894130,-0.298793,-0.335816,-0.244168,-0.228734,-0.350240,-0.098579,-0.251705};

vector<float> BL_shift_8SiPM_17_run57_HV56 = {0,-0.645652,-0.096661,-0.069783,-0.012360,-0.018906,-0.135682,0.098434,-0.122733};
vector<float> BL_shift_8SiPM_17_run54_HV57 = {0,-0.664145,-0.127923,-0.102607,-0.036950,-0.054989,-0.175302,0.061189,-0.136929};
vector<float> BL_shift_8SiPM_17_run55_HV58 = {0,-0.702267,-0.164472,-0.159815,-0.051495,-0.082261,-0.225175,0.033823,-0.142871};
vector<float> BL_shift_8SiPM_17_run56_HV59 = {0,-0.777759,-0.181139,-0.217059,-0.103472,-0.136753,-0.275943,0.000302,-0.198927};
vector<float> BL_shift_8SiPM_17_run53_HV60 = {0,-0.853458,-0.264629,-0.332783,-0.177919,-0.201986,-0.356012,-0.061629,-0.212546};
vector<float> BL_shift_8SiPM_17_run63_HV61 = {0,-1.055880,-0.352358,-0.540187,-0.297888,-0.362430,-0.560380,-0.210956,-0.348794};

vector<float> BL_shift_40SiPM_1_sw4_run38 = {0,-1.269104,-0.613759,-0.406221,-0.664414,-0.747090,-0.523321,-0.540161,-0.660863};
vector<float> BL_shift_40SiPM_1_sw2_run39 = {0,-1.030745,-0.384704,-0.228097,-0.361703,-0.572680,-0.333203,-0.330794,-0.49887};
vector<float> BL_shift_40SiPM_1_sw2_run41 = {0,-1.009915,-0.391341,-0.184119,-0.341918,-0.508204,-0.340291,-0.292745,-0.489331};
vector<float> BL_shift_40SiPM_1_sw1_run42 = {0,-0.853001,-0.283740,-0.073848,-0.203796,-0.349942,-0.193949,-0.172469,-0.338451};
vector<float> BL_shift_40SiPM_1_sw1_run43 = {0,-0.885457,-0.373082,-0.202642,-0.309442,-0.450317,-0.205849,-0.285526,-0.422478};

vector<float> BL_shift_40SiPM_2_sw4_run44 = {0,-1.000876,-0.454525,-0.354981,-0.455208,-0.715616,-0.494147,-0.480527,-0.688425};
vector<float> BL_shift_40SiPM_2_sw2_run47 = {0,-0.852014,-0.283640,-0.281441,-0.185205,-0.313035,-0.418639,-0.178663,-0.371231};
vector<float> BL_shift_40SiPM_2_sw1_run50 = {0,-0.751363,-0.200432,-0.143382,-0.098068,-0.177394,-0.272392,-0.088074,-0.218430};

vector<float> BL_shift_dummy(9,0);

vector<float> BL_shift = BL_shift_8SiPM_17_run56_HV59;

// clibrated sum
int amp_array_printRate = 1000;

int wavesPrintRate = 1000;
int sumWOMAPrintRate = 1000000;
int sumWOMBPrintRate = 1000000;
int ch0PrintRate = 1000000;
int trigPrintRate = 1000000;//100
int signalPrintRate = 100000;//100
double coef = 2.5 / (4096 * 10);
string WCHU ("AB"), WCAlexander ("CD");

//External Variables - mostly definded in main.C
extern string WCVersion;
extern int runNr;
extern float horizontal;
extern float vertical;
extern float angle;
extern int pdgID;
extern float energy;
extern int isSP;
extern int mp;
extern int safPMT2;
extern int safPMT1;
extern int safSiPM;
extern int trackL;

void read(TString _inFileList, TString _inDataFolder, TString _outFile){
  /*
  Function that is used by main(). read() does most of the analysis of the data. The raw .bin files 
  are read, converted into root-histograms and then the analysis calls are done. Most of the longer
  functions are defined in alanysis.C.
  The read() function the saves all the events (event by event) of that particular run to a root tree
  which is the saved in /runs/runName/out.root.
  */

  /*Formerly used in TB17 analysis when calculation CFD. Not used in TB18 */
  TF1* fTrigFit = new TF1("fTrigFit","gaus");
  fTrigFit->SetParameter(0,800);
  fTrigFit->SetParameter(2,1);
  fTrigFit->SetLineWidth(1);
  
  
  
  /*Create root-file and root-tree for data*/
  TFile *rootFile = new TFile( _outFile, "RECREATE");
  if (rootFile->IsZombie()) {
    cout << "PROBLEM with the initialization of the output ROOT ntuple "
    << _outFile << ": check that the path is correct!!!"
    << endl;
    exit(-1);
  }
  TTree *tree = new TTree("T", "USBWC Data Tree");
  TTree::SetBranchStyle(0);

  /*Declare & define the variables that are to be saved in the root-tree or that are used during the analysis.*/
  Int_t EventNumber=-999;
  Int_t LastEventNumber=-999;
  Float_t SamplingPeriod = -999;
  Double_t EpochTime = -999;
  Int_t Year = -999;
  Int_t Month = -999;
  Int_t Day = -999;
  Int_t Hour = -999;
  Int_t Minute = -999;
  Int_t Second = -999;
  Int_t Millisecond = -999;
  Float_t trigT = -999;//t_trig = (t0+t1+t2+t3)/4
  Float_t tPMT1 = -999;
  Float_t tPMT2 = -999;
  Float_t tPMT2i = -999;
  Float_t tSUMp = -999;
  Float_t tSUMm = -999;
  Float_t trigTp = -999;//t_trig' = [(t0+t1)-(t2+t3)]/4
  Float_t t0t1 = -999;//t0t1 = [(t0-t1)]
  Float_t t2t3 = -999;//t2t3 = [(t2-t3)]
  Int_t isVeto = -999; //variable to define veto, 1 if veto, 0 if not, -999 if undefined
  Int_t isTrig = -999;
  Float_t tsumWOMA_invCFD = -999;
  Float_t tsumWOMB_invCFD = -999;
  Float_t tsumWOMA_invCFD_wrtTrig = -999;
  Float_t tsumWOMB_invCFD_wrtTrig = -999;
  Int_t isLastEvt = -999;
  Int_t isGoodSignal_5 = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  Int_t WOMID[16];  //1=A, 2=B, 3=C, 4=D
  float PE_WOM1, PE_WOM2;

  // clibrated sum
  float amp_array; // clibrated sum 
  Float_t t_amp_array_invCFD = -999;
  Float_t t_amp_array_invCFD_wrtTrig = -999;

  std::vector<float> amp(16,-999);
  std::vector<float> amp_max(16,-999);
  std::vector<float> amp_inRange(16,-999);
  std::vector<float> amp_BL(16,-999);
  std::vector<float> max(16,-999);
  std::vector<float> min(16,-999);
  Float_t t[16];
  Float_t tSiPM[16];

  Float_t BL_lower[16];//store baseline for 16 channels for 0-75ns range
  Float_t BL_RMS_lower[16];//store rms of baseline for 16 channels for 0-75ns range
  Float_t BL_Chi2_lower[16];//store chi2/dof of baseline-fit for 16 channels for 0-75ns range
  Float_t BL_pValue_lower[16];
  Float_t BL_upper[16];//store baseline for 16 channels for 220-320ns range
  Float_t BL_RMS_upper[16];//store rms of baseline for 16 channels for 220-320ns range
  Float_t BL_Chi2_upper[16];//store chi2/dof of baseline-fit for 16 channels for 220-320ns range
  Float_t BL_pValue_upper[16];

  Float_t BL_used[16];
  Float_t BL_Chi2_used[16];
  Float_t BL_pValue_used[16];

  float BL_output[4];//array used for output getBL-function

  int nPeaks = 4; // maximum number of peaks to be stored by peakfinder; has to be set also when creating branch
  Double_t peakX[16][nPeaks];
  Double_t peakY[16][nPeaks];

  float Integral[16];
  float Integral_mVns[16];
  float Integral_0_300[16]; //array used to store Integral of signal from 0 to 300ns
  float Integral_inRange[16];
  int NumberOfBins;
  Int_t EventIDsamIndex[16];
  Int_t FirstCellToPlotsamIndex[16];
  std::vector<TH1F*> hChSum;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChSum_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChSum.push_back(h);
  }
  std::vector<TH1F*> hChShift;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChShift_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChShift.push_back(h);
  }
  std::vector<TH1F> hChtemp;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChtemp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChtemp.push_back(h);
  }
  std::vector<TH1F> hChShift_temp;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChShift_temp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChShift_temp.push_back(h);
  }
  Short_t amplValues[16][1024];
  TH1F hCh("hCh","dummy;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
  TString plotSaveFolder  = _inDataFolder;
  plotSaveFolder.ReplaceAll("data","runs");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(4,4);
  TCanvas csumWOMA("csumWOMA","csumWOMA",1000,700);
  csumWOMA.Divide(4,2);
  TCanvas csumWOMB("csumWOMB","csumWOMB",1000,700);
  csumWOMB.Divide(3,3);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);

  // clibrated sum
  TCanvas C_amp_array("C_amp_array","C_amp_array",1000,700);
  C_amp_array.Divide(3,3);

  /*Create branches in the root-tree for the data.*/
  tree->Branch("EventNumber",&EventNumber, "EventNumber/I");
  tree->Branch("SamplingPeriod", &SamplingPeriod,  "SamplingPeriod/F");
  tree->Branch("EpochTime",&EpochTime, "EpochTime/D");
  tree->Branch("Year",&Year, "Year/I");
  tree->Branch("Month",&Month, "Month/I");
  tree->Branch("Day",&Day, "Day/I");
  tree->Branch("Hour",&Hour, "Hour/I");
  tree->Branch("Minute",&Minute, "Minute/I");
  tree->Branch("Second",&Second, "Second_/I");
  tree->Branch("Millisecond",&Millisecond, "Millisecond/I");
  tree->Branch("trigT",&trigT, "trigT/F");
  tree->Branch("tPMT1",&tPMT1, "tPMT1/F");
  tree->Branch("tPMT2",&tPMT2, "tPMT2/F");
  tree->Branch("tPMT2i",&tPMT2i, "tPMT2i/F");
  tree->Branch("tSUMp",&tSUMp, "tSUMp/F");
  tree->Branch("tSUMm",&tSUMm, "tSUMm/F");
  tree->Branch("runNr",&runNr, "runNr/I");//run number in google table
  tree->Branch("horiz",&horizontal,"horiz/F");// horizontal position of the box units: [cm]
  tree->Branch("vert",&vertical,"vert/F");//vertical position of the box, units: [cm]
  tree->Branch("angle",&angle,"angle/F");
  tree->Branch("pdgID",&pdgID,"pdgID/I");
  // tree->Branch("WOMID",WOMID,"WOMID[nCh]/I");
  tree->Branch("energy",&energy,"energy/F");
  tree->Branch("isSP",&isSP,"isSP/I");
  tree->Branch("mp",&mp,"mp/I");
  tree->Branch("safPMT2",&safPMT2,"safPMT2/I");//solid angle factor
  tree->Branch("safPMT1",&safPMT1,"safPMT1/I");//solid angle factor
  tree->Branch("safSiPM",&safSiPM,"safSiPM/I");//solid angle factor
  tree->Branch("trackL",&trackL,"trackL/I");//track length
  tree->Branch("isLastEvt",&isLastEvt,"isLastEvt/I");
  tree->Branch("trigGate",&trigGate,"trigGate/F");
  tree->Branch("trigTp",&trigTp, "trigTp/F");
  tree->Branch("t0t1",&t0t1, "t0t1/F");//t0t1 = [(t0-t1)]
  tree->Branch("t2t3",&t2t3, "t2t3/F");
  tree->Branch("isVeto",&isVeto,"isVeto/I");
  tree->Branch("isTrig",&isTrig,"isTrig/I");
  tree->Branch("tsumWOMA_invCFD",&tsumWOMA_invCFD,"tsumWOMA_invCFD/F");
  tree->Branch("tsumWOMB_invCFD",&tsumWOMB_invCFD,"tsumWOMB_invCFD/F");
  tree->Branch("tsumWOMA_invCFD_wrtTrig",&tsumWOMA_invCFD_wrtTrig,"tsumWOMA_invCFD_wrtTrig/F");
  tree->Branch("tsumWOMB_invCFD_wrtTrig",&tsumWOMB_invCFD_wrtTrig,"tsumWOMB_invCFD_wrtTrig/F");  
  tree->Branch("isGoodSignal_5",&isGoodSignal_5,"isGoodSignal_5/I");
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("WOMID",WOMID,"WOMID[nCh]/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  tree->Branch("amp",amp.data(), "amp[nCh]/F");
  tree->Branch("amp_max",amp_max.data(), "amp_max[nCh]/F");
  tree->Branch("amp_inRange",amp_inRange.data(), "amp_inRange[nCh]/F");
  tree->Branch("amp_BL",amp_BL.data(), "amp_BL[nCh]/F");
  tree->Branch("max",max.data(), "max[nCh]/F");
  tree->Branch("min",min.data(), "min[nCh]/F");
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("tSiPM", tSiPM, "tSiPM[nCh]/F");
  tree->Branch("PE_WOM1",&PE_WOM1, "PE_WOM1/F");
  tree->Branch("PE_WOM2",&PE_WOM2, "PE_WOM2/F");
  tree->Branch("BL_lower", BL_lower, "BL_lower[nCh]/F");
  tree->Branch("BL_RMS_lower", BL_RMS_lower, "BL_RMS_lower[nCh]/F");
  tree->Branch("BL_Chi2_lower", BL_Chi2_lower, "BL_Chi2_lower[nCh]/F");
  tree->Branch("BL_pValue_lower", BL_pValue_lower, "BL_pValue_lower[nCh]/F");
  tree->Branch("BL_upper", BL_upper, "BL_upper[nCh]/F");
  tree->Branch("BL_RMS_upper", BL_RMS_upper, "BL_RMS_upper[nCh]/F");
  tree->Branch("BL_Chi2_upper", BL_Chi2_upper, "BL_Chi2_upper[nCh]/F");
  tree->Branch("BL_pValue_upper", BL_pValue_upper, "BL_pValue_upper[nCh]/F");
  tree->Branch("BL_used", BL_used, "BL_used[nCh]/F");
  tree->Branch("BL_Chi2_used", BL_Chi2_used, "BL_Chi2_used[nCh]/F");
  tree->Branch("BL_pValue_used", BL_pValue_used, "BL_pValue_used[nCh]/F");
  tree->Branch("peakX",peakX,"peakX[nCh][4]/D");
  tree->Branch("peakY",peakY,"peakY[nCh][4]/D");
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral_inRange", Integral_inRange, "Integral_inRange[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("Integral_mVns", Integral_mVns, "Integral_mVns[nCh]/F");
  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

  // clibrated sum
  tree->Branch("t_amp_array_invCFD",&t_amp_array_invCFD,"t_amp_array_invCFD/F");
  tree->Branch("t_amp_array_invCFD_wrtTrig",&t_amp_array_invCFD_wrtTrig,"t_amp_array_invCFD_wrtTrig/F");
  tree->Branch("amp_array",&amp_array, "amp_array/F");

  /*Start reading the raw data from .bin files.*/
  int nitem = 1;
  ifstream inList;
  TString fileName;
  inList.open(_inFileList);
  assert(inList.is_open());

  // clibrated sum
  int amp_array_PrintStatus=-1;

  int wavePrintStatus=-1;
  int sumWOMAPrintStatus=-1;
  int sumWOMBPrintStatus=-1;
  int ch0PrintStatus=-1;
  int trigPrintStatus=-1;
  int signalPrintStatus=-1;
  while(inList >> fileName){
    fileName = _inDataFolder + fileName;
    cout << endl;
    cout << fileName << endl;
    FILE* pFILE = fopen(fileName.Data(),"rb");
    if (pFILE==NULL) {
      fputs ("File error",stderr); 
      assert(0);
    }
    fseek (pFILE , 0 , SEEK_END);
    int totFileSizeByte = ftell (pFILE);
    rewind (pFILE);
    cout<<"totFileSizeByte = "<<totFileSizeByte<<endl;
    int size_of_header;
    /*During 2018 testbeam measurements two WaveCatchers were used. One from the Berlin group
    and one from Alexander from Geneva. As these two Wavecatchers had two different versions
    there are two types of raw data files that have different header lengths.*/
    if (WCVersion == WCHU){
      size_of_header = 328;
    }
    else if (WCVersion == WCAlexander){
      size_of_header = 327;
    }
    char header[size_of_header];
    nitem=fread(header,1,size_of_header,pFILE);

    cout << "Header:\n" << header << endl;

    char* word;
    word = strtok(header," \n");
    while(word != NULL){
      if(strcmp("ACQUIRED:",word) == 0){
        word = strtok(NULL, " \n");
        nActiveCh = atoi(word);
        break;
      }
      word = strtok(NULL, " \n");
    }

    if(nActiveCh>9){
      cout << endl;
      char dummy;
      nitem=fread(&dummy,1,1,pFILE);
    }

    int whileCounter = 0;
    /*Loop over events. Events are processed and analysed one by one in order.*/
    while(nitem>0){ //event loop
      std::vector<TObject*> eventTrash;
      whileCounter++;
      nitem = fread (&EventNumber, sizeof(int), 1, pFILE);
      nitem = fread (&EpochTime, sizeof(double), 1, pFILE);
      nitem = fread (&Year, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Month, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Day, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Hour, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Minute, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Second, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Millisecond, sizeof(unsigned int), 1, pFILE);
      if (WCVersion == WCHU){
        nitem = fread (&nCh ,sizeof(unsigned int),1,pFILE); // since V2.8.14 the number of stored channels is written for each event
      }
      else if (WCVersion == WCAlexander){
        nCh = 16;
      }
	


      if(EventNumber%100==0){
        printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber,Year,Month,Day,Hour,Minute,Second,Millisecond,nActiveCh,nCh);
      }

      float	MeasuredBaseline[16];
      float	AmplitudeValue[16];
      float	ComputedCharge[16];
      float	RiseTimeInstant[16];
      float	FallTimeInstant[16];
      float	RawTriggerRate[16];
      float floatR=-1;

      /*Loop over individual channels. For each event the data from every channel is 
      processed and analysed one by one in order*/
      for(int i = 0;i<nCh;i++){
        nitem = fread (&ChannelNr[i], sizeof(int), 1, pFILE);
        nitem = fread (&EventIDsamIndex[i], sizeof(int), 1, pFILE);
        nitem = fread (&FirstCellToPlotsamIndex[i],sizeof(int),1,pFILE);
        nitem = fread (&floatR,1,4,pFILE); MeasuredBaseline[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); AmplitudeValue[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); ComputedCharge[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); RiseTimeInstant[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); FallTimeInstant[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); RawTriggerRate[i] = floatR;
        ChannelNr[i]=i;

        /*
        __ Set WOMID _________________________________________________________
        The labeling of the WOMs in the box was done using the letters A,B,C,D. For convinience these letters are here replaced by the numbers 1-4 which is stored in the root-tree for every channel and every event.
        */
        if (WCVersion == WCAlexander){
          if (i <= 6){ WOMID[i] = 3; }
          else if (i >= 7 && i <= 14){ WOMID[i] = 4; }
        }
        else {
          if (i <= 6){ WOMID[i] = 1; }
          else if (i >= 7 && i <= 14){ WOMID[i] = 2; }
        }

        TString title("");
        title.Form("ch %d, ev %d",i,EventNumber);
        hCh.Reset();
        hCh.SetTitle(title);

        /*
        __ Waveform Histogram _______________________________________________
        Writing the signal amplitude values into the root-histogram hCh.
        */
        if (i == 15){
          for(int j = 0;j<1024;j++){
            nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
            hCh.SetBinContent(j+1,-(amplValues[i][j]*coef*1000));
          }
        }
        else {
          for(int j = 0;j<1024;j++){
            nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
            hCh.SetBinContent(j+1,(amplValues[i][j]*coef*1000));
          }
        }

        /*The error of each value in each bin is set to 0.5 mV.*/
        for(int j=1;j<=hCh.GetXaxis()->GetNbins();j++){
          hCh.SetBinError(j,0.5);
        }

        /*Analysis if the event/signal starts.*/
        max[i] = hCh.GetMaximum();
        min[i] = hCh.GetMinimum();


        /*Saving the histogram of that event into a temporary histogram hChtemp. These histograms are available outside of the channel-loop. If analysis using the signals/events of multiple channels needs to be done, this can be accomplished by using hChtemp after the channel-loop.*/
        hChtemp.at(i) = hCh;
	  
        /*
        __ Baseline Fit _______________________________________________________
        Calculate baseline values infront and after the triggered signal
        Triggered signal is expected in the range fromm 100 to 150 ns
        */
        BL_fit(&hChtemp.at(i), BL_output, 0.0, 75.0);
        BL_lower[i] = BL_output[0];
        BL_RMS_lower[i] = BL_output[1];
        BL_Chi2_lower[i] = BL_output[2];
        BL_pValue_lower[i] = BL_output[3];
        BL_fit(&hChtemp.at(i), BL_output, 220.0, 320.0);
        BL_upper[i] = BL_output[0];
        BL_RMS_upper[i] = BL_output[1];
        BL_Chi2_upper[i] = BL_output[2];
        BL_pValue_upper[i] = BL_output[3];

        if (BL_Chi2_upper[i] <= BL_Chi2_lower[i])
        {
          BL_used[i] = BL_upper[i];
          BL_Chi2_used[i] = BL_Chi2_upper[i];
          BL_pValue_used[i] = BL_pValue_upper[i];
        }
        else
        {
          BL_used[i] = BL_lower[i];
          BL_pValue_used[i] = BL_pValue_lower[i];
          BL_Chi2_used[i] = BL_Chi2_lower[i];
        }

        /*
        __ Peakfinder _________________________________________________________
        Implemented to search double-muon-event candiates
        Set maximum number of peaks stored in beginning of script -> nPeaks
        peakX/Yarray[nCh][nPeaks] stores peak coordinates as branches in tree
        Switch on/off with pfON
        -> when off:  set peakX/Yarray[nCh][nPeaks] to zero
        */
        gErrorIgnoreLevel = kError; // suppress root terminal output 

        bool pfON = false;
        if (i<15) {pfON = false;} // switch on/off peakfinder 
        int sigma = 10; // sigma of searched peaks
        Double_t thrPF = 0.1; // peakfinder threshold
        TPolyMarker pm; // store polymarker showing peak position, print later
        peakfinder(&hCh,0,130, nPeaks, sigma, thrPF, peakX[i], peakY[i], &pm, pfON);

        gErrorIgnoreLevel = kUnset; // return to normal terminal output

        // baseline-correct Y-values and convert to units of p.e.
        if (pfON)
        {
          for (int j = 0; j < nPeaks; ++j)
          {
            peakY[i][j] = amp2pe(peakY[i][j], calib_amp[i], BL_used[i]);
          }
        }

        // printf("X: %d %f %f %f %f \n",i,peakX[i][0],peakX[i][1],peakX[i][2],peakX[i][3]);
        // printf("Y: %d %f %f %f %f \n",i,peakY[i][0],peakY[i][1],peakY[i][2],peakY[i][3]);

        /*
        __ Max. Amplitude in Range __________________________________________
        Record maximum amplitude in range before expected signal (100-130 ns)
        */
        // if (i!=2)
        // {
        //   amp_inRange[i] = max_inRange(&hCh,100,150);
        //   // convert p.e. and BL-correct
        //   amp_inRange[i] = amp2pe(amp_inRange[i], calib_amp[i], BL_shift[i]);
        // }
        amp_inRange[i] = max_inRange(&hCh,145,165);
        // convert p.e. and BL-correct
        amp_inRange[i] = amp2pe(amp_inRange[i], calib_amp[i], BL_shift[i]);
        
        Integral_inRange[i] = integral(&hCh,145,170,BL_shift[i]) / g_8SiPM_17_constBL_Charge_HV59[i];
        /*
        __ CFD _____________________________________________________________
        Setting the signal time by using a constant fraction disriminator method.
        The SiPM and the trigger sinals are handled differently using different thresholds.
        */
        if (i == 0){ //trigger
          t[i] = CDF(&hCh,0.5);
        }
        else { //SiPMs
          t[i] = CFD2(&hCh,0.35);
          if (t[i] < 95){
            t[i] = CFDinvert2(&hCh,0.35);
          }
        }



        /*
        __Print Raw Data to .txt ______________________________________________
        Select channel. Prints histogram raw data in a two column text file
        */

        // if (i==4 && BL_chi2[4]<1.7 && BL_chi2[4]>0.7)
        // if (i==4)
        // {
        //   TString histDataName;
        //   histDataName.Form("Ch%d_hist_data.txt",i);
        //   TString path2hist_data;
        //   path2hist_data.Form("%s/%s",(const char*)plotSaveFolder,(const char*)histDataName);
        //   FILE * histOut;
        //   histOut = fopen(path2hist_data,"a"); // produces overhead, maybe put this infront of loop
          
        //   Int_t nbins_x = hCh.GetNbinsX(); // bins at k==0 and k==nbins_x seem to have BinContent==0
        //   for (Int_t k=1; k<=nbins_x; k++)
        //   {
        //     fprintf(histOut,"%.4f %.8f\n",
        //     hCh.GetBinLowEdge(k)+hCh.GetBinWidth(k)/2,
        //     hCh.GetBinContent(k));
        //   }
        //   fclose(histOut);
        // }
      if(EventNumber%sumWOMAPrintRate==0&&i<7){
          csumWOMA.cd(i+1);
          hCh.DrawCopy();
      }
      if(EventNumber%sumWOMBPrintRate==0 && i>6){
         csumWOMB.cd(i-6);
         hCh.DrawCopy();
      }
      // clibrated sum
      // if(EventNumber%amp_array_printRate==0 && (i>0 && i!=8)){
      if(EventNumber%amp_array_printRate==0 && (i>0)){
         C_amp_array.cd(i);
         hCh.DrawCopy();
      }

        /*
        __ Integral & Amplitude ________________________________________
        There are several definitions of the integral of a signal used here. Those are:
        - Integral_0_300: Integration over the entire time window (~320ns)
        - Integral: Integration over a smaller time window (~50ns) relative to the trigger
        Additionally the number of p.e. is now calculated using the amplitude
        and the calibration factors in the calib_amp-vactor. The function 'PE' calculates the amplitude of the signal, subtracts the better BL value and divides by the calibration factor.
        */
        Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-0.0*1024*SP);
        if (BL_Chi2_upper[i] <= BL_Chi2_lower[i]){
        	Integral[i] = Integrate_50ns(&hCh, BL_upper[i]);
        	amp[i] = PE(&hCh,calib_amp.at(i),BL_upper[i], 180.0, 230.0);
          // amp[i] = PE(&hCh,calib_amp_AB.at(i),BL_upper[i], t[i]-20, t[i]+30);
          // BL_used[i] = BL_upper[i];
          // BL_Chi2_used[i] = BL_Chi2_upper[i];
          // BL_pValue_used[i] = BL_pValue_upper[i];
          hCh.GetXaxis()->SetRange((t[i]-20)/SP,(t[i]+30)/SP);
          int max_bin = hCh.GetMaximumBin();
          amp_max[i] = (hCh.GetBinContent(max_bin)-BL_upper[i])/calib_amp.at(i);
          hCh.GetXaxis()->SetRange(0,1024);
        }
        else{
        	Integral[i] = Integrate_50ns(&hCh, BL_lower[i]);
         	amp[i] = PE(&hCh,calib_amp.at(i),BL_lower[i], 180.0, 230.0);
          // amp[i] = PE(&hCh,calib_amp_AB.at(i),BL_lower[i], t[i]-20, t[i]+30);
          // BL_used[i] = BL_lower[i];
          // BL_Chi2_used[i] = BL_Chi2_lower[i];
          // BL_pValue_used[i] = BL_pValue_lower[i];
          hCh.GetXaxis()->SetRange((t[i]-20)/SP,(t[i]+30)/SP);
          int max_bin = hCh.GetMaximumBin();
          amp_max[i] = (hCh.GetBinContent(max_bin)-BL_lower[i])/calib_amp.at(i);
          hCh.GetXaxis()->SetRange(0,1024);
        }

        hCh.GetXaxis()->SetRange(0./SP,75./SP);
        int max_bin = hCh.GetMaximumBin();
        amp_BL[i] = hCh.GetBinContent(max_bin);
        hCh.GetXaxis()->SetRange(0,1024);
      // End of loop over inividual channels

        /*
        __ Printing Wafevorms ____________________________________________
        The signals for events can be printed to a .pdf file called waves.pdf. The rate at which the events are drawn to waves.pdf is set via the variable wavesPrintRate. Additional requirements can be set in the if-statement to look at specific events only.
        The entire if-statement so far also plots lines at the found signal maximum, the corresponding integration limit, as well as the BL values to each of the histograms.
        */
        if(EventNumber%wavesPrintRate==0){
          cWaves.cd(1+4*(i%4)+(i)/4);
          hCh.DrawCopy();
          hCh.GetXaxis()->SetRange((t[i]-20)/SP,(t[i]+30)/SP);
          int max_bin = hCh.GetMaximumBin();
          int lower_bin = max_bin - 20.0/SP;
          int upper_bin = max_bin + 30.0/SP;
          // double x = h->GetXaxis()->GetBinCenter(binmax);
          float max_time = hCh.GetXaxis()->GetBinCenter(max_bin);
          float lower_time = hCh.GetXaxis()->GetBinCenter(lower_bin);
          float upper_time = hCh.GetXaxis()->GetBinCenter(upper_bin);
          hCh.GetXaxis()->SetRange(0,1024);
          TLine* ln = new TLine(max_time,-2000,max_time,2000);
          TLine* ln2 = new TLine(lower_time,-2000,lower_time,2000);
          TLine* ln3 = new TLine(upper_time,-2000,upper_time,2000);
          TLine* ln4 = new TLine(0,BL_lower[i],75,BL_lower[i]);
          TLine* ln5 = new TLine(220,BL_upper[i],320,BL_upper[i]);
          TText *text = new TText(.5,.5,Form("%f %f",BL_lower[i],BL_upper[i]));
          TText *text2 = new TText(max_time+3,amp[i]+1,Form("%f %f %f",amp[i], amp_max[i], BL_used[i]));
          ln->SetLineColor(2);
          ln2->SetLineColor(3);
          ln3->SetLineColor(3);
          ln4->SetLineColor(2);
          ln5->SetLineColor(2);
          ln->Draw("same");
          ln2->Draw("same");
          ln3->Draw("same");
          ln4->Draw("same");
          ln5->Draw("same");
          text->Draw("same");
          text2->Draw("same");
          if (pfON){pm.Draw();} // print peakfinders polymarker
        }
      }

      /*
      __ Number of P.E. _____________________________________________________
      Calculate & save the number of p.e. for an entire WOM. Note: for the 1st WOM of each WC one channel was not recorded. Thus, there are only 7 values from 8. The result for the WOM is therefore scaled up by 8/7 to make the numbers comparable.
      */
      PE_WOM1 = 8/7*(amp[0]+amp[1]+amp[2]+amp[3]+amp[4]+amp[5]+amp[6]);
      PE_WOM2 = (amp[7]+amp[8]+amp[9]+amp[10]+amp[11]+amp[12]+amp[13]+amp[14]);
      trigT = t[15];
      for (int i=0; i<=14; i++){
        tSiPM[i] = t[i] - trigT;
        /*
        if (tSiPM[i+7] < -66){
          t[i+7] = CDFinvert(&hChtemp.at(i+7),0.33);
          tSiPM[i+7] = t[i+7] - trigT;
        }
        */
      }

      /*
      __ Calibrated SUM SIGNAL of all channels ______________
      */

      // clibrated sum
      TH1F hSum("hSum","Calibrated Sum ;ns;Amplitude in npe",1024,-0.5*SP,1023.5*SP);
      for(int hSumIndex=1;hSumIndex<9;hSumIndex++){
        // if (hSumIndex!=2)
        // {
        //   TF1* f_const = new TF1("f_const","pol0",0,320);
        // f_const->SetParameter(0,BL_shift[hSumIndex]);

        // hChtemp.at(hSumIndex).Add(f_const, -1);
        // hChtemp.at(hSumIndex).Scale(1.0/calib_amp.at(hSumIndex));

        // hSum.Add(&hChtemp.at(hSumIndex),1);
        // }
        TF1* f_const = new TF1("f_const","pol0",0,320);
        f_const->SetParameter(0,BL_shift[hSumIndex]);

        hChtemp.at(hSumIndex).Add(f_const, -1);
        hChtemp.at(hSumIndex).Scale(1.0/calib_amp.at(hSumIndex));

        hSum.Add(&hChtemp.at(hSumIndex),1);
        
      }

      amp_array = max_inRange(&hSum, 130, 180.0);
      // amp_array = 8.0/7.0 * max_inRange(&hSum, 100.0, 150.0);
      C_amp_array.cd(9);
      hSum.DrawCopy();

      t_amp_array_invCFD = CFDinvert2(&hSum,0.4);
      t_amp_array_invCFD_wrtTrig = trigT-t_amp_array_invCFD;

      /* Filling Sum Histograms for WOM A and B and determine time Resolution */
      /*
      TH1F hSumA("hSumA","Sum A;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
      for(int hSumIndexA=0;hSumIndexA<7;hSumIndexA++){
        hSumA.Add(&hChtemp.at(hSumIndexA),1);
      }
      csumWOMA.cd(8);
      hSumA.DrawCopy();

      tsumWOMA_invCFD = CFDinvert2(&hSumA,0.4);
      tsumWOMA_invCFD_wrtTrig = trigT-tsumWOMA_invCFD;
*/
      // TH1F hSumB("hSumB","Sum B;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
      // for(int hSumIndexB=7;hSumIndexB<15;hSumIndexB++){
      //   TF1* f_const = new TF1("f_const","pol0",0,320);
      //   f_const->SetParameter(0,BL_used[hSumIndexB]);

      //   hChtemp.at(hSumIndexB).Add(f_const, -1);
      //   hChtemp.at(hSumIndexB).Scale(1.0/calib_amp.at(hSumIndexB));

      //   hSumB.Add(&hChtemp.at(hSumIndexB),1);
      // }

      // PE_WOM2 = PE(&hSumB,1,0, 100.0, 150.0);
      // csumWOMB.cd(9);
      // hSumB.DrawCopy();

      // tsumWOMB_invCFD = CFDinvert2(&hSumB,0.4);
      // tsumWOMB_invCFD_wrtTrig = trigT-tsumWOMB_invCFD;
/*
      */
      /* end */

      /*temporary for cfdScan_
      if (PE_WOM2 >= 5){    //Switch PE_WOM of you want WOM 1 or 2
        for(int i=0;i<N_CFD_points;i++){
          v_hTimeRes.at(i)->Fill(trigT-CFD2(&hSumB,0.01*i)); // Switch hSumA or B or CFD2 or CFDinvert2
        }
      }
      *///temporary end
      
      /*
      trigGate = abs(*(std::max_element(t,t+4))-*(std::min_element(t,t+4)));  
      
      if(max[0]<1240&&max[1]<1240&&max[2]<1240&&max[3]<1240&&isVeto==0){
        isTrig=1;
        if(isTrig&&BL[0]<1.1&&BL[1]<1.1&&BL[2]<1.1&&BL[3]<1.1){
          isTrig=1;
          if(trigT<140&&trigT>90&&trigGate<10){
            isTrig=1;
          }
          else isTrig=0;
        }
        else isTrig=0;
      }
      else isTrig=0;
      if(isTrig==1){
        int shift = (int)((140-trigT)/SP);
        for(int j=0;j<(int)hChtemp.size();j++){
          hChSum.at(j)->Add(&hChtemp.at(j),1);
          hChShift_temp.at(j).Reset();
          for(int bin=1;bin<=hCh.GetXaxis()->GetNbins()-shift;bin++){
            hChShift_temp.at(j).SetBinContent(shift+bin,hChtemp.at(j).GetBinContent(bin));
          }
          hChShift.at(j)->Add(&hChShift_temp.at(j),1);
        }
      }
      
      if(isTrig==1&&max[5]<1240)isGoodSignal_5=1;
      else isGoodSignal_5=0;
      */

      /*Saving the plotted signals/events to a new page in the .pdf file.*/
      if(EventNumber%wavesPrintRate==0) {
        if(wavePrintStatus<0){
          cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");
          wavePrintStatus=0;
        }
        else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
      }
      if(EventNumber%sumWOMAPrintRate==0){
        if(sumWOMAPrintStatus<0){
          csumWOMA.Print((TString)(plotSaveFolder+"/sumWOMA.pdf("),"pdf");
          sumWOMAPrintStatus=0;
        }
        else csumWOMA.Print((TString)(plotSaveFolder+"/sumWOMA.pdf"),"pdf");
      }
      if(EventNumber%sumWOMBPrintRate==0){
        if(sumWOMBPrintStatus<0){
           csumWOMB.Print((TString)(plotSaveFolder+"/sumWOMB.pdf("),"pdf");
           sumWOMBPrintStatus=0;
        }
        else csumWOMB.Print((TString)(plotSaveFolder+"/sumWOMB.pdf"),"pdf");
      }
      if(EventNumber%trigPrintRate==0){
        if(trigPrintStatus<0){
          cTrig.Print((TString)(plotSaveFolder+"/trig.pdf("),"pdf");
          trigPrintStatus=0;
        }
        else cTrig.Print((TString)(plotSaveFolder+"/trig.pdf"),"pdf");
      }
      if(EventNumber%signalPrintRate==0){
        if(signalPrintStatus<0){
          cSignal.Print((TString)(plotSaveFolder+"/signal.pdf("),"pdf");
          signalPrintStatus=0;
        }
        else cSignal.Print((TString)(plotSaveFolder+"/signal.pdf"),"pdf");
      }

      // clibrated sum
      if(EventNumber%amp_array_printRate==0){
        if(amp_array_PrintStatus<0){
           C_amp_array.Print((TString)(plotSaveFolder+"/amp_array.pdf("),"pdf");
           amp_array_PrintStatus=0;
        }
        else C_amp_array.Print((TString)(plotSaveFolder+"/amp_array.pdf"),"pdf");
      }

      /*Writing the data for that event to the tree.*/
      tree->Fill();
    }
    fclose(pFILE);
  }

  /*Clearing objects and saving files.*/
  inList.close();
  cWaves.Clear();
  cWaves.Print((TString)(plotSaveFolder+"/waves.pdf)"),"pdf");
  csumWOMA.Clear();
  csumWOMA.Print((TString)(plotSaveFolder+"/sumWOMA.pdf)"),"pdf");
  csumWOMB.Clear();
  csumWOMB.Print((TString)(plotSaveFolder+"/sumWOMB.pdf)"),"pdf");
  cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf)"),"pdf");
  cTrig.Print((TString)(plotSaveFolder+"/trig.pdf)"),"pdf");
  cSignal.Print((TString)(plotSaveFolder+"/signal.pdf)"),"pdf");

  // clibrated sum
  C_amp_array.Clear();
  C_amp_array.Print((TString)(plotSaveFolder+"/amp_array.pdf)"),"pdf");

  /* temporary for cfdScan
  TH1F* sigma_timeRes_Fit = new TH1F("sigma_timeRes_Fit",";CFD threshold;Timeresolution, ns",N_CFD_points,0,0.01*N_CFD_points);
  TF1* fGaus = new TF1("fGaus","gaus",45,60);
  fGaus->SetLineWidth(1);
  TCanvas c_TimeRes("c_TimeRes");
  for(int i=0;i<N_CFD_points;i++){
    c_TimeRes.cd();
    v_hTimeRes.at(i)->Draw();
    v_hTimeRes.at(i)->Fit("fGaus","R");
    sigma_timeRes_Fit->SetBinContent(i+1,fGaus->GetParameter(2));
    sigma_timeRes_Fit->SetBinError(i+1,fGaus->GetParError(2));
    TString name("");
    name.Form("hTimeRes%d.png",i);
    c_TimeRes.Print((TString)(plotSaveFolder+"/"+name));
  }
  c_TimeRes.cd();
  sigma_timeRes_Fit->SetMarkerSize(1);
  sigma_timeRes_Fit->SetMarkerStyle(20);
  sigma_timeRes_Fit->SetLineColor(kRed);
  //c_TimeRes.Print((TString)("/home/maximilian/Dokumente/Arbeit/CERNTestBeam2018/Auswertung/CERN-TestBeam-2018/cfd-scan/"+name));
  c_TimeRes.Write();
  */// temporary end

  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}