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
int wavesPrintRate = 100;
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
  Int_t isLastEvt = -999;
  Int_t isGoodSignal_5 = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  Int_t WOMID[16];  //1=A, 2=B, 3=C, 4=D
  std::vector<float> amp(16,-999);
  std::vector<float> max(16,-999);
  std::vector<float> min(16,-999);
  Float_t t[16];
  Float_t tSiPM[16];
  Float_t BL[16];//store baseline for 16 channels
  Float_t BL_RMS[16];//store rms of baseline for 16 channels
  float BL_output[2];//array used for output getBL-function
  float Integral_0_300[16];//array used to store Integral of signal from 0 to 300ns
  float Integral[16];
  float Integral_Correction;
  float Integral_mVns[16];
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
  TString plotSaveFolder  = _outFile;
  plotSaveFolder.ReplaceAll("out.root","");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(4,4);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);

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
  // tree->Branch("WOMID",&WOMID,"WOMID/I");
  tree->Branch("WOMID",WOMID, "WOMID[nCh]/I");
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
  tree->Branch("isGoodSignal_5",&isGoodSignal_5,"isGoodSignal_5/I");
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  tree->Branch("amp",amp.data(), "amp[nCh]/F");
  tree->Branch("max",max.data(), "max[nCh]/F");
  tree->Branch("min",min.data(), "min[nCh]/F");
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("tSiPM", tSiPM, "tSiPM[nCh]/F");
  tree->Branch("BL", BL, "BL[nCh]/F");
  tree->Branch("BL_RMS", BL_RMS, "BL_RMS[nCh]/F");
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("Integral_Correction",&Integral_Correction, "Integral_Correction/F");
  tree->Branch("Integral_mVns", Integral_mVns, "Integral_mVns[nCh]/F");
  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

  /*Start reading the raw data from .bin files.*/
  int nitem = 1;
  ifstream inList;
  TString fileName;
  inList.open(_inFileList);
  assert(inList.is_open());

  int wavePrintStatus=-1;
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

        /*The labeling of the WOMs in the box was done using the letters A,B,C,D. For convinience
        these letters are here replaced by the numbers 1-4 which is stored in the root-tree for
        every channel and every event.*/
        if (WCVersion == WCAlexander){
          if (i <= 6){ WOMID[i] = 1; }
          else if (i >= 7 && i <= 14){ WOMID[i] = 2; }
        }
        else {
          if (i <= 6){ WOMID[i] = 3; }
          else if (i >= 7 && i <= 14){ WOMID[i] = 4; }
        }

        TString title("");
        title.Form("ch %d, ev %d",i,EventNumber);
        hCh.Reset();
        hCh.SetTitle(title);

        /*Writing the signal amplitude values into the root-histogram hCh.*/
        for(int j = 0;j<1024;j++){
          nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
          hCh.SetBinContent(j+1,(amplValues[i][j]*coef*1000));
        }

        /*Analysis if the event/signal starts.*/
        max[i] = hCh.GetMaximum();
        min[i] = hCh.GetMinimum();
        amp[i] = hCh.GetMaximum();
	  
        getBL(&hCh, BL_output,0,30);
        BL[i] = BL_output[0];
        BL_RMS[i] = BL_output[1];

        amp[i] = hCh.GetMaximum() - BL[i];
	  
        /*The error of each value in each bin is here set to the root-mean-square of the
        baseline.*/
        for(int j = 1;j <= hCh.GetXaxis()->GetNbins();j++){
          hCh.SetBinError(j,BL_RMS[i]);
        }

        /*Saving the histogram of that event into a temporary histogram hChtemp.
        These histograms are available outside of the channel-loop. If analysis using
        the signals/events of multiple channels needs to be done, this can be accomplished
        by using hChtemp after the channel-loop.*/
        hChtemp.at(i) = hCh;

        /*Setting the signal time by using a constant fraction disriminator method.
        The SiPM and the trigger sinals are handled differently using different thresholds.*/
        if (i == 15){ //trigger
          t[i] = CDF(&hCh,fTrigFit,0.5);
        }
        else { //SiPMs
          t[i] = CDF(&hCh,fTrigFit,0.1);
        }

        /*The signals for events can be printed to a .pdf file called waves.pdf. The rate at
        which the events are drawn to waves.pdf is set via the variable wavesPrintRate. Additional
        requirements can be set in the if-statement to look at specific events only.*/
        if(EventNumber%wavesPrintRate==0){
          cWaves.cd(1+4*(i%4)+(i)/4);
          hCh.DrawCopy();
          TLine* ln = new TLine(t[i],-2000,t[i],2000); //draw red vertical line at signal time
          ln->SetLineColor(2);
          ln->Draw("same");
        }

        /*There are several definitions of the integral of a signal used here. Those are:
        - Integral_0_300: Integration over the entire time window (~320ns)
        - Integral: Integration over a smaller time window (~50ns) relative to the trigger*/
        Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-BL[i]*1024*SP);

      }

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
      if(EventNumber%wavesPrintRate==0){
        if(wavePrintStatus<0){
          cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");
          wavePrintStatus=0;
        }
        else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
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

      /*Writing the data for that event to the tree.*/
      tree->Fill();
    }
  }

  /*Clearing objects and saving files.*/
  inList.close();
  cWaves.Clear();
  cWaves.Print((TString)(plotSaveFolder+"/waves.pdf)"),"pdf");
  cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf)"),"pdf");
  cTrig.Print((TString)(plotSaveFolder+"/trig.pdf)"),"pdf");
  cSignal.Print((TString)(plotSaveFolder+"/signal.pdf)"),"pdf");
  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}