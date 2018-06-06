#include "Plot2.h"

using namespace std;

//function: "main"                                                                                                          
void Plot3(){

  //directory where ROOT file is stored
  TString dir = "condor/finishedJobs/";

  //ROOT file path for accesing variables to be plotted
  TString histPath1 = "DQMData/Run 1/Phase2OuterTracker/Run summary/Stubs/";
  TString histPath2 = "DQMData/Run 1/Phase2OuterTracker/Run summary/Stubs/NStubs/";
  TString treePath = "Phase2PixelStubs/eventTree";

  vector<TString> filename1;
  //filename1.push_back("DQM_TTbar_PU200_932_OT613_200_IT4025_deadFPix2Pos_noTrackingChange.root");
  filename1.push_back("step4_DQM_opt8s4l.root");
  filename1.push_back("step4_DQM_opt6s3l.root");
  filename1.push_back("step4_DQM_opt7s4l.root");
  filename1.push_back("step4_DQM_opt8s3l.root");
  filename1.push_back("DQM_TTbar_PU200_932_OT613_200_IT4025_deadFPix2Pos_noTrackingChange.root");

  vector<TString> filename2;
  //filename2.push_back("DQM_opt8s4l_noTrackingChange.root");
  filename2.push_back("Phase2PixelStubs_opt8s4l_stubFix.root");
  filename2.push_back("Phase2PixelStubs_opt6s3l_stubFix.root");
  filename2.push_back("Phase2PixelStubs_opt7s4l_stubFix.root");
  filename2.push_back("Phase2PixelStubs_opt8s3l_stubFix.root");
  filename2.push_back("DQM_opt8s4l_noTrackingChange.root");

  vector<TString> geomID;
  //geomID.push_back("8s4l_deadDisc");
  geomID.push_back("8s4l");
  geomID.push_back("6s3l");
  geomID.push_back("7s4l");
  geomID.push_back("8s3l");
  geomID.push_back("8s4l_deadFPix2Pos_noTrackingChange");

  vector<Leaf> leaves;
  Leaf nstub,eta,barrel,ring,disc,discP1,discP2,discP3,discP4,discP5,discN1,discN2,discN3,discN4,discN5,disc_BW,disc_FW;
  nstub.SetValues("nstub","Number of Stubs","Events","Number of Stubs Per Event",50,6000,23000);leaves.push_back(nstub);
  eta.SetValues("Stub_Eta","Stub #eta","Number of Stubs", "Stub #eta");leaves.push_back(eta);
  barrel.SetValues("NStubs_Barrel","Barrel Layer","Number of Stubs","Stubs Per Barrel Layer");leaves.push_back(barrel);
  ring.SetValues("NStubs_Endcap_Ring","Endcap Ring","Number of Stubs","Stubs Per Endcap Ring");leaves.push_back(ring);
  disc.SetValues("NStubs_Endcap_Disc","Endcap Disc","Number of Stubs","Stubs Per Endcap Disc");leaves.push_back(disc);
  discP1.SetValues("NStubs_Disc+1","Endcap Ring (Disc +1)","Number of Stubs","Stubs Per Endcap Disc (+1)");leaves.push_back(discP1);
  discP2.SetValues("NStubs_Disc+2","Endcap Ring (Disc +2)","Number of Stubs","Stubs Per Endcap Disc (+2)");leaves.push_back(discP2);
  discP3.SetValues("NStubs_Disc+3","Endcap Ring (Disc +3)","Number of Stubs","Stubs Per Endcap Disc (+3)");leaves.push_back(discP3);
  discP4.SetValues("NStubs_Disc+4","Endcap Ring (Disc +4)","Number of Stubs","Stubs Per Endcap Disc (+4)");leaves.push_back(discP4);
  discP5.SetValues("NStubs_Disc+5","Endcap Ring (Disc +5)","Number of Stubs","Stubs Per Endcap Disc (+5)");leaves.push_back(discP5);
  discN1.SetValues("NStubs_Disc-1","Endcap Ring (Disc -1)","Number of Stubs","Stubs Per Endcap Disc (-1)");leaves.push_back(discN1);
  discN2.SetValues("NStubs_Disc-2","Endcap Ring (Disc -2)","Number of Stubs","Stubs Per Endcap Disc (-2)");leaves.push_back(discN2);
  discN3.SetValues("NStubs_Disc-3","Endcap Ring (Disc -3)","Number of Stubs","Stubs Per Endcap Disc (-3)");leaves.push_back(discN3);
  discN4.SetValues("NStubs_Disc-4","Endcap Ring (Disc -4)","Number of Stubs","Stubs Per Endcap Disc (-4)");leaves.push_back(discN4);
  discN5.SetValues("NStubs_Disc-5","Endcap Ring (Disc -5)","Number of Stubs","Stubs Per Endcap Disc (-5)");leaves.push_back(discN5);
  disc_FW.SetValues("NStubs_Endcap_Disc_Fw","Forward Endcap Disc","Number of Stubs","Stubs Per Endcap Disc (z>0)");leaves.push_back(disc_FW);
  disc_BW.SetValues("NStubs_Endcap_Disc_Bw","Backward Endcap Disc","Number of Stubs","Stubs Per Endcap Disc (z<0)");leaves.push_back(disc_BW);

  bool scale = true;
  bool log = false;

  int fileNum = filename1.size();
  int leafNum = leaves.size();

  //Tree and hist arrays
  TTree *rootTreeArr[fileNum];
  TH1D *histArr[fileNum][leafNum]; 

  //Don't display canvas when drawing
  gROOT->SetBatch(kTRUE);

  //Get or create histograms to be plotted and make single plots for each
  for (int i = 0; i < fileNum; ++i) {
    //Open ROOT Files
    TFile *f1 = TFile::Open(dir+filename1[i]);
    TFile *f2 = TFile::Open(dir+filename2[i]);

    for (int j = 0; j < leafNum; ++j) {

      if (leaves[j].GetName() == "nstub") histArr[i][j] = CreateHistogram1(rootTreeArr[i],histArr[i][j],dir+filename2[i],treePath,geomID[i],leaves[j]);
      else if (leaves[j].GetName() == "Stub_Eta") histArr[i][j] = (TH1D*)f1->Get(histPath1+leaves[j].GetName());
      else histArr[i][j] = (TH1D*)f1->Get(histPath2+leaves[j].GetName());

      histArr[i][j]->Sumw2();

      SinglePlot(histArr[i][j],geomID[i],leaves[j],scale,log);
    }
  }

  //Make comparison and ratio plots for all variables (Compares all other histograms with the one in the first file)
  for (int i = 0; i < leafNum; ++i) {
    for (int j = 1; j < fileNum; ++j) {
      ComparisonPlot(histArr[0][i],histArr[j][i],leaves[i],geomID[0],geomID[j],scale,log);
      RatioPlot(histArr[0][i],histArr[j][i],leaves[i],geomID[0],geomID[j],scale,log);
    }
  }
}

