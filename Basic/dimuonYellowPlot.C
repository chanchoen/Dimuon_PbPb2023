//
// dimuonYellowPlot.C
// Plot the dimuon invariant mass from CMS data
// input: oniaTree
// you have to provide the location+filename of the input oniaTree
//


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <Riostream.h>

#include <TInterpreter.h>
#include <TROOT.h>
#include <TSystem.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1F.h>

#include <TLatex.h>
#include <TLegend.h>

#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include "CMS/tdrstyle.C"
#include "CMS/CMS_lumi.C"

#endif

void dimuonYellowPlot(bool isAA=true,
		      const char* inFileAA = "/eos/cms/store/group/phys_heavyions/dileptons/Data2023/Oniatrees/inputfile.root",
		      bool addLogo = false,
		      float ptMu=4)
{
  cout << "dimuonYellowPlot: Starting macro dimuonYellowPlot" << endl;
  cout << "dimuonYellowPlot: Setting styles..." << endl;

  setTDRStyle();
 
  // binning of the histo
  double bins[100000];
  bins[0] = 0.2;
  int nBins = 0;
  for (int i=1; bins[i-1]<200; i++) {
    bins[i] = bins[i-1]*1.015;//1.04;//1.015;//
    nBins++;
  }

  cout << "dimuonYellowPlot: making mass histogram with " << nBins << " bins" << endl;
  TH1F *phAA = new TH1F("phAA","phMain;m_{#mu#mu} (GeV/c^{2});Events/(GeV/c^{2})",nBins,bins);
  // reading the input oniaTree
  cout << "dimuonYellowPlot: Opening input file ..."<< endl;
  TFile *infAA, *infPA, *infPP, *infXX;
  if (isAA){
    infAA = new TFile(Form("%s",inFileAA),"READ");
    // each muons is (glb&Trk), trigger = HLT_HIL1DoubleMu0_HighQ
    TTree *tree=(TTree*)infAA->Get("hionia/myTree");
    short reco_QQ_mumi_idx[1000];
    short reco_QQ_mupl_idx[1000];
    short reco_QQ_type[1000];
    short reco_QQ_sign[1000];
    long reco_QQ_trig[1000];
    short reco_QQ_size;
    short reco_mu_charge[1000];
    vector<float> *reco_QQ_4mom_m = new vector<float>;
    vector<float> *reco_QQ_4mom_pt = new vector<float>;
    vector<float> *reco_QQ_4mom_eta = new vector<float>;
    
    vector<float> *reco_mu_4mom_m = new vector<float>;
    vector<float> *reco_mu_4mom_pt = new vector<float>;
    vector<float> *reco_mu_4mom_eta = new vector<float>;

    tree->SetBranchAddress("Reco_QQ_mumi_idx", reco_QQ_mumi_idx);
    tree->SetBranchAddress("Reco_QQ_mupl_idx", reco_QQ_mupl_idx);
    tree->SetBranchAddress("Reco_QQ_size", &reco_QQ_size);
    tree->SetBranchAddress("Reco_QQ_trig", &reco_QQ_trig);
    tree->SetBranchAddress("Reco_QQ_type", &reco_QQ_type);
    tree->SetBranchAddress("Reco_QQ_sign", &reco_QQ_sign);
    tree->SetBranchAddress("Reco_QQ_4mom_m", &reco_QQ_4mom_m);
    tree->SetBranchAddress("Reco_QQ_4mom_pt", &reco_QQ_4mom_pt);
    tree->SetBranchAddress("Reco_QQ_4mom_eta", &reco_QQ_4mom_eta);
    tree->SetBranchAddress("Reco_mu_4mom_m", &reco_mu_4mom_m);
    tree->SetBranchAddress("Reco_mu_4mom_pt", &reco_mu_4mom_pt);
    tree->SetBranchAddress("Reco_mu_4mom_eta", &reco_mu_4mom_eta);
    tree->SetBranchAddress("Reco_mu_charge", reco_mu_charge);
    
    int nEvent = tree->GetEntries();
    for(int evt = 0; evt < nEvent; evt++){
      tree->GetEntry(evt);        
      for(int irqq = 0; irqq < reco_QQ_size; irqq++){
        float qq_M = reco_QQ_4mom_m->at(irqq);
        float qq_pt = reco_QQ_4mom_pt->at(irqq);
        float qq_eta = reco_QQ_4mom_eta->at(irqq);
        float mupl_pt, mumi_pt;
        if (reco_mu_charge[irqq] == 1) {
          mupl_pt = reco_mu_4mom_pt->at(reco_QQ_mupl_idx[irqq]);
        }
        if (reco_mu_charge[irqq] == -1){
          mumi_pt = reco_mu_4mom_pt->at(reco_QQ_mumi_idx[irqq]);
        }
        // cout<< qq_M<<endl;
        if((reco_QQ_trig[irqq]&9) == 9  && reco_QQ_sign[irqq] == 0 && qq_pt > ptMu){
          phAA->Fill(qq_M);
        }
      }
    }
  }

  //------------------------------------- DRAWING SECTION !!!!!
  TH1 *phAxis = new TH1D("phAxis","phAxis;m_{#mu#mu} (GeV/c^{2});Events/(GeV/c^{2})",1,2,200);
  phAxis->GetYaxis()->SetRangeUser(0.1,2e7);
  phAxis->SetDirectory(0);

  TCanvas* yellowPlot = new TCanvas("yellowPlot","yellowPlot");
  phAxis->Draw();  
  gPad->SetLogx();
  gPad->SetLogy();
 
  phAA->SetLineWidth(2);
  phAA->SetLineColor(kBlue+1);
  phAA->SetFillColor(kYellow);
  phAA->Draw("HIST same");
 
  gPad->RedrawAxis();

  //
  // Options to be used with the CMS_lumi macro (for what options are available see CMS_lumi.h and CMS_lumi.C)
  //
  //
  // iPeriod options: 99 for pPb 5.02 TeV, 101 for PbPb 2011, 102 for pp 2013, 104 for pp 2015, 105 for PbPb 2015
  //
    int iPeriod      = 104;   
    if(isAA) iPeriod = 101;
    lumiTextOffset   = 0.3; // default 0.28

  // Call the CMS_lumi macro to draw:
  // CMS preliminary, aligned on the right and justified (iPos=33, third argument)
  // iPeriod: integrated luminosity (drawn on top left, out of frame, or use lumiTextOffset)
  // center of mass energy (drawn on top right, out of frame, or use lumiTextOffset)
  
  CMS_lumi(yellowPlot,iPeriod,33);

  float H = gPad->GetWh();
  float W = gPad->GetWw();
  float l = gPad->GetLeftMargin();
  float t = gPad->GetTopMargin();
  float r = gPad->GetRightMargin();
  float b = gPad->GetBottomMargin();
  float posX_ =   l + 0.045*(1-l-r)*W/H;
  float posY_ = 1-t - 0.045*(1-t-b);
  float xl_0 = posX_;
  float yl_0 = posY_ - 0.1;
  float xl_1 = posX_ + 0.1*H/W;
  float yl_1 = posY_;
  TASImage* CMS_logo = new TASImage("CMS-Color.gif");
  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
  if(addLogo)
    {
      pad_logo->Draw();
      pad_logo->cd();
      CMS_logo->Draw("X");
      pad_logo->Modified();
    }
  yellowPlot->cd();
  
  // Booking the TLatex class and setting the parameters
  // CMS Preliminary (not used anymore) left as a dummy
  TLatex* latex = new TLatex(0.62,0.88,"CMS Preliminary");
  latex->SetNDC();
  /*
  latex->DrawLatex(0.24,0.88,"#rho, #omega");
  latex->DrawLatex(0.32,0.88,"#phi");
  */
  latex->DrawLatex(0.22,0.8,"J/#psi");
  // tex = new TLatex(0.38,0.94,"J/#psi");
  latex->DrawLatex(0.3,0.7,"#psi(2S)");
  //  tex = new TLatex(0.45,0.82,"#psi(2S)");
  latex->DrawLatex(0.45,0.7,"#varUpsilon(1,2,3S)");
  //  tex = new TLatex(0.59,0.76,"#varUpsilon(1,2,3S)");
  latex->DrawLatex(0.8,0.5,"Z");
  // tex = new TLatex(0.84,0.57,"Z");

  latex->DrawLatex(0.7,0.7,Form("p_{T}^{#mu#mu} > %.0f GeV/c",ptMu));

  yellowPlot->SaveAs(Form("%s_ptMu%.0f_isAA%d_hasLogo%d.pdf",yellowPlot->GetTitle(),ptMu,isAA,addLogo));
  yellowPlot->SaveAs(Form("%s_ptMu%.0f_isAA%d_hasLogo%d.png",yellowPlot->GetTitle(),ptMu,isAA,addLogo));
  //-----------------------------------
  // TCanvas* pc1 = new TCanvas("pc1","pc1");
  // phMain->Divide(phComp);
  // phMain->Draw();
  return;
}
