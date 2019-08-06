#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"
#include "THStack.h"
#include "TColor.h"


int checkggRatio(){

  bool doMHT = false;

  const unsigned nbins = doMHT ? 48 : 28;

  TFile *fvbf = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/amagnan/output_lighttree_180718_preFiring/MC_Powheg-VBFHtoinv-mH125.root");
  TFile *fgg = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/amagnan/output_lighttree_180718_preFiring/MC_Powheg-GluGluHtoinv-mH125.root");

  TH2F *hvbf = new TH2F("hvbf",";E_{T}^{miss} (GeV);M_{jj} (GeV);events",nbins,130,130+10*nbins,30,1000,4000);
  if (doMHT) hvbf->GetXaxis()->SetTitle("H_{T}^{miss} (GeV)");
  TH2F *hgg = new TH2F("hgg",";E_{T}^{miss} (GeV);M_{jj} (GeV);events",nbins,130,130+10*nbins,30,1000,4000);
  if (doMHT) hgg->GetXaxis()->SetTitle("H_{T}^{miss} (GeV)");
  TH2F *hr1 = new TH2F("hr1",";E_{T}^{miss} (GeV);M_{jj} (GeV);ggH/VBFH",nbins,130,130+10*nbins,30,1000,4000);
  if (doMHT) hr1->GetXaxis()->SetTitle("H_{T}^{miss} (GeV)");

  TCanvas *c1 = new TCanvas("c1", "myPlots",0,67,600,600);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  gStyle->SetPalette(1);

  //std::string lcut1 = "total_weight_lepveto*(nvetomuons==0&&nvetoelectrons==0&&nvetotaus==0&&n_jets_csv2medium==0&&jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1000 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>130 && nloosephotons==0 && dijet_dphi<1.5 && fourjetsmetnomu_mindphi>0.5)";
  std::string lcut1;
  if (!doMHT) lcut1 = "weight_xsection*(nvetomuons==0&&nvetoelectrons==0&&nvetotaus==0&&n_jets_csv2medium==0&&jet1_eta*jet2_eta<0 && abs(jet1_eta)<5 && abs(jet2_eta)<5 && dijet_M>1000 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>130 && nloosephotons==0 && dijet_dphi<1.5 && fourjetsmetnomu_mindphi>0.5)";
  else lcut1 = "weight_xsection*(nvetomuons==0&&nvetoelectrons==0&&nvetotaus==0&&n_jets_csv2medium==0&&jet1_eta*jet2_eta<0 && abs(jet1_eta)<5 && abs(jet2_eta)<5 && dijet_M>1000 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>0 && mht>130 && nloosephotons==0 && dijet_dphi<1.5 && fourjetsmetnomu_mindphi>0.5)";

  fvbf->cd();
  TTree *tvbf = (TTree*)gDirectory->Get("LightTree");
  fgg->cd();
  TTree *tgg = (TTree*)gDirectory->Get("LightTree");

  c1->cd();
  if (!doMHT){
    tvbf->Draw("dijet_M:metnomuons>>hvbf",lcut1.c_str());
    tgg->Draw("dijet_M:metnomuons>>hgg",lcut1.c_str());
  }
  else {
    tvbf->Draw("dijet_M:mht>>hvbf",lcut1.c_str());
    tgg->Draw("dijet_M:mht>>hgg",lcut1.c_str());
  }

  for (int i(1);i<nbins;++i){ 
    for (int j(1);j<31;++j){ 
      hr1->SetBinContent(i,j,hgg->Integral(i,nbins,j,31)/hvbf->Integral(i,nbins,j,31));
    }
  }

  std::cout << hgg->Integral(13,nbins,4,31) << " " << hvbf->Integral(13,nbins,4,31) << " " << hr1->GetBinContent(1,1) << std::endl;

  hr1->GetZaxis()->SetTitleOffset(1.5);
  hr1->Draw("colz");

  c1->Update();
  c1->Print(doMHT?"ggHoverVBFvsMjjvsMHT.pdf":"ggHoverVBFvsMjjvsMET.pdf");

  TFile *outfile;
  if (doMHT) outfile = new TFile("ggHoverVBFvsMjjvsMHT.root","RECREATE");
  else outfile = new TFile("ggHoverVBFvsMjjvsMET.root","RECREATE");
  outfile->cd();
  hr1->Write();
  outfile->Write();

  return 0;
}//main
