#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"

#include "TDRStyle.h"


int plotCutOptim(){


  bool doPU = true;
  //std::string histo = "mht30";
  std::string histo = "metnolep";
  //std::string histo = "Mjjdetajj";
  //std::string histo = "Mjjdphijj";
  //std::string histo = "dphijj";
  //std::string histo = "detajj";
  //std::string histo = "Mjj";
  //std::string histo = "Jet1_pt";
  //std::string histo = "Jet2_pt";
  bool do2D = histo.find("Mjjd")!=histo.npos;

  //std::string plotSuf = "_MjjdetajjCheck";
  //std::string plotSuf = "_mhtTuning";
  std::string plotSuf = "_metmhtCheck";

  unsigned ngrid = histo.find("metnolep")!=histo.npos ||histo.find("mht30")!=histo.npos  ? 30:histo.find("Mjj")!=histo.npos ? 40 : histo.find("dphi")!=histo.npos ? 36 : 30;
  TFile *fin = TFile::Open("Mjj1300deta4met0dphi18mht250/plots0e0mu.root");
  if (!fin) return 1;
  else {
    std::cout << " File found: pointer " << fin << std::endl << std::flush;
    std::cout << fin->GetName() << std::endl;
  }
  fin->cd();

  std::string lSigName = "Signal_"+histo+"_";
  if (doPU) lSigName += "200PU";
  else lSigName += "noPU";
  std::string lBkgName = "Background_"+histo+"_";
  if (doPU) lBkgName += "200PU";
  else lBkgName += "noPU";

  TH1F *num = 0;
  TH1F *den = 0;
  TH2F *num2D = 0;
  TH2F *den2D = 0;

  bool intUp = histo.find("dphijj")==histo.npos;

  if (do2D){
    std::cout << " -- 2D histo!" << std::endl;
    num2D = (TH2F*)(gDirectory->Get(lSigName.c_str()))->Clone("signal");
    den2D = (TH2F*)(gDirectory->Get(lBkgName.c_str()))->Clone("background");
  }
  else {
    std::cout << " -- 1D histo!" << std::endl;
    num = (TH1F*)(gDirectory->Get(lSigName.c_str()))->Clone("signal");
    den = (TH1F*)(gDirectory->Get(lBkgName.c_str()))->Clone("background");
  }
  
  if (do2D && (!num2D || !den2D)) return 1;
  else if (!do2D && (!num || !den)) return 1;

  std::cout << " Signal histo: " << lSigName << " " << (do2D?num2D->Integral():num->Integral()) << std::endl
	    << " Background histo: " << lBkgName << " " << (do2D?den2D->Integral():den->Integral()) << std::endl;

  SetTdrStyle();

  TCanvas *myc = new TCanvas("myc","myc");
  TCanvas *mycR = new TCanvas("mycR","mycR");

  gStyle->SetOptStat(0);


  if (!do2D){
    myc->cd();
    //den->Rebin(5);
    den->SetLineColor(1);
    den->SetMarkerColor(1);
    den->SetMarkerStyle(22);
    den->GetXaxis()->SetNdivisions(505);
    den->Draw("hist");
    //num->Rebin(5);
    num->SetLineColor(2);
    num->SetMarkerColor(2);
    num->SetMarkerStyle(20);
    num->Draw("PEsame");
    den->Draw("same");
  }
  else {
    num2D->GetYaxis()->SetRangeUser(500,5000);
    myc->Divide(2,1);
    myc->cd(1);
    gPad->SetRightMargin(0.15);
    den2D->Draw("colz");
    myc->cd(2);
    num2D->Draw("colz");
  }
  myc->Update();
  myc->Print(("OPTIM/"+histo+plotSuf+".pdf").c_str());
  
  mycR->cd();
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);
  if (do2D) mycR->SetRightMargin(0.15);

  TH1F *ratio = 0;
  TH1F *ratioS = 0;
  TH2F *ratio2D = 0;
  TH2F *ratioS2D = 0;

  if (do2D){
    ratio2D = (TH2F*)num2D->Clone("ratio2D");
    ratioS2D = (TH2F*)num2D->Clone("ratioS2D");
  }
  else {
    ratio = (TH1F*)num->Clone("ratio");
    ratioS = (TH1F*)num->Clone("ratioS");
    if (!ratio || !ratioS) return 1;
  }
  double maxR = 0;
  double maxS = 0;
  int maxRibin = 0;
  int maxRjbin = 0;
  int maxi = do2D?ratio2D->GetNbinsX()+1:ratio->GetNbinsX()+1;
  int maxj = do2D?ratio2D->GetNbinsY()+1:1;
  for (int ibin(0); ibin<maxi; ibin++){
    for (int jbin(0); jbin<maxj; jbin++){
      
      double errS=0,errB=0;
      double signal = 0;
      double bkg = 0;
      if (!do2D){
	signal = intUp?num->IntegralAndError(ibin,ratio->GetNbinsX()+1,errS):num->IntegralAndError(0,ibin,errS);
	bkg = intUp?den->IntegralAndError(ibin,ratio->GetNbinsX()+1,errB):den->IntegralAndError(0,ibin,errB);
      }
      else {
	signal = intUp?num2D->Integral(ibin,ratio2D->GetNbinsX()+1,jbin,ratio2D->GetNbinsY()+1):num2D->Integral(0,ibin,jbin,ratio2D->GetNbinsY()+1);
	bkg = intUp?den2D->Integral(ibin,ratio2D->GetNbinsX()+1,jbin,ratio2D->GetNbinsY()+1):den2D->Integral(0,ibin,jbin,ratio2D->GetNbinsY()+1);
      }
      
      if (bkg<=0) continue;
      double asimov = sqrt(2*((signal+bkg)*log(1+signal/bkg)-signal));
      double errSig = sqrt(pow(1/sqrt(bkg)*errS,2)+pow(0.5*signal*errB*pow(bkg,-3./2),2));
      if (do2D) std::cout << ibin << " " << jbin << " " << ratio2D->GetXaxis()->GetBinLowEdge(ibin) << " " << ratio2D->GetYaxis()->GetBinLowEdge(jbin) << " " << signal << " " << bkg << " " << asimov << std::endl;
      else {
	std::cout << ibin << " " << ratio->GetXaxis()->GetBinLowEdge(ibin) << " " << signal << " " << bkg << " " << asimov << std::endl;
      }

      if (!do2D){
	ratio->SetBinContent(ibin,asimov);
	ratioS->SetBinContent(ibin,signal/sqrt(bkg));
	ratio->SetBinError(ibin,errSig);
	ratioS->SetBinError(ibin,errSig);
      } else {
	ratio2D->SetBinContent(ibin,jbin,asimov);
	ratioS2D->SetBinContent(ibin,jbin,signal/sqrt(bkg));
      }
      if (asimov>maxR) {
	maxR = asimov;
	maxRibin = ibin;
	maxRjbin = jbin;
      }
      if (signal/sqrt(bkg)>maxS) {
	maxS = signal/sqrt(bkg);
      }
    }
  }

  if (do2D) {
    std::cout << " Maximum Asimov found at bin: " << maxRibin << " " << maxRjbin << " " << ratio2D->GetXaxis()->GetBinLowEdge(maxRibin) << " " << ratio2D->GetYaxis()->GetBinLowEdge(maxRjbin) << " " << maxR << std::endl;
  }

  if (!do2D){
    ratio->SetMaximum(maxR*1.1);
    ratio->GetYaxis()->SetTitle("Asimov sig.");
    ratio->SetTitle("");
    ratio->SetLineColor(3);
    ratio->SetFillColor(5);
    ratio->GetXaxis()->SetNdivisions(510);
    ratio->Draw("hist");
    
    ratioS->SetLineColor(2);
    ratioS->SetMarkerColor(2);
    ratioS->SetMarkerStyle(20);
    ratioS->Draw("PEsame");
  
    TPad *grid = new TPad("grid","",0,0,1,1); 
    grid->Draw();
    grid->cd();
    grid->SetGrid();
    grid->SetFillStyle(4000); 
    
    TH2 *hgrid = new TH2C("hgrid","",ngrid,ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),1,0,ratio->GetMaximum());
    hgrid->Draw();
    hgrid->GetXaxis()->SetNdivisions(ngrid);
    hgrid->GetYaxis()->SetNdivisions(1);
    hgrid->GetYaxis()->SetLabelOffset(999.); 
    hgrid->GetXaxis()->SetLabelOffset(999.); 
  } else {
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    ratio2D->GetZaxis()->SetTitle("Asimov sig.");
    ratio2D->SetTitle("");
    ratio2D->Draw("colz"); 
  }

  mycR->Update();
  mycR->Print(("OPTIM/AsimovSig"+histo+plotSuf+".pdf").c_str());


  return 0;
}//main
