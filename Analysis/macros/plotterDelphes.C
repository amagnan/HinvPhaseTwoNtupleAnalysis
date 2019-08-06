#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>

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

#include "tdrstyle.C"
#include "CMS_lumi.C"

float getggRatio(const float mjj, const float met, const float mht){

  TFile *file = TFile::Open(met<1?"ggHoverVBFvsMjjvsMHT.root":"ggHoverVBFvsMjjvsMET.root");
  file->cd();
  TH2F *htmp = (TH2F*)gDirectory->Get("hr1");


  int binx = htmp->GetXaxis()->FindBin(met<1?mht+0.1 : met+0.1);
  int biny = htmp->GetYaxis()->FindBin(mjj+0.1);

  if (binx==0) binx=1;
  if (binx==htmp->GetXaxis()->GetNbins()+1) binx=htmp->GetXaxis()->GetNbins();
  if (biny==0) biny=1;
  if (biny==htmp->GetYaxis()->GetNbins()+1) biny=htmp->GetYaxis()->GetNbins();


  float result = htmp->GetBinContent(binx,biny);

  std::cout << " Check: " << met << " " << mht << " " << mjj << " " << result << std::endl;

  return result;
}

void SetHistStyle(TH1 *histPlot,const std::string & histX){

  /*histPlot->SetFillStyle(1001);
  histPlot->SetLineStyle(0);
  histPlot->SetMarkerStyle(20);
  histPlot->GetXaxis()->SetLabelFont(42);
  histPlot->GetXaxis()->SetLabelOffset(0.007);
  histPlot->GetXaxis()->SetLabelSize(0.05);
  histPlot->GetXaxis()->SetTitleSize(0.06);
  histPlot->GetXaxis()->SetTitleOffset(0.9);
  histPlot->GetXaxis()->SetTitleFont(42);
  histPlot->GetYaxis()->SetLabelFont(42);
  histPlot->GetYaxis()->SetLabelOffset(0.007);
  histPlot->GetYaxis()->SetLabelSize(0.05);
  histPlot->GetYaxis()->SetTitleSize(0.06);
  histPlot->GetYaxis()->SetTitleOffset(1.25);
  histPlot->GetYaxis()->SetTitleFont(42);
  histPlot->GetZaxis()->SetLabelFont(42);
  histPlot->GetZaxis()->SetLabelOffset(0.007);
  histPlot->GetZaxis()->SetLabelSize(0.05);
  histPlot->GetZaxis()->SetTitleSize(0.06);
  histPlot->GetZaxis()->SetTitleFont(42);*/
  histPlot->SetTitle(histX.c_str());
}

  void getNgen(const std::string & paramfile,
	     std::map<std::string,unsigned> & aMap){

  std::ifstream lInput;
  lInput.open(paramfile.c_str());
  if (!lInput.is_open()){
    std::cout << "Could not open input param file " << paramfile << std::endl;
    exit(1);
  }
  while(!lInput.eof()){
    std::string tmp;
    unsigned Ngen = 0;
    lInput>>tmp>>Ngen;
    
    if (Ngen>0) aMap.insert(std::pair<std::string,unsigned>(tmp,Ngen));
  }
  std::cout << " -- Found " << aMap.size() << " samples in file " << paramfile << std::endl;

};


double getXS(const std::string & process){
  //in pb...
  if (process.find("VBFH")!=process.npos) return 4.278;
  else if (process.find("QCD_Mdijet-1000toInf")!=process.npos) return 1135;
  else if (process.find("EWKZ2Jets_ZToLL")!=process.npos) return 4.38;
  else if (process.find("EWKZ2Jets_ZToNuNu")!=process.npos) return 11.05;
  else if (process.find("EWKWMinus2Jets")!=process.npos) return 23.6;
  else if (process.find("EWKWPlus2Jets")!=process.npos) return 30.37;
  else if (process.find("ST_s-channel")!=process.npos) return 11.14;
  else if (process.find("ST_tW_top")!=process.npos) return 45.06;
  else if (process.find("ST_tW_antitop")!=process.npos) return 45.02;
  else if (process.find("ST_tch_top")!=process.npos) return 48.03;
  else if (process.find("ST_tch_antitop")!=process.npos) return 29.2;
  else if (process.find("TT")!=process.npos) return 864.4;
  else if (process.find("W0JetsToLNu")!=process.npos) return 47340.;
  else if (process.find("W1JetsToLNu")!=process.npos) return 10370.;
  else if (process.find("W2JetsToLNu")!=process.npos) return 2965.;
  else if (process.find("W3JetsToLNu")!=process.npos) return 1891;
  else if (process.find("ZJetsToNuNu_HT-100To200")!=process.npos) return 304.5;
  else if (process.find("ZJetsToNuNu_HT-200To400")!=process.npos) return 85.88;
  else if (process.find("ZJetsToNuNu_HT-400To600")!=process.npos) return 12.35;
  else if (process.find("ZJetsToNuNu_HT-600To800")!=process.npos) return 3.029;
  else if (process.find("ZJetsToNuNu_HT-800To1200")!=process.npos) return 1.400;
  else if (process.find("ZJetsToNuNu_HT-1200To2500")!=process.npos) return 0.345;
  else if (process.find("ZJetsToNuNu_HT-2500ToInf")!=process.npos) return 0.0069;
  else if (process.find("DYJetsToLL_M-50_HT-70to100")!=process.npos) return 161.6;
  else if (process.find("DYJetsToLL_M-50_HT-100to200")!=process.npos) return 150.;
  else if (process.find("DYJetsToLL_M-50_HT-200to400")!=process.npos) return 32.95;
  else if (process.find("DYJetsToLL_M-50_HT-400to600")!=process.npos) return 3.911;
  else if (process.find("DYJetsToLL_M-50_HT-600to800")!=process.npos) return 0.8301;
  else if (process.find("DYJetsToLL_M-50_HT-800to1200")!=process.npos) return 0.3852;
  else if (process.find("DYJetsToLL_M-50_HT-1200to2500")!=process.npos) return 0.08874;
  else if (process.find("DYJetsToLL_M-50_HT-2500toInf")!=process.npos) return 0.001755;
  else {
    std::cout << " Cross section not found for process " << process << " Please add !" << std::endl;
    exit(1);
  }
  return 0;
  
}



int getRegion(const int I, const std::string & prodDate, const std::string & plotbasedir, const double & mhtcut, const double & metcut, const double & mjjcut, const double & detajjcut, const double & dphijjcut, const double & ggFfrac, const bool doAllPlots=false, const unsigned doJes = 0)
{
  double lumi = 3000;//in fb-1...

  setTDRStyle();

  //Set global variables for CMS style
  writeExtraText = true;       // if extra text
  cmsText     = "CMS Phase-2";
  extraText  = "Simulation Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "3000 fb^{-1} (14 TeV)";
  relPosX = 0.28;

  std::map<std::string,unsigned> nEvtMap;


  std::string paramFile = "/afs/cern.ch/work/a/amagnan/UPSGAna/HinvPhaseTwoNtupleAnalysis/Analysis/scripts/ParamFile_"+prodDate+".dat";
  getNgen(paramFile,nEvtMap);
  std::map<std::string,unsigned>::iterator  nEvtMapIter = nEvtMap.begin();
  for (; nEvtMapIter!=nEvtMap.end(); ++nEvtMapIter){
    std::cout << nEvtMapIter->first << " " << nEvtMapIter->second << std::endl;
  }

  const unsigned nProc = 9+12;
  const unsigned nmaxvar = 11;
  const unsigned nvar = doAllPlots?(I>1?nmaxvar:nmaxvar-1):1;
  const unsigned nsel = 6;
  const unsigned nSamples = 29;
  const unsigned npu = 1;

  //std::string histName[nvar] = {"metnolep","Mjj"};
  //std::string histNameShort[nvar] = {"metnolep","Mjj"};
  //std::string binRange[nvar] = {"(41,190,600)","(45,1500,6000)"};
  //std::string histLabelX[nvar] = {";E_{T}^{miss} (GeV)",";M_{jj} (GeV)"};

  std::string histName[nmaxvar] = {(I==3?"Mee":I==5?"Mmumu":I==2?"ele_mt":I==4?"mu_mt":"metnolep"),"detajj" ,"mht30","dphijj","Mjj","Jet1_pt","Jet2_pt","Jet1_eta","Jet2_eta","jetmetnolepmindphi","metnolep"};
  std::string histNameShort[nmaxvar] = {(I==3?"Mee":I==5?"Mmumu":I==2?"ele_mt":I==4?"mu_mt":"metnolep"),"detajj" ,"mht30","dphijj","Mjj","Jet1pt","Jet2pt","Jet1eta","Jet2eta","jetmetnolepmindphi","metnolep"};
  std::string binRange[nmaxvar] = {((I==3 || I==5)?"(30,60,120)":(I==2 || I==4)?"(20,0,160)":"(20,190,600)"),"(21,4.,9.25)","(20,100,2500)","(18,0,3.24)","(20,2500,5000)","(25,80,500)","(25,40,300)","(20,-5,5)","(20,-5,5)","(15,0.5,3.1416)","(20,190,600)"};
  std::string histLabelX[nmaxvar] = {(I==3?";M_{ee} (GeV);Events":I==5?";M_{#mu#mu} (GeV);Events":I==2?";m_{T}^{ele} (GeV);Events":I==4?";m_{T}^{#mu} (GeV);Events":";E_{T}^{miss} (GeV);Events"),";#Delta#eta_{jj};Events" ,";H_{T}^{miss} (GeV);Events",";#Delta#phi_{jj};Events",";M_{jj} (GeV);Events",";p_{T}^{j1} (GeV);Events",";p_{T}^{j2} (GeV);Events",";#eta^{j1};Events",";#eta^{j2};Events",";Min#Delta#phi(jet,E_{T}^{miss});Events",";E_{T}^{miss} (GeV);Events"};

  //std::string histName[nvar] = {(I==3?"Mee":I==5?"Mmumu":I==2?"ele_mt":I==4?"mu_mt":"metnolep"),((I==3 || I==2)?"ele1_pt":(I==5 || I==4)?"mu1_pt":"GenLep1_pt"),((I==3 || I==2)?"ele1_eta":(I==5 || I==4)?"mu1_eta":"GenLep1_eta")};
  //std::string binRange[nvar] = {((I==3 || I==5)?"(60,60,120)":(I==2 || I==4)?"(60,0,300)":"(100,0,600)"),"(100,0,100)","(80,-4,4)"};
  //std::string histLabelX[nvar] = {(I==3?";M_{ee} (GeV)":I==5?";M_{#mu#mu} (GeV)":I==2?";m_{T}^{ele} (GeV)":I==4?";m_{T}^{#mu} (GeV)":";E_{T}^{miss} (GeV)"),((I==3 || I==2)?";p_{T}^{ele} (GeV)":(I==5 || I==4)?";p_{T}^{#mu} (GeV)":";p_{T}^{genLep} (GeV)"),((I==3 || I==2)?";#eta^{ele}":(I==5 || I==4)?";#eta^{#mu}":";#eta^{genLep}")};


  std::string ProcessNames[nSamples] = {
    "VBFH",
    "QCD_Mdijet-1000toInf",
    "ST_s-channel","ST_tch_top","ST_tW_top","ST_tW_antitop","TT",
    "EWKZ2Jets_ZToLL",
    "DYJetsToLL_M-50_HT-70to100",
    "DYJetsToLL_M-50_HT-100to200",
    "DYJetsToLL_M-50_HT-200to400",
    "DYJetsToLL_M-50_HT-400to600",
    "DYJetsToLL_M-50_HT-600to800",
    "DYJetsToLL_M-50_HT-800to1200",
    "DYJetsToLL_M-50_HT-1200to2500",
    "DYJetsToLL_M-50_HT-2500toInf",
    "EWKZ2Jets_ZToNuNu",
    "ZJetsToNuNu_HT-100To200",
    "ZJetsToNuNu_HT-200To400",
    "ZJetsToNuNu_HT-400To600",
    "ZJetsToNuNu_HT-600To800",
    "ZJetsToNuNu_HT-800To1200",
    "ZJetsToNuNu_HT-1200To2500",
    "EWKWMinus2Jets","EWKWPlus2Jets",
    "W0JetsToLNu","W1JetsToLNu","W2JetsToLNu","W3JetsToLNu"
  };
    
  std::string PU[npu] = {"200PU"};
  //2016 selection
  //std::string baseSelection = "Mjj>1300 && metnolep>250 && detajj>4 && dphijj<1.5 && jetmetnolepmindphi>0.5 && nmediumbjets==0 && Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && ntaus==0";
  std::ostringstream baseSelection;
  baseSelection << "Mjj>" << mjjcut << " && metnolep>" << metcut << " && mht30>" << mhtcut << " && detajj>" << detajjcut << " && dphijj<" << dphijjcut << " && jetmetnolepmindphi>0.5 && Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && ntaus==0 && nmediumbjets==0";

  std::string selection[nsel] = 
    {"",
     //" && (nlooseEle==0 || (nlooseEle==1 && TMath::Abs(ele1_eta)>2.5))  && (nlooseMu==0 || (nlooseMu==1 && TMath::Abs(mu1_eta)>2.4)) && ntightGamma==0",
     " && nlooseEle==0 && nlooseMu==0 && ntightGamma==0",
     " && ntightEle==1 && nlooseMu==0 && ele1_pt>40 && met>60 && ele_mt<160",// && TMath::Abs(ele1_eta)<2.5",
     " && nlooseEle==2 && ntightEle>=1 && ele1_pt>40 && nlooseMu==0 && Mee>60 && Mee<120",// && TMath::Abs(ele1_eta)<2.5",
     " && nlooseEle==0 && ntightMu==1 && mu_mt<160",// && TMath::Abs(mu1_eta)<2.4",
     " && nlooseEle==0 && nlooseMu==2 && ntightMu>=1 && Mmumu>60 && Mmumu<120"// && TMath::Abs(mu1_eta)<2.4"
    };

  std::string Regions[nsel] = {"","0e0mu","1e0mu","2e0mu","0e1mu","0e2mu"};
  enum  {VBFH, EWK_Zll,EWKZee,EWKZmumu,EWKZtautau, QCD_Zll,QCDZee,QCDZmumu,QCDZtautau,EWK_Z_nunu, QCD_Z_nunu, EWK_W_lnu,EWKWenu,EWKWmunu,EWKWtaunu, QCD_W_lnu,QCDWenu,QCDWmunu,QCDWtaunu,Top, QCD};
  std::string Yields_names [nProc] = {"VBFH", "Zll(EWK)","Zee(EWK)","Zmumu(EWK)","Ztautau(EWK)", "Zll(QCD)","Zee(QCD)","Zmumu(QCD)","Ztautau(QCD)","Z->nunu(EWK)","Z->nunu(QCD)", "W->lnu(EWK)","Wenu(EWK)","Wmunu(EWK)","Wtaunu(EWK)","W->lnu(QCD)","Wenu(QCD)","Wmunu(QCD)","Wtaunu(QCD)","Top","QCD"};
  
  std::string extraSel = "";
  //extraSel = " && (TMath::Abs(GenLep1_pdgid)==13 || (TMath::Abs(GenLep1_pdgid)==15 && TMath::Abs(GenLep2_pdgid)==13)) && GenLep1_pt>10 && TMath::Abs(GenLep1_eta)<2.4";
  //extraSel = " && TMath::Abs(GenLep1_pdgid)==15 && TMath::Abs(GenLep2_pdgid)==0 && GenLep1_pt>20 && TMath::Abs(GenLep1_eta)<2.3";

  for (unsigned iS(0); iS<nsel; ++iS){
    selection[iS] = baseSelection.str()+selection[iS]+extraSel;
  }


  TH1 *histPlot=0,*dataPlot=0;
  TH2 *histPlot2D=0,*dataPlot2D=0;
    
  TTree* LightTree = NULL;
  
  TCanvas *c1 = new TCanvas("c1", "myPlots",0,67,600,600);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-102.5,-10.38415,847.5,69.4939);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  double leg_xl = 0.62, leg_xr = 0.948, leg_yb = 0.62, leg_yt = 0.948 ;
  TLegend* leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt);
          
  std::ofstream Yield((plotbasedir+"/Yield_R"+Regions[I]+".txt").c_str());
  std::ofstream Datacard((plotbasedir+"/Datacard_R"+Regions[I]+".txt").c_str());
  Yield<<"Process "<<Regions[I]<<std::endl;
  //Yield<<"Format: N_sel, N_gen, w, w*N_sel, w*N_sel/N_gen"<<std::endl;
  TFile *outfile = TFile::Open((plotbasedir+"/plots"+Regions[I]+".root").c_str(),"RECREATE"); 
  if (!outfile) return 1;
  for (unsigned pu = 0; pu<npu; pu++)
    {//loop on PU
      double Process_Yields[nProc];
      double Process_Yields_err[nProc];
      Yield<<PU[pu]<<std::endl;
      for (unsigned yields = 0; yields<nProc;yields++)
	{
	  Process_Yields[yields] = 0;
	  Process_Yields_err[yields] = 0;
	}
      for (unsigned variable = 0; variable <nvar; variable++)
	{//loop on variables
	  bool do2D = histName[variable].find(":")!=histName[variable].npos;

	  THStack *A = new THStack("test","test");
	  TH1 *background = 0;
	  TH2 *background2D = 0;

	  TFile *fin[nSamples];

	  for (unsigned i = 0;i < nSamples; i++)
	    {
	      if (nEvtMap.find(ProcessNames[i])==nEvtMap.end()) {
		std::cout << " Sample " << ProcessNames[i] << " not found, skipping." << std::endl;
		continue;
	      }
	      fin[i] = 0;
	      if (doJes>0 && i==0) {
		if (doJes==1) fin[i] = TFile::Open(("/afs/cern.ch/work/a/amagnan/public/UPSGAna/"+prodDate+"/JESUP/200PU/"+ProcessNames[i]+"/HistosFile_"+ProcessNames[i]+"_"+PU[pu]+".root").c_str());
		else if (doJes==2) fin[i] = TFile::Open(("/afs/cern.ch/work/a/amagnan/public/UPSGAna/"+prodDate+"/JESDOWN/200PU/"+ProcessNames[i]+"/HistosFile_"+ProcessNames[i]+"_"+PU[pu]+".root").c_str());
		else fin[i] = TFile::Open(("/afs/cern.ch/work/a/amagnan/public/UPSGAna/"+prodDate+"/200PU/"+ProcessNames[i]+"/HistosFile_"+ProcessNames[i]+"_"+PU[pu]+".root").c_str());
	      }
	      else {
		fin[i] = TFile::Open(("/afs/cern.ch/work/a/amagnan/public/UPSGAna/"+prodDate+"/"+ProcessNames[i]+"_"+PU[pu]+".root").c_str());
	      }
	      //try second option if already merged
	      //if (!fin[i]) fin[i] = TFile::Open(("/afs/cern.ch/work/a/amagnan/public/UPSGAna/"+prodDate+"/"+ProcessNames[i]+"_"+PU[pu]+"/HistosFile_"+ProcessNames[i]+"_"+PU[pu]+".root").c_str());
	      if (!fin[i]) {
		continue;
	      }
	      fin[i]->GetObject("LightTree",LightTree);
	      if (LightTree!=NULL){

		//std::cout << LightTree << " " << LightTree->GetEntries() << std::endl;

		//if (variable==0) std::cout << " Nevts = " << nEvtMap[ProcessNames[i]] << std::endl;
		
		double weight = getXS(ProcessNames[i])*lumi*1000/nEvtMap[ProcessNames[i]];
		if (weight > 1000) {
		  weight = 0;
		  if (variable==0)std::cout << " -- Weight too large, ignoring this sample... setting weight to 0" << std::endl;
		}
		
		if (variable==0) std::cout << " -- Sample " << ProcessNames[i] << " xs " << getXS(ProcessNames[i]) << " weight " << weight << std::endl;
		
		std::ostringstream lcut;
		lcut << weight << "*(" << selection[I] << ")";
		//doesn't work for vetos...
		//lcut << weight << "*TMath::Power(0,ntaus)*TMath::Power(0,nmediumbjets)*(" << selection[I] << ")";
		//syst: *(1+2.5*nvetotaus)

		if (i==0)
		  {
		    //signal process      
		    int counters1  = LightTree->Draw((histName[variable]+">>SignalPlot"+PU[pu]+binRange[variable]).c_str(),lcut.str().c_str());
		    if (!do2D) dataPlot= (TH1*)gDirectory->Get(("SignalPlot"+PU[pu]).c_str());
		    else dataPlot2D= (TH2*)gDirectory->Get(("SignalPlot"+PU[pu]).c_str());
		    
		    //			for (int j = 0; j<dataPlot->GetNbinsX()+1;j++)
		    //{
		    //  dataPlot->SetBinContent(j,dataPlot->GetBinContent(j)*weight);
		    
		    //}
		    if (!do2D){
		      dataPlot->SetLineColor(kBlack);
		      double integral = dataPlot->Integral(0,dataPlot->GetNbinsX()+1);
		      if (variable==0 && pu==0) leg->AddEntry(dataPlot,"VBF H(125)","L");
		      std::cout << " Check yield: " 
			//<< 1.0*weight*LightTree->GetEntries(selection[I].c_str()) << " " 
				<< dataPlot->Integral() << " " << integral << std::endl;
		      if (variable==0) {
			Process_Yields[VBFH]+= integral;
			double eff = integral/weight*1./nEvtMap[ProcessNames[i]];
			Process_Yields_err[VBFH]+= eff*(1-eff)*nEvtMap[ProcessNames[i]]*weight*weight;
		      }
		    }
		    if (doJes>0) return 1;
		  }
		else
		  {
		    //backgrounds     
		    int counters2  = LightTree->Draw((histName[variable]+">>"+ProcessNames[i]+"_"+PU[pu]+binRange[variable]).c_str(),lcut.str().c_str());
		    if (!do2D) histPlot= (TH1*)gDirectory->Get((ProcessNames[i]+"_"+PU[pu]).c_str());
		    else histPlot2D= (TH2*)gDirectory->Get((ProcessNames[i]+"_"+PU[pu]).c_str());
		    
		    if (!do2D){
		      double integral = histPlot->Integral(0,histPlot->GetNbinsX()+1);
		      double eff = (weight>0) ? integral/weight*1./nEvtMap[ProcessNames[i]]: 0;
		      double interrsq = eff*(1-eff)*nEvtMap[ProcessNames[i]]*weight*weight;
		      double ytot = (weight>0) ? integral/weight : 0;
		      double yee = 0,ymumu=0,ytautau=0;
		      double yenu=0,ymunu=0,ytaunu=0;
		      if (ProcessNames[i].find("ToLL")!=ProcessNames[i].npos){
			yee = LightTree->GetEntries((selection[I]+"&& TMath::Abs(GenLep1_pdgid)==11 &&  TMath::Abs(GenLep2_pdgid)==11").c_str());
			ymumu = LightTree->GetEntries((selection[I]+"&& TMath::Abs(GenLep1_pdgid)==13 &&  TMath::Abs(GenLep2_pdgid)==13").c_str());
			ytautau = LightTree->GetEntries((selection[I]+"&& TMath::Abs(GenLep1_pdgid)==15 &&  TMath::Abs(GenLep2_pdgid)==15").c_str());
			//rescale by total
			double yemutau = yee+ymumu+ytautau;
			yee *= yemutau>0 ? ytot/yemutau : 0;
			ymumu *= yemutau>0 ? ytot/yemutau : 0;
			ytautau *= yemutau>0 ? ytot/yemutau : 0;
		      } else if (ProcessNames[i].find("EWKW")!=ProcessNames[i].npos || ProcessNames[i].find("JetsToLNu")!=ProcessNames[i].npos){
			yenu = LightTree->GetEntries((selection[I]+"&& (TMath::Abs(GenLep1_pdgid)==11 || (TMath::Abs(GenLep1_pdgid)==15 && TMath::Abs(GenLep2_pdgid)==11))").c_str());
			ymunu = LightTree->GetEntries((selection[I]+"&& (TMath::Abs(GenLep1_pdgid)==13 || (TMath::Abs(GenLep1_pdgid)==15 && TMath::Abs(GenLep2_pdgid)==13))").c_str());
			ytaunu = LightTree->GetEntries((selection[I]+"&& TMath::Abs(GenLep1_pdgid)==15 &&  TMath::Abs(GenLep2_pdgid)==0").c_str());
			//rescale by total
			double yemutau = yenu+ymunu+ytaunu;
			yenu *= yemutau>0 ? ytot/yemutau : 0;
			ymunu *= yemutau>0 ? ytot/yemutau : 0;
			ytaunu *= yemutau>0 ? ytot/yemutau : 0;
		      }

		      if (variable==0) std::cout << " Check yield: " 
					 //<< 1.0*weight*LightTree->GetEntries(selection[I].c_str()) << " " 
					 //<< histPlot->Integral() 
						 << " " << integral 
						 << " ee+mm+tt " <<(yee+ymumu+ytautau)*weight
						 << " e+m+t " << (yenu+ymunu+ytaunu)*weight
						 << std::endl;
		      
		      if (ProcessNames[i].find("QCD")!=ProcessNames[i].npos){
			histPlot->SetFillColor(kYellow);
			histPlot->SetLineColor(kYellow);
			if (variable==0 && pu==0) leg->AddEntry(histPlot,"QCD","F");
			if (variable==0) {
			  Process_Yields[QCD]+= integral;
			  Process_Yields_err[QCD]+= interrsq;
			}
		      }
		      
		      else if (ProcessNames[i].find("EWKZ2Jets_ZToLL")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#91ABC4"));
			histPlot->SetLineColor(TColor::GetColor("#91ABC4"));
			if (variable==0 && pu==0 && I!=1) leg->AddEntry(histPlot,"EWK Zll","F");
			if (variable==0) {
			  Process_Yields[EWK_Zll]+= integral;
			  Process_Yields[EWKZee] += weight*yee;
			  Process_Yields[EWKZmumu] += weight*ymumu;
			  Process_Yields[EWKZtautau] += weight*ytautau;
			  Process_Yields_err[EWK_Zll]+= interrsq;
			  Process_Yields_err[EWKZee] += ytot>0? interrsq*yee/ytot : 0;
			  Process_Yields_err[EWKZmumu] += ytot>0? interrsq*ymumu/ytot : 0;
			  Process_Yields_err[EWKZtautau] += ytot>0? interrsq*ytautau/ytot : 0;
			  
			}
		      }
		      
		      else if (ProcessNames[i].find("DYJetsToLL")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#9A9EAB"));
			histPlot->SetLineColor(TColor::GetColor("#9A9EAB"));
			if (ProcessNames[i].find("HT-70")!=ProcessNames[i].npos && variable==0 && pu==0 && I!=1) leg->AddEntry(histPlot,"QCD Zll","F");
			if (variable==0) {
			  Process_Yields[QCD_Zll]+= integral;
			  Process_Yields[QCDZee] += weight*yee;
			  Process_Yields[QCDZmumu] += weight*ymumu;
			  Process_Yields[QCDZtautau] += weight*ytautau;
			  Process_Yields_err[QCD_Zll]+= interrsq;
			  Process_Yields_err[QCDZee] += ytot>0? interrsq*yee/ytot : 0;
			  Process_Yields_err[QCDZmumu] += ytot>0? interrsq*ymumu/ytot : 0;
			  Process_Yields_err[QCDZtautau] += ytot>0? interrsq*ytautau/ytot : 0;
			}
			
		      }
		      
		      else if (ProcessNames[i].find("ST")!=ProcessNames[i].npos || ProcessNames[i].find("TT")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#CF3721"));
			histPlot->SetLineColor(TColor::GetColor("#CF3721"));
			if (ProcessNames[i].find("ST_s-channel")!=ProcessNames[i].npos && variable==0 && pu==0) leg->AddEntry(histPlot,"Top","F");
			if (variable==0) {
			  Process_Yields[Top]+= integral;
			  Process_Yields_err[Top]+= interrsq;
			}
		      }
		      
		      else if (ProcessNames[i].find("EWKW")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#0066CC"));
			histPlot->SetLineColor(TColor::GetColor("#0066CC"));
			if (ProcessNames[i].find("Plus")!=ProcessNames[i].npos && variable==0 && pu==0) leg->AddEntry(histPlot,"EWK W+jets#rightarrow l#nu","F");
			if (variable==0) {
			  Process_Yields[EWK_W_lnu]+= integral;
			  Process_Yields[EWKWenu] += weight*yenu;
			  Process_Yields[EWKWmunu] += weight*ymunu;
			  Process_Yields[EWKWtaunu] += weight*ytaunu;
			  Process_Yields_err[EWK_W_lnu]+= interrsq;
			  Process_Yields_err[EWKWenu] += ytot>0? interrsq*yenu/ytot : 0;
			  Process_Yields_err[EWKWmunu] += ytot>0? interrsq*ymunu/ytot : 0;
			  Process_Yields_err[EWKWtaunu] += ytot>0? interrsq*ytaunu/ytot : 0;
			}
		      }
		      else if (ProcessNames[i].find("EWKZ2Jets_ZToNuNu")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#00CCCC"));
			histPlot->SetLineColor(TColor::GetColor("#00CCCC"));
			if (variable==0 && pu==0) leg->AddEntry(histPlot,"EWK Z#rightarrow#nu#nu","F");
			if (variable==0) {
			  Process_Yields[EWK_Z_nunu]+= integral;
			  Process_Yields_err[EWK_Z_nunu]+= interrsq;
			}
		      }
		      else if (ProcessNames[i].find("JetsToLNu")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#E19D07"));
			histPlot->SetLineColor(TColor::GetColor("#E19D07"));
			if (ProcessNames[i].find("0JetsToLNu")!=ProcessNames[i].npos && variable==0 && pu==0) leg->AddEntry(histPlot,"QCD W+jets#rightarrow l#nu","F");
			if (variable==0) {
			  Process_Yields[QCD_W_lnu]+= integral;
			  Process_Yields[QCDWenu] += weight*yenu;
			  Process_Yields[QCDWmunu] += weight*ymunu;
			  Process_Yields[QCDWtaunu] += weight*ytaunu;
			  Process_Yields_err[QCD_W_lnu] += interrsq;
			  Process_Yields_err[QCDWenu] += ytot>0? interrsq*yenu/ytot : 0;
			  Process_Yields_err[QCDWmunu] += ytot>0? interrsq*ymunu/ytot : 0;
			  Process_Yields_err[QCDWtaunu] += ytot>0? interrsq*ytaunu/ytot : 0;
			}
		      }
		      else if (ProcessNames[i].find("ZJetsToNuNu")!=ProcessNames[i].npos){
			histPlot->SetFillColor(TColor::GetColor("#4D975D"));
			histPlot->SetLineColor(TColor::GetColor("#4D975D"));
			if (ProcessNames[i].find("100To200")!=ProcessNames[i].npos && variable==0 && pu==0) leg->AddEntry(histPlot,"QCD Z#rightarrow#nu#nu","F");
			if (variable==0) {
			  Process_Yields[QCD_Z_nunu]+= integral;
			  Process_Yields_err[QCD_Z_nunu]+= interrsq;
			}
		      }
		      
		      //for (int j = 0; j<histPlot->GetNbinsX()+1;j++)
		      //{
		      //  histPlot->SetBinContent(j,histPlot->GetBinContent(j)*weight);
		      
		      //}
		      SetHistStyle(histPlot,histLabelX[variable]);
		    
		      A->Add(histPlot);
		      if (i==1) background = (TH1*) histPlot->Clone("totalBkg");
		      else background->Add(histPlot);
		    } else {
		      if (i==1) background2D = (TH2*) histPlot2D->Clone("totalBkg");
		      else background2D->Add(histPlot2D);
		    }
                    
                    
		  }//if bkg
	      }//ifLT exists
	      else {
		std::cout << "LightTree not found !" << std::endl;
	      }
	    }//loop on samples
	  
	  if (!dataPlot || (!do2D && !A)) return 1;
          
	  if (!do2D){
	    //dataPlot->SetFillColor(kBlack);
	    dataPlot->SetLineColor(kBlack);
	    dataPlot->SetLineWidth(3);
	    dataPlot->SetFillStyle(0);
	    //dataPlot->SetMarkerStyle(20);
	    
	    SetHistStyle(dataPlot,histLabelX[variable]);
	    
	    //dataPlot->Draw("e1 goff");
	    c1->Clear();
	    c1->cd();
	    A->Draw("hist goff");
	    SetHistStyle((TH1*)A->GetHistogram(),histLabelX[variable]);
	    if (histName[variable].find("jetmetnolepmindphi")!=histName[variable].npos) ((TH1*)A->GetHistogram())->GetYaxis()->SetRangeUser(0,9500);
	    
	    dataPlot->Draw("same hist goff");
	    gPad->RedrawAxis();
	    leg->Draw("goff");
	    
	    CMS_lumi(c1,0,0);


	    /*char buf[100];
	    sprintf(buf,"%3.0f  fb^{-1} (14 TeV)",lumi);
	    
	    TLatex *   tex = new TLatex(0.95,0.96,buf);
	    tex->SetNDC();
	    tex->SetTextAlign(31);
	    //tex->SetTextFont(42);
	    tex->SetTextSize(0.035);
	    tex->SetLineWidth(2);
	    tex->Draw("goff");
	    tex = new TLatex(0.15,0.96,("CMS Delphes Simulation "+PU[pu]).c_str());
	    tex->SetNDC();
	    //tex->SetTextFont(52);
	    tex->SetTextSize(0.035);
	    tex->SetLineWidth(2);
	    tex->Draw("goff");*/
	  }

	  outfile->cd();
	  if (!do2D){
	    dataPlot->Write(("Signal_"+histNameShort[variable]+"_"+PU[pu]).c_str());
	    background->Write(("Background_"+histNameShort[variable]+"_"+PU[pu]).c_str());
	    c1->Update();
	    c1->Print((plotbasedir+"/plots"+Regions[I]+"/"+histNameShort[variable]+"_"+PU[pu]+"_auto.pdf").c_str());
	    c1->Print((plotbasedir+"/plots"+Regions[I]+"/"+histNameShort[variable]+"_"+PU[pu]+"_auto.png").c_str());
	    c1->Print((plotbasedir+"/plots"+Regions[I]+"/"+histNameShort[variable]+"_"+PU[pu]+"_auto.C").c_str());
	    c1->Print((plotbasedir+"/plots"+Regions[I]+"/"+histNameShort[variable]+"_"+PU[pu]+"_auto.root").c_str());
	  }
	  else {
	    dataPlot2D->Write(("Signal_"+histNameShort[variable]+"_"+PU[pu]).c_str());
	    background2D->Write(("Background_"+histNameShort[variable]+"_"+PU[pu]).c_str());
	  }
          
	  //cout<<"save"<<std::endl;
	  A->Delete();
	  for (unsigned i = 0;i < nSamples; i++){
	    fin[i]->Close();
	  }

	}//loop on variables
      for (unsigned yields = 0;yields<nProc; yields++){
	Yield<<Yields_names[yields]<<" "<<Process_Yields[yields]<<" \\pm " << sqrt(Process_Yields_err[yields]) 
	     << " " << sqrt(Process_Yields_err[yields])/Process_Yields[yields] << std::endl;
	std::cout <<Yields_names[yields]<<" "<<Process_Yields[yields]<<" \\pm " << sqrt(Process_Yields_err[yields])
		  << " " << sqrt(Process_Yields_err[yields])/Process_Yields[yields] <<std::endl;
      }
      Yield<<"  "<<std::endl;
      std::cout <<"  "<<std::endl;

      Datacard << "rate \t" ;
      if (I==1) {
	Datacard << Process_Yields[VBFH] << "\t"
		 << Process_Yields[VBFH]*ggFfrac << "\t " 
		 << Process_Yields[QCD_Z_nunu] << "\t"
		 << Process_Yields[EWK_Z_nunu] << "\t";
      }
      
      //Factor 0.5 for We due to lower eff of medium electrons.
      double factWe = I==1 ? 0.5 : 1;
      Datacard	 << Process_Yields[QCD_Zll] << "\t"
		 << Process_Yields[EWK_Zll] << "\t"
		 << Process_Yields[QCDWmunu] << "\t"
		 << Process_Yields[EWKWmunu] << "\t"
		 << Process_Yields[QCDWenu]*factWe << "\t"
		 << Process_Yields[EWKWenu]*factWe << "\t"
		 << Process_Yields[QCDWtaunu] << "\t"
		 << Process_Yields[EWKWtaunu] << "\t"
		 << Process_Yields[Top] << "\t"
		 << Process_Yields[QCD]
		 << std::endl;

      Datacard << "------------" << std::endl;
      
      if (I==1){
	Datacard << setprecision(4) << "CMS_VBFHinv_qqH_norm \t lnN \t"
		 << 1+sqrt(Process_Yields_err[VBFH])/Process_Yields[VBFH] << "\t" 
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;

	Datacard << setprecision(4) << "CMS_VBFHinv_zvv_qcd_norm \t lnN \t" 
		 << "-" << "\t"
		 << "-" << "\t" 
		 << 1+sqrt(Process_Yields_err[QCD_Z_nunu])/Process_Yields[QCD_Z_nunu] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;


	Datacard << setprecision(4) << "CMS_VBFHinv_zvv_ewk_norm \t lnN \t" 
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[EWK_Z_nunu])/Process_Yields[EWK_Z_nunu] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;
      }

      if (I==3 || I==5){
	if (I==3) Datacard << setprecision(4) << "CMS_VBFHinv_zee_qcd_norm \t lnN \t" ;
	else Datacard << setprecision(4) << "CMS_VBFHinv_zmumu_qcd_norm \t lnN \t" ;
	Datacard << 1+sqrt(Process_Yields_err[QCD_Zll])/Process_Yields[QCD_Zll] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;


	if (I==3) Datacard << setprecision(4) << "CMS_VBFHinv_zee_ewk_norm \t lnN \t" ;
	else Datacard << setprecision(4) << "CMS_VBFHinv_zmumu_ewk_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[EWK_Zll])/Process_Yields[EWK_Zll] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;
      }

      if (I==1 || I==4 || I==5){
	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_wmu_qcd_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else if (I==4) Datacard << setprecision(4) << "CMS_VBFHinv_WMCR_wmu_qcd_norm \t lnN \t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_ZMMCR_wmu_qcd_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[QCDWmunu])/Process_Yields[QCDWmunu] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;

	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_wmu_ewk_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else if (I==4) Datacard << setprecision(4) << "CMS_VBFHinv_WMCR_wmu_ewk_norm \t lnN \t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_ZMMCR_wmu_ewk_norm \t lnN \t";
	Datacard<< "-" << "\t"
		<< "-" << "\t"
		<< "-" << "\t"
		<< 1+sqrt(Process_Yields_err[EWKWmunu])/Process_Yields[EWKWmunu] << "\t"
		<< "-" << "\t"
		<< "-" << "\t"
		<< "-" << "\t"
		<< "-" << "\t"
		<< "-" << "\t"
		<< "-" 
		<< std::endl;
      }
      if (I==1 || I==2 || I==3){

	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_wel_qcd_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else if (I==2) Datacard << setprecision(4) << "CMS_VBFHinv_WECR_wel_qcd_norm \t lnN \t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_ZEECR_wel_qcd_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[QCDWenu])/Process_Yields[QCDWenu] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;

	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_wel_ewk_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else if (I==2) Datacard << setprecision(4) << "CMS_VBFHinv_WECR_wel_ewk_norm \t lnN \t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_ZEECR_wel_ewk_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[EWKWenu])/Process_Yields[EWKWenu] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;
      }
      if (I==1 || I==2){
	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_wtau_qcd_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_WECR_wtau_qcd_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[QCDWtaunu])/Process_Yields[QCDWtaunu] << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;
      
	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_wtau_ewk_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_WECR_wtau_ewk_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[EWKWtaunu])/Process_Yields[EWKWtaunu] << "\t"
		 << "-" << "\t"
		 << "-" 
		 << std::endl;
      }

      if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_top_norm \t lnN \t"
			 << "-" << "\t"
			 << "-" << "\t"
			 << "-" << "\t"
			 << "-" << "\t";
      else if (I==2) Datacard << setprecision(4) << "CMS_VBFHinv_WECR_top_norm \t lnN \t";
      else if (I==3) Datacard << setprecision(4) << "CMS_VBFHinv_ZEECR_top_norm \t lnN \t";
      else if (I==4) Datacard << setprecision(4) << "CMS_VBFHinv_WMCR_top_norm \t lnN \t";
      else Datacard << setprecision(4) << "CMS_VBFHinv_ZMMCR_top_norm \t lnN \t";
      Datacard << "-" << "\t"
	       << "-" << "\t"
	       << "-" << "\t"
	       << "-" << "\t"
	       << "-" << "\t"
	       << "-" << "\t"
	       << "-" << "\t"
	       << "-" << "\t"
	       << 1+sqrt(Process_Yields_err[Top])/Process_Yields[Top] << "\t"
	       << "-" 
	       << std::endl;
      
      if (I==1 || I==2 || I==4){
	if (I==1) Datacard << setprecision(4) << "CMS_VBFHinv_SR_qcd_norm \t lnN \t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t"
			   << "-" << "\t";
	else if (I==2) Datacard << setprecision(4) << "CMS_VBFHinv_WECR_qcd_norm \t lnN \t";
	else Datacard << setprecision(4) << "CMS_VBFHinv_WMCR_qcd_norm \t lnN \t";
	Datacard << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << "-" << "\t"
		 << 1+sqrt(Process_Yields_err[QCD])/Process_Yields[QCD]
		 << std::endl;
      }


    }//loop on PU
  Yield.close();
  Datacard.close();
  outfile->Close();

  return 0;
 
}//getRegion method

int plotterDelphes(){
  //std::string Regions[nsel] = {"","0e0mu","1e0mu","2e0mu","0e1mu","0e2mu"};

  int errcode=0;

  std::string prodDate = "180809";//default
  //std::string prodDate = "180917";//smear45%
  //std::string prodDate = "181008";//run2 MET
  std::ostringstream plotbasedir;

  //plotbasedir = "Mjj1000deta4mht170dphi18";
  //ggF fraction = 0.3
  //for (int regions = 1; regions<6; regions++){
  //errcode=getRegion(regions,prodDate,plotbasedir,170,0,1000,4,1.8,0.3);
  //}
  const unsigned nMet = 1;
  //int met[nMet]={150,180,190,200,210,220,250,300,350,400};
  //int met[nMet]={360,380};
  int met[nMet]={190};
  const unsigned nMht = 1;
  //int mht[nMht]={160,200,250,300,350};
  //int mht[nMht]={160,200,300,400,500,600,700,800};
  int mht[nMht]={100};
  const unsigned nMjj = 1;//16;
  int mjj[nMjj] = {2500};
  //int mjj[nMjj] = {
  //1000,1500,
  //2000,2200,2400,2500,2600,2800,
  //3000,3500,
  //4000};



  for (unsigned im(0); im<nMet;++im){//loop on met
    for (unsigned in(0); in<nMht;++in){//loop on mht
      for (unsigned jm(0); jm<nMjj; ++jm){//loop on mjj
	//if (met[im]!=150 && met[im]!=170 && met[im]!=210 && mjj[jm]<1850) continue;
	plotbasedir.str("");
	plotbasedir << "Mjj" << mjj[jm] << "deta4met" << met[im] << "dphi18mht" << mht[in];
      //plotbasedir << "Mjj" << mjj[jm] << "deta4met0dphi18mht" << mht[im];
      //ggF fraction = 0.25
	if (system(("./prepareDirs.sh "+plotbasedir.str()).c_str())) return 1;

	float ggRatio = getggRatio(mjj[jm],met[im],mht[in]);
	std::cout << " ggRatio: " << ggRatio << std::endl;
	for (int regions = 1; regions<2; regions++){
	  errcode=getRegion(regions,prodDate,plotbasedir.str(),mht[in],met[im],mjj[jm],4,1.8,ggRatio,true);
	}

      }//loop on mjj
    }//loop on mht
  }//loop on met
  
  return errcode;
  
}
