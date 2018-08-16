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

float getggRatio(const float mjj, const float met){

  TFile *file = TFile::Open("ggHoverVBFvsMjjvsMHT.root");
  file->cd();
  TH2F *htmp = (TH2F*)gDirectory->Get("hr1");


  int binx = htmp->GetXaxis()->FindBin(met+0.1);
  int biny = htmp->GetYaxis()->FindBin(mjj+0.1);

  float result = htmp->GetBinContent(binx,biny);

  std::cout << " Check: " << met << " " << mjj << " " << result << std::endl;

  return result;
}

