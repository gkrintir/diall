
//================================================================================================
//
// Perform fits to recoil against Z->ll events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <list>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include <math.h>
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"

#include "TLorentzVector.h"

#include "TClass.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooVoigtian.h>
#include <RooGaussian.h>
#include <RooHistPdf.h>
#include <RooFormulaVar.h>

#include <unordered_map>

using namespace RooFit;
using namespace std;                       


struct FileInfo
{
    std::string nEvents;
    std::string crossSection;
};

const std::unordered_map<std::string,FileInfo> eventWeights = {
  {"uParaZllPt",
        {
	    "try",
            "try"
        }
    },
};

RooRealVar x ("x", "x", -400, 400); // changed the axis range, we only need 800 for zll_pt going to 500 GeV.
RooRealVar g_w ("g_w", "width Gaus", 10., 0., 100., "GeV");	//40
RooRealVar gamma_Z0 ("gamma_Z0_U", "Z0 width", 2.3, 0, 100, "GeV");	//20
RooRealVar v_m ("v_m", "v_m",0,-10.,10.);

RooVoigtian *voigt;
RooFitResult *result;
RooAddPdf *model;

double f;
double efwhm;
//TCanvas *c1 = new TCanvas ("c1", "c1", 800, 800);

//--------------------------------------------------------------------------------------------------
Double_t sigmaFunc(Double_t *x, Double_t *par) {
  // par[0]: quadratic coefficient
  // par[1]: linear coefficient
  // par[2]: constant term
  
  Double_t a  = par[0];
  Double_t b  = par[1];
  Double_t c  = par[2];
  Double_t d  = par[3];
    
  return a*x[0]*x[0] + b*x[0] + c;
}

//--------------------------------------------------------------------------------------------------
// function to describe relative fraction in a double Gaussian based on 
// functions for sigma0, sigma1, and sigma2
Double_t frac2Func(Double_t *x, Double_t *par) {
  // par[0..3]:  sigma0
  // par[4..7]:  sigma1
  // par[8..11]: sigma2
  
  TF1 s0("_s0",sigmaFunc,0,7000,4); s0.SetParameters(par[0],par[1],par[2],par[3]);
  TF1 s1("_s1",sigmaFunc,0,7000,4); s1.SetParameters(par[4],par[5],par[6],par[7]);
  TF1 s2("_s2",sigmaFunc,0,7000,4); s2.SetParameters(par[8],par[9],par[10],par[11]);
  
  return (s0.Eval(x[0]) - s1.Eval(x[0]))/(s2.Eval(x[0]) - s1.Eval(x[0]));
}

/*
//--------------------------------------------------------------------------------------------------
Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[2];
  df[0] = 1;
  df[1] = x;
  Double_t err2 = df[0]*df[0]*(fs->GetCovarianceMatrix()[0][0]) 
                  + df[1]*df[1]*(fs->GetCovarianceMatrix()[1][1]) 
		  + 2.0*df[0]*df[1]*(fs->GetCovarianceMatrix()[0][1]);
  assert(err2>=0);
  return sqrt(err2);
}

//--------------------------------------------------------------------------------------------------
Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[4];
  Double_t a  = fcn->GetParameter(0);
  Double_t b  = fcn->GetParameter(1);
  Double_t c  = fcn->GetParameter(2);
  
  df[0] = x*x;
  df[1] = x;
  df[2] = 1;
  
  Double_t err2=0;
  for(Int_t i=0; i<3; i++) {
    err2 += df[i]*df[i]*(fs->GetCovarianceMatrix()[i][i]);
    for(Int_t j=i+1; j<3; j++) {
      err2 += 2.0*df[i]*df[j]*(fs->GetCovarianceMatrix()[i][j]);
    }
  }
  assert(err2>=0);
  return sqrt(err2);
}
*/
// perform fit of recoil component
void performFit(vector<TH1*> hv, vector<TH1*> hbkgv, 
                const Int_t model, const Bool_t sigOnly,
                const char *plabel, const char *xlabel,
                Double_t *meanArr,   Double_t *meanErrArr,
                Double_t *sigma0Arr, Double_t *sigma0ErrArr,
                Double_t *sigma1Arr, Double_t *sigma1ErrArr,
                Double_t *sigma2Arr, Double_t *sigma2ErrArr,
                Double_t *sigma3Arr, Double_t *sigma3ErrArr,
                Double_t *frac2Arr,  Double_t *frac2ErrArr,
                Double_t *frac3Arr,  Double_t *frac3ErrArr, TString outputDir) ;

void 
identifyToken (std::set<char> & token, const std::string& tokenize ) ;

void 
removeWhiteSpaces( std::string &strg ) ;

/*
void
ConstructModel(RooDataHist Hist,RooDataHist *bkg_hist, bool BKGSubtract) {

  f=0;
  efwhm=0;

  v_m.setVal(Hist.mean(x) );
  v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));

  voigt =new RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w);
  
  
  if(BKGSubtract) {
    RooHistPdf *bkg_pdf = new RooHistPdf("bkg_pdf","bkg_pdf",RooArgSet(x),*bkg_hist);
    RooRealVar *lAbkgFrac =new RooRealVar("AbkgFrac","AbkgFrac",0.5,0.,1.);
    RooFormulaVar * sigbkgFrac= new RooFormulaVar("bkgfrac","@0",RooArgSet(*lAbkgFrac));
    model = new RooAddPdf("modelSB","modelSB",*voigt,*bkg_pdf,*sigbkgFrac);
    result = model->fitTo (Hist, RooFit::Minimizer("Minuit2","migrad"),RooFit::Strategy(2), RooFit::SumW2Error (kFALSE), RooFit::Save (kTRUE), RooFit::PrintLevel (-1));	// -1 verbose
                         
  } else {
      result = voigt->fitTo (Hist, RooFit::Minimizer("Minuit2","migrad"),RooFit::Strategy(2), RooFit::SumW2Error (kFALSE), RooFit::Save (kTRUE), RooFit::PrintLevel (-1));	// -1 verbose //
  }

  //if(result->status()!=0) voigt=0;



  //Get the FWHM
  double sigma = g_w.getVal ();
  double gamma = gamma_Z0.getVal ();
  double esigma = g_w.getError ();
  double egamma = gamma_Z0.getError ();

  double Vsg = result->correlation (g_w, gamma_Z0);
  double Vgs = result->correlation (gamma_Z0, g_w);
  double Vss = result->correlation (g_w, g_w);
  double Vgg = result->correlation (gamma_Z0, gamma_Z0);
  cout << "correlacion Vgs " << Vgs << " y correlacion Vsg" << Vsg << endl;
  f = FWHM (sigma, gamma);
  efwhm = FWHMError_fixed (sigma, gamma, esigma, egamma, Vss, Vsg, Vgs, Vgg);

  return;

}
*/

  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;

  Double_t ptbins[] = {30,40,50,60,70,80,90,100};
  Int_t nbins = sizeof(ptbins)/sizeof(Double_t)-1;
  //Int_t nbins = 1 ;

//=== MAIN MACRO ================================================================================================= 

void fitRecoil(TString Datafilename="/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/AnaResults.root",  // input ntuple for data
	       TString Bkgfilename="/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/out/AnaResults_1.root",  // input ntuple for bkg
	       Int_t   pfu1model=1,   // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
	       TString outputDir="./" // output directory
	       )
{

  TFile *file;

  //Data
  file = new TFile(Datafilename);
  if (file->IsZombie()) gSystem->Abort();
  TList * metPFRawData = dynamic_cast<TList*>(file->Get("metPFRaw"));
  TTree * metPerformanceInfoData = dynamic_cast<TTree *>(metPFRawData->FindObject("METPerformanceInf"));

  //Bkg.
  Bool_t sigOnly=0; // to be added in the config. file !!

  file = new TFile(Bkgfilename);
  if (file->IsZombie()) gSystem->Abort();
  TList * metPFRawBgr = dynamic_cast<TList*>(file->Get("metPFRaw_ee"));
  TTree * metPerformanceInfoBgr = dynamic_cast<TTree *>(metPFRawBgr->FindObject("METPerformanceInf"));

  
  TString condition="(xsec)*(puWeight)";//"(weighttotal)*(channel=="+dileptonch +")"; //to be done!!

  //Init variables to read from the input TTrees
  RooRealVar MassZCands("MassZCands", "MassZCands", 0, 1000);
  RooRealVar PtZCands("PtZCands", "PtZCands", 0, 200);
  RooRealVar EtaZCands("EtaZCands", "EtaZCands", -5., 5.);
  RooRealVar uParaZllPt("uParaZllPt", "uParaZllPt", -50, 50);

  //std::string cut("PtZCands>40 && PtZCands<60");
  std::map<std::string, std::string > matEffHistosToPdgIdsToSelect;
  //histos1DToPdgIdsToSelect_[pdgNames[iPart]].push_back(histos1D);

  ifstream fin ("Config.dat");
  //fin.open( "Config.dat" );
  std::string line, val;
  //fin >> name >> val;
  if (!fin.is_open())
    cout<<"Could not open file\n";
  else {
    do {
      getline(fin,line);
      int count = 0;
      if (fin) {
        stringstream s(line); 
        while (s.good()) { 
          if (s >> val){
	    if (count==1)
	      std::cout<<val.data()<<std::endl;
	    count++;
	  }
        } 
      } 
    } while (fin);
    fin.close();
  }
  /*
  std::vector<std::string> vctr_particles;
  std::set<char> token;
  identifyToken(token, val);
  std::istringstream particles_toDraw(val);
  std::string s_particle;
  while (std::getline(particles_toDraw, s_particle, (*token.begin()))) {
    removeWhiteSpaces(s_particle);
    vctr_particles.push_back(s_particle);
    std::cout<<s_particle.data()<<std::endl;
  }
  */
  /*
  fNPar = (int) val;
  for(int i=0; i!=fNPar; ++i) {
    fin >> name;
    fin >> val;
    fin >> err;
    fParNam[i] = name;
    fParVal[i] = val;
    fParErr[i] = err;
  }
  */
  fin.close();

  for (int i =0 ; i<1; i++) {
    std::vector<TH1*> hPFu1v, hPFu1Bkgv;
    for(Int_t ibin=0; ibin<nbins; ibin++) {
      TString cut(Form("PtZCands>%f && PtZCands<%f", ptbins[ibin], ptbins[ibin+1]));
      RooDataSet *data = new RooDataSet(Form("data_%d", ibin), "data", RooArgSet(uParaZllPt, PtZCands), Import(*metPerformanceInfoData), Cut(cut.Data()));
      RooDataSet *bkg = new RooDataSet(Form("bkg_%d", ibin), "bkg", RooArgSet(uParaZllPt, PtZCands), Import(*metPerformanceInfoBgr), Cut(cut.Data()));
      RooDataHist data_binned(Form("data_%d", ibin),"binned version of data",RooArgSet(uParaZllPt), *data) ;
      //TH1 *h_data = data_binned.createHistogram("uParaZllPt", 20);
      hPFu1v.push_back(data_binned.createHistogram("uParaZllPt", 20));   
      
      RooDataHist bkg_binned(Form("bkg_%d", ibin),"binned version of bkg",RooArgSet(uParaZllPt), *bkg) ;
      //TH1 *h_bkg = bkg_binned.createHistogram("uParaZllPt", 20);
      hPFu1Bkgv.push_back( bkg_binned.createHistogram("uParaZllPt", 20));   
      
    }
    
    std::string processName = "uParaZllPt";
    auto it = eventWeights.find(processName);
    if (it==eventWeights.end())
      {
	gSystem->Error("", "no event weight information available for process name '") ;
	gSystem->Abort();
      }
    else
      {
	printf("%s", it->second.crossSection.data());
      }
    
    std::cout<<hPFu1Bkgv.size()<<std::endl;
    
    TString formulaPFu1mean("pol1");
    
    TFitResultPtr fitresPFu1mean;   TF1 *fcnPFu1mean   = new TF1("fcnPFu1mean",formulaPFu1mean,0,7000);
    TFitResultPtr fitresPFu1sigma1; TF1 *fcnPFu1sigma1 = new TF1("fcnPFu1sigma1",sigmaFunc,0,7000,3);
    TFitResultPtr fitresPFu1sigma2; TF1 *fcnPFu1sigma2 = new TF1("fcnPFu1sigma2",sigmaFunc,0,7000,3);  
    TFitResultPtr fitresPFu1sigma0; TF1 *fcnPFu1sigma0 = new TF1("fcnPFu1sigma0",sigmaFunc,0,7000,3);    
    TFitResultPtr fitresPFu1frac2;  TF1 *fcnPFu1frac2  = new TF1("fcnPFu1frac2",frac2Func,0,7000,12);
    
    // Control the parameters of the fitting functions here
    
    //   fcnPFu1mean->SetParameter(0,4.47);fcnPFu1mean->SetParLimits(0,4.45,4.5);
    //   fcnPFu1mean->SetParameter(1,-0.964);fcnPFu1mean->SetParLimits(1,-1,-0.95);
    
    fcnPFu1sigma1->SetParameter(0,-5e-5); //fcnPFu1sigma1->SetParLimits(0,-5e-2,-200e-5);
    fcnPFu1sigma1->SetParameter(1,0.25); fcnPFu1sigma1->SetParLimits(1,-1,1);
    fcnPFu1sigma1->SetParameter(2,15);   fcnPFu1sigma1->SetParLimits(2,0,20);
    
    fcnPFu1sigma2->SetParameter(0,-2e-3);  //fcnPFu1sigma2->SetParLimits(0,-2e-3,-1e-5);
    fcnPFu1sigma2->SetParameter(1,0.5);  fcnPFu1sigma2->SetParLimits(1,-1,3);
    fcnPFu1sigma2->SetParameter(2,16);   fcnPFu1sigma2->SetParLimits(2,5,25);
    fcnPFu1sigma0->SetParameter(0,-2e-5);// fcnPFu1sigma0->SetParLimits(0,-1e-3,0);
    fcnPFu1sigma0->SetParameter(1,0.07);  fcnPFu1sigma0->SetParLimits(1,-1,1);
    fcnPFu1sigma0->SetParameter(2,14);     fcnPFu1sigma0->SetParLimits(2,0,25); 
    
    Double_t xval[nbins], xerr[nbins];
    for(Int_t ibin=0; ibin<nbins; ibin++) {
      xval[ibin] = 0.5*(ptbins[ibin+1]+ptbins[ibin]);
      xerr[ibin] = 0.5*(ptbins[ibin+1]-ptbins[ibin]);
    }
    
    // ------- Arrays and graphs to store fit results -----------
    TGraphErrors *grPFu1mean=0;   Double_t pfu1Mean[nbins],   pfu1MeanErr[nbins];
    TGraphErrors *grPFu1sigma0=0; Double_t pfu1Sigma0[nbins], pfu1Sigma0Err[nbins];
    TGraphErrors *grPFu1sigma1=0; Double_t pfu1Sigma1[nbins], pfu1Sigma1Err[nbins];
    TGraphErrors *grPFu1sigma2=0; Double_t pfu1Sigma2[nbins], pfu1Sigma2Err[nbins];
    TGraphErrors *grPFu1sigma3=0; Double_t pfu1Sigma3[nbins], pfu1Sigma3Err[nbins];
    TGraphErrors *grPFu1frac2=0;  Double_t pfu1Frac2[nbins],  pfu1Frac2Err[nbins];  
    TGraphErrors *grPFu1frac3=0;  Double_t pfu1Frac3[nbins],  pfu1Frac3Err[nbins];
    
    TGraphErrors *grPFu2mean=0;   Double_t pfu2Mean[nbins],   pfu2MeanErr[nbins];
    TGraphErrors *grPFu2sigma0=0; Double_t pfu2Sigma0[nbins], pfu2Sigma0Err[nbins];
    TGraphErrors *grPFu2sigma1=0; Double_t pfu2Sigma1[nbins], pfu2Sigma1Err[nbins];
    TGraphErrors *grPFu2sigma2=0; Double_t pfu2Sigma2[nbins], pfu2Sigma2Err[nbins];
    TGraphErrors *grPFu2sigma3=0; Double_t pfu2Sigma3[nbins], pfu2Sigma3Err[nbins];
    TGraphErrors *grPFu2frac2=0;  Double_t pfu2Frac2[nbins],  pfu2Frac2Err[nbins];  
    TGraphErrors *grPFu2frac3=0;  Double_t pfu2Frac3[nbins],  pfu2Frac3Err[nbins];
    
    //TCanvas *c = new TCanvas("c","c",800,600);
    
    performFit(hPFu1v, hPFu1Bkgv, pfu1model, sigOnly,
	       "pfu1", "PF u_{1} [GeV]",
	       pfu1Mean,   pfu1MeanErr,
	       pfu1Sigma0, pfu1Sigma0Err,
	       pfu1Sigma1, pfu1Sigma1Err,
	       pfu1Sigma2, pfu1Sigma2Err,
	       pfu1Sigma3, pfu1Sigma3Err,
	       pfu1Frac2,  pfu1Frac2Err,
	       pfu1Frac3,  pfu1Frac3Err, outputDir);
    
    
    
    //--------------------------------------------------------------------------------------------------------------
    // Make plots, requires external utility to draw plots
    // I commented out the stuff that has external dependence for now.
    // If you want plots to print, which is helpful for checking the fitted functions,
    // it needs to be updated to work without it.
    //==============================================================================================================  
    
    char fitparam[100];
    char chi2ndf[50]; 
    
    
    TPaveText *tb = new TPaveText(0.65,0.87,0.95,0.65,"NDC");
    tb->SetTextColor(kBlack);
    tb->SetFillStyle(0);
    tb->SetBorderSize(0);
    
    // Plotting u1 vs. dilepton pT
    grPFu1mean = new TGraphErrors(nbins,xval,pfu1Mean,xerr,pfu1MeanErr);
    grPFu1mean->SetName("grPFu1mean");
    fitresPFu1mean = grPFu1mean->Fit("fcnPFu1mean","QMRN0FBSE");
    sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu1mean->GetChisquare())/(fcnPFu1mean->GetNDF()));
    tb->AddText(chi2ndf);
    
    grPFu1mean->SetTitle("");
    grPFu1mean->GetXaxis()->SetTitle("p_{T}(ll) [GeV/c]");
    grPFu1mean->GetYaxis()->SetTitle("#mu(u_{2}) [GeV]");
    grPFu1mean->SetMarkerColor(kBlack);
    grPFu1mean->SetMarkerStyle(kOpenCircle);
    grPFu1mean->Draw();
    fcnPFu1mean->SetLineColor(kRed);
    fcnPFu1mean->Draw("same");
    
    sprintf(fitparam,"p_{0} = %.3f #pm %.3f",fcnPFu1mean->GetParameter(0),fcnPFu1mean->GetParError(0)); tb->AddText(fitparam);
    sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu1mean->GetParameter(1),fcnPFu1mean->GetParError(1)); tb->AddText(fitparam);
    sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu1mean->GetParameter(2),fcnPFu1mean->GetParError(2)); tb->AddText(fitparam);
    tb->Draw("same");
    
    //c->SaveAs(TString(outputDir+"/plots/"+"pfu1mean.png"));
    //c->Clear();
    //tb->Clear();
    
    /*
    // u1, sigmas
    grPFu1sigma1 = new TGraphErrors(nbins,xval,pfu1Sigma1,xerr,pfu1Sigma1Err);  
    grPFu1sigma1->SetName("grPFu1sigma1");
    fitresPFu1sigma1 = grPFu1sigma1->Fit("fcnPFu1sigma1","QMRN0SE");
    sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu1sigma1->GetChisquare())/(fcnPFu1sigma1->GetNDF())); tb->AddText(chi2ndf);
    tb->SetX1NDC(0.21);
    tb->SetX2NDC(0.51);
    grPFu1sigma1->SetTitle("");
    grPFu1sigma1->GetXaxis()->SetTitle("p_{T}(ll) [GeV/c]");
    grPFu1sigma1->GetYaxis()->SetTitle("#sigma_{1}(u_{1}) [GeV]");
    grPFu1sigma1->SetMarkerColor(kBlack);
    grPFu1sigma1->SetMarkerStyle(kOpenCircle);
    grPFu1sigma1->Draw();
    fcnPFu1sigma1->SetLineColor(kRed);
    fcnPFu1sigma1->Draw("same");
    sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu1sigma1->GetParameter(0)),1e5*(fcnPFu1sigma1->GetParError(0))); tb->AddText(fitparam);
    sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu1sigma1->GetParameter(1),fcnPFu1sigma1->GetParError(1)); tb->AddText(fitparam);
    sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu1sigma1->GetParameter(2),fcnPFu1sigma1->GetParError(2)); tb->AddText(fitparam);
    tb->Draw("same");
    //c->SaveAs(TString(outputDir+"/plots/"+"pfu1sigma1.png"));
    //c->Clear();
    tb->Clear();
    */
    /*
      TString samplephys14; TString variablename; TString xvariable; TString tchannel; bool drawchi2=false; bool WantBKGSubtract=false;
      TString variablenamepng=variablename;
      variablenamepng.ReplaceAll("/","over");
      
      TString DestFolder; // Different folder for background subtracted or not subtracted files
      if(WantBKGSubtract){
      DestFolder="BKG_Subtraction";
      }
      else {
      DestFolder="Not_BKG_Subtraction";
      }
      
      TH1::SetDefaultSumw2() ;
      
      
      cout << "sample phys " << samplephys14 << endl;
      TString folder = "DY";
  if (samplephys14.Contains ("TT"))
    folder = "TTbar";
  if (samplephys14.Contains ("GJet"))
    folder = "Gamma";
  if(samplephys14.Contains ("QCD"))
    folder = "QCD";
  if(samplephys14.Contains("Data"))
   folder ="Data"; 
   if (samplephys14.Contains("Pseudo"))
   folder="Pseudo";
 

   cout << " folder   " << folder << "  -   " << "Destfolder " << DestFolder << endl;
  
  
  TString titley = "";
  if (variablename == "pfmetx")
    titley = "#sigma(MET_{x}) GeV";
  if (variablename == "pfmety")
    titley = "#sigma(MET_{y}) GeV";




  TCanvas *c1 = new TCanvas ("c1", "c1", 800, 800);

  c1->SetTickx ();
  c1->SetTicky ();
  c1->SetFillColor (kWhite);
  c1->SetFillStyle (0);
  c1->SetRightMargin (0.05);
  c1->SetTopMargin (0.08);




  TFile filephys14 (samplephys14);



  TTree *treephys14 = (TTree *) filephys14.Get ("METtree");




  std::vector < TH1F * >resolution;
  resolution.clear ();




  int limitdown (0), limitup (0);
  TString strlimitup = "0";
  TString strlimitdown = "0";

  int tempsizexarray = 0;

  if (xvariable == "nVert")
    tempsizexarray = 6;
  if (xvariable == "met_sumEt")
    tempsizexarray = 6;
  if (xvariable == "zll_pt")
    tempsizexarray = 10; //6


  const int sizexarray=tempsizexarray;
  TH1F *histonvertex=new TH1F("histonvertex","histonvertex",50,0,50); // added later
  TH1F *histozll_pt=new TH1F("histozll_pt","histozll_pt",100,0,1200); // added later
  TH1F *histomet_sumEt=new TH1F("histomet_sumEt","histomet_sumEt",100,0,4); // added later
  TH1F *histomet_uPara_zllzll_pt=new TH1F("histomet_uPara_zllzll_pt","histomet_uPara_zllzll_pt",50,-300,300); // added later
  TH1F *histomet_uPerp_zll=new TH1F("histomet_uPerp_zll","histomet_uPerp_zll",50,-300,300); // added later
  //TH1F *met_uPara_zllresponse1=new TH1F("met_uPara_zllresponse1","met_uPara_zll",50,-300,300); // added later
  //TH1F *zll_ptresponse1=new TH1F("zll_ptresponse1","zll_pt",100,0,100); // added later
  TH1F *histoscale=new TH1F("histoscale","histoscale",50,-50,50); // added later
  
  
  TString dileptonch="";
  if (tchannel=="MuMu")  dileptonch="1";
  if (tchannel=="EE" || tchannel=="Gamma") dileptonch="0";
  
  
  // Plot inclusive distributions of the main variables
  TString condition="(xsec)*(puWeight)";//"(weighttotal)*(channel=="+dileptonch +")";
  if (samplephys14.Contains("Data")) condition="";
  //cout << condition.Data() << endl;
  if (!samplephys14.Contains("raw") && !samplephys14.Contains("jes_")){ treephys14->Draw ("nVert >> histonvertex", condition.Data());
  treephys14->Draw ("zll_pt >> histozll_pt", condition.Data());
  treephys14->Draw ("met_sumEt/1000 >> histomet_sumEt", condition.Data());
  treephys14->Draw("met_uPara_zll/zll_pt >> histoscale" ) ;
  treephys14->Draw ("met_uPara_zll+zll_pt >> histomet_uPara_zllzll_pt", condition.Data());
  treephys14->Draw ("met_uPerp_zll >> histomet_uPerp_zll",condition.Data());//condition.Data());//, condition.Data());
   }
          
  histonvertex->Draw();
  histonvertex->GetXaxis ()->SetTitle ("Number of Vertices");
  if (!samplephys14.Contains("jes_"))  c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/nVert_inclusive_.png");
  c1->SetLogy();
  histozll_pt->Draw();
  histozll_pt->GetXaxis ()->SetTitle ("zll_pt [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/"+ folder + "/"+tchannel +"/zll_pt_inclusive_.png");
  histomet_sumEt->Draw();
  histomet_sumEt->GetXaxis ()->SetTitle("sumE_{T} [TeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_sumEt_inclusive_.png");
  histomet_uPara_zllzll_pt->Draw();
  histomet_uPara_zllzll_pt->GetXaxis ()->SetTitle ("u_{||}+zll_pt [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/"+ folder + "/"+tchannel +"/met_uPara_zllzll_pt_inclusive_.png");
  histomet_uPerp_zll->Draw();
  histomet_uPerp_zll->GetXaxis ()->SetTitle ("u_{#perp}   [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_uPerp_zll_inclusive_.png");
  c1->SetLogy(0);
  histoscale->Draw();
  histoscale->GetXaxis ()->SetTitle ("u_{#perp}   [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_uPara_zlloverzll_pt.png");


  Double_t tgraphx[sizexarray], tgraphy[sizexarray], etgraphy[sizexarray],
    etgraphx[sizexarray], tgraphxchi2[sizexarray], tgraphychi2[sizexarray];
    
    



  for (int index = 0; index < sizexarray; index++)
    {
      if (xvariable == "nVert")
	limitup = (index + 1) * 5;
      if (xvariable == "met_sumEt")
	limitup = (index + 1) * 200;
      if (xvariable == "zll_pt")
	limitup = (index + 1) * 12;
      strlimitup = Form ("%d", limitup);
cout << variablename << endl;

cout << "index " << index << endl;
cout << "limit down " << limitdown << endl;
cout << "limit up " << limitup << endl;
      if(variablenamepng.Contains("over"))
        resolution.push_back (new TH1F (Form ("resx%d", index), " ", 200, -200, 200)); // changing the x axis for met_uPara_zlloverzll_pt, so that models are more visible
      else
        resolution.push_back (new TH1F (Form ("resx%d", index), " ", 200, -400, 400));
      


     TString totalnevents="1";
     TString lumi="40.3";
     TH1F * h0 = (TH1F*)filephys14.Get("Count");
     totalnevents=NToString(h0->Integral());
         
      TString conditionbkg="";
      

      if (xvariable == "nVert") condition = "(" + xvariable + "==" + strlimitup +")";
      else condition = "(" + xvariable + "<" + strlimitup + ")*(" + xvariable + ">" + strlimitdown + ")";
      
      conditionbkg=condition;
      if(!samplephys14.Contains("Data")) condition=condition+"*((xsec)*(puWeight)*("+lumi+"))/"+totalnevents;   
      
      cout << "condition data"  << condition.Data() << endl;
      treephys14->Draw (variablename + ">>" + TString (resolution[index]->GetName ()),			condition.Data(), "sames");

    //  c1->Print("~/www/"+TString (resolution[index]->GetName ())+".png");     
      //double m =  resolution[index]->GetMean (); (unused)
      //double um = resolution[index]->GetMean () - resolution[index]->GetRMS ();
      //double uM = resolution[index]->GetMean () + resolution[index]->GetRMS ();



      ////////
      
      RooDataHist Hist ("Hist", "Hist", x,			(TH1 *) resolution[index]->Clone ());

      RooDataHist *bkg_histogram=0;

      if(WantBKGSubtract) {
	TFile *file_ ;
	
	if (tchannel=="Gamma") file_=TFile::Open("QCD_BKG_Train.root");
	else file_=TFile::Open("TTJets13TeV.root");
	
	TString totalnbkg="1";
	TH1F * h1 = (TH1F*)file_->Get("Count");
		totalnbkg=NToString(h1->Integral());
  conditionbkg=conditionbkg+"*((xsec)*(puWeight)*("+lumi+"))/("+totalnbkg+")"; 
	TTree *treephys14bkg = (TTree *) file_->Get ("METtree");
	int bkgbin(0);
	if(variablenamepng.Contains("over"))bkgbin = 5;
	else bkgbin = 400;
	TH1F *h_ = new TH1F("h_"," ", 200, -bkgbin, bkgbin);
	cout << "condition bkg " << conditionbkg << endl;
	treephys14bkg->Draw (variablename + ">>" + TString (h_->GetName ()),	conditionbkg.Data(),"sames");
	bkg_histogram= new RooDataHist("bkg_histogram","bkg_histogram",x,h_);

      }

      // construct the voightian model
      // fit the Hist Dataset also
      // fill f and efwhm that are the parameter of the voightian
      ConstructModel(Hist, bkg_histogram, WantBKGSubtract);
                       
      //if (f/2.3 < 5) continue;
      
      RooPlot *xFrame=x.frame();
      Hist.plotOn (xFrame);

      TString titlexfit = "";
      if (variablename == "pfmetx")
	titlexfit = "MET_{x} [GeV]";
      if (variablename == "pfmety")
	titlexfit = "MET_{y} [GeV]";
      xFrame->SetXTitle (titlexfit);

      int color=kBlack;
      if (xvariable == "nVert")
	color = kRed;
      if (xvariable == "met_sumEt")
	color = kGreen+3;
      if (xvariable == "zll_pt")
	color = kBlue;
      

      cout << "plot made " << endl;
      c1->cd();
      
      //Hist.plotOn(xFrame2);
      //model->plotOn(xFrame2,RooFit::LineColor(kBlack));
      if ( WantBKGSubtract  ){
        model->plotOn(xFrame);
        model->plotOn(xFrame, Components("bkg_pdf"), LineColor(kRed), LineStyle(kDashed), FillColor(kRed), DrawOption("F"));
        model->plotOn(xFrame, Components("voigt"), LineColor(kGreen), LineStyle(kDashed), FillColor(kGreen+1), DrawOption("L"));
      }  
      else {
	Hist.plotOn(xFrame);
        voigt->plotOn(xFrame, RooFit::FillColor(kGray), VisualizeError(*result,1), RooFit::Components(*voigt)); // 1 sigma band in gray
        voigt->plotOn(xFrame, RooFit::LineColor(color));
      }                             
      TString histoname = resolution[index]->GetName ();
      xFrame->GetYaxis()->SetRangeUser(1,200);
      xFrame->GetXaxis() ->SetRangeUser(-200,200);
      xFrame->Draw();
      c1->SetLogy();
      if (!samplephys14.Contains("jes_"))      c1->Print ("~/www/"+DestFolder+"/METModel/" + folder + "/" + tchannel +"/" + histoname + "_" +	variablenamepng + "_vs_" + xvariable + ".png");
      c1->SetLogy(0);

      //c1->Print ("~/www/"+DestFolder+"/METFits/" + folder + "/" + tchannel +"/" + histoname + "_" +	variablenamepng + "_vs_" + xvariable + ".png");

      //Print chi2/dof value

      Double_t chi2 = xFrame->chiSquare ();	//"voigt", "Hist", 3);
      //cout << "chi2 = " << chi2 << endl;

      tgraphx[index] = index;

      if (xvariable == "nVert")
	tgraphx[index] = limitup;
      if (xvariable == "met_sumEt")
	tgraphx[index] = limitup * 0.001;	//For the x axis to be in TEV
      if (xvariable == "zll_pt")
	tgraphx[index] = limitup;

      tgraphxchi2[index] = tgraphx[index];
      tgraphychi2[index] = chi2;

      if (chi2 != chi2 || chi2 >= 100)
    	tgraphychi2[index] = -0.2;
      tgraphy[index] = f / 2.3548;
      if ((variablename == "met_uPara_zll_raw/zll_pt") || (variablename =="met_uPara_zll/zll_pt")|| (variablename =="met_uPara_zll_down/zll_pt")|| (variablename =="met_uPara_zll_up/zll_pt")){
	    tgraphy[index] = -resolution[index]->GetMean ();
	    

// 	    if(tgraphy[index] > 1.0){
// 	       TString condition2="(weighttotal)*(channel=="+dileptonch +")";
//                treephys14->Draw ("met_uPara_zll >> met_uPara_zllresponse1", condition.Data());
//                treephys14->Draw ("zll_pt >> zll_ptresponse1", condition.Data(), "sames");
// 	       TString NumberStr;          // string which will contain the result
//                ostringstream convert;   // stream used for the conversion
//                convert << index;      // insert the textual representation of 'Number' in the characters in the stream
//                NumberStr = convert.str(); 
// 	       met_uPara_zllresponse1->Draw("hist");
//                met_uPara_zllresponse1->GetXaxis ()->SetTitle ("u_{|| }   [GeV]");
//                c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_uPara_zll_response1_"+NumberStr+".png");
//                c1->SetLogy(0);
// 	       zll_ptresponse1->Draw("hist");
//                zll_ptresponse1->GetXaxis ()->SetTitle ("zll_pt   [GeV]");
//                c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/zll_pt_response1_"+NumberStr+".png");
//                c1->SetLogy(0);
// 	    }
	    //cout << index << "  and mean: " << -resolution[index]->GetMean () << endl;
      }
      etgraphy[index] = efwhm / 2.3548;
      if (variablename == "met_uPara_zll/zll_pt" || variablename == "met_uPara_zll_raw/zll_pt"|| variablename == "met_uPara_zll_up/zll_pt"|| variablename == "met_uPara_zll_down/zll_pt")
	    {
	    etgraphy[index] = resolution[index]->GetMeanError ();
      etgraphx[index] = 0;}


      //Set limit down
      limitdown = limitup;
      strlimitdown = Form ("%d", limitdown);


    }




  TGraph *gr =  new TGraphErrors (sizexarray, tgraphx, tgraphy, etgraphx, etgraphy);
  gr->SetMarkerColor (4);
  gr->SetMarkerStyle (21);

  TGraph *grchi2 = new TGraph (sizexarray, tgraphxchi2, tgraphychi2);
  grchi2->SetMarkerColor (2);
  grchi2->SetMarkerStyle (34);
  
  
  if (xvariable == "met_sumEt")
    {
      gr->GetXaxis ()->SetTitle ("sumE_{T} [TeV]");
      grchi2->GetXaxis ()->SetTitle ("sumE_{T} [TeV]");
    }
  if (xvariable == "nVert")
    {
      gr->GetXaxis ()->SetTitle ("Number of Vertices");
      grchi2->GetXaxis ()->SetTitle ("Number of Vertices");
    }
  if (xvariable == "zll_pt")
    {
      gr->GetXaxis ()->SetTitle ("zll_pt [GeV]");
      grchi2->GetXaxis ()->SetTitle ("zll_pt [GeV]");
    }



  gr->GetYaxis ()->SetTitle (titley);
  if (variablename == "met_uPara_zll" )
    gr->GetYaxis ()->SetTitle ("#sigma(u_{||}) [GeV]");
  if (variablename == "met_uPara_zll_raw")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{|| raw}) [GeV]");
  if (variablename == "met_uPara_zll+zll_pt")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{||} +zll_pt) [GeV]");
  if (variablename == "met_uPara_zll/zll_pt")
    gr->GetYaxis ()->SetTitle ("-<u_{||}> /zll_pt ");
  if (variablename == "met_uPara_zll_raw/zll_pt")
    gr->GetYaxis ()->SetTitle ("-<u_{|| raw}> /zll_pt ");
  if (variablename == "met_uPara_zll_raw+zll_pt")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{|| raw} +zll_pt) [GeV]");
  if (variablename == "met_uPerp_zll")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{#perp}  ) [GeV]");
  if (variablename == "met_uPerp_zll_raw")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{#perp raw}  ) [GeV]");
  

  if (drawchi2) {
  grchi2->GetYaxis ()->SetTitle ("#Chi^{2}");
  grchi2->Draw ("AP");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/" + tchannel + "/" +  variablenamepng  + "_vs_" +	     xvariable + "_chi2.png");
  c1->Clear (); }
  
  TString direction="";
  if (variablename.Contains("_up")) direction="_up_";
  if (variablename.Contains("_down")) direction="_down_";
  TFile f2 (DestFolder+folder+ direction+"tgraphs_"+samplephys14, "UPDATE");
  gr->Write (tchannel+"_"+variablenamepng + "_vs_" + xvariable);


  c1->Update ();

  TLegend *leg;
  leg = new TLegend (0.60, 0.6, 0.91, 0.81);
  leg->SetFillStyle (0);
  leg->SetBorderSize (0);
  leg->SetTextSize (0.04);
  leg->SetTextFont (42);



  leg->SetFillColor (0);





  // leg->Draw ();
  TLatex l1;
  l1.SetTextAlign (12);
  l1.SetTextSize (0.04);
  l1.SetNDC ();
  l1.DrawLatex (0.155, 0.98, "CMS Preliminary, #sqrt{s} = 13 TeV");

  
  

  gr->Draw ("AP");
  if (variablename!="met_uPara_zll/zll_pt" && variablename!="met_uPara_zll_raw/zll_pt"&& variablename!="met_uPara_zll_down/zll_pt"&& variablename!="met_uPara_zll_up/zll_pt") gr->GetYaxis()->SetRangeUser(0,70);
  else gr->GetYaxis()->SetRangeUser(0.8,1.2);
  c1->Update();
  
  
  if (variablename == "met_uPara_zll/zll_pt" || variablename=="met_uPara_zll_raw/zll_pt" || variablename=="met_uPara_zll_up/zll_pt" || variablename=="met_uPara_zll_down/zll_pt")
    {
      TLine *lineR =  new TLine ( gr->GetHistogram ()->GetXaxis ()->GetXmin (), 1, gr->GetHistogram ()->GetXaxis ()->GetXmax (), 1);
      lineR->SetLineColor (kBlue + 1);
      lineR->SetLineWidth (2);
      lineR->SetLineStyle (2);
      lineR->Draw ();
    }

  else
    {
      
      TLine *lineR =	new TLine (gr->GetHistogram ()->GetXaxis ()->GetXmin (), 0,  gr->GetHistogram ()->GetXaxis ()->GetXmax (), 0);
      lineR->SetLineColor (kBlue + 1);
      lineR->SetLineWidth (2);
      lineR->SetLineStyle (2);
      lineR->Draw ();


    }

  

   if (!samplephys14.Contains("jes_"))  c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/" + tchannel + "/" + variablenamepng  + "_vs_" +	     xvariable + ".png");

  */
  }
}



//--------------------------------------------------------------------------------------------------
void performFit(vector<TH1*> hv, vector<TH1*> hbkgv, 
                const Int_t model, const Bool_t sigOnly,
                const char *plabel, const char *xlabel,
		Double_t *meanArr,   Double_t *meanErrArr,
		Double_t *sigma0Arr, Double_t *sigma0ErrArr,
		Double_t *sigma1Arr, Double_t *sigma1ErrArr,
		Double_t *sigma2Arr, Double_t *sigma2ErrArr,
		Double_t *sigma3Arr, Double_t *sigma3ErrArr,
		Double_t *frac2Arr,  Double_t *frac2ErrArr,
		Double_t *frac3Arr,  Double_t *frac3ErrArr, 
        TString outputDir
) {
  char pname[50];
  char ylabel[50];
  char binlabel[50];
  char nsigtext[50];
  char nbkgtext[50];
  
  char meantext[50];
  char sig0text[50];
  char sig1text[50];
  char sig2text[50];
  char sig3text[50];
  
  
  for(Int_t ibin=0; ibin<nbins; ibin++) {
    if (hv.size()!=nbins || hbkgv.size()!=nbins)gSystem->Abort();

    RooRealVar u("u","u",hv[ibin]->GetXaxis()->GetXmin(),hv[ibin]->GetXaxis()->GetXmax());
    u.setBins(100);
    RooDataHist dataHist("dataHist","dataHist",RooArgSet(u),hv[ibin]);
    
    //
    // Set up background histogram templates
    //
    RooDataHist bkgHist("bkgHist","bkgHist",RooArgSet(u),hbkgv[ibin]);
    RooHistPdf bkg("bkg","bkg",u,bkgHist,0);
    
    //
    // Set up fit parameters
    //
    RooRealVar mean("mean","mean",
                    hv[ibin]->GetMean(),
                    //hv[ibin]->GetXaxis()->GetXmin(),
                    //hv[ibin]->GetXaxis()->GetXmax());
		    hv[ibin]->GetMean()-1.,
		    hv[ibin]->GetMean()+1.);
    RooRealVar mean1("mean1","mean1",
                    hv[ibin]->GetMean(),
                    //hv[ibin]->GetXaxis()->GetXmin(),
                    //hv[ibin]->GetXaxis()->GetXmax());
		    hv[ibin]->GetMean()-1.,
		    hv[ibin]->GetMean()+1.);
    RooRealVar mean2("mean2","mean2",
                    hv[ibin]->GetMean()-2,
                    //hv[ibin]->GetXaxis()->GetXmin(),
                    //hv[ibin]->GetXaxis()->GetXmax());
		    hv[ibin]->GetMean()-3.,
		     hv[ibin]->GetMean());
    RooRealVar sigma1("sigma1","sigma1",0.5*(hv[ibin]->GetRMS()),0,1.5*(hv[ibin]->GetRMS()));
    RooRealVar sigma2("sigma2","sigma2",0.5*(hv[ibin]->GetRMS()),0,2.0*(hv[ibin]->GetRMS()));
    RooRealVar sigma3("sigma3","sigma3",10.*(hv[ibin]->GetRMS()),0,200); 
    RooRealVar frac2("frac2","frac2",0.5,0.1,0.9);
    RooRealVar frac3("frac3","frac3",0.05,0,0.15);
    RooGaussian gauss1("gauss1","gauss1",u,mean1,sigma1);
    RooGaussian gauss2("gauss2","gauss2",u,mean2,sigma2);
    RooGaussian gauss3("gauss3","gauss3",u,mean,sigma3);
    RooVoigtian voight("voigt", "Voigtian", u, mean, gamma_Z0, g_w);

    
    if(ibin>0) {
      
      // set the limits of the widths and means
      mean.setVal(meanArr[ibin-1]);
      sigma1.setMin(0);
      sigma1.setMax(1.5*(hv[ibin-1]->GetRMS()));
      sigma1.setVal(0.3*(hv[ibin-1]->GetRMS()));
      sigma2.setMin(0);
      sigma2.setMax(1.8*(hv[ibin-1]->GetRMS()));
      sigma2.setVal(1.0*(hv[ibin-1]->GetRMS()));
      
      //             mean.setVal(meanArr[ibin-1]);
      //       sigma1.setMin(0.1*(sigma0Arr[ibin-1]));
      //       sigma1.setMax(1.5*(sigma0Arr[ibin-1]));
      //       sigma1.setVal(0.65*(sigma0Arr[ibin-1]));
      //       sigma2.setMin(0.1*(sigma0Arr[ibin-1]));
      //       sigma2.setMax(1.8*(sigma0Arr[ibin-1]));
      //       sigma2.setVal(1.2*(sigma0Arr[ibin-1]));
      
    }
    
    //
    // Define formula for overall width (sigma0)
    //
    char formula[100];
    RooArgList params;
    if(model==1) {
      sprintf(formula,"sigma1");
      params.add(sigma1);
    
    } else if(model==2) {
      sprintf(formula,"(1.-frac2)*sigma1 + frac2*sigma2");
      params.add(frac2);
      params.add(sigma1);
      params.add(sigma2);

    } else if(model==3) {
      sprintf(formula,"(1.-frac2-frac3)*sigma1 + frac2*sigma2 + frac3*sigma3");
      params.add(frac2);
      params.add(frac3);
      params.add(sigma1);
      params.add(sigma2);
      params.add(sigma3);
    }       
    RooFormulaVar sigma0("sigma0",formula,params);
    
    //
    // Construct fit model
    //
    
    RooArgList shapes;
    if(model>=3) shapes.add(gauss3);
    if(model>=2) shapes.add(gauss2);
    shapes.add(gauss1);
    //shapes.add(voight);
    
  
    RooArgList fracs;
    if(model>=3) fracs.add(frac3);
    if(model>=2) fracs.add(frac2);
    
    RooAddPdf sig("sig","sig",shapes,fracs);
    
    RooArgList parts;
    parts.add(sig);
    if(!sigOnly) parts.add(bkg);
    
    RooArgList yields;
    RooRealVar nsig("nsig","nsig",0.98*(hv[ibin]->Integral()),0.,hv[ibin]->Integral());
    yields.add(nsig);
    RooRealVar nbkg("nbkg","nbkg",0.01*(hv[ibin]->Integral()),0.,0.50*(hv[ibin]->Integral()));
    if(!sigOnly) yields.add(nbkg);
    else         nbkg.setVal(0);
    
    RooAddPdf modelpdf("modelpdf","modelpdf",parts,yields);
        
    //
    // Perform fit
    //
    RooFitResult *fitResult=0;
    fitResult = modelpdf.fitTo(dataHist,
                               //RooFit::Minos(),
			       RooFit::Strategy(2),
	                       RooFit::Save());
    
    if(sigma1.getVal() > sigma2.getVal()) {
      Double_t wide = sigma1.getVal();
      Double_t thin = sigma2.getVal();
      sigma1.setVal(thin);
      sigma2.setVal(wide);
      frac2.setVal(1.0-frac2.getVal());
      fitResult = modelpdf.fitTo(dataHist,
                                 //RooFit::Minos(),
			         RooFit::Strategy(2),
	                         RooFit::Save());
    }
    
    meanArr[ibin]      = mean.getVal();
    meanErrArr[ibin]   = mean.getError();
    sigma1Arr[ibin]    = sigma1.getVal();
    sigma1ErrArr[ibin] = sigma1.getError();
    if(model>=2) {
      sigma0Arr[ibin]    = sigma0.getVal();
      sigma0ErrArr[ibin] = sigma0.getPropagatedError(*fitResult);
      sigma2Arr[ibin]    = sigma2.getVal();
      sigma2ErrArr[ibin] = sigma2.getError();
      frac2Arr[ibin]     = frac2.getVal();
      frac2ErrArr[ibin]  = frac2.getError();
    }
    if(model>=3) {    
      sigma3Arr[ibin]    = sigma3.getVal();
      sigma3ErrArr[ibin] = sigma3.getError();
      frac3Arr[ibin]     = frac3.getVal();
      frac3ErrArr[ibin]  = frac3.getError();
    }
    
    //
    // Plot fit results
    //
    RooPlot *frame = u.frame(Bins(100));
    dataHist.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelpdf.plotOn(frame);
    if(!sigOnly) modelpdf.plotOn(frame,Components("bkg"),LineStyle(kDotted),LineColor(kMagenta+2));
    if(model>=2) sig.plotOn(frame,Components("gauss1"),LineStyle(kDashed),LineColor(kRed));
    if(model>=2) sig.plotOn(frame,Components("gauss2"),LineStyle(kDashed),LineColor(kCyan+2));
    if(model>=3) sig.plotOn(frame,Components("gauss3"),LineStyle(kDashed),LineColor(kGreen+2));
    
    sprintf(pname,"%sfit_%i",plabel,ibin);
    sprintf(ylabel,"Events / %.1f GeV",hv[ibin]->GetBinWidth(1));
    sprintf(binlabel,"%i < p_{T} < %i",(Int_t)ptbins[ibin],(Int_t)ptbins[ibin+1]);    
    if(sigOnly) {
      sprintf(nsigtext,"N_{evts} = %i",(Int_t)hv[ibin]->Integral());
    } else {
      sprintf(nsigtext,"N_{sig} = %.1f #pm %.1f",nsig.getVal(),nsig.getError());
      sprintf(nbkgtext,"N_{bkg} = %.1f #pm %.1f",nbkg.getVal(),nbkg.getError());
    }
    sprintf(meantext,"#mu = %.1f #pm %.1f",meanArr[ibin],meanErrArr[ibin]);
    sprintf(sig1text,"#sigma = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);
    if(model>=2) {
      sprintf(sig0text,"#sigma = %.1f #pm %.1f",sigma0Arr[ibin],sigma0ErrArr[ibin]);
      sprintf(sig1text,"#sigma_{1} = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);          
      sprintf(sig2text,"#sigma_{2} = %.1f #pm %.1f",sigma2Arr[ibin],sigma2ErrArr[ibin]);
    }
    if(model>3)
      sprintf(sig3text,"#sigma_{3} = %.1f #pm %.1f",sigma3Arr[ibin],sigma3ErrArr[ibin]);
    


    TPaveText *tb = new TPaveText(0.21,0.60,0.51,0.85,"NDC");
    tb->SetTextColor(kBlack);
    tb->SetFillStyle(0);
    tb->SetBorderSize(0);
    tb->AddText(binlabel);
    if(sigOnly) tb->AddText(nsigtext);
    if(model==1) {tb->AddText(meantext); tb->AddText(sig1text);}
    if(model==2) {tb->AddText(meantext); tb->AddText(sig0text); tb->AddText(sig1text); tb->AddText(sig2text);}
    if(model==3) {tb->AddText(meantext); tb->AddText(sig0text); tb->AddText(sig1text); tb->AddText(sig2text); tb->AddText(sig3text);}

    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(xlabel);
    frame->GetYaxis()->SetTitle(ylabel);
    frame->Draw();
    tb->Draw("same");
    //c->SaveAs(TString(outputDir+"/plots/"+pname+".png"));
    
    sprintf(pname,"%sfitlog_%i",plabel,ibin);
    //c->SetLogy();
    //c->SaveAs(TString(outputDir+"/plots/"+pname+".png"));
    
    //c->Clear();
    //c->SetLogy(0);
    delete tb;

  }
}

void 
identifyToken (std::set<char> & token, const std::string& tokenize ) 
{
    //Function to automatically identify token used by the user in the string literal 'tokenize'.
    //Token could be everything other a number, i.e. everything with ASCII code outside the [48,57] range
    //other than a letter, i.e. everything with ASCII code outside the [65,90] and [97,122] range
    //and other than an underscore, i.e. with ASCII code different than 95
    for(unsigned int count =0; count < tokenize.size(); count++)
      if (( (int) tokenize[count]>=48 && (int) tokenize[count]<=57 ) ||
	  ( (int) tokenize[count]>=65 && (int) tokenize[count]<=90 ) || 
	  ( (int) tokenize[count]>=97 && (int) tokenize[count]<=122 ) ||
	  (int) tokenize[count]==95 ) 
      
	continue;
      else
	if ( (int) tokenize[count]!=32 )
	  token.insert(tokenize[count]);

    if (token.size()>1) throw "Please define one delimeter at most!";
}

void 
removeWhiteSpaces( std::string &strg ) 
{  
    //Function to remove white spaces that accidentally have been added by the user in the string literal
    int i;
    i = strg.find("  ");
  
    if ( i > -1 )
    {
        removeWhiteSpaces (  strg.replace( i,2,"" ) ); //recursive calling of the function itself
    }
    else if ( i == -1 )
    {
        i = strg.find(" ");
	if ( i > -1 )
	{
	    strg.erase(0,1);
	}
    }
}

/*
double
FWHM (double sigma, double gamma)
{

  double f_g = 2 * sigma * sqrt (2 * log (2));
  double f_l = 2 * gamma;

  return 0.5346 * 2 * gamma + sqrt (0.2166 * f_l * f_l + f_g * f_g);
}



double
FWHMError (double sigma, double gamma, double esigma, double egamma,
	   double Vss, double Vsg, double Vgs, double Vgg)
{


  double a = 0.5346;
  double b = 0.2166;
  double ef_g = 2 * esigma * sqrt (2 * log (2));
  double ef_l = 2 * egamma;

  double dg =
    2 * a + 4 * b * gamma / sqrt (4 * b * pow (gamma, 2) +
				  4 * pow (sigma, 2) * log (2));

  double ds =
    (sigma * log (4)) / sqrt (b * pow (gamma, 2) + pow (sigma, 2) * log (2));
  
  double p1 = ef_l * ef_l * Vgg * dg;
  double p2 = ef_g * ef_l * Vsg * dg * ds;	//identical (should be)
  double p3 = ef_g * ef_l * Vgs * dg * ds;
  double p4 = ef_g * ef_g * Vss * ds;

  return sqrt (abs (p1) + abs (p2) + abs (p3) + abs (p4));

}


double
FWHMError_fixed (double sigma, double gamma, double esigma, double egamma,
	   double Vss, double Vsg, double Vgs, double Vgg)
{
  // Vss = correlation(sigma, sigma)
  // Vsg = correlation(sigma, gamma)
  // etc
  double a = 0.5346;
  double b = 0.2166;
  double c = 2 * sqrt( 2*log(2) );
  double f_g = c * sigma;
  double f_l = 2 * gamma;
  double sq = sqrt( b * pow(f_l, 2) + pow(f_g, 2) );
  
  // Partial derivatives of f_voigtian w.r.t sigma and gamma
  // f = a * f_l + sqrt( b * f_l^2 + f_g^2 )
  double dfds = c * ( f_g / sq ) ;
  double dfdg = 2 * ( a + b * f_l / sq ) ;
  
  // esigma * esigma * pow( Vss, 2 ) gives covariance(sigma, sigma) etc
  double p1 = dfds * dfds * esigma * esigma * pow( Vss, 2 );
  double p2 = dfds * dfdg * esigma * egamma * pow( Vsg, 2 );
  double p3 = dfdg * dfds * egamma * esigma * pow( Vgs, 2 );
  double p4 = dfdg * dfdg * egamma * egamma * pow( Vgg, 2 );

  return sqrt ( p1 + p2 + p3 + p4 );
}
TString
NToString(Float_t type) 
{ 
   TString name; name.Form("%f",type); 
      cout << name << endl;   
         return name;
         
            } 
            
*/            
