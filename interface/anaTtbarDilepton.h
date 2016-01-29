#ifndef anaTtbarDilepton_h
#define anaTtbarDilepton_h

#include "TString.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TClonesArray.h"

#include "UserCode/diall/interface/anaBaseTask.h"

#include "UserCode/diall/interface/triggerMap.h"
#include "UserCode/diall/interface/lwJetContainer.h"
#include "UserCode/diall/interface/particleBase.h"
#include "UserCode/diall/interface/diParticle.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooVoigtian.h>
#include <RooHistPdf.h>
#include <RooFormulaVar.h>

#include <utility>
#include <vector>

class anaTtbarDilepton : public anaBaseTask {
   
public:
  enum metType {
    kGen   = 0,
    kGenEm = 1,
    kPFRaw = 2,
    kVS    = 3,
    kPuppi = 4,
    kCS    = 5
  };
  
   anaTtbarDilepton() {;}
   anaTtbarDilepton(const char *name, const char *title);
   virtual ~anaTtbarDilepton() {;}
   void Exec(Option_t *option="");
   void ConstructModel(RooDataHist Hist, RooDataHist *bkg_hist, bool BKGSubtract);
   std::pair<double, double> compMETProjU(diParticle* zP4, double metPx, double metPy, int& errorFlag, int count);
   std::pair<double, double> compHadronicRecoilProjU(diParticle* zP4, TLorentzVector MET, int& errorFlag, int count);
   void CreateOutputObjects();
   void FillDileptonArray(const Int_t nLeptonLead, const Int_t nLeptonSublead);
   double FWHM (double sigma, double gamma);
   double FWHMError (double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg);
   double FWHMError_fixed (double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg);
   
   TString NToString(Float_t type);
   void SetCheckPid(Bool_t b)          { fCheckPid = b; }
   void SetDileptonChargeSign(Int_t b) { fDileptonChargeSign = b; }
   void SetMetType(metType t)          { fMetType = t; }
   void SetMinPt(Float_t m)            { fMinPt = m; }
   void SetTriggerMapName(TString name){ fTriggerMapName = name ; }
   void SetRecoLeptonLeadName(TString name) { fRecoLeptonLeadName = name ; }
   void SetRecoLeptonSubleadName(TString name)  { fRecoLeptonSubleadName = name ; }
   void SetGenLeptonName(TString name)  { fGenLeptonName = name ; }
   void SetGenJetsName(TString name)   { fGenJetsName = name ; }     
   void SetRecoJetsName(TString name)  { fRecoJetsName = name ; } 
   void SetParticlesName(TString name) { fParticlesName = name ; }
   
 protected:
   Bool_t            CheckPid(particleBase *p);

   Bool_t            fCheckPid;             //check if candidates are really muons (for simulation)
   metType           fMetType;              //matching type (defines where to store)
   Float_t           fMinPt;                //minimum pT of particles
   TString           fTriggerMapName;       //name of triggers
   TString           fRecoLeptonLeadName;   //name of reco leptons
   TString           fRecoLeptonSubleadName;//name of reco leptons
   TString           fGenLeptonName;        //name of gen leptons
   TString           fGenJetsName;          //name of generated jet container
   TString           fRecoJetsName;         //name of reconstrcuted jet container
   TString           fParticlesName;        //name of pf (particle flow) Particles
   Int_t             fDileptonChargeSign;   //charge sign(+1 or -1) of dilepton candidates
   TString           fDileptonName;         //name of dilepton candidates
   TString           fTopsName;             //name of top quark candidates

   triggerMap       *fTriggerMap;           //!trigger map
   TClonesArray     *fRecoLeptonLead;       //!reco lepton array
   TClonesArray     *fRecoLeptonSublead;    //!reco lepton array
   TClonesArray     *fGenLepton;            //!gen lepton arra
   lwJetContainer   *fGenJets;              //!jet container
   lwJetContainer   *fRecoJets;             //!jet container
   TClonesArray     *fParticles;            //!pfParticle array
   TClonesArray     *fDilepton;             //!dilepton candidates 
   TClonesArray     *fTops;                 //!Top quark candidates container

   Int_t            fDileptonCands;
   TH1F             *fNEvents;            //!# selected events
   TTree            *fAnaTtbarDileptonInfo;
   
   std::vector<Float_t>   fEventWeight;         //!# selected b-tagged jets in event
   std::vector<Int_t>     fNRecoLeptonLead;        //!# selected electrons in event
   std::vector<Int_t>     fNRecoLeptonSublead;     //!# selected muons in event
   std::vector<Int_t>     fNJetsIncl;          //!# selected jets in event
   std::vector<Int_t>     fNBJetsIncl;         //!# selected b-tagged jets in event
   std::vector<Float_t>   fBJetsDiscr;         //!# selected b-tagged jets in event
   
   std::vector<Float_t>   fHT;                 //!# HT varialbe
   std::vector<Float_t>   fMassDilepton;            //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetPt;          //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonPt;       //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonAbsEta;   //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonEta;      //!# dilepton mass varialbe
   std::vector<Int_t>     fNJets;              //!# dilepton mass varialbe

   
   std::vector<Double_t>  fDeltaPhi_Rec2jets;  //!# HT varialbe
   std::vector<Float_t>   fMassDilepton_Rec2jets;   //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetPt_Rec2jets; //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetAbsEta_Rec2jets; //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetEta_Rec2jets;    //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonPt_Rec2jets;  //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonAbsEta_Rec2jets;//!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonEta_Rec2jets;   //!# dilepton mass varialbe

   TH2F                 *fh2MetCentPtMin[10];   //!MET vs centrality for various min pt cuts
   
   //ClassDef(anaTtbarDilepton,2)
};
#endif
