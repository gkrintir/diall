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
   virtual ~anaTtbarDilepton();
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
   void SetMinDilepMass(Float_t m)     { fMinMass = m; }
   void SetMinMET(Float_t m)           { fMinMET = m; }
   void SetTriggerMapName(TString name){ fTriggerMapName = name ; }
   void SetTriggerPath(TString name){ fTriggerPath = name ; }
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
   Float_t           fMinPt;                //minimum pT of PFparticles
   Float_t           fMinMass;              //minimum dilepton mass
   Float_t           fMinMET;                //minimum MET
   TString           fTriggerMapName;       //name of triggers
   TString           fTriggerPath;          //name of triggers
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
   std::vector<lwJet*> fRecoRemovalJets;    //!jet container
   TClonesArray     *fParticles;            //!pfParticle array
   TClonesArray     *fDilepton;             //!dilepton candidates 
   TClonesArray     *fTops;                 //!Top quark candidates container

   Int_t            fDileptonCands;
   TH1F             *fNEvents;            //!# selected events
   TTree            *fAnaTtbarDileptonInfo;

   std::vector<Int_t>     fTriggerFired;        //!# addtional trigger of interest 
   std::vector<Float_t>   fEventWeight;         //!# generator weight
   std::vector<Float_t>   fProcessXSection;     //!# generator x-section
   std::vector<Int_t>     fNRecoLeptonLead;     //!# selected electrons in event
   std::vector<Int_t>     fNRecoLeptonSublead;  //!# selected muons in event
   std::vector<Int_t>     fNJetsIncl;           //!# selected jets in event
   std::vector<Int_t>     fNBJetsIncl;          //!# selected b-tagged jets in event
   std::vector<Float_t>   fBJetsDiscr;          //!# selected b-tagged jets in event
   std::vector<Float_t>   fNDileptons;          //!# dilepton mass varialbe
   std::vector<Float_t>   fMassDilepton_new;    //!# dilepton mass varialbe
   std::vector<Float_t>   fSignDilepton_new;    //!# dilepton mass varialbe
   std::vector<Float_t>   fEtaDilepton_new;     //!# dilepton mass varialbe
   std::vector<Float_t>   fPhiDilepton_new;     //!# dilepton mass varialbe
   std::vector<Float_t>   fPtDilepton_new;      //!# dilepton mass varialbe
   std::vector<Float_t>   fDeltaPhi_new;
   std::vector<Float_t>   fLeadJetPt_new;           //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetEta_new;
   std::vector<Float_t>   fLeadRecoLeptonPt_new;    //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonAbsEta_new;//!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonEta_new;   //!# dilepton mass varialbe
   //std::vector<Float_t>   fMETAbs;              //!# dilepton mass varialbe
   //std::vector<Float_t>   fMETPhi;              //!# dilepton mas
   
   //Dilepton-related variables
   std::vector<Float_t>   fHT;                  //!# HT varialbe
   std::vector<Float_t>   fMassDilepton;        //!# dilepton mass varialbe
   std::vector<Float_t>   fEtaDilepton;         //!# dilepton mass varialbe
   std::vector<Float_t>   fPhiDilepton;         //!# dilepton mass varialbe
   std::vector<Float_t>   fPtDilepton;          //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetPt;           //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonPt;    //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonAbsEta;//!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonEta;   //!# dilepton mass varialbe
   std::vector<Float_t>   fMETAbs;              //!# dilepton mass varialbe
   std::vector<Float_t>   fMETPhi;              //!# dilepton mas
   std::vector<Float_t>   fMETAbsNoHF;              //!# dilepton mass varialbe
   std::vector<Float_t>   fMETPhiNoHF;              //!# dilepton mas   
   std::vector<Int_t>     fNJets;               //!# dilepton mass varialbe

   //Muon-related variables
   std::vector<Float_t>   fMuonPtLead;         // !# muon varialbe
   std::vector<Float_t>   fMuonEtaLead;        // !# muon varialbe
   std::vector<Float_t>   fMuonIsoLead;
   //Muon-related variables
   std::vector<Float_t>   fMuonPtSublead;      // !# muon varialbe
   std::vector<Float_t>   fMuonEtaSublead;     // !# muon varialbe
   std::vector<Float_t>   fMuonIsoSublead;
   
   //Electron-related variables
   std::vector<Float_t>   fElePtLead;         // !# electron varialbe
   std::vector<Float_t>   fEleEtaLead;        // !# electron varialbe
   std::vector<Float_t>   fEleIsoLead;        // !# electron varialbe
   std::vector<Float_t>   fEleIso03Lead;      // !# electron varialbe
   std::vector<Float_t>   fEleIso04Lead;      // !# electron varialbe
   std::vector<Float_t>   fdEtaAtVtxLead;     // !# electron varialbe
   std::vector<Float_t>   fdPhiAtVtxLead;     // !# electron varialbe
   std::vector<Float_t>   fSigmaIEtaIEtaLead; // !# electron varialbe
   std::vector<Float_t>   fHoverELead;        // !# electron varialbe
   std::vector<Float_t>   fDxyLead;           // !# electron varialbe
   std::vector<Float_t>   fDzLead;            // !# electron varialbe
   std::vector<Float_t>   fEoverPInvLead;     // !# electron varialbe
   std::vector<Int_t>     fMissHitsLead;      // !# electron varialbe
   std::vector<Int_t>     fConversionVetoLead;// !# electron varialbe
   
   //Electron-related variables
   std::vector<Float_t>   fElePtSublead;        // !# electron varialbe
   std::vector<Float_t>   fEleEtaSublead;        // !# electron varialbe
   std::vector<Float_t>   fEleIsoSublead;        // !# electron varialbe
   std::vector<Float_t>   fEleIso03Sublead;      // !# electron varialbe
   std::vector<Float_t>   fEleIso04Sublead;      // !# electron varialbe
   std::vector<Float_t>   fdEtaAtVtxSublead;     // !# electron varialbe
   std::vector<Float_t>   fdPhiAtVtxSublead;     // !# electron varialbe
   std::vector<Float_t>   fSigmaIEtaIEtaSublead; // !# electron varialbe
   std::vector<Float_t>   fHoverESublead;        // !# electron varialbe
   std::vector<Float_t>   fDxySublead;           // !# electron varialbe
   std::vector<Float_t>   fDzSublead;            // !# electron varialbe
   std::vector<Float_t>   fEoverPInvSublead;     // !# electron varialbe
   std::vector<Int_t>     fMissHitsSublead;      // !# electron varialbe
   std::vector<Int_t>     fConversionVetoSublead;// !# electron varialbe

   
   std::vector<Double_t>  fDeltaPhi_Rec2jets;       //!# HT varialbe
   std::vector<Float_t>   fMassDilepton_Rec2jets;   //!# dilepton mass varialbe
   std::vector<Float_t>   fEtaDilepton_Rec2jets;   //!# dilepton mass varialbe
   std::vector<Float_t>   fPtDilepton_Rec2jets;   //!# dilepton mass varialbe
   std::vector<Float_t>   fPhiDilepton_Rec2jets;   //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetPt_Rec2jets;      //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetAbsEta_Rec2jets;  //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetEta_Rec2jets;     //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonPt_Rec2jets;    //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonAbsEta_Rec2jets;//!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonEta_Rec2jets;   //!# dilepton mass varialbe

   std::vector<Double_t>  fDeltaPhi_Rec2jets_METCut; //!# HT varialbe
   std::vector<Float_t>   fMassDilepton_Rec2jets_METCut;   //!# dilepton mass varialbe
   std::vector<Float_t>   fEtaDilepton_Rec2jets_METCut;   //!# dilepton mass varialbe
   std::vector<Float_t>   fPtDilepton_Rec2jets_METCut;   //!# dilepton mass varialbe
   std::vector<Float_t>   fPhiDilepton_Rec2jets_METCut;   //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetPt_Rec2jets_METCut;      //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetAbsEta_Rec2jets_METCut;  //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadJetEta_Rec2jets_METCut;     //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonPt_Rec2jets_METCut;    //!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonAbsEta_Rec2jets_METCut;//!# dilepton mass varialbe
   std::vector<Float_t>   fLeadRecoLeptonEta_Rec2jets_METCut;   //!# dilepton mass varialbe
   
   TH2F                 *fh2MetCentPtMin[10];   //!MET vs centrality for various min pt cuts
   
   //ClassDef(anaTtbarDilepton,2)
};
#endif
