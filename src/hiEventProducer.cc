//
// producer of container with heavy-ion event properties
//

#include "UserCode/diall/interface/hiEventProducer.h"

#include <ostream>
using namespace std;

ClassImp(hiEventProducer)

//__________________________________________________________
hiEventProducer::hiEventProducer() :
inputBase("hiEventProducer"),
  fhiEventContName("hiEventContainer"),
  fhiEventContainer(0),
  fRun(-1),
  fEvt(-1),
  fLumi(-1),
  fWeight(-1),
  fVx(0.),
  fVy(0.),
  fVz(0.),
  fHiBin(0),
  fHiHF(0),
  fNPV(0),
  fHLT_AK4PFJet80_Eta5p1_v1(0),
  fHLT_HISinglePhoton30_Eta3p5_v1(0),
  fHLT_FullTrack24ForPPRef_v1(0),
  fHLT_FullTrack45ForPPRef_v1(0),
  fcollisionEvtSel(0),
  fHBHENoise(0)
{
  //default constructor
}

//__________________________________________________________
hiEventProducer::hiEventProducer(const char *name) :
  inputBase(name),
  fhiEventContName("hiEventContainer"),
  fhiEventContainer(0),
  fRun(-1),
  fEvt(-1),
  fLumi(-1),
  fWeight(-1),
  fVx(0.),
  fVy(0.),
  fVz(0.),
  fHiBin(0),
  fHiHF(0),
  fNPV(0),
  fHLT_AK4PFJet80_Eta5p1_v1(0),
  fHLT_HISinglePhoton30_Eta3p5_v1(0),
  fHLT_FullTrack24ForPPRef_v1(0),
  fHLT_FullTrack45ForPPRef_v1(0),
  fcollisionEvtSel(0),
  fHBHENoise(0)
{
  //standard constructor
}

//__________________________________________________________
void hiEventProducer::SetInput(TChain *chain) {

  inputBase::SetInput(chain);
  //Init();
  
}

//__________________________________________________________
Bool_t hiEventProducer::Init() {

  if(!inputBase::Init()) return kFALSE;

  if(fInputMode==hiForest) {
    fChain->SetBranchAddress("run",   &fRun,   &b_run);
    fChain->SetBranchAddress("evt",   &fEvt,   &b_evt);
    fChain->SetBranchAddress("lumi",  &fLumi,  &b_lumi);
    if(fChain->GetBranch("weight"))
      fChain->SetBranchAddress("weight",&fWeight,&b_weight);
    fChain->SetBranchAddress("vx",    &fVx,    &b_vx);
    fChain->SetBranchAddress("vy",    &fVy,    &b_vy);
    fChain->SetBranchAddress("vz",    &fVz,    &b_vz);
    fChain->SetBranchAddress("hiBin", &fHiBin, &b_hiBin);
    fChain->SetBranchAddress("hiHF",  &fHiHF,  &b_hiHF);
    fChain->Print();
    if(fChain->GetBranch("nPU"))
      fChain->SetBranchAddress("nPU",&fNPV,&b_nPU);
    if(fChain->GetBranch("HLT_AK4PFJet80_Eta5p1_v1"))
      fChain->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1",&fHLT_AK4PFJet80_Eta5p1_v1,&b_HLT_AK4PFJet80_Eta5p1_v1);
    if(fChain->GetBranch("HLT_HISinglePhoton30_Eta3p5_v1"))
      fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta3p5_v1",&fHLT_HISinglePhoton30_Eta3p5_v1,&b_HLT_HISinglePhoton30_Eta3p5_v1);
    if(fChain->GetBranch("HLT_FullTrack24ForPPRef_v1"))
      fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v1",&fHLT_FullTrack24ForPPRef_v1,&b_HLT_FullTrack24ForPPRef_v1);
    if(fChain->GetBranch("HLT_FullTrack45ForPPRef_v1"))
      fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v1",&fHLT_FullTrack45ForPPRef_v1,&b_HLT_FullTrack45ForPPRef_v1);
    if(fChain->GetBranch("PAcollisionEventSelection"))
      fChain->SetBranchAddress("PAcollisionEventSelection",&fcollisionEvtSel,&b_collisionEvtSel);
    if(fChain->GetBranch("pHBHENoiseFilterResultProducer"))
      fChain->SetBranchAddress("pHBHENoiseFilterResultProducer",&fHBHENoise,&b_HBHENoise);

    fInit = kTRUE;
  }
  return kTRUE;
}
//__________________________________________________________
Bool_t hiEventProducer::InitEventObjects() {
  
  //Create event objects
  if(!fEventObjects) {
    Printf("%s: fEventObjects does not exist. Cannot store output",GetName());
    return kFALSE;
  } else {
    if(!fEventObjects->FindObject(fhiEventContName)) {
      fhiEventContainer = new hiEventContainer(fhiEventContName);
      fEventObjects->Add(fhiEventContainer);
    }
  }
  return kTRUE;
}

//__________________________________________________________
Bool_t hiEventProducer::Run(Long64_t entry) {

  //run analysis
  Long64_t centry = LoadTree(entry);
  if(centry<0) return kFALSE;

  if(!InitEventObjects()) return kFALSE;
    
  fhiEventContainer->SetRun(fRun);
  fhiEventContainer->SetEvent(fEvt);
  fhiEventContainer->SetLumi(fLumi);
  fhiEventContainer->SetWeight(fWeight);
  fhiEventContainer->SetVx(fVx);
  fhiEventContainer->SetVy(fVy);
  fhiEventContainer->SetVz(fVz);
  fhiEventContainer->SetHiBin(fHiBin);
  fhiEventContainer->SetHiHF(fHiHF);
  fhiEventContainer->SetNPV(fNPV);
  fhiEventContainer->SetPFJet80(fHLT_AK4PFJet80_Eta5p1_v1);
  fhiEventContainer->SetPhoton30(fHLT_HISinglePhoton30_Eta3p5_v1);
  fhiEventContainer->SetTrk24(fHLT_FullTrack24ForPPRef_v1);
  fhiEventContainer->SetTrk45(fHLT_FullTrack45ForPPRef_v1);
  fhiEventContainer->SetColl(fcollisionEvtSel);
  fhiEventContainer->SetHBHENoise(fHBHENoise);

  return kTRUE; 
}

//__________________________________________________________
Long64_t hiEventProducer::LoadTree(Long64_t entry) {

  //overloaded LoadTree function
  if(!fChain) {
    Printf("fChain doesn't exist");
    return -5;
  }

  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Init();
    //Printf("%lld fCurrent: %d",entry,fCurrent);
  }
 
  // fChain->SetMakeClass(1);
 
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) {
    Printf("hiEventProducer: centry smaller than 0");
    return centry;  
  }
  
  fChain->GetEntry(entry);
  
  return centry;  
}
