//
// producer of muon candidates
//

#include "UserCode/diall/interface/lwMuonProducer.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/lwMuon.h"


ClassImp(lwMuonProducer)

//__________________________________________________________
lwMuonProducer::lwMuonProducer() :
inputBase("lwMuonProducer"),
  flwMuonsRecoName("lwMuonsReco"),
  flwMuonsReco(0x0),
  flwMuonsGeneName("lwMuonsGene"),
  flwMuonsGene(0x0),
  fMuons(),
  fPtMin(10.),
  fMaxEtaAbs(2.1),
  fMaxTrkChi2(4.),
  fMaxGlbChi2(10.),
  fMinNMuHits(0),
  fMinMS(1),
  fMaxDxy(0.2),//3.),
  fMaxDz(0.5),//15.),
  fMaxTrkDxy(0.3),//3.),
  fMaxTrkDz(20.),
  fMinNPixHits(0),
  fMinTrkLWM(5)
{
  //default constructor
}

//__________________________________________________________
lwMuonProducer::lwMuonProducer(const char *name) :
  inputBase(name),
  flwMuonsRecoName("lwMuonsReco"),
  flwMuonsReco(0x0),
  flwMuonsGeneName("lwMuonsGene"),
  flwMuonsGene(0x0),
  fMuons(),
  fPtMin(18.),
  fMaxEtaAbs(2.4),
  fMaxTrkChi2(4.),
  fMaxGlbChi2(10.),
  fMinNMuHits(0),
  fMinMS(1),
  fMaxDxy(0.2),//3.),
  fMaxDz(0.5),//15.),
  fMaxTrkDxy(0.3),//3.),
  fMaxTrkDz(20.),
  fMinNPixHits(0),
  fMinTrkLWM(5)
{
  //standard constructor
}

//__________________________________________________________
Bool_t lwMuonProducer::Init() {

  if(!inputBase::Init()) return kFALSE;
  
  if(fInputMode==hiForest) {
    // Gen Info
    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("nMC", 1);
    fChain->SetBranchStatus("mc*", 1);
    if (fChain->GetBranch("nMC"))
      fChain->SetBranchAddress("nMC", &fMuons.Gen_nMC, &fMuons.b_Gen_nMC);
    if (fChain->GetBranch("mcPID"))
      fChain->SetBranchAddress("mcPID", &fMuons.Gen_pid, &fMuons.b_Gen_pid);
    if (fChain->GetBranch("mcStatus"))
      fChain->SetBranchAddress("mcStatus", &fMuons.Gen_status, &fMuons.b_Gen_status);
    if (fChain->GetBranch("mcMass"))
      fChain->SetBranchAddress("mcMass", &fMuons.Gen_mass, &fMuons.b_Gen_mass);
     if (fChain->GetBranch("mcE"))
      fChain->SetBranchAddress("mcE", &fMuons.Gen_e, &fMuons.b_Gen_e);
    if (fChain->GetBranch("mcPt"))
      fChain->SetBranchAddress("mcPt", &fMuons.Gen_pt, &fMuons.b_Gen_pt);
    if (fChain->GetBranch("mcEta"))
      fChain->SetBranchAddress("mcEta", &fMuons.Gen_eta, &fMuons.b_Gen_eta);
    if (fChain->GetBranch("mcPhi"))
      fChain->SetBranchAddress("mcPhi", &fMuons.Gen_phi, &fMuons.b_Gen_phi);
    if (fChain->GetBranch("mcMomPID"))
      fChain->SetBranchAddress("mcMomPID", &fMuons.Gen_mompid, &fMuons.b_Gen_mompid);
    if (fChain->GetBranch("mcMommass"))
      fChain->SetBranchAddress("mcMommass", &fMuons.Gen_mommass, &fMuons.b_Gen_mommass);
    if (fChain->GetBranch("mcMomE"))
      fChain->SetBranchAddress("mcMomE", &fMuons.Gen_mome, &fMuons.b_Gen_mome);
    if (fChain->GetBranch("mcMomPt"))
      fChain->SetBranchAddress("mcMomPt", &fMuons.Gen_mompt, &fMuons.b_Gen_mompt);
    if (fChain->GetBranch("mcMomEta"))
      fChain->SetBranchAddress("mcMomEta", &fMuons.Gen_mometa, &fMuons.b_Gen_mometa);
    if (fChain->GetBranch("mcMomPhi"))
      fChain->SetBranchAddress("mcMomPhi", &fMuons.Gen_momphi, &fMuons.b_Gen_momphi);
     if (fChain->GetBranch("mcGMomPID"))
      fChain->SetBranchAddress("mcGMomPID", &fMuons.Gen_Gmompid, &fMuons.b_Gen_Gmompid);

    // Reco Info
    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("nMu", 1);
    fChain->SetBranchStatus("mu*", 1);
    if (fChain->GetBranch("nMu"))
      fChain->SetBranchAddress("nMu", &fMuons.Glb_nptl, &fMuons.b_Glb_nptl);
    if (fChain->GetBranch("muCharge"))
      fChain->SetBranchAddress("muCharge", &fMuons.Glb_charge, &fMuons.b_Glb_charge);
    if (fChain->GetBranch("muPt"))
      fChain->SetBranchAddress("muPt", &fMuons.Glb_pt, &fMuons.b_Glb_pt);
    if (fChain->GetBranch("muEta"))
      fChain->SetBranchAddress("muEta", &fMuons.Glb_eta, &fMuons.b_Glb_eta);
    if (fChain->GetBranch("muPhi"))
      fChain->SetBranchAddress("muPhi", &fMuons.Glb_phi, &fMuons.b_Glb_phi);
    if (fChain->GetBranch("muD0"))
      fChain->SetBranchAddress("muD0", &fMuons.Glb_dxy, &fMuons.b_Glb_dxy);
    if (fChain->GetBranch("muDz"))
      fChain->SetBranchAddress("muDz", &fMuons.Glb_dz, &fMuons.b_Glb_dz);
    if (fChain->GetBranch("muPixelHits"))
      fChain->SetBranchAddress("muPixelHits", &fMuons.Glb_nValPixHits, &fMuons.b_Glb_nValPixHits);
    //if (fChain->GetBranch("Glb_nValTrkHits"))
    //fChain->SetBranchAddress("Glb_nValTrkHits", &fMuons.Glb_nValTrkHits, &fMuons.b_Glb_nValTrkHits);
    if (fChain->GetBranch("muMuonHits"))
      fChain->SetBranchAddress("muMuonHits", &fMuons.Glb_nValMuHits, &fMuons.b_Glb_nValMuHits);
    if (fChain->GetBranch("muTrkQuality"))
      fChain->SetBranchAddress("muTrkQuality", &fMuons.Glb_TrkQuality, &fMuons.b_Glb_TrkQuality);
    if (fChain->GetBranch("muIsGood"))
      fChain->SetBranchAddress("muIsGood", &fMuons.Glb_isGood, &fMuons.b_Glb_isGood);
    if (fChain->GetBranch("muChi2NDF"))
      fChain->SetBranchAddress("muChi2NDF", &fMuons.Glb_glbChi2_ndof, &fMuons.b_Glb_glbChi2_ndof);
    if (fChain->GetBranch("muStations"))
      fChain->SetBranchAddress("muStations", &fMuons.Glb_nMatchedStations, &fMuons.b_Glb_nMatchedStations);
    if (fChain->GetBranch("muInnerD0"))
      fChain->SetBranchAddress("muInnerD0", &fMuons.Glb_TrkDxy, &fMuons.b_Glb_TrkDxy);
    if (fChain->GetBranch("muInnerDz"))
      fChain->SetBranchAddress("muInnerDz", &fMuons.Glb_TrkDz, &fMuons.b_Glb_TrkDz);
    if (fChain->GetBranch("muPixelLayers"))
      fChain->SetBranchAddress("muPixelLayers", &fMuons.Glb_pixLayerWMeas, &fMuons.b_Glb_pixLayerWMeas);
    if (fChain->GetBranch("muTrkLayers"))
      fChain->SetBranchAddress("muTrkLayers", &fMuons.Glb_TrkLayerWMeas, &fMuons.b_Glb_TrkLayerWMeas);
    if (fChain->GetBranch("muPFChIso"))
      fChain->SetBranchAddress("muPFChIso", &fMuons.Glb_pfChIso, &fMuons.b_Glb_pfChIso);
    if (fChain->GetBranch("muPFPhoIso"))
      fChain->SetBranchAddress("muPFPhoIso", &fMuons.Glb_pfPhoIso, &fMuons.b_Glb_pfPhoIso);
    if (fChain->GetBranch("muPFNeuIso"))
      fChain->SetBranchAddress("muPFNeuIso", &fMuons.Glb_pfNeuIso, &fMuons.b_Glb_pfNeuIso);
    if (fChain->GetBranch("muPFPUIso"))
      fChain->SetBranchAddress("muPFPUIso", &fMuons.Glb_pfPUIso, &fMuons.b_Glb_pfPUIso);
    
    fInit = kTRUE;
  }
  return kTRUE;
}

//__________________________________________________________
Bool_t lwMuonProducer::InitEventObjects() {

  //Create event objects
  if(!fEventObjects) {
    Printf("%s: fEventObjects does not exist. Cannot store output",GetName());
    return kFALSE;
  } else {
    if(!fEventObjects->FindObject(flwMuonsRecoName)) {
      flwMuonsReco = new TClonesArray("lwMuon");
      flwMuonsReco->SetName(flwMuonsRecoName);
      fEventObjects->Add(flwMuonsReco);
    }
    if(!fEventObjects->FindObject(flwMuonsGeneName) && !flwMuonsGeneName.IsNull()) {
      flwMuonsGene = new TClonesArray("genParticle");
      flwMuonsGene->SetName(flwMuonsGeneName);
      fEventObjects->Add(flwMuonsGene);
    }
  }
  
  return kTRUE;
}

//__________________________________________________________
Bool_t lwMuonProducer::Run(Long64_t entry) {

  //overloaded run funtion
  Long64_t centry = LoadTree(entry);
  if(centry<0) return kFALSE;

  if(!InitEventObjects()) return kFALSE;
  
  //clear arrays
  flwMuonsReco->Delete();
  if(flwMuonsGene) flwMuonsGene->Delete();

  //reconstructed muons
  Int_t muCount = 0;
  for(Int_t i = 0; i<fMuons.Glb_nptl; i++) {
    if(!AcceptMuon(i)) continue;
    lwMuon *mu = new lwMuon(fMuons.Glb_pt->at(i),
                            fMuons.Glb_eta->at(i),
                            fMuons.Glb_phi->at(i),
                            0,
                            i);
    mu->SetCharge(fMuons.Glb_charge->at(i));
    (*flwMuonsReco)[muCount] = mu;
    ++muCount;
  }
  flwMuonsReco->Sort();
  //Printf("%d reconstructed muons",muCount);

  //generated muons
  if(flwMuonsGene) {
    muCount = 0;
    for(Int_t i = 0; i<fMuons.Gen_nMC; i++) {
      genParticle *mu = new genParticle(fMuons.Gen_pt->at(i),
                                        fMuons.Gen_eta->at(i),
                                        fMuons.Gen_phi->at(i),
                                        0,
                                        i);
      mu->SetCharge(fMuons.Gen_pid->at(i)/abs(fMuons.Gen_pid->at(i)));
      mu->SetPID(fMuons.Gen_pid->at(i));
      mu->SetPIDMom(fMuons.Gen_mompid->at(i));
      mu->SetPIDGMom(fMuons.Gen_Gmompid->at(i));
      (*flwMuonsGene)[muCount] = mu;
      ++muCount;
    }
    flwMuonsGene->Sort();
    //Printf("%d generated muons",muCount);
  }
  
  return kTRUE;
}

//__________________________________________________________
Bool_t lwMuonProducer::AcceptMuon(Int_t i) {

  //muon quality selection 
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  //plus https://github.com/CmsHI/quickZMacros/blob/master/ggHistos.C#L4
  if(!(fMuons.Glb_isGood->at(i)))                       return kFALSE;
  else if((fMuons.Glb_pt->at(i))<fPtMin)                return kFALSE;
  else if((fabs(fMuons.Glb_eta->at(i)))>fMaxEtaAbs)     return kFALSE;
  //else if((fMuons.Glb_trkChi2_ndof->at(i))>fMaxTrkChi2) return kFALSE;
  else if((fMuons.Glb_glbChi2_ndof->at(i))>fMaxGlbChi2) return kFALSE;
  else if((fMuons.Glb_nValMuHits->at(i))<fMinNMuHits)   return kFALSE;
  else if((fMuons.Glb_nMatchedStations->at(i))<fMinMS)  return kFALSE;
  //else if((fMuons.Glb_dxy->at(i))>fMaxDxy)            return kFALSE;
  //else if((fMuons.Glb_dz->at(i))>fMaxDz)              return kFALSE;
  else if((fabs(fMuons.Glb_TrkDxy->at(i)))>fMaxDxy)     return kFALSE;
  else if((fabs(fMuons.Glb_TrkDz->at(i)))>fMaxDz)       return kFALSE;
  else if((fMuons.Glb_nValPixHits->at(i))<fMinNPixHits) return kFALSE;
  else if((fMuons.Glb_TrkLayerWMeas->at(i))<fMinTrkLWM) return kFALSE;
  else if (!(fMuons.Glb_TrkQuality->at(i)))             return kFALSE; 
  else return kTRUE;
}

//__________________________________________________________
Long64_t lwMuonProducer::LoadTree(Long64_t entry) {

  //overloaded LoadTree function 
  if(!fChain) {
    Printf("fChain doesn't exist");
    return -5;
  }
  
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Init();
    Printf("%lld fCurrent: %d",entry,fCurrent);
  }

  //  fChain->SetMakeClass(1);
 
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) {
    Printf("hiEventProducer: centry smaller than 0");
    return centry;  
  }
  
  fChain->GetEntry(entry);

  return centry;
}
