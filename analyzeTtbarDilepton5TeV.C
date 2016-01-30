#include "UserCode/diall/interface/hiEventProducer.h"
#include "UserCode/diall/interface/lwMuonProducer.h"
#include "UserCode/diall/interface/lwElectronProducer.h"
#include "UserCode/diall/interface/pfParticleProducer.h"
#include "UserCode/diall/interface/lwJetContainer.h"
#include "UserCode/diall/interface/lwJetFromForestProducer.h"
#include "UserCode/diall/interface/triggerProducer.h"
#include "UserCode/diall/interface/anaBaseTask.h"
#include "UserCode/diall/interface/anaTtbarDilepton.h"
#include "UserCode/diall/interface/anaJetMatching.h"
#include "UserCode/diall/interface/anaJetQA.h"
#include "UserCode/diall/interface/anaRhoProducer.h"
#include "UserCode/diall/interface/anaMET.h"
#include "UserCode/diall/interface/anaPFCandidates.h"
#include "UserCode/diall/interface/anaPuppiProducer.h"
#include "UserCode/diall/interface/anaPuppiParticles.h"

#include <TList.h>
#include <TChain.h>
#include <TFile.h>

//#include "CollectFiles.C"

#include <iostream>

using namespace std;

void analyzeTtbarDilepton5TeV(std::vector<std::string> urls, const char *outname = "AnaResults.root", Long64_t nentries = 20, Int_t firstF = -1, Int_t lastF = -1, Int_t firstEvent = 0, int isData = 1) {


  
  //Printf("anaFile: %d",anaFile);
  size_t firstFile = 0;
  size_t lastFile = urls.size();

  if(firstF>-1) {
    firstFile = (size_t)firstF;
    lastFile = (size_t)lastF;
  }
  std::cout << "firstFile: " << firstFile << "  lastFile: " << lastFile << std::endl;

  Int_t lastEvent = nentries;
  if(firstEvent>0) {
    lastEvent = firstEvent + nentries;
  }
  std::cout << "firstEvent: " << firstEvent << std::endl;
  std::cout << "isData: " << isData << std::endl;  

  //add files to chain
  TChain *chain = NULL;
  chain = new TChain("hiEvtAnalyzer/HiTree");
  for(size_t i=firstFile; i<lastFile; i++) chain->Add(urls[i].c_str());
  Printf("hltTree done");


  TChain *hiEvtTree = new TChain("hltanalysis/HltTree");
  for(size_t i=firstFile; i<lastFile; i++) hiEvtTree->Add(urls[i].c_str());
  chain->AddFriend(hiEvtTree);
  Printf("hiTree done");


  TChain *skimTree = new TChain("skimanalysis/HltTree");
  for(size_t i=firstFile; i<lastFile; i++) skimTree->Add(urls[i].c_str());
  chain->AddFriend(skimTree);
  Printf("skimTree done");

  TChain *pfTree = new TChain("pfcandAnalyzer/pfTree");
  for(size_t i=firstFile; i<lastFile; i++) pfTree->Add(urls[i].c_str());
  chain->AddFriend(pfTree);
  Printf("pfTree done");
  
  TChain *muTree = new TChain("ggHiNtuplizer/EventTree");
  for(size_t i=firstFile; i<lastFile; i++) muTree->Add(urls[i].c_str());
  chain->AddFriend(muTree);
  Printf("muTree done");

  TChain *pfJetTree = new TChain("ak4PFJetAnalyzer/t");
  for(size_t i=firstFile; i<lastFile; i++) pfJetTree->Add(urls[i].c_str());
  chain->AddFriend(pfJetTree);
  Printf("pfJetTree done");

  
  // TChain *genTree = new TChain("HiGenParticleAna/hi");
  // for(size_t i=firstFile; i<lastFile; i++) genTree->Add(urls[i].c_str());
  // chain->AddFriend(genTree);
  // Printf("genTree done");
  
  TList *fEventObjects = new TList();


  //---------------------------------------------------------------
  // producers
  //
  hiEventProducer *p_evt = new hiEventProducer("hiEvtProd");
  p_evt->SetInput(chain);
  p_evt->SetHIEventContName("hiEventContainer");
  p_evt->SetEventObjects(fEventObjects);

  triggerProducer *p_trg = new triggerProducer("trigProd");
  p_trg->SetInput(hiEvtTree);
  p_trg->SetTriggerMapName("triggerMap");
  p_trg->SetEventObjects(fEventObjects);
  
  pfParticleProducer *p_pf = new pfParticleProducer("pfPartProd");
  p_pf->SetInput(chain);
  p_pf->SetpfParticlesName("pfParticles");
  p_pf->SetEventObjects(fEventObjects);

  lwMuonProducer *p_mu = new lwMuonProducer("lwMuonProd");
  p_mu->SetInput(chain);
  p_mu->SetlwMuonsRecoName("lwMuonsReco");
  p_mu->SetlwMuonsGeneName("lwMuonsGene");
  p_mu->SetEventObjects(fEventObjects);

  lwElectronProducer *p_ele = new lwElectronProducer("lwElectronProd");
  p_ele->SetInput(chain);
  p_ele->SetlwElectronsRecoName("lwElectronsReco");
  //p_ele->SetlwElectronsGeneName("lwElectronsGene");
  p_ele->SetEventObjects(fEventObjects);


  lwJetFromForestProducer *p_pfJet = new lwJetFromForestProducer("lwJetForestProdPF");
  p_pfJet->SetInput(chain);
  p_pfJet->SetJetContName("akt4PF");
  p_pfJet->DoPFJetID(true);
  p_pfJet->SetGenJetContName("akt4Gen");
  p_pfJet->SetEventObjects(fEventObjects);
  p_pfJet->SetRadius(0.4);
  
  //---------------------------------------------------------------
  //analysis modules
  //

  //handler to which all modules will be added
  anaBaseTask *handler = new anaBaseTask("handler","handler");

  
  anaTtbarDilepton *anaEMu = new anaTtbarDilepton("anaEMu","anaEMu");
  anaEMu->ConnectEventObject(fEventObjects);
  anaEMu->SetHiEvtName("hiEventContainer");
  //anaEMu->SetTriggerMapName("triggerMap");
  anaEMu->SetParticlesName("pfParticles");
  anaEMu->SetRecoLeptonLeadName("lwElectronsReco");
  anaEMu->SetRecoLeptonSubleadName("lwMuonsReco");
  anaEMu->SetGenLeptonName("lwMuonsGene");
  anaEMu->SetRecoJetsName("akt4PF");
  anaEMu->SetGenJetsName("akt4Gen");
  //anaEMu->SetMetType(anaTtbarDilepton::kPFRaw);
  handler->Add(anaEMu);

  anaTtbarDilepton *anaMuMu = new anaTtbarDilepton("anaMuMu","anaMuMu");
  anaMuMu->ConnectEventObject(fEventObjects);
  anaMuMu->SetHiEvtName("hiEventContainer");
  //anaMuMu->SetTriggerMapName("triggerMap");
  anaMuMu->SetParticlesName("pfParticles");
  anaMuMu->SetRecoLeptonLeadName("lwMuonsReco");
  anaMuMu->SetRecoLeptonSubleadName("lwMuonsReco");
  anaMuMu->SetGenLeptonName("lwMuonsGene");
  anaMuMu->SetRecoJetsName("akt4PF");
  anaMuMu->SetGenJetsName("akt4Gen");
  //anaMuMu->SetMetType(anaTtbarDilepton::kPFRaw);
  handler->Add(anaMuMu);
  
  anaTtbarDilepton *anaEleEle = new anaTtbarDilepton("anaEleEle","anaEleEle");
  anaEleEle->ConnectEventObject(fEventObjects);
  anaEleEle->SetHiEvtName("hiEventContainer");
  //anaEleEle->SetTriggerMapName("triggerMap");
  anaEleEle->SetParticlesName("pfParticles");
  anaEleEle->SetRecoLeptonLeadName("lwElectronsReco");
  anaEleEle->SetRecoLeptonSubleadName("lwElectronsReco");
  anaEleEle->SetGenLeptonName("lwMuonsGene");
  anaEleEle->SetRecoJetsName("akt4PF");
  anaEleEle->SetGenJetsName("akt4Gen");
  //anaEleEle->SetMetType(anaTtbarDilepton::kPFRaw);
  handler->Add(anaEleEle);
  

  /*
  //PF-GEN matching
  anaJetMatching *matchingGenPFJet = new anaJetMatching("matchingGenPFJet","matchingGenPFJet");
  matchingGenPFJet->ConnectEventObject(fEventObjects);
  matchingGenPFJet->SetHiEvtName("");
  matchingGenPFJet->SetJetsNameBase("akt4Gen");
  matchingGenPFJet->SetJetsNameTag("akt4PF");
  matchingGenPFJet->SetNCentBins(1);
  handler->Add(matchingGenPFJet);  
  
  anaPFvsCaloJet *anaGenPFJet = new anaPFvsCaloJet("anaGenVsPFJet","anaGenVsPFJet");
  anaGenPFJet->ConnectEventObject(fEventObjects);
  anaGenPFJet->SetHiEvtName("");
  anaGenPFJet->SetPFJetsName("akt4Gen");
  anaGenPFJet->SetCaloJetsName("akt4PF");
  handler->Add(anaGenPFJet);
  */
  
  //---------------------------------------------------------------
  //Event loop
  //---------------------------------------------------------------	  
  Long64_t entries_tot =  chain->GetEntries(); 

  std::cout << entries_tot << std::endl;
  if(nentries<0) lastEvent = chain->GetEntries();  

  // Long64_t nentries = 20;//chain->GetEntriesFast();

  Printf("nentries: %lld  tot: %lld",nentries,entries_tot);
  for (Long64_t jentry=firstEvent; jentry<lastEvent; ++jentry) { 
    //  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) Printf("Processing event %d  %d",(int)(jentry), (int)(lastEvent));
    //Run producers
    if (jentry>=entries_tot) break;
    // Printf("produce hiEvent");
    if(!p_evt->Run(jentry));   //hi event properties

    //if(!p_trg->Run(jentry));
    if(!p_mu->Run(jentry)); //lepton collection
    if(!p_ele->Run(jentry)); //lepton collection
    if(!p_pf->Run(jentry));    //pf particles
    if(!p_pfJet->Run(jentry)); //jets
    
    //Execute all analysis tasks
    handler->ExecuteTask();
  }
  fEventObjects->Print();

  TFile *out = new TFile(outname,"RECREATE");
  TList *tasks = handler->GetListOfTasks();
  TIter next(tasks);
  anaBaseTask *obj;
  while ((obj = dynamic_cast<anaBaseTask*>(next()) ))
    if(obj->GetOutput()) obj->GetOutput()->Write(obj->GetName(),TObject::kSingleKey);
  out->Write();
  out->Close();

}
