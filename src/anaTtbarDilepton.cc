
#include "UserCode/diall/interface/anaTtbarDilepton.h"
#include "UserCode/diall/interface/diParticle.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/lwJetContainer.h"
#include "UserCode/diall/interface/particleBase.h"
#include "UserCode/diall/interface/pfParticle.h"


#include "TLorentzVector.h"
#include "TMath.h"

#include "TClass.h"

anaTtbarDilepton::anaTtbarDilepton(const char *name, const char *title) 
:anaBaseTask(name,title),
  fCheckPid(kFALSE),
  fMetType(),      
  fMinPt(0.), 
  fTriggerMapName(""),
  fRecoLeptonLeadName(""),
  fRecoLeptonSubleadName(""),  
  fGenLeptonName(""),
  fGenJetsName(""),
  fRecoJetsName(""),
  fParticlesName(""),
  fDileptonChargeSign(-1),
  fDileptonName(""),
  fTopsName(""),
  fTriggerMap(0x0),
  fRecoLeptonLead(0x0), 
  fRecoLeptonSublead(0x0),     
  fGenLepton(0x0),
  fGenJets(0x0),
  fRecoJets(0x0),
  fParticles(0x0),
  fDilepton(0x0),
  fTops(0x0),
  fDileptonCands(0),
  fNEvents(),
  fAnaTtbarDileptonInfo(0x0)

{

  //for(Int_t j = 0; j<10; ++j)
  // fh2MetCentPtMin[j] = 0;
  
}


//----------------------------------------------------------
void anaTtbarDilepton::Exec(Option_t * /*option*/)
{

   //printf("anaTtbarDilepton executing\n");
   //if(!SelectEvent()) return;

   TLorentzVector met;
   fDileptonCands = 0; 
   if(!fInitOutput) CreateOutputObjects();

   fNEvents->Fill(1);

   //Get objects from event

   //Get event container
   fHiEvent = dynamic_cast<hiEventContainer*>(fEventObjects->FindObject(fEvtName.Data()));
   if(!fHiEvent) { Printf("No %s HiEventContainer found", fEvtName.Data()); return; }

   fEventWeight.push_back(fHiEvent -> GetWeight());
   
   //Get triggers
   if(!fTriggerMap && !fTriggerMapName.IsNull()) {
     fTriggerMap = dynamic_cast<triggerMap*>(fEventObjects->FindObject(fTriggerMapName.Data()));
   }
   //if(!fTriggerMap) { Printf("No %s TriggerMap found", fTriggerMapName.Data());  }

   //if(!fTriggerMap->TriggerFired("HLT_HIL3Mu15_v1")) return;
   //std::cout<< fTriggerMap->TriggerFired("HLT_HIL3Mu15_v1") <<std::endl;
   //fTriggerMap->PrintTriggers();
      
   fNEvents->Fill(2); //only for now!
   fNEvents->Fill(3); //only for now!

   //Get leading lepton
   if(!fRecoLeptonLead && !fRecoLeptonLeadName.IsNull()) {
     fRecoLeptonLead = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fRecoLeptonLeadName.Data()));
   }
   if(!fRecoLeptonLead) { Printf("No %s RecoLeptonLead found", fRecoLeptonLeadName.Data()); return; }

   const Int_t nRecoLeptonLead = fRecoLeptonLead->GetEntriesFast();
   fNRecoLeptonLead.push_back(nRecoLeptonLead);

   //Get subleading lepton
   if(!fRecoLeptonSublead && !fRecoLeptonSubleadName.IsNull()) {
     fRecoLeptonSublead = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fRecoLeptonSubleadName.Data()));
   }
   if(!fRecoLeptonSublead) { Printf("No %s RecoLeptonSublead found", fRecoLeptonSubleadName.Data()); return; } 

   const Int_t nRecoLeptonSublead = fRecoLeptonSublead->GetEntriesFast();
   fNRecoLeptonSublead.push_back(nRecoLeptonSublead);
   
   //Get gen lepton
   if(!fGenLepton && !fGenLeptonName.IsNull()) {
     fGenLepton = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fGenLeptonName.Data()));
   }
   if(!fGenLepton) { Printf("No %s GenLepton found", fGenLeptonName.Data()); return; } 
   //const Int_t nGenLeptons = fGenLepton->GetEntriesFast();
   //std::cout << "nleptons "<< nGenLeptons << std::endl;
   
   //Make array for emu candidates
   if(!fEventObjects->FindObject(fDileptonName) && !fDilepton) {
      fDilepton = new TClonesArray("diParticle");
      fDilepton->SetName(fDileptonName);
      fEventObjects->Add(fDilepton);
    }
   if(fDilepton) fDilepton->Delete();


   //Make array for top quark candidates
   //if(!fEventObjects->FindObject(fTopsName) && !fTops) {
   //   fTops = new TClonesArray("diParticle");
   //   fTops->SetName(fTopsName);
   //   fEventObjects->Add(fTops);
   // }
   //if(fTops) fTops->Delete();
   
   //Get reco jets
   if(!fRecoJets && !fRecoJetsName.IsNull())
     fRecoJets = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fRecoJetsName.Data()));
   if(!fRecoJets) { Printf("No %s RecoJets found", fRecoJetsName.Data()); return; }
   const Int_t nRecoJets= fRecoJets->GetNJets();
   fNJetsIncl.push_back(nRecoJets);
   
   //Get gen jets
   if(!fGenJets && !fGenJetsName.IsNull())
     fGenJets = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fGenJetsName.Data()));
   if(!fGenJets) { Printf("No %s GenJets found", fGenJetsName.Data()); return; }
   const Int_t nGenJets= fGenJets->GetNJets();
     
   for (int i = 0; i < nGenJets; i++) {
     lwJet * jet = fGenJets->GetJet(i);
     std::cout<<"disc " <<jet->Pt()<<std::endl;
   }
   
   FillDileptonArray(nRecoLeptonLead, nRecoLeptonSublead);


   Int_t nbjets = 0;
   for (int i = 0; i < fRecoJets->GetNJets(); i++) {
     lwJet * jet = fRecoJets->GetJet(i);
     fBJetsDiscr.push_back(jet->GetCsvSimpleDiscr());
     if(jet->GetCsvSimpleDiscr()>0.9)nbjets++;
     //std::cout<<"disc" <<jet->GetCsvSimpleDiscr()<<std::endl;
   }
   fNBJetsIncl.push_back(nbjets);
   
   //Get particles from which MET will be calculated
   if(!fParticles && !fParticlesName.IsNull()) {
     //fEventObjects->Print();
     fParticles = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fParticlesName.Data()));
     //if(!fParticles) {
       //check if in jet branch
       //lwJetContainer *jetsCont = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fParticlesName.Data()));
       //if(jetsCont) fParticles = jetsCont->GetJets();
     //}
   }
   
   if(!fParticles) {
     Printf("%s: WARNING: Couldn't locate %s branch",GetName(),fParticlesName.Data());
     return;
   }
   
   /*   
   //Double_t cent = 5.;//fHiEvent->GetCentrality();
   TLorentzVector p4(0.,0.,0.,0.);
   //Double_t sumEt = 0.;

   const Int_t nptmin = 10;
   //Double_t ptarr[nptmin] {0.,1.,2.,3.,10.,20.,30.,40.,50.,60.};
   TLorentzVector r4[nptmin];
   for(Int_t j = 0; j<nptmin; ++j)
     r4[j].SetPtEtaPhiM(0.,0.,0.,0.);
   for (int i = 0; i < fParticles->GetEntriesFast(); i++) {
     particleBase *p = static_cast<particleBase*>(fParticles->At(i));
     if(!p) {
       Printf("%s ERROR: couldn't get particle",GetName());
       continue;
     }
     
     TLorentzVector l;
     if(fMetType==kGen || fMetType==kPFRaw) {
     //  if(p->Pt() < fMinPt)
       //  continue;
       l = p->GetLorentzVector();
       printf("!! mpika edp %f\n", fMinPt);

     }
     else if(fMetType==kVS) {
       pfParticle *pf = dynamic_cast<pfParticle*>(p);
       if(!pf) {
         Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
         return;
       }
       l.SetPtEtaPhiM(pf->PtVS(),pf->Eta(),pf->Phi(),pf->M());
     }
     else if(fMetType==kPuppi) {
       pfParticle *pf = dynamic_cast<pfParticle*>(p);
       if(!pf) {
         Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
         return;
       }
       if (pf->GetPuppiWeight()!=0) printf("weight %f\n", pf->GetPuppiWeight());
       l = pf->GetPuppiWeight()*p->GetLorentzVector();
     }
     
     for(Int_t j = 0; j<nptmin; ++j) {
       if(l.Pt()>ptarr[j])
         r4[j]+=l;
     }
    
     if(l.Pt() < fMinPt) continue;
     fh3PtEtaPhi->Fill(l.Pt(),l.Eta(),l.Phi());
     p4+=l;
     sumEt+=l.Et();
 
   }//particle loop

   met = -p4;
   //fh2MetCent->Fill(cent,met.Pt());
   //fh2SumEtCent->Fill(cent,sumEt);

   //Int_t nhists = TMath::Min(nptmin,10);
   //for(Int_t j = 0; j<nhists; ++j) {
   //  TLorentzVector met2 = -r4[j];
   //  fh2MetCentPtMin[j]->Fill(cent,met2.Pt());
   //}
   Printf("muon loop %d\n", fLeptonSubleads->GetEntriesFast()-1);
     */

   
   //int error;
   //fuParaZllPt->Fill(fParticles->GetEntriesFast());
   Bool_t filled = false; Bool_t filled_Rec2jets = false;
   for(int i = 0; i<fDileptonCands; ++i) 
   { 
       diParticle *pPart = (diParticle*)fDilepton->At(i);
       particleBase *pPart1 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(0));
       particleBase *pPart2 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(1));
	  
       bool dr_flag = false; 
       for (int ijet = 0; ijet < fRecoJets->GetNJets(); ijet++) 
       {
	   lwJet * jet = fRecoJets->GetJet(i);
	   if(jet->DeltaR(pPart1)<0.4 || jet->DeltaR(pPart2)<0.4) dr_flag = true;
	   
       }
       if(!dr_flag)
       {
	 float ht = 0;
	 if(pPart->M()>20)
	 {
	     std::cout<< " meta "<< nbjets<<std::endl;
	     fMassDilepton.push_back(pPart->M());
	     if (pPart1->Pt()>pPart2->Pt()) 
	     {
	         fLeadRecoLeptonPt.push_back(pPart1->Pt());
		 fLeadRecoLeptonAbsEta.push_back(TMath::Abs(pPart1->Eta()));
		 fLeadRecoLeptonEta.push_back(pPart1->Eta());
	     }
	     else 
	     {
	         fLeadRecoLeptonPt.push_back(pPart2->Pt());
		 fLeadRecoLeptonAbsEta.push_back(TMath::Abs(pPart2->Eta()));
		 fLeadRecoLeptonEta.push_back(pPart2->Eta());
	     }
	     if(!filled)
	     {
	         fNEvents->Fill(4);
		 for (int jet = 0; jet < fRecoJets->GetNJets(); jet++) 
		 {
		     float pt = fRecoJets->GetJet(jet)->Pt();
		     if (jet==0)
		     {
		         fLeadJetPt.push_back(pt);
		     }
		     ht += fRecoJets->GetJet(jet)->Pt();
		 }
		 fHT.push_back(ht);
		 fNJets.push_back(nRecoJets);
		 if (nRecoJets>=2)fNEvents->Fill(5);
		 if (nbjets>=1)fNEvents->Fill(6);
		 switch(nbjets)
		 {
		     case 1:
		       fNEvents->Fill(7);
		       break;
		     case 2:
		       fNEvents->Fill(8);
		       break;
		 }
		 filled=true;
	     }
	     if(nRecoJets>=2)
	     {
	         if(!filled_Rec2jets)
		 {
		     fLeadJetPt_Rec2jets.push_back(fRecoJets->GetJet(0)->Pt());
		     fLeadJetAbsEta_Rec2jets.push_back(TMath::Abs(fRecoJets->GetJet(0)->Eta()));
		     fLeadJetEta_Rec2jets.push_back(fRecoJets->GetJet(0)->Eta());
		     filled_Rec2jets=true;
		 }
		    
		 fDeltaPhi_Rec2jets.push_back(TMath::Abs(pPart1->DeltaPhi(pPart2))/TMath::Pi());
		 fMassDilepton_Rec2jets.push_back(pPart->M());
		 if (pPart1->Pt()>pPart2->Pt()) 
		 {
		     fLeadRecoLeptonPt_Rec2jets.push_back(pPart1->Pt());
		     fLeadRecoLeptonAbsEta_Rec2jets.push_back(TMath::Abs(pPart1->Eta()));
		     fLeadRecoLeptonEta_Rec2jets.push_back(pPart1->Eta());
		 }
		 else 
		 {
		     fLeadRecoLeptonPt_Rec2jets.push_back(pPart2->Pt());
		     fLeadRecoLeptonAbsEta_Rec2jets.push_back(TMath::Abs(pPart2->Eta()));
		     fLeadRecoLeptonEta_Rec2jets.push_back(pPart2->Eta());
		 }
	     }
	 }
       }
   }
   
     fAnaTtbarDileptonInfo->Fill();
     fEventWeight.clear();
     fNRecoLeptonLead.clear();       
     fNRecoLeptonSublead.clear();             
     fNJetsIncl.clear();          
     fNBJetsIncl.clear();  
     fBJetsDiscr.clear();
     fHT.clear();                 
     fMassDilepton.clear();            
     fLeadJetPt.clear();          
     fLeadRecoLeptonPt.clear();       
     fLeadRecoLeptonAbsEta.clear();   
     fLeadRecoLeptonEta.clear();      
     fNJets.clear();              
   
     fDeltaPhi_Rec2jets.clear();  
     fMassDilepton_Rec2jets.clear();   
     fLeadJetPt_Rec2jets.clear(); 
     fLeadJetAbsEta_Rec2jets.clear();
     fLeadJetEta_Rec2jets.clear();   
     fLeadRecoLeptonPt_Rec2jets.clear();
     fLeadRecoLeptonAbsEta_Rec2jets.clear();
     fLeadRecoLeptonEta_Rec2jets.clear();   

     //}
     //if (pPart->Pt()>20 && pPart->M()>60 && pPart->M()<120)
     //std::pair<double, double> u = compHadronicRecoilProjU(pPart,met, error, count);
     //if (pPart->Pt()>48 && pPart->Pt()<60)
     //fuParaZllPt->Fill(std::get<1>(u)+pPart->Pt());
     //for (int k = 0; k < fRecoJets->GetNJets(); k++) {
     // acos(cos(muPhi[0]-jtphi[0]))>0.2 && acos(cos(muPhi[0]-jtphi[1]))>0.2 && acos(cos(elePhi[0]-jtphi[0]))>0.2 && acos(cos(elePhi[0]-jtphi[1]))>0.2
     //	 }

     //if(pPart->M()>20)
     
     //printf("%f", pPart->Pt());
   
}
 

//----------------------------------------------------------                                                                         
bool anaTtbarDilepton::CheckPid(particleBase *p) {
  //check generated particle ID                                                                                            
  genParticle *gp = dynamic_cast<genParticle*>(p);
  if(!gp) {std::cout<< "edo generated! " << std::endl;return kFALSE;}
  if(abs(gp->GetPID())==13 || abs(gp->GetPID())==11) return kTRUE;
  return kFALSE;
}

std::pair<double, double> 
anaTtbarDilepton::compMETProjU(diParticle* zP4, double metPx, double metPy, int& errorFlag, int count)
{
  if ( zP4->Pt() == 0. ) {
    Warning ("compMEtProjU", " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!");
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4->Px();
  double qY = zP4->Py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  double uX = -metPx;
  double uY = -metPy;
  uX -= qX;
  uY -= qY;
  
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;

  return std::pair<double, double>(u1,u2);
}

std::pair<double, double> 
anaTtbarDilepton::compHadronicRecoilProjU(diParticle* zP4, TLorentzVector MET, int& errorFlag, int count)
{
  if ( zP4->Pt() == 0. ) {
    Warning ("compMEtProjU", " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!");
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4->Px();
  double qY = zP4->Py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  TLorentzVector zBoson = zP4->GetLorentzVector();
  TLorentzVector hadronicRecoil = -(MET+zBoson);

  double uX = hadronicRecoil.Px();
  double uY = hadronicRecoil.Py();
  uX -= qX;
  uY -= qY;
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;

  return std::pair<double, double>(u1,u2);
}

//----------------------------------------------------------
void 
anaTtbarDilepton::CreateOutputObjects() {

  anaBaseTask::CreateOutputObjects();

  if(!fOutput) {
    Printf("anaTtbarDilepton: fOutput not present");
    return;
  }
  
  fNEvents = new TH1F("fNEvents", "Entries per selection step", 8, 1,9);
  TString  step[8] = {"Initial", "Trigger", "PVfilter", "Dilepton", "#req 2jets","#req 1 b-tags","= 1 b-tag", "=2btag"};
  for (int i=1; i < fNEvents->GetNbinsX()+1; i++) 
    fNEvents->GetXaxis()->SetBinLabel(i,step[i-1]);

  fOutput->Add(fNEvents);

  fAnaTtbarDileptonInfo = new TTree("anaTtbarDileptonInfo","Tree with info related to ttbar emu analysis");
  fAnaTtbarDileptonInfo->Branch("fEventWeight", &fEventWeight);
  fAnaTtbarDileptonInfo->Branch("fNRecoLeptonLead", &fNRecoLeptonLead);
  fAnaTtbarDileptonInfo->Branch("fNRecoLeptonSublead",&fNRecoLeptonSublead);
  fAnaTtbarDileptonInfo->Branch("fNJetsIncl", &fNJetsIncl);
  fAnaTtbarDileptonInfo->Branch("fNBJetsIncl",&fNBJetsIncl);
  fAnaTtbarDileptonInfo->Branch("fBJetsDiscr",&fBJetsDiscr);

  fAnaTtbarDileptonInfo->Branch("fHT", &fHT);          
  fAnaTtbarDileptonInfo->Branch("fMassDilepton",&fMassDilepton);
  fAnaTtbarDileptonInfo->Branch("fLeadJetPt",&fLeadJetPt);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonPt",&fLeadRecoLeptonPt);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonAbsEta",&fLeadRecoLeptonAbsEta);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonEta",&fLeadRecoLeptonEta);     
  fAnaTtbarDileptonInfo->Branch("fNJets",&fNJets);

  fAnaTtbarDileptonInfo->Branch("fDeltaPhi_Rec2jets", &fDeltaPhi_Rec2jets);          
  fAnaTtbarDileptonInfo->Branch("fMassDilepton_Rec2jets",&fMassDilepton_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadJetPt_Rec2jets",&fLeadJetPt_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadJetAbsEta_Rec2jets",&fLeadJetAbsEta_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonPt_Rec2jets",&fLeadRecoLeptonPt_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonAbsEta_Rec2jets",&fLeadRecoLeptonAbsEta_Rec2jets);     
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonEta_Rec2jets",&fLeadRecoLeptonEta_Rec2jets);
  
  fOutput->Add(fAnaTtbarDileptonInfo);


  //for(Int_t j = 0; j<10; ++j) {
  //fh2MetCentPtMin[j] = new TH2F(Form("fh2MetCentPtMin%d",j),"fh2MetCent;centrality;MET",100,0,100,500,0,1000.);
  //fOutput->Add(fh2MetCentPtMin[j]);
  //}
    
  
}
//----------------------------------------------------------
void 
anaTtbarDilepton::FillDileptonArray(const Int_t nRecoLeptonLead, const Int_t nRecoLeptonSublead) {
  
  diParticle *pPart(0);
  for (int i = 0; i < nRecoLeptonLead; i++) {
     particleBase *leptonLead = static_cast<particleBase*>(fRecoLeptonLead->At(i));
     if(!leptonLead) {
       Printf("%s ERROR: couldn't get leading lepton",GetName());
       continue;
     }
     //CheckPid(leptonLead);
     int init = 0;
     if (fRecoLeptonLeadName.EqualTo(fRecoLeptonSubleadName))
       init = i + 1;

     //if(fCheckPid)
     //if(!CheckPid(leptonSublead1)) continue;
     for (int j = init; j < nRecoLeptonSublead; j++) {
       particleBase *leptonSublead = static_cast<particleBase*>(fRecoLeptonSublead->At(j));
       if(!leptonSublead) {
         Printf("%s ERROR: couldn't get subleading lepton",GetName());
         continue;
       }
       //CheckPid(leptonSublead);
       TLorentzVector l1 = leptonLead->GetLorentzVector();
       TLorentzVector l2 = leptonSublead->GetLorentzVector();
       TLorentzVector dilepton = l1 + l2;
       
       //dilepton pair should be of opposite sign
       if(leptonLead->GetCharge()*leptonSublead->GetCharge()==fDileptonChargeSign) {

         //if(fCheckPid)
	 //if(!CheckPid(mu2)) continue;
                   
         //fh3CentPtInvMass->Fill(cent,dimu.Pt(),dimu.M());
         
         //Store dilepton pair candidates in event
         if(fDilepton) {
           pPart = new ((*fDilepton)[fDileptonCands])
             diParticle(dilepton.Pt(),
                        dilepton.Eta(),
                        dilepton.Phi(),
                        dilepton.M(),
                        11); //dummy PDG
           pPart->SetCharge(0);
           pPart->AddParticle(leptonLead);
           pPart->AddParticle(leptonSublead);
           fDileptonCands++;
	 }
       }              
     }//subleading  loop
   }//leading lepton loop
}
