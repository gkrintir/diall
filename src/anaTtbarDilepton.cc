
#include "UserCode/diall/interface/anaTtbarDilepton.h"
#include "UserCode/diall/interface/diParticle.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/lwElectron.h"
#include "UserCode/diall/interface/lwMuon.h"
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
  fMinMass(20.),
  fMinMET(30.),
  fTriggerMapName(""),
  fTriggerPath(""),
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


}

//----------------------------------------------------------
anaTtbarDilepton::~anaTtbarDilepton()
{

  //fRecoRemovalJets->Clear();
  
}
//----------------------------------------------------------
void anaTtbarDilepton::Exec(Option_t * /*option*/)
{

   //printf("anaTtbarDilepton executing\n");
   //if(!SelectEvent()) return;


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
   if(!fTriggerMap) { Printf("No %s TriggerMap found", fTriggerMapName.Data());  }

   //if(!fTriggerMap->TriggerFired("HLT_HIL3Mu15_v1") || !fTriggerMap->TriggerFired("HLT_HIL3Mu7_NHitQ15_v1")) return;
   //if(!fTriggerMap->TriggerFired("HLT_HIL3Mu15_v1")) return; 
   //if(!fTriggerMap->TriggerFired("HLT_HISinglePhoton10_Eta3p1_v1")) return;
   //if(!fTriggerMap->TriggerFired("HLT_HISinglePhoton15_Eta1p5_v1")) return;
   if (!fTriggerPath.IsNull()) if(!fTriggerMap->TriggerFired(fTriggerPath.Data())) return;
   //std::cout<< fTriggerMap->TriggerFired("HLT_HIL3Mu15_v1") <<std::endl;
   //fTriggerMap->PrintTriggers();
     
   fNEvents->Fill(2); //only for now!
   fNEvents->Fill(3); //only for now!

   //fTriggerFired.push_back(fTriggerMap->TriggerFired("HLT_HISinglePhoton20_Eta3p1_v1_Prescl"));

   //Get leading lepton
   if(!fRecoLeptonLead && !fRecoLeptonLeadName.IsNull()) {
     fRecoLeptonLead = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fRecoLeptonLeadName.Data()));
     if(!fRecoLeptonLead) { Printf("No %s RecoLeptonLead found", fRecoLeptonLeadName.Data()); return; }
   }


   const Int_t nRecoLeptonLead = fRecoLeptonLead->GetEntriesFast();
   fNRecoLeptonLead.push_back(nRecoLeptonLead);

   //Get subleading lepton
   if(!fRecoLeptonSublead && !fRecoLeptonSubleadName.IsNull()) {
     fRecoLeptonSublead = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fRecoLeptonSubleadName.Data()));
     if(!fRecoLeptonSublead) { Printf("No %s RecoLeptonSublead found", fRecoLeptonSubleadName.Data()); return; } 
   }


   const Int_t nRecoLeptonSublead = fRecoLeptonSublead->GetEntriesFast();
   fNRecoLeptonSublead.push_back(nRecoLeptonSublead);
   
   //Get gen lepton
   if(!fGenLepton && !fGenLeptonName.IsNull()) {
     fGenLepton = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fGenLeptonName.Data()));
     if(!fGenLepton) { Printf("No %s GenLepton found", fGenLeptonName.Data()); return; } 
   }
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
   if(!fRecoJets && !fRecoJetsName.IsNull()) {
     fRecoJets = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fRecoJetsName.Data()));
   if(!fRecoJets) { Printf("No %s RecoJets found", fRecoJetsName.Data()); return; 
   }}
   
   //Get gen jets
   if(!fGenJets && !fGenJetsName.IsNull()){
     fGenJets = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fGenJetsName.Data()));
     if(!fGenJets) { Printf("No %s GenJets found", fGenJetsName.Data()); return; }}
   //const Int_t nGenJets= fGenJets->GetNJets();
     
   //for (int i = 0; i < nGenJets; i++) {
   //  lwJet * jet = fGenJets->GetJet(i);
   //  std::cout<<"disc " <<jet->Pt()<<std::endl;
   // }
   FillDileptonArray(nRecoLeptonLead, nRecoLeptonSublead);
   fDilepton->Sort();
   fNDileptons.push_back(fDileptonCands);


   if (fDileptonCands>0)
   { 
       diParticle *pPart = (diParticle*)fDilepton->At(0);
       particleBase *pPart1 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(0));
       particleBase *pPart2 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(1));
       fSignDilepton_new.push_back(pPart->GetCharge());
       fMassDilepton_new.push_back(pPart->M());
       fEtaDilepton_new.push_back(pPart->Eta());
       fPtDilepton_new.push_back(pPart->Pt());
       fPhiDilepton_new.push_back(pPart->Phi());
       fDeltaPhi_new.push_back(TMath::Abs(pPart1->DeltaPhi(pPart2))/TMath::Pi());
       
       if (pPart1->Pt()>pPart2->Pt()) 
       {
	   fLeadRecoLeptonPt_new.push_back(pPart1->Pt());
	   fLeadRecoLeptonAbsEta_new.push_back(TMath::Abs(pPart1->Eta()));
	   fLeadRecoLeptonEta_new.push_back(pPart1->Eta());
       }
       else 
       {
	   fLeadRecoLeptonPt_new.push_back(pPart2->Pt());
	   fLeadRecoLeptonAbsEta_new.push_back(TMath::Abs(pPart2->Eta()));
	   fLeadRecoLeptonEta_new.push_back(pPart2->Eta());
       }
       if(pPart1->GetId()==11)
       {
	   lwElectron *electronLead = dynamic_cast<lwElectron*>(pPart1);
	   fElePtLead.push_back(pPart1->Pt());
	   fEleEtaLead.push_back(pPart1->Eta());
	   double iso, deposit ;
	   deposit =  electronLead->GetPFNeuIso() +  electronLead->GetPFPhoIso() - electronLead->GetEffAreaTimesRho() ;
	   iso = (electronLead->GetPFChIso() + TMath::Max (0.0, deposit ))/electronLead->Pt();

	   fEleIsoLead.push_back( iso );
	   deposit =  electronLead->GetPFNeuIso03() +  electronLead->GetPFPhoIso03() - electronLead->GetEffAreaTimesRho() ;
	   iso = (electronLead->GetPFChIso03() + TMath::Max (0.0, deposit ))/electronLead->Pt() ;
	   fEleIso03Lead.push_back( iso );
	   deposit =  electronLead->GetPFNeuIso04() +  electronLead->GetPFPhoIso04() - electronLead->GetEffAreaTimesRho() ;
	   iso = (electronLead->GetPFChIso04() + TMath::Max (0.0, deposit ))/electronLead->Pt() ;
	   fEleIso04Lead.push_back( iso  );
	   fdEtaAtVtxLead.push_back(electronLead->GetdEta());
	   fdPhiAtVtxLead.push_back(electronLead->GetdPhi());   
	   fSigmaIEtaIEtaLead.push_back(electronLead->GetSigmaIEtaIEta());
	   fHoverELead.push_back(electronLead->GetHoverE());       
	   fDxyLead.push_back(electronLead->GetD0());          
	   fDzLead.push_back(electronLead->GetDZ());           
	   fEoverPInvLead.push_back(electronLead->GetEoverPInv());    
	   fMissHitsLead.push_back(electronLead->GetMissHits());     
	   fConversionVetoLead.push_back(electronLead->GetConversionVeto());
       }
       if(pPart2->GetId()==11)
       {
	   lwElectron *electronSublead = dynamic_cast<lwElectron*>(pPart2);
	   fElePtSublead.push_back(pPart2->Pt());
	   fEleEtaSublead.push_back(pPart2->Eta());
	   double iso, deposit ;
	   deposit =  electronSublead->GetPFNeuIso() +  electronSublead->GetPFPhoIso() - electronSublead->GetEffAreaTimesRho() ;
	   iso = (electronSublead->GetPFChIso() + TMath::Max (0.0, deposit ))/electronSublead->Pt();

	   fEleIsoSublead.push_back( iso );
	   deposit =  electronSublead->GetPFNeuIso03() +  electronSublead->GetPFPhoIso03() - electronSublead->GetEffAreaTimesRho() ;
	   iso = (electronSublead->GetPFChIso03() + TMath::Max (0.0, deposit ))/electronSublead->Pt() ;
	   fEleIso03Sublead.push_back( iso );
	   deposit =  electronSublead->GetPFNeuIso04() +  electronSublead->GetPFPhoIso04() - electronSublead->GetEffAreaTimesRho() ;
	   iso = (electronSublead->GetPFChIso04() + TMath::Max (0.0, deposit ))/electronSublead->Pt() ;
	   fEleIso04Sublead.push_back( iso  );
	   
	   fdEtaAtVtxSublead.push_back(electronSublead->GetdEta());
	   fdPhiAtVtxSublead.push_back(electronSublead->GetdPhi());   
	   fSigmaIEtaIEtaSublead.push_back(electronSublead->GetSigmaIEtaIEta());
	   fHoverESublead.push_back(electronSublead->GetHoverE());       
	   fDxySublead.push_back(electronSublead->GetD0());          
	   fDzSublead.push_back(electronSublead->GetDZ());           
	   fEoverPInvSublead.push_back(electronSublead->GetEoverPInv());    
	   fMissHitsSublead.push_back(electronSublead->GetMissHits());     
	   fConversionVetoSublead.push_back(electronSublead->GetConversionVeto());
       }
       if(pPart1->GetId()==13)
       {
	   lwMuon *muonLead = dynamic_cast<lwMuon*>(pPart1);
	   fMuonPtLead.push_back(pPart1->Pt());
	   fMuonEtaLead.push_back(pPart1->Eta());
	   double iso, deposit ;
	   deposit =  muonLead->GetPFNeuIso() +  muonLead->GetPFPhoIso() - 0.5*muonLead->GetPFPUIso() ;
	   iso = (muonLead->GetPFChIso() + TMath::Max (0.0, deposit ))/muonLead->Pt();

	   fMuonIsoLead.push_back( iso );
       }
       if(pPart2->GetId()==13)
       {
	   lwMuon *muonSublead = dynamic_cast<lwMuon*>(pPart2);
	   fMuonPtSublead.push_back(pPart2->Pt());
	   fMuonEtaSublead.push_back(pPart2->Eta());
	   double iso, deposit ;
	   deposit =  muonSublead->GetPFNeuIso() +  muonSublead->GetPFPhoIso() - 0.5*muonSublead->GetPFPUIso() ;
	   iso = (muonSublead->GetPFChIso() + TMath::Max (0.0, deposit ))/muonSublead->Pt();
	   fMuonIsoSublead.push_back( iso );
	   
       }
   }
   else
   {
       fSignDilepton_new.push_back(-999);
       fMassDilepton_new.push_back(-999);
       fEtaDilepton_new.push_back(-999);
       fPtDilepton_new.push_back(-999);
       fPhiDilepton_new.push_back(-999);
       fDeltaPhi_new.push_back(-999);

       fLeadRecoLeptonPt_new.push_back(-999);
       fLeadRecoLeptonAbsEta_new.push_back(-999);
       fLeadRecoLeptonEta_new.push_back(-999);

       fElePtLead.push_back(-999);
       fEleEtaLead.push_back(-999);
       fEleIsoLead.push_back(-999);
       fEleIso03Lead.push_back(-999);
       fEleIso04Lead.push_back(-999);
       fdEtaAtVtxLead.push_back(-999);
       fdPhiAtVtxLead.push_back(-999);
       fSigmaIEtaIEtaLead.push_back(-999);
       fHoverELead.push_back(-999);
       fDxyLead.push_back(-999);
       fDzLead.push_back(-999);
       fEoverPInvLead.push_back(-999);
       fMissHitsLead.push_back(-999);
       fConversionVetoLead.push_back(-999);

       fElePtSublead.push_back(-999);
       fEleEtaSublead.push_back(-999);
       fEleIsoSublead.push_back(-999);
       fEleIso03Sublead.push_back(-999);
       fEleIso04Sublead.push_back(-999);
       fdEtaAtVtxSublead.push_back(-999);
       fdPhiAtVtxSublead.push_back(-999);
       fSigmaIEtaIEtaSublead.push_back(-999);
       fHoverESublead.push_back(-999);
       fDxySublead.push_back(-999);
       fDzSublead.push_back(-999);
       fEoverPInvSublead.push_back(-999);
       fMissHitsSublead.push_back(-999);
       fConversionVetoSublead.push_back(-999);
       
       fMuonPtLead.push_back(-999);
       fMuonEtaLead.push_back(-999);
       fMuonIsoLead.push_back(-999);
       
       fMuonPtSublead.push_back(-999);
       fMuonEtaSublead.push_back(-999);
       fMuonIsoSublead.push_back(-999);
	
   }
   
   for(int i = 0; i<fDileptonCands; ++i) 
    { 
	diParticle *pPart = (diParticle*)fDilepton->At(i);
	particleBase *pPart1 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(0));
	particleBase *pPart2 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(1));

	for (int ijet = 0; ijet < fRecoJets->GetNJets(); ijet++) 
	{
	    lwJet * jet = dynamic_cast<lwJet*>(fRecoJets->GetJet(ijet));
	    if(jet->DeltaR(pPart1)>0.4 && jet->DeltaR(pPart2)>0.4) 
	      fRecoRemovalJets.push_back(jet);
	}
    }

    const Int_t nRecoJets= fRecoRemovalJets.size();
    fNJetsIncl.push_back(nRecoJets);

    Int_t nbjets = 0;
    if (fRecoRemovalJets.size() > 0 )
    {
        fLeadJetPt_new.push_back(fRecoRemovalJets[0]->Pt());
	fLeadJetEta_new.push_back(fRecoRemovalJets[0]->Eta());
	for (unsigned int i = 0; i < fRecoRemovalJets.size(); i++) 
	{
	    lwJet * jet = dynamic_cast<lwJet*>(fRecoRemovalJets[i]);
	    fBJetsDiscr.push_back(jet->GetCsvV2Discr());
	    if(jet->GetCsvV2Discr()>0.9) ++nbjets;
	}
    }
    else 
    {
        fLeadJetPt_new.push_back(-999);
	fLeadJetEta_new.push_back(-999);
	fBJetsDiscr.push_back(-999);
    }

   fNBJetsIncl.push_back(nbjets);

   
   //Get particles from which MET will be calculated
   if(!fParticles && !fParticlesName.IsNull()) {
     fParticles = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fParticlesName.Data()));
     if(!fParticles) {
     Printf("%s: WARNING: Couldn't locate %s branch",GetName(),fParticlesName.Data());
     return;
     }
   }
   

   TLorentzVector met(0.,0.,0.,0.); 
   TLorentzVector p4(0.,0.,0.,0.);
   Double_t sumEt = 0.;
   
   TLorentzVector metNoHF(0.,0.,0.,0.);
   TLorentzVector p4NoHF(0.,0.,0.,0.);
   Double_t sumEtNoHF = 0.;

   for (int i = 0; i < fParticles->GetEntriesFast(); i++) {
     particleBase *p = static_cast<particleBase*>(fParticles->At(i));
     if(!p) {
       Printf("%s ERROR: couldn't get particle",GetName());
       continue;
     }
     
     TLorentzVector l;
     if(fMetType==kGen || fMetType==kPFRaw) 
     {
         //  if(p->Pt() < fMinPt)
         //  continue;
         l = p->GetLorentzVector();
     }
     else if(fMetType==kVS) 
     {
         pfParticle *pf = dynamic_cast<pfParticle*>(p);
	 if(!pf) 
	 {
	     Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
	     return;
	 }
	 l.SetPtEtaPhiM(pf->PtVS(),pf->Eta(),pf->Phi(),pf->M());
     }
     else if(fMetType==kPuppi) 
     {
         pfParticle *pf = dynamic_cast<pfParticle*>(p);
	 if(!pf) 
	 {
	     Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
	     return;
	 }
	 if (pf->GetPuppiWeight()!=0) printf("weight %f\n", pf->GetPuppiWeight());
	 l = pf->GetPuppiWeight()*p->GetLorentzVector();
     }
     
     //if(l.Pt() < fMinPt) continue;
     p4+=l;
     sumEt+=l.Et();
     if(abs(l.Eta()) <= 3.){
       p4NoHF+=l;
       sumEtNoHF+=l.Et();
     }
     
   }//particle loop
   
   met = -p4; 
   fMETAbs.push_back(met.Pt());
   fMETPhi.push_back(met.Phi());
   metNoHF = -p4;
   fMETAbsNoHF.push_back(metNoHF.Pt());
   fMETPhiNoHF.push_back(metNoHF.Phi());

   /*
   //int error;
   //fuParaZllPt->Fill(fParticles->GetEntriesFast());
   Bool_t filled = false; Bool_t filled_Rec2jets = false; Bool_t filled_Rec2jets_METCut = false;
   for(int i = 0; i<fDileptonCands; ++i) 
   { 
       diParticle *pPart = (diParticle*)fDilepton->At(i);
       particleBase *pPart1 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(0));
       particleBase *pPart2 = dynamic_cast<particleBase*>(pPart->GetDecayParticles()->At(1));

       bool dr_flag = false; 
       
       if(!dr_flag)
       {
	 float ht = 0;
	 if(pPart->M()>fMinMass)
	 {
	     fMassDilepton.push_back(pPart->M());
	     fEtaDilepton.push_back(pPart->Eta());
	     //if (pPart->Pt() < 20 ) std::cout << pPart1->Pt() << " 2nd "<< pPart2->Pt() << std::endl;
	     fPtDilepton.push_back(pPart->Pt());
	     fPhiDilepton.push_back(pPart->Phi());
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
	     //std::cout<< pPart1->GetId() <<std::endl;
	     if(pPart1->GetId()==11)
	     {
	         lwElectron *electronLead = dynamic_cast<lwElectron*>(pPart1);
		 fElePtLead.push_back(pPart1->Pt());
		 fEleEtaLead.push_back(pPart1->Eta());
		 fdEtaAtVtxLead.push_back(electronLead->GetdEta());
		 fdPhiAtVtxLead.push_back(electronLead->GetdPhi());   
		 fSigmaIEtaIEtaLead.push_back(electronLead->GetSigmaIEtaIEta());
		 fHoverELead.push_back(electronLead->GetHoverE());       
		 fDxyLead.push_back(electronLead->GetD0());          
		 fDzLead.push_back(electronLead->GetDZ());           
		 fEoverPInvLead.push_back(electronLead->GetEoverPInv());    
		 fMissHitsLead.push_back(electronLead->GetMissHits());     
		 fConversionVetoLead.push_back(electronLead->GetConversionVeto());
	     }
	     if(pPart2->GetId()==11)
	     {
		 lwElectron *electronSublead = dynamic_cast<lwElectron*>(pPart2);
		 fElePtSublead.push_back(pPart2->Pt());
		 fEleEtaSublead.push_back(pPart2->Eta());
		 fdEtaAtVtxSublead.push_back(electronSublead->GetdEta());
		 fdPhiAtVtxSublead.push_back(electronSublead->GetdPhi());   
		 fSigmaIEtaIEtaSublead.push_back(electronSublead->GetSigmaIEtaIEta());
		 fHoverESublead.push_back(electronSublead->GetHoverE());       
		 fDxySublead.push_back(electronSublead->GetD0());          
		 fDzSublead.push_back(electronSublead->GetDZ());           
		 fEoverPInvSublead.push_back(electronSublead->GetEoverPInv());    
		 fMissHitsSublead.push_back(electronSublead->GetMissHits());     
		 fConversionVetoSublead.push_back(electronSublead->GetConversionVeto());
	     }
	     TLorentzVector met(0.,0.,0.,0.); 

	     if(!filled)
	     { 
	         fNEvents->Fill(4);
		 
		 for (unsigned int ijet = 0; ijet < fRecoRemovalJets.size(); ijet++) 
		 {
		     float pt = fRecoRemovalJets[ijet]->Pt();
		     if (ijet==0)
		     {
		         fLeadJetPt.push_back(pt);
		     }
		     ht += fRecoRemovalJets[ijet]->Pt();
		 }
		 if (fRecoRemovalJets.size() == 0 ) fLeadJetPt.push_back(0.);
		 fHT.push_back(ht);
		 fNJets.push_back(nRecoJets);
		 if (nRecoJets>=2)fNEvents->Fill(5);
		 
		 //if (nbjets>=1)fNEvents->Fill(6);
		// switch(nbjets)
		// {
		//     case 1:
		//       fNEvents->Fill(7);
		//       break;
		//     case 2:
		//       fNEvents->Fill(8);
		//       break;
		// }

		 //Double_t cent = 5.;//fHiEvent->GetCentrality();
		 TLorentzVector p4(0.,0.,0.,0.);
		 Double_t sumEt = 0.;
		 
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
		     //printf("!! mpika edp %f\n", fMinPt);
		     
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
		   
		   //for(Int_t j = 0; j<nptmin; ++j) {
		   //  if(l.Pt()>ptarr[j])
		   //   r4[j]+=l;
		   // }
		   
		   //if(l.Pt() < fMinPt) continue;
		   if(abs(l.Eta()) <= 3.){
		     p4+=l;
		     sumEt+=l.Et();
		   }
		   
		 }//particle loop

		 met = -p4; 
		 fMETAbs.push_back(met.Pt());
		 fMETPhi.push_back(met.Phi());
		 if (met.Pt()>fMinMET) fNEvents->Fill(6);
		 filled=true;
		 
	     }
	     if(nRecoJets>=2)
	     {
	         
	         if(!filled_Rec2jets)
		 {
		     fLeadJetPt_Rec2jets.push_back(fRecoRemovalJets[0]->Pt());
		     fLeadJetAbsEta_Rec2jets.push_back(TMath::Abs(fRecoRemovalJets[0]->Eta()));
		     fLeadJetEta_Rec2jets.push_back(fRecoRemovalJets[0]->Eta());
		     filled_Rec2jets=true;
		 }
		    
		 fDeltaPhi_Rec2jets.push_back(TMath::Abs(pPart1->DeltaPhi(pPart2))/TMath::Pi());
		 fMassDilepton_Rec2jets.push_back(pPart->M());
		 fEtaDilepton_Rec2jets.push_back(pPart->Eta());
		 fPtDilepton_Rec2jets.push_back(pPart->Pt());
		 fPhiDilepton_Rec2jets.push_back(pPart->Phi());
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
		 if(met.Pt()>fMinMET)
		 {
		     if(!filled_Rec2jets_METCut)
		     {
		         fLeadJetPt_Rec2jets_METCut.push_back(fRecoRemovalJets[0]->Pt());
			 fLeadJetAbsEta_Rec2jets_METCut.push_back(TMath::Abs(fRecoRemovalJets[0]->Eta()));
			 fLeadJetEta_Rec2jets_METCut.push_back(fRecoRemovalJets[0]->Eta());
			 filled_Rec2jets_METCut=true;
		     }
		    
		     fDeltaPhi_Rec2jets_METCut.push_back(TMath::Abs(pPart1->DeltaPhi(pPart2))/TMath::Pi());
		     fMassDilepton_Rec2jets_METCut.push_back(pPart->M());
		     fEtaDilepton_Rec2jets_METCut.push_back(pPart->Eta());
		     fPtDilepton_Rec2jets_METCut.push_back(pPart->Pt());
		     fPhiDilepton_Rec2jets_METCut.push_back(pPart->Phi());
		     if (pPart1->Pt()>pPart2->Pt()) 
		     {
		         fLeadRecoLeptonPt_Rec2jets_METCut.push_back(pPart1->Pt());
			 fLeadRecoLeptonAbsEta_Rec2jets_METCut.push_back(TMath::Abs(pPart1->Eta()));
			 fLeadRecoLeptonEta_Rec2jets_METCut.push_back(pPart1->Eta());
		     }
		     else 
		     {
		         fLeadRecoLeptonPt_Rec2jets_METCut.push_back(pPart2->Pt());
			 fLeadRecoLeptonAbsEta_Rec2jets_METCut.push_back(TMath::Abs(pPart2->Eta()));
			 fLeadRecoLeptonEta_Rec2jets_METCut.push_back(pPart2->Eta());
		     }
		 }
	     }
	   }
       }
   }
   */
     fAnaTtbarDileptonInfo->Fill();
     fTriggerFired.clear();
     fRecoRemovalJets.clear();
     fEventWeight.clear();
     fNRecoLeptonLead.clear();       
     fNRecoLeptonSublead.clear();             
     fNJetsIncl.clear();          
     fNBJetsIncl.clear();  
     fBJetsDiscr.clear();
     fHT.clear();        
     fNDileptons.clear();
         
     fMassDilepton.clear();            
     fEtaDilepton.clear();
     fPtDilepton.clear();
     fPhiDilepton.clear();
     fLeadJetPt.clear();          
     fLeadRecoLeptonPt.clear();       
     fLeadRecoLeptonAbsEta.clear();   
     fLeadRecoLeptonEta.clear();      
     
     fMassDilepton_new.clear();            
     fSignDilepton_new.clear();
     fEtaDilepton_new.clear();
     fPtDilepton_new.clear();
     fPhiDilepton_new.clear();
     fLeadJetPt_new.clear();          
     fLeadJetEta_new.clear();
     fLeadRecoLeptonPt_new.clear();       
     fLeadRecoLeptonAbsEta_new.clear();   
     fLeadRecoLeptonEta_new.clear();      
     fDeltaPhi_new.clear();

     fMETAbs.clear();
     fMETPhi.clear();
     fMETAbsNoHF.clear();
     fMETPhiNoHF.clear();
     fNJets.clear();              
     
     fElePtLead.clear();
     fEleEtaLead.clear();
     fEleIsoLead.clear();
     fEleIso03Lead.clear();
     fEleIso04Lead.clear();
     fdEtaAtVtxLead.clear();
     fdPhiAtVtxLead.clear();
     fSigmaIEtaIEtaLead.clear();
     fHoverELead.clear();
     fDxyLead.clear();
     fDzLead.clear();
     fEoverPInvLead.clear();
     fMissHitsLead.clear();
     fConversionVetoLead.clear();

     fElePtSublead.clear();
     fEleEtaSublead.clear();
     fEleIsoSublead.clear();
     fEleIso03Sublead.clear();
     fEleIso04Sublead.clear();
     fdEtaAtVtxSublead.clear();
     fdPhiAtVtxSublead.clear();
     fSigmaIEtaIEtaSublead.clear();
     fHoverESublead.clear();
     fDxySublead.clear();
     fDzSublead.clear();
     fEoverPInvSublead.clear();
     fMissHitsSublead.clear();
     fConversionVetoSublead.clear();
     
     fDeltaPhi_Rec2jets.clear();  
     fMassDilepton_Rec2jets.clear();   
     fEtaDilepton_Rec2jets.clear(); 
     fPtDilepton_Rec2jets.clear(); 
     fPhiDilepton_Rec2jets.clear();
     fLeadJetPt_Rec2jets.clear(); 
     fLeadJetAbsEta_Rec2jets.clear();
     fLeadJetEta_Rec2jets.clear();   
     fLeadRecoLeptonPt_Rec2jets.clear();
     fLeadRecoLeptonAbsEta_Rec2jets.clear();
     fLeadRecoLeptonEta_Rec2jets.clear();  

     fDeltaPhi_Rec2jets_METCut.clear();  
     fMassDilepton_Rec2jets_METCut.clear();   
     fEtaDilepton_Rec2jets_METCut.clear(); 
     fPtDilepton_Rec2jets_METCut.clear(); 
     fPhiDilepton_Rec2jets_METCut.clear();
     fLeadJetPt_Rec2jets_METCut.clear(); 
     fLeadJetAbsEta_Rec2jets_METCut.clear();
     fLeadJetEta_Rec2jets_METCut.clear();   
     fLeadRecoLeptonPt_Rec2jets_METCut.clear();
     fLeadRecoLeptonAbsEta_Rec2jets_METCut.clear();
     fLeadRecoLeptonEta_Rec2jets_METCut.clear(); 
     

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
  
  fNEvents = new TH1F("fNEvents", "Entries per selection step", 9, 1, 10 );
  TString  step[9] = {"Initial", "Trigger", "PVfilter", "Dilepton", "#req 2jets", "MET > 30", "#req 1 b-tags","= 1 b-tag", "=2btag"};
  for (int i=1; i < fNEvents->GetNbinsX()+1; i++) 
    fNEvents->GetXaxis()->SetBinLabel(i,step[i-1]);

  fOutput->Add(fNEvents);

  fAnaTtbarDileptonInfo = new TTree("anaTtbarDileptonInfo","Tree with info related to ttbar emu analysis");
  fAnaTtbarDileptonInfo->Branch("fTriggerFired", &fTriggerFired);
  fAnaTtbarDileptonInfo->Branch("fEventWeight", &fEventWeight);
  fAnaTtbarDileptonInfo->Branch("fNRecoLeptonLead", &fNRecoLeptonLead);
  fAnaTtbarDileptonInfo->Branch("fNRecoLeptonSublead",&fNRecoLeptonSublead);
  fAnaTtbarDileptonInfo->Branch("fNJetsIncl", &fNJetsIncl);
  fAnaTtbarDileptonInfo->Branch("fNBJetsIncl",&fNBJetsIncl);
  fAnaTtbarDileptonInfo->Branch("fBJetsDiscr",&fBJetsDiscr);

  fAnaTtbarDileptonInfo->Branch("fHT", &fHT);          
  fAnaTtbarDileptonInfo->Branch("fNDileptons",&fNDileptons);
  fAnaTtbarDileptonInfo->Branch("fSignDilepton",&fSignDilepton_new);
  fAnaTtbarDileptonInfo->Branch("fMassDilepton",&fMassDilepton_new);
  fAnaTtbarDileptonInfo->Branch("fEtaDilepton",&fEtaDilepton_new);
  fAnaTtbarDileptonInfo->Branch("fPtDilepton",&fPtDilepton_new);
  fAnaTtbarDileptonInfo->Branch("fPhiDilepton",&fPhiDilepton_new);
  fAnaTtbarDileptonInfo->Branch("fDeltaPhiDilepton",&fDeltaPhi_new);
  fAnaTtbarDileptonInfo->Branch("fLeadJetPt",&fLeadJetPt_new);
  fAnaTtbarDileptonInfo->Branch("fLeadJetEta",&fLeadJetEta_new);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonPt",&fLeadRecoLeptonPt_new);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonAbsEta",&fLeadRecoLeptonAbsEta_new);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonEta",&fLeadRecoLeptonEta_new);     
  fAnaTtbarDileptonInfo->Branch("fMETAbs",&fMETAbs); 
  fAnaTtbarDileptonInfo->Branch("fMETPhi",&fMETPhi);
  fAnaTtbarDileptonInfo->Branch("fMETAbsNoHF",&fMETAbsNoHF);
  fAnaTtbarDileptonInfo->Branch("fMETPhiNoHF",&fMETPhiNoHF);
  fAnaTtbarDileptonInfo->Branch("fNJets",&fNJets);

  fAnaTtbarDileptonInfo->Branch("fElePtLead", &fElePtLead);
  fAnaTtbarDileptonInfo->Branch("fEleEtaLead", &fEleEtaLead);
  fAnaTtbarDileptonInfo->Branch("fEleIsoLead", &fEleIsoLead);
  fAnaTtbarDileptonInfo->Branch("fEleIso03Lead", &fEleIso03Lead);
  fAnaTtbarDileptonInfo->Branch("fEleIso04Lead", &fEleIso04Lead);
  fAnaTtbarDileptonInfo->Branch("fdEtaAtVtxLead", &fdEtaAtVtxLead);
  fAnaTtbarDileptonInfo->Branch("fdPhiAtVtxLead", &fdPhiAtVtxLead);
  fAnaTtbarDileptonInfo->Branch("fSigmaIEtaIEtaLead", &fSigmaIEtaIEtaLead);
  fAnaTtbarDileptonInfo->Branch("fHoverELead", &fHoverELead);
  fAnaTtbarDileptonInfo->Branch("fDxyLead", &fDxyLead);
  fAnaTtbarDileptonInfo->Branch("fDzLead", &fDzLead);
  fAnaTtbarDileptonInfo->Branch("fEoverPInvLead", &fEoverPInvLead);
  fAnaTtbarDileptonInfo->Branch("fMissHitsLead", &fMissHitsLead);
  fAnaTtbarDileptonInfo->Branch("fConversionVetoLead", &fConversionVetoLead);

  fAnaTtbarDileptonInfo->Branch("fElePtSublead", &fElePtSublead);
  fAnaTtbarDileptonInfo->Branch("fEleEtaSublead", &fEleEtaSublead);
  fAnaTtbarDileptonInfo->Branch("fEleIsoSublead", &fEleIsoSublead);
  fAnaTtbarDileptonInfo->Branch("fEleIso03Sublead", &fEleIso03Sublead);
  fAnaTtbarDileptonInfo->Branch("fEleIso04Sublead", &fEleIso04Sublead);
  fAnaTtbarDileptonInfo->Branch("fdEtaAtVtxSublead", &fdEtaAtVtxSublead);
  fAnaTtbarDileptonInfo->Branch("fdPhiAtVtxSublead", &fdPhiAtVtxSublead);
  fAnaTtbarDileptonInfo->Branch("fSigmaIEtaIEtaSublead", &fSigmaIEtaIEtaSublead);
  fAnaTtbarDileptonInfo->Branch("fHoverESublead", &fHoverESublead);
  fAnaTtbarDileptonInfo->Branch("fDxySublead", &fDxySublead);
  fAnaTtbarDileptonInfo->Branch("fDzSublead", &fDzSublead);
  fAnaTtbarDileptonInfo->Branch("fEoverPInvSublead", &fEoverPInvSublead);
  fAnaTtbarDileptonInfo->Branch("fMissHitsSublead", &fMissHitsSublead);
  fAnaTtbarDileptonInfo->Branch("fConversionVetoSublead", &fConversionVetoSublead);
  
  fAnaTtbarDileptonInfo->Branch("fMuonPtLead", &fMuonPtLead);
  fAnaTtbarDileptonInfo->Branch("fMuonEtaLead", &fMuonEtaLead);
  fAnaTtbarDileptonInfo->Branch("fMuonIsoLead", &fMuonIsoLead);

  fAnaTtbarDileptonInfo->Branch("fMuonPtSublead", &fMuonPtSublead);
  fAnaTtbarDileptonInfo->Branch("fMuonEtaSublead", &fMuonEtaSublead);
  fAnaTtbarDileptonInfo->Branch("fMuonIsoSublead", &fMuonIsoSublead);

  
  fAnaTtbarDileptonInfo->Branch("fDeltaPhi_Rec2jets", &fDeltaPhi_Rec2jets);          
  fAnaTtbarDileptonInfo->Branch("fMassDilepton_Rec2jets",&fMassDilepton_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fEtaDilepton_Rec2jets",&fEtaDilepton_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fPtDilepton_Rec2jets",&fPtDilepton_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fPhiDilepton_Rec2jets",&fPhiDilepton_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadJetPt_Rec2jets",&fLeadJetPt_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadJetAbsEta_Rec2jets",&fLeadJetAbsEta_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonPt_Rec2jets",&fLeadRecoLeptonPt_Rec2jets);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonAbsEta_Rec2jets",&fLeadRecoLeptonAbsEta_Rec2jets);     
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonEta_Rec2jets",&fLeadRecoLeptonEta_Rec2jets);
  
  fAnaTtbarDileptonInfo->Branch("fDeltaPhi_Rec2jets_METCut", &fDeltaPhi_Rec2jets_METCut);          
  fAnaTtbarDileptonInfo->Branch("fMassDilepton_Rec2jets_METCut",&fMassDilepton_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fEtaDilepton_Rec2jets_METCut",&fEtaDilepton_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fPtDilepton_Rec2jets_METCut",&fPtDilepton_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fPhiDilepton_Rec2jets_METCut",&fPhiDilepton_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fLeadJetPt_Rec2jets_METCut",&fLeadJetPt_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fLeadJetAbsEta_Rec2jets_METCut",&fLeadJetAbsEta_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonPt_Rec2jets_METCut",&fLeadRecoLeptonPt_Rec2jets_METCut);
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonAbsEta_Rec2jets_METCut",&fLeadRecoLeptonAbsEta_Rec2jets_METCut);     
  fAnaTtbarDileptonInfo->Branch("fLeadRecoLeptonEta_Rec2jets_METCut",&fLeadRecoLeptonEta_Rec2jets_METCut);

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
     if (fRecoLeptonLeadName.Contains("Muon")) leptonLead->SetId(13);
     else leptonLead->SetId(11);
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
       if (fRecoLeptonSubleadName.Contains("Muon")) leptonSublead->SetId(13);
       else leptonSublead->SetId(11);

       TLorentzVector l1 = leptonLead->GetLorentzVector();
       TLorentzVector l2 = leptonSublead->GetLorentzVector();
       TLorentzVector dilepton = l1 + l2;
       
       //dilepton pair should be of opposite sign
       //if(leptonLead->GetCharge()*leptonSublead->GetCharge()==fDileptonChargeSign) {

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
                        0); //dummy PDG
           pPart->SetCharge(leptonLead->GetCharge()*leptonSublead->GetCharge());
           pPart->AddParticle(leptonLead);
           pPart->AddParticle(leptonSublead);
           fDileptonCands++;
	 }
	 //}//sign              
     }//subleading  loop
   }//leading lepton loop
}
