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

#include <TH1D.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>

//#include "CollectFiles.C"

#include <iostream>

using namespace std;

void countPrescalesForTtbarDilepton5TeV(std::vector<std::string> urls, const char *outname = "AnaResults.root", Long64_t nentries = 20, Int_t firstF = -1, Int_t lastF = -1, Int_t firstEvent = 0, int isData = 1) {


  
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


  //---------------------------------------------------------------
  //Event loop
  //---------------------------------------------------------------	  
  TFile *out = new TFile("test_runnumber.root","RECREATE");
  Long64_t entries_tot =  chain->GetEntries(); 
  Int_t run, lum;// presc_val_1, presc_val_2, presc_val_3, presc_val_4, presc_val_5, presc_val_6;
  Int_t presc_val[6];
  std::vector<std::string> presc_name = {"HLT_HISinglePhoton10_Eta1p5_v1_Prescl", "HLT_HISinglePhoton10_Eta3p1_v1_Prescl"
					 ,"HLT_HISinglePhoton15_Eta1p5_v1_Prescl", "HLT_HISinglePhoton15_Eta3p1_v1_Prescl"
					 ,"HLT_HISinglePhoton20_Eta1p5_v1_Prescl", "HLT_HISinglePhoton20_Eta3p1_v1_Prescl"};
  std::set<Int_t> runNumber;
  std::set<Int_t>::iterator it;
  std::pair<std::set<Int_t>::iterator,Bool_t> check;
  TList * tList;
  tList = new TList(); tList->SetOwner(); tList->SetName("Presc_vs_Lumi");
  std::cout << entries_tot << std::endl;
  chain->SetBranchStatus("Run", 1);
  chain->SetBranchStatus("LumiBlock", 1);

  for(unsigned int i=0; i<presc_name.size(); i++)
  {
    chain->SetBranchStatus(presc_name[i].data(), 1);
    chain->SetBranchAddress(presc_name[i].data(), &presc_val[i]);
  }

  chain->SetBranchAddress("Run", &run);
  chain->SetBranchAddress("LumiBlock", &lum);

  //chain->SetBranchAddress("HLT_HISinglePhoton20_Eta1p5_v1_Prescl", &presc);
  //chain->SetBranchAddress("HLT_HISinglePhoton20_Eta3p1_v1_Prescl", &presc);
  if(nentries<0) lastEvent = chain->GetEntries();  

  // Long64_t nentries = 20;//chain->GetEntriesFast();
  Printf("nentries: %lld  tot: %lld", nentries,entries_tot);
  std::cout<< runNumber.size()<<std::endl;

  for (Long64_t jentry=firstEvent; jentry<lastEvent; ++jentry) 
  { 
      //Run producers
      if (jentry>=entries_tot) break;
      if (jentry%1000==0) Printf("Processing event %d  %d",(int)(jentry), (int)(lastEvent));
      chain->GetEntry(jentry);

      //std::cout<< runNumber.size()<<std::endl;
      check = runNumber.insert(run);
      //std::cout<< runNumber.size()<<std::endl;

      //std::cout<< check.second <<std::endl;
      std::string name(std::to_string(run));

      if (check.second==kTRUE) 
      {
	  TDirectory * first_level = out->mkdir(name.data());
	  first_level->cd();
	  for(unsigned int i=0; i<presc_name.size(); i++)
	  {
	      //TH2D * hist = new TH2D(Form("%s_%s", name.data(), presc_name[i].data()), Form("%s_%s", name.data(), presc_name[i].data()), 10001, -0.5, 10000.5, 10001, -0.5, 10000.5);
	    TH1D * hist = new TH1D(Form("%s_%s", name.data(), presc_name[i].data()), Form("%s_%s", name.data(), presc_name[i].data()), 10001, -0.5, 10000.5);//, 10001, -0.5, 10000.5);
	      tList->Add(hist);
	      hist->SetBinContent(lum, presc_val[i]);
	  }
      }
      else
      {
	
	//it=check.first;
	for(unsigned int i=0; i<presc_name.size(); i++)
	{
	    TH1D * hist = dynamic_cast<TH1D*>(tList->FindObject(Form("%s_%s", name.data(), presc_name[i].data())));
	    hist->SetBinContent(lum, presc_val[i]);
	}
      }
      
  }

  tList->ls();

  out->Write();
  out->Close();
}




