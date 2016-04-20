#include "UserCode/diall/analyzeTtbarDilepton5TeV.C"

//#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TString.h"
#include <iostream>

using namespace std;

//
// MAIN METHOD
//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  //  AutoLibraryLoader::enable();
  FWLiteEnabler::enable();

  //check arguments
  if ( argc < 2 ) 
    {
      std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
      return 0;
    }

  // int anaFile = 1;
  int firstFile = 0;
  int lastFile  = 1;
  Int_t firstEvent = 0;
  int isData = 1;
  
  std::cout << "Have " << argc << " arguments:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout << argv[i] << std::endl;
    switch(argc)
      {
      
      case 3: 
	isData = atoi(argv[2]);
	break;
      case 5:
        isData = atoi(argv[2]);
        firstFile = atoi(argv[3]);
        lastFile = atoi(argv[4]);
        break;
      case 6:
	isData = atoi(argv[2]);
	firstFile = atoi(argv[3]);
	lastFile = atoi(argv[4]);
	firstEvent = atoi(argv[5]);
	break;
      }
  }
  // read configuration
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("config");
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");

  for (std::vector<std::string>::const_iterator i = urls.begin(); i != urls.end(); ++i)
    std::cout << *i << std::endl;

  //Output
  std::string outname = "AnaResults.root";
  // std::string outname = runProcess.getParameter<std::string>("output");
  //Events
  int maxEvts = runProcess.getParameter<int>("maxEvents");
  //Trigger Object
  TString triggerPath = runProcess.getParameter<string>("triggerPath");
  //RECO Muon Object
  float muonPtCut = runProcess.getParameter<double>("muonPtCut");
  float muonEtaCut = runProcess.getParameter<double>("muonEtaCut");
  //RECO Electron Object
  float electronPtCut = runProcess.getParameter<double>("electronPtCut");
  float electronEtaCut = runProcess.getParameter<double>("electronEtaCut");
  int electronVetoID = runProcess.getParameter<int>("electronVetoID");
  int electronLooseID = runProcess.getParameter<int>("electronLooseID");
  int electronMediumID = runProcess.getParameter<int>("electronMediumID");
  int electronTightID = runProcess.getParameter<int>("electronTightID");
  //RECO Jet Object
  float jetPtCut = runProcess.getParameter<double>("jetPtCut");
  float jetEtaCut = runProcess.getParameter<double>("jetEtaCut");

  //Analysis Objects
  float dilepMassCut = runProcess.getParameter<double>("dilepMassCut");
  int dilepSign = runProcess.getParameter<int>("dilepSign");
  float metCut = runProcess.getParameter<double>("metCut");

  std::cout<< maxEvts<< firstFile << lastFile << firstEvent<< isData<< std::endl;//firstEvent << isData << std::endl;

  analyzeTtbarDilepton5TeV(urls,outname.c_str(),maxEvts,firstFile,lastFile,firstEvent,isData,
			   triggerPath,
			   muonPtCut, muonEtaCut,
			   electronPtCut, electronEtaCut, electronVetoID, electronLooseID, electronMediumID, electronTightID, 
			   jetPtCut, jetEtaCut, 
			   dilepMassCut, dilepSign, metCut);
  
  cout << "Results have been stored in " << outname << endl;

}
