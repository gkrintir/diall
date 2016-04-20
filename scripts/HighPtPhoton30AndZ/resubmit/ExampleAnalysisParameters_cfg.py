import FWCore.ParameterSet.Config as cms
from UserCode.diall.storeTools_cff import *
import FWCore.ParameterSet.VarParsing as VarParsing

## parse some command line arguments
#options = VarParsing.VarParsing ('analysis')

#options.register ('anaFile',
#                  1, # default value
#                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
#                  VarParsing.VarParsing.varType.int,       # string, int, or float
#                  "Number of files to process (-1 for all)")

#options.anaFile = 1

#options.parseArguments()

#print('anaFile: ',options.anaFile)

tt_PbPb=("HiForest_0.root","/store/cmst3/group/hintt/CMSSW_7_5_8_patch2/TTbar/Forest/merge/")
DY_PbPb=("HiForest_0.root","/store/cmst3/group/hintt/CMSSW_7_5_8_patch2/DY/Forest/merge/")
W_PbPb =("HiForest_0.root","/store/cmst3/group/hintt/CMSSW_7_5_8_patch2/W/Forest/merge/")
WW_PbPb=("HiForest_0.root","/store/cmst3/group/hintt/CMSSW_7_5_8_patch2/WW/Forest/merge/")

centralityRequirements={"inc":[0,200],
                        "0to20":[0,40],
                        "20to50":[40,100],
                        "50to80":[100,160],
                        "80to100":[160,200]}


#sample=tt_PbPb
#sample=DY_PbPb
#sample=W_PbPb
sample=WW_PbPb 


FSimfile = open('/home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/HighPtPhoton30AndZ_Run2015E-PromptReco-v1_v1.data', 'r')
FSimlist = cms.vstring()
FSimlist.extend( [line.strip() for line in FSimfile.read().splitlines()] )


centralityBins=centralityRequirements["inc"]
#centralityBins=centralityRequirements["0to20"]
#centralityBins=centralityRequirements["20to50"]
#centralityBins=centralityRequirements["50to80"]
#centralityBins=centralityRequirements["80to100"]

config = cms.PSet(
    output = cms.string('cen_%dto%d_%s'%(centralityBins[0],centralityBins[1],sample[0])),
    input  = FSimlist,
        #cms.vstring( 
        #fillFromStore(sample[1])
        #'root://eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch2/TTbar/Forest/HiForest_999.root'
        #'root://cmsxrootd.fnal.gov//store/user/velicanu/Merged/Pythia8_Z30mumuJet_pthat30Norm_TuneCUETP8M1_5020GeV_cff_ppFOREST_PrivMC_v10/0.root'
        #'/storage/data/cms/store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_merged.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_2.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_3.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_4.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_5.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_6.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_7.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_8.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_9.root',
        #'root://cmsxrootd.fnal.gov//store/user/gkrintir/TopHI/HighPtLowerPhotons_Run2015E-PromptReco-v1_v1/HighPtLowerPhotons/crab_TopHI/160111_200108/0000/HiForestAOD_10.root'
        #), #fillFromStore(sample[1]) 
    maxEvents = cms.int32(-1),#-1),
    minCentrality = cms.int32(centralityBins[0]),
    maxCentrality = cms.int32(centralityBins[1]),
    triggerPath = cms.string('HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1'),
    muonPtCut = cms.double(18.),
    muonEtaCut = cms.double(2.4),
    electronPtCut = cms.double(30.),
    electronEtaCut = cms.double(2.4),
    electronVetoID = cms.int32(0),
    electronLooseID = cms.int32(1), 
    electronMediumID = cms.int32(0),
    electronTightID = cms.int32(0),
    jetPtCut = cms.double(25.),
    jetEtaCut = cms.double(3.),
    dilepMassCut = cms.double(20.),
    dilepSign = cms.int32(-1), 
    metCut = cms.double(30.)
   # anaFile = cms.int32(options.anaFile)
)
