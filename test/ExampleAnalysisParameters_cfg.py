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

tt_PbPb=("ZMMpp_v0.root","/store/group/phys_heavyions/velicanu/forest/Run2015E/SingleMuHighPt/RAW-RECO/ZMM-PromptReco-v1/Merged/")
DY_PbPb=("DY_PbPb.root","/store/cmst3/group/hintt/CMSSW/DY/NTUPLE")
QCDEM_PbPb=("QCDEM_PbPb.root","/store/cmst3/group/hintt/CMSSW/QCDEMriched_80to250/NTUPLE")
QCDDJ_PbPb=("QCDDJ_PbPb.root","/store/cmst3/user/mverweij/jetsPbPb/Run2Prep/Dijet80CMSSW753p1/v6/PyquenUnquenched_Dijet_NcollFilt_pthat80_740pre8_MCHI1_74_V4_GEN-SIM_v3/crab_HiForestDijet80Run2Fullv6/151020_135458")
HCALMeth2=("HCALM2.root","/store/user/istaslis/PyquenUnquenched_Dijet_NcollFilt_pthat80_753patch1_RECODEBUG_FOREST/")

centralityRequirements={"inc":[0,200],
                        "0to20":[0,40],
                        "20to50":[40,100],
                        "50to80":[100,160],
                        "80to100":[160,200]}


sample=tt_PbPb
#sample=DY_PbPb
#sample=QCDEM_PbPb
#sample=HCALMeth2 
#sample = QCDDJ_PbPb

centralityBins=centralityRequirements["inc"]
#centralityBins=centralityRequirements["0to20"]
#centralityBins=centralityRequirements["20to50"]
#centralityBins=centralityRequirements["50to80"]
#centralityBins=centralityRequirements["80to100"]

config = cms.PSet(
    output = cms.string('cen_%dto%d_%s'%(centralityBins[0],centralityBins[1],sample[0])),
    input  = cms.vstring( fillFromStore(sample[1]) ),
    maxEvents = cms.int32(-1),#-1),
    minCentrality = cms.int32(centralityBins[0]),
    maxCentrality = cms.int32(centralityBins[1])#,
   # anaFile = cms.int32(options.anaFile)
)
