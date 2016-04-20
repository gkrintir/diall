#!/bin/bash

source /nfs/home/fynu/gkrintiras/.profile
export X509_USER_PROXY=/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/proxyforprod
export HOME=/nfs/home/fynu/gkrintiras/
export OUTDIR=/nfs/scratch/fynu/gkrintiras/merge

cp /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/Muons/runMergingForest.C .
cp /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/Muons/mergeFileMerger.C .

args=("$@")

echo Number of arguments in the batch job: $#
echo 1st argument: ${args[0]}
echo 2nd argument: ${args[1]}
echo 3nd argument: ${args[2]}
echo 4nd argument: ${args[3]}
echo 5nd argument: ${args[4]}

root -b -q 'runMergingForest.C("'${args[0]}'", '${args[2]}', '${args[3]}', "'${args[4]}'")'
export OUTPUT=${args[4]}

gridhost=130.104.133.253:2811
gridpath=/pnfs/cism.ucl.ac.be/data/cms/sca06/users/gkrintir/TopHI/SingleMuHighPt_Run2015E-PromptReco-v1_v1/SingleMuHighPt/crab_TopHI/160126_203339/0000/First_Merge

cp ${args[4]} $OUTDIR/$OUTPUT

rm HiForest_*.root


