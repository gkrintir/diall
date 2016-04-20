#!/bin/bash
source /nfs/home/fynu/gkrintiras/.profile
export X509_USER_PROXY=/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/proxyforprod
export OUTDIR=/home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/test_condor_output/WW/test1
cd /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/
eval `scram r -sh`
cd -
cp /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/WW/ExampleAnalysisParametersSS_cfg.py .

args=("$@")

echo Number of arguments in the batch job: $#
echo 1st argument: ${args[0]}
echo 2nd argument: ${args[1]}

runTtbarDilepton5TeV ${args[0]} ${args[1]} ${args[2]} ${args[3]}
export OUTPUT=AnaResults_${args[4]}.root
cp AnaResults.root $OUTDIR/$OUTPUT
rm AnaResults.root