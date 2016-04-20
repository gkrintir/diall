#!/bin/bash
source /nfs/home/fynu/gkrintiras/.profile
export X509_USER_PROXY=/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/proxyforprod
export OUTDIR=/home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/test_condor_output
cd /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/
eval `scram r -sh`
cd -
cp /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/ExampleAnalysisParameters_cfg.py .
runTtbarDilepton5TeV $@ $@@
export OUTPUT=AnaResults_1.root
cp AnaResults.root $OUTDIR/$OUTPUT
rm AnaResults.root
