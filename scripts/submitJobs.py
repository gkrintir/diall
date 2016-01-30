#!/usr/bin/env python
import os,sys
import optparse
import commands
import time

import FWCore.ParameterSet.Config as cms
#from storeTools_cff import *
#python scripts/submitJobs.py -q 1nh -j 2 -n 10000 -f 1 -o /afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/first_result_HighPtMuon

# python submitJobs.py -q 1nd -j 50 -n 10000 --proxy proxyforprod

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,dest='queue'  ,help='batch queue'          ,default='1nd')
parser.add_option('-j', '--jobs'       ,dest='jobs'   ,help='number of jobs'       ,default=1,   type=int)
parser.add_option('-f', '--files'      ,dest='files'  ,help='files per job'        ,default=5,   type=int)
parser.add_option('-n', '--nevts'      ,dest='nevts'  ,help='number of events/job' ,default=-1,  type=int)
parser.add_option(      '--proxy'      ,dest='proxy'  ,help='proxy to be used'     ,default=None, type='string')
parser.add_option('-o', '--output'     ,dest='output' ,help='output directory'     ,default='')
(opt, args) = parser.parse_args()

#prepare working directory
print(os.getcwd())
workBase=os.getcwd()
cmsswBase=os.environ['CMSSW_BASE']
outDir=opt.output

jobsBase='%s/FARM%s'%(workBase,time.time())
os.system('mkdir -p %s'%jobsBase)
os.system('cp %s/src/UserCode/diall/test/ExampleAnalysisParameters_cfg.py %s' % (cmsswBase,jobsBase))

#init a new proxy if none has been passed
if opt.proxy is None:
    print 'Initiating a new proxy'
    os.system('voms-proxy-init --voms cms --valid 72:00')
    os.system('cp /tmp/x509up_u`id -u \`whoami\`` %s/proxyforprod' % workBase)
    print 'Production proxy is now available in %s/proxyforprod (valid for 72h)' % workBase

#loop over the required number of jobs
for n in xrange(1,opt.jobs+1):

    #create a wrapper for standalone job
    scriptFile = open('%s/runJob_%d.sh'%(jobsBase,n), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('export X509_USER_PROXY=%s/proxyforprod\n' % workBase)
    scriptFile.write('export OUTDIR=%s\n' % outDir)
    scriptFile.write('cd %s/src\n'%cmsswBase)
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd -\n')
    scriptFile.write('cp %s/ExampleAnalysisParameters_cfg.py .\n' % jobsBase)
    #scriptFile.write('runTtbarEMuData5TeV ExampleAnalysisParameters_cfg.py %d %d\n' % ((n-1)*opt.files,(n-1)*opt.files+opt.files))
    scriptFile.write('runTtbarDilepton5TeV ExampleAnalysisParameters_cfg.py %d %d %d %d\n' % (0, (1-1)*opt.files, (1-1)*opt.files+opt.files, (n-1)*opt.nevts ))
    #scriptFile.write('cmsMkdir $OUTDIR\n')
    scriptFile.write('export OUTPUT=AnaResults_%d.root\n' % n)
    scriptFile.write('cp AnaResults.root $OUTDIR/$OUTPUT\n')
    scriptFile.write('rm AnaResults.root\n')
    scriptFile.close()

    #prepare to run it
    os.system('chmod u+rwx %s/runJob_%d.sh'%(jobsBase,n))

    if opt.queue=='':
        print 'Job #%d will run locally' % n
        #os.system('%s/runJob_%d.sh' % (jobsBase,n) )
    else:
        print 'Job #%d will run remotely' % n
        os.system("bsub -q %s -R \"swp>1000 && pool>30000\" -J tt%d \'%s/runJob_%d.sh\'" % (opt.queue,n,jobsBase,n) )
        
