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

jobsBase='%s/%s'%(workBase,opt.output)
#os.system('mkdir -p %s'%jobsBase)
#os.system('cp %s/src/UserCode/diall/test/ExampleAnalysisParameters_cfg.py %s' % (cmsswBase,jobsBase))

#create a wrapper for standalone job
scriptFile = open('%s/diall.condor'%(jobsBase), 'w')
scriptFile.write('executable = runJob_1.sh \n')
scriptFile.write('universe = vanilla \n')
scriptFile.write('requirements   = (CMSFARM =?= TRUE) \n')
scriptFile.write('output = log/log$(Process).out \n')
scriptFile.write('error = log/log$(Process).error  \n')
scriptFile.write('log = log/log$(Process).log \n')
scriptFile.write('should_transfer_files = YES \n')
scriptFile.write('when_to_transfer_output = ON_EXIT \n')

#loop over the required number of jobs
for n in xrange(1,opt.jobs+1):

    scriptFile.write('arguments = /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/tW_antitop/ExampleAnalysisParameters_cfg.py %d %d %d %d\n' % ( 1, (n-1)*(opt.files-1), n*(opt.files-1), n ))
    #scriptFile.write('arguments = /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/TTbar/ExampleAnalysisParameters_cfg.py %d %d %d %d %d\n' % ( 0, 0, 1, (n-1)*opt.nevts, n ))
    scriptFile.write('queue \n')

scriptFile.close()
