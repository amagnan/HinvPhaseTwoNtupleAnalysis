#!/usr/bin/env python
import subprocess
import os,sys
import optparse
import commands
import math
import random

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'  ,    dest='queue'             , help='batch queue'             , default='2nd')
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-n', '--numEvts' ,  dest='numEvts'  , help='Number of events to process' , default=0)

parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

pulist = ['200PU']

#inputlistdir = 'filelists/180807/split/'
inputlistdir = 'filelists/180809/'

runDir=os.getcwd()

for pu in pulist :
    os.system('echo " - pu %s"'%pu)

    for inputlist in os.listdir(inputlistdir):
        lastchar=inputlist.rfind("_%s.dat"%(pu))
        print "%s %d"%(inputlist,lastchar)
        sample=inputlist[0:lastchar]
        os.system('echo " -- sample: %s"'%(sample))
        outDir='%s/%s_%s'%(opt.out,sample,pu)
        os.system('mkdir -p %s'%outDir)

    #wrapper
        scriptFile = open('%s/runTreeJob.sh'%(outDir), 'w')
        scriptFile.write('#!/bin/bash\n')
        scriptFile.write('export BASEDIR=$PWD \n')
        scriptFile.write('export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n')
        scriptFile.write('source $VO_CMS_SW_DIR/cmsset_default.sh \n')
        scriptFile.write('cd /afs/cern.ch/work/a/amagnan/CMSSW_9_3_5/src/PhaseTwoAnalysis/delphesInterface/ \n')
        scriptFile.write('eval `scramv1 runtime -sh` \n')
        scriptFile.write('source env.sh\n')
        scriptFile.write('cd $BASEDIR\n')
        scriptFile.write('pwd\n')
        scriptFile.write('echo "running with input filelist: %s/%s/%s"\n'%(os.getcwd(),inputlistdir,inputlist))
        scriptFile.write('echo $ROOTSYS \n')
        #run light tree maker
        scriptFile.write('/afs/cern.ch/work/a/amagnan/UPSGAna/HinvPhaseTwoNtupleAnalysis/Analysis/bin/simpleTree ./ %s %s %s/%s %d | tee runTree.log\n'%(sample,pu,os.getcwd(),inputlistdir,opt.numEvts))
        scriptFile.write('cp runTree.log %s/runTree.log\n'%(outDir))
        scriptFile.write('ls -lh * \n')
        scriptFile.write('mv HistosFile* %s/\n'%(outDir))
        scriptFile.write('echo " -- ntuple to lighttree done"\n')
        scriptFile.write('echo " -- All done"\n')
        scriptFile.close()
        
        #get dir from filelist
        #tmpDir=subprocess.check_output(['head','-n 1','%s/%s.dat'%(inputlistdir,sample)])
        #lastchar=tmpDir.rfind("/")
        #sampleDir=tmpDir[0:lastchar]
        #print "Sample dir is: %s"%sampleDir

        os.system('chmod u+rwx %s/runTreeJob.sh'%outDir)

        condorFile = open('%s/condorSubmitTree.sub'%(outDir), 'w')
        condorFile.write('universe = vanilla\n')
        condorFile.write('+JobFlavour = "nextweek"\n')
        condorFile.write('Executable = %s/runTreeJob.sh\n'%outDir)
        condorFile.write('Output = %s/condorTree.out\n'%outDir)
        condorFile.write('Error = %s/condorTree.err\n'%outDir)
        condorFile.write('Log = %s/condorTree.log\n'%outDir)
        condorFile.write('Queue 1\n')
        condorFile.close()

    #submit
        if opt.nosubmit :
            os.system('echo condor_submit %s/condorSubmitTree.sub'%(outDir))
            #os.system('echo bsub -q %s %s/runTreeJob.sh'%(opt.queue,outDir))
            #os.system('cat %s/runTreeJob.sh'%(outDir))
            #os.system('cat %s/condorSubmitTree.sub'%(outDir))
        else: 
            #os.system("bsub -q %s \'%s/runTreeJob.sh\'"%(opt.queue,outDir))
            os.system('condor_submit %s/condorSubmitTree.sub'%(outDir))
