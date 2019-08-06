#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)

parser.add_option('-q', '--queue' ,  dest='queue'  , help='batch queue' , default='1nd')
parser.add_option('-d', '--outdir' ,  dest='outdir'  , help='output directory' , default='LooseVBFsel')
parser.add_option('-p', '--process' ,  dest='process'  , help='process short name' , default='')
parser.add_option('-l', '--filelist' ,  dest='filelist'  , help='directory path with filelists' , default='filelists_20171221/')
parser.add_option('-n', '--numEvts' ,  dest='numEvts'  , help='Number of events to process' , default=0)

parser.add_option('-S', '--no-submit',action="store_true",  dest='nosubmit', help='Do not submit batch job.')
(opt, args) = parser.parse_args()


#workdir='/afs/cern.ch/work/a/amagnan/UPSGAna/'



#pulist = ['noPU','PU200']
pulist = ['200PU']

for pu in pulist :

    #processList=['DYToLL-M-50_0J','DYToLL-M-50_1J','DYToLL-M-50_2J','DYToLL-M-50_3J','EWKWMinus2Jets_WToLNu_M-50','EWKWPlus2Jets_WToLNu_M-50','EWKZ2Jets_ZToLL_M-50','EWKZ2Jets_ZToNuNu','QCD_Mdijet-1000toInf','ST_tW_DR_14TeV_top','ST_tW_DR_14TeV_antitop','ST_tch_14TeV_antitop','ST_tch_14TeV_top','TT_TuneCUETP8M2T4','VBFH','WToLNu_0J','WToLNu_1J','WToLNu_2J','WToLNu_3J','ZJetsToNuNu_HT-100To200','ZJetsToNuNu_HT-1200To2500','ZJetsToNuNu_HT-200To400','ZJetsToNuNu_HT-400To600','ZJetsToNuNu_HT-600To800','ZJetsToNuNu_HT-800To1200']
    processList=['VBFH']


    if len(opt.process)>0:
        processList=[opt.process]


    for myproc in processList :
        outDir='%s/%s/%s/'%(opt.outdir,pu,myproc)
        os.system('mkdir -p %s'%outDir)

        scriptFile = open('%s/runJob.sh'%(outDir), 'w')
        scriptFile.write('#!/bin/bash\n')
        scriptFile.write('cd /afs/cern.ch/work/a/amagnan/CMSSW_9_3_5/src/\n')
        scriptFile.write('cmsenv\n')
        scriptFile.write('cd -\n')
        scriptFile.write('mkdir -p %s\n'%(opt.filelist))
        scriptFile.write('cp %s/%s/* %s/\n'%(os.getcwd(),opt.filelist,opt.filelist))
        scriptFile.write('ls *\n')
    #os.system('./bin/simpleTree %s %s'%(outDir,myproc))
        scriptFile.write('%s/bin/simpleTree %s %s %s %s %d 1 | tee %s/runJob.log\n'%(os.getcwd(),outDir,myproc,pu,opt.filelist,opt.numEvts,outDir))
        scriptFile.write('echo "All done"\n')
        scriptFile.close()
 
  #submit                                                                                                                                           
        os.system('chmod u+rwx %s/runJob.sh'%outDir)
        if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(opt.queue,outDir))
        else: os.system("bsub -q %s \'%s/runJob.sh\'"%(opt.queue,outDir))
