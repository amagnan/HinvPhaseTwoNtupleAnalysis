#!/bin/sh

#

if (( "$#" != "3" ))
    then
    echo $# $*
    echo "Input parameter needed: <filedir> <prod> <puSuffix>"
    echo "./mergeOutputs.sh /afs/cern.ch/work/a/amagnan/public/UPSGAna/ 180716 200PU"
     exit
fi

BASEDIR=$1
PROD=$2
PU=$3
EOSDIR=/eos/cms/store/user/amagnan/DelphesNtuples/Trees/$PROD

INPUTDIR=$BASEDIR$PROD

for sample in DYJetsToLL_M-50_HT-70to100 DYJetsToLL_M-50_HT-100to200 DYJetsToLL_M-50_HT-200to400 DYJetsToLL_M-50_HT-400to600 DYJetsToLL_M-50_HT-600to800 DYJetsToLL_M-50_HT-800to1200 DYJetsToLL_M-50_HT-1200to2500 DYJetsToLL_M-50_HT-2500toInf EWKWMinus2Jets EWKWPlus2Jets EWKZ2Jets_ZToLL EWKZ2Jets_ZToNuNu QCD_Mdijet-1000toInf ST_s-channel ST_tch_top ST_tW_antitop ST_tW_top TT VBFH W0JetsToLNu  W1JetsToLNu W2JetsToLNu W3JetsToLNu ZJetsToNuNu_HT-100To200 ZJetsToNuNu_HT-200To400 ZJetsToNuNu_HT-400To600 ZJetsToNuNu_HT-600To800 ZJetsToNuNu_HT-800To1200 ZJetsToNuNu_HT-1200To2500
do
    echo " - Processing $sample"
    ls $INPUTDIR/${sample}_${PU}.root
    if (( "$?" == 0 )); then
	echo " - File exists already, checking more files to hadd."
	ls $INPUTDIR/${sample}*/Hist*.root
	if (( "$?" == 0 )); then
	    echo " - found more to hadd."
	    mv $INPUTDIR/${sample}_${PU}.root /tmp/amagnan/${sample}_${PU}_prev.root
	    hadd $INPUTDIR/${sample}_${PU}.root /tmp/amagnan/${sample}_${PU}_prev.root $INPUTDIR/${sample}*/Hist*.root
	else
	    echo " - no more to hadd, just clean and cp to eos."
	fi
    else
	hadd $INPUTDIR/${sample}_${PU}.root $INPUTDIR/${sample}*/Hist*.root
    fi
    
    if (( "$?" == 0 )); then
	echo " - success: deleting individual files."
	rm $INPUTDIR/${sample}*/Hist*.root
	rm /tmp/amagnan/${sample}_${PU}_prev.root
	eos mkdir -p $EOSDIR
	eos cp $INPUTDIR/${sample}_${PU}.root $EOSDIR/${sample}_${PU}.root
    else
	echo " -hadd failed: please check."
    fi
    
done
