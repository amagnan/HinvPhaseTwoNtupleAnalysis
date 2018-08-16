#!/bin/sh


if (( "$#" != "3" ))
    then
    echo $# $*
    echo "Input parameter needed: <filedir> <prod> <log>"
    echo "./makeParamFile.sh /afs/cern.ch/work/a/amagnan/public/UPSGAna/ 180716 runTree.log"
    exit
fi

BASEDIR=$1
PROD=$2
LOGFILE=$3

INPUTDIR=$1/$2

rm ParamFile_${PROD}.dat


for sample in DYJetsToLL_M-50_HT-70to100 DYJetsToLL_M-50_HT-100to200 DYJetsToLL_M-50_HT-200to400 DYJetsToLL_M-50_HT-400to600 DYJetsToLL_M-50_HT-600to800 DYJetsToLL_M-50_HT-800to1200 DYJetsToLL_M-50_HT-1200to2500 DYJetsToLL_M-50_HT-2500toInf EWKWMinus2Jets EWKWPlus2Jets EWKZ2Jets_ZToLL EWKZ2Jets_ZToNuNu QCD_Mdijet-1000toInf ST_s-channel ST_tch_top ST_tW_antitop ST_tW_top TT VBFH W0JetsToLNu  W1JetsToLNu W2JetsToLNu W3JetsToLNu ZJetsToNuNu_HT-100To200 ZJetsToNuNu_HT-200To400 ZJetsToNuNu_HT-400To600 ZJetsToNuNu_HT-600To800 ZJetsToNuNu_HT-800To1200 ZJetsToNuNu_HT-1200To2500
#for sample in `ls $INPUTDIR`
do
    Nevts=0
    echo " - Processing $sample"
    for subsample in `ls -d $INPUTDIR/${sample}_*_200PU`
    do
	ls $subsample/runTree.log
	if (( "$?" == 0 )); then
	    tmpevt=`grep "\-\- Processing" $subsample/$LOGFILE | awk '{print $3}'`
	    echo "----"$subsample" "$tmpevt
	    let Nevts=$Nevts+tmpevt
	fi
    done
    echo $sample" "$Nevts" " >> ParamFile_${PROD}.dat
done
