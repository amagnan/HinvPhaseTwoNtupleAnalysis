#!/bin/sh

#

if (( "$#" != "4" ))
then
    echo $# $*
    echo "Input parameter needed: <filedir1> <filedir2> <prod> <puSuffix>"
    echo "./makeFilelists.sh /eos/cms/store/user/amagnan/DelphesNtuples/200PU/ /eos/cms/store/group/phys_higgs/future/amagnan/ 180804 200PU"
    exit
fi

BASEDIR1=$1
BASEDIR2=$2
PROD=$3
PU=$4

FILEPATH1=root://eoscms.cern.ch/$BASEDIR1
FILEPATH2=root://eoscms.cern.ch/$BASEDIR2
OUTPUTDIR=/afs/cern.ch/work/a/amagnan/UPSGAna/HinvPhaseTwoNtupleAnalysis/Analysis/filelists/$PROD/
mkdir -p $OUTPUTDIR

for sample in DYJetsToLL_M-50_HT-70to100 DYJetsToLL_M-50_HT-100to200 DYJetsToLL_M-50_HT-200to400 DYJetsToLL_M-50_HT-400to600 DYJetsToLL_M-50_HT-600to800 DYJetsToLL_M-50_HT-800to1200 DYJetsToLL_M-50_HT-1200to2500 DYJetsToLL_M-50_HT-2500toInf EWKWMinus2Jets EWKWPlus2Jets EWKZ2Jets_ZToLL EWKZ2Jets_ZToNuNu QCD_Mdijet-1000toInf ST_s-channel ST_tch_top ST_tW_antitop ST_tW_top TT VBFH W0JetsToLNu  W1JetsToLNu W2JetsToLNu W3JetsToLNu ZJetsToNuNu_HT-100To200 ZJetsToNuNu_HT-200To400 ZJetsToNuNu_HT-400To600 ZJetsToNuNu_HT-600To800 ZJetsToNuNu_HT-800To1200 ZJetsToNuNu_HT-1200To2500
do
    echo " - Processing $sample"
    rm $OUTPUTDIR/${sample}_${PU}.dat
    BASEDIR=$BASEDIR1
    FILEPATH=$FILEPATH1
    if [ "$sample" == "EWKWMinus2Jets" ] || [ "$sample" == "EWKWPlus2Jets" ] || [ "$sample" == "W2JetsToLNu" ]; then
	BASEDIR=$BASEDIR2
	FILEPATH=$FILEPATH2
    else 
	BASEDIR=$BASEDIR1
	FILEPATH=$FILEPATH1
    fi
    echo $BASEDIR
    #eos ls $BASEDIR/ | grep $sample
    for file in `eos ls $BASEDIR | grep $sample`;
    do
	eos ls $BASEDIR/$file
	if (( "$?" != 0 )); then
	    echo " - File not found !" 
	else
	    echo "$FILEPATH$file" >> $OUTPUTDIR/${sample}_${PU}.dat
	fi
    done
    wc -l $OUTPUTDIR/${sample}_${PU}.dat
    mkdir $OUTPUTDIR/split
    split -l 5 -d $OUTPUTDIR/${sample}_${PU}.dat $OUTPUTDIR/split/${sample}_
    for splitsample in $OUTPUTDIR/split/${sample}_*; 
    do 
	mv $splitsample ${splitsample}_${PU}.dat
    done
    
done
