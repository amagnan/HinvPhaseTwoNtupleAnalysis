#!/bin/sh


if (( "$#" != "3" ))
    then
    echo $# $*
    echo "Input parameter needed: <mjj> <met> <mht>"
    echo "./correctGluonYield.sh 1000 130 100"
    exit
fi

mjj=$1
met=$2
mht=$3

root -b -q 
result=`root -l -b -q "getggRatio.C++($mjj,$met,$mht)" | grep "Check" | awk '{print $5}'`

echo "Script output: $met $mht $mjj $result"


#tmpstr1=`head -n 1 Mjj${mjj}deta4met${met}dphi18mht100/Datacard_R0e0mu.txt`
#tmpstr2=`head -n 1 Mjj${mjj}deta4met${met}dphi18mht100/Datacard_R0e0mu.txt | awk -v result="$result" '{$3=$2*result}{print $0}'`
tmpstr1=`head -n 1 Mjj${mjj}deta4met${met}dphi18mht${mht}/Datacard_R0e0mu.txt`
tmpstr2=`head -n 1 Mjj${mjj}deta4met${met}dphi18mht${mht}/Datacard_R0e0mu.txt | awk -v result="$result" '{$3=$2*result}{print $0}'`

echo $tmpstr1
echo $tmpstr2

#sed "0,/$tmpstr1/{s/$tmpstr1/$tmpstr2/}" Mjj${mjj}deta4met${met}dphi18mht100/Datacard_R0e0mu.txt > tmp 
#mv tmp Mjj${mjj}deta4met${met}dphi18mht100/Datacard_R0e0mu.txt
sed "0,/$tmpstr1/{s/$tmpstr1/$tmpstr2/}" Mjj${mjj}deta4met${met}dphi18mht${mht}/Datacard_R0e0mu.txt > tmp 
mv tmp Mjj${mjj}deta4met${met}dphi18mht${mht}/Datacard_R0e0mu.txt
