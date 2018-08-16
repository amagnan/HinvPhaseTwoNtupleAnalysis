#!/bin/sh


if (( "$#" != "1" ))
    then
    echo $# $*
    echo "Input parameter needed: <filedir>"
    echo "./prepareDirs.sh Mjj800deta2met100"
    exit
fi

plotdir=$1
mkdir -p $plotdir/plots0e0mu  $plotdir/plots0e1mu $plotdir/plots0e2mu $plotdir/plots1e0mu $plotdir/plots2e0mu
