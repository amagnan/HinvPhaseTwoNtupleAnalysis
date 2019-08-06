#!/bin/sh


if (( "$#" != "1" ))
    then
    echo $# $*
    echo "Input parameter needed: <filedir>"
    echo "for met in 190; do for mht in 100; do for mjj in 2500;do ./prepareDirs.sh METSMEARED/Mjj\${mjj}deta4met\${met}dphi18mht\${mht}; done; done; done"
    exit 1
fi

plotdir=$1
echo "-- Preparing directory structure for $plotdir"
mkdir -p $plotdir/plots0e0mu  $plotdir/plots0e1mu $plotdir/plots0e2mu $plotdir/plots1e0mu $plotdir/plots2e0mu
