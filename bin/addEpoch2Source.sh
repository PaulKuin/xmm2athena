#!/bin/bash
# 
# How to use:
#  addEpoch2SourceSUSS.sh XMM-OM-SUSS5.0.fits <output_file_name>
# where the input  is the SUSS or UVOTSSC file with two extensions #1 and #2 for the
#   SRCLIST and SUMMARY 
#
# The OUTPUT is just the SRCLIST with added column "obs_epoch"
# 
# command line test:
#java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe  in=XMM-OM-SUSS5.0.fits#2 ifmt=fits out=tmp.fits ofmt=fits-basic cmd="addcol obs_epoch  mjdToDecYear(0.5*(MJD_START+MJD_END));keepcols 'OBSID obs_epoch'"
#
echo "addEpoch2Source.sh"
#
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe  \
in=$1\#2 ifmt=fits out=tmp2.fits ofmt=fits-basic \
cmd="addcol -units "y" obs_epoch  mjdToDecYear(0.5*(MJD_START+MJD_END));keepcols 'OBSID obs_epoch'"
#
if [ "$2" != "" ]; then 
   OUT=$2
else 
   OUT="withObsEpoch_$1"
fi
# 
# command line test:
#java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 in1=$IN\#1 ifmt1=fits in2=tmp.fits ifmt2=fits out=$IN_with_epoch.fits ofmt=fits-basic matcher=exact values1="OBSID" values2="OBSID" join=all1 find=best1  fixcols=dups suffix1=  suffix2=_2
#
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 \
in1=$1\#1 ifmt1=fits in2=tmp2.fits ifmt2=fits \
out=tmp1.fits ofmt=fits-basic \
matcher=exact values1="OBSID" values2="OBSID" \
join=all1 find=best1  fixcols=dups suffix1=  suffix2=_2 
#
# edit the extra OBSID out
#
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe  \
in=tmp1.fits ifmt=fits out=$OUT ofmt=fits-basic \
cmd="delcols 'OBSID_2 GroupID GroupSize'"
#
# recombining with the SUMMARY extension from the SUSS has not been done 
# e.g., fappend $1#2 $OUT
# clean up
#
rm tmp1.fits tmp2.fits
