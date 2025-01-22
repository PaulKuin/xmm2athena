#!/bin/bash
#
IN=$1
echo input file $IN
split -l 500 $IN
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti in=@xaa ifmt=fits multi=false ofmt=fits out=mergeda.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti in=@xab ifmt=fits ofmt=fits out=mergedb.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xac out=mergedc.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xad out=mergedd.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xae out=mergede.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xaf out=mergedf.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xag out=mergedg.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xah out=mergedh.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xai out=mergedi.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xaj out=mergedj.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xak out=mergedk.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xal out=mergedl.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xam out=mergedm.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xan out=mergedn.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xao out=mergedo.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xap out=mergedp.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xaq out=mergedq.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xar out=mergedr.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xas out=mergeds.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xat out=mergedt.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xau out=mergedu.fits
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti ifmt=fits ofmt=fits in=@xav out=mergedv.fits

ls -1 merged?.fits > merged.lis
#
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti multi=true ifmt=fits ofmt=fits in=@merged.lis out=merged_all1.fits

java -jar /Users/kuin/bin/topcat-full.jar -stilts tmultin nin=22 \
ifmt1=fits in1=mergeda.fits \
ifmt2=fits in2=mergedb.fits \
ifmt3=fits in3=mergedc.fits \
ifmt4=fits in4=mergedd.fits \
ifmt5=fits in5=mergede.fits \
ifmt6=fits in6=mergedf.fits \
ifmt7=fits in7=mergedg.fits \
ifmt8=fits in8=mergedh.fits \
ifmt9=fits in9=mergedi.fits \
ifmt10=fits in10=mergedj.fits \
ifmt11=fits in11=mergedk.fits \
ifmt12=fits in12=mergedl.fits \
ifmt13=fits in13=mergedm.fits \
ifmt14=fits in14=mergedn.fits \
ifmt15=fits in15=mergedo.fits \
ifmt16=fits in16=mergedp.fits \
ifmt17=fits in17=mergedq.fits \
ifmt18=fits in18=mergedr.fits \
ifmt19=fits in19=mergeds.fits \
ifmt20=fits in20=mergedt.fits \
ifmt21=fits in21=mergedu.fits \
ifmt22=fits in22=mergedv.fits \
out=merged_all2.fits ofmt=fits
