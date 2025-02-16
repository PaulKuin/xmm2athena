echo Run this in the test directory 
echo make a test set, e.g.:
echo make_testset.sh suss_6.1.fit uses "287.80 < RA < 287.83"
echo
#alias topcat="java -jar /Users/kuin/bin/topcat-full.jar"  
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe   in=$1    ifmt=fits  out=test_cat.fit ofmt=fits  cmd='select (RA>287.80&&RA<287.83)'
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 in1=$1"#2" ifmt1=fits in2=test_cat.fit ifmt2=fits out=test_cat_summary.fits matcher=exact+exact values1='N_SUMMARY OBSID' values2='N_SUMMARY OBSID' find=best1 join=all2 fixcols=none
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$1"#2" ifmt=fits out=suss_summary.fit  
echo 
python add_summary_to_srclist.py
echo
echo Done
