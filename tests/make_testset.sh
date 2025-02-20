echo Run this in the test directory 
echo make a test set, e.g.:
echo make_testset.sh suss_6.1.fit uses "287.80 < RA < 287.83"
echo
#alias topcat="java -jar /Users/kuin/bin/topcat-full.jar"  
date
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe   in=$1  ifmt=fits  out=test_cat.fit ofmt=fits  cmd='select (RA>287.80&&RA<287.83)'
date
echo 
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 in1=$1"#2" ifmt1=fits in2=test_cat.fit ifmt2=fits out=test_cat_summary.fits icmd2='keepcols "N_SUMMARY OBSID"' matcher=exact+exact values1='N_SUMMARY OBSID' values2='N_SUMMARY OBSID' find=best1 join=all2  ocmd='delcols "N_SUMMARY_2 OBSID_2 Separation"' ocmd='colmeta -name N_SUMMARY "N_SUMMARY_1";colmeta -name OBSID "OBSID_1"'
date
echo
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$1"#2" ifmt=fits out=suss_summary.fits  
echo 
#java -jar /Users/kuin/bin/topcat-full.jar -stilts tmulti in=test_cat.fits in=test_cat_summary.fits out=test_cat.fit
#rm test_cat.fits 
python add_summary_to_srclist.py
echo
echo Done
