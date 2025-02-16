echo 
echo make a test set, e.g.:
echo make_testset.sh suss_6.1.fit 287.80 287.83
echo
#alias topcat="java -jar /Users/kuin/bin/topcat-full.jar"  
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$1 ifmt=fits out=test_cat.fit ofmt=fits cmd='select (RA>287.80&&RA<287.83)'
echo Done
