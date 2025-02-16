echo 
echo make a test set, e.g.:
echo make_testset.sh suss_6.1.fit 287.80 287.83
echo
stilts in=$1 out=test_cat.fit icmd="select (RA>$2&&RA<$3)"
