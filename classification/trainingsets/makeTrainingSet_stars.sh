echo
echo make the training set for stars using the single source catalogue and 
echo criteria on the proper motion and error
echo
echo  - use the existing training set for now
echo 
echo match SUSS to single_source_catalog1 (sources with a PM match in Gaia)
echo
date
#alias topcat="java -jar /Users/kuin/bin/topcat-full.jar"  
topcat -stilts tpipe in=singlesourcecat1.fits infmt=fits out=ssctmp.fits addcol="pm=sqrt(pmra**2+pmdec**2);pm_error=sqrt(pmra_error**+pmdec_error**2)" 
echo
topcat -stilts in=ssctmp.fits infmt=fits out=stars_training_set.fits icmd='select pm/pm_error>30'
# it would be nice to also limit the quality_flag set to 0, but the filters used are unknown
#topcat -stilts in=ssctmp.fits infmt=fits out=stars_training_set.fits icmd='select pm/pm_error>30&&{}_QUALITY_FLAG==0'
echo
echo      STARS training set done      
date