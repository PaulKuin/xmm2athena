echo 
echo match SUSS to SDSS to get spectroscopically identified QSOs and galaxies for the training sets
echo
echo - use the existing training sets for now 
date
topcat -stilts cdsskymatch cdstable=V/154/sdss16 in=singlesourcecat2.fits out=sdss_matches.fits find=best ra=RAJ2000Ep2000 dec=DEJ2000Ep2000 radius=3 
echo      QSO training set       
topcat -stilts in=sdss_matches.fits ifmt=fits out=qso_training_set.fits icmd='select (spcl=="QSO")&&(Q==3)'  ocmd=delcols "GroupID GroupSize"
echo   galaxy training set
topcat -stilts in=sdss_matches.fits ifmt=fits out=galaxy_training_set.fits icmd='select (spcl=="GALAXY")&&(Q==3)'  ocmd=delcols "GroupID GroupSize"
date                 
