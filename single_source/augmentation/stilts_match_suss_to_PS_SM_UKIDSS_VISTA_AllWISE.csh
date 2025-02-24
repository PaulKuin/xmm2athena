#!/bin/csh
echo 
echo   match SUSS single source catalog "singlesourcecat2.fits" to PS1, SM, UKIDSS, VISTA, WISE
echo
date
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/349/ps1 in=singlesourcecat2.fits out=SUSS_PS.fits find=each ra=ra2000Ep dec=dec2000Ep radius=3 
echo   SUSS_PS.fits Pan-Starrs done
date
echo 
#java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/358/smss in=singlesourcecat2.fits out=SUSS_SMdr1.1.fits find=each  ra=ra2000Ep dec=dec2000Ep radius=3 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/379/smssdr4 in=singlesourcecat2.fits out=SUSS_SMdr4.fits find=each  ra=ra2000Ep dec=dec2000Ep radius=3 
echo   SUSS_SM.fits    Sky Mapper done
date                 
echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/319/gcs9 in=singlesourcecat2.fits out=SUSS_GCS.fits find=each  ra=ra2000Ep dec=dec2000Ep  radius=3 
echo        UKIDSS GCS  done
date                 
echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/319/las9 in=singlesourcecat2.fits out=SUSS_LAS.fits find=each  ra=ra2000Ep dec=dec2000Ep  radius=3 
echo        UKIDSS LAS  done
date                 
echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/319/dxs9 in=singlesourcecat2.fits out=SUSS_DXS.fits find=each  ra=ra2000Ep dec=dec2000Ep  radius=3 
echo        UKIDSS DXS
date                 
echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/316/gps6 in=singlesourcecat2.fits out=SUSS_GCS.fits find=each ra=ra2000Ep dec=dec2000Ep radius=3 
echo        UKIDSS GPS 
date    
echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/343/viking2 in=singlesourcecat2.fits out=SUSS_VIK.fits find=each  ra=ra2000Ep dec=dec2000Ep radius=3 
echo        VISTA VIKING DR2 
date                   
echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/367/vhs_dr5 in=singlesourcecat2.fits out=SUSS_VHS.fits find=each  ra=ra2000Ep dec=dec2000Ep  radius=3 
echo        VISTA  VHS DR5
date                 
 echo 
java -jar /Users/kuin/bin/topcat-full.jar  -stilts cdsskymatch cdstable=II/328/allwise in=singlesourcecat2.fits out=SUSS_WIS.fits find=each  ra=ra2000Ep dec=dec2000Ep radius=3 
echo        AllWISE
date                 
           
	
