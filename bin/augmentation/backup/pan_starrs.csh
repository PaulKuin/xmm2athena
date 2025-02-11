date
#topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS1.fits icmd='select (N_SUMMARY>0)&&(N_SUMMARY<=1000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& /dev/null
date
#topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS2.fits icmd='select (N_SUMMARY>1000)&&(N_SUMMARY<=2000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& /dev/null
date
#topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS3.fits icmd='select (N_SUMMARY>2000)&&(N_SUMMARY<=3000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& /dev/null
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS4.fits icmd='select (N_SUMMARY>3000)&&(N_SUMMARY<=4000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS4.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS5.fits icmd='select (N_SUMMARY>4000)&&(N_SUMMARY<=5000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS5.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS6.fits icmd='select (N_SUMMARY>5000)&&(N_SUMMARY<=6000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS6.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS7.fits icmd='select (N_SUMMARY>6000)&&(N_SUMMARY<=7000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS7.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS8.fits icmd='select (N_SUMMARY>7000)&&(N_SUMMARY<=8000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS8.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS9.fits icmd='select (N_SUMMARY>8000)&&(N_SUMMARY<=9000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS9.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS10.fits icmd='select (N_SUMMARY>9000)&&(N_SUMMARY<=10000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS10.log
date
topcat -stilts coneskymatch in=XMM-OM-SUSS5.0.fits.gz ifmt=fits out=SUSS_PS11.fits icmd='select (N_SUMMARY>10000)&&(N_SUMMARY<=11000)&&(DEC>=-30.0)' find=each ra=RA dec=DEC servicetype=cone serviceurl='http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=PS1DR1&' sr=0.000833333 verb=3 parallel=5 >& PS11.log
date
