#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out2.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>1000)&&(N_SUMMARY<=2000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out3.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>2000)&&(N_SUMMARY<=3000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out4.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>3000)&&(N_SUMMARY<=4000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out5.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>4000)&&(N_SUMMARY<=5000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out6.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>5000)&&(N_SUMMARY<=6000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out7.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>6000)&&(N_SUMMARY<=7000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out8.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>7000)&&(N_SUMMARY<=8000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
#date
#topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out9.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>8000)&&(N_SUMMARY<=9000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
#keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
#colmeta -name ra_x_2 "ra_x";\
#colmeta -name dec_x_2 "dec_x";\
#'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out10.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>9000)&&(N_SUMMARY<=10000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out11.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>10000)&&(N_SUMMARY<=11000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts tcat in="out1.fits out2.fits out3.fits out4.fits out5.fits out6.fits out7.fits out8.fits out9.fits out10.fits out11.fits" out=out.fits 
