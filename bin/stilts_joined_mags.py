#!python3
#
# 
#  stilts_joined_mags 
#   - remove columns not needed for the classification in the input catalogue
#   - add mag columns from version input catalogue which was matched to SDSS13 
#   - add columns with weighted value the magnitude
#   - remove the instrumental magnitudes once used to compute the mean
#  csh scripts in ~/bin
#  python script in /Users/data/catalogs/
#  data files in /Users/data/catalogs subdirectories: see below 
#  output data cover the SED of the UV-optical-IR for SUSS sources, with added Gaia PM and parallax
#
# NPMK 2023/ XMM2Athena WG7

import os
from astropy.table import Table
from astropy.io import fits
import numpy as np
import numpy.ma as ma

HOME=os.getenv('HOME')
BIN=HOME+'/bin'
WORK=os.getcwd()

delta_err = 1e-5 # minimum error (added to all errors)
fix_blank_error = True
blank_err_fill = 0.05

# from Atrium SUSS single source catalogue augmented with PanStars,SkyMapper, UKIDS, VISTA and ALLWISE
IN = "/Users/data/catalogs/atrium/XMMOM_SUSS5.0_Sources_v0_1_aug.fits"
# IN2 is an intermediate product by merging SUSS and SDSS16 in TopCat using CDS crossmatch
IN2 = "/Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_xSDSSxGaia_v0.1.csv"
# output file (both FITS and CSV)
OUT1 = "/Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_v0_1_traininginput.fits"
OUT = "/Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_v0_1_traininginput.csv"
TEMP1=WORK+'/temp1.fit'
#
def merge_band2(t, newname, mag1,err1,mag2,err2):
    """
    merge_band2:
    
    input parameters:
    t : astropy.table.Table 
    newname: string
      the new column name to be added. the derived error column has prepended e_
    (mag1,err1,mag2,err2): strings (for 2 bands)
       mag? will be the column name in t of a magnitude
       err? will be the corresponding error column name
    delta_err: real
       extra error term, mainly to avoid zero error   
    
    output:   
    compute the weighted magnitude and error for two input bands which have masked 
    values. The default method will only combine those rows with no masked value, thus
    not including the rows with one or two valid magnitudes. This program will treat 
    that case, so that the resulting columns have a value for each row that contains 
    at least one measurement.
    The resulting magnitude and error columns will be added to the table t.
    """
    newcol = np.copy(t[mag1]) # initialise new column with properties of mag1
    errcol = np.copy(t[err1]) # ditto for err1
    # now find number of inputs in each row (or use mask)
    q1 = np.isfinite(t[mag1])
    q2 = np.isfinite(t[mag2])
    if fix_blank_error:  
       p1 = np.isfinite(t[mag1]) & np.isnan(t[err1])        # has mag, but error isnan
       t[err1][p1] = blank_err_fill
       p2 = np.isfinite(t[mag2]) & np.isnan(t[err2])        # has mag, but error isnan
       t[err2][p2] = blank_err_fill
    qnum = np.asarray(q1,dtype=int) + np.asarray(q2,dtype=int)
    # now add the values from the second columns where first has no value
    q = (qnum == 1) & q2
    newcol[q] = t[mag2][q]
    errcol[q] = t[err2][q] 
    # now compute weighted mean and add (here we add delta_err to avoid NULLs)
    q = (qnum == 2)
    e1 = t[err1][q]*t[err1][q] + delta_err
    e2 = t[err2][q]*t[err1][q] + delta_err
    errcol[q] = 1./( 1./e1 + 1./e2 )
    newcol[q] = errcol[q]*( t[mag1][q]/e1 + t[mag2][q]/e2 )
    errcol[q] = np.sqrt(errcol[q])
    tt.add_columns([newcol,errcol],names=[newname,'e_'+newname])
    return t

def merge_band3(t, newname, mag1,err1,mag2,err2,mag3,err3):
    """
    merge_band3:
    
    input parameters:
    t : astropy.table.Table 
    newname: string
      the new column name to be added. the derived error column has prepended e_
    (mag1,err1,mag2,err2,mag3,err3): strings (for 3 bands)
       mag? will be the column name in t of a magnitude
       err? will be the corresponding error column name
    delta_err: real
       extra error term, mainly to avoid zero error   
    
    output:   
    compute the weighted magnitude and error for three input bands which have masked 
    values. The default method will only combine those rows with no masked value, thus
    not including the rows with one or two valid magnitudes. This program will treat 
    that case, so that the resulting columns have a value for each row that contains 
    at least one measurement.
    The resulting magnitude and error columns will be added to the table t.
    """
    from numpy import arange, sqrt,asarray
    newcol = np.copy(t[mag1]) # initialise new column with properties of mag1
    errcol = np.copy(t[err1]) # ditto for err1
    # now find number of inputs in each row (or use mask)
    q1 = np.isfinite(t[mag1])
    q2 = np.isfinite(t[mag2])
    q3 = np.isfinite(t[mag3])
    if fix_blank_error:  
       p1 = np.isfinite(t[mag1]) & np.isnan(t[err1])        # has mag, but error isnan
       t[err1][p1] = blank_err_fill
       p2 = np.isfinite(t[mag2]) & np.isnan(t[err2])        # has mag, but error isnan
       t[err2][p2] = blank_err_fill
       p3 = np.isfinite(t[mag3]) & np.isnan(t[err3])        # has mag, but error isnan
       t[err3][p3] = blank_err_fill

    qnum = np.asarray(q1,dtype=int) + np.asarray(q2,dtype=int)+ np.asarray(q3,dtype=int)
    # now add the values from the second columns where first has no value, same for third
    q = (qnum == 1) & q2
    newcol[q] = t[mag2][q]
    errcol[q] = t[err2][q]
    q = (qnum == 1) & q3
    newcol[q] = t[mag3][q]
    errcol[q] = t[err3][q]
    # now the case of three values: weighted mean and error; add delta_err
    q = qnum == 3 
    e1 = t[err1][q]*t[err1][q] + delta_err
    e2 = t[err2][q]*t[err2][q] + delta_err
    e3 = t[err3][q]*t[err3][q] + delta_err
    errcol[q] = 1./(1./e1 + 1./e2 + 1./e3)
    newcol[q] = errcol[q]*( t[mag1][q]/e1 + t[mag2][q]/e2 + t[mag3][q]/e3 )
    errcol[q] = np.sqrt(errcol[q])
    
    # now the case of only two values: split over 1 and 2, 1 and 3, or 2 and 3
    q = (qnum == 2) & q1 & q2
    e1 = t[err1][q]*t[err1][q] + delta_err
    e2 = t[err2][q]*t[err2][q] + delta_err
    #e3 = t[err3][q]*t[err3][q] + delta_err
    errcol[q] = 1./( 1./e1 + 1./e2 )
    newcol[q] = errcol[q]*( t[mag1][q]/e1 + t[mag2][q]/e2 )
    errcol[q] = np.sqrt(errcol[q])
    q = (qnum == 2) & q1 & q3
    e1 = t[err1][q]*t[err1][q] + delta_err
    #e2 = t[err2][q]*t[err2][q] + delta_err
    e3 = t[err3][q]*t[err3][q] + delta_err
    errcol[q] = 1./(1./e1 + 1./e3)
    newcol[q] = errcol[q]*( t[mag1][q]/e1 + t[mag3][q]/e3 )
    errcol[q] = np.sqrt(errcol[q])
    q = (qnum == 2) & q2 & q3
    #e1 = t[err1][q]*t[err1][q] + delta_err
    e2 = t[err2][q]*t[err2][q] + delta_err
    e3 = t[err3][q]*t[err3][q] + delta_err
    errcol[q] = 1./( 1./e2 + 1./e3 )
    newcol[q] = errcol[q]*( t[mag2][q]/e2 +t[mag3][q]/e3 )
    errcol[q] = np.sqrt(errcol[q])
        
    tt.add_columns([newcol,errcol],names=[newname,'e_'+newname])
    return t

def merge_band5(t, newname, mag1,err1, mag2,err2, mag3,err3, mag4,err4, mag5,err5):
    """
    merge_band5:
    
    input parameters:
    t : astropy.table.Table 
    newname: string
      the new column name to be added. the derived error column has prepended e_
    (mag1,err1,mag2,err2,mag3,err3,mag4,err4,mag5,err5): strings (for 5 bands)
       mag? will be the column name in t of a magnitude
       err? will be the corresponding error column name
    delta_err: real, environment variable
       extra error term, mainly to avoid zero error   
    
    output:   
    compute the weighted magnitude and error for three input bands which have masked 
    values. The default method will only combine those rows with no masked value, thus
    not including the rows with one or two valid magnitudes. This program will treat 
    that case, so that the resulting columns have a value for each row that contains 
    at least one measurement.
    The resulting magnitude and error columns will be added to the table t.
    """
    from numpy import arange, sqrt,asarray
    newcol = np.copy(t[mag1]) # initialise new column with properties of mag1
    errcol = np.copy(t[err1]) # ditto for err1
    # now find number of inputs in each row (or use mask)
    q1 = np.isfinite(t[mag1])
    q2 = np.isfinite(t[mag2])
    q3 = np.isfinite(t[mag3])
    q4 = np.isfinite(t[mag4])
    q5 = np.isfinite(t[mag5])
    if fix_blank_error:  
       p1 = np.isfinite(t[mag1]) & np.isnan(t[err1])        # has mag, but error isnan
       t[err1][p1] = blank_err_fill
       p2 = np.isfinite(t[mag2]) & np.isnan(t[err2])        # has mag, but error isnan
       t[err2][p2] = blank_err_fill
       p3 = np.isfinite(t[mag3]) & np.isnan(t[err3])        # has mag, but error isnan
       t[err3][p3] = blank_err_fill
       p4 = np.isfinite(t[mag4]) & np.isnan(t[err4])        # has mag, but error isnan
       t[err4][p4] = blank_err_fill
       p5 = np.isfinite(t[mag5]) & np.isnan(t[err5])        # has mag, but error isnan
       t[err5][p5] = blank_err_fill
       
    qnum = np.asarray(q1,dtype=int) + np.asarray(q2,dtype=int)+ np.asarray(q3,dtype=int)\
     + np.asarray(q4,dtype=int)+ np.asarray(q5,dtype=int)
    mm = [mag1,mag2,mag3,mag4,mag5] 
    ee = [err1,err2,err3,err4,err5]
    qq = [q1,q2,q3,q4,q5]
    # single value in row
    for k in np.arange(1,5):
        q = (qnum == 1) & qq[k]
        newcol[q] = t[mm[k]][q]
        errcol[q] = t[mm[k]][q]
    # two values in row 
    for k1 in np.arange(5):
        for k2 in arange(k1+1,5):
            q = (qnum == 2) & qq[k1] & qq[k2]
            errcol[q] = 1./( 1./(t[ee[k1]][q]*t[ee[k1]][q]+delta_err) + 1./(t[ee[k2]][q]*t[ee[k2]][q]+delta_err) )
            newcol[q] = errcol[q]*( t[mm[k1]][q]/(t[ee[k1]][q]*t[ee[k1]][q]+delta_err) + t[mm[k2]][q]/(t[ee[k2]][q]*t[ee[k2]][q]+delta_err) )
            errcol[q] = np.sqrt(errcol[q])
    # three values in row
    for k1 in np.arange(5):
        for k2 in arange(k1+1,5):
            for k3 in arange(k2+1,5):
                q = (qnum == 3) & qq[k1] & qq[k2] &qq[k3]
                errcol[q] = 1./( 1./(t[ee[k1]][q]*t[ee[k1]][q]+delta_err) + \
                                 1./(t[ee[k2]][q]*t[ee[k2]][q]+delta_err) + \
                                 1./(t[ee[k3]][q]*t[ee[k3]][q]+delta_err)  )
                newcol[q] = errcol[q]*( t[mm[k1]][q]/(t[ee[k1]][q]*t[ee[k1]][q]+delta_err) + \
                                        t[mm[k2]][q]/(t[ee[k2]][q]*t[ee[k2]][q]+delta_err) +
                                        t[mm[k3]][q]/(t[ee[k3]][q]*t[ee[k3]][q]+delta_err) )
                errcol[q] = np.sqrt(errcol[q])
    # four values in row
    for k1 in np.arange(5):
       for k2 in arange(k1+1,5):
           for k3 in arange(k2+1,5):
               for k4 in arange(k3+1,5):
                   q = (qnum == 4) & qq[k1] & qq[k2] & qq[k3] & qq[k4]
                   errcol[q] = 1./( 1./(t[ee[k1]][q]*t[ee[k1]][q]+delta_err) + \
                                    1./(t[ee[k2]][q]*t[ee[k2]][q]+delta_err) + \
                                    1./(t[ee[k4]][q]*t[ee[k4]][q]+delta_err) + \
                                    1./(t[ee[k3]][q]*t[ee[k3]][q]+delta_err)  )
                   newcol[q] = errcol[q]*( t[mm[k1]][q]/(t[ee[k1]][q]*t[ee[k1]][q]+delta_err) + \
                                           t[mm[k2]][q]/(t[ee[k2]][q]*t[ee[k2]][q]+delta_err) + \
                                           t[mm[k4]][q]/(t[ee[k4]][q]*t[ee[k4]][q]+delta_err) + \
                                          t[mm[k3]][q]/(t[ee[k3]][q]*t[ee[k3]][q]+delta_err) )
                   errcol[q] = np.sqrt(errcol[q])
    # five values in row
    q = qnum == 5
    errcol[q] = 1./( 1./(t[ee[1]][q]*t[ee[1]][q]+delta_err) + \
                     1./(t[ee[2]][q]*t[ee[2]][q]+delta_err) + \
                     1./(t[ee[4]][q]*t[ee[4]][q]+delta_err) + \
                     1./(t[ee[0]][q]*t[ee[0]][q]+delta_err) + \
                     1./(t[ee[3]][q]*t[ee[3]][q]+delta_err)  )
    newcol[q] = errcol[q]*( t[mm[1]][q]/(t[ee[1]][q]*t[ee[1]][q]+delta_err) + \
                            t[mm[2]][q]/(t[ee[2]][q]*t[ee[2]][q]+delta_err) + \
                            t[mm[3]][q]/(t[ee[3]][q]*t[ee[3]][q]+delta_err) + \
                            t[mm[4]][q]/(t[ee[4]][q]*t[ee[4]][q]+delta_err) + \
                            t[mm[0]][q]/(t[ee[0]][q]*t[ee[0]][q]+delta_err) )
    errcol[q] = np.sqrt(errcol[q])
    tt.add_columns([newcol,errcol],names=[newname,'e_'+newname])
    return t
    
# create a new work file from IN and IN2
command=BIN+"/stilts_joined_mags.sh "+TEMP1
status = os.system(command)
# read the catalogue
tt = Table.read(TEMP1)
t = ma.asarray(tt)
# tt = Table.read(TEMP1)
# t = np.copy(tt)
print (f"now joining data of each band from various catalogues into a single column")
#
# gmag (add PS, SM, SDSS)
merge_band3( t, 'gmag_M', 'PS_gmag', 'PS_e_gmag', 'SM_gmag', 'SM_e_gmag', 'gmag', 'e_gmag')
print (f" gmag (add PS, SM, SDSS) done \n")
#
# gPKmag Kron and Petrosian (PS, SM)
merge_band2( t, 'gPKmag_M','PS_gKmag', 'PS_e_gKmag', 'SM_gPmag', 'SM_e_gPmag') 
print (f" gPKmag Kron and Petrosian \(PS, SM\) done \n")
#
# umag (add OM, SM, SDSS) even though they have different band centres ...
merge_band3( t, 'umag_M', 'SM_umag', 'SM_e_gmag', 'U_AB_MAG', 'U_AB_MAG_ERR', 'umag', 'e_umag')
print (f" umag (add OM+SM) done\n")
#
# uPmag = SM_uPmag (i.e., not renamed)
#
# vmag (add OM and SM)
merge_band2( t,  'vmag_M', 'V_AB_MAG', 'V_AB_MAG_ERR',  'SM_vmag', 'SM_e_vmag') 
print (f" vmag (add OM and SM) done\n")
#
# vPmag = SM_vPamg (i.e., not renamed)
#
# rmag (add PS, SM, SDSS)
merge_band3( t, 'rmag_M', 'PS_rmag', 'PS_e_rmag', 'SM_rmag', 'SM_e_rmag', 'rmag', 'e_rmag')
print (f" rmag (add PS, SM, SDSS) done\n")
#
# rPKmag (add PS Kron and SM Petrosian)
merge_band2( t, 'rPKmag_M', 'PS_rKmag', 'PS_e_rKmag', 'SM_rPmag', 'SM_e_rPmag') 
print (f" rPKmag (add PS Kron and SM Petrosian) done\n")
#
# imag (add PS, SM, SDSS)
merge_band3( t,  'imag_M', 'PS_imag', 'PS_e_imag', 'SM_imag', 'SM_e_imag', 'imag', 'e_imag')
print (f" imag (add PS, SM, SDSS) done\n")
# 
# iPKmag (add PS Kron and SM Petrosian)
merge_band2( t,  'iPKmag_M', 'PS_iKmag', 'PS_e_iKmag', 'SM_iPmag', 'SM_e_iPmag')
print (f" iPKmag (add PS Kron and SM Petrosian) done\n")
#
# zmag (add PS, SM, UK, VI, SDSS)
merge_band5( t,  'zmag_M', 'PS_zmag', 'PS_e_zmag', 'SM_zmag', 'SM_e_zmag', 'UK_zmag', 'UK_e_zmag', 'VI_zmag', 'VI_e_zmag', 'zmag', 'e_zmag')
print (f" zmag (add PS, SM, UK, VI, SDSS) done\n")
#
# zPKmag (add PS Kron and SM, VI Petrosian)
merge_band3( t, 'zPKmag_M', 'PS_zKmag', 'PS_e_zKmag', 'SM_zPmag', 'SM_e_zPmag', 'UK_zPmag', 'UK_e_zPmag') 
print (f" zPKmag (add PS Kron and SM, VI Petrosian) done\n")
#
# ymag (add PS, UK, VI)
merge_band3( t, 'ymag_M', 'PS_ymag', 'PS_e_ymag', 'UK_ymag', 'UK_e_ymag', 'VI_ymag', 'VI_e_ymag')
print (f" ymag (add PS, UK, VI) done\n")
#
# yPKmag (add PS Kron and SM, VI Petrosian)
merge_band2( t,  'yPKmag_M', 'PS_yKmag', 'PS_e_yKmag','UK_yPmag', 'UK_e_yPmag')
print (f" yPKmag (add PS Kron and SM, VI Petrosian) done\n")
#
# jmag (add UK, VI)
merge_band2( t, 'jmag_M',  'UK_jmag', 'UK_e_jmag', 'VI_jmag', 'VI_e_jmag')
print (f" jmag (add UK, VI) done\n")
#
# jPmag = UK_jPmag (i.e., not renamed)
#
# hmag (add UK, VI)
merge_band2( t, 'hmag_M', 'UK_hmag', 'UK_e_hmag', 'VI_hmag', 'VI_e_hmag')
print (f" hmag \(add UK, VI\) done\n")
#
# hPmag = UK_hPmag (i.e., not renamed)
#
# kmag (add UK, VI)
merge_band2( t, 'kmag_M', 'UK_kmag', 'UK_e_kmag', 'VI_kmag', 'VI_e_kmag')
print (f" kmag (add UK, VI) done\n-------------")
#
# kPmag = UK_kPmag (i.e., not renamed)
#
print ("discarding t -- working on tt")
t = ''
tt.colnames
#
# no joins for WISE mag's WI_W1mag,WI_W2mag,WI_W3mag,WI_W4mag
# no joins for UVOT uvw2, uvm2, uvw1, b
#
# Question: should we add old catalogue with B-magnitudes like USNO-B1 ? Probably not.
#
#   remove catalogue mag columns no longer needed (keep  uvw1, gmag, zmag, WI_w1mag)
#   define colour terms 
#   remove merged mag columns no longer needed for classification
#
#t.remove_columns('PS_* SM_?mag SM_e_?mag UK_?mag UK_e_?mag VI_* '; \ 
tt.add_column(tt['UVW2_AB_MAG']-tt['UVM2_AB_MAG'], name="uvw2_uvm2" )
tt.add_column(tt['UVM2_AB_MAG']-tt['UVW1_AB_MAG'], name="uvm2_uvw2" )  
tt.add_column(tt['UVW1_AB_MAG']-tt['U_AB_MAG'], name="uvw1_u" )  
print ("added UVW2-UVM2, UVM2-UVW1, and UVW1-U")
tt.add_column(tt['B_AB_MAG']-tt['V_AB_MAG'], name="b_v" )  
tt.add_column(tt['gmag_M']-tt['rmag_M'], name="g_r" )  
tt.add_column(tt['rmag_M']-tt['imag_M'], name="r_i" )  
tt.add_column(tt['imag_M']-tt['zmag_M'], name="i_z" )  
print ("added B-V, g-r, r-i, i-z")
tt.add_column(tt['zmag_M']-tt['ymag_M'], name="z_y" )  
tt.add_column(tt['jmag_M']-tt['zmag_M'], name="j_z" )  
tt.add_column(tt['hmag_M']-tt['jmag_M'], name="h_j" )  
tt.add_column(tt['kmag_M']-tt['hmag_M'], name="k_h" )  
print ("added z-y, j-z, h-j, k-h")
tt.add_column(tt['WI_W1mag']-tt['zmag_M'], name="W1_z" )  
tt.add_column(tt['WI_W2mag']-tt['WI_W1mag'], name="W2_W1" )  
tt.add_column(tt['WI_W3mag']-tt['WI_W2mag'], name="W3_W2" )  
tt.add_column(tt['WI_W4mag']-tt['WI_W3mag'], name="W4_W3" )  
print ("added WISE W1-z,W2-W1, W3-W2,W4-W3")
tt.add_column(tt['gPKmag_M']-tt['SM_uPmag'], name="gPK_uP" )  
tt.add_column(tt['rPKmag_M']-tt['gPKmag_M'], name="rPK_gPK" )  
tt.add_column(tt['iPKmag_M']-tt['rPKmag_M'], name="iPK_rPK" )  
print ("added Petrosian/Kron gPK-uP, rPK-gPK, iPK-rPK")
tt.add_column(tt['gPKmag_M']-tt['gmag_M'], name="gPK_g" )  
tt.add_column(tt['rPKmag_M']-tt['rmag_M'], name="rPK_r" )  
tt.add_column(tt['iPKmag_M']-tt['imag_M'], name="iPK_i" )  
print ("added Petrosian/Kron to default gP/K-g, rP/K-r, iP/K-i, next: zPK-z")

tt.add_column(tt['zPKmag_M']-tt['zmag_M'], name="zPK_z" )  
 
#t.remove_columns( imag, jmag, hmag, kmag, WI_e_*, WI_W4mag, WI_W3mag, WI_W2mag )
print (f"writing {OUT1}")
tt.write(OUT1,overwrite=True)
#t.write(OUT,format='csv',overwrite=True)

