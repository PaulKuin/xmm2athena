#python3

#
# 
# Code for adding the observing epoch to the SUSS[SOURCES]:
#
#    mjd2epoch(mjd)
#    delta_time_tdb_minus_tt(ra,dec,obstimes)
#    apply_proper_motion(ra,dec,ref_epoch,pmra,pmde,obs_epoch)
#    degrees2sexagesimal(ra,dec,)
#    bary_centric_correction_to_SUMMARY(
#      suss = '/data/catalogs/XMM-OM-SUSS5.0.fits',
#      output='/data/catalogs/suss5.0_summary_added_tdb_minus_tt')
#    _add_epoch_col()
#
# Code for matching SUSS[SOURCES] per observation to the Gaia DR3 high 
# proper motion sources in the nearby field. Match after precessing the
# Gaia DR3 to the SUSS epoch:
#
#    gaia_cat(root='.',o1=0,o2=10)
#       (calls several shell scripts, including one provided by http://www.arches-fp7.eu)
#
# Fixes:
#    fix_corruped(table,fill=None,nozeros=True)
#    _fix_corrupted_cols(col,err="--",fill=np.nan,nozeros=True)
#    fix_coordinates(infile,outfile)
#
# Add class for each training set:
#
#    add_class(cat, classx={'STAR':0,'QSO':1,'GALAXY':2,},
#      star_training="/Users/data/catalogs/work/SUSS5_gaiadr3-parallaxoErr_gt_10_25Jan23.fits",
#      QSO_training="/Users/data/catalogs/suss_gaia_epic/QSO_v2.csv",
#      gal_training="/Users/data/catalogs/suss_gaia_epic/galaxy_v2.csv")
#

from numba import jit
import numpy as np
import numpy.ma as ma


def mjd2epoch(mjd):
    from astropy import time
    t = time.Time(mjd,format='mjd')
    return t.jyear


def delta_time_tdb_minus_tt(ra,dec,obstimes):
    '''  (SUSS SUMMARY table)
    Given the source/obsid position and (average) observation time in TT, 
    compute the time correction from Terrestrial Time to Barycentric time
    (TDB). This does not provide the fully relativistic TCB time which 
    corrects for the time as seen by an observer outside the local solar 
    system gravitational field. Accuracy is better than the time resolution 
    of XMM OM SUSS. 
    
    input parameters:
    
    ra : float 
       right ascension (ICRS,J2000) in deg
    dec : float
       declination (ICRS,J2000) in deg
    obsstimes : array 
       times in MJD based on TT
       
    output: 
       delta times as a Time object
       
    2021-05-17 npmk; based on astropy documentation
    '''
    import numpy as np
    from astropy import time, coordinates as coord, units as u
    
    pointing = coord.SkyCoord(ra,dec,unit=(u.deg,u.deg),frame='icrs')
    location = coord.EarthLocation.of_site('greenwich')
    tijden = time.Time(obstimes, format='mjd', scale='tt', location=location)
    #print (pointing)
    ltt_baryc = tijden.light_travel_time(pointing,location=location,
         kind='barycentric' )
    #ltt_helio = time.light_travel_time(pointing, kind='heliocentric' )
    # note: heliocenter will move during Earth orbit 
    corrtimes = tijden.tdb + ltt_baryc
    delta_time_tdb_minus_tt = corrtimes.mjd - tijden.mjd
    return delta_time_tdb_minus_tt


#@jit    
def apply_proper_motion(ra_eDR3,dec_eDR3,ref_epoch,pmra,pmde,obs_epoch):
    ''' (SUSS SOURCES table)
    Sky position at 'epoch', given Gaia eDR3 parameters
    
    input parameters:
    ra_eDR3,dec_eDR3 (ICRS) in degrees 
    ref_epoch (yr) the epoch of ra_eDR3, dec_eDR3 
    pmra, pmde proper motion in mas/yr
    obs_epoch (yr) epoch of the observation
    
    output 
    ra,dec (deg) at epoch 
    '''
    # delta in deg due to PM in RA and Dec *mas-> s*1000, s -> deg /3600. 
    dra  = (obs_epoch-ref_epoch)*pmra/3.6e6
    ddec = (obs_epoch-ref_epoch)*pmde/3.6e6
    # add the PM corrected seconds 
    ra_epoch  = ra_eDR3  + dra
    dec_epoch = dec_eDR3 + ddec
    # ensure 0=< ra_epoch <360
    ra_epoch[ra_epoch <  0.0 ] += 360.0
    ra_epoch[ra_epoch >= 360.] -= 360.0
    # ensure -90 =< dec_epoch =< 90
    #if (dec_epoch < -90.0) | (dec_epoch > +90.0): 
    #   raise ValueError(
    #   f"The declination is outside the range after applying the PMDE correction.")
    return ra_epoch, dec_epoch
    

def degrees2sexagesimal(ra,dec,as_string=False):
       ''' 
   simple code to convert RA, Dec from decimal degrees to sexigesimal
   
   input ra, dec in decimal degrees
   output ra,dec in sexagesimal strings
       '''
       import numpy as np
   
       ra = np.mod(ra,360.0) # 0<= ra < 360, positive
       rah=int(ra/15.)
       ram=int( (ra/15.-rah)*60.)
       ras=(ra/15.-rah-ram/60.)*3600
       newra="%02i:%02i:%05.2f" % (rah,np.abs(ram),np.abs(ras)) 
   
       sdec = np.sign(dec)
       ded=int(dec)
       dem=int((dec-ded)*60.)
       des=np.abs((dec-ded-dem/60.)*3600)
       ded=sdec*np.abs(ded)
       dem=np.abs(dem)
       newdec="%0i:%02i:%04.1f"%(ded,dem,des) 
       if as_string: 
          return newra,newdec
       else: 
          return rah,ram,ras,ded,dem,des

def bary_centric_correction_to_SUMMARY(
      suss = '/data/catalogs/XMM-OM-SUSS5.0.fits',
      output='/data/catalogs/suss5.0_summary_added_tdb_minus_tt'):
   from astropy.io import fits,ascii as ioascii
   from astropy import coordinates,units as u,constants
   from astropy.table import Table
   import xmm2athena as xmm

   fsuss = fits.open(suss,'update',memmap=True)
   summary = fsuss['SUMMARY']
   mjd = 0.5*(summary.data['MJD_START'] + summary.data['MJD_END'])
   ra = summary.data['RA_PNT']
   dec = summary.data['DEC_PNT']
   tab = Table(summary.data)

   # compute the barycentric time correction and write to a tmp file
   res = []
   for x,y,z in zip(ra,dec,mjd):
      res.append(xmm.delta_time_tdb_minus_tt(x,y,z)*86400.)
   tdb_minus_tt = np.asarray(res,dtype=float)

   # add the barycentric correction to the SUMMARY Table and write the new table 
   tab.add_column(tdb_minus_tt,name='TDB_MINUS_TT')
   tab.write(output+'.fits')

def _add_epoch_col():
   # the main "SOURCES" part of the SUSS does not include theepoch of observation. 
   # The epoch is in the 'SUMMARY" part. But each source observation is determined 
   # by spatial location and time of observation so the epoch should be added 
   # The epoch is needed to correct for proper motion of nearby sources, and is 
   # also useful for variability, and for practical reasons adding it to the 
   # 'SOURCES' part of SUSS is a good thing. 
   #
   # input file is a join of SUSS SUMMARY and Gaia DR3 
   #
   # returns input table, list of mid-time of observation obs_epoch derived from summary 
   #  and the positions corrected for proper motion using the input file
   # 
   # 
   from astropy.io import fits,ascii as ioascii
   from astropy import coordinates,units as u,constants
   from astropy.table import Table, vstack

   d2r = np.pi/180.0
   # this catalog is all gaia sources within 20' radius of suss obsids, with SUSS 
   #  Summary table data added
   dr3suss = fits.open('/data/catalogs/gaia.mpgt30.susssummary.fits')
   g = Table(dr3suss[1].data) # gaia source with PM >30mas within 20' of SUSS5.0 (RA_PNT,DEC_PNT)
   
   g.sort("RAJ2000")
   obs_epoch = []
   pos_epoch = []
   for k in np.arange(len(g)):
      gk = g[k]
      obs_ep =  mjd2epoch(0.5*(gk['MJD_START']+gk['MJD_END']))
      obs_epoch.append( obs_ep )
      # RA_ICRS, DE_ICRS are DR3 positions at epoch 2016.0
      ra_ep,dec_ep = apply_proper_motion(gk['RA_ICRS'],gk['DE_ICRS'],2016.0,gk['pmRA'],
           gk['pmDE'], obs_ep)
      pos_epoch.append([ ra_ep,dec_ep ] )    
   obs_epoch = np.array(obs_epoch)
   pos_epoch = np.array(pos_epoch)
   #
   # returns the input table g and the new columns to be added after review 
   #   obs_epoch is the positions at the observation epoch
   #   pos_epoch is the Gaia DR3 epoch of 2016 
   #
   return g,obs_epoch,pos_epoch
   
       
def gaia_cat(root='.',o1=0,o2=10):
    """
    make a SUSS catalogue file based on Gaia DR3 astrometry (JUNE 13 2022)
    
    (A) retrieve high PM sources from gaiadr3 => gaiadr3_highpm.fits 
        This was completed 2022-06-14 using the gaiadr3_lite, PM > 25 mas/yr
    
    (1) use obsid to find pointing RA,DEC, obs_epoch
        restrict processing using [o1:o2] range
    
    (B) loop over obsids :
        - extract cone from gaiadr3_highpm.fits centre pointing, radius 15 arc min > gtemp
        - extract suss data for obsid
        - compute the 2000 epoch for gtemp positions
        - change gtemp epoch from 2016.0 to suss using stilts function 
        - edit catalog upload columns if needed
        - xarches[cds xmatch] job to match them by position after upload (retain all suss records)
    
    #(2) for each obsid retrieve the Gaia entries ( minimally
        RA(J2000), DE(J2000), Position error, PMRA, PMDE, Gmag, epoch)
        in a radius of 11 arcminutes around pointing from GAIA eDR3
    
    #(3) apply proper motion correction to put data on the obs_epoch 
        of the SUSS data *** copy RA,DEC to new columns; filter on 
        'parallax' 'pm' since some sources don't have. 
        The apply PM correction on sources with PM, plx 
        
    
    #(4) iterate over the obsids and write Gaia for the SUSS catalog. 
    
    (C) assume all suss sources without high PM match coordinates have coordinate near
        epoch 2000 within 0.5 arcsec. 
        still by obsid: 
        - xarches [xmatch cds] match all records using epoch 2000 positions to gaiadr3  

    (D) merge the matched records 
    
    (E) clean up
    
    (5) write out a restricted columns catalog table [ra,dec,perr,pmra,pmde,gmag] 

    using conesearch in astroquery
    
    """
    
    import os
    import numpy as np 
    from astropy.table import Table, vstack
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    from astropy.io import fits, ascii as ioascii
    
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source_lite"
    radius = 15./60.  # u.Quantity(15./60., u.deg)

    gaia_epoch='2016.0'  # for DR3
    omcat = "/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits#1"
    groot1 = "/Users/data/catalogs/gaia_xmm/"    # temp products in adir, matched, archive
       # in archive are the by-obsid extractions of gaia eDR3 with PM corrections to the XMM epoch
       # in adir these have been corrected for sources with no plx/PM 
       # in matched we want to output from matching for each obsid of XMM-OM and Gaia from adir
    groot2 =  "/Users/data/catalogs/suss_gaia_epic/"  # catalog products
    fh = fits.open(groot2+'suss5.0_obsids_pointing.fits')
    obsid_pos_epoch_list = fh[1].data
    obsids  = obsid_pos_epoch_list['obsid'.upper()]
    ra_pnt  = obsid_pos_epoch_list['ra_pnt'.upper()]
    dec_pnt = obsid_pos_epoch_list['dec_pnt'.upper()]
    epochs  = mjd2epoch( obsid_pos_epoch_list['MJD_START'] )
    fh.close()
    kk = 0
    
    for ra_, dec_, obs_epoch, obsid in zip(ra_pnt[o1:o2],dec_pnt[o1:o2],epochs[o1:o2], obsids[o1:o2]): 
        print (f"\nnew record={kk+o1}, {ra_}, {dec_}, {obs_epoch}, {obsid} \n")
        kk+=1
        # (B)
               # select from high PM gaia catalog a subset gtemp
        gtemp = "gtemp.fit" 
        infile = groot2+"gaiadr3_highpm.fits"
        command= f"stilts_gaiadr3_for_obsid.sh {ra_} {dec_} {radius} {infile} {gtemp} "
        status = os.system(command)
        print ("status:",status," -- ",command)      
               # select from the suss a subset omtemp
        omtemp = "suss5_obsid.fit"
        command= f"stilts_suss_select_obsid.sh {obsid} '{omcat}' {omtemp} "
        status = os.system(command)
        print (status," -- ",command)      
              # add epoch 2000 coordinates to gtemp
              # add epoch obs_epoch coordinates to gtemp
        command= f"stilts_add_positions_at_epochs.sh {obs_epoch} {gtemp} "
        status = os.system(command)
        print (status," -- ",command)      
              # cds xmatch or xarches
        output1 = f'/Users/data/catalogs/xmerg_temp.fit'      
        output = f'/Users/data/catalogs/june/xmerg_{obsid}.fit'      
        command= f"xarches_sub.sh {output} "
        status = os.system(command)
        print (status," -- ",command)
        status = os.system("rm gtemp.fit suss5_obsid.fit")
               # drop gaia hig PM columns gaiadr3_earlier
       
        # after the initial series of calls, some failed files with just headers were seen
        # those were then processed seperately by identifying their position in the 
        # suss obsids pointing file and rerunning this program
# following this program the 10600+ files were concatenated using 'stilts_cat_files.sh'

# since the resulting positions posRA,posDec were matched accordong to the observed epoch, final 
# J2000, Epoch2000 positions will be adopted using the observed positions for the 
# objects that have no match in the highPM gaiadr3 file, after that we update the 
# positions of the sources with highPM match (and high probability of a good match 
# proba_AB > proba_A_B -- with the GaiaDR3 position. The column names will be 
# {RA,DEC}_EPJ2000   
# add matching column 4match=="OBSID*10,000,000+SRCNUM"
# add column PM, ePM
# remove columns : ePos_Epoch, posRA, posDec,
# suss_*_SRCDIST, suss_RA_HMS, suss_DEC_DMS, suss_*_RATE, suss_*_RATE_ERR, suss*VEGA*
# suss_POSERR ?
# gaiadr3*
# do xmatch cds to gaiaed3 epoch2000 for whole catalog keep photometry+error, astrometry
# do xmatch cds to sdss9 j2000 ep2000
# 
        
    # next stage: combine all the files suss5 + new position columns; drop gaiadr3 columns
    # copy the suss ra,dec for rows *without* parallax into the ra2000ep, dec2000ep columns
    #   note: these are the gaiadr3 matched sources
    # xmatch to gaiadr3, sdss9, sdss12, ps2, ...
    
def fix_coordinates(infile,outfile):
    """
    2022-08-20 fin sept 4. 
    
    the input being the SUSS matched to highest PM gaia, with astrometry parameters
    rest catalogue objects at observed position. XArches added the probability of a 
    good match prob_ab to the input file which is being used here. 
    
    -- the adjusted positions are marked J2000 but are at the observed epoch
    -- make new positions col RAJ2000Ep2000, DecJ2000Ep2000 by copying the observed 
       positions for all from the observed epoch except for the ones with astrometry.
    -- outfile must be of type <name>.fits   
     
    """    
    from astropy.table import Table
    from astropy.io import fits
    import numpy as np

    groot2 =  "/Users/data/catalogs/suss_gaia_epic/"  # catalog products
    fh = fits.open(groot2+'suss5.0_obsids_pointing.fits') # list of pointing for each 
                                                          # observation
    obsid_pos_epoch_list = fh[1].data
    obsids  = obsid_pos_epoch_list['obsid'.upper()]
    epochs  = mjd2epoch( obsid_pos_epoch_list['MJD_START'] )
    fh.close()

    f = fits.open(infile)
    t = Table(f[1].data)
    ra = t['posRA']
    ra.name='RAJ2000Ep2000'
    dec = t['posDec']
    dec.name='DecJ2000Ep2000'
    suss_epoch = []
    for k in range(len(t)):
        obsid = t[k]['suss_OBSID']
        xep = np.array([epochs[np.where(obsid == obsids)[0][0]]])
        suss_epoch.append(xep[0])
        if np.isfinite(t[k]['gaiadr3_parallax']):
            prob_ab = np.array([t[k]['proba_AB']])
            xra = np.array([ra[k]])
            xde = np.array([dec[k]])
            xpmra = np.array([t[k]['gaiadr3_pmra']])
            xpmde = np.array([t[k]['gaiadr3_pmdec']])
            if np.isfinite(prob_ab) :
                #print (k,xra[0],xde[0],xpmra[0],xpmde[0],xep[0],prob_ab[0])
                if prob_ab > 0.5:
                    yra,yde = apply_proper_motion(xra,xde,2000,xpmra,xpmde,xep)
                    ra[k] = yra[0]
                    dec[k] = yde[0]
                else:   # keep the SUSS position as it is not a good match
                    yra = t[k]['suss_RA']    
                    yde = t[k]['suss_DEC']
                    print (f"prob_ab < 0.5: {k}  - {ra[k]},{dec[k]}")
        
    t.add_column(suss_epoch,name="suss_epoch")
    f.close()
    t.write(outfile,overwrite=True)
    return t, ra, dec, suss_epoch
    
    
def add_class(cat, classx={'STAR':0,'QSO':1,'GALAXY':2,}, 
#  row='spCl', 
  star_training="/Users/data/catalogs/work/SUSS5_gaiadr3-parallaxoErr_gt_10_25Jan23.fits",
  QSO_training="/Users/data/catalogs/suss_gaia_epic/QSO_v2.csv",
  gal_training="/Users/data/catalogs/suss_gaia_epic/galaxy_v2.csv"):
    """
    input cat: file to add class column to according to classes 
          best input format csv
    
    see Hugo auto_classes.py for more complicated selections across many catalogues
    
    """
    from astropy.table import Table
    import numpy.ma as ma
    # load stars training set
    input_table = Table.read(cat)
    input_table['class'] = np.nan
    input_table = ma.asarray(input_table)

    star_set = Table.read(star_training)
    starids = star_set['IAUNAME_1']
    QSO_set = Table.read(QSO_training)
    qsoids = QSO_set['IAUNAME']
    gal_set = Table.read(gal_training)
    galids = gal_set['IAUNAME']

    dum, a1_idx, a2_adx = np.intersect1d(input_table['IAUNAME'],starids,return_indices=True)
    input_table['class'][a1_idx] = classx['STAR']
    dum, a1_idx, a2_adx = np.intersect1d(input_table['IAUNAME'],qsoids,return_indices=True)    
    input_table['class'][a1_idx] = classx['QSO']
    dum, a1_idx, a2_adx = np.intersect1d(input_table['IAUNAME'],galids,return_indices=True)    
    input_table['class'][a1_idx] = classx['GALAXY']
    
    #for C in classes.keys():
    #    input_table['class'][input_table[row] == C] = classes[C]
    
    return Table(input_table)
    
def _fix_corrupted_cols(col,err="--",fill=np.nan,nozeros=True):
    """
    Input: masked array column which has floating point values but is in string 
           format with some erronous row values, like '--', here in parm 'err' 
           parameter fill: fill value 
           
    Output: masked array column with erronous values removed and turned into 
           string format
           
    Example:  
        t = Table.read(file)
        c1 = t['U_AB_MAG_ERR']   # get string column
        o1 = xmm._fix_corrupted_cols(c1)  # fix column (fails if not string)
        t[o1.name] = o1  # replace column in t
        t.write(fileo)
           
    """   
    from astropy.table import MaskedColumn    
    import numpy as np
    name = col.name
    xa = ma.asarray(col)
    xa.fill_value = fill
    xav = xa.data
    xam = xa.mask
    try: # skip if number
        q = np.where(xav == err)
        # set mask 
        xam[q] = True
        # set values
        xav[q] = fill
    except:
        pass    
    if nozeros:
       xav = np.asarray(xav,dtype=float)
       xam[xav == 0] = True
       xav[xav == 0] = fill
    return MaskedColumn(data=xav,mask=xam,dtype=float,name=name) 

def fix_corruped(table,fill=None,nozeros=True):
    """
    the columns in the csv file that are being read in as strings or have 
    zeros
    
    Use fill=None then  astropy.ascii.write output in csv format 
      (writing fits has 1e20 in all masked fields). 
      
    example:
      from astropy.table import Table
      from cats import xmm2athena as xmm
      file = 'XMMOM_SUSS5.0_Sources_v0_1.csv'
      t = Table.read(file)  
      tnew = xmm.fix_corruped(t,fill=None) 
      tnew.write("XMMOM_SUSS5.0_Sources_v0.1.csv",format='csv')  
    """
    names=np.array(['UVW2_SRCDIST', 'UVM2_SRCDIST', 'UVW1_SRCDIST', 'U_SRCDIST',
       'B_SRCDIST', 'V_SRCDIST', 'UVW2_SIGNIF', 'UVM2_SIGNIF',
       'UVW1_SIGNIF', 'U_SIGNIF', 'B_SIGNIF', 'V_SIGNIF', 'UVW2_AB_FLUX',
       'UVW2_AB_FLUX_ERR', 'UVM2_AB_FLUX', 'UVM2_AB_FLUX_ERR',
       'UVW1_AB_FLUX', 'UVW1_AB_FLUX_ERR', 'U_AB_FLUX', 'U_AB_FLUX_ERR',
       'B_AB_FLUX', 'B_AB_FLUX_ERR', 'V_AB_FLUX', 'V_AB_FLUX_ERR',
       'UVW2_AB_MAG', 'UVW2_AB_MAG_ERR', 'UVM2_AB_MAG', 'UVM2_AB_MAG_ERR',
       'UVW1_AB_MAG', 'UVW1_AB_MAG_ERR', 'U_AB_MAG', 'U_AB_MAG_ERR',
       'B_AB_MAG', 'B_AB_MAG_ERR', 'V_AB_MAG', 'V_AB_MAG_ERR',
       'UVW2_AB_FLUX_MAX', 'UVW2_AB_MAG_MAX', 'UVM2_AB_FLUX_MAX',
       'UVM2_AB_MAG_MAX', 'UVW1_AB_FLUX_MAX', 'UVW1_AB_MAG_MAX',
       'U_AB_FLUX_MAX', 'U_AB_MAG_MAX', 'B_AB_FLUX_MAX', 'B_AB_MAG_MAX',
       'V_AB_FLUX_MAX', 'V_AB_MAG_MAX'], dtype='<U29')
    for x in names:
       col = table[ x ]
       fixd = _fix_corrupted_cols(col,fill=fill,nozeros=nozeros)
       table[x] = fixd      
    return table      
 


