
# code for updating / adding to the SUSS catalog SUMMARY and SOURCE 

#
# preprocessing
#
# same as make_file_variable below
catversion = 16
# this is where the code expects the topcat-full.jar executable
TOPCATPATH = "/Users/kuin/bin" 
"""
    print('match')
    # match 
    in1 = "allsources.txt"
    catin1 = "allsources_1.txt"
    catin2 = inputcat
    outcat = f"SimbadxSUSS5_variable_sources_v{catversion}.fits"
    xgaiavar = "xgaiavar.fits"
    xsimbad = "xsimbad.fits"
"""
# output catalogue for the variable sources
#
outcat = f"SimbadxSUSS5_variable_sources_v{catversion}.fits"
#
# *** inputs for make_file_variable( ) :
#
inputdir = '/Volumes/DATA11/data/catalogs/suss_gaia_epic/v3/'

# the single source catalogue input
inputcat = "XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits"

# filter to process
band = 'uvw1' # one of uvw2,uvm2,uvw1,u,b,v

# minimum required chi-squared reduced for inclusion (will reject data with quality # 0)
chi2_red_min = 5.0

# select only records from the single source catalogue which have good quality ? 
onlyqualzero = False

# set the minimum required source distance for inclusion in the variability 
min_srcdist  = 6.0

# minimum number of valid data points for inclusion into result
minnumber = 5

# name of file listing bad obsids
bad_obsidfile = "remove_obsidfile_here.txt"

# verbosity
chatter = 2

# names of RA, Dec columns at Epoch 2000, ICRS
RA2000 = "RAJ2000Ep2000" 
DE2000 = "DEJ2000Ep2000"

# make plots 
make_plots = True

# directory for putting the light curves fits files and the plots, relative to inputdir
var_lc = '../var_lc/'



##########################################################################################

#from numba import jit
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


#  #@jit nopython=False
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
    ra_epoch = np.asarray(ra_epoch)

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
       newdec="%02i:%02i:%04.1f"%(ded,dem,des) 
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


def obsid_epoch(file):
   """
   The 4XMM does not have a separate SUMMARY as the XMM-OM SUSS, so create one
   This is needed for the Gaia eDR3 first cut to areas around OBSIDs
   using the TOPCAT match to 20 arcmin from the mean 4XMM ra,dec of the OBSIDs
   
   this is followed with a match to the 4XMM within 2.5" for eDR3 and 4XMM RA,Dec
   records: 20,119 with PM > 75mas/yr; 2561 obsids.; matched obsid-OBS_ID=> 1990
     but many single observation, so remove those 404 groups left with multiple 
     observations: 4xmm-eDR3-PMgt75-obsid-cleaned_multipleobs.fits 
   """
   from astropy.io import fits,ascii as ioascii
   from astropy.table import Table, vstack
   # file was written from 4XMM after adding epoch and renaming ra,dec to *_4xmm   
   #file = '/Users/kuin/4xmm_obsids.fits'
   f = fits.open(file)
   f.info()
   t = Table(f[1].data)
   tg = t.group_by('OBS_ID')
   g = open('/Users/data/catalogs/4XMM_obsids.txt','w')
   N = len(tg.groups)
   for k in range(N):
       g.write(f"{tg.groups[k]['OBS_ID'][0]} {tg.groups[k]['RA_4xmm'].mean()} {tg.groups[k]['DEC_4xmm'].mean()} {tg.groups[k]['epoch_obs'].mean()} \n")
   g.close()


def xmatch_it():
   '''
   USE: out = xmatch_it()
   
   Using TopCat:
   
   first all sources in Gaia EDR3 with PM > 30 were selected and saved to file.
   
   next the selected EDR3 soures were matched with the obsid pointings in the Summary 
   of the SUSS5.0 within a range of 20' (1200"). That is the file we will work 
   from.  
   
   This program takes 1000x more time as it is not vectorised. See later version for DR3.
    
   '''
   from astropy.io import fits,ascii as ioascii
   from astropy import coordinates,units as u,constants
   from astropy.table import Table, vstack
   #import xmm2athena as xmm

   d2r = np.pi/180.0

   #edr3suss = fits.open('/data/catalogs/EDR3highPM-SUSS5match.fit')
   # this catalog is all gaia sources within 20' radius of suss obsids, with SUSS 
   #  Summary table data added
   edr3suss = fits.open('/data/catalogs/gaia.mpgt30.susssummary.fits')
   g = Table(edr3suss[1].data) # gaia source with PM >30mas within 20' of SUSS5.0 (RA_PNT,DEC_PNT)
   g.sort("RAJ2000")
   raproj = g['RA_ICRS']*np.cos(g['DEC_ICRS']*d2r)
   obs_epoch = []
   pos_epoch = []
   for k in np.arange(len(g)):
      gk = g[k]
      obs_ep =  mjd2epoch(0.5*(gk['MJD_START']+gk['MJD_END']))
      obs_epoch.append( obs_ep )
      ra_ep,dec_ep = apply_proper_motion(gk['RA_ICRS'],gk['DE_ICRS'],2000.0,gk['pmRA'],
           gk['pmDE'], obs_ep)
      pos_epoch.append([ ra_ep,dec_ep ] )    

   # the following is the SUSS srclist table
   suss = '/data/catalogs/XMM-OM-SUSS5.0.fits'
   fsuss = fits.open(suss,'update',memmap=True)
   srclist = fsuss['SRCLIST']
   s = Table(srclist.data)
   s.sort('RA') # to speed up matching
   
   # method:
   # in g calculate for each row, obs_epoch of SUSS obsid, 
   # then calculate gaia source position equinox J2000 at observed epoch using PM 
   # finally, match the gaia source position to SUSS positions to within 0.7" and 1.5"
   # in order to judge how different they were (not much difference so went with 1.5")
   # for each match, write data to an ascii file
   
   output = 'xmatch_GaiaEDR3_PM.gt.30_SUSS5.0.dat'
   # pre-process epochs and proper motion corrected positions 
   g = np.array(g)
   outputlist = _match(g, s, raproj, obs_epoch, pos_epoch, d2r )
   return outputlist

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
   # Note: see the script called addEpoch2SourceSUSS.sh
   #
   from astropy.io import fits,ascii as ioascii
   from astropy import coordinates,units as u,constants
   from astropy.table import Table, vstack

   d2r = np.pi/180.0
   # this catalog is all gaia sources within 20' radius of suss obsids, with SUSS 
   #  Summary table data added
   edr3suss = fits.open('/data/catalogs/gaia.mpgt30.susssummary.fits')
   g = Table(edr3suss[1].data) # gaia source with PM >30mas within 20' of SUSS5.0 (RA_PNT,DEC_PNT)
   g.sort("RAJ2000")
   #raproj = g['RA_ICRS']*np.cos(g['DE_ICRS']*d2r)
   obs_epoch = []
   pos_epoch = []
   for k in np.arange(len(g)):
      gk = g[k]
      obs_ep =  mjd2epoch(0.5*(gk['MJD_START']+gk['MJD_END']))
      obs_epoch.append( obs_ep )
      # RA_ICRS, DE_ICRS are EDR3 positions at epoch 2016.0
      ra_ep,dec_ep = apply_proper_motion(gk['RA_ICRS'],gk['DE_ICRS'],2016.0,gk['pmRA'],
           gk['pmDE'], obs_ep)
      pos_epoch.append([ ra_ep,dec_ep ] )    
   obs_epoch = np.array(obs_epoch)
   pos_epoch = np.array(pos_epoch)
   # g.add_column(obs_epoch,name='obs_epoch',)
   # g.add_column(pos_epoch[0],name='RA_OBS_EP')
   # g.add_column(pos_epoch[1],name='DE_OBS_EP')
   return g,obs_epoch,pos_epoch


def write(outputlist):   
   xo = open( output, 'w')
   for r in outputlist: 
       xo.write(r,'\n')
   xo.close()
   print (f"outputlist written to {output}")
   restab = Table(outputlist)
   restab.write('xmatch_GaiaEDR3_PM.gt.30_SUSS5.0.fits',overwrite=True)
   
# @jit
def _match(g, s, raproj, obs_ep, pos_epoch, d2r):
   # WG7 
   # too slow - need to  use TopCat/
   #d2r = np.pi/180.0
   radius = 1.0/3600.
   radius2 = radius*radius
   tablemade = False  
   outputlist = []
   
   for k in np.arange(100):
   #for k in arange(len(g)):
      gk = g[k]

      obs_epoch = obs_ep[k] #mjd2epoch(0.5*(gk['MJD_START']+gk['MJD_END']))
      ra_ep,dec_ep = pos_epoch[k] # apply_proper_motion(gk['RAJ2000'],gk['DEJ2000'],2000.0,gk['pmRA'], # gk['pmDE'], obs_epoch)
      
      q = (np.abs(raproj-ra_ep*np.cos(dec_ep*d2r)) < radius) & (np.abs(s['DEC']-dec_ep) < radius)
      if len(q) > 0:
        ss = s[q] 
        distance2 = ((raproj[q] - ra_ep*np.cos(dec_ep*d2r))**2 + (ss['DEC'] - dec_ep)**2)
        rq = distance2 < radius2
        if len(rq) > 0: 
           result = ss[rq]
           dist = np.sqrt(distance2[rq])
           result.add_column(dist,name='distmatch')
           result.sort('distmatch')
           #if tablemade:
           #   tmp = restab
           #   restab = vstack([tmp, result],join_type='exact') 
           #else: 
           #   restab = result
           #   tablemade = True
           m = 1
           for r in result: 
              outputlist.append([gk, r, [ra_ep, dec_ep, obs_epoch, m]])
              m = m+1
              #xo.write(f"{gk['EDR3Name']} {gk['Source']} {gk['RAJ2000']} {gk['DEJ2000']}  {r['IAUNAME']} {gk['OBSID']} {r['N_SUMMARY']} {r['RA']} {r['DEC']}   {ra_ep} {dec_ep}  {r['distmatch']}\n")  # i.e., result
   return outputlist   
      

def skydistance(ra1,dec1,ra2,dec2):
   """
   input in degrees, output in arcsec 
   (except at poles) 
   """
   import numpy as np
   x1 = ra1 / np.cos(np.pi*dec1/180.0)
   x2 = ra2 / np.cos(np.pi*dec2/180.0)
   d = (x1-x2)**2 + (dec1-dec2)**2
   return np.sqrt(d)*3600.0
   
def match_negative_refcat(infile='4xmmdr11_obsids+refcat.fit', maxdist=240.0,
       matchup='neg_refcat_matched_to_positive_refcat-v2.txt',
       nomatchup='neg_refcat_no_matches-v2.txt', dk = 78,
       match_swift='neg_refcat_match_swift-v2.txt',
       cross_match='neg_refcat_all_matches-v2.txt',
       not_any_match='neg_refcat_not_any_match-v2.txt'):
    """
    XMM SAS WG4 
    here we input the xxmdr11_obsid+refcat.fit file, 
    sort on ra, 
    find the indices for the negative refcat entries.
    then for each of those records, search nearby entries in ra,dec -> skydistance 
       search within dk records
       if within maxdist arcsec, then go to matchup list, 
       if within none within maxdist, then put negative refcat record in nomatchup list 
       - also check if matches swift position 
    """   
    from astropy.io import fits, ascii as ioascii
    from astropy.table import Table
    nmatch = 0
    f = fits.open(infile)
    f.info()
    match_out = open(matchup,'w')
    nomatch_out = open(nomatchup,'w')  
    swmatch = open(match_swift,'w')
    allmatch = open(cross_match,'w')
    without_any_match = open(not_any_match,'w')
    t = Table(f[1].data)
    t.sort("RA")
    negrefcat = np.where(t['REFCAT'] < 0 )[0]
    
    for k in negrefcat[:]:
        #print (60*"=!")
        #print (k)
        #print (f"{t[k]}\ncompare to:\n")
        k1 = np.max([(k - dk), 0]) 
        k2 = np.min([(k + dk), len(t)])
        small_t = t[k1:k2]
        small_t = small_t[ small_t['REFCAT'] > 0 ] # only the positive ones
        if len(small_t) <= 0 : 
            raise RuntimeError(f"number of entries {dk} around record {t[k]} is too small ")
        #print (small_t)    
        found_match=[]
        ra1 = t['RA'][k]
        de1 = t['DEC'][k]    
        for h in small_t:
            ra2 = h['RA']
            de2 = h['DEC']
            dd = skydistance(ra1,de1,ra2,de2)   
            #print (f"distance ({ra1},{de1}) ({ra2},{de2}) = {dd}")
            if dd < maxdist:
                found_match.append( (*h.as_void(),"\n"))
        if len(found_match) > 0:
            nmatch += 1
            match_out.write(f"{nmatch} {t[k].as_void()} \n")
            #for rec in found_match:
            #    match_out.write(f"{nmatch} match {rec.data} \n")
            match_out.write(f"{nmatch} matched_good_EPIC: {found_match} \n")
        else:    
            nmatch += 1
            nomatch_out.write(f"{nmatch}  {t[k].as_void()} \n")
        N = check_swift(ra1,de1)
        if N > 0: 
            swmatch.write( f"{nmatch} N_SwiftObsids={N}  {t[k].as_void()} \n")
            if len(found_match) > 0:
               allmatch.write(f"{nmatch} {t[k].as_void()}\n")
            else:
               without_any_match.write(f"{nmatch} {t[k].as_void()}\n")  
    match_out.close()
    nomatch_out.close()
    swmatch.close()
    allmatch.close()
    without_any_match.close()

def check_swift(ra,dec,radius=6.0):
    """
    XMM SAS study, compare 
    
    check if a source is in the Swift XRT archive
    within a 6 arcmin radius
    """
    from swifttools.swift_too import Swift_ObsQuery
    import numpy as np
    from astropy import coordinates,units
    query = Swift_ObsQuery()
    query.ra, query.dec = ra, dec
    pos = coordinates.SkyCoord(ra,dec,frame=coordinates.ICRS,unit='deg')
    query.skycoord = pos
    query.radius = radius / 60.
    xx = query.submit()
    N = len(query)
    if xx and (N > 0):
        return N
    else:
        return -1

def astrometry_match(obsid_bad,obsid_good,root='/Users/data/XMM/',workdir=None,download=True):
    """
    XMM SAS WG4 astrometry 
    
    Scripting the XMM work
    
    (1) create a new working directory (workdir='xxx')
    (2) download the data from the archive (needs obsids)
    (3) untar/unzip OBSMLI files
    (4) edit required headers (remove kwds in bad, rename cols in good)
    (5) display the images + regions 
    (6) initialise SAS environment
    (7) run catcorr on obsid_bad OBSML1 data using the obsid_good one
    (8) display the header of the, possibly fixed, obsid_bad
    2022-02-01 NPMK
    """
    #import numpy as np
    import os
    #from astropy.io import fits
    #from matplotlib import 
    
    if workdir is None: 
        raise RuntimeError(f"you need to specify a work directory containing the good and \nbad obsids starting from {root}.\n")
    else:
        #test for present 
        os.system(f"mkdir 0277 {root}{workdir}")
    
    if download:
        print ("downloading the obsid data")
        command = f'cd {root}{workdir};curl -o files.tar "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsid_bad}";tar xf files.tar'
        if not os.system(command): 
            raise RuntimeError(f"problem downloading {obsid_bad}") 
           
        command = f'cd {root}{workdir};curl -o files.tar "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsid_good}";tar xf files.tar'
        if not os.system(command): 
            raise RuntimeError(f"problem downloading {obsid_good}") 
    
    bad_obsmli_dir = f"{root}{workdir}/{obsid_bad}/pps/"
    good_obsmli_dir = f"{root}{workdir}/{obsid_good}/pps/"
    bad_obsmli_ = f"{bad_obsmli_dir}/P{obsid_bad}EPX000OBSMLI0000."
    good_obsmli_ = f"{good_obsmli_dir}/P{obsid_good}EPX000OBSMLI0000."
    
    # rename FTZ to fit.gz
    os.system(f"mv {bad_obsmli_}FTZ {bad_obsmli_}fit.gz;gunzip {bad_obsmli_}fit.gz")
    os.system(f"mv {good_obsmli_}FTZ {good_obsmli_}fit.gz;gunzip {good_obsmli_}fit.gz")
    
    fix_bad_(f"{bad_obsmli_}fit")
    goodfile = fix_good_(f"{good_obsmli_}fit",root=bad_obsmli_dir)
    print (f"\ncatalog file is now {goodfile}\n")
    
    img_bad  = f"{bad_obsmli_dir}/P{obsid_bad}EPX000OIMAGE8000.FTZ"
    img_good = f"{good_obsmli_dir}/pps/P{obsid_good}EPX000OIMAGE8000.FTZ"
    
    fix_fitsimg(img_bad)
    fix_fitsimg(img_good)

    srcreg_bad  = f"{bad_obsmli_dir}/P{obsid_bad}EPX000REGION0000.ASC"
    srcreg_good = f"{good_obsmli_dir}/pps/P{obsid_good}EPX000REHGION000.ASC"
    
    display_img(img_bad,srcreg_bad,img_good,srcreg_good)
    
    try:
        odfdir = f"{root}{workdir}/{obsid_bad}/odf/"
        print (f"\ntrying to untar/zip {odfdir}{obsid_bad}.tar.gz \n   ")
        os.system(f"cd {odfdir};tar xzvf {odfdir}{obsid_bad}.tar.gz")
        a8 = os.listdir(f"{odfdir}")
        for x8 in a8:
           if '.TAR' in x8:
               os.system(f"cd {odfdir} ;tar xf {odfdir}{x8}")
    except:
        pass    
    
    print (f"\nrunning catcorr on {obsid_bad} using {obsid_good}")
    run_catcorr(bad_obsmli_,goodfile,bad_obsmli_dir,mingood=2,syserr=1.0,maxoffset=10.0,V=7)
    
    # check the header to see if the astrometry was successfully corrected
    
    # done
    
def run_catcorr(bad_obsmli_,good_obsmli_,bad_obsmli_dir,mingood=2,syserr=1.0,maxoffset=10.0,V=4):
    """
    XMM SAS 
    
    initialise_sas(workdir) in calling program
    then run catcorr
    """
    import os
    odfdir = bad_obsmli_dir.split('/pps/')[0] + "/odf"
    command = f"cd {odfdir};"+\
        f"export SAS_ODF={odfdir};cifbuild;export SAS_CCF=`pwd`/ccf.cif;"+\
        f"odfingest;echo odfingest completed;export SAS_ODF=`pwd`/`ls -1 *SUM.SAS`;sasversion;"+\
        f"/Users/scisoft/sas/xmmsas_20210317_1624/bin/catcorr  srclistset={bad_obsmli_}fit "+\
        f"catset={good_obsmli_} mingood={mingood} syserr={syserr} maxoffset={maxoffset} -V {V};"
    print (command)    
    os.system(command)    
    # or write to a file, and use, e.g., execvpe() ut how to pass env? example needed
    
    
def display_img(bad, badreg, good, goodreg):
    # XMM SAS  stub
    print ("\nTBD Display image\n")    
    
def fix_fitsimg(file):
    """
    XMM SAS 
    
    add keyword RADECSYS FK5
    """
    import os
    file = decompress_ftz(file)
    os.system(f"fthedit {file} operation=add keyword=RADECSYS value='FK5'")
    
def decompress_ftz(file): 
    # XMM SAS 
    
    import os
    f1 = file.split('.')  
    if 'FTZ' in f1:
       f1 = file.split('.FTZ')[0]
       os.system(f"mv {file} {f1}.fit.gz;gunzip {f1}.fit.gz"  )
       file = f"{f1}.fit.gz"
    return file
    
def fix_bad_(file):
    """
    WG4 XMM SAS: study aspect corrections 
    
    check POSCORROK = False
    remove the keywords that were input in previous call of catcorr
    """
    import os
    file = decompress_ftz(file)
    tempfile = '_fix1_'
    f = open (tempfile,"w")
    f.write(f"-RAOFFSET \n")
    f.write(f"-DEOFFSET \n")
    f.write(f"-ROT_CORR \n")
    f.write(f"-RAOFFERR \n")
    f.write(f"-DEOFFERR \n")
    f.write(f"-ROT_ERR \n")
    f.write(f"-POFFSET \n")
    f.write(f"-POFFERR \n")
    f.write(f"-LIK_HOOD \n")
    #f.write(f"-LIKHOOD0 \n")
    f.write(f"-NMATCHES \n")
    f.write(f"-POSCOROK \n")
    f.write(f"-ASTCORR \n")
    f.write(f"-REFCAT \n")
    f.write(f"COMMENT attempt to remove keywords astrometry")
    f.close() 
    os.system(f"fthedit {file} @{tempfile}; rm {tempfile}")   
    # remove column SYSERRCC ? 

def fix_good_(file,root='.'):
    """
    import os
    check POSCORROK = True
    rename columns RA_CORR to CAT_RA, DEC_CORR to CAT_DEC, and RADEC_ERR to CAT_RADEC_ERR
    (possibly rename EXTNAME to EPXMATCH)
    """
    import os
    from astropy.io import fits
    from astropy.table import Table
    
    #file = decompress_ftz(file)
    f = fits.open(file)
    t = Table(f['SRCLIST'].data)
    ra = t['RA_CORR']
    ra.name = 'CAT_RA'
    de = t['DEC_CORR']
    de.name='CAT_DEC'
    err = t['RADEC_ERR']
    err.name='CAT_RADEC_ERR'
    tt = Table([ra,de,err])
    tt.write(output=root+'epiccat.fit',format='fits',overwrite=True)
    f.close()
    f = fits.open(root+'epiccat.fit',mode='update')
    f[1].header['EXTNAME']='EPICCAT'
    f.flush()
    f.close()
    return root+'epiccat.fit'    
    
def gaia_cat(root='.',o1=0,o2=10):
    """
    make a catalogue file based on Gaia DR3 (JUNE 13 2022)
    
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
 #   groot3 =  groot1+"adir/"
 #   groot4 =  groot1+"matched/"
    
    # (1)
    # ra_pnt, dec_pnt, epochs, obsids = pnt_from_obsid()  this was the XMM EPIC list now the ones that were not included:
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
              
        ''' eDR3 version 
           print (f"processing {obsid} pointing=({ra_},{dec_}) with epoch={obs_epoch} ...")
           gin = groot3+"gaiaedr3_for_"+f"{obsid:010}"+f"_a.vot"
   
           coord = SkyCoord(ra=ra_,dec=dec_,unit=(u.degree,u.degree),frame='icrs')
           j = Gaia.cone_search_async(coord, radius)
           jr = j.get_results()
           rep = jr['ref_epoch'][0]  # Problem here 
           n = len(jr['ra'])
           jr.add_column(jr['ra'], index=6, name=f"ra{gaia_epoch}" )   
           jr.add_column(jr['dec'],index=9, name=f"dec{gaia_epoch}") 
           
           raq, decq = apply_proper_motion(jr['ra'], jr['dec'],rep, jr['pmra'],jr['pmdec'],obs_epoch)
           q = np.isfinite(raq)
           jr['ra'][q] = raq[q]
           jr['dec'][q] = decq[q]
           jr['ref_epoch'] += obs_epoch - rep 
           
           # update the epoch of ra,dec
           jr.write(f"{gin}",format='votable',table_id=f"{obsid}",overwrite=True) 
           print (f"wrote table {gin}")
           
              # match to SUSS using bash script
           stemp = 'temp_'+f"{obsid:010}"+'.fits'  # SUSS 
           omtemp = groot4+stemp
           omfinal = groot4+"suss5.0_gaiaedr3_"+f"{obsid:010}"+".fits"
           command = f"stilts_match_suss_gaia.sh "+f"{obsid:010}"
           #command = f"stilts_suss_select_obsid.sh {obsid} '{omcat}#1' {omtemp} "
           status = os.system(command)
           print (status," -- ",command)      
           #command = f"stilts_match_by_obsid.sh {omtemp} {gin}  {omfinal} "
           #status = os.system(command)
           #print (status," -- ",command)
       '''
       
       
       
"""              this worked also to match the OBSID string to column:
       obsid_short=int(obsid)
       #extract the obsid data from SUSS
       cmd = f"select matches(OBSID,padWithZeros({obsid_short},10))"
       print (cmd)
       command=f'stilts tpipe in={SUSS} out={suss_obsid_tmp.fit}  cmd={cmd}'  
       print (command)     
       status = os.system(command)
       print (status)
       go = False

======

groot = "/Users/data/catalogs/gaia_xmm/"
obslist = "/Users/data/catalogs/gaia_xmm/suss5.0_obsids_pointing.fits"
suss = "/Users/data/catalogs/"
alias stilts = "java -jar /Users/kuin/bin/topcat-full.jar -stilts "

fix the error in the renamed column labels:
 rename columns use something like:

cd /Users/data/catalogs/gaia_xmm/
a = os.listdir('/Users/data/catalogs/gaia_xmm/')
for x in a:
    os.system(f"cp {x} test.vot")
    command=f"java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=test.vot"+\
    f" cmd='colmeta -name RA2015.5  $8;colmeta -name DEC2015.5  $10' out=testout.vot"
    os.system(command)
    os.system(f"mv testout.vot {x}")

from astropy.io.votable import parse_single_table, from_table, writeto
from astropy.table import Table
for x in a:
    xout = x.split('.')[0]+"_a.vot"
    t = parse_single_table(x) 
    data = t.array
    plx = data['parallax']
    data['ra'][plx.mask] = data['RA2015.5'][plx.mask]
    data['dec'][plx.mask] = data['DEC2015.5'][plx.mask]
    tab = Table (data)
    tab.rename_columns((tab.colnames[7],tab.colnames[9]),('RA2015.5','DEC2015.5'))
    votable = from_table(tab)
    writeto(votable,xout)
    print (f"completed {xout}")
    
    
========

def make_suss_gaia():   # gaia-eDR3 code 
    from astropy.io import fits, ascii as ioascii
    from astropy.table import Table

    groot1 = "/Users/data/catalogs/gaia_xmm/"    # temp products in adir, matched, archive
       # in archive are the by-obsid extractions of gaia eDR3 with PM corrections to the XMM epoch
       # in adir these have been corrected for sources with no plx/PM 
       # in matched we want to output from matching for each obsid of XMM-OM and Gaia from adir
    groot2 =  "/Users/data/catalogs/suss_gaia_epic/"  # catalog products
    groot3 =  groot1+"adir/"
    groot4 =  groot1+"matched/"
    obslist = "/Users/data/catalogs/suss_gaia_epic/suss5.0_obsids_pointing.fits" #list of OBSIDs
    omcat = "/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits"
    stilts = f"java -jar /Users/kuin/bin/topcat-full.jar -stilts "
    
    f = fits.open(obslist)
    #f.info()
    obs = f['SUMMARY'].data["OBSID"] # named array
    f.close()
    
    doall = False
    
    if doall:
       f = fits.open(omcat)

    for o in obs:
        stemp = 'temp_'+o+'.fits'  # SUSS 
     if doall:  
          q = f[1].data["OBSID"] == o

          t = Table(f[1].data[q])
          t.keep_columns(["IAUNAME","N_SUMMARY","OBSID","SRCNUM","RA","DEC",])
          stemp = 'temp_'+o+'.fits'  # SUSS 
          gtemp = groot+"gaiaedr3_for_"+o+".vot"
          mout = groot+"suss_gaia_match_"+o+".vot"
          mfinal = groot+
          t.write(stemp)
          command = "{stilts} tmatch2 in1="+stemp+" in2="+gtemp+"  out="+mout+\
           ' matcher=sky values1="RA/deg DEC/deg" values2="ra dec" '+\
           ' params=1.0 find=best1 '

        gin = groot3+"gaiaedr3_for_"+o+"_a.vot"
        stemp = 'temp_'+o+'.fits'  # SUSS 
        omtemp = groot4+stemp
        omfinal = groot4+"suss5.0_gaiaedr3_"+o+".fits"

        shortobsid = np.int(o)
        command = f"SHORTOBSID=`echo {shortobsid} | bc`; "+f"{stilts} tpipe in="+omcat+"#1  out="+omtemp+" cmd=\'"+\
          "select matches(OBSID,padWithZeros('${SHORTOBSID}',10))';echo $SHORTOBSID"
        status = os.system(command)
        print (status," -- ",command)
      
      #stilts tmatch2 in1=suss_0000110101.fits in2=gaiaedr3_for_0000110101_a.vot \
      # out=suss5.0_gaiaedr3_0000110101.fits  matcher=sky values1="RA DEC" values2="ra dec" params=1.0 find=best1
        command = "{stilts} tmatch2 in1="+omtemp+" in2="+gin+"  out="+omfinal+\
           ' matcher=sky values1="RA DEC" values2="ra dec" params=1.0 find=best1 '
        status = os.system(command)
        print (status," -- ",command)
========== 
 alternatively
 
 for o in obs[100:3000]:
    ...:     os.system(f"stilts_match_suss_gaia.sh {o}")

   
    for o in obs[:10]:
      gin = groot3+"gaiaedr3_for_"+o+"_a.vot"
      stemp = 'temp_'+o+'.fits'  # SUSS 
      omtemp = groot4+stemp
      omfinal = groot4+"suss5.0_gaiaedr3_"+o+".fits"
      command = f"stilts_suss_select_obsid.sh {o} '{omcat}#1' {omtemp} "
      status = os.system(command)
      print (status," -- ",command)      
      command = f"stilts_match_by_obsid.sh {omtemp} {gin}  {omfinal} "
      status = os.system(command)
      print (status," -- ",command)
       
     
=============       
       
    query = "SELECT gaia_source.ra,gaia_source.ra_error,"\
    "gaia_source.dec,gaia_source.dec_error,gaia_source.pmra,gaia_source.pmra_error,"\
    "gaia_source.pmdec,gaia_source.pmdec_error,"\
    "gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.bp_rp,"\ 
    "DISTANCE(POINT(ra_pnt,dec_pnt), POINT(ra,dec) ) AS ang_sep"\
    "FROM gaiaedr3.gaia_source "\
    "WHERE  1 = CONTAINS( POINT( ra_pnt, dec_pnt), CIRCLE(ra, dec, 17./60. )) "\
    "ORDER BY ang_sep",
    job = Gaia.launch_job_async( query, dump_to_file=True, 
        output_file=f"gaia_edr3_for_{obsid}_{ra_pnt}_{dec_pnt}.vot", 
    output_format='votable')
    
    res = job.get_results()
    
    #file = decompress_ftz(file)
    f = fits.open(output_file)
    t = Table(f['SRCLIST'].data)
    ra = t['RA_CORR']
    ra.name = 'CAT_RA'
    de = t['DEC_CORR']
    de.name='CAT_DEC'
    err = t['RADEC_ERR']
    err.name='CAT_RADEC_ERR'
    tt = Table([ra,de,err])
    tt.write(output=root+'epiccat.fit',format='fits',overwrite=True)
    f.close()
    f = fits.open(root+'epiccat.fit',mode='update')
    f[1].header['EXTNAME']='GAIAeDR3'
    f.flush()
    f.close()
    return root+'epiccat.fit'    
"""
    
def pnt_from_obsid(obsid=None, obsidlist='/users/data/catalogs/suss_gaia_epic/4xmmdr11_obslist.fit'):
    """
    obsid is 10 bytes long text
    return all rows if obsid == None
    """
    from astropy.table import Table
    from astropy.io import fits
    ft = fits.getdata(obsidlist, 1 )
    #t = Table.read(obsidlist,hdu=1)
    t = Table(ft) 
    if obsid != None:
        q = t['OBS_ID'] == obsid
        tq = t[q]   
        if len(tq) != 1: 
           print (f"ERROR: obsid {obsid} was not found in {file}\n")
    else:
        tq = t      
    # pointing, position angle
    epoch = mjd2epoch(0.5*(tq['MJD_START']+tq['MJD_STOP']))
    return tq['RA'],tq['DEC'],epoch,tq['OBS_ID']

    
def fix_coordinates(infile,outfile):
    """
    2022-08-20 fin sept 4. 
    
    the input being the SUSS matched to highest PM gaia, with astrometry parameters
    rest catalogue objects at observed position
    -- the adjusted positions are marked J2000 but are in fact at the observed epoch
    -- make new positions col RAJ2000Ep2000, DecJ2000Ep2000 by copying the observed 
       positions for all from the observed epoch except for the ones with astrometry.
    -- outfile must be of type <name>.fits   
    """    
    from astropy.table import Table
    from astropy.io import fits
    import numpy as np

    groot2 =  "/Users/data/catalogs/suss_gaia_epic/"  # catalog products
    fh = fits.open(groot2+'suss5.0_obsids_pointing.fits')
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
    #for k in range(2000,4000,1):
        obsid = t[k]['suss_OBSID']
        xep = np.array([epochs[np.where(obsid == obsids)[0][0]]])
        suss_epoch.append(xep[0])
        if np.isfinite(t[k]['gaiadr3_parallax']):
            #print (t[k])
            #obsid = t[k]['suss_OBSID']
            #print (f"OBSID={obsid}")
            prob_ab = np.array([t[k]['proba_AB']])
            xra = np.array([ra[k]])
            xde = np.array([dec[k]])
            #xep = np.array([epochs[np.where(obsid == obsids)[0][0]]])
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
        
    #t.add_columns([ra,dec],indexes=[1,2],)
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
 
    
''' process the catalogs -- notes 

   Processing : locating the high proper motion sources by matching with Gaia eDR3


%load_ext autoreload
autoreload 2

You might also find the Table hstack() function useful. This will let you add 
all the columns of one table into another table. 

  https://docs.astropy.org/en/stable/table/operations.html#id6

In your case this might look like:

from astropy.table import hstack
tab = hstack([tab, c1])

note: way too slow!
   
   # writing to a FITS table : Column definitions 
   cols = summary.columns
   c1 = fits.Column(name='TDB_MINUS_TT', array=col_tdb_minus,unit='s', format='E')
   c2 =  ioascii.read('/data/catalogs/suss5.0_summary_col_tdb_minus_tt.dat')
   cols.add_col(c2['col_tdb_minus'])
   
   # make a new fits bintable extension after reading in the ascii table and 
   # defining the columns in cols, the data are in tab (TBD)

note: getting the Gaia EDR3 with high PM (full sky)

SELECT TOP 12001000 gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,
gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_over_error,gaia_source.pm,
gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,
gaia_source.ruwe,gaia_source.phot_g_mean_flux,gaia_source.phot_g_mean_flux_over_error,
gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_flux_over_error,
gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_flux_over_error,
gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.dr2_radial_velocity,
gaia_source.dr2_radial_velocity_error,gaia_source.l,gaia_source.b
FROM gaiaedr3.gaia_source 
WHERE (gaiaedr3.gaia_source.pm>=25)

===

14 June 2022 https://gea.esac.esa.int/archive/ log in 
gaiadr3 from lite catalogue all sources with PM>25

SELECT TOP 30000000 * FROM gaiadr3.gaia_source_lite
WHERE (gaiadr3.gaia_source_lite.pmra*gaiadr3.gaia_source_lite.pmra+gaiadr3.gaia_source_lite.pmdec*gaiadr3.gaia_source_lite.pmdec) > 625

=> gaiadr3_highpm
===

2 Mar 2023 https://gea.esac.esa.int/archive/ log in 

:retrieve all gaia with PM > 2mas/yr Epoch 2016, for matching our source to gaiadr3 ID. 

SELECT TOP 30100100 gaia_source.designation,gaia_source.source_id,gaia_source.ra,
gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,
gaia_source.parallax_over_error,gaia_source.pm,gaia_source.pmra,gaia_source.pmra_error,
gaia_source.pmdec,gaia_source.pmdec_error FROM gaiadr3.gaia_source WHERE (gaiadr3.gaia_source.pm>=10.0)

=> gaiadr3_pmgt2

==== from notebook 1
28 june 23

New match of the new SUSS single source result V2 which was rematched to the previous 
version match in order to keep the better Gaia match from the Xmatch site.   
(TBD file name used)XMMOM_SUSS5.0_Sources_v0_1xGaiaDR3.fits
(2)XMM-OM-SUSS5.0_singlerecs_v2.fits
(3)XMMOM_highpm_probaAB_gt_half.fits

The reason to do it again is that the initial match of the high PM sources came back with 
Gaia source_ids. 
When I added further Gaia DR3 columns to the whole SUSS, I used a match on position, but
some pf the high PM sources were then mismatched to another source. So first we split off 
the high PM sources, and match them using the Gaia source_id. However, we need to choose 
how to do that. We have the downloaded subset of Gaia DR3 with pm > 10, but that has 
not all the required columns. The other way is to match using TAP on the whole 
gaiadr3.gaia_source catalogue.

Using TopCat we split  XMM-OM-SUSS5.0_singlerecs_v2xhighpm_probabGThalf.fits into 
(1) temp_pm_sources.fits (66,557 sources)  
(2) temp_nopm_sources.fits (5,896,514 sources)

match (1) using TAP on gaiadr3.gaia_source : 
====
26 June 2023 NPMK

New version April-June 2023

- new code for SUSS6 produced to add proper motion effects (sent to Simon)
[ a] match2gaiadr3_positions2Epoch
[ b] addEpoch2Source.sh (previously developed for this project)
The code works by sorting the input first for observations done 
within a range of Epochs (typically 0.1 year); then matching the 
selected subset to Gaia DR3 sources with PM > 25 mas/year which 
were placed on the relevant Epoch. Matching is done using STILTS.
 
This is more efficient but less accurate then the method used for this project 
using the Xmatch website, but that requires uploading and downloading the data.
Xmatch advantage is the consideration of how crowded the field is on the 
possibility of mis-match. 

- software to produce the single-source catalogue in two steps, 
   - match sources based on their IAUNAME
   - match results on PM, and plx_over_error being the same
   
The software was changed to derive new statistics, namely chi-squared, skew. 
Also, two fields were added capturing OBSIDS and Epochs of the input. 
The latter two are the same nearby objects unless at the edge of an observation. 
The magnitude is now the median value, while the extreme is the minimum magnitude 
in the input, so maximum flux. 

The single source processing thus produces two files, one with significant PM and 
one where the objects are not reported in Gaia DR3 as having PM>25 mas/yr. Only the 
first set has a Gaia DR3 source_id match. 

- In the single source catalogue there are sources with plx < 0 or with large 
positional uncertainties. Examination of the histograms of these data with separation 
of source and Gaia DR3 match showed requiring plx > 0 is consistent with the 
expected drop off of the number of matches with separation, while the excluded 
set with plx < 0 did not fall off with separation. 
[positional errors? ]

The matched sources with plx>0 are considered to be of higher quality.

==== 
30 Jun 2023 TopCat 

part1: match only pmcat.gaiadr3_source_id from TAP_UPLOAD (temp_pm_sources.fits)

SELECT pmcat.*, refcat.source_id, 
refcat.pm,refcat.parallax_over_error,
refcat.phot_g_mean_flux_over_error,
refcat.phot_g_mean_mag,
refcat.phot_bp_mean_flux_over_error,
refcat.phot_bp_mean_mag,
refcat.phot_rp_mean_flux_over_error,
refcat.phot_rp_mean_mag,
refcat.bp_rp,
refcat.bp_g,
refcat.g_rp,
refcat.radial_velocity,
refcat.radial_velocity_error,
refcat.phot_variable_flag,
refcat.classprob_dsc_combmod_galaxy,
refcat.classprob_dsc_combmod_quasar,
refcat.classprob_dsc_combmod_star,
refcat.non_single_star
FROM gaiadr3.gaia_source AS refcat
RIGHT OUTER JOIN TAP_UPLOAD.t5 AS pmcat
ON refcat.source_id=pmcat.gaiadr3_source_id


This works only if using the source_id in the upload temp_pm_sources.fits. 
Afterwards, match locally with the full input temp_pm_sources using Topcat
XMMOM_SUSS5.0_SourcesV2part1.fits => 66,557 sources

*** NEXT part2: the sources with no PM match. 

unfortunately, for matching this with position error is not possible (I would 
need a local copy of Gaia DR3 for STILTS)

So I did a TopCat match CDSxmatch on position 'best' with radius 0.9 arcsec for 
temp_nopm_sources.fits with gaiadr3; only keep columns like part1 above. 


Further selection is needed: some of the sources in the high PM set (1) are likely 
not good matches. After plotting the histograms it became clear that selecting 
on the parallax_over_error also removed sources with bad positional match:
-- Keep only those with parallax_over_error > 0 which removes 

There remain now some sources with not any Gaia DR3 match 
5,896,514 in temp_nopm_sources.fits, and 4,534,565 matches in XMMOM_SUSS5.0_SourcesV2part2.fits, so 
4,534,565
--------- -
1,361,949 sources do _not_ have a match in GaiaDR3 within 1.5 arcsec.  

Proposed pmflag:
pmflag : 1 for all matches in part1 with parallax_ovser_error > 3 based on position error per object
         2 for matched in part1 otherwise ; also based on position error per object
         3 for part2 with match parallax_over_error > 3 with position within 1.5 arcsec 
         4 for part2 otherwise with position within 1.5 arcsec
         5 for no Gaia DR3 match found within 1.5 arcsec 
         
         
merging  to  XMMOM_SUSS5.0_SourcesV2part1+2   (5,963,071 sources)     
=====
reporting July 5, 2023 

starting with (3) suss_gaia_epic/XMMOM_SUSS5.0_Sources_v2part1+2.fits:

I added a new column 'chi2red' defined as 

median([UVW2_CHISQ/UVW2_NOBS , UVM2_CHISQ/UVM2_NOBS, UVW1_CHISQ/UVW1_NOBS. U_CHISQ/U_NOBS, B_CHISQ/B_NOBS, V_CHISQ/V_NOBS])

This distribution as a log(chi2red peaks at a value of around 0.5. I define therefore a subset with 
possibly variability of interest as chi2red > 2, which has 214,229 rows/sources, that is 4% of all sources. 
I save this as
 (4) suss_gaia_epic/source_v2_chi2redGT2.fits

****** After discussion with Mat, he thinks that we should select on the extreme value of the reduce ChiSq. 

****** On hindsight (July 18) I think we need to select also those with the higher number of observations in anyone of the filters, 
       because doing more (automated) analysis with, e.g., wavelets or power spectrum needs a decent number of points to fit. 

changed definition to chisq/(N-1)

Next I use TopCat match (4) to (5) XMM_SUSS5.0_Sources_v0_1_traininginput_v2.fits using 
IUANAME then omitting from the result  *var*, colours, columns, and leave the merged 
magnitudes per colour. The result is then matched to 
(6) XMMOM_SUSS5.0_Sources_v0_1_classifications.fits  which adds the classifications to the 
result being: (7) source_v2_chi2redGT2_aux_class.fits (214,211 sources)

The result was finally matched by J2000,Epoch 2000 positions (1.0" radius) to SIMBAD using TopCat 
-> (8) source_v2_chi2redGT2_aux_class_SIMBAD.fits  (14,963 rows) 

    #     NOBS(max)
w2  5051  36
m2  5902  38
w1 11336  59
u   7951  67
b   4663  41
v   4226  58

The non-SIMBAD has longer NOBS series

Notice that we did not rerun the classification at this time. The main change is in the calculation
of Chi-Squareds. 

July 18, updated these notes. 

In the past week we selected some sources to investigate how their data supposted the chi-squared indication
of variability. One source showed a clearly decreasing flux over the years, with the data also 
spread fairly wide around the trend. Other sources showed large changes in flux without a clear 
trend. The issue is partly that the data have been taken in a rather random way over time. 

The examination was rather time consuming using TopCat since the source information in (8) was used 
on the original SUSS5.0+Epoch. I wrote a program that prepares a fits file that combines 
all for one source: 
make_for_many_variable(pos=None, poserr=None, iauname=None, recnums=[0,5], chatter=1)



   '''
"""   
# obsolete because does not work (matches whole exposures):

# 11000 last one done -12000-13000 yellow 14000-16000 white now
def make_for_many_variable(pos=None, poserr=None, iauname=None, recnums=[], chatter=1):
    
    given the single source catalogue selection of vaiable sources, this will 
    create a fits file with just the SUSS observations of a given source.
    
    Note that the iauname has a single space in it before the 'J'
    
    position input 'pos' must be as an astropy coordinate object
    
    recnum will pull the source of the 'recnum' record after sorting on IAUNAME
    
    output file name is based on iauname, the class from our classification, and number
    of UV+u records.
    
    For each filter another extension, and add MJD. 
    
    fecit npmk june 2023 
    

    from astropy import coordinates, units as u, constants as c
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table, vstack

    indir = '/Users/data/catalogs/suss_gaia_epic/'
    suss = 'XMM-OM-SUSS5.0ep.fits'
    varfile = 'source_v2_chi2redGT2_aux_class_SIMBAD.fits' # 14,963 rows 
    fsuss = fits.open(indir+suss)
    fvars = fits.open(indir+varfile)
    t = Table(fsuss[1].data)
    v = Table(fvars[1].data)
    v.sort('IAUNAME')
    v.add_index('IAUNAME')
    print (f"the total number of sources is now {len(v)}")

    for recnum in range(recnums[0],recnums[1]):
        
        # what to look for
        if iauname != None: 
           recnum = v['IAUNAME'] == iauname
    
        if pos != None:
           ra = pos.ra.deg
           dec = pos.dec.deg
       # NEED recnum from match to v   
       
        if pos == None and iauname == None:
            if recnum < 0:
                raise IOError("give position pos(astropy.coordinates) or iauname as input")
            
        # new: select records based on OBSIDS, EPOCHS in v => get all sources in image
        rec = v[recnum]
        t.sort('OBSID')
    
        ob = rec['OBSIDS'].split('_')
        ep = rec['EPOCHS'].split('_')  
        # remove duplicates
        epstr = []
        for o9,e9 in zip(ob,ep): epstr.append( f"{o9}{e9[:10]}" )
        uniobep = np.unique(epstr)
        ob = []
        ep = []
        for o9e9 in uniobep:
            ob.append( o9e9[:10] )
            ep.append( np.array(o9e9[10:],dtype=float) )
        # locate records in SUSS    
        qq = []  # list of matches    
        for ob1, ep1 in zip(ob,ep):
             ts = t['OBSID'] == ob1  # boolean list
             tee = np.abs(t['obs_epoch'] - np.float(ep1)) < 3e-5 # bool sub-list
             qq.append(ts & tee)      
        # select records SUSS for that image ; next find the right source
        t1s = [ ]
        for k in np.arange(len(qq)):
            t1s.append(t[qq[k]])
        t2 = vstack(t1s)          
       
        # position desired source: 
        acheck = np.isnan(rec['pm'])   
        if acheck | (rec['pm'] < 43):
            poserr = 1.5/3600.
            if not acheck: poserr = (1.5 + rec['pm']*0.03)/3600.
            # define Epoch 2000 position of source        
            ra  = rec['RAJ2000Ep2000']
            dec = rec['DecJ2000Ep2000']
            pos = coordinates.SkyCoord(ra,dec,unit=(u.deg,u.deg),frame='icrs')
            dra = np.abs(t2['RA']-pos.ra.deg) < poserr
            dde = np.abs(t2['DEC']-pos.dec.deg) < poserr
            t1 = t2[dra & dde]
        else:    
            poserr = 1.8/3600. 
            # need to precess to the epochs of observation
            ra  =  rec['RAJ2000Ep2000']
            dec =  rec['DecJ2000Ep2000']
            pmra = rec['pm']  *0.5 # need ['pmra'] in v catalogue
            pmde = rec['pm']  *0.5 # need ['pmde']
            ra_epoch=[]
            dec_epoch=[]
            t1s = []
            for obsep in t2['obs_epoch']:  # for each star
                ra9, dec9 = apply_proper_motion(ra,dec,2000.0,pmra,pmde, obsep )
                dra = np.abs(t2['RA'] -ra9) < poserr
                dde = np.abs(t2['DEC']-dec9) < poserr
                if (dra&dde).sum() > 0:
                    t1s.append(t2[dra & dde])
            t1 = vstack(t1s)    

#  write output file
               
        nobs = rec['UVW2_NOBS']+rec['UVM2_NOBS']+rec['UVW1_NOBS']+rec['U_NOBS']
        class1 = rec['class_prediction']
        name = rec['IAUNAME'].replace(" ","_")   # replace space with underscore
        outfilename=indir+f"../var/var_{name}_{class1}_{int(nobs)}."
    
        if chatter > 0: print (recnum, outfilename)
    
        hdu = fits.PrimaryHDU()
        hdulist = fits.HDUList(hdu)
        hdulist[0].header['COMMENT']="observations of a single source from XMM-OM SUSS"
    
        # create bintable HDU from records 
        hdu1 = fits.BinTableHDU( t1, name=rec['IAUNAME'])
        hdu1.header['iauname'] = rec['IAUNAME']
        hdu1.header['RA']      = (rec['RAJ2000Ep2000'],'J2000ep2000')
        hdu1.header['Dec']     = (rec['DecJ2000Ep2000'],'J2000ep2000')
        hdu1.header['class']   = (rec['class_prediction'],'predicted class')
        if np.isfinite(rec['chi2red']): hdu1.header['chi2red'] = (rec['chi2red'],'reduced chi-squared, median from six filters')
        hdu1.header['main_id'] = (rec['main_id'],'SIMBAD')
        hdu1.header['maintype'] = (rec['main_type'],'SIMBAD')
        if np.isfinite(rec['redshift']): hdu1.header['redshift'] = rec['redshift']
        hdu1.header['sp_type'] = rec['sp_type']
        hdu1.header['w2_nobs'] = rec['UVW2_NOBS']
        hdu1.header['m2_nobs'] = rec['UVM2_NOBS']
        hdu1.header['w1_nobs'] = rec['UVW1_NOBS']
        hdu1.header['u_nobs']  = rec['U_NOBS']
        if np.isfinite(rec['UVW2_AB_MAG']): hdu1.header['w2_abmag'] = (rec['UVW2_AB_MAG'],'median mag')
        if np.isfinite(rec['UVM2_AB_MAG']): hdu1.header['m2_abmag'] = (rec['UVM2_AB_MAG'],'median mag')
        if np.isfinite(rec['UVW1_AB_MAG']): hdu1.header['w1_abmag'] = (rec['UVW1_AB_MAG'],'median mag') 
        if np.isfinite(rec['U_AB_MAG']): hdu1.header['u_abmag'] = (rec['U_AB_MAG'],'median mag')
        if np.isfinite(rec['Gmag_gaia']): hdu1.header['Gmag']    = (rec['Gmag_gaia'],'Gaia DR3 Gmag')
        if np.isfinite(rec['pm']): hdu1.header['PM'] = (rec['pm'], 'Gaia DR3 proper motion mas/yr')
        hdulist.append(hdu1)
    
        # create extension with variable single source record ? SOURCE_INFO
        hdu2 = fits.BinTableHDU(Table(rec),name="varinfo")
        hdulist.append(hdu2)
        # create extension for each filter: e.g., LC_UU (filter names W2,M2,W1,UU,VV,BB)
    
        hdulist.writeto(outfilename+'fits', overwrite=True)
        
        
   December 5, 2023
   
I matched the subset which has gaiadr3_source_id values to the Gaia DR3 lite catalogue 
in order to retrieve the lost columns for Gaia magnitude & colours Gmag, G-BP, G-RP, BP-RP
Previously, the rows of V3/XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2_srcnum.fits
with gaiadr3_source_id > 0 were selected for the input file (UPLOAD.t1):
        
SELECT * 
   FROM TAP_UPLOAD.t1 INNER JOIN gaiadr3.gaia_source_lite 
   ON TAP_UPLOAD.t1.gaiadr3_source_id = gaiadr3.gaia_source_lite.source_id    
   
This used the TAP Query interface in TopCat with ivo://esavo/gaia/tap  on table 
gaiadr3.gaiasource_lite.    

The other rows, i.e., with no gaiadr3 source id were saved. 66557 sources had match with PM
correction. Now the other 5 896 493 sources. Match to Gaia DR3 2016 using positions (RA,DEC) 
from SUSS using the CDS interface I/355/gaiadr3

three parts
gaia_tmp3b  source with PM  - match by gaiadr3 source id
gaia_tmp4b  source in Gaia, small PM  - match by position (ra,dec) best within 1.0"
gaia_tmp5b  source not in Gaia 

after merging and renaming the files: 

-rw-r--r--  1 kuin  staff    362989440 Dec  5 18:45 XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b_part1.fits
66,577 sources with PM > 25 mas/yr
-rw-r--r--  1 kuin  staff  22573393920 Dec  5 18:51 XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b_part2.fits
4,149,511 sources with PM < 25 mas/yr and no match to Gaia DR3 with angDist < 1 arcsec
-rw-r--r--  1 kuin  staff   9402304320 Dec  5 21:36 XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b_part3.fits
-rw-r--r--  1 kuin  staff  32516568000 Dec  6 00:09 XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b.fits

Problems found with SUSS5.0 

1. some srcnum contain two entries for the same sources, e.g., 2429739 what is that? 
   single image , summed image ? the magnitudes are different 
   origin is the SUSS5.0 
   
   wrote script to make sure in light curve just a single time/mag per OBSID.
   
2. rudy-5 or summed image issue ? 
   in:XMMOMV2_J004210.09411529.3_qso_uvw1_n72.fits OBSID:0202230301 srcnum:112546 
   and many more in that file
   
   no, also here two clusters in RA,DEC -- both variable, but at different mean mag.
   
3. angular distance for match SUSS5 to Gaia DR3

2024-01-23 

Need for matching SUSS within more than 1" to Gaia DR3. Mat pointed out that 
there is a possibility that for the bright stars the OM positions are not as accurate as
1 arcsec. So I match the part3 file mentioned above to Gaia DR3 within 3.0".  It looks 
like there are good matches up to 2.0 arcsec, after which we may pick up some bad ones. 

XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b_part3A ... .fits has 476,795 sources 
that match with a Gaia source to within 3 arcsec angDist 
XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b_part3B ... .fits has 1,270,007 
sources that can find no match to Gaia within 3". 


2024-02-06 

Classification rework. 

Selected a new set of 30 input parameters. Better result. Mat came up with an expression 
for differencing Gaia G mag to ground-based magnituds with sensitivity to the spatial 
extent of a source. It starts with a filter. 

gaia_extended =  ((BII>10.0)||(BII<-10.0))&&(gmag_M<19.0)&&!(gaiadr3_ra<500.0) ? 
  (20.0-gmag_M) : ((gmag_M > 0.0) ? (gaia_Gmag-gmag_M) : (gaia_Gmag-V_AB_MAG))

where BII is the galactic lattitude.  


The claxbol software (classify_new.py) when running the optimisation needs to be 
interrupted (CNTL-C) which leads to incomplete classifications file. 
Also, the API has been updated requiring replacement of product() with prod(). 

Priors: The selection for the priors has been done by crossmatching with SIMBAD, and then
determining the proportion of stars, and non-stars, being 69% and 31%. SIMBAD selection 
of QSOs and galaxies then gives 6% for QSOs and 25% for galaxies. 
In version 0.1 of the classification, we used as priors [0.65,10,25]% and got as a
result proportions for stars, QSO, galaxies = [0.78,0.04,0.18]. 
Alternatively, using the proportions of QSO and galaxies from SDSS16 which are the 
training set, we would get [0.69, 0.13, 0.18]. 
Since there is uncertainty in the proportions of stars, QSO, and galaxies, a range of 
values for the priors have been explored. 

2024-02-22 Do a rerun with only quality_flag == 1 data for variable sources catalogue

make_file_variable (onlyqualzero=True)

2024-05-21 and later

Need to redo the earlier classification with gaia_Gmag - WISE-W1 and gaia_extended, 
using gaiadr3_source_id>0  instead of gaiadr3_ra< 500. 

Also create 

  Gaia_G_WISE_W1=((BII>10.0)||(BII<-10.0))&&
  (WI_W1mag<16.0)&&
  !(gaiadr3_ra>0)?(20.0-WI_W1mag):(Gaia_Gmag<15.6)&&
  !(WI_W1mag>0.0)?(Gaia_Gmag-17.1):(gaia_Gmag-WI_W1mag)
  
    Gaia_G_WISE_W1=((BII>10.0)||(BII<-10.0))&&(WI_W1mag<16.0)&&!(gaiadr3_ra>0)?(20.0-WI_W1mag):(Gaia_Gmag<15.6)&&!(WI_W1mag>0.0)?(Gaia_Gmag-17.1):(gaia_Gmag-WI_W1mag)

and thus, 

  gaia_extended =  ((BII>10.0)||(BII<-10.0))&&
  (gmag_M<19.0)&&
  !(gaiadr3_ra>0)?(20.0-gmag_M):((gmag_M>0.0)?(gaia_Gmag-gmag_M):(gaia_Gmag-V_AB_MAG))

   gaia_extended=  ((BII>10.0)||(BII<-10.0))&&(gmag_M<19.0)&&!(gaiadr3_ra>0)?(20.0-gmag_M):((gmag_M>0.0)?(gaia_Gmag-gmag_M):(gaia_Gmag-V_AB_MAG))

where BII is the galactic lattitude.  gaiadr3_ra is being used, since there are some 
Gaia sources with an source_is but without any photometry. 

Starting now once more from 
   XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b.fits (Jan 24, 2024)
   
create new input for classification in htranin directory
    XMMOM_SUSS5.0_Sources_4claxbol_v6.csv
    XMMOM_SUSS5.0_Sources_4claxbol_v6.in   
classification is found in classificationcode/
-rw-r--r--  1 kuin  staff         288  5 Jun 14:50 classification_set_v6.metrics
-rw-r--r--  1 kuin  staff  1592182757  5 Jun 14:50 classification_set_v6.csv

merged clssification_set_v6 [IAUNAME] to subset of XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b.fits
(not VAR3, OBSID, SRCNAM, SRCDIST, SIGNIF, FLUX), and thus retains some basic aux/gaia data.
this gave a catalogue in directory v3:
====>
   XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits
The classification report based on the input classification sets remains virtually 
the same as before, but due to including upper limits which are needed since 
the Gaia Gmag only goes down to ~19.5 magnitudes, and making that update for 
sources outside the Galactic plane, we find now some previously classified as stars 
are being better represented as QSO or galaxy class. 
   
   star          QSO        galaxy
A : NpC0=3809450, NpC1=516515, NpC2=1637085
   previously (no upper limits)
B : NpC0=4136491, NpC1=370822, NpC2=1455737  

The ratio's of change are star:0.92, QSO:1.39, galaxy:1.12
and the relative fractions are stars:0.639, QSO:0.087, galaxy:0.274

And this show a dramatic increase in QSO's. We therefore now test against 
external catalogues as before using the SIMBAD object types. 
  in TopCat x SIMBAD and the catalogue match ra,dec J2000, ep2000, delta_r=2.5 arcsec. 
  -> 239,227 rows
  write to 
  /Volumes/DATA11/data/catalogs/suss_gaia_epic/v3/XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6xSIMBAD.fits
 
If we treat the SIMBAD main_type=Star,QSO,Galaxy as truth, we get the following
matrix

true->	        	star	QSO		galaxy	retrieved 
predicted star		83003	268		4173	89.3%
predicted QSO		6282	6404	1076	94.4%
predicted galaxy	3479	109		12817	70.9%
TruePos.Rate		94.8%	46.5%	78%

The low t.p.r. for QSOs is due to the large number of misclassified stars.
Indeed, 2118 have a parallax in SIMBAD/Gaia. Also, we took the training sets 
from the SDSS which is partly in SIMBAD, so this is a bit biased, but, for 
example, for true star, predicted QSO, and not in training set, we get 6278
sources, only 7 less than in the table above. 
We need to discriminate between stars and QSO better. 
The outlier histogram for this set as compared to the full catalogue shows 
 a different shape. Whereas the full catalogue has outlier peak around 9, 
 this set is flattish between 8 and 19, so there is a poor match to the class.
 Restricting outlier to < 15 leaves still 4626 misclassified sources. 

"""

def make_file_variable( band='uvw1', minnumber=10, maxnumber=1000, 
         chi2_red_min=5., 
         inputdir='/Volumes/DATA11/data/catalogs/suss_gaia_epic/v3/',
            ##inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6xSIMBAD.fits.gz",
         inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits",
         selectrecs=[0,60000000.], plotit=False, onlyqualzero=False,
         min_srcdist=6.0,
         bad_obsidfile="remove_obsidfile_here.txt",
         chatter=2):
    """
    
    given the single source catalogue (SSCenh) which includes NOBS and 
    Chi-squared, this will create a fits file with just the SUSS 
    observations of a given source.
    
    Note that the iauname has a single space in it before the 'J'
    
    minnumber requires that <band>_NOBS >= minnumber for processing
    
    #position input 'pos' must be as an astropy coordinate object
    recnum will pull the source of the 'recnum' record after selecting 
    
    output file name is based on 
    IAUNAME2, <band>, <band>_NOBS, <band>_CHI2, Class 
    iauname2 is based on position at J2000 epoch 2000, while the 
    iauname is that of one record in SSCenh.  
    
    output: header pulls data from SSCenh, and data records from SUSS.
    For each filter another extension, and add an MJD column. 
    
    fecit npmk June-October 2023
    dec 11, 2023: updated lc filename (+ in dec), input catalogue, output directory
    
    Jan 2, 2024: start adding screening for double OBSID entries: keep only the one with 
       the longest exposure time/smallest error
       
      Variables 
    for input combine  
    1.XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2...   
    2.classification_set_v4.csv
    match(XMMOMSUSS5,classificationv4) -> 
    3.XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v4.fits
    
    Feb 22: add option to restrict to quality == 0 (good data)
    June/July 
       change criterium chisq red to chisq/(Nobs-1) > 5(chi2_red_min)
       only count using point source data 
       modify plot 
    
    Aug 2024: added reading in for the bad obsids 
    
    """

    from astropy import coordinates, units as u, constants as c
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table, vstack
    import os

    brief = inputdir+"../var_lc/brief.txt" # for adding filename, ra, dec, IAUNAME
    brieffh = open(brief,"a")
    brieffh.write(f"IAUNAME  RA Dec  filename\n")

    indir = inputdir 
    suss =  '../../XMM-OM-SUSS5.0ep.fits'
    varfile = inputcat # "XMM-OM-SUSS5.0ep_singlerecs_v2_srcnum_aux_stage2_v2b.fits"
    fsuss = fits.open(indir+suss)
    print ("SUSS opened")
    # 
    t = Table(fsuss[1].data)   # this is the SUSS5.0 file + obs_epoch
    thead = t.colnames
    t.sort('SRCNUM')
    t.add_index('SRCNUM')
    print ("SUSS sorted by SRCNUM")
    #
    bad_obsids = []    
    if os.access(indir+bad_obsidfile,os.F_OK):
       with open(indir+bad_obsidfile) as gf:
           xx = gf.readlines()
       for obs in xx:
          bad_obsids.append(obs.split("\n")[0])
       print (f"Read in {len(bad_obsids)} bad obsids")    

    fvars = fits.open(indir+varfile)
    print (f"Single Source Catalogue {varfile} opened")
    #
    vv = fvars[1].data[ np.isfinite(fvars[1].data[f"{band.upper()}_CHISQ"]) ] 
    print (f"restricted to finite chi2 ")
    v = Table(vv)   # this is the list of variable sources in band
    print (f"cast SSC into Table at start length = {len(v)}")
    fvars.close()
    
    filt = f"{band.upper()}_NOBS"
    print (f"using {filt} for selection")
    
    #v.sort(filt)
    #v.add_index(filt)
    #print (f"sort SSC by filter")
    chi2red = v[f'{band.upper()}_CHISQ']/(v[filt]-1)  #20240701 changed denom from N to N-1
    
    # keep only records with quality=0 in the single source catalogue (subset)
    if onlyqualzero:
        sscqual0  = v[f"{band.upper()}_QUALITY_FLAG"] == 0
        #sussqual0 = t[f"{band.upper()}_QUALITY_FLAG"] == 0
    else:
        sscqual0  = True   
        #sussqual0 = True
    
    print (f"\t all sources = {len(v)} ")
    qb = (v[filt] >= minnumber) & sscqual0   
    allsrcgt10 = qb.sum()
    print (f"\tN > {minnumber} quality=0, number of SSC sources = {qb.sum()}")

    qb = (v[filt] >= minnumber) & (chi2red > chi2_red_min) & sscqual0
    print (f"\tchisq reduced > {chi2_red_min} and number of source = {qb.sum()}")
    ch2rgt6 = qb.sum()
    
    v = v[qb]
    recnums = [0,len(v)]
    
    # ensure the top limit is constrained 
    #selectrecs[1] = int(np.min([selectrecs[1],len(v)]))
    #v = v[selectrecs[0]:selectrecs[1]]
    selectrecs[1] = len(v)
    ngood = len(v)
    print (f"the remaining number of sources has been restricted to {len(v)}")
    recnums[1] = np.min([recnums[1], len(v)])
           
    """
    June 2024
    So far, we have selected transients based on Nobs and Chi2-reduced, but not selected 
    on whether the source was detected in the OM processing as extended or pointlike. 
    That is important because extended sources are measured using isophotal photometry, 
    and that depends on the exposure time, where the longer exposures exhibit the fainter 
    outer regions of galaxies. So in the next part we collect the relevant observations 
    for each source, but then deselect those which have less than 10 pointlike 
    observations. 
    
    25 July 2024 forced chi2 to use only qual=0 data 
    band, all $N(qual=0) > 10$ & 4N(qual=0)>10$ + $\chi^2_red(qual=0) > 6$ \\
    uvw2 &   4066 &   48 \\
    uvm2 &   3533 &   24 \\
    uvw1 &  23095 &  430 \\
    u    &  13549 &  275 \\
    b    &  14399 &  163 \\
    v    &  15970 &  168 \\
    all  &  74612 & 1108 \\
    
    29Aug 2024: using XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fit... and selecting 
    UVW2_CHISQ/UVW2_NOBS>5&UVW2_SRCDIST>6&UVW2_QUALITY_FLAG==0&UVW2_NOBS>5|UVM2_CHISQ/UVM2_NOBS>5&UVM2_SRCDIST>6&UVM2_QUALITY_FLAG==0&UVM2_NOBS>5|UVW1_CHISQ/UVW1_NOBS>5&UVW1_SRCDIST>6&UVW1_QUALITY_FLAG==0&UVW1_NOBS>5|U_CHISQ/U_NOBS>5&U_SRCDIST>6&U_QUALITY_FLAG==0&U_NOBS>5|B_CHISQ/B_NOBS>5&B_SRCDIST>6&B_QUALITY_FLAG==0&B_NOBS>5|V_CHISQ/V_NOBS>5&V_SRCDIST>6&V_QUALITY_FLAG==0&V_NOBS>5
    I obtain 10,632 sources
    But.. SimbadxSUSS5_variable_sources.fits (V8) only has 3448 records. Writing to 
    SimbadxSUSS5_variable_sources_alternate.fits in var_lc
    My code recalculates Chi2red and then the number of sources with chi2red>5 drops to 
    3448 sources
    
    
    """
    flog = open(inputdir+f"{band}_stats","w")
    flog.write(f"{band} allsrc>10={allsrcgt10} and ch2r>6={ch2rgt6}\n")
    flog.write(f"srcn  band  n_zero_qual  counter chi2red \n")
    counter = 0
    
    for recnum in range(recnums[0],recnums[1]): # process each variable source of <band> 
        print (f"\n\tStarting search for records in  recnum= {recnum} out of {recnums[0]},{recnums[1]}")

        # select records based on SRCNUM in v => get all sources in image
        
        rec = v[recnum]
        srcn1 = rec['SRCNUMS'].split('_')
        epoc1 = rec['EPOCHS'].split('_')
        vobsi1 = rec['OBSIDS'].split('_')

        # remove duplicates (obsid,epoch)
        epstr = []
        for o9 in srcn1: epstr.append( f"{o9}" )
        uniobep = np.unique(epstr)
        
        srcn = []
        for s9 in uniobep:
            if s9 == "":
               continue
            else:   
               srcn.append( s9 )
        if len(srcn) > 0:       
            print (f"the source id numbers are: {srcn}")   
            srcn = np.array(srcn,dtype=np.int64)    
        else:
            print ("WARNING:  There are no valid sources found for recnum={recnum}.")
            continue
        
        # select all (that is SRCNUM) records in original SUSS for that source (includes also other sources)
        
        ta = t[ t['SRCNUM'] == srcn[0] ] 
        ta2 = ta[ np.isfinite(ta [f'{band.upper()}_AB_MAG']) ]
        # check the OBSID is in the vobsi1 list from rec   22 aug 2024
        qobsid = False
        trecnums = []
        for vob in vobsi1:
            for tk9, trec9 in enumerate(ta2):
                if trec9['OBSID'] == vob:
                    trecnums.append(tk9)
                    qobsid = True
        if len(trecnums) == 0:
            continue            
        t2 = ta2[ trecnums ]
        
        if len(srcn) > 1:    
            t1s = [t2]  # list of Table matches    
            for s9 in srcn[1:]:
                 ts = t[t['SRCNUM'] == s9]   
                 print (f"for SRCNUM={s9}, len = {len(ts)}")  # <---- Table() it ?
                 t1s.append( ts ) 
            t2 = vstack(t1s)          # SUSS exposures
        else:
            print (f"only one srcnum with table of lenght {len(t2)}")
        
        # keep cols from SUSS ['obs_epoch',f"{band}_*" ? use stilts .. 
        n1 = ['obs_epoch', 'IAUNAME', 'N_SUMMARY', 'OBSID', 'SRCNUM', 'RA',\
           'DEC', 'RA_HMS', 'DEC_DMS', 'POSERR', 'LII', 'BII', 'N_OBSID']
        n2 =['_SRCDIST','_SIGNIF','_RATE','_RATE_ERR','_AB_FLUX','_AB_FLUX_ERR',\
        '_AB_MAG','_AB_MAG_ERR','_VEGA_MAG','_VEGA_MAG_ERR','_MAJOR_AXIS',\
        '_MINOR_AXIS','_POSANG','_QUALITY_FLAG','_QUALITY_FLAG_ST',\
        '_EXTENDED_FLAG','_SKY_IMAGE'] 
        for nx in n2:
           n1.append(f"{band.upper()}{nx}")

    # 8-8-2024 removed sources with SRCDIST < 6"
        if np.max(t2[f"{band.upper()}_SRCDIST"]) <= min_srcdist :
            print (f"skipping source with srcdist < {min_srcdist}")
            continue 
           
    # 2024-01-02 scan for duplicate OBSIDs and keep best  ; also remove bad obsids
        if chatter > 0: print (f"1938 calling fix_multiple_obsids_in_SUSS t2")
        #t2 = fix_multiple_obsids_in_SUSS5(t2,band,chatter=chatter)  
        t2 = single_obsid_for_srcnum(t2, band,bad_obsids=bad_obsids, chatter=chatter)
        
        # filter point-like sources and quality=0
        qxt = t2[f"{band.upper()}_EXTENDED_FLAG"] == 0
        # 2024-07-24     
        zero_qual_pnt = qxt & (t2[f"{band.upper()}_QUALITY_FLAG"]==0)
        n_zero_qual = zero_qual_pnt.sum()
             
        # recalculate chi-squared-reduced using n_zero_qual filter
        chi2r = 100
        from cats.single_source_cat_v2 import stats
        if (n_zero_qual >= minnumber):
           t3 = t2[zero_qual_pnt]
           medmag = np.median(t3[f"{band.upper()}_AB_MAG"])
           mag = t3[f"{band.upper()}_AB_MAG"] 
           err = t3[f"{band.upper()}_AB_MAG_ERR"]
           statn, statmedmag, ncChisq, statsd, statsk, statV3  = stats(mag,err=err,syserr=0.02) 
           chi2r = ncChisq/(statn-1)
           print ("    magnitudes t3= ",t3[f"{band.upper()}_AB_MAG"],"   median= ",medmag)
           print (f"1952 =>  re-checking if chi2red < 6: {chi2r}")
           print (f"{srcn}  {band}  {n_zero_qual}  {counter}  {chi2r}\n")
        
        flog.write(f"{srcn}  {band}  {n_zero_qual}  {counter}  {chi2r}\n")
        if (n_zero_qual < minnumber) | (chi2r < chi2_red_min): 
            print (f"skipping {srcn} {band} {n_zero_qual} {chi2r}\n")
            continue
        else:
            counter+=1
            
#  write output file
               
        #nobs = rec['UVW2_NOBS']+rec['UVM2_NOBS']+rec['UVW1_NOBS']+rec['U_NOBS']
        nobs = rec[f"{band.upper()}_NOBS"]
        cl1 = np.array(['star','qso','galaxy'])
        class1 = cl1[rec['predicted_class']]
        
        # create IAUname2 for Epoch 2000, ICRS
        ra2  = rec[RA2000] 
        dec2 = rec[DE2000]
        newra,newdec = degrees2sexagesimal(ra2,dec2,as_string=True)
        newra = newra.replace(":","")
        newdec = newdec.replace(":","")
        decsign = '+'
        if dec2 < 0 :
            decsign = ''  # error in dec deg when less than 10- no leading zero
        if rec['IAUNAME'] == "":
           name = f"XMMOM_J{newra:8.1f}{newdec:+7}"  # 2024-02-21 added formatting
           print (f"No IAUNAME found for this row {rec} ")
           brieffh.close()
           raise RuntimeError(f"ERROR - no column IAUNAME in input file found for\n{rec}\ntried to compose name = {name}")
        else:   
           name = rec['IAUNAME'].replace(" ","_")  
           # replace space with underscore 2024-02-22 back to iauname to have traceability
        outfilename=indir+f"../var_lc/{name}_{class1}_{band.lower()}_n{int(nobs)}"
        
        brieffh.write(f"{name:23s}   {ra2:10.5f} {dec2:10.5f}   {outfilename}.fits\n")
        print(f"{name:23s}   {ra2:10.5f} {dec2:10.5f}   {outfilename}.fits\n")
       
    
        if chatter > 0: print ('working  on .... ',outfilename)
    
        hdu = fits.PrimaryHDU()
        hdulist = fits.HDUList(hdu)
        hdulist[0].header['COMMENT']=f"{band} observations of a single source from XMM-OM SUSS"
    
        # create bintable HDU from records 
        print (f"t2 type - { type(t2) }  - length - {len(t2)}")
        extname = name
        hdu1 = fits.BinTableHDU( Table(t2)[n1] ,name="lightcurve")
        hdu1.header['EXTNAME'] = name
        
        hdu1.header['iauname'] = rec['IAUNAME']
        hdu1.header['RA']      = (rec[RA2000],'J2000ep/ICRS')
        hdu1.header['Dec']     = (rec[DE2000],'J2000ep/ICRS')
        hdu1.header['class']   = (rec['predicted_class'],'predicted class')
        #if np.isfinite(rec['chi2red']): hdu1.header['chi2red'] = (rec['chi2red'],'reduced chi-squared, median from six filters')
        #hdu1.header['main_id'] = (rec['main_id'],'SIMBAD')
        #hdu1.header['maintype'] = (rec['main_type'],'SIMBAD')
        #if np.isfinite(rec['redshift']): hdu1.header['redshift'] = rec['redshift']
        #hdu1.header['sp_type'] = rec['sp_type']
        hdu1.header['w2_nobs'] = rec['UVW2_NOBS']
        hdu1.header['m2_nobs'] = rec['UVM2_NOBS']
        hdu1.header['w1_nobs'] = rec['UVW1_NOBS']
        hdu1.header['u_nobs']  = rec['U_NOBS']
        hdu1.header['b_nobs']  = rec['B_NOBS']
        hdu1.header['v_nobs']  = rec['V_NOBS']
        if np.isfinite(rec['UVW2_AB_MAG']): hdu1.header['w2_abmag'] = (rec['UVW2_AB_MAG'],'median mag')
        if np.isfinite(rec['UVM2_AB_MAG']): hdu1.header['m2_abmag'] = (rec['UVM2_AB_MAG'],'median mag')
        if np.isfinite(rec['UVW1_AB_MAG']): hdu1.header['w1_abmag'] = (rec['UVW1_AB_MAG'],'median mag') 
        if np.isfinite(rec['U_AB_MAG']): hdu1.header['u_abmag'] = (rec['U_AB_MAG'],'median mag')
        #if np.isfinite(rec['Gmag_gaia']): hdu1.header['Gmag']    = (rec['Gmag_gaia'],'Gaia DR3 Gmag')
        #if np.isfinite(rec['pm']): hdu1.header['PM'] = (rec['pm'], 'Gaia DR3 proper motion mas/yr')
        
        hdulist.append(hdu1)
    
        # create extension with variable single source record ? SOURCE_INFO
        hdu2 = fits.BinTableHDU(Table(rec),name="varinfo")
        hdulist.append(hdu2)
        
        # create extension for each filter: e.g., LC_UU (filter names W2,M2,W1,UU,VV,BB)
        if chatter > 0: print (f"saving {outfilename}.fits")
        hdulist.writeto(outfilename+'.fits', overwrite=True)  
        
        if plotit:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(t2['obs_epoch'],t2[f"{band.upper()}_AB_MAG"] )
    flog.close()
    brieffh.close()
         

def single_obsid_for_srcnum(t, filt, bad_obsids=[],chatter=0):
    """
    In SUSS5 for each OBSIDs all exposures should have been summed into one row. 
    Unfortunately, for some sources there are multiple rows. 
    Assuming that one row is the summed data, the other not, and discriminating 
    on <filt>_AB_MAG_ERR we keep the row with the lowest mag error. 
    
    Requirement: input Table is for a single SRCNUM and FILTER filt. 
    
    2024-08-14: retain only the good obsids
    
    """    
    if chatter > 2: print (f"only retaining unique obsids for table filter={filt}")
    from astropy.table import table
    
    obsids = np.unique(t['OBSID'])
    ind_obsid = []
    
    for k in range(len(obsids)):
        a = obsids[k]
 
        q = t['OBSID'] == a      # q is/are the index/indices in t 
        m = len(t['OBSID'][q])
        
        if a in bad_obsids:   # skip the bad obsid list
            print (f"skipping obsid {a}")
            m = 0         
            # not adding this record
        elif m == 1:
            iq = np.arange(len(t))[q]
            ind_obsid.append(iq[0])
        elif m > 1:
            tt = t[q]
            errors = tt[f'{filt.upper()}_AB_MAG_ERR']
            #print (20*"= = ",k,f"\n{tt}\n{errors}\n+++++++")
            qs = np.min(errors) == errors # true for index
            iq = np.arange(len(q))[q]
            if len(iq[qs]) == 0:   # failsafe : select any 
               ind_obsid.append(iq[0])
            else:
               if chatter > 2: 
                   print (iq)
                   print (iq[qs])
               ind_obsid.append(iq[qs][0])
        else:
           raise RuntimeError(f"length of OBSIDs is not allowed to be {m}")
    if chatter > 2: print (f" -> {ind_obsid}") 
              
    tout = t[ind_obsid]
    return tout
    

def remove_problematic_extended_and_not_variables(extendedProblemFile):
    """
    The OM may have observed sources and sometimes classified them as extended, sometimes
    as not. This may affect for example compact galaxies depending on the exposure time
    or filter. These need to be identified and filtered out of the db of variable sources. 
    
    
    """
    import os
    from astropy.io import fits
    from astropy.table import table
    
    files = os.listdir('/Volumes/DATA11/data/catalogs/suss_gaia_epic/var_lc/')

    with open(extendedProblemFile,'w') as pf: 
        for file in files:
            band = file.split('_')
            if band[0] != "XMMOMV2":
                continue
            f = fits.open(file)
            t = Table(f[1].data)
            f.close()
             
            ext = t[f"{band[3].upper()}_EXTENDED_FLAG"]
            n_extended = ext.sum()
            if (len(ext) == n_extended) or (n_extended == 0) :
                pf.write(f"{file:30s} {len(ext):04i} {n_extended:04i}     OK\n") 
            else:
                pf.write(f"{file:30s} {len(ext):04i} {n_extended:0.4} NOT_OK\n")
                
                
def plot_lc(lcfile,simbad=False,minpnts=10,noplot=False):
    """
    Make light curve plots of the variables found using the var_lc data files 
    
    Example of use:
    
    cd suss_gaia_epic/var_lc
    ls -1 *.fits > allfiles.txt
    vi allfiles.txt <= add first row "name"
    t = Table.read('allfiles.txt',format='ascii')
    from cats import xmm2athena as xmm2
    for ff in t['name'][0:10]:
         xmm2.plot_lc(ff) 

    """
    import os
    from astropy.io import fits
    from astropy.table import Table
    import matplotlib.pyplot as plt   
    from astroquery.skyview import SkyView 
    import astropy.coordinates as coord
    from astropy import units
    from matplotlib.ticker import MaxNLocator
    
    fn = lcfile
    #if mixed: print (f"plot_lc processing {fn}")
    band = lcfile.split('_')
    with fits.open(lcfile) as f:
        t = Table(f[1].data)
        ext = t[f"{band[3].upper()}_EXTENDED_FLAG"]
        colr = {'UVW2':'purple','UVM2':'blue','UVW1':'deepskyblue','U':'g','B':'y','V':'r'}[f"{band[3].upper()}"]
        q1 = ext == 1
        q0 = ext == 0
        npoint = q0.sum()
        
        if npoint < minpnts:
           print (f"number of points = {npoint}: not enough points for lc plot in {lcfile}")
           return
        
        mark = 'd' 
        if ext[0] == 1: mark = 'x'
        x = t['obs_epoch']
        y = t[f"{band[3].upper()}_AB_MAG"]
        e = t[f"{band[3].upper()}_AB_MAG_ERR"]
        # approximate position for getting DSS plot
        ra = t["RA"].mean()
        dec = t["DEC"].mean()
        radius = 3.0*units.arcmin
        
        ingaia = f['VARINFO'].data['gaiadr3_ra'] > -2.  
        
        #lc_class = classify_light_curve(x, y, e, ax3)
        pa = f[2].data
        sed = make_sed(pa)
        """
        ax1: SED
        the labels on the SED x-axis are crowding each other. so plot log10 wavelength, 
        or slant the axis labels (or both)
        """
 
        
        fig = plt.figure()
        ax  = fig.add_axes([0.13,0.55,0.62,0.35])  # light curve in a filter
        ax1 = fig.add_axes([0.1,0.10,0.42,0.35])  # sed
        ax2 = fig.add_axes([0.75,0.45,0.24,0.40],frameon=False) # text
        #ax3 = fig.add_axes([0.55,0.10,0.40,0.35])  # for PSD plot

        #ax3.tick_params(top=False,bottom=False,left=False,right=True,labelleft=False,labelbottom=True)
        
        if q0.sum() > 0:
              ax.errorbar(x[q0],y[q0],e[q0],fmt=mark,color=colr,label=f"{band[3].upper()}")
        if q1.sum() > 0:    
              ax.errorbar(x[q1],y[q1],e[q1],fmt=mark,color=colr,label=f"{band[3].upper()}-extended", alpha=0.25)
        #ax.legend(fontsize=6)
        ax.set_title(fn)
        ax.set_xlabel('obs epoch (yr)')
        ax.set_ylabel(f'OM {band[3].upper()} (AB mag)')
        dy = 0.1*np.abs(np.max(y)-np.min(y))
        ax.set_ylim(np.max(y)+dy,np.min(y)-dy)
 
        lc_class, class2, FAP, period, minmax, ierup, idimming = classify_light_curve(x, y, e, fig, noplot)
        
        if sed != None:
            if not noplot: 
               t1 = int(np.min(sed[0]))
               t2 = int(np.max(sed[0])+1)
               x_ax1 = ax1.axes.get_xaxis()
               x_ax1.set_major_locator(MaxNLocator(integer=True))
               ax1.plot(sed[0],sed[1],ls="",marker=mark,color=colr)
               #ax1.set_xlim(t1,t2)
            ax1.set_ylabel('SED (AB mag)')
            ax1.set_xlabel('log10(wavelength in nm)')
            ax1.invert_yaxis()
            #ax1.semilogx()
            
        if not noplot: ax2.plot([0,1],[0,1],'w')
        #ax2.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        #ax2.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
        ax2.axis("off")
        
        cl0 = {0:'star',1:'QSO',2:'galaxy'}
        cl1 = pa['predicted_class'][0]
        #print ("cl0=",cl0,"cl1=",cl1)
        cl = cl0[cl1]
        s9 = ""
        if "main_type" in t.colnames:
           s9 = f"SIMBAD main type = {t['main_type']} "
        if cl1 == 0 :  #star
            s10 = ""    
            s11 = "" 
            if ingaia:
                s10 = " in GaiaDR3:\n"
                Teff=f"{pa['Teff'][0]:.0f}"
                logg=f"{pa['logg'][0]:.2f}"
                try:
                    fe=  f"{pa['[Fe/H]'][0]:.2f}"
                except:    
                    fe=  f"{pa['FeAbundance'][0]:.2f}"
                dist=f"{pa['Dist'][0]:.0f}"
                AG=  f"{pa['AG'][0]:.2f}"
                if len(Teff)>0: s10+=" Teff = "+Teff+"\n"
                if len(logg)>0: s10+=" log g = "+logg+"\n"
                if len(fe)>0  : s10+=" Fe/H = "+fe+"\n"
                if len(dist)>0: s10+=" Distance = "+dist+" pc\n"
                if len(AG)>0:   s10+=" A(G) = "+AG+" mag"
                s10+=s9+"\n"
            if FAP < 0.03:
                s11 +=f"\n FAP = {FAP:.3f}\n Period={period/365.25:.3f} d"  
            ax2.text(0.02,0.15,
                f"RA={pa[RA2000][0]:10.5f} deg\nDec={pa[DE2000][0]:10.5f} deg\n"+
                    f"classification={cl}\nfilter={band[3].upper()}\n"+s10+s11
                    #f"={cl}\n{lc_class}\nfilter={band[3].upper()}\n"+s10+s11
                  ,ha='left',fontsize=8,
                  )
        else:           
           ax2.text(0.02,0.3,f"RA={pa[RA2000][0]:10.5f} deg\nDec={pa[DE2000][0]:10.5f} deg\n"+
                    f"classification={cl}\nfilter={band[3].upper()}\n{s9}\n",ha='left',fontsize=9,
                    #f"={cl}\n{lc_class}\nfilter={band[3].upper()}\n{s9}\n",ha='left',fontsize=9,
                  )
        """          
        print ('FAP=',FAP)
        if FAP < 0.030 : 
            print ('setting labels ax3')        
            ax3.set_xlabel(r'frequency (yr$^{-1}$)') 
            ax3.set_ylabel('normalised power')
        else: 
            ax3.remove                    
            plt.show()
        """    
        if FAP > 0.03:
           ax3 = fig.add_axes([0.55,0.10,0.40,0.35]) # for DSS plot
           ax3.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
           pos = coord.SkyCoord(f"{ra} {dec}", unit=(units.deg, units.deg) )
           path = SkyView.get_images(position=pos,survey=['DSS'],scaling='log',grid=True,gridlabels=True, radius=radius)
           img = path[0]
           ax3.imshow(img[0].data)
           SkyView.clear_cache()
           
        plt.savefig(fn.split('fits')[0]+'pdf')
        print (f"saved ",fn.split('fits')[0]+'pdf')
        plt.close()
        msf = "%60s %6s %s %.2f"%(fn.ljust(60),cl,class2,FAP)
        msf2= "%55s %10.4f %9.3f %2i %2i"%(fn.ljust(55),FAP, period/365.26, ierup, idimming)
        command=f"echo '{msf2}' >> lc_summary.in"
        print (command)
        os.system(command)
        
def make_sed(t): 
    """ 
    typically data is an array with photometry from the SUSS SSC variables
    
    """
    import numpy as np
    bands = {'uvw2':'UVW2_AB_MAG','uvm2':'UVM2_AB_MAG','uvw1':'UVW1_AB_MAG',
      'OM-u':'U_AB_MAG',  'U':'umag_M', 'OM-b':'B_AB_MAG', 
      'OM-v':'V_AB_MAG', 'V':'vmag_M', 'Gaia-G':'Gmag',  'g':'gmag_M',
      'r':'rmag_M', 'i':'imag_M',   'z':'zmag_M',  'y':'ymag_M', 
        'J':'jmag_M',   'H':'hmag_M', 'K':'kmag_M',  'WISE-W1':'WI_W1mag',
         'WISE-W2':'WI_W2mag', 'WISE-W3':'WI_W3mag', 'WISE-W4':'WI_W4mag',
         }
    bandinv = {'UVW2_AB_MAG':'uvw2','UVM2_AB_MAG':'uvm2','UVW1_AB_MAG':'uvw1',
      'U_AB_MAG':'OM-u', 'umag_M':'U', 'B_AB_MAG':'OM-b', 
      'V_AB_MAG':'OM-v', 'vmag_M':'V', 'Gmag':'Gaia-G',  'gmag_M':'g',
      'rmag_M':'r', 'imag_M':'i',   'zmag_M':'z',  'ymag_M':'y', 
        'jmag_M':'J',   'hmag_M':'H', 'kmag_M':'K',  'WI_W1mag':'WISE-W1',
         'WI_W2mag':'WISE-W2', 'WI_W3mag':'WISE-W3', 'WI_W4mag':'WISE-W4',
         }
    wcentral = {'uvw2':193, 'uvm2':225, 'uvw1':260,"OM-u":350.,"OM-b":433,"OM-v":540,
    'U':359, 'B':440  , 'V':550.  ,'Gaia-G':650  ,'g':476  ,'r':620  ,'i':762 ,
    'z':912.  ,'y': 960 ,'J':1250  ,'H':1635  ,'K': 2200 ,'WISE-W1': 3400. ,
    'WISE-W2':4600  ,'WISE-W3': 12000. ,'WISE-W4': 22000 }

    # sort wcentral by wl
    
    # find the ones that have SED data    
    has_data = []
    for key in list(bands.values()):
        if t[key] > 0: 
            has_data.append(bandinv[key]) 
    
    # require more than 2 photometry points
    if len(has_data) < 2: 
        return None
    x1 = []
    y1 = []
    for key in has_data:
        #print ("key:",key)
        x1.append(wcentral[key])
        ykey = bands[key]
        #print ("ykey = ",ykey)
        y1.append(t[ykey][0])   
    x1 = np.array(x1) 
    y1 = np.array(y1)
    x1 = np.log10(x1)
    return x1, y1 

                                                              

def fix_multiples(cat="v3/XMM-OM-SUSS5.0ep_singlerecs_v2_7291gal_7294qso.fits",
    indir='/Volumes/DATA11/data/catalogs/suss_gaia_epic/',cols=["SRCNUMS","OBSIDS","EPOCHS"]):
    """
    clean up the OBSIDS, EPOCHS, SRCNUMS columns by removing duplicates
    assumes only one table in file
    
    2024-01-02 Decided not to do that for the main enhances/auxil +gaia file
    
    """
    from astropy.table import Table
    t = Table.read(indir+cat)
    print (f'finished reading {cat}')
    for name in col: 
        c = t[name]
        for r in c: 
            xxx = r.split("_")
            r = np.unique(xxx) 
            # assumes this updates c value
        t[name] = c
        # assume this updates t
    t.write(indir+cat+".1")
    # rename to originan file if OK after check
     

def classify_light_curve(time, flux, err, fig, noplot):
    # 'flux, err are AB mags
    # for plot_lc calls 2024-02-16
    import numpy as np
    from astropy import units as u
    from astropy.timeseries import LombScargle
    from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap
    from astropy import units as u
    
    # offset time by launch time (obs_epoch) 
    time = (time - 1995.0) # * u.yr  # unit = years
    flux = flux # * u.mag
    err = err # * u.mag
    freq_min, freq_max = 1./60., 1/0.0015 #   np.array([60.,0.0015]) * u.yr

    # Compute Lomb-Scargle periodogram (astropy 5.3.4 does not support all)
    
    frequency, power = LombScargle(time, flux, err).autopower(
        minimum_frequency=freq_min,
        maximum_frequency=freq_max,
        #fit_mean=True,
        #nterms=3
        )
    
    LS = lomb_scargle(time, flux, dy=err, omega=frequency)
    # Bootstrap to estimate significance

    significance = lomb_scargle_bootstrap(time, flux, dy=err,
                 N_bootstraps=100, random_state=0,omega=frequency)
    
    
    max_power_freq = frequency[np.argmax(power)] 
    period = 1./max_power_freq
    try:
       max_significance_freq = FAP = frequency[np.argmax(significance)]
       print (f"max power freq at={max_power_freq}\tperiod={period}")
       print (f"max significance freq at {max_significance_freq}")
    except:
       max_significance_freq = FAP = 99.
       print (f"2423 skipping Lomb-Scargle ")
       pass   
    
    #  astropy v6 ?  FAP = LS.false_alarm_probability(power.max(),method='bootstrap')   
    
    # reporting 
    
    class1 = ""
    class2 = ""
    eclipse = ""
    period = -1.
    minmax = [np.min(flux),np.max(flux)]
    
    if (FAP <= 0.03) and (not noplot):
        print (f"FAP={FAP}...plotting freq.-power ")
        ax = fig.add_axes([0.55,0.10,0.40,0.35])  # for PSD plot
        ax.tick_params(top=False,bottom=False,left=False,right=True,labelleft=False,labelbottom=True)
        ax.plot(frequency, power)
        ax.legend(title=f"@freq:{max_power_freq:.2f}:\nnormalised FAP={FAP:.2f}")
        ax.set_xlabel(r'frequency (yr$^{-1}$)') 
        ax.set_ylabel('normalised power')
    
        # for FAP < 0.01 assume periodic 
    if FAP < 0.01: 
            class1 = f"periodic: P = {period:.4f} yr\n"
            class2 = f"periodic\tP = {period:.4f} yr"
        #if FAP > 0.9: 
        #else:
        #    class1 = f"likely non-periodic\n" 
        #    class2 = f"non-periodic" 
        
    # test for eruptive : 10% or less have 5 sigma above median flux  and 3 mag 
    print (f"eruptive? {np.median(flux)} - 5x {flux.std()}   ")
    q = flux < (np.median(flux) - np.max([5.0 * flux.std(), 1.0]) )  
    print ("q.sum=",q.sum()," flux[q]=",flux[q]," median flux=",np.median(flux))
    if (q.sum() > 0) and (q.sum() < 0.1 * len(flux)): 
        erup = f"probably eruptive"
        ierup = 1
    else:
        erup = "" #f"likely non-eruptive"  
        ierup = 0
        
    # test for eclipse/absorption (flux is magnitudes)
    q = flux > np.median(flux)+5.0 * (0.03 + err) # (0.03 mag for systematic error)
    if len(q) > 0:
        eclipse = f"possible dimming"   
        idiming = 1
    else:
        idiming = 0    

    classification = f"{class1}{erup}\n{eclipse}"
    classif2 = f"{class2}\t{erup}\t{eclipse}"

    return classification, classif2, FAP, period, minmax, ierup, idiming


def concatenate_variable_sussfiles(
       indir="/Volumes/DATA11/data/catalogs/suss_gaia_epic/var_lc/lc_fits/"):
    """
    We created a fits file with the enhanced SUSS records for each variable. 
    For easy searching and stats we want to concatenate the first extensions 
    of those files. 
    """
    import os
    from astropy.table import Table
    files = os.listdir(indir)
    t0 = Table()
    for file in files:
        if ".fits" in file:
            t=Table.read(file,ext=1)
            t0 = vstack([t0,t])
    t0.write(f"{indir}/../concatenated_variable_sussfiles.fits")
 
def preprocessing(
        minnumber=10,  
        chi2_red_min=5., 
        inputdir='/Volumes/DATA11/data/catalogs/suss_gaia_epic/v3/',
        inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6_UVqual=0.fits",
        onlyqualzero=False,
        min_srcdist=6.0,
        chatter=2,
        minpnts=5, 
        threesigma=5.0,
        catversion=16,
    ): 
    """
    first extract the set of light curves for selected sources
    
    The selection criteria are:
       1. data quality good: <band>_quality_flag=0
       2. more than 5 good quality observations
       3. distance to nearest source > 6"
       4. chi-squared-reduced > 5
       
    mkdir var_lc
    
    do in ipython:
    from cats import xmm2athena as xmm2
    xmm2.make_file_variable(band='uvw2',chi2_red_min=5.,minnumber=5,min_srcdist=6.)
    xmm2.make_file_variable(band='uvm2',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='uvw1',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='u',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='b',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='v',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
   
    found 3448 sources (w2:263,m2:150,w1:1278,u:702,b:506,v:549)  
    
    where def make_file_variable( band='uvw1', minnumber=10, maxnumber=1000, 
         chi2_red_min=5., 
         inputdir='/Volumes/DATA11/data/catalogs/suss_gaia_epic/v3/',
            ##inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6xSIMBAD.fits.gz",
         inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits",
         selectrecs=[0,60000000.], plotit=False, onlyqualzero=False,
         min_srcdist=6.0,
         chatter=2):
   
    in shell do: 
    cd var_lc
    ls -1  *.fits > allsources.txt
    
    in topcat edit allsources.txt:
    topcat: new col; test=split(col1,"_"); name=concat(test[0]," ",test[1])
    rename col1 filename
    
    Select all sources in 
        inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits",
    which match allsources.txt (but need to replace _ with " " in name
        match with all from allsources.txt: name , filename 
        => SimbadxSUSS5_variable_sources_v8.fits
    use "best match for each table 1 row" and "all from 1"     
    
    match this to CDS/GaiaVarSumm to find matches to Gaia variables 
    add column GaiaVar after using a match to Gaia Var in topcat and retaining the 
        object variable classification, and Gmagmean for selection as GaiaVar
    add column SimbadVar for the match to Simbad variables (see below for the selection)
    --> 268 vars (out 3448) 7.7%
    
    match this to Simbad 
    
    select from Simbad/Main_type the variable types identified by Ada
    main_type=="QSO"|main_type=="AGN"|main_type=="BLLAC"|main_type=="Blazar"|
    main_type=="SN"|main_type=="**"|main_type=="EB*"|main_type=="Cepheid"|
    main_type=="AeBe"|main_type=="BYDra"|main_type=="Eruptive"|main_type=="XRB"|
    main_type=="LPV"|main_type=="Mira"|main_type=="OrionV"|main_type=="Nova"|
    main_type=="Pulsating*"|main_type=="Pulsar"|main_type=="Variable"|
    main_type=="HMXB"|main_type=="RRLyr"|main_type=="Rotation"|
    main_type=="SpBinaries"|main_type=="Symbiotic*"|main_type=="YSO"|
           or:
main_type=="QSO"|main_type=="AGN"|main_type=="BLLAC"|main_type=="Blazar"|main_type=="SN"|main_type=="**"|main_type=="EB*"|main_type=="Cepheid"|main_type=="AeBe"|main_type=="BYDra"|main_type=="Eruptive"|main_type=="XRB"|main_type=="LPV"|main_type=="Mira"|main_type=="OrionV"|main_type=="Nova"|main_type=="Pulsating*"|main_type=="Pulsar"|main_type=="Variable"|main_type=="HMXB"|main_type=="RRLyr"|main_type=="Rotation"|main_type=="SpBinaries"|main_type=="Symbiotic*"|main_type=="YSO"

    
    -> 174 variables (out of 3448) much overlap with GaiaVar
    add SimbadVar, GaiaVar boolean columns
    
    save the catalogue file:
        => SimbadxSUSS5_variable_sources_vXX.fits

    """
    import os
    from cats import xmm2athena as xmm2

    if not os.access(f"{inputdir}/../var_lc",os.F_OK):
       print (f"WARNING ould not locate the var_lc directory, making one below {inputdir}")
       os.system(f"mkdir {inputdir}/../var_lc")
    indir = f"{inputdir}/../var_lc/"
    os.system(f"cd {indir}")
    for band in ["uvw2","uvm2","uvw1","u","b","v"]:
        xmm2.make_file_variable(band=band,
            chi2_red_min=chi2_red_min,
            minnumber=minnumber,
            min_srcdist=min_srcdist,
            inputdir=inputdir,
            inputcat=inputcat,
            )
    print (f"writing the list of X*.fits files as allsources.txt")        
    os.system(f"ls -1 X*.fits > allsources.txt")
    #
    # make catalogue for selection - match allsources.txt to inputcat
    #
    TOPCATPATH = "/Users/kuin/bin" 
    print('match')
    # match 
    in1 = "allsources.txt"
    catin1 = "allsources_1.txt"
    catin2 = inputcat
    #outcat = f"SimbadxSUSS5_variable_sources_v{catversion}.fits"
    xgaiavar = "xgaiavar.fits"
    xsimbad = "xsimbad.fits"
    
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"addcol test split(col1,'_');addcol name concat(test[0],' ',test[1]);addcol filename col1 "+\
    f"in={in1} out={catin1}" 
    print('2600 ',command)
    status = os.system(command)
    if status > 0: 
        print(f"2600 ERROR in {command} ")
       
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"tmatch2 in1={catin1} ifmt1=ascii  in2={catin2} ifmt2=fits "+\
    f"out={outcat} ofmt=fits-basic  matcher=exact "+\
    f"values1='name' values2='IAUNAME'  "+\
    f" params={match_par}   join=all1 find=best1 fixcols=dups suffix1=  suffix2=suss  "+\
    f" ocmd=delcols 'GroupID';delcols 'GroupSize' "
    print('2600 ',command)
    status = os.system(command)
    if status > 0: 
        print(f"2600 FOUT!ERROR in {command} ")
        
    # match to GaiaVar
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"cdsskymatch cdstable=simbad in={catin2} ifmt=fits out={xsimbad} "+\
    f"icmd='keepcols  IAUNAME gaiadr3_ra gaiadr3_dec' "+\
    f"ocmd='keepcols main_type' radius=2.0 "+\
    f"ra=gaiadr3 dec=gaiadr3_dec find=best "
    
    # match to Simbad
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"cdsskymatch cdstable=I/358/varisum in={catin2} ifmt=fits out={xgaiavar} "+\
    f"icmd='keepcols  IAUNAME gaiadr3_ra gaiadr3_dec' "+\
    f"ocmd='keepcols Gmagmean' radius=2.0 "+\
    f"ra=gaiadr3 dec=gaiadr3_dec find=best "
    
    # merge together 
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"tmatch2 in1={catin1} ifmt1=fits  in2={xgaiavar} ifmt2=fits "+\
    f"out={outcat} ofmt=fits-basic  matcher=exact "+\
    f"values1='name' values2='IAUNAME'  "+\
    f" params={match_par}   join=all1 find=best1 fixcols=dups suffix1=  suffix2=suss  "+\
    f" ocmd=delcols 'GroupID';delcols 'GroupSize' "
    
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"tmatch2 in1={catin1} ifmt1=ascii  in2={xsimbad} ifmt2=fits "+\
    f"out={outcat} ofmt=fits-basic  matcher=exact "+\
    f"values1='name' values2='IAUNAME'  "+\
    f" params={match_par}   join=all1 find=best1 fixcols=dups suffix1=  suffix2=suss  "+\
    f" ocmd=delcols 'GroupID';delcols 'GroupSize' "
    
    # create new col GaiaVar and SimbadVar
    # create IAUNAME from filename in file :
    # topcat: new col; test=split(col1,"_")
    # name=concat(test[0]," ",substring(test[1],0,9),substring(test[1],10,17))
    
    
    
def post_process_variables1(       
    #indir="/Volumes/DATA11/data/catalogs/suss_gaia_epic/var_lc/",
    indir = inputdir+var_lc,
    minpnts=minnumber, 
    threesigma=chi2_red_min,
    ):
    """
    first extract the set of light curves for selected sources
    
    The selection criteria are:
       1. data quality good: <band>_quality_flag=0
       2. more than 5 good quality observations
       3. distance to nearest source > 6"
       4. chi-squared-reduced > 5
       
    mkdir var_lc
    
    do in ipython:
    from cats import xmm2athena as xmm2
    xmm2.make_file_variable(band='uvw2',chi2_red_min=5.,minnumber=5,min_srcdist=6.)
    xmm2.make_file_variable(band='uvm2',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='uvw1',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='u',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='b',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
    xmm2.make_file_variable(band='v',chi2_red_min=5.,minnumber=5,min_srcdist=6.,chatter=3)
   
    found 3448 sources (w2:263,m2:150,w1:1278,u:702,b:506,v:549)  
    
    where def make_file_variable( band='uvw1', minnumber=10, maxnumber=1000, 
         chi2_red_min=5., 
         inputdir='/Volumes/DATA11/data/catalogs/suss_gaia_epic/v3/',
            ##inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6xSIMBAD.fits.gz",
         inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits",
         selectrecs=[0,60000000.], plotit=False, onlyqualzero=False,
         min_srcdist=6.0,
         chatter=2):
   
    in shell do: 
    cd var_lc
    ls -1  *.fits > allsources.txt
    
    in topcat edit allsources.txt:
    topcat: new col; test=split(col1,"_"); name=concat(test[0]," ",test[1])
    rename col1 filename
    
    Select all sources in 
        inputcat="XMM-OM-SUSS5.0ep_singlerecs_v4_classified_v6.fits",
    which match allsources.txt (but need to replace _ with " " in name
        match with all from allsources.txt: name , filename 
        => SimbadxSUSS5_variable_sources_v8.fits
    use "best match for each table 1 row" and "all from 1"     
    
    match this to CDS/GaiaVarSumm to find matches to Gaia variables 
    add column GaiaVar after using a match to Gaia Var in topcat and retaining the 
        object variable classification, and Gmagmean for selection as GaiaVar
    add column SimbadVar for the match to Simbad variables (see below for the selection)
    --> 268 vars (out 3448) 7.7%
    
    match this to Simbad 
    
    select from Simbad/Main_type the variable types identified by Ada
    main_type=="QSO"|main_type=="AGN"|main_type=="BLLAC"|main_type=="Blazar"|
    main_type=="SN"|main_type=="**"|main_type=="EB*"|main_type=="Cepheid"|
    main_type=="AeBe"|main_type=="BYDra"|main_type=="Eruptive"|main_type=="XRB"|
    main_type=="LPV"|main_type=="Mira"|main_type=="OrionV"|main_type=="Nova"|
    main_type=="Pulsating*"|main_type=="Pulsar"|main_type=="Variable"|
    main_type=="HMXB"|main_type=="RRLyr"|main_type=="Rotation"|
    main_type=="SpBinaries"|main_type=="Symbiotic*"|main_type=="YSO"|
           or:
main_type=="QSO"|main_type=="AGN"|main_type=="BLLAC"|main_type=="Blazar"|main_type=="SN"|main_type=="**"|main_type=="EB*"|main_type=="Cepheid"|main_type=="AeBe"|main_type=="BYDra"|main_type=="Eruptive"|main_type=="XRB"|main_type=="LPV"|main_type=="Mira"|main_type=="OrionV"|main_type=="Nova"|main_type=="Pulsating*"|main_type=="Pulsar"|main_type=="Variable"|main_type=="HMXB"|main_type=="RRLyr"|main_type=="Rotation"|main_type=="SpBinaries"|main_type=="Symbiotic*"|main_type=="YSO"

    
    -> 174 variables (out of 3448) much overlap with GaiaVar
    add SimbadVar, GaiaVar boolean columns
    
    save the catalogue file:
        => SimbadxSUSS5_variable_sources_vXX.fits
  - - - -      
    Remove obsids which have maximum brightness (~3sigma from median or > 1mag) in 
    multiple sources or have a dimming in excess (~3sigma or > 0.2 mag)in multiple 
    sources for the same obsid. (see code below) 
      ... 5 sigme does not find any multiple sources with same obsid
    
    -> write a file with bad obsids. 
    -> rerun the above with the bad obsids input file   
       the copy of the input file (temp_file.fits) has been updated  
    -> remove the XMM*.fits files 
       check the bad_obsids.txt file
    now we have hopefully culled the obsids which were systematically high or low. 
    
    -> rerun make_file_variable() with indir = './' and inputcat='temp_file.fits'
    
    Update the above fits files per source, then the following:    
        
    for each source, 
       for each of its OBSIDs 
          find the number of sources with that obsid -> list
       take the maximum number of sources within that list 
    => source, ra, dec, number of sources 
          
    compare all sources to those flagged as Gaia Variables to determine maximum 
    number of sources acceptable within an obsid
    
    remove the sources which exceed the maximum number 
    
    Remove obsids which have > 20 Var Srcs of which < 20% are Gaia Variable
    
    rerun the above => v9 
    
    make a list of obsids where there is a significant brightening, and then 
    make a histogram of the number of sources with brightning per obsid.
    
    
    ========================
    further post-processing: 
    
    match to nearby galaxies and remove sources which are inside the galaxy 
    
    match to bright nebulae and remove sources 
    
    Paired correlated lightcurves: might be because of scattered light features
    limit search for source correlation to other source that have a matched obsid
    output source, matched source, obsid with the high point? 
    
    More general approach:
    for all sources 
       for each obsid 
           compute the correlation with random 2 other sources within a 50(tbd) arcsec radius
    -> source, obsid, correlation
    sort by obsid       
    for high correlations in an obsid, 
       select its sources 
          select for each source the obsids
              remove one obsid at a time and compute the correlation/chi-squared
              a significant? drop will flag the obsid which causes the flare
          find out which of the same obsid causes a flare in all sources
          remove such an obsid from all light curves (reprocess)
          
              
    """
    from astropy.table import Table
    from astropy.io import fits
    import os
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    
    file1 = "SimbadxSUSS5_variable_sources.fits"
    file2 = "temp_file.fits"
    os.system(f"cp {indir}{file1} {file2}")
    with fits.open(indir+file1) as f1:
       t0 = Table(f1[1].data)
       gaiavar=t0['GaiaVar']
       results = []
       for row in t0:
           source = row['IAUNAME']
           filename = row['filename'].split("/")[-1]
           obsids = row['OBSIDS'].split("_")
           ra = row['RAJ2000Ep2000']
           de = row['DEJ2000Ep2000']
           results.append([source, filename, obsids,ra,de])

 # - - - - 
   
    # test for eruptive : have 5 sigma (staterr+syserr) above median flux 
    # test for dimming  : have 5 sigma (staterr+syserr) above median flux 
    obspeak=[]
    syserr = 0.02
    with open("peak_in_obsid_list.txt","w") as bo:
        for res in results:
            filen = res[1]
            band = filen.split("_")[3]
            #obses = res[2]            
            with fits.open(filen) as gd:
                obses = gd[1].data['OBSID']
                mag   = gd[1].data[f"{band.upper()}_AB_MAG"]
                magerr   = gd[1].data[f"{band.upper()}_AB_MAG_ERR"]
                print (f"magnitudes = {mag}")
            # 23 Aug - adding error in the mean as proxy for that in the median
            # 27 aug - removing the maximum and minimum points for the test of outliers
            medmag=np.median(mag)
            minmag=np.min(mag)
            maxmag=np.max(mag)
            magdel=[np.abs(medmag-minmag),np.abs(maxmag-medmag)]
            if len(mag)>4:
               stdmean = mag[(mag<maxmag)&(mag>minmag)].std()
            else:
               stdmean = mag.std()
            mval = np.abs(  (mag - np.median(mag)) / (threesigma * (magerr+syserr+stdmean) )  ) 
            q = mval > 1.
            print (f" median {np.median(mag)}  {threesigma} sigma {threesigma * (magerr+syserr)}",
                f"is {mval} > 1 ?" )
            if q.sum() > 0:
                obsmem = obses, mag, q
                print (f"found peaks in {filen}, {mag} {obses} {q}\n")
                if q.sum() == 1:
                    bo.write(f"{obses[q][0]}\n")
                    obspeak.append(f"{obses[q][0]}")
                else:    
                   for o in obses[q]:
                      bo.write(f"{o}\n")
                      obspeak.append(o)
    uniq_obspeak=np.unique(obspeak) 
    obspeak = np.array(obspeak)                     
# - - - -
    obspeakstat=[]
    min_peaks = 2   # for N > 5 test 23 Aug
    for o1 in uniq_obspeak:
        no1 = 0
        q= o1 == obspeak
        no1 = q.sum()
        obspeakstat.append([o1,no1])
    obspeakstat = np.array(obspeakstat,dtype=int)
    if len(obspeakstat) > 1:    
       npeaks = obspeakstat[:,1]
       q = np.where(npeaks > min_peaks)           
       baddies = obspeakstat[q,:][0]
       badobs = baddies[:,0]
       
       # remove / update the peaky obsids from the file 
       update_cataloguefile(file2,badobs)
                   
    else:
       npeaks=0
       baddies= []
       badobs = []   

    if len(badobs) > 0:
        with open("remove_jumpy.txt","w") as jump:
            for k in badobs:
                jump.write(f"{k}\n")      

 # - - - - SECOND RUN now that the inputcat has been updated (to file2) 
    for band_ in ["uvw2","uvm2","uvw1","u","b","v"]:
        make_file_variable( band=band_, minnumber=10, maxnumber=1000, 
         chi2_red_min=5., inputdir="./",inputcat=file2,min_srcdist=6.0,
         bad_obsidfile="remove_obsids.txt",chatter=2)

    """
    for band_ in ["uvw2","uvm2","uvw1","u","b","v"]:
        xmm2.make_file_variable( band=band_, minnumber=minpnts, maxnumber=1000, 
         chi2_red_min=5., inputdir="./",inputcat=file2,min_srcdist=6.0,
         bad_obsidfile="remove_obsids.txt",chatter=2)
    """

 #    cull obsids with many sources, but not enough GaiaVar sources
    
    with fits.open(file2) as f1:
       t0 = Table(f1[1].data)
       gaiavar=t0['GaiaVar']
       results = []
       for row in t0:
           source = row['IAUNAME']
           filename = row['filename']
           obsids = row['OBSIDS'].split("_")
           ra = row['RAJ2000Ep2000']
           de = row['DEJ2000Ep2000']
           results.append([source, filename, obsids,ra,de])
    obs = []

    for x in results:
        oo = x[2]       
        for o in oo:
            if o == "":
               continue
            else:   
               obs.append(o)
    uniq_obsids = np.unique(obs) 
    print (f"number of unique obsids is {len(uniq_obsids)}")
# - - - -  
    # find for each source 
    source_results = []
    for x in results:
        source = x[0]
        obs = x[2]
        n_src_list = [] 
        # for each obsid find matching sources 
        # keep track how many gaiavar sources in each obsid
        n_gaiavar_list = []
        for o in obs:
            n_src = 0       # count for a given obsid
            n_gaia_srcs = 0
            for r,gv in zip(results,gaiavar):  # not excluding the originating source 
                if o in r[2]:
                    n_src += 1
                    if gv: n_gaia_srcs += 1 
            n_src_list.append(n_src)
            n_gaiavar_list.append(n_gaia_srcs)
        source_results.append([source,np.max(n_src_list),n_src_list,n_gaiavar_list]) 

    n1=[]   # list the maximum number of sources in any of the obsids 
    for x in source_results: 
        n1.append(x[1])        
    n_all = np.array(n1)        
    
    # select the subset of sources which are GaiaVar 
    n_var = np.array(n1)[gaiavar] 
    
    # for all variable, select max number gaia var
    n1=[]
    for x in source_results:
       xm = 0
       try:
           xm = np.max(x[3])
       except:
           print (f"error max gaiaVar: \n{x}\nxm={xm}")    
       n1.append(xm)
    n_gaia = np.array(n1)                                
 # - - - - - figure();plot(n_all,n_gaia/n_all,'+b',)
          
    bad_obsid = []  
    # define the limits for exclusion 
    n_src_limit = 20
    n_var_ratio_limit = 0.20  
    # find sources with an obsid which is bad because too many sources with few/no GaiaVars
    Nsrc = len(results)   # loop over source
    arebad = (n_all < n_src_limit) & (n_gaia/n_all > n_var_ratio_limit)
    bad_res = []
    bad_src = []
    n_gaia_4these = []
    for k in range(Nsrc):
        bad_res.append(results[k])
        bad_src.append(source_results[k])
        n_gaia_4these.append(n_gaia[k])
    # now we need to filter out the bad obsids from the good ones     
    for o in uniq_obsids:  # for each OBSID check if good or bad
        for k in range(len(bad_res)):
            res    = bad_res[k]
            srcres = bad_src[k]   
            ngaia1 = n_gaia_4these[k]
            obsid1 = res[2]   # list of obsids 
            n_src1 = srcres[2]  # matching list of sources 
            for obs,nnn in zip(obsid1,n_src1):
                 if o == obs: # match outer loop
                     # print (f"{o}  x={nnn}  y={ngaia1/nnn} {ngaia1}")
                     if (nnn > n_src_limit) & (ngaia1/nnn < n_var_ratio_limit):  # bad selection
                         bad_obsid.append(o)
    uniq_bad_obsids=np.unique(bad_obsid)
    
    good_obsids = []
    for obs in uniq_obsids:
        if obs in uniq_bad_obsids:
            continue
        else: good_obsids.append(obs)    
    print (f"number of bad obsids (too many sources, not enough Gaia Var's) is {len(uniq_bad_obsids)}")    
# - - - -    plot(n_all,n_gaia/n_all,'+g',)
      
    with open(indir+"good_obsids.txt","w") as gf:
        for obs in good_obsids:
           gf.write(f"{obs}\n")
    with open(indir+"bad_obsids.txt","w") as gf:
        for obs in uniq_bad_obsids:
           gf.write(f"{obs}\n")
                           
    # remove / update the obsids with many sources above expected from the file 
    file3 = f"step3_{file2}"
    os.system(f"cp {file2}  {file3}")
    update_cataloguefile(f"step3_{file2}",uniq_bad_obsids)
 # - - - -                   
    # move the results from step 1 out of the way
    os.system("mkdir var_lc_step1/;mv X*.fits var_lc_step1")
    # rerun with corrected obsids 
    
    for band_ in ["uvw2","uvm2","uvw1","u","b","v"]:
         make_file_variable( band=band_, minnumber=minpnts, maxnumber=1000, 
         chi2_red_min=5., inputdir="./",inputcat=file3,min_srcdist=6.0,
         bad_obsidfile="remove_obsids.txt",chatter=2)
    """     
    for band_ in ["uvw2","uvm2","uvw1","u","b","v"]:
        xmm2.make_file_variable( band=band_, minnumber=minpnts, maxnumber=1000, 
         chi2_red_min=5., inputdir="./",inputcat=file3,min_srcdist=6.0,
         bad_obsidfile="remove_obsids.txt",chatter=2)
    """ 
    # create light curve plots (may need to adjust minpnts)
    os.system("ls -1 X*.fits > newfiles.txt")
    f2 = open("newfiles.txt")
    files = f2.readlines()
    f2.close()
    for i, lcfile in enumerate(files[0:]):
        file=lcfile.split("\n")[0]
        print (f"==============\nfile {i}  {file}")
        trynumber = 1
        while trynumber < 4:
            try:
                plot_lc(file,simbad=False,minpnts=minpnts,noplot=False)
            # xmm2.plot_lc(lcfile.split("\n")[0],simbad=False,minpnts=minpnts,noplot=False)
                trynumber=10
            except:
                plot_lc(file,simbad=False,minpnts=minpnts,noplot=False)
                trynumber+=1
                pass          
    
    # repeat earlier code with 'file3' to recompute n_all, n_gaia (not implemented here)
     
"""
  Code to extract the summary/index file to the light curve directory

$ ls -1 > allfiles__.txt

python:
from astropy.io import fits
from astropy.table import Table

g = open('allfiles__.txt')
fnames = g.readlines()
g.close()
# for each lightcurve add a row (with multiple rows for a single source possible)
name=fnames[0].split('\n')[0]
x = fits.open("lc_fits/"+name)
x1 = f[1].data
t1 = Table(x1)
x = fits.open("lc_fits/"+name)
x2 = x[2].data
tt = Table(x2)
N=len(fnames))
for k in np.arange(2,N):
    name=fnames[k].split('\n')[0]
    x.close()
    x = fits.open("lc_fits/"+name)
    x2 = x[2].data
    t2 = Table(x2)
    d2 = dict(t2)
    tt.add_row(vals=d2)
tt.write("SUSS_variables.fits")

unames=[]
for k in np.arange(N):
    unames.append(fnames[k][:21])
 
# find a single match to each in unames  
snames=[]
for k in unames:
    for m in fnames:
        if k == m[:21]: 
            snames.append(m)
            break

name=snames[0].split('\n')[0]
x = fits.open("lc_fits/"+name)
x2 = x[2].data
tt = Table(x2)
N=len(snames)
for k in np.arange(1,N):
    name=snames[k].split('\n')[0]
    x.close()
    x = fits.open("lc_fits/"+name)
    x2 = x[2].data
    t2 = Table(x2)
    d2 = dict(t2)
    tt.add_row(vals=d2)
tt.write("SUSS_variables_single.fits")

"""                
                             
    return n_all,n_gaia, bad_obsid, good_obsids, results, source_results
    
 # end post processing                               

 # end post processing                               

 # end post processing                               


def update_cataloguefile(file2,badobsids):
    """
    This procedure is to update the catalogue by removing 
    obsids in the three columns called OBSIDS, SRCNUMS, EPOCHS.
    
    """
    from astropy.io import fits
    
    print (f"updating {file2}")
    with fits.open(file2,"update") as f1:
         o1in = f1[1].data['OBSIDS']
         s1in = f1[1].data['SRCNUMS']
         e1in = f1[1].data['EPOCHS']
           # o1 is a list of obsids as a long string for each source 
           # need to iterate over the sources 
         for k1 in range(len(o1in)):  
             o1 = o1in[k1].split("_")
             s1 = s1in[k1].split("_")
             e1 = e1in[k1].split("_")
             
             o3 = ""  # to fill with updated lists
             s3 = ""
             e3 = ""
             
             for o2, s2, e2 in zip(o1,s1,e1):
                 # o2, s2, and e2 are strings
                 ok = True
                 for bk in badobsids:
                     if (len(o2)>0):
                         if (int(bk) == int(o2)): 
                             ok = False
                 if ok: 
                    o3 += f"{o2:10}_"
                    s3 += f"{s2}_"
                    e3 += f"{e2:10}_"        
             o1in[k1] = o3[:-1]       
             s1in[k1] = s3[:-1]       
             e1in[k1] = e3[:-1]       
         #update f1 
         f1[1].data['OBSIDS']  = o1in
         f1[1].data['SRCNUMS'] = s1in
         f1[1].data['EPOCHS']  = e1in

    print (f"finished update")
    """

fix for problem with IAUNAME missing by creating that from the filename file

    t = Table.read("newfiles.txt",format="ascii")    
    
    f = open("newiauname.txt","w")

In [164]: for x in t["filename"]:
     ...:     x1 = x.split("_")
     ...:     if (x1[1][10] == "+") | (x1[1][10] =="-"):
     ...:         f.write (f"{x1[0].strip()}_{x1[1][:9]}{x1[1][10:17]:7}\n" )
     ...:     else:
     ...:         f.write (f"{x1[0].strip()}_{x1[1][:9]}+{x1[1][10:17]:7}\n" )
     ...: 

In [165]: f.close()

    
    #hist(n_all,bins=40,color='orange',label='all sources')
    #hist(n_var,bins=40,color='b',label='Gaia Variables')
    #xlabel("for each source: max number of OBSIDs")     
     
    # next a try to remove sources; but in the end, we will remove obsids  
    qhigh = n_all > 150
    goodsrc=[]
    for k in range(len(results)):
        if not qhigh[k]:
            goodsrc.append(results[k])
    print (f"remaining good sources with Nmax < 150 is {len(goodsrc)}")
    highsrc = []
    for k in range(len(results)):
        if qhigh[k]:
            highsrc.append(results[k]) 
    print (f"in GaiaVars there are {len(n_var[n_var > 150])} sources with >150 other sources in an obsid")
  
    # make file to move the high a_all sources
    os.system(f"mkdir ./high_sourcenumbers")
    with open(indir+"move_high_srcnumber_data","w") as outf:
        for k in range(len(n_all)): 
            if n_all[k] > 150:
                file = results[k][1]
                file2 = file.split('.fits')[0]+'.pdf'
                outf.write(f"mv {indir}{file}  {indir}high_sourcenumbers/{file}\n") 
                outf.write(f"mv {indir}{file2}  {indir}high_sourcenumbers/{file2}\n") 
    #os.system(f"chmod a+x {indir}move_high_srcnumbers;{indir}/high_source_numbers/")  

    
    # screen for bright nebula sources 
    # lynds catalogue 7009: 
    lynds = fits.getdata('7009.fit',ext=1)
    lyRA  = lynds.field('_RAJ2000')
    lyDE  = lynds.field('_DEJ2000')
    diam1  = lynds.field('Diam1')
    diam2  = lynds.field('Diam2')
    diam=np.max([diam1,diam2],axis=0)
                    
    in_a_cloud = []                                 
    for k in range(len(lyRA)):
        # look for sources in each lynds field of bright nebulosity
        c1 = SkyCoord(lyRA[k]*u.deg, lyDE[k]*u.deg, frame='icrs')
        cldiam = diam[k]*u.arcmin
        for kk in range(len(n_all)):
           srcra = results[k][3]
           srcde = results[k][4]
           c2 = SkyCoord(srcra*u.deg, srcde*u.deg, frame='icrs')
           sep = c1.separation(c2)
           if sep < cldiam:
              in_a_cloud.append(results[k])
    # result : in_a_cloud=[]     -- no sources inside a bright cloud     
          
    """     
 

"""    
from astropy import coordinates as coord

a = coord1[0][0].split('.')    
ra = f"{a[0][:-4]}:{a[0][-4:-2]}:{a[0][-2:]}.{a[1]}"

b = coord1[0][1].split(".")
dec = f"{int(b[0][:-4]):02}:{b[0][-4:-2]}:{b[0][-2:]}.{b[1]}"

pos = coord.SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg) )

 
 
"""    



def fix_multiple_obsids_in_SUSS5(t,filt,chatter=0):
    """
    In SUSS5 for each OBSIDs all exposures should have been summed into one row. 
    Unfortunately, for some sources there are two rows. Assuming that one row is the 
    summed data, the other not, and discriminating on <filt>_AB_MAG_ERR we remove 
    the lesser boy and return the corrected Table. 
    
    Requirement: input Table is for a single SRCNUM and FILTER filt. 
    """
    # check type t is Table...
    t.sort('OBSID')
    bad_rows = [] # list of bad row numbers
    #row_a = t[0]
    k1 = 0
    kays = [0]
    bad_rows = [-1]
    for k in np.arange(1,len(t)):
        row_a =  t[k1]
        if chatter > 1: print (f"k = {k} , k1 = {k1}")
        err_a = row_a[filt.upper()+'_AB_MAG_ERR']
        row_b = t[k]
        err_b = row_b[filt.upper()+'_AB_MAG_ERR']
        if chatter > 1: print (f"{row_a['OBSID']} {row_b['OBSID']}  {err_a} {err_b}  ")
        if row_a['OBSID'] == row_b['OBSID']:
            # find the bad one (was: assume only pairs) 
            kays.append(k)
            if chatter > 1: print (f"SAME OBSIDS: ")
            s1 = row_a[f'{filt.upper()}_AB_MAG_ERR'] < row_b[f'{filt.upper()}_AB_MAG_ERR'] 
            if s1:
               if bad_rows[-1] not in kays: #k-1:  # more than 2 in a row
                  bad_rows.append(k1)
            else:  # error row_b < err row_a
               k1 = k
               if bad_rows[-1] not in kays: 
                  bad_rows.append(k)      
            if chatter > 1: print (bad_rows)    
        # reset
        else:  # no OBSID match 
           k1 = k 
           kays = [k1]             
    # now remove the bad rows
    a = bad_rows.pop(0)
    if chatter > 1 : print (f"BAD ROWS are:{bad_rows}")
    if len(bad_rows) > 0:
        if chatter > 0: print (f"2be removed {len(bad_rows)} erronous OBSIDs")
        t.remove_rows(bad_rows)   
    return t

################# end 