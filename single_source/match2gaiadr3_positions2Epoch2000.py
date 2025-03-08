#!python3

"""
match2gaiadr3_positions2Epoch2000.py : 

This program does the initial  processing of the SUSS/UVOTSSC tables to 
add ObsEpoch, positions at Epoch 2000, and Epoch 2016 (the Gaia DR3 epoch)

Prior to executing this, we retrieved the subset from Gaia with higher PMs (2 Mar 2023) 
from the Gaia archive, https://gea.esac.esa.int/archive/ (have to log in). 
We retrieved all gaia with PM > 10 mas/yr Epoch 2016, for matching our source to 
gaiadr3 source, listing the matching Gaia ID.  

SELECT TOP 30100100 gaia_source.designation,gaia_source.source_id,gaia_source.ra,
gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,
gaia_source.parallax_over_error,gaia_source.pm,gaia_source.pmra,gaia_source.pmra_error,
gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,
gaia_source.bp_g,gaia_source.g_rp,gaia_source.non_single_star
FROM gaiadr3.gaia_source 
WHERE (gaiadr3.gaia_source.pm>=25.0)

=> gaiadr3_pmgt25.fits 

    match the SUSS/UVOTSSC catalogue to Gaia DR3 astrometry for  obs_epoch in catalogue
    and then compute the positions for Epochs in list, typically the Gaia epoch, and 
    epoch 2000. 
    
    Method:
    
    (A) independently retrieve the higher PM sources from gaiadr3 => gaiadr3_pmgt10.fits 
        This was completed 2023-march2  
        (using the gaiadr3_lite, PM > 25 mas/yr > gaiadr3_highpm.fits 
         2 March 2023 for > 10mas/yr => gaiadr3_pmgt10.fits )
        18 Feb 2025 back to a 25 mas/yr 
           
    (B) The retrieved catalogue is large,  and to speed up things, it has been subsetted
        using the class split_catalogue -> split_catalogue.
    
    (1) processing: done for a limited range in Epoch, make list of obsids
        use obsid list to find pointing RA,DEC
        restrict processing using that range
    
    (C) loop over epochs involves :
        - extract the tiles near the observed image pointing and on those select
          with radius > 20 arcmin to make a list gtemp
        - extract suss data for obsid 
        - compute the 2000 epoch for gtemp positions
        - change gtemp epoch from 2016.0 to suss using stilts function 
        - edit catalog upload columns if needed
        - %xarches[cds xmatch] job to match them by position after upload (retain all suss records)
        - match using skymatcherr 
        
    (3) apply proper motion correction to put data on the obs_epoch 
        of the SUSS data *** copy RA,DEC to new columns; filter on 
        'parallax' 'pm' since some sources don't have. 
        The apply PM correction on sources with PM, plx 
            
    (D) assume all suss sources without high PM match coordinates have coordinate near
        epoch 2016/2000 within 0.5 arcsec. 

    (D) merge the matched records -> finalise
    
    (E) clean up while doing this
    
    (5) write out a restricted columns catalog table [ra,dec,perr,pmra,pmde,gmag] 
    
Requires: astropy cdshealpix     

NPMK, April 2023, UVOTSSC2 processing (about 2 Epochs per hour on MacBook Pro 
    2019, 2.3GHz 8-core i9,32GB 2667MHz  memory) 
    SUSS5.0 two epochs take 440s system time.

Method: use_healpix works, the other, dividing into RA,Dec parcels makes 
   duplicates in output and does not cover the poles well. 
   
   
"""

import os
import numpy as np
import time
#import numpy.ma as ma
from astropy.table import Table, vstack
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
#    from astroquery.gaia import Gaia
from astropy.io import fits, ascii as ioascii
from cdshealpix import nested

########################
#  Global definitions  #
########################

ROOT = '/Volumes/DATA11/'  # data disk
TOPCATPATH = "/Users/kuin/bin" # /disk/xray16/npmk/bin
work = os.getcwd()+"/"

gaiacat  = ROOT+'data/catalogs/gaiadr3_pmgt25.fits' #  /disk/xray16/npmk/ ???
gaiapath = ROOT+'data/catalogs/gaiadr3_PM/' # store here the split up cat.
gaiaHP   = ROOT+'data/catalogs/gaiadr3_healpix.fits' # gaiacat with healpix col.
gaiaObsidTable = work+'/gaiadr3_obsid_tiles.fits'
GaiaEpoch = 2016.0

process_epoch_range = [2005.0, 2024.0] #[2007.7,2007.9]# [2005.0,2022.2]
process_epoch_step  = 0.1 # if making smaller, edit filenames which assume 0.1yr

catalog             = 'SUSS' #'test2' # options are: 'UVOTSSC2', 'SUSS', 'test'

RADIUS        = 20./60 # circular area in arcsec to include in search for PM objects
gaia_temp_cat = work+"/gaiadr3_selected_4epoch.fits" 
obs_cat       = work+'/epoch_slice_cat.fits'
epochfilesdir = work #ROOT+'data/catalogs/tests/'

use_healpix = True
chatter     = 0
do_epoch    = True #True if running for first time on chunks
prints2file = False #TBD not implemented

if catalog == 'SUSS':
   #rootobs    = ROOT+'data/catalogs/suss_gaia_epic/'
   #chunks     = ['XMM-OM-SUSS5.0.fits']
   rootobs    = ROOT+'data/catalogs/suss6.2/'
   chunks     = ['XMM-OM-SUSS6.2.fits']
   match_par  = 2.5  # see line 471: typical positional error in Gaia DR3 match is 2.5 
elif catalog == 'UVOTSSC2':
   rootobs    = ROOT+'data/catalogs/uvotssc2/' # '/disk/xray16/npmk/uvotssc2/'
   chunks     = []
   for k in range(24):
       chunks.append(f'ovotsso{k:02}.fits')
   match_par  = 0.45   # see line 471: typical positional error 
elif catalog == 'test':   
   rootobs    = ROOT+'data/catalogs/uvotssc2/'
   chunks     = ['test_sources.fits']
   match_par  = 0.45   # see line 471: typical positional error 
elif catalog == 'test2':  # do not change "test2" as it is hardcoded in this module
   rootobs    = ROOT+'data/catalogs/tests/'
   chunks     = ['test_cat.fit']
   match_par  = 2.5   # see line 471: typical positional error 
      
else: raise IOError(f'Error : the catalog {catalog} is not defined in the code.')


#### END GLOBAL DEFINITIONS


class finalise(object):
   """
   take the products from the earlier processing which are files by epoch
   that contain an extension for matches to the Gaia DR3 PM>10mas/yr data, 
   and those with no match found there. 
   
   The nomatch records need to be matched first to Vizier using cdsskymatch
   in STILTS, so we have the Gaia DR3 designation and sourceid for those that
   have a match in Gaia but no PM of consequence. 
   
   Thereafter, populate the positions for epoch 2016 and 2000, and combine
   all these table bits into one, split over sky location (TBD)
   
   Finally, add the appropriate summary table bits from the original catalogue.
   
   currently still editing (@npmkuin,25-march-2023)
   """
   
   def __init__(self,epochfilesdir=epochfilesdir):
       import os
       print (f"finalise step1a on processing the catalogue: fix positions to Ep2000 & add some GaiaDR3 data")
       self.indir = epochfilesdir
       self.matchout = work+'/matchout_all.fits' # renamed cat_out/cat_out_all_pm.fits
       self.nomatch = work+'/nomatch_all.fits'   # renamed cat_out/cat_out_all_nopm.fits
       self.outfile = work+'/matched2gaiadr3.fits' # final product
       # further match using CDS and 2.5" -> cat_out_nopmxgdr3.Vizier1355.fits
       self.infiles = []
       self.infiles_match =""
       self.infiles_nomatch = ""
       x = os.listdir(self.indir)
       for a in x:
           if a[:8] == 'cat_out_':  
               self.infiles.append(a)
               if a[-7:-5] == "-1":
                   self.infiles_match   += a+" "
               elif a[-7:-5] == "-2":    
                   self.infiles_nomatch += a+" "
               else:
                   raise RuntimeError(f"match2gaiadr2_ *finalise.init(): I found file {a}!!! but there should not be a cat_out_ file that is not -1 or -2 type.")    
       if (self.infiles_match=="" )| (self.infiles_nomatch==""):
           raise RuntimeError(f"match2gaiadr3_* finalise.init: no files found infiles_match={self.infiles_match} infiles_nomatch={self.infiles_nomatch}")        
       self.cleanup = "rm  epoch_slice_cat.fits gaiadr3_obsid_tiles.fits gaiadr3_selected_4epoch.fits tmp_aaa_1.fit"


   def stilts_process(self):
       # match a nomatch table to gaiadr3 using cdsskymatch 
       # retain the relevant columns
       self.stilts_merge()
       print ("\n merged the cat_out_* files \n")
       self.match2gaiadr3()
       print ("\nmade match to gaiadr3 by cdsskymatch\n")
       # now the nomatch cat lacks some columns so, we match allowing them to be added
       # different columns in 1 = "designation source_id ra2016gaia ra_error dec2016gaia dec_error "+\
       # " parallax parallax_over_error pm pmra pmra_error pmdec pmdec_error healpix "+\
       # " raObs decObs obsEp  "
       import os
       #command = f"/bin/csh ~/bin/stilts_fin_match2gaia.sh {self.matchout} {self.nomatch} {self.outfile}"
       #status=os.system(command,30*"-=")   
       #if status > 0: print(f"finalise.stilts_process 224 Error in {command} ")
       
       print (f"\nthere are now: a file matched with gaiadr3_pmgt25 tmp_{self.matchout},\n")
       print (f"one for nomatch with pm>25 tmp_{self.nomatch} and \none all matched to gaiadr3 from cds: sussxgaiadr3_ep2000.fits \n")
       
       print ("deleting temporary files")
       status=os.system("rm tmp_aaa* tmp_gaia_pm_m*")
       status=os.system(self.cleanup)
       
       
   def stilts_merge(self):
       # merge the list of tables using stilts tcat     
       import os
       # since the previous processing sometimes add some columns for GroupSize and GroupID, 
       # for now we make sure the match extension is present, and keep columns 1-138 only
       #     - this was solved - 
       """
       from astropy.io import fits       # obsolete - replaced by scripts
       from astropy.table import Table
       filmat = ""
       filnomat = ""
       for fil in self.infiles_match:
           filout  = fil[:-5]+"-1.fits"
           filout1 = "tmpfil-1.fits"  # temp first extension (may have duplicate col names!
           filout2 = fil[:-5]+"-2.fits"
           x = fits.open(fil)
           try:
               t = Table (x["PMMATCH"].data)
               filmat += filout+" "
               t.write(filout1,format='fits',overwrite=True)
               command=f'java -jar ~/bin/topcat-full.jar -stilts tpipe in={filout1} out={filout}  cmd="keepcols $1-$138"'
               status = os.system(command)
               if status > 0: print(f"finalise.stilts_merge 229 Error in {command} ")
           except: 
               print (f"finalise.stilts_merge 231 - failure with {fil}[PMMATCH] ")    
           try:
               t2 = Table( x["NOMATCH"].data)
               filnomat += filout2
               t2.write(filout2,format='fits',overwrite=True)
           except:
               print (f"finalise.stilts_merge 237 no NOMATCH extension in {fil}")    
           print (f"writing {filout}")
           
       command1 = f"java -jar ~/bin/topcat-full.jar -stilts tcat in='{filmat}' "+\
         f" ifmt=fits out={self.matchout}  ofmt=fits lazy=True"  
       print ("\nfinalise.merge match merge :\n"+command1+"\n")  
       status=os.system(command1)
       if status > 0: print(f"finalise.stilts_merge 249 Error in {command1} ")

       command2 = f"java -jar ~/bin/topcat-full.jar -stilts tcat in='{filnomat}' "+\
         f" ifmt=fits out={self.nomatch}  ofmt=fits lazy=True "               
       print (f"in finalise.stilt_merge commands: \n{command2}\n")
       status=os.system(command2)
       if status > 0: print(f"finalise.stilts_merge 251 Error in {command2} ")    
       """
       
       command1 = f"java -jar ~/bin/topcat-full.jar -stilts tcat in='{self.infiles_match}' "+\
         f" ifmt=fits out={self.matchout}  ofmt=fits lazy=True"  #+\
         #f" icmd='delcols GroupID' icmd='delcols GroupSize' icmd='delcols Separation' icmd='delcols 2000Ep' "+\
       print (f"in finalise.stilt_merge commands: \n{command1}\n")
       status=os.system(command1)
       if status > 0: 
           print(f"finalise.stilts_merge 263 Error in :\n{command1}\n\n ")
           raise RuntimeError()
          
                 
       command2 = f"java -jar ~/bin/topcat-full.jar -stilts tcat in='{self.infiles_nomatch}' "+\
         f" ifmt=fits out={self.nomatch}  ofmt=fits lazy=True " # +\
         #" icmd='delcols GroupID' icmd='delcols GroupSize' icmd='delcols Separation' icmd='delcols 2000Ep' "                
       print (f"in finalise.stilt_merge commands: \n{command2}\n")     
       status=os.system(command2)
       if status > 0: 
           print(f"\nfinalise.stilts_merge 271 Error in \n{command2}\n\n ")
           raise RuntimeError()
       
       
   def match2gaiadr3(self,):
       # 
       # since cds renames columns, we match it all again, so that we have a consistent 
       # set of columns. Note that a few of the earlier matches with high PM can switch 
       # to another Gaia source. 
       #
       import os
       
       command4 = f"/bin/csh ~/bin/stilts_fin_match2gaia.sh {self.matchout} {self.nomatch} {self.outfile}"
       print (f"finalise.match2gaiadr3 -- call stilts_fin_match2gaia.sh: \n{command4}\n")
       status = os.system(command4)
       if status > 0: print(f"finalise.match2gaiadr3 285 Error in {command4} ")
       
       
# end class finalise()       
       
   
class healpix_catalogue(object):

    def __init__(self, catalog, out="healpix_catalogue.fits", ext=1, hpcat=True):
        import numpy as np
        from astropy.io import fits
        from astropy.table import Table
        from cdshealpix import nested
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        self.cat = fits.open(catalog)
        self.ext=ext
        self.depth=7
        self.hpcat = hpcat
        self.t = Table.read(self.cat)
        if self.hpcat:
            self.colnames = self.t.colnames
            if not 'healpix' in self.colnames:
                raise RuntimeError(f"the input catalogue does not include a healpix column")
            else:
                self.healpix = self.t['healpix']    


    def add_healpixno(self):
        ra = self.t['ra'] 
        dec = self.t['dec'] 
        pos = SkyCoord(ra, dec , frame='icrs')
        self.ipix = nested.skycoord_to_healpix(pos,self.depth)
        self.t.add_column(self.ipix,name='healpix')
        self.t.sort('healpix')
        self.t.write(gaiaHP,overwrite=True)
       
    def get_records(self,pos):
        # central tile
        cpix = nested.skycoord_to_healpix(pos,self.depth)
        pixes = nested.neighbours(cpix,self.depth)
        #print (f"nested neighbours are {pixes}")
        return pixes
        """
        recs = []
        for p in pixes[0]:
            r = self.t[self.healpix == p]
            #print (f"number of sources in {p} is {len(r)}")
            for x in r:
                recs.append(x)
        return recs 
        """  
        
    def get_header(self):
        return self.t.colnames    

class split_catalogue(object):  # OBSOLETE
   """
   OBSOLETE 
   
   Split the catalogue into blocks of given size using ra, dec   
   provide method to query only central and nearby blocks

   example: 
       cat = split_catalogue( gaiacat, gaiapath, )
       cat.create()
       filenames = cat.selectNames( 210.0, -48.4 )
       
       filenames: list of paths to the sections of gaiacat around pointing 210.0,-48.4
       
       Problem near poles. use healpix!
       
       Duplication (same source found 24, 48, or 72 times)
          
   """
   def __init__(self, catalog, outdir, ext=1):
       import numpy as np
       from astropy.io import fits
       from astropy.table import Table
       self.cat = fits.open(catalog)
       # where to put the results, base name
       self.path = outdir+'/cat_'
       # blocks in RA (x1:x2) and Dec (y1:y2)
       self.x1 = np.arange(360)
       self.x2 = self.x1+1.
       self.y1 = np.arange(-80,80,1)
       self.y2 = self.y1 + 10.
       self.top = 80.0
       self.bottom = -80.0
       self.fntop = filename = self.path+f"_top.fits"
       self.fnbottom = self.path+f"_bottom.fits"
       self.ext=ext
       
   def create(self,):    
       # split the catalogue (sort RA, split, sort Dec, split )
       t = Table(self.cat[self.ext].data)
       t.sort('ra')
       for a,b in zip(self.x1,self.x2):
           qra = (t['ra'] >= a) & (t['ra'] < b)
           tra = t[qra]
           tra.sort('dec')
           for c,d in zip(self.y1,self.y2):
               qdec = (tra['dec'] >= c) & (tra['dec'] < d)
               tdec = tra[qdec]
               filename = self.path+f"{a:003}"+f"_{c:003}.fits"
               print (f"writing {filename}")
               try:
                   tdec.write(filename,overwrite=False)
               except: pass    
       # add a cap to the poles of the set to make sure of coverage        
       t.sort('dec')
       qdec = t['dec'] > self.top
       tdec = t[qdec]
       filename = self.fntop
       print (f"writing {self.fntop}")
       tdec.write(self.fntop,overwrite=False)
       qdec = t['dec'] < self.bottom
       tdec = t[qdec]
       filename = self.fnbottom
       print (f"writing {self.fnbottom}")
       try:
           tdec.write(self.fnbottom,overwrite=False)
       except: pass    
                              
   def selectNames(self,ra_pnt,dec_pnt):
       i1 = np.searchsorted(1.0*self.x1,ra_pnt +0.5,side='right') - 1
       j1 = np.searchsorted(1.0*self.y1,dec_pnt+0.999, side='right') - 1
       filenames = []
       for i in [i1-1,i1]:
          for j in [j1-1,j1]:
              if i < 0: i = i+360
              if i > 359: i = i-360
              if (j >= 0) & (j<18): 
                  a = self.x1[i]
                  b = self.y1[j]
                  fname = self.path+f"{a:003}"+f"_{b:003}.fits"
                  filenames.append(fname)
       if dec_pnt > self.top:    filenames.append(self.fntop )
       if dec_pnt < self.bottom: filenames.append(self.fnbottom)         
       return filenames
       
   def __exit__(self):
       self.cat.close()    
       

class pointing_of_observation(object):
    """
    Find pointing RA,Dec 
    """
    def __init__(self, catalogue, extension,multi=False,chunk=None):
        from astropy.table import Table
        from astropy.io import fits
        import numpy as np
        
        if not multi:
            t = fits.getdata(catalogue, ext=extension)
            t = Table(t)
            t.sort('OBSID')
            self.obsid = t['OBSID']
            self.ra =    t['RA_PNT']
            self.dec =   t['DEC_PNT']
            self.file = catalogue
            self.ext = extension
        else:
            #  run over the bit in chunk
            obs = []
            ra = []
            de = []
            for bit in chunk:
                #print (f"input {catalogue+bit}[{extension}]")
                t = fits.getdata(catalogue+bit, ext=extension)
                t = Table(t)
                obsids = t['OBSID']
                ras = t['RA_PNT'] 
                decs = t['DEC_PNT']
                for a,b,c in zip(obsids, ras, decs):
                    obs.append(a)
                    ra.append(b)
                    de.append(c)
            print (f"len obs:{len(obs)}; len ra:{len(ra)}; dec:{len(de)}  ")        
            self.obsid = np.asarray(obs)
            self.ra = np.asarray(ra)
            self.dec= np.asarray(de)
            self.file = f"multiple files {chunk} ext="
            self.ext = extension
        self.multi = multi


    def obspnt(self,obsid_in):
        """
    usage:
     p = pointing_of_observation(catalogue, 2)
     or 
     p = pointing_of_observation(catalogue, 2, multi=True, chunk=chunks)
     ra_pnt, dec_pnt = p.obspnt(obsid )
     finds the pointing for given obsid in the file
        """
        q = self.obsid == obsid_in
        if q.sum() == 0:
           raise runtimeError(f"the required OBSID pointing is not there {self.file}[{self.extension}]")
        ra1 = self.ra[q]
        dec1 = self.dec[q] 
        #print (f"pointing {ra1[0]},{dec1[0]}")
        if self.multi:
            return ra1[0], dec1[0]
        else:   
            return ra1.data[0], dec1.data[0]
            
            
    def get_data(self):
        print (f"table of poiting of observation\n ")
        for a,b,c in zip(self.obsid,self.ra,self.dec): 
            print (f"{a} {b} {c}")
        
                

def match2gaia_cat(root='.',catin=None,cat_matched=None,cat_nomatch=None,
        gaia_epoch=2016.0, obsEpochRange=[2005.0,2005.1],toEpochs=[2000.0],
        gaia_cat='gaiadr3_pmgt10.fits',outputf='merged_output.fits' ):
    """
    match the SUSS/UVOTSSC catalogue to Gaia DR3 astrometry for  obs_epoch in catalogue
    and then compute the positions for Epochs in list, typically the Gaia epoch, and 
    epoch 2000. 
    
    Input Parameters
    
    root : path
       path to the root directory for the catalogues
       work will be in the current directory   
    
    o1,o2 : int
       process records of input catin numbered from o1 to o2
    
    catin : str
       filename of input om/uvot catalogue (subset meant to be restricted in ObsEpoch)
       table in extension #1
    
    cat_matched, cat_nomatch: str
       filenames for the output of the records matched with Gaia PM and those without
       
    gaia_epoch: float 
       epoch of Gaia PM cat
       
    obsEpochRange: list
       range of Epochs in catin
       
    toEpochs: float 
       The epoch for which positions need to be added (in addition to the Gaia Epoch)
       
    outputf: float
       merged output file for this obsEpochRange   
    
    """
    
    print (f"\nstarting match2gaia_cat ... ")
    work=os.getcwd()
    
    if obsEpochRange[1]-obsEpochRange[0] > 1.0 :
       raise IOError(f"The range of the Epochs in the input file {catin} is more than a year! Aborting!")
    midEpoch = 0.5*(obsEpochRange[1]+obsEpochRange[0])
    gaiaout = "tmp_aaa_1.fits"
    matchout = "tmp_aaa_2.fits"
    nomatchout = "tmp_aaa_3.fits"
    gcatin = gaia_cat
    outputf1 = outputf[:-5]+"-1.fits"
    outputf2 = outputf[:-5]+"-2.fits"

    print ('1 precess')
    # precess gaiadr3_pmxxx to midEpoch -> (raObs,decObs):
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tpipe "+\
    f"in={gcatin} ifmt=fits out={gaiaout} ofmt=fits "+\
    f'cmd="addcol epxxxx epochProp({midEpoch}-{gaia_epoch},ra,dec,parallax,pmra,pmdec,0.)" '+\
    f"cmd='addcol raObs pick(epxxxx,0)[0];addcol decObs pick(epxxxx,1)[0]' "+\
    f"cmd='delcols 'epxxxx';addcol obsEp {midEpoch}'" 
    #print (command)
    status = os.system(command)
    if status > 0: print(f"560 Error in {command} ")

    print('2 match')
    # match 
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"tmatch2 in1={catin} ifmt1=fits  in2={gaiaout} ifmt2=fits "+\
    f"out={matchout} ofmt=fits-basic  matcher=skyerr "+\
    f"values1='RA DEC POSERR' values2='ra dec ra_error*0.001'  "+\
    f" params={match_par}   join=1and2 find=best2 fixcols=dups suffix1=  suffix2=2016gaia  "+\
    f" ocmd='delcols 'GroupID';delcols 'GroupSize' '"
#  f"ocmd='addcol ra2016gaia ra';'addcol dec2016gaia dec' "+\
#  f"ocmd=delcols 'ra dec raObs decObs obsEp' "
    #print('465 ',command)
    status = os.system(command)
    if status > 0: print(f"574 Error in {command} ")

    print ('3 nomatch')
    # nomatch
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"tmatch2 in1={catin} ifmt1=fits in2={gaiaout} ifmt2=fits  "+\
    f"out={nomatchout} ofmt=fits-basic  matcher=skyerr "+\
    f"values1='RA DEC POSERR'  values2='ra dec ra_error*0.001' "+\
    f" params={match_par}  join=1not2 find=best2 fixcols=dups suffix1=  suffix2=_2  "

    status = os.system(command)
    #print ('583 ',command)
    if status > 0: print(f"584 Error in {command} ")

    print ('4 matchout - add epoch position')
# edit the extra stuff out
#
#java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe  \
#in=tmp1.fits ifmt=fits out=$OUT ofmt=fits-basic \
#cmd="delcols 'OBSID_2 GroupID GroupSize'"
#   
    # add to cat_matched (matchout) the positions at gaia_epoch (already there) and toEpochs
    # precess gaiadr3_pmxxx to midEpoch -> (raObs,decObs):
    for toEpoch in toEpochs:
        matchout2="tmpx_"+matchout
        ra_ep = "ra"+str(toEpoch).split('.')[0]+'Ep'
        de_ep = "dec"+str(toEpoch).split('.')[0]+'Ep'
        to_ep = str(toEpoch).split('.')[0]+'Ep'
        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tpipe "+\
        f"in={matchout} ifmt=fits out={matchout2} ofmt=fits-basic "+\
        f'cmd="addcol epxxxx epochProp({toEpoch}-{gaia_epoch},ra,dec,parallax,pmra,pmdec,0.)" '+\
        f'cmd="addcol {ra_ep} pick(epxxxx,0)[0];addcol {de_ep} pick(epxxxx,1)[0]" '+\
        f"cmd='delcols 'epxxxx';addcol {to_ep} {toEpoch}' "
        print ("605 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
        
        os.system(f"mv {matchout2} {matchout}")
        print (f"output moved to {matchout}\n")

    print ('5 nomatchout copy RA,Dec to other epoch columns')
    # add to cat_nomatch a copy of positions for gaia_epoch and toEpochs
    for toEpoch in toEpochs:
        nomatchout2="tmpx"+nomatchout
        ra_ep = "ra"+str(toEpoch).split('.')[0]+'Ep'
        de_ep = "dec"+str(toEpoch).split('.')[0]+'Ep'
        to_ep = str(toEpoch).split('.')[0]+'Ep'
        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
        f" tpipe in={nomatchout} ifmt=fits out={nomatchout2} ofmt=fits-basic "+\
        f"cmd='addcol {ra_ep} RA;addcol {de_ep} DEC;addcol {to_ep} {toEpoch}' " 
        print (command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")

        os.system(f"mv {nomatchout2} {nomatchout}")
        print (f"output moved to {nomatchout}\n")

    #print ('6 merge into fits file with 2 extensions') obsolete
    #
    #  add fits extension names
    #
    try:
        add_extname(f"{matchout}","PMMATCH")
    except:
        print (f" match2gaiadr3_ 631 failed to update extname for {toEpoch} PM match ") 
        pass
    try:    
        add_extname(f"{nomatchout}","NOMATCH")
    except: 
        print (f" match2gaiadr3_ 631 failed to update extname for {toEpoch} No PM match ") 
        pass    
    # merge matchout and nomatchout obsolete
    #   write separate files for PM matches and non-matches
    status = os.system(f"mv {matchout} {outputf1}")
    status = os.system(f"mv {nomatchout} {outputf2}")
    """ obsolete
    command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
    f"tmulti in='{matchout} {nomatchout}' ifmt=fits "+\
    f" out={outputf} ofmt=fits-basic"  
    print (f'303 {command}')
    status = os.system(command)
    if status > 0: print(f"238 Error in {command} ")
    print (f"306 end matching output to {outputf}\n")
    """
    
def add_extname(file,extname):
    from astropy.io import fits
    with fits.open(file,mode="update") as f:
       f[1].name=extname
       f.flush()
     
    
def editObscat(rootobs,work,TOPCATPATH):
# assume the catalogue has two extensions, SRCLIST   and SUMMARY

#  edit input catalog, e.g., - add obsEpoch to catalog  if missing
  for bit in chunks:
    print (f"editing - add epoch to {bit} ")
    h = fits.open(rootobs+bit)
    temp1 = work+'/_ori_'+bit
    temp2 = work+'/tmpobs.fits'
    temp3 = work+'/tmpobs2.fits'
    temp4 = work+'/summary.fits'
    xtname = 'SOURCES'
    if (catalog == "SUSS") | (catalog == "test2"): xtname = 1 # 'SRCLIST' or 'Joined'
    if not 'obs_epoch' in h[xtname].data.names:
        # passed check not yet edited... 
        
        command=f"cp {rootobs+bit} {temp1}"  # copy SUSS to temp1
        print ("\neditObscat:561 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
        
        # copy the summary extension to its own file  SUSS-SUMMARY to temp4
        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tpipe "+\
        f"in={temp1+'#2'} ifmt=fits out={temp4} ofmt=fits-basic"
        print ("\n331 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
       
        # edit the dates in the uvotssc2 input file to MJD
        if (catalog == 'UVOTSSC2'):  # copy cat-SRCxxx to temp2, add MJD  
           command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tpipe "+\
    f"in={temp1+'#2'} ifmt=fits out={temp2} ofmt=fits-basic "+\
    f"cmd='addcol -units d MJD_START isoToMjd(DATE_MIN)' "+\
    f"cmd='addcol -units d MJD_END isoToMjd(DATE_MAX)' "+\
    f"cmd='delcols DATE_M*' "
        else: 
           command=f"cp {temp4} {temp2}" # copy SUSS-SRCxxx to temp2
        print ("\neditObscat587 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
                
        print (f"\nnow recombine with summary ")
        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
         f"tmulti in='{temp1+'#1'} {temp2}' ifmt=fits out={temp3} ofmt=fits-basic"  
        print ("\n349 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
        
        # add the epoch 
        command=f"addEpoch2Source.sh {temp3} {temp1}"
        print ("355 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
                
        #print (f"\n recombine result once more with summary ")
        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tmulti "+\
         f"in='{temp1} {temp4}' ifmt=fits out={rootobs+bit} ofmt=fits-basic"  
        print ("\n362 ",command)
        status = os.system(command)
        if status > 0: print(f"Error in {command} ")
    h.close()
    
def getEpochObscat(ep1,ep2,rootobs,obs_cat):    
    """
    write obs_cat for the requested Epoch and return list of OBSIDs
    """
    print ("368 ",ep1, ep2, obs_cat)
    keepcol = ['IAUNAME', 'OBSID', 'SRCNUM', 'RA', 'DEC', 'POSERR', 'obs_epoch']
    t = []
    ndata = 0
    obsids = []
    gotdat = True
    for bit in chunks: 
        print (f"processing {bit}")
        d = fits.open(rootobs+bit)
        td = Table(d[1].data)
        if not ((catalog == "SUSS") |(catalog == "test2")) :
            td.add_column(np.max([td['RA_ERR'],td['DEC_ERR']]),name='POSERR')
            dd = td[keepcol]
        else: dd= td    
        dd.sort('obs_epoch')
        q = (dd['obs_epoch'] >= ep1) & (dd['obs_epoch'] < ep2)
        print ("380 total number of records in epoch found is ",q.sum())
        ndata += q.sum()
        if q.sum() > 0:
           t.append(dd[q])
           print (f"adding {q.sum()} epoch records for {bit} ")
        d.close()   
    
    if len(t) > 0:       
       tab = vstack(t,join_type='exact')    
       tab.write(obs_cat,overwrite=True)
       print (f"388 now the vstacked catalog {obs_cat} has been made.\n")
    else:
       print (f"390 no data in epoch range [{ep1},{ep2}]\n")
       gotdat = False
    
    # now process the Gaiadr3 to obtain a proper subset for matching, subset name is 
#    gaia_temp_cat = work+"/gaiadr3_selected_4epoch.fits" 

    # list of unique obsids in this epoch slice    
    print (f"397 sorting for unique list of obsids")
    if gotdat:
        ox = tab['OBSID']
        ox.sort()
        x = ox[0]
        obsids = [x]
        for k in ox:
           if k != x:
               obsids.append(k)
               x = k
    print (f"found {len(obsids)} unique obsids") 
    print (f"found {ndata} rows in this epoch")  
    return obsids
    
 
def  processGaiaPM(rootobs,chunks=None, xcat=None, obsids=None, chatter=1):
    """
    old processing based on using RA, Dec tiles
    """
    print (f"412 processing Gaiadr3_pm catalogue")
    t1 = time.time()
    # select gaia subset for these obsids 
    tempfile = work+'processgaiapm_tmp.fits'
    gaia_sel = []
    n_obsid = 0
    gfiles=[]
    allrecs = []
    if catalog == 'UVOTSSC2':
        pointing = pointing_of_observation(rootobs,2,multi=True,chunk=chunks)
    for bit in chunks:
        n_bit = 0
        if catalog != 'UVOTSSC2':
            pointing = pointing_of_observation(rootobs+bit,2)
        for ox in obsids:
            ra_pnt, dec_pnt = pointing.obspnt(ox)
            if np.isfinite(ra_pnt):
                n_obsid += 1
                if not use_healpix:
                    filelist = xcat.selectNames(ra_pnt,dec_pnt)
                    #print (f"427 processGaiaPM pnt=[{ra_pnt},{dec_pnt}]")
                    for gaiafil in filelist:
                        gfil = gaiafil.split('.')[0].split('/')[-1]
                        tmp = work+f"/tmp_{gfil}_{n_bit}_{n_obsid:00005}.fits"
                        n_bit += 1
                        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
                      f' tpipe in={gaiafil} ifmt=fits out={tmp} ofmt=fits '+\
                      f' cmd="select skyDistanceDegrees(ra,dec,{ra_pnt},{dec_pnt})<{RADIUS}"'
                        if chatter > 3: print (command)
                        status = os.system(command)
                        if status > 0: print(f"416 Error in {command} ")
                        # remove empty data files
                        ft = Table.read(tmp)
                        if len(ft) > 0:
                            gfiles.append(tmp)
                            #if chatter > 2: print (f"adding {tmp} -- {n_obsid}:{ox}  {n_bit}:{bit}")
                        else:
                            #if chatter > 2: print (f"no data in match for {tmp} - skip")  
                            os.system(f"rm {tmp}")  
    if not use_healpix:       
      # combine gfiles stuff    
      #print (f"487 {len(gfiles)} gfiles: {gfiles}  \nn_obsid = {n_obsid}\n")
      print (f"573 number of files = {len(gfiles)}\nn_obsid = {n_obsid}\n")
      xxx = ""
      gaia_k = 0
      # sum 50 at a time:
      nt = int(len(gfiles)/50) 
      for nk in np.arange(nt+1):
        n1 = 50 * nk
        n2 = np.min([50 * (nk + 1),len(gfiles)])
        #print (f"summing nk={nk}; n1,n2=[{n1},{n2}]")
        if gaia_k > 0:
           os.system(f"cp {gaia_temp_cat} {tempfile}") 
           xxx = tempfile+" "
        for fil in gfiles[n1:n2]: 
            xxx = xxx+fil+" " # this is a list of fits files to stack
            if chatter > 4: print (f"files to stack: {xxx}")
        command = f"java -jar {TOPCATPATH}/topcat-full.jar -stilts "+\
            f"tcat in='{xxx}' ifmt=fits out={gaia_temp_cat} ofmt=fits-basic"  
        if chatter > 3: print (f"{command}\n")
        status = os.system(command)
        if status > 0: print(f"454 Error in {command} ")
        gaia_k = 1
        
      if chatter > 4: 
        print (f"the full file list is: {gfiles}")
        print (f"check that also the last entries have been processed into the {gaia_temp_cat}")         
      t2 = time.time()
      print (f"Elapsed time subsetting Gaia = {t2-t1}\n")
 
def  processGaiaPM2(rootobs,chunks=None, xcat=None, obsids=None, chatter=1):
    """
    Processing based on Healpix tiles
    """
    print (f"412 processing Gaiadr3_pm catalogue")
    t1 = time.time()
    n_obsid = 0
    tiles = []  # list of tiles
    obsidTiles = []
    if catalog == 'UVOTSSC2':
        pointing = pointing_of_observation(rootobs,2,multi=True,chunk=chunks)
    for bit in chunks:
        n_bit = 0
        if catalog != 'UVOTSSC2':
            pointing = pointing_of_observation(rootobs+bit,2)
        for ox in obsids:
            ra_pnt, dec_pnt = pointing.obspnt(ox)
            if np.isfinite(ra_pnt):
                row = [ox]
                n_obsid += 1
                if type(ra_pnt) == Quantity:
                        pos = SkyCoord(ra_pnt,dec_pnt,frame='icrs')
                else:
                        pos = SkyCoord(ra_pnt*u.deg,dec_pnt*u.deg,frame='icrs')
                pixes = xcat.get_records(pos)[0]
                for p in pixes:                    
                    row.append(p)
                    tiles.append(p)
                if len(pixes) < 9:
                    for p in range(len(pixes)-9):
                        row.append(None)
                        tiles.append(None)
        obsidTiles.append(row)        
    onames=['obsid','pix1','pix2','pix3','pix4','pix5','pix6','pix7','pix8','pix9']   
    obsidTiles = np.asarray(obsidTiles)
    #print (obsidTiles.shape) 
    obsidTable = Table(data=obsidTiles,names=onames,)    
    obsidTiles = None
    # now get list of unique tiles in tiles    
    print (f"before removal multiples number of tiles is {len(tiles)}")           
    tiles.sort()
    yr =[]
    yr.append(tiles[0])
    a = tiles[0]
    for b in tiles[1:]:
         if b != a:
             yr.append(b)
             a = b
    tiles = yr  # unique tiles 
    print (f"after removal multiples number of tiles is {len(tiles)}")       
    # 
    # get records in xcat.t for unique tiles and precess them -> new table tp
    #
    allrecs = []
    for p in tiles:
        x = xcat.t[xcat.healpix == p]
        allrecs.append(x)  
    tab = vstack(allrecs) 
    # for further processing we need both the table connecting obsid to tiles and
    # the table of all tiles to be precessed
    #
    obsidTable.write(gaiaObsidTable,overwrite=True) 
    tab.write(gaia_temp_cat,overwrite=True)
    t2 = time.time()
    print (f"Elapsed time subsetting Gaia = {t2-t1}\n")

  
def mainx():
    # initialise xcat
    print (f"initialising the Gaia catalogue bits")
    if use_healpix:
        try:  
            xcat = healpix_catalogue(gaiaHP, None, ext=1 )
        except:
            xcat = healpix_catalogue(gaiacat, gaiaHP, ext=1 )
            use_helpix = False
    else:
        xcat = split_catalogue(gaiacat, gaiapath, ext=1 )
    
    print (f"adding epoch to input cat?")    
    if do_epoch: 
        editObscat(rootobs,work,TOPCATPATH)
        print (f"\n \033[1m{'done editing input catalogue by adding obs_epoch'}\033[0m\n")
    # possibly, check that it worked and they all have obs_epoch

    # set list of chunks (ep1,ep2) given as parameter input, e.g. all the 'hour' files
    epoch_array = np.arange(process_epoch_range[0],process_epoch_range[1],process_epoch_step)
    # so, we want to get epochs from epoch_array[k] to epoch_array[k]+process_epoch_step
    n_epochs = len(epoch_array) - 1

    # collect all the data for epoch number k and strip out every column
    #    except IAUNAME, OBSID, SRCNUM, RA, DEC, obs_epoch

    merged_output_file_list=[]   

    for k in np.arange(n_epochs):
        ep1 = epoch_array[k]
        ep2 = ep1 + process_epoch_step
        t1 = time.time()

      # process our catalog to get a subset for this epoch
      # obs_cat = work+'/epoch_slice_cat.fits'
    
        obsids = getEpochObscat(ep1,ep2,rootobs,obs_cat)    
        if len(obsids) == 0: 
            continue    
  #
  # at this point I have a problem - match against the whole Gaia pm>10 catalogue 
  # or first make a subset of that catalogue for all the observed fields and then 
  # match ? 
  # matching one obsid takes several minutes. We have hundreds of thousand 
  # obsids so of order 1e5+ s goes into this subsetting ... 
  # split the gaia catalog in pieces and then limit the search to those for each pointing 
  #
        if use_healpix:
            processGaiaPM2(rootobs,chunks=chunks, xcat=xcat, obsids=obsids,chatter=1) 
            # was (chunks,rootobs,obsids,gaiacat,gaia_temp_cat,TOPCATPATH,RADIUS,xcat)

        else:
            processGaiaPM(rootobs,chunks=chunks, xcat=xcat, obsids=obsids, chatter=1) 
            # output is catalogue for all gaia tiles for this epoch

  # now we can do the matching 
        t2 = time.time()
        merged_out = f"cat_out_{ep1:06.1f}.fits"
        match2gaia_cat(root='.', catin=obs_cat ,cat_matched='matched.fits',\
          cat_nomatch='nomatch.fits', gaia_epoch=GaiaEpoch, \
          obsEpochRange=[ep1,ep2],toEpochs=[2000.0], gaia_cat=gaia_temp_cat,\
          outputf=merged_out )
     
        print (f"for epoch {ep1} the results is in file {merged_out} \n removing temp files\n")
        os.system(f"rm -f tmp_cat_*.fits")
        merged_output_file_list.append(merged_out)
        t3 = time.time()
        print (f"elapsed time matching to Gaia = {t3-t1}")

  #  clean up any tmp files left
        #try:
        #   os.system(f"rm  tmp_cat_1*.fits ")     
        #   os.system(f"rm  tmp_cat_2*.fits ")     
        #   os.system(f"rm  tmp_cat_*.fits")     
        #except: pass

   
