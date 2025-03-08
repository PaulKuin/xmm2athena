#   python3
#
# for merging the high pm sources
# this program is stage 2 of the single source cat making 
# the single_source_cat_v2.py (May 2023) should first run
#   also further editing of the failed masks is needed.
# only needs to be run on matched Gaia dr3-UVOT sources
#
# made by using the diff from the v1 files
#
#    @npmkuin, 2023, Copyright 2023, license 3-clause BSD style (see github.com/PaulKuin)
#
#  This version (V2) computes the chiSquared values for variability.
#  stage 2, June 19, 2023, npmk
# March 2025, npmk
# added cluster analysis for high PM sources in order to identify which sources have 
# same PM within a limit, and are in the same area on the sky.

import os
import numpy as np
import numpy.ma as ma
import astropy, numpy
from astropy.io import fits, ascii as ioascii
from astropy.table import Table, Column, MaskedColumn

# updates: return list of obsids, epochs [can be used to search for original records)
#  return median, skew, etc.
# not yet: ABmag upper limits from SUMMARY 

########################################################################################
# globals
# define globals:

# output filename : see fileio() how it is being derived
fix_duplicate = False # do not set! needs a lot of eval to get it to work - 
# the start/stop times are for agregated images for one obsid, with often 
# multiple image snapshots being merged 
colsofinterest = ["IAUNAME","SRCNUM","obs_epoch","OBSIDS","EPOCHS",'UVW2_ABMAG','UVM2_ABMAG','UVW1_ABMAG','U_ABMAG',"RA_ERR","DEC_ERR","UVW2_NOBS","UVM2_NOBS","UVW1_NOBS","U_NOBS"]

iauname="IAUNAME" 
bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent

catalog = "SUSS" #'testsuss' # 'SUSS' 'UVOTSSC2'
outtable = None
matched2gaia = "sussxgaiadr3_ep2000.fits"
Rpm_limit = 10 # if PM/error < Rpm_limit the record will not be merged
chatter = 0


if catalog == 'SUSS':
    rootobs = '/Volumes/DATA11/data/catalogs/test_2/'
    chunk = "sussxgaiadr3_ep2000_singlerecs.fits" 
    pm_cds = "PM_cds"  # proper motion 
    Rplx = "RPlx"  # Parallax over Error
    onedone = False # if the sources without large pm have been written already
    hipm_chunk = "sussxgaiadr3_ep2000_singlerecs_hipm.fits"
elif catalog == 'UVOTSSC2':
    rootobs = '/Users/kuin/pymodules/cats/' # /Users/data/catalogs/uvotssc2/' 
    chunk =  "TBD x gaia" # 'cat_out_all_pm.fits' 
elif catalog == 'test':
    rootobs = '/Users/data/catalogs/uvotssc2/'
    chunk = 'test_sources.fits'
elif catalog == 'testsuss':
    rootobs = '/Volumes/DATA11/data/catalogs/test/'
    chunk = "sussxgaiadr3_ep2000_singlerecs.csv"  
    pm_cds = "PM_cds"  # proper motion 
    Rplx = "RPlx"  # Parallax over Error

tmax = int(5e8) # mximum number of records to read in (adjust for test)

print (f"checking proper motion in col {pm_cds} and using parallax over error column {Rplx}")

################ END GLOBAL SECTION ######################################################


def make_chunks(chunk, dn=3020000):
    """
      the UVOTSSC input file is >54GB and too large to have in memory, so we need to process in 
      parts: split up the file by iauname - it has 159,845,628 records, so 53 slices 
      of 3M records
    """
    chunklist = []
    f = fits.open(chunk,memmap=True)
    tab = f[1].data
    nrow = len(tab)
    print (f"{nrow}")
    nch = nrow//dn
    n1=0
    n2=dn
    go1 = True
    testname = tab[iauname]
    idiot = 0
    while go1:
      for k in np.arange(nch+1): 
      
         print (f"{n1} - {n2}   {testname[n1]}-{testname[n2]}")
         outf = f"chunk_{k:02}.fit"
         t1 = testname[n2]
         st = True
         while st:
            n2+=1 
            if testname[n2] != t1:
                st = False 
         tchunk = Table(tab[n1:n2])
         tchunk.write(outf)
         chunklist.append(outf)
         n1=n2
         n2=n1+dn
         if n2 >= nrow:
             outf = f"chunk_{k+1:02}.fit"
             tchunk = Table(tab[n1:nrow])
             tchunk.write(outf)
             chunklist.append(outf)
             go1 = False
             break
         idiot+=1 
         if idiot > 60: 
             raise RuntimeError()  
              
    # 
    f.close()
    del f[1].data
    # now process them 
    for chunk in chunklist:
       # use ftools to append the summary extension to the chunck for further processing
       summ = "uvotssc2_in_summary.fits"
       command = f"ftappend {summ} {chunk} "
       os.system(command)
       # make single object catalogue
       mainsub(chunk)
       
def chunklist():
    chunklist = []
    nch = 54
    for k in np.arange(nch+1):
       chunklist.append( f"chunk_{k:02}_singlerecs_v2.csv" )
    return chunklist
           
   
def make_new(tx):
    return tx
    #   OBSOLETE:
    # remove and add columns to main table
    tx.remove_columns(['N_SUMMARY','N_OBSID'])
    for  b in bands:
        tx.remove_columns(["OBSID","obs_epoch"])
        tx.add_columns( [None,None,None,None],names=['IAUNAME2','EPOCHS','OBSIDS','SRCNUMS']  ) 
        # NOTE: REFID = number of filters observed
    #tx.add_columns([None,None],names=['OBSIDS', 'EPOCHS'])           
    return tx         

def stats(array, err=[None], syserr=0.005, nobs=None, chi2=None, sd=None, skew=None):
    """
    input parameters for merging second stage
    ----------------
       array: 1D array
       err: 1D array 
          elements correspond to array elements
       syserr: float
          syserr is the non-ramdom systematic error estimate  (default 0.005) 
    
    Median:
       sort the 1D array in order
       select the middle element(s)

    chi-squared:
       one value from mean

    skewness:
    
    
    kurtosis:
    

    Returns:   
       
       return values of 
         the number of points, 
         median, 
         chi-squared, 
         standard deviation
         skewness
         variability measure V3
         
         In addition, a triangle-shaped filter used in wavelet analysis can determine 
         flaring (fast rise, slow decay). Probably needs 32+ data points. See the 
         paper by Marcus Aschwanden et al, 1998. 
         (check if there is a paper by Curtis Saxton on this also).

    @npmkuin, 2023, Copyright 2023, license 3-clause BSD style (see github.com/PaulKuin)
    note the input array, err can have null elements
    
    """
    import numpy as np
    from scipy.stats import kurtosis
    from scipy.stats import skew
    
    # first remove nan values and if no elements left, return, otherwise continue
    #n1 = len(array) input can be a float
    a = np.asarray(array)
    q = np.isnan(a)   # to screen out missing values 
    a = a[q == False]
    na = len(a)
    if len(a) < 1: 
        return 0, np.nan, np.nan, np.nan, np.nan, np.nan
    if err.all() == None:    
       e = syserr
       var = np.ones(len(a))*(0.1+e)**2 # assumed photometric error 0.1 mag if missing
       # should multiply var with N/(N-1)
       var = var*na/(na-1)
    else:
       err = err[q == False]
       e = np.asarray(err)+syserr
       var = e*e
       # should multiply var with N/(N-1)
       if na >1: var = var*na/(na-1)
       else: var = None
    if a.ndim > 1: 
       raise IOError(f"median: ithe input is not 1D")
    n = len(a.flatten())  # number for single filter
    a.sort()
    e.sort()
    if n == 1:
       return 1, a[0], np.NaN, e[0], np.nan, np.nan
    if n/2*2 == n:
       median = 0.5*( a[int((n-1)/2)]+ a[int((n+1)/2)] ) 
       sd = 0.5*( e[int((n-1)/2)]+ e[int((n+1)/2)]  )
    else:
       median = a[int(n/2)]     
       sd = e[int(n/2)]  
    # expected standard deviation sd if normally distributed would depend on the median error      
    p = a.mean()
    chisq = (a-p)*(a-p) 
    chisq = chisq.sum() / p 
    # non-central chi-squared distribution (Wolfram)
    ncChisq = 0.5*(a-p)*(a-p)/var
    ncChisq = 2.*ncChisq.sum() 
    # variability3 (x_i - x_k)/(3err_k) with k at median error sd 
    V3 = (a - median)/(3*sd)
    V3 = np.max(V3)
    #kur = kurtosis(a,nan_policy='omit')
    sk = skew(a,nan_policy='omit')
    
    return n, median, ncChisq, sd, sk, V3   
      

def fileio(infile,outstub="_stg2",outdir=""):
    #
    # use the input file name infile to construct the output filename
    # open the output file
    # read in the table from the input file
    # return the table and the output file handle
    a = infile.rsplit('.')
    instub = ""
    ft = a[-1]  # file extension 
    for x in a[:-1]: 
        instub += x+"."
    instub = instub[:-1]    
    print (f"reading {infile} ...\n")
    
    if (ft.lower() == 'csv'):
       print (f"reading input csv file ")
       t = ioascii.read(infile)  # table object
    elif (ft[:3] == 'fit'):
       print (f"reading input fits file")
       f = fits.open(infile)
       x = f[1].data #[:tmax]  # LIMITED TABLE FOR TESTS
       
       #summ = f[2].data # summary - has obstimes and exposures of obsids
       t = Table(x)   
       #summ = Table(summ)
       x = ''
       f.close()
    else:
       print("problem")
       raise RuntimeError("Fatal Error: *** input file was not read ***")   

    #print (f"the crashes happen next, so print all the keys to t: {t.colnames}")
    #t.add_column(np.sqrt(t["gaiadr3_pmra"]*t["gaiadr3_pmra"]+t["gaiadr3_pmdec"]*t["gaiadr3_pmdec"]),name="gaiadr3_pm")
    if len(t) < 1:
        raise IOError("the catalogue table loading failed.\n")  
    
    print (f"total number of input sources before merging high PM ones is {len(t)}")    
    # now just limit the processing to high PM entries   (tx) and leave the rest (ty)  
    try:     
        tx = t[ t[pm_cds] > 20. ]  
        ty = t[ (np.isnan(t[pm_cds]) ) | (t[pm_cds] <= 20.)  ]
    except:
        print ("279 problem reading table -- colnames are \n {t.colnames}")    
        exit            

    outfile=instub+outstub+"a.csv"
    if not onedone: 
       outf = open(outdir+outfile,'w')  # file handle for all records that do not need another merge
    else:   
       outf = "outf_296"   
    outfile2=instub+outstub+"b.csv"
    outf2 = open(outdir+outfile2,'w')  # for new merges 
    outerr = open(outdir+instub+"_err.csv",'w')  #was:file handle weird ones with same PM but different parallax_over_error
    # 
    # we now need to group the sources according to close position and PM
    #
    # distance: best is distance < 0.001xPM arcminutes
    # delta PM < 1e-2 
    # 
    # method: index of group after sort by PM, then sort by position, then check membership 
    #
    # after that we can call up the record/row selection by that.
    #
    #tx.add_index(pm_cds)
    sources = tx[pm_cds]
    sources = np.unique(sources)
    
    print (f"NUMBER OF UNIQUE SOURCES : low PM:{len(ty)} high PM: {len(sources)}\n")
    try:
        del f[1].data
        del f[2].data
    except: pass
    return tx, outf, sources, ty, outf2, outerr

   

def evaluate_extended_nature(tsrc):
    # compare the extended flags for long exposures, and consistency in the 
    # various bands 
    # input for each observation in target
    #
    # are more than half of observations in a filter extended ?
    #  check with the same PA?
    #  combine for all bands ?
    #  limit to with sufficient exposure
    #
    # => was the observation trailing? 
    # => transient?
    # => consistent?  
    #
    # simplest solution
    #    obsids = tsrc['OBSID']
    #    tsumm = summ.loc[obsids] # find records for obsids - needs first to summ.add_index('obsid')
        # refer to exposure by tsumm['EXPOSURE_UVW2'] etc. 
    # need check all
    bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent
    base = {'UVW2':255,'UVM2':255,'UVW1':255,'U':255,'B':255,'V':255}
    for band in bands:
       if (catalog == 'SUSS') | (catalog == "testsuss"):
           exflag = tsrc[band+'_EXTENDED_FLAG']
       elif (catalog == 'UVOTSSC2') | (catalog == "test"):
           exflag = tsrc[band+'_EXTENDED']
       q = exflag != 255
       if np.sum(q) > 0:
           base[band]= int(np.mean(exflag[q]) > 0.5) # more than half say extended
    # so valid bands have base not equal to 255 
        # TBD:
    # updates to the base solution:
    # case 1: single obs. in one band: pass the info
    # case 2: multiple obs, but after weeding only 1 : pass merged info
    # case 3: all bands extended in long exposures : OK extended
    # case 4: all bands not extended in long exposures: OK point source
    # case X: something else ...
     
    return base
    

def create_csv_output_record(trow,outf,col,nc):
    # create the output record to write to a csv file 
    # input is a table object for one source after the magnitudes have been 
    # combined, and the extended nature has been evaluated
    #
    #print     ( f"in create_csv_output 357 : type:{type(trow)} ... values ...\n{trow},\n---" )
    for k in range(nc-1): 
       outf.write( f"{trow[col[k]]}," )
    outf.write( f"{trow[col[nc-1]]}\n" )

    
# The set flags are ranked as follows (rank:flag bit) (bad->worse->verybad order):
f_rank = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 6:2, 7:0, 8:3, 9:10,10:9, 11:11}
inverse_f_rank ={0:7,1:5,2:6,3:8,4:2,5:1,6:3,7:4,8:0,9:10,10:9,11:11}
f_rk = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 7:0, 9:10,10:9, 11:11} # without f_carry

def qual_st_merge(qual_st,err=0.02):
    # produce a qual_ft (derive the integer flag from it    
    #print (f"input QUALITY_ST= {qual_st} ")
    
    # The set flags are ranked as follows (rank:flag bit) (bad->worse->verybad order):
    #f_rank = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 6:2, 7:0, 8:3, 9:10,10:9, 11:11}
    #inverse_f_rank ={0:7,1:5,2:6,3:8,4:2,5:1,6:3,7:4,8:0,9:10,10:9,11:11}
    #f_rk = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 7:0, 9:10,10:9, 11:11} # without f_carry
    # if one of the following flags is set in any of the source detections, carry it
    f_carry = [6,8] # original flags, i.e., [1,3] in ordered flags
    
    qualout = ""
    qual_st = np.asarray(qual_st)
    if qual_st.size == 1:  # nothing to decide
        return qual_st,  ranked_qual_st2qual(qual_st)
    else: 
       nq = 12  # number of bits
       qual = "FFFFFFFFFFFF"     # start with first one
       # check each flag is not set, select the lowest set
       if isinstance(qual_st,np.flexible):
           print (f"\nWARNING 352 (char) np.flexible : is {qual_st} just a string? \n\n")
           return qual_st,  ranked_qual_st2qual(qual_st)
       else:
          # convert to ranked and logicals
          ranked_quals = aa =  ranked_qual_st2qual(qual_st) #np.array(aa)
          if chatter > 3: 
              print (f"qual_st_merge 383 {60*'-'}\noriginal {qual_st}\ninteger flag values are: {aa}")    
              print (f"qual_st_merge385 summed rows: {aa}")
          qa = ma.min(aa) == aa
          k = ma.where(qa)[0]
          qualout = qual_st[k][0]
          if qualout.strip() == 'T':
              raise RuntimeError(f"398 return value error {qualout}\n input was {qual_st} err={err}")   
             
          return qualout , ranked_quals       
          
def qual_stflag2decimal(flag_st):
    if flag_st == "            ": return -2147483648
    out = 0
    for k in range(12):
        if flag_st[k] == "T":
           out += 2**k
    return out       

def qual_dec2qual_stflag(decimalflag):
    if decimalflag == -2147483648 : return "            "
    if decimalflag == -999 : return "            "
    if decimalflag == -32768 : return "            "
    a = str(bin(decimalflag))[2:]
    n = len(a)
    out = ""
    for k in range(12-n):
        out+='F'
    for k in a:
        if k == "0":
           out+='F'
        else: out+="T"    
    a = ""
    for k in range(11,-1,-1):
        a += out[k]        
    return a  


#f_rk = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 7:0, 9:10,10:9, 11:11} # without f_carry

def ranked_qual_st2qual(qual_st):
    # convert string qual to number qual after ranking with f_rk 
    if type(qual_st) == str: 
       if not (("F" in qual_st) | ("T" in qual_st)):
          print (f"single_source_cat_merge_v2.ranked_qual_st2qual426 qual_st input is {qual_st} in error")
       print (f"single_source_cat_merge_v2.ranked_qual_st2qual426 expected an array of quality ")   
    aa = []
    for k in range(qual_st.size):  # all rows
        x = ma.asarray(qual_st[k])
        if ~x.mask:
           bb = 0
           for j in f_rk.keys():
                s = str(x)
                if s[f_rk[j]] == 'T' :
                    bb += 2**(f_rk[j])
           aa.append(bb)
        else:
           aa.append(None)      # was -9999 instead of None
    return ma.masked_values(aa,None) # ditto
 
########## clustering 


def cluster_sources(ra, dec, pm_ra, pm_dec, id, pm_tolerance=0.001, pm_distance_factor=0.01, output_stuff=True,chatter=0):
    """
    Cluster sources by similar proper motion and positional proximity with dynamic distance threshold.
    Includes source IDs in the output.
    
    Parameters:
    ra, dec: arrays of right ascension and declination in degrees
    pm_ra, pm_dec: arrays of proper motion in RA and Dec directions (e.g., mas/yr)
    id: array of unique identifiers for each source (e.g., strings or integers)
    pm_tolerance: tolerance for matching proper motions (e.g., mas/yr)
    pm_distance_factor: factor to compute max separation (arcmin) as pm_distance_factor * PM
    output_stuff: optional 
    
    Returns:
    List of clusters, where each cluster is a dictionary with 'ids' (list of source IDs) and 'indices' (list of source indices)
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from sklearn.cluster import DBSCAN
    from sklearn.preprocessing import StandardScaler
    from sklearn.neighbors import BallTree
    import time
 
    start_time = time.time()
    n_sources = len(ra)
    if len(id) != n_sources:
        raise ValueError("Length of id array must match length of ra, dec, pm_ra, pm_dec")
    if chatter > 2: print(f"Processing {n_sources} sources...")

    # Step 1: Compute total proper motion for each source
    pm_total = np.sqrt(pm_ra**2 + pm_dec**2)  # Total PM in mas/yr

    # Step 2: Cluster by proper motion using DBSCAN
    pm_features = np.vstack((pm_ra, pm_dec)).T  # Shape: (n_sources, 2)
    pm_features_scaled = StandardScaler().fit_transform(pm_features)
    
    # Tune eps based on pm_tolerance
    pm_eps = pm_tolerance / np.std(pm_features_scaled[:, 0])  # Approximate scaling
    pm_dbscan = DBSCAN(eps=pm_eps, min_samples=2).fit(pm_features_scaled)
    pm_labels = pm_dbscan.labels_
    
    pm_cluster_count = len(set(pm_labels)) - (1 if -1 in pm_labels else 0)
    noise_count = np.sum(pm_labels == -1)
    if chatter > 2: print(f"Proper motion clustering done in {time.time() - start_time:.2f} seconds")
    if chatter > 2: print(f"Found {pm_cluster_count} PM clusters, {noise_count} sources labeled as noise")

    # Step 3: For each PM cluster, cluster by position with dynamic threshold
    clusters = []
    unique_pm_labels = set(pm_labels) - {-1}  # Ignore noise (-1 label)
    
    for pm_label in unique_pm_labels:
        pm_cluster_indices = np.where(pm_labels == pm_label)[0]
        pm_cluster_size = len(pm_cluster_indices)
        if pm_cluster_size < 2:
            continue
        
        # Compute representative PM for this cluster (median for robustness)
        cluster_pm = np.median(pm_total[pm_cluster_indices])
        max_separation = pm_distance_factor * cluster_pm  # Arcminutes
        max_separation_rad = (max_separation / 60.0) * (np.pi / 180.0)  # Convert to radians for BallTree
        
        if chatter > 2: 
            print(f"PM Cluster {pm_label}: {pm_cluster_size} sources, Median PM: {cluster_pm:.4f} mas/yr, "
              f"Max separation: {max_separation:.2f} arcmin")
        
        # Get positions and IDs for this PM cluster
        pm_cluster_ra = ra[pm_cluster_indices]
        pm_cluster_dec = dec[pm_cluster_indices]
        pm_cluster_ids = id[pm_cluster_indices]
        
        # Convert positions to SkyCoord for accurate distance calculations
        coords = SkyCoord(ra=pm_cluster_ra*u.deg, dec=pm_cluster_dec*u.deg, frame='icrs')
        
        # Convert positions to radians for BallTree
        ra_rad = coords.ra.to(u.radian).value
        dec_rad = coords.dec.to(u.radian).value
        pos_features = np.vstack((ra_rad, dec_rad)).T
        
        # Build BallTree for efficient neighbor search
        tree = BallTree(pos_features, metric='haversine')
        
       # Find neighbors within the dynamic max_separation
        neighbors = tree.query_radius(pos_features, r=max_separation_rad)
        
        # Build clusters from neighbors (connected components)
        parent = list(range(pm_cluster_size))
        
        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        def union(x, y):
            parent[find(x)] = find(y)
        
        # Union all pairs of neighbors
        for i in range(len(neighbors)):
            for j in neighbors[i]:
                if i != j:
                    union(i, j)
        
        # Group by connected components
        cluster_dict = {}
        for i in range(pm_cluster_size):
            root = find(i)
            if root not in cluster_dict:
                cluster_dict[root] = {'ids': [], 'indices': []}
            cluster_dict[root]['ids'].append(pm_cluster_ids[i].item() if isinstance(pm_cluster_ids[i], np.ndarray) else pm_cluster_ids[i])
            cluster_dict[root]['indices'].append(pm_cluster_indices[i])
        
        # Add clusters with 2 or more sources to the final list
        for root, cluster in cluster_dict.items():
            if len(cluster['ids']) >= 2 :
                clusters.append(cluster)
                if chatter > 2: print(f"  Positional cluster: {len(cluster['ids'])} sources, IDs: {cluster['ids'][:5]}...")

    if chatter > 2:  print(f"Total clustering done in {time.time() - start_time:.2f} seconds")
    if chatter > 2:  print(f"Found {len(clusters)} final clusters")

    # Optionally save clusters to a file
    if output_stuff:
        cluster_data = {
            "n_sources": n_sources,
            "pm_tolerance": pm_tolerance,
            "pm_distance_factor": pm_distance_factor,
            "clusters": [
                {
                    "ids": cluster['ids'],
                    "indices": cluster['indices'],
                    "size": len(cluster['ids']),
                    "median_pm": float(np.median(pm_total[cluster['indices']])),
                    "max_separation": float(pm_distance_factor * np.median(pm_total[cluster['indices']]))
                }
                for cluster in clusters
            ]
        }

    return clusters, cluster_data

 
 

########## main ()  


def mainsub2(chunk,chatter=0):
    # get list of unique sources, output file handle, table, summary
    if onedone:   # limiting processing to just the hipm data (saves lots of time!)
       chunk = hipm_chunk
    tx, outf, sources, ty, outf2, outerr = fileio(rootobs+chunk)
    
    low_pm_nrows = len(ty)
    high_pm_rows = len(tx)
    
    # edit the columns 
    t = make_new(tx) # remove and add columns to the Table
    #print(t[:20],"\n")
    
    # write column HEADERS to ascii/csv file
    col = t.colnames
    nc = len(col)
    for k in range(nc-1): 
        if not onedone: outf.write( f"{col[k]}," )
        outf2.write( f"{col[k]}," )
    if not onedone: outf.write(f"{col[nc-1]}\n")  
    outf2.write(f"{col[nc-1]}\n")  

    #print("evaluating ",col," sources\n")
    if not onedone:
        k9 = 0    # counter records

        # first with the sources with low PM
        for trow in ty:
            print ("call create_csv_  484")
            create_csv_output_record(trow,outf,col,nc)
            k9 += 1
        outf.close()
        print (f"wrote {k9} records in {outf} for the low-pm sources ...")
        
    # reset
    k9 = 0
# now process the high pm rows
    print (f"646 processing {len(tx)} high PM records")
    #
    # clean out high pm records that have bad PM values (PM/error < 10) > remove from tx
    #
    for trow in tx:
        if chatter > 2: print (f"650 clean out PM/error < {Rpm_limit} ")
        Rpm = trow[pm_cds]/np.sqrt( trow['e_pmRA']**2 + trow['e_pmDE']**2)
        if Rpm < Rpm_limit:
            trow['mismatchedPM'] = 2
            create_csv_output_record(trow,outf2,col,nc)
            k9 += 1
    print (f"removed and flag to {k9} records ") 
    print (f"658 check that tx updated too # mismatched is {len(tx[tx['mismatchedPM']==2])}")       
    #        
    # remove these records from further consideration       
    q = tx['mismatchedPM'] != 2 
    print (f"661 now processing remaining {q.sum()} high PM records")
    #
    # TBD clean out high PM records that have too large a mismatch Gmag and V ?
    #
    high_pm_rows = q.sum()
    print (f"processing {high_pm_rows} of high PM records\n")
    print (f"now check for sources with multiple records with same PM (which have likely moved) ")
    #
    # do a clustering analysis to identify those sources that cluster in PM and Ep 2000 coordinates
    #
    ra = tx['ra2000Ep'][q]
    dec = tx['dec2000Ep'][q]
    pm_ra = tx['pmRA_cds'][q]
    pm_dec = tx['pmDE'][q]
    id = tx['IAUNAME'][q]
    clusters, cluster_data = cluster_sources(ra, dec, pm_ra, pm_dec, id, pm_tolerance=0.001, pm_distance_factor=0.01)

    # these records (IAUNAME/ids) are clusters 
    srcids = []
    for cl in cluster_data['clusters']:
        for id1 in cl['ids']:
            srcids.append(id1)
    print (f"of the high PM records, {len(srcids)} have multiple observation IAUNAME IDs\n")
    
    #t = Table(tx[1].data)    
    notinit = []
    for src3 in tx['IAUNAME']:
        if not src3 in srcids:
           notinit.append(src3)
    print (f"and {len(notinit)} are single sources")       

    # write the hi PM records that do not cluster
    for name in notinit:
        trow = tx[np.where(tx['IAUNAME'] == name)]
        create_csv_output_record(trow[0],outf2,col,nc)
        k9 += 1
    print (f"wrote so far {k9} records in {outf2} for the single sources ... ")
    #       
    # process the clusters    
    #
    #tx.remove_index(pm_cds)
    tx.add_index(iauname)

    for cluster in clusters:   
        
        # find the cluster members
        indices = cluster['indices']
        names   = cluster['ids']
        tab_src = tx.loc[names]        
        k9+=1                       
        tab_ma = ma.asarray(tab_src)
        nrow = tab_ma.size
                
        if type(tab_src) != astropy.table.row.Row: # is Table object ?
            new_row = tab_src[0] # on new row per iauname  
        else: 
            new_row = tab_src  # is Table Row
            
        # how to deal with now merging records with multiple IAUNAME: take the first one    
        firstname = new_row['IAUNAME']    
            
        # now check if there are double entries and remove those based on obs_epoch and iauname (src)
        if nrow > 1:
            tab_src = Table(tab_src)
            idup = [0]
                            
        if nrow == 1:  # meaning no matches within 20' distance 
                raise RuntimeError ("call create_csv_  539 nrow==1  should not occur!!")
                
        else:    
            base = evaluate_extended_nature(tab_src)
            obsidlist = tab_src['OBSIDS'][0]
            epochlist = str(tab_src['EPOCHS'][0])
            srcnumlist = str(tab_src['SRCNUMS'][0])
            for ix in np.arange(1,nrow):
                obsidlist += "_"+tab_src['OBSIDS'][ix]   
                epochlist += "_"+str(tab_src['EPOCHS'][ix])
                srcnumlist+= "_"+str(tab_src['SRCNUMS'][ix])
            new_row['OBSIDS'] = obsidlist
            new_row['EPOCHS'] = epochlist
            new_row['SRCNUMS']= srcnumlist                
        #
         #   print (f"=== === === === === \n 746 tab.ma.size={nrow}\n new_row is\n{new_row};\ntype is {type(new_row)}\n === === === === === 491\n")
         #   print (f"OBSIDS, EPOCHS:\n{new_row['OBSIDS']}\n{new_row['EPOCHS']}")    
              
            check = True    
            for band in bands:   # loop over all filters - merge into one 
                nobs_all = 0
                qf = tab_ma[band+"_QUALITY_FLAG_ST"]
                if chatter > 2: print (f"751 checking {band}: quality_flag_st qf = {qf}\n qf.data={qf.data}523\n")
                qfaa = qf.data
                qfok = []
                for x6 in qfaa:
                    if len(x6) > 6: 
                        qfok.append(x6)
                for x7 in tab_ma[band+"_NOBS"].data:
                    nobs_all += x7
                    if chatter > 3: print (f"759  for {x7} nobs_all = {nobs_all} cc")
                if chatter > 2: print (f"757 qfok = {qfok}")        
                if nrow > 1:
                    # QUALITY   
                    if (catalog == 'SUSS') | (catalog == "testsuss"):
                        qf = qfok
                        if len(qf) > 1:
                            if chatter > 2: print (f"759 qf={qf}  err=\n{tab_src[band+'_AB_MAG_ERR']}")
                            if check:
                                qual_st_out,ranked_quals = qual_st_merge(qf,err=tab_src[band+"_AB_MAG_ERR"])
                                qual_out = qual_stflag2decimal(qual_st_out)
                        elif len(qf) == 1:
                                qual_st_out, ranked_quals, qual_out = qf[0],qf[0],qual_stflag2decimal(qf[0])
                        else:
                                qual_st_out,ranked_quals, qual_out = "", "", -2147483648
                 
                    elif (catalog == 'UVOTSSC2') | (catalog == "test"):
                        qual_st_out,ranked_quals = qual_st_merge(qf,err=tab_src[band+"_ABMAG_ERR"])
                        qual_out = qual_stflag2decimal(qual_st_out)
                        ix = qual_st_out == -999
                        #print (f"x: {type(ix)}  {ix}")
                        if ix:
                            if chatter > 2: print (f" qual {qual_out}, st:{qual_st_out}\n")
                            qual_st_out = "" 
                            qual_out = np.nan
                    #if len(qual_st_out) != 12:
                    #    raise RuntimeError(f"448 ERROR: {k9} the quality flag is wrong {qual_st_out}")
                    new_row[band+'_QUALITY_FLAG_ST'] = qual_st_out
                    new_row[band+'_QUALITY_FLAG'] = qual_out
                    # SIGNIFICANCE  - pick best value - not sure if in SUSS
                    if check & (catalog == 'UVOTSSC2') | (catalog == 'test'):
                        signif = tab_src[band+"_SIGNIF"]
                        signif = signif.data
                        qsignif = signif[np.isfinite(signif)]
                        if len(qsignif) == 1: 
                            new_row[band+"_SIGNIF"]= qsignif
                        elif len(qsignif) > 1:     
                            new_row[band+"_SIGNIF"]= np.max(qsignif)
                        else: 
                            new_row[band+"_SIGNIF"]= signif[0]  
                            
          # the other case, one row input, is taken case of by the version 1 
          #  med_mag, max_mag, sigma_mag, varFlag_mag, mean_mag = \
          #     combine_a_set_of_(tab_src[band+"_AB_MAG"],
          #        tab_src[band+"_AB_MAG_ERR"],N=5,
          #        qual=new_row[band+"_QUALITY_FLAG"])
            if check & (catalog == 'SUSS') | (catalog == "testsuss"):
                magx = tab_src[band+"_AB_MAG"]
                errx = tab_src[band+"_AB_MAG_ERR"]
                #print (f"magx={magx}, errx={errx}")
                #
                #  note: if we are merging data that have already multiple observations 
                #        included, we should consider expanding those to do stats
                #        WARNING that is not done here (although we keep track on the
                #        nobs_all 
                #
                nObs, med_mag, chisq, sigma_mag, skew, var3 = stats(magx,err=errx,syserr=0.005)    
                qmag = magx > 5
                if len(magx[qmag]) > 0:
                    min_mag = np.min(magx[qmag])
                else:
                    magx = np.nan
                    min_mag = np.nan    
                new_row[band+'_AB_MAG'] = med_mag
                new_row[band+'_AB_MAG_ERR'] = sigma_mag
                new_row[band+'_AB_MAG_MIN'] = min_mag
                new_row[band+'_NOBS'] = nobs_all # nObs - use to test variability
                if nObs > 1:
                    new_row[band+'_CHI2RED'] = chisq/nObs
                else:
                    new_row[band+'_CHI2RED'] = np.nan                  
                new_row[band+'_CHISQ'] = chisq
                new_row[band+'_SKEW'] = skew
                new_row['mismatchedPM'] = 3

            elif (catalog == 'UVOTSSC2') | (catalog == "test"):    
                magx = tab_src[band+"_ABMAG"]
                errx = tab_src[band+"_ABMAG_ERR"]
                nObs, med_mag, chisq, sigma_mag, skew, var3 = stats(magx,err=errx,syserr=0.005)    
                min_mag = np.min(magx)
                new_row[band+'_ABMAG'] = med_mag
                new_row[band+'_ABMAG_ERR'] = sigma_mag
                new_row[band+'_ABMAG_MIN'] = min_mag
                new_row[band+'_VAR3'] = var3   # variability on 3-sigma 
                new_row[band+'_NOBS'] = nObs
                #new_row[band+'_CHI2RED'] = chisq/nObs
                new_row[band+'_CHISQ'] = chisq
                new_row[band+'_SKEW'] = skew

            if check & (catalog == 'SUSS') | (catalog == "testsuss"):
                new_row[band+'_EXTENDED_FLAG'] = base[band]  
            elif (catalog == 'UVOTSSC2') | (catalog == "test"):
                new_row[band+'_EXTENDED'] = base[band]  

        if chatter > 0: print ("call create_csv_  834 ")
        create_csv_output_record(new_row,outf2,col,nc)  

    outf2.close()
    outerr.close()

    #
    #  join the two files
    #
    command = f"cat {outf} {outf2} > {instub+outstub+'.fits'}"
    try:
        if onedone == False:
            os.system(command)
    except:
        pass

######### ######## ####### END



    