#   python3
#  take a catalogue like SUSS5 where there are multiple observations of a single source
#
#  scan the catalogue and create initially a lookup table for the records for each source
#  combine the photometry in the multiple observations as follows:
#  For each colour/flux
#    - compute median, sigma and range  
#    - compute error on median value
#    - populate a variability column as range / median
#  For other data columns
#    - Find a consistent distance/parallax/PM
#    - determine a consistent source extended flag/value
#    - compare classifiers from other catallogues.
#
#  Augmented: which columns & criteria to include in definition of class ?
#
#  WARNING : output format of maked float columns converts to string and there can be 
#    some entries with '--' (need to find out why)
#    There is a python routine in cats/xmm2athena.py which converts to a floating 
#    point array called: _fix_coordinates(column)
#
# original npmk, October-November 2022
# 
# updates for v2:
# calculate the median instead of the middle of the magnitudes, ie, sort the values, 
#   take the the one at index halfway
# calculate the chi-squared for variability 
#    for version 2.1 also calculate chisq-reduced for quality=0 data
# report the brightest magnitude
# report the OBSIDs, epochs, and SRCNUM as a merged list
#    remaining :
#  ? report upper limits if observed same source in other filter (from summary)  
# merged extended & variability flags needed ? 
#
# to do: fail calculating mag_min, due to failed masks, so -999 is always lowest
#
# 2025-02-15: added code to select only qual=0 points to calculate chisq, but keep the 
#   value using all data in chisq_all 
#    @npmkuin, 2023, Copyright 2023-2025, license 3-clause BSD style (see github.com/PaulKuin)

__version__ = 2.1
########################################################################################
# globals
import os
import numpy as np
import numpy.ma as ma
import astropy, numpy
from astropy.io import fits, ascii as ioascii
from astropy.table import Table, Column, MaskedColumn

# updates: return list of obsids, epochs [can be used to search for original records)
#  return median, skew, etc.
# not yet: ABmag upper limits from SUMMARY 

# define globals:

fix_duplicate = False # do not set! needs a lot of eval to get it to work - 
# the start/stop times are for agregated images for one obsid, with often 
# multiple image snapshots being merged 
colsofinterest = ["IAUNAME","obs_epoch","OBSIDS","EPOCHS",'UVW2_ABMAG','UVM2_ABMAG','UVW1_ABMAG','U_ABMAG',"RA_ERR","DEC_ERR","UVW2_NOBS","UVM2_NOBS","UVW1_NOBS","U_NOBS"]

iauname="IAUNAME" 
bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent
#catalog = 'UVOTSSC2'
catalog = "SUSS" #'testsuss' #'SUSS'
outtable = None

if catalog == 'SUSS':
    rootobs = '/Volumes/DATA11/data/catalogs/test_2/'
    chunk = "sussxgaiadr3_ep2000.fits"   
    sussdir= '/Volumes/DATA11/data/catalogs/suss6./'
    suss='XMM-OM-SUSS6.2.fits'
    summary_file = "suss_summary.fits"
elif catalog == 'UVOTSSC2':
    rootobs = '/Users/data/catalogs/uvotssc2/'
    chunk =  'uvotssc2_in_photometry.fits'  
    #'uvotssc2_in_photometry.fits'  'test_main_in500k.fits' 
    # combine the files in  one with fix_uvotssc2.combine_input_files
elif catalog == 'test':
    rootobs = '/Users/data/catalogs/uvotssc2/'
    chunk = 'test_sources.fits'
elif catalog == 'testsuss':
    rootobs = '/Volumes/DATA11/data/catalogs/test/'
    chunk = "sussxgaiadr3_ep2000.fits"   
    summary_file = "test_cat_summary.fits"

#output file name is by appending "_singlerecs" to input file name (see fileio() ).
tmax = 500000000 # maximum number of records to read in for test

################ END GLOBAL SECTION ######################################################


def make_chunks(chunk, dn=3020000):
    """
      the uvotssc input file is >54GB and too large to have in memory, so we need to process in 
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
         outf = f"schunk_{k:02}.fit"
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
    
   
def make_new(tx):
    # remove and add columns to main table
    tx.remove_columns(['N_SUMMARY','N_OBSID'])
    for  b in bands:
        if (catalog == 'SUSS') | (catalog == "testsuss"): 
            tx.remove_columns([ b+'_RATE',b+'_RATE_ERR',b+'_SKY_IMAGE'])
            tx.remove_columns([ b+'_VEGA_MAG',b+'_VEGA_MAG_ERR'])
            tx.remove_columns([ b+"_MAJOR_AXIS",b+"_MINOR_AXIS",b+"_POSANG"])
            tx.add_columns([0.,0.,0.,None,0.],
                names=[b+"_CHI2RED",b+"_CHISQ",b+"_NOBS",b+'_AB_MAG_MIN',b+'_SKEW'])  
        elif (catalog == 'UVOTSSC2') | (catalog == "test"):        
            tx.add_columns([0.,0.,0.,None,0.],
               names=[b+"_CHI2RED",b+"_CHISQ",b+"_NOBS",b+'_ABMAG_MIN',b+'_SKEW'])  
#               names=[b+"_CHISQ",b+"_NOBS",b+'_ABMAG_MIN',b+'_VAR3',b+'_SKEW'])  
            tx.remove_columns([b+"_FLUX",b+"_FLUX_ERR"])
    tx.add_columns([None,None,None],names=['OBSIDS', 'EPOCHS','SRCNUMS'])
           
    return tx         

def stats(array, err=[None], syserr=0.005):
    """
    input parameters
    ----------------
       array: 1D array
          require all a > 0
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
    
    """
    import numpy as np
    from scipy.stats import kurtosis
    from scipy.stats import skew
    
    # first remove nan values and if no elements left, return, otherwise continue
    #n1 = len(array) input can be a float
    a = np.asarray(array)
    q = (np.isnan(a)) | (a <= 0.)  # added a <= 0 2023-06-16
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
       var = var*na/(na-1)
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
    ncChisq = (a-p)*(a-p)/var
    ncChisq = ncChisq.sum() 
    # variability3 (x_i - x_k)/(3err_k) with k at median error sd 
    V3 = (a - median)/(3*sd)
    V3 = np.max(V3)
    #kur = kurtosis(a,nan_policy='omit')
    sk = skew(a,nan_policy='omit')
    
    return n, median, ncChisq, sd, sk, V3   
      

def fileio(infile,outstub="_singlerecs",outdir=""):
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
    
    if (ft.lower == 'csv'):
       t = ioascii.read(infile)  # table object
    elif (ft[:3] == 'fit'):
       f = fits.open(infile)
       x = f[1].data #[:tmax]  # LIMITED TABLE FOR TESTS
       if len(f) == 3:
           summ = f[2].data # summary - has obstimes and exposures of obsids
       else:
           with fits.open(summary_file) as sumfile:
               sumx = sumfile[1].data
       t = Table(x)   
       summ = Table(sumx)
       x = ''
       f.close()
    if len(t) < 1:
        raise IOError("the catalogue table loading failed.\n")   
        
    if catalog == 'UVOTSSC2': #add the same kind of column as in SUSS for the processing
        # add the logical quality columns
        bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent
        for band in bands:
            decimalflag = t[band+"_QUALITY_FLAG"]
            out = []
            for flag in decimalflag:
               out.append(qual_dec2qual_stflag(flag))
            t.add_column(out,name=band+"_QUALITY_FLAG_ST")
            
    outfile=instub+outstub+".csv"
    outf = open(outdir+outfile,'w')  # file handle
    # 
    # we use the astropy to flag data that have the same IAUNAME
    # after that we can call up the record/row selection by that.
    #
    t.add_index(iauname)
    sources = t[iauname]
    sources = np.unique(sources)
    print (f"NUMBER OF UNIQUE SOURCES = {len(sources)}")
    summ.add_index("OBSID") # to get a list of exposure times 
    # to find the table with rows that corresponds to a source, 
    # tab_src = t.loc[sources[k]] 
    #   the type is eather a astropy.table.row.Row for a single record, or
    #   astropy.table.table.Table if multiple records
    try:
        del f[1].data
        del f[2].data
    except: pass
    return t, outf, sources, summ

   

def evaluate_extended_nature(tsrc,summ):
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
    for k in range(nc-1): 
       #print     ( f"{trow[col[k]]}," )
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
          #print (f"qual_st_merge409 {60*'-'}\noriginal {qual_st}\n\ninteger flag values are: {aa}")    
          #print (f"qual_st_merge410 summed rows: {aa}")
          qa = ma.min(aa) == aa
          k = ma.where(qa)[0]
          qualout = qual_st[k][0]
          #   print (f"387 {qualout}, len -> {len(qualout)}")
          # qualout sets the same flags as in SUSS5.0
          # force extended and bright star 
          
          #print (f"390 qualout is {qualout}, len -> {len(qualout)}, {type(qualout)}\n")        
          #for k in range(len(qual_st)):  # all rows
          #  if ~qual_st.mask[k]:
          #    rerank = False
          #    if qual_st[k][6] == 'T':
          #        x = list(qualout)
          #        x[6] = 'T'
          #        qua = ""
          #        qualout=qua.join(x)
          #        rerank = True
          #    if (qual_st[k][8] == 'T') and (err[k] < 0.03):   
          #        x = list(qualout)
          #        x[8] = 'T'
          #        qua = ""
          #        qualout=qua.join(x)
          #        rerank = True
          #if rerank:
          #    aa =  ranked_qual_st2qual(ma.asarray([qualout]))
          #    ranked_quals = aa[0]   
          #print (f"qualout to return is {qualout}, len -> {len(qualout)}\n")     
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
 



########## main ()  


def mainsub(chunk):
    # get list of unique sources [iauname], output file handle, table, summary
    t1, outf, sources, summ = fileio(rootobs+chunk)
    # edit the columns 
    t = make_new(t1) # remove and add columns to the Table
    #print(t[:20],"\n")
    
    # write column HEADERS
    col = t.colnames
    nc = len(col)
    for k in range(nc-1): 
        outf.write( f"{col[k]}," )
    outf.write(f"{col[nc-1]}\n")  
    outfits = rootobs+chunk[:-4]+"_srcnum.fits"     

    #print("evaluating ",col," sources\n")
    k9 = 0    # counter records

    for src in sources[:]:   
        if int(k9/5000)*5000 == k9:
            print (f"451 src number={k9}, source={src} "); 
        k9+=1
        tab_src = t.loc[src]
        tab_ma = ma.asarray(tab_src)
        nrow = tab_ma.size
            
        # now check if there are doubles and remove those based on obs_epoch and iauname (src)
        if nrow > 1:
            tab_src = Table(tab_src)
            tab_src.sort('obs_epoch')
            #print (f"473 nrow={nrow}, {tab_ma.size},  tab_src={tab_src} ;--> {type(tab_src)} ; tab_ma={tab_ma}")   
            oep = tab_src[0]['obs_epoch']
            idup = [0]
            for itn in np.arange(1,nrow):
                if tab_src[itn]['obs_epoch'] != oep:
                    idup.append(itn)
                    oep = tab_src[itn]['obs_epoch']
            if len(idup) < nrow:
                #print (f"*********\n duplicate records found nrow={nrow}, rows={idup}:\n")
                #print(f"\n ********* original: {tab_src[colsofinterest]} \n ******* corrected: {tab_src[idup][colsofinterest]}\n")
                if fix_duplicate:  # test found match but different observed filters
                   tab_src = tab_src[idup]
                   nrow = len(idup)
                   tab_ma = ma.asarray(tab_src)
                   # test found new_row['obs_epoch'] not to be adjusted? Output corruption
                
        if type(tab_src) != astropy.table.row.Row: # is Table object ?
            new_row = tab_src[0] # on new row per iauname  
        else: 
            new_row = tab_src  # is Table Row
            
        base = evaluate_extended_nature(tab_src,summ)
        if nrow == 1:
            new_row['OBSIDS'] = f"{tab_src['OBSID']}"  #tab_src['OBSID']
            new_row['EPOCHS'] = f"{tab_src['obs_epoch']:.5f}" #tab_src['obs_epoch']
            new_row['SRCNUMS'] = f"{tab_src['SRCNUM']}" #
            
        else:
            obsidlist = f"{tab_src['OBSID'][0]}"
            epochlist = f"{tab_src['obs_epoch'][0]:.5f}" #str(tab_src['obs_epoch'][0])
            srcnumlist = f"{tab_src['SRCNUM'][0]}" 
            for ix in np.arange(1,nrow):
                obsidlist += f"_{tab_src['OBSID'][ix]}"  
                epochlist += f"_{tab_src['obs_epoch'][ix]:.5f}"
                srcnumlist += f"_{tab_src['SRCNUM'][ix]}"
            new_row['OBSIDS'] = obsidlist
            new_row['EPOCHS'] = epochlist
            new_row['SRCNUMS'] = srcnumlist
                                
        for band in bands:
            qf = tab_ma[band+"_QUALITY_FLAG_ST"]
            qf2 = tab_ma[band+"_QUALITY_FLAG"]
            qfaa = qf.data
            qfa2 = qf2.data
            qf_q = qfaa != ''
            qf2_q= qfa2 != ''
            qual_out = 0
            #print (f"546 qfaa -- {type(qfaa)}  {qfaa} selected: {qf[qf_q]} {len(qf[qf_q])} ")
             
            if nrow > 1:
                # QUALITY   
                if (catalog == 'SUSS') | (catalog == "testsuss"):
                    qf = qf[qf_q]
                    if len(qf) > 1:
                        qual_st_out,ranked_quals = qual_st_merge(qf,err=tab_src[band+"_AB_MAG_ERR"])
                        qual_out = qual_stflag2decimal(qual_st_out)
                    elif len(qf) == 1:
                        #print (f"qual_stflag2decimal input =|{qf}|")
                        qual_st_out, ranked_quals, qual_out = qf[0], qf[0], qual_stflag2decimal(qf[0])   
                    else:
                        qual_st_out,ranked_quals, qual_out = "", "", -2147483648 # need to define band+quality_flag as a float     
                elif (catalog == 'UVOTSSC2') | (catalog == 'test'):
                    qual_st_out,ranked_quals = qual_st_merge(qf,err=tab_src[band+"_ABMAG_ERR"])
                    qual_out = qual_stflag2decimal(qual_st_out)
                    ix = qual_st_out == -999
                    #print (f"x: {type(ix)}  {ix}")
                    if ix:
                       print (f" qual {qual_out}, st:{qual_st_out}\n")
                       qual_st_out = "" 
                       qual_out = np.nan
                #print (f"566  qual_st_out={qual_st_out},  qual_out={qual_out} *****")       
                #if len(qual_st_out) != 12:
                #    raise RuntimeError(f"448 ERROR: {k9} the quality flag is wrong {qual_st_out}")
                new_row[band+'_QUALITY_FLAG_ST'] = qual_st_out
                new_row[band+'_QUALITY_FLAG'] = qual_out
                # SIGNIFICANCE  - pick best value - not sure if in SUSS
                if catalog == 'UVOTSSC2':
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
            if (catalog == 'SUSS') | (catalog == "testsuss"):
                magx = tab_src[band+"_AB_MAG"]
                errx = tab_src[band+"_AB_MAG_ERR"]
                #print (f"magx={magx},\nerrx={errx}\n = = = = = next ")
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
            elif (catalog == 'UVOTSSC2') | (catalog == "test"):    
                magx = tab_src[band+"_ABMAG"]
                errx = tab_src[band+"_ABMAG_ERR"]
                nObs, med_mag, chisq, sigma_mag, skew, var3 = stats(magx,err=errx,syserr=0.005)    
                qmag = magx > 5 
                if len(magx[qmag]) > 0: 
                    min_mag = np.min(magx[qmag])
                else: min_mag = np.nan
                new_row[band+'_ABMAG'] = med_mag
                new_row[band+'_ABMAG_ERR'] = sigma_mag
                new_row[band+'_ABMAG_MIN'] = min_mag
            #new_row[band+'_VAR3'] = var3   # variability on 3-sigma 
            if nrow >1:
                if (qual_out == 0) and (len(qf) > 1):  # there is a good detection amongst list
                    qchi = qf2[qf2_q] == 0  # index data with qual==0
                    nq0 = len(qchi)
                    if np.isnan(magx).all():
                       chisq = np.nan
                       skew  = np.nan
                    else:
                       nObs, med_mag, chisq, sigma_mag, skew, var3 = stats(magx[qchi],err=errx[qchi],syserr=0.005)   
                    if nq0>1: 
                        new_row[band+'_CHI2RED'] = chisq /(nq0-1)  # this should only be based on the qual==0 data points
                   
                new_row[band+'_CHISQ'] = chisq   # this is based on any quality 
                new_row[band+'_SKEW'] = skew
                new_row[band+'_NOBS'] = nObs
            if (catalog == 'SUSS') | (catalog == "testsuss"):
                new_row[band+'_EXTENDED_FLAG'] = base[band]  
            elif (catalog == 'UVOTSSC2') | (catalog == "test"):
                new_row[band+'_EXTENDED'] = base[band]  
            # edit masks:    
            
        # perhaps change mask value here for some cols
        create_csv_output_record(new_row,outf,col,nc)  
    outf.close()




######### ######## ####### END END END END 
