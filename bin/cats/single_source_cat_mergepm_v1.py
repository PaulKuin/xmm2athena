#   python3
#  refine the SUSS5.0 single source catalogue which has multiple observations of a 
#  single source of high proper motion
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
# npmk, October-November 2022
#
#  Note: this version computes the VerFlag number for variability, which in version 2
#        has been replaced with athe chiSquared measure.
#

# globals
import numpy as np
import numpy.ma as ma
import numpy, astropy
from astropy.io import fits, ascii as ioascii
from astropy.table import Table, Column, MaskedColumn

"""
copy of single source cat .py with edits to merge rows with same PM

"""
bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent
   

def make_new(tx):
    # remove and add columns to main table
    return tx         


def fileio(infile,outstub="_fixpm",outdir=""):
    #
    # use the input file name infile to construct the output filename
    # open the output file
    # read in the table from the input file
    # return the table and the output file handle
    import numpy as np
    a = infile.rsplit('.')
    instub = ""
    ft = a[-1]
    for x in a[:-1]: 
        instub += x+"."
    instub = instub[:-1]    
    print (f"reading {infile} ...\n")
    if (ft.lower() == 'csv'):
       t = ioascii.read(infile)  # table object
    if len(t) < 1:
        raise IOError("the catalogue table loading failed.\n")   
    outfile=instub+outstub+".csv"
    outf = open(outdir+outfile,'w')  # file handle
    print (f"output is written to {outdir+outfile} ")
    # header first:
    col = t.colnames
    nc = len(col)
    for k in range(nc-1): 
       outf.write( f"{col[k]}," )
    outf.write(f"{col[nc-1]}\n")       
    # 
    # we use the astropy to flag data that have the same PM
    # first split the catalog into two parts and use here the one with PM
    #
    # print (f"writing {} records without PM")
    # after that we can call up the record/row selection by PM.
    # select records with PM value 
    print (f"now processing {len(t)} records with PM to reduce duplication")    
    t.add_index("gaiadr3_pm")
    sources = t["gaiadr3_pm"]
    import numpy
    sources = numpy.unique(sources)
    print (f"there are {len(sources)} unique sources with PM to be processed.")
    #   the type is eather a astropy.table.row.Row for a single record, or
    #   astropy.table.table.Table if multiple records
    return t, outf, sources

   
def combine_a_set_of_(mag,err,qual=None,N=3):
    #
    # input is an array for one filter where some values are missing
    # determine mean, min, max, and 1-sigma
    # if mean-min or max-mean is larger than N * RMS(sigma,err): set variable flag
    #
    # return mean_mag, quality, max_mag, sigma_mag, varFlag_mag
    mag = ma.asarray(mag)
    varFlag = ""
    #print (f"mag line 215 {mag}, type is {type(mag)}\n")
    # now check if this data has no qual set - then this filter was not observed
    # use_row = qual > 0
       
    q = np.isfinite(mag) 
    N = len(mag[q])   
    if N < 2:  # only one good value or None
       if N == 0:   # # again, filter not observed or only once - no need to average, etc.
          return None, None, None, varFlag, None
    
       med_mag = mean_mag = min_mag = max_mag = mag[q]
       sigma_mag = err[q]
       return med_mag, max_mag, sigma_mag, varFlag, mean_mag
       
    # N > 1
    mean_mag = mag[q].mean()
    sigma_mag = mag[q].std()
    er=np.min(err[q])
    #error2 = er*er+sigma_mag*sigma_mag # errors squared on N array elements  
    min_mag = mag[q].min()
    max_mag = mag[q].max()
    med_mag = 0.5*(max_mag+min_mag)

    # variability: 
    #   (a) find data with lower half error region "best data"
    #   (b) compute mean mag and standard deviation for that
    #   (c) set varFlag if there are data outside N*sigma from mean if error on data is 
    #       small enough , i.e., 3-sigma from mean+3 sigma
    es = ma.asarray( np.sort(err[q]) )
    nk = np.int(len(es)/2)
    q2 = es <= es[nk] # (a)
    if ma.sum(q2) <= 1:   # only one left to test 
        varFlag = ""
    else:    
        mean2 = mag[q][q2].mean()  # (b)
        sigma2 = mag[q][q2].std()  # (b)
        a = ma.abs(mean2 - mag[q][q2]) > N * (ma.sqrt(sigma2) + err[q][q2])
        varFlag = np.int(ma.any(len(a) > 0))
    
    return med_mag, max_mag, sigma_mag, varFlag, mean_mag


def evaluate_extended_nature(tsrc):
    # compare the extended flags for long exposures, and consistency in the 
    # various bands 
    # input for each observation in target
    #
    # is it extended in all observations?
    #  with the same PA?
    # in all bands ?
    #  with sufficient exposure
    #
    # => was the observation trailing? 
    # => transient?
    # => consistent?  
    #
    # simplest solution
    #tsumm = summ.loc[obsids] # find records for obsids 
        # refer to exposure by tsumm['EXPOSURE_UVW2'] etc. 
    # need check all
    bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent
    base = {'UVW2':255,'UVM2':255,'UVW1':255,'U':255,'B':255,'V':255}
    for band in bands:
       exflag = tsrc[band+'_EXTENDED_FLAG']
       q = exflag != 255
       if np.sum(q) > 0:
           base[band]= np.int(np.mean(exflag[q]) > 0.5) # more than half say extended
    # so valid bands have base not equal to 255 
        # TBD: ...
    # updates to the base solution:
    # case 1: single obs. in one band: pass the info
    # case 2: multiple obs, but after weeding only 1 : pass merged info
    # case 3: all bands extended in long exposures : OK extended
    # case 4: all bands not extended in long exposures: OK point source
    # case X: something else ...
     
    return base
    

def create_output_record(trow,outf):
    # create the output record to write to a csv file 
    # input is a table object for one source after the magnitudes have been 
    # combined, and the extended nature has been evaluated
    #
    #nc = len(col)
    #print (f"line 286 trow colnames:\n{col}\n")
    for k in range(nc-1): 
       #print     ( f"{trow[col[k]]}," )
       outf.write( f"{trow[col[k]]}," )
    outf.write( f"{trow[col[nc-1]]}\n" )

    
# The set flags are ranked as follows (rank:flag bit) (bad->worse->verybad order):
f_rank = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 6:2, 7:0, 8:3, 9:10,10:9, 11:11}
inverse_f_rank ={0:7,1:5,2:6,3:8,4:2,5:1,6:3,7:4,8:0,9:10,10:9,11:11}
f_rk = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 7:0, 9:10,10:9, 11:11} # without f_carry

def qual_st_merge(qual_st,err=0.02):
    # produce qoalout, from qual_ft (derive the ranked flag)    
    # output ranked_quals is the numerical quality flag, ranked using f_rk defined below
    #
    # The set flags are ranked as follows (rank:flag bit) (bad->worse->verybad order):
    #f_rank = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 6:2, 7:0, 8:3, 9:10,10:9, 11:11}
    #inverse_f_rank ={0:7,1:5,2:6,3:8,4:2,5:1,6:3,7:4,8:0,9:10,10:9,11:11}
    #f_rk = {0:8, 1:5, 2:4, 3:6, 4:7, 5:1, 7:0, 9:10,10:9, 11:11} # without f_carry
    # if one of the following flags is set in any of the source detections, carry it
    f_carry = [6,8] # original flags, i.e., [1,3] in ordered flags
    qualout = ""
       
    if len(qual_st) == 1:  # nothing to decide
            return qual_st,  ranked_qual_st2qual(qual_st)
    else: 
       nq = 12  # number of bits
       qual = "FFFFFFFFFFFF"     # start with first one
       # check each flag is not set, select the lowest set
       if isinstance(qual_st,np.flexible):
           print (f"\nWARNING 208 (char) np.flexible : is {qual_st} just a string? \n\n")
           return qual_st,  ranked_qual_st2qual(qual_st)
       else:
          print (f"prior to ranking {qual_st}")
          ranked_quals = aa =  ranked_qual_st2qual(qual_st) #np.array(aa)
          #print (f"{60*'-'}\noriginal {qual_st}\n\ninteger flag values are: {aa}")    
          #print (f"summed rows: {aa}")
          qa = ma.min(aa) == aa
          k = ma.where(qa)[0]
          #print (f"218 qa={qa}, k={k} {qual_st[~qual_st.mask]} and qual_st[k] = {qual_st[k]}")
          qualout = qual_st[k][0]
          #print (f"225 qualout = {qualout}")   
 
          # qualout sets the same flags as in SUSS5.0
          # force extended and bright star 
          
          #for k in range(len(qual_st)):  # all rows
          #    if ~qual_st.mask[k]:
          #       if qual_st[k][6] == 'T':
          #           x = list(qualout)
          #           x[6] = 'T'
          #           qualout = ""
          #           for i in x: qualout += i
          #       if (qual_st[k][8] == 'T') and (err[k] < 0.03):   
          #           x = list(qualout)
          #           x[8] = 'T'
          #           qualout = ""
          #           for i in x: qualout += i
                  
          #print (f"qualout to return is {qualout}, len -> {len(qualout)}\n")        
          return qualout , ranked_quals       
          
def qual_stflag2decimal(flag_st):
    if flag_st == "            ": return -2147483648
    out = 0
    for k in range(12):
        if flag_st[k] == "T":
           out += 2**k
    return out       

def ranked_qual_st2qual(qual_st):
    # convert string qual to number qual after ranking with f_rk 
    aa = []
    for k in range(len(qual_st)):  # all rows
        x = ma.asarray(qual_st[k])
        if ~x.mask:
           bb = 0
           for j in f_rk.keys():
                s = str(x)
                if s[f_rk[j]] == 'T' :
                    bb += 2**(f_rk[j])
           aa.append(bb)
        else:
           aa.append(-9999)      
    return ma.masked_values(aa,-9999) 
  
########## main ()  ##########
  
    
sussdir= '/Users/data/catalogs/suss_gaia_epic/'
suss='XMMOM_SUSS5.0_Sources_GaiaDR3motion.part1.csv'
t, outf, sources = fileio(sussdir+suss)
#t = make_new( t1) # replace make_record inside loop
#print(t,"\n")
col = t.colnames
nc = len(col)
#print("evaluating ",col," sources\n")
k9 = 0    

for src in sources[:]:   
    k9 += 1
    mergethem = True 
    tab_src = t.loc[src]
    tab_ma = ma.asarray(tab_src)
    if tab_ma.size == 1:
        new_row = tab_ma
        n_row = 1
    elif tab_ma.size > 1:
        n_row = tab_ma.size
        new_row = tab_ma[0] 
    # screen low PM out
        ra = tab_ma['RA']
        dec = tab_ma['DEC']
        pm = tab_ma['gaiadr3_pm']*1e-3*21.5
        for k in range(1,len(ra)):
            if ma.sqrt( (ra[k]-ra[0])**2 + (dec[k]-dec[0])**2 ) > pm.mean(): 
                mergethem == False
    else: 
        raise RuntimeError(f"empy data ? {k9}: {tab_ma}")#f"ERROR: the type of the selection is not recognised.")    

    if n_row > 1 and mergethem:
      base = evaluate_extended_nature(tab_src)
      for band in bands:
    
         qf = tab_ma[band+"_QUALITY_FLAG_ST"]
         print (f"304 (main) {band} k9={k9} qf={qf}")

         if n_row > 1:  # must be multiple elements 
            if np.all(qf.mask): 
               qual_st_out, ranked_quals = qfx[0],""
            else:    
               qual_st_out,ranked_quals = qual_st_merge(qf,err=tab_src[band+"_AB_MAG_ERR"])
               new_row[band+'_QUALITY_FLAG_ST'] = qual_st_out
               new_row[band+'_QUALITY_FLAG'] = qual_stflag2decimal(qual_st_out)
        # the other case, one row input, is taken case of by the earlier copy 

         med_mag, max_mag, sigma_mag, varFlag_mag, mean_mag = \
               combine_a_set_of_(tab_src[band+"_AB_MAG"],
               tab_src[band+"_AB_MAG_ERR"],N=5,
               qual=new_row[band+"_QUALITY_FLAG"])
         new_row[band+'_AB_MAG'] = med_mag
         new_row[band+'_AB_MAG_ERR'] = sigma_mag
         new_row[band+'_AB_MAG_MAX'] = max_mag
         new_row[band+'_VAR'] = varFlag_mag
         new_row[band+'_EXTENDED_FLAG'] = base[band]  
                
         # same for flux
         med_flx, max_flx, sigma_flx, varFlag_flx, mean_flx = \
           combine_a_set_of_(tab_src[band+"_AB_FLUX"],tab_src[band+"_AB_FLUX_ERR"],
             N=5,qual=new_row[band+"_QUALITY_FLAG"])
         new_row[band+'_AB_FLUX'] = med_flx
         new_row[band+'_AB_FLUX_ERR'] = sigma_flx
         new_row[band+'_AB_FLUX_MAX'] = max_flx
        
         # the significances are not updated, but could insert the minimum 
         # of all valid entries
    if mergethem:    
        create_output_record(new_row,outf) 
    else:
        for row in new_row:
            create_output_record(row,outf) 
    
outf.close()
print (f"output file completed {outf}")

    