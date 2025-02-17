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
# npmk, October-November 2022


# globals
import numpy as np
import numpy.ma as ma
import astropy, numpy
from astropy.io import fits, ascii as ioascii
from astropy.table import Table, Column, MaskedColumn

"""
column mapping 

[     'IAUNAME',
 'N_SUMMARY',
 'OBSID',      => obsid list?
      'SRCNUM',
      'RA',        => epoch 2000 J2000 ? GRAB FROM HIGH PM FILE BASED ON OBSID+SRCID
      'DEC',       => epoch 2000 J2000 ?
      'POSERR',   
      'LII',
      'BII',
 'N_OBSID',      ??
       'UVW2_SIGNIF',
       'UVM2_SIGNIF',
       'UVW1_SIGNIF',
       'U_SIGNIF',
       'B_SIGNIF',
       'V_SIGNIF',
 'UVW2_RATE',       ??
 'UVW2_RATE_ERR',   ??
 'UVM2_RATE',       ??
 'UVM2_RATE_ERR',   ??
 'UVW1_RATE',       ??
 'UVW1_RATE_ERR',   ??
 'U_RATE',          ??
 'U_RATE_ERR',      ??
 'B_RATE',          ??
 'B_RATE_ERR',      ??
 'V_RATE',          ??
 'V_RATE_ERR',      ??
       'UVW2_AB_FLUX',
       'UVW2_AB_FLUX_ERR',
       'UVM2_AB_FLUX',
       'UVM2_AB_FLUX_ERR',
       'UVW1_AB_FLUX',
       'UVW1_AB_FLUX_ERR',
       'U_AB_FLUX',
       'U_AB_FLUX_ERR',
       'B_AB_FLUX',
       'B_AB_FLUX_ERR',
       'V_AB_FLUX',
       'V_AB_FLUX_ERR',
       
       'UVW2_AB_MAG',       MEDIAN
       'UVW2_AB_MAG_ERR',   COMPUTE WEIGHTED ? 1-SIGMA ERROR
#       'UVW2_AB_MAG_MIN'
       'UVW2_AB_MAG_MAX'
       
       'UVM2_AB_MAG',
       'UVM2_AB_MAG_ERR',
       'UVM2_AB_MAG_MIN'
       'UVM2_AB_MAG_MAX'
       
       'UVW1_AB_MAG',
       'UVW1_AB_MAG_ERR',
       'UVW1_AB_MAG_MIN'
       'UVW1_AB_MAG_MAX'
       
       'U_AB_MAG',
       'U_AB_MAG_ERR',
       'U_AB_MAG_MIN'
       'U_AB_MAG_MAX'
       
       'B_AB_MAG',
       'B_AB_MAG_ERR',
       'V_AB_MAG',
       'V_AB_MAG_ERR',
       
       WHAT TO DO ABOUT THESE ?
 'UVW2_MAJOR_AXIS',
 'UVW2_MINOR_AXIS',
 'UVM2_MINOR_AXIS',
 'UVM2_MAJOR_AXIS',
 'UVW1_MAJOR_AXIS',
 'UVW1_MINOR_AXIS',
 'U_MAJOR_AXIS',
 'U_MINOR_AXIS',
 'B_MAJOR_AXIS',
 'B_MINOR_AXIS',
 'V_MAJOR_AXIS', 
 'V_MINOR_AXIS',
 'UVW2_POSANG',
 'UVM2_POSANG',
 'UVW1_POSANG',
 'U_POSANG',
 'B_POSANG',
 'V_POSANG',
              => BEST VALUES WEIGHTED BY EXPOSURE SOMEHOW
 'UVW2_QUALITY_FLAG',
 'UVM2_QUALITY_FLAG',
 'UVW1_QUALITY_FLAG',
 'U_QUALITY_FLAG',
 'B_QUALITY_FLAG',
 'V_QUALITY_FLAG',
 
 'UVW2_QUALITY_FLAG_ST',
 'UVM2_QUALITY_FLAG_ST',
 'UVW1_QUALITY_FLAG_ST',
 'U_QUALITY_FLAG_ST',
 'B_QUALITY_FLAG_ST',
 'V_QUALITY_FLAG_ST',
             => SELECT BEST QUALITY 
 'UVW2_EXTENDED_FLAG',
 'UVM2_EXTENDED_FLAG',
 'UVW1_EXTENDED_FLAG',
 'U_EXTENDED_FLAG',
 'B_EXTENDED_FLAG',
 'V_EXTENDED_FLAG',
             => NEW MERGED EXTENDED FLAGS
 'UVW2_SKY_IMAGE',
 'UVM2_SKY_IMAGE',
 'UVW1_SKY_IMAGE',
 'U_SKY_IMAGE',
 'B_SKY_IMAGE',
 'V_SKY_IMAGE'
 

                 'VAR_FLAG'
 ]

"""
bands=['UVW2','UVM2','UVW1','U','B','V']  # I think White is absent
   
def make_new(tx):
    # remove and add columns to main table
    tx.remove_columns(['N_SUMMARY','N_OBSID'])
    for  b in bands:
        tx.remove_columns([ b+'_RATE',b+'_RATE_ERR',b+'_SKY_IMAGE'])
        tx.remove_columns([ b+'_VEGA_MAG',b+'_VEGA_MAG_ERR'])
        tx.remove_columns([b+"_MAJOR_AXIS",b+"_MINOR_AXIS",b+"_POSANG"])
        tx.add_columns([0.,0.,None],
           names=[b+'_AB_FLUX_MAX',b+'_AB_MAG_MAX',b+'_VAR'])  
    return tx         


def fileio(infile,outstub="_singlerecs_v2",outdir=""):
    #
    # use the input file name infile to construct the output filename
    # open the output file
    # read in the table from the input file
    # return the table and the output file handle
    a = infile.rsplit('.')
    instub = ""
    ft = a[-1]
    for x in a[:-1]: 
        instub += x+"."
    instub = instub[:-1]    
    print (f"reading {infile} ...\n")
    if (ft.lower == 'csv'):
       t = ioascii.read(infile)  # table object
    elif (ft[:3] == 'fit'):
       f = fits.open(infile)
       x = f[1].data #[:tmax]  # LIMITED TABLE FOR TESTS
       summ = f[2].data # summary - has obstimes and exposures of obsids
       t = Table(x)   
       summ = Table(summ)
       x = ''
       f.close()
    if len(t) < 1:
        raise IOError("the catalogue table loading failed.\n")   
    outfile=instub+outstub+".csv"
    outf = open(outdir+outfile,'w')  # file handle

    # 
    # we use the astropy to flag data that have the same IAUNAME
    # after that we can call up the record/row selection by that.
    #
    t.add_index("IAUNAME")
    sources = t["IAUNAME"]
    sources = np.unique(sources)
    summ.add_index("OBSID") # to get a list of exposure times 
    # to find the table with rows that corresponds to a source, 
    # tab_src = t.loc[sources[k]] 
    #   the type is eather a astropy.table.row.Row for a single record, or
    #   astropy.table.table.Table if multiple records
    return t, outf, sources, summ

   
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
    es = np.asarray( np.sort(err[q]) )
    nk = np.int(len(es)/2)
    q2 = es <= es[nk] # (a)
    if np.sum(q2) <= 1:   # only one left to test 
        varFlag = ""
    else:    
        mean2 = mag[q][q2].mean()  # (b)
        sigma2 = mag[q][q2].std()  # (b)
        a = np.abs(mean2 - mag[q][q2]) > N * (np.sqrt(sigma2) + err[q][q2])
        varFlag = np.int(np.any(len(a) > 0))
    
    return med_mag, max_mag, sigma_mag, varFlag, mean_mag


def evaluate_extended_nature(tsrc,summ):
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
    obsids = tsrc['OBSID']
    tsumm = summ.loc[obsids] # find records for obsids 
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
        # TBD:
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
    #trow = Table(tablerow)
    #trow.remove_columns(['OBSID'])
    #col = trow.colnames
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
         #print (f"{60*'-'}\noriginal {qual_st}\n\ninteger flag values are: {aa}")    
          ranked_quals = aa =  ranked_qual_st2qual(qual_st) #np.array(aa)
          #print (f"summed rows: {aa}")
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
           aa.append(-9999)      
    return ma.masked_values(aa,-9999) 
 
########## main ()  
  
    
sussdir= '/Users/data/catalogs/suss_gaia_epic/'
suss='XMM-OM-SUSS5.0.fits'
tmax = 25000 # mximum number of records to read in
t1, outf, sources, summ = fileio(sussdir+suss)
t = make_new( t1) # replace make_record inside loop
#print(t,"\n")
col = t.colnames
nc = len(col)
for k in range(nc-1): 
       outf.write( f"{col[k]}," )
outf.write(f"{col[nc-1]}\n")       
#print("evaluating ",col," sources\n")
k9 = 0    

for src in sources[:]:   # do just 100 sources
#for src in sources[:5000]:   # do just 100 sources
    #print (f"number={k9}, source={src}"); 
    k9+=1
    tab_src = t.loc[src]
    tab_ma = ma.asarray(tab_src)
    nrow = tab_ma.size
    if type(tab_src) != astropy.table.row.Row: # is Table object
        new_row = tab_src[0]
    else: 
        new_row = tab_src    
        
    #  obsids = tab_src['OBSID']
    #  new_row = make_record( tab_src )
    #print (f"300 new_row is\n{new_row}; type is {type(new_row)}\n")
    base = evaluate_extended_nature(tab_src,summ)
    for band in bands:
    
        qf = tab_ma[band+"_QUALITY_FLAG_ST"]
        
        if nrow > 1: #not isinstance(qf, np.flexible):  # must be multiple elements 
            qual_st_out,ranked_quals = qual_st_merge(qf,err=tab_src[band+"_AB_MAG_ERR"])
            if len(qual_st_out) != 12:
                raise RuntimeError(f"448 ERROR: {k9} the quality flag is wrong {qual_st_out}")
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
        
    create_output_record(new_row,outf)     
outf.close()


    