#!python3

# npmk March 2023

TOPCATPATH = "/Users/kuin/bin"
import os
import numpy as np

class fix_uvotssc(object):

    def __init__(self,):
        import numpy as np, os
        from astropy.io import fits
        self.filters=['UVW2','UVM2','UVW1','U','B','V']
        self.ABmVega_zp = {'UVW2':19.11-17.38, 'UVM2':18.54-16.85,\
            'UVW1':18.95-17.44, 'U':19.36-18.34, 'B':18.98-19.11, \
            'V':17.88-17.89}
        self.catalog = ""

    def _fix_phot(self, band,):
    
        catin = self.catalog 
        catout = 'tmp_bbb.fits'
        vega2ab = self.ABmVega_zp[band]
        command=f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tpipe "+\
           f"in={catin} ifmt=fits out={catout} ofmt=fits "+\
           f"cmd='addcol {band}_VEGAMAG {band}_ABMAG' "+\
           f"cmd='delcols {band}_ABMAG' "+\
           f"cmd='addcol {band}_ABMAG {vega2ab}+{band}_VEGAMAG' "
        print (command)
        status = os.system(command)
        if status > 0: 
            print(f"21 Error in {command} ")
        os.system(f'mv {catout} {catin}')
        
    def fix_catalog(self, catalogue):
        from astropy.io import fits
        
        self.catalog = catalogue
        for band in self.filters:
            self._fix_phot(band)
        f = fits.open(self.catalog,'update')
        hdr = f[1].header
        hdr['COMMENT'] = f" Updated Vega and AB magnitudes from level6 "
        f.flush()
        f.close()
        
#    def __exit__(self):
        
        
            
def photometry_fix():
    
    obx = fix_uvotssc()
    
    # photometry
    for k in range(24):
        cat = '/Users/data/catalogs/uvotssc2/' + f'ovotsso{k:02}.fits'
        obx.fix_catalog(cat)    
        print (f"fixed {cat} magnitudes AB, Vega")
        
# def finish( ):        
    # IAUNAME  UVOTSSC2 Jxxxxxx.x-xxxxxx  is TBD
    
    # MERGE NEARBY SOURCES 
    
    # SELECT ONLY SOURCES WITH AT LEAST ONE DETECTION OF 5 SIGMA  
        
        
def combine_input_files(in_cat='uvotssc2_in',chatter=0,just1=False,just2=False):
    """
    combine the epochal files into one big file for input to the single_source_cat_v2
    
    """
    import numpy as np
    import os
    root1 = '/Users/data/catalogs/uvotssc2/' 
    filenames=[]
    for k in np.arange(24): filenames.append(f'ovotsso{k:02}.fits')
    
    if not just2:
        command = f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tcat "
        for k in np.arange(len(filenames)):
            command += f"in={filenames[k]}#1 "
        command += f" ifmt=fits out={in_cat}_photometry.fits ofmt=fits-basic"
        if chatter > 3: print (f"{command}\n")
        status = os.system(command)
        if status > 0: print(f"photometry Error in {command} ")

    print (30*"=-")
    if not just1:
        command = f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tcat "
        for k in np.arange(len(filenames)):
            command += f"in={filenames[k]}.gz#2 "
        command += f" ifmt=fits out={in_cat}_summary.fits ofmt=fits-basic"

        if chatter > 3: print (f"{command}\n")
        status = os.system(command)
        if status > 0: print(f"summary Error in {command} ")

    print (30*"=-")

        
        
def combine_epochal_files(rootdir,out_cat="cat_out_all",chatter=0):
    """
    combine the epochal files into one big file for input to the single_source_cat_v2
    
    """
    import numpy as np
    import os
    
    # cd to the directory
    
    process_epoch_range = [2005.0,2022.1]# [2005.0,2022.2]
    process_epoch_step = 0.1
    epoch_array = np.arange(process_epoch_range[0],process_epoch_range[1],process_epoch_step)
    # so, we want to get epochs from epoch_array[k] to epoch_array[k]+process_epoch_step
    n_epochs = len(epoch_array) - 1
    
    filenames=[]
    for k in np.arange(n_epochs): filenames.append(f"cat_out_{epoch_array[k]:06.1f}.fits")
    k = 'cat_out_2007.7.fits'
    filenames.pop(filenames.index(k))

    command = f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tcat "
    for k in np.arange(len(filenames)):
        command += f"in={filenames[k]}#1 "
    command += f" ifmt=fits out={out_cat}_pm.fits ofmt=fits-basic"
            
    if chatter > 3: print (f"{command}\n")
    status = os.system(command)
    if status > 0: print(f"454 Error in {command} ")

    print (30*"=-")

    command = f"java -jar {TOPCATPATH}/topcat-full.jar -stilts tcat "
    for k in np.arange(len(filenames)): 
        command += f"in={filenames[k]}#2 "
    command += f" ifmt=fits out={out_cat}_nopm.fits ofmt=fits-basic"  
            
    if chatter > 3: print (f"{command}\n")
    status = os.system(command)
    if status > 0: print(f"454 Error in {command} ")
   
    print (30*"=-")
   