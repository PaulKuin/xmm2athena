#!python3

# NPMK, 2023-05-16 (npkuin @ gmail) distance nearby sources UVOTSSC/SUSS

"""
distance_nearby_source.py <infile> <outfile>


"""

HEALPIX = {"scheme":"nested", "depth":17}
nside = 2**HEALPIX["depth"]

#import optparse
import numpy as np
#import numpy.ma as ma
from cdshealpix import nested
from astropy.io import fits
from astropy.table import Table
import os, sys 
from astropy.coordinates import SkyCoord
import astropy.units as u


class healpix_stuff(object):

    def __init__(self, infile, out="", ext=1):
        import numpy as np
        from astropy.io import fits
        from astropy.table import Table
        from cdshealpix import nested
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        fcat = Table.read(infile,ext)
        self.t = fcat['IAUNAME','RA','DEC']
        self.ext=ext
        self.depth=HEALPIX["depth"]
        #self.t = Table.read(self.cat)

    def add_healpixno(self):
        ra = self.t['RA'] 
        dec = self.t['DEC'] 
        self.pos = SkyCoord(ra, dec , unit='deg', frame='icrs')
        self.ipix = nested.skycoord_to_healpix(self.pos,self.depth)
        self.t.add_column(self.ipix,name='healpix')
        self.t.sort('healpix')
       
    def get_neighbours(self):
        pixes = nested.neighbours(self.ipix,self.depth)
        return pixes
        
    def get_records(self,ipix):
        xloc = self.t['healpix'] == ipix 
        return self.t[xloc]    
        
    def get_header(self):
        return self.t.colnames   
        
    def __exit__(self):
        self.cat.close()
        del self.cat
        

def distance(ra1,dec1,ra2,dec2):
    c1 = SkyCoord(ra1*u.deg,dec1*u.deg,frame='icrs')
    c2 = SkyCoord(ra2*u.deg,dec2*u.deg,frame='icrs')
    sep = c1.separation(c2)
    return sep.arcsecond  
    
def dist(pos1, ra2,dec2):
    pos2 = SkyCoord(ra2*u.deg,dec2*u.deg,frame='icrs')
    sep = pos1.separation(pos2)
    return sep.arcsecond    
                
        
 ####################################################       
        
if __name__ == '__main__': 
    #usage = "usage: %prog infile outfile"
    #epilog = "required input is input catalogue and output catalogue"
    #parser = optparse.OptionParser(usage=usage,epilog=epilog)

    infile = sys.argv[1]
    outfile = sys.argv[2]

else: 

    #raise RuntimeError ("need to run from command line ")
    infile = input("input file")
    outfile = input("output file")

print (f"input file = {infile}\noutput file = {outfile}")

hp = healpix_stuff(infile, out=outfile, ext=1)

# first compute the healpix coordinate for each source

hp.add_healpixno()


# secondly determine for each source healpix the neighboring healpix

pixes = hp.get_neighbours()

separations = []
nrec = len(pixes[:,0])
for k in np.arange(nrec):   # each object
    match1 = []
    for k1 in range(9): 
        inxpix = hp.get_records( pixes[k,k1] ) 
        for x in inxpix:
            match1.append(x)
    tsep = []    
    
    for mt in match1:
        ra2, dec2 = mt['RA'],mt['DEC']
        dis1 = dist (hp.pos, ra2, dec2)
        qdis = dis1 > 1e-3
        if qdis.sum() > 0:
           tsep.append(dis1[qdis])
        adis = np.min(np.asarray(tsep).flatten())
    separations.append(adis)

hp.t.add_column(adis*u.arcsec,name='nearest_distance')

# write out table  

hp.t.write(outfile)






