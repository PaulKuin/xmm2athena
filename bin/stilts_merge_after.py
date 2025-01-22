#python
from astropy.io import fits,ascii as ioasci
pwd
g = open('column_header_names.lis')
names=[]
for x in g:
    names.append(x.split('\n')[0] )
g.close()
f = fits.open('suss5.0_gaiaedr3_0690744101.fits')

from astropy.table import Table

l1=open('list1.txt')
l2=open('list2.txt')
mkdir table

l1n = l1.readlines()
l1.close()

k = 0
for fnam in l1n:
     fnam = fnam.split("\n")[0]
     f = fits.open(fnam)
     d = f[1].data
     t = Table(d)
     s = t[names]
     out = "table/suss5."+fnam.split('.')[1]+".txt"
     print (f"{k} - processing {fnam}; output is {out}")
     ioasci.write(s,output=out,names=names)
     f.close()
     k += 1
     
l2n = l2.readlines()
l2.close()
k = 0
for fnam in l2n:
     fnam = fnam.split("\n")[0]
     f = fits.open(fnam)
     d = f[1].data
     t = Table(d)
     s = t[names]
     out = "table/suss5."+fnam.split('.')[1]+".txt"
     print (f"{k} - processing {fnam}; output is {out}")
     ioasci.write(s,output=out,names=names)
     f.close()
     k+=1
