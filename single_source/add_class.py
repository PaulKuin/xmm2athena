#!python
#
# add class column values after running add_class.csh has created class0, class1, class2
#
from astropy.table import Table
import numpy as np
file='temp4python.fits'
t = Table.read(file)
t.add_column(np.nan,name='class')
t['class'][t['class0'] == 0] = 0
t['class'][t['class1'] == 1] = 1
t['class'][t['class2'] == 2] = 2
# remove some that are in both star and QSO
t['class'][(t['class0'] == 0) & (t['class2']==1)] = np.nan
# remove some that are in both star and galaxy 
t['class'][(t['class0'] == 0) & (t['class2']==2)] = np.nan
t.write(file,overwrite=True)            