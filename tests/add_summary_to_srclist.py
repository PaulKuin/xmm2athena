#!python
file1 = "test_cat.fit"
file2 = "test_cat_summary.fits"
from astropy.io import fits
f1 = fits.open(file1)
f2 = fits.open(file2)
f1.append(f2[1])
f1[1].name="SRCLIST"
f1[2].name="SUMMARY"
f1.writeto(file1,overwrite=True)