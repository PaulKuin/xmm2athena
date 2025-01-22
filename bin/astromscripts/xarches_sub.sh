#!/bin/csh 
# {output} 
#  using ficed file names for xarches
###############################################################################
# Name:
# Description:
# Input files:
# Output files:
# main: xmm2athena.gaia_cat()
###############################################################################
set OUT=$1

# Start server first
#echo first execute "archesxmatch.bash i"

# upload files
archesxmatch.bash p gtemp.fit
archesxmatch.bash p gtemp.fit
archesxmatch.bash p suss5_obsid.fit
archesxmatch.bash x test1_arches.txt
# list and download result
archesxmatch.bash l
archesxmatch.bash g xmerg.fits

# rename 
mv xmerg.fits $OUT

# clean 
archesxmatch.bash r gtemp.fit
archesxmatch.bash r suss5_obsid.fit
archesxmatch.bash r xmerg.fits
archesxmatch.bash r gaia.vot
archesxmatch.bash r suss.vot
