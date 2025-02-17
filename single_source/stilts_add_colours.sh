#!/bin/csh 
#
#set IN = /Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_xSDSSxGaia_v0.1.csv
#set OUT = /Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_xSDSSxGaia_v0.2.csv
#
echo adding colours
#
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe \
  in=/Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_xSDSSxGaia_v0_1.csv \
  out=/Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_xSDSSxGaia_v0_2.csv \
  ifmt='csv' ofmt='csv' \ 
  cmd="addcol uvw2_uvw1 UVW2_AB_MAG-UVW1_AB_MAG; \
  addcol uvw1_u UVW1_AB_MAG-U_AB_MAG;\
  addcol b_v B_AB_MAG-V_AB_MAG; \
  addcol g_r gmag-rmag; \
  addcol r_i rmag-imag;\
  addcol i_z imag-zmag "