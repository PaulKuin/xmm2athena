#!/bin/csh
#
#  remove some columns from the SUSS+Gaia cds matched file
#
set IN=$1
#
date
echo 
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe \
in=$IN ifmt=csv out=single_source_cat_trimmed.fits ofmt=fits \
cmd='delcols "RAdeg DEdeg errHalfMaj errHalfMin errPosAng SolID RandomI e_RAdeg e_DEdeg RADEcor RAPlxcor RApmRAcor RApmDEcor DEPlxcor DEpmRAcor DEpmDEcor PlxpmRAcor PlxpmDEcor pmRApmDEcor NAL NAC NgAL NbAL gofAL chi2AL epsi sepsi Solved APF nueff pscol e_pscol RApscolCorr DEpscolCorr PlxpscolCorr pmRApscolCorr pmDEpscolCorr MatchObsA Nper amax MatchObs IPDgofha IPDgofhp IPDfmp IPDfow o_Gmag FG e_FG RFG o_BPmag FBP e_FBP RFBP o_RPmag FRP e_FRP RFRP E(BP/RP) NBPcont NBPblend NRPcont NRPblend Mode n_RV o_RV o_RVd RVNper RVS/N RVgof RVTdur RVamp RVtempTeff RVtemplogg RVtemp[Fe/H] Vatmparam vbroad e_Vbroad o_Vbroad e_GRVSmag o_GRVSmag RVSS/N PQSO PGal PSS b_Teff_cds B_Teff_cdsa b_logg_cds B_logg_cdsa b_[Fe/H]_cds B_[Fe/H]_cdsa b_Dist_cds B_Dist_cdsa A0 b_A0_cds B_A0_cdsa b_AG_cds B_AG_cdsa b_E(BP-RP)_cds B_E(BP-RP)_cdsa Lib RAJ2000 DEJ2000 e_RAJ2000 e_DEJ2000 RADEcorJ2000 angDist"'
echo
date