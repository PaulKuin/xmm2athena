#!/bin/csh
echo
echo make auxiliary data ready for test data set
echo this is a restricted data set, so not all surveys are used
echo
echo run the match with CDS catalogues
echo
~/bin/stilts_match_suss_to_PS_SM_UKIDDS_VISTA_AllWISE.csh
echo
echo now start merging with VISTA
echo
/bin/csh ~/bin/vista_merge.csh
echo
echo we don't need to run UKIRT_merge.csh, UKIRT_augment.csh and only part of 
echo Vega_to_AB.csh
echo
echo next is concat_data_v2.csh though this needs to skip the UKIDSS part
/bin/csh ~/bin/concat_data_v2.csh
echo
