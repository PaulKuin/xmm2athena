#/bin/bash
""""
  merging table fits files 
  works only if number of columns in tables match exactly

""""
split -l 2000 list1.txt list1_
echo "written `ls list1_??` - now processing aa to ae "
ftmerge @list1_aa merged_all3.fits columns=@column_header_names.lis clobber=yes chatter=5
echo merged_all3.fits >> list1_ab
ftmerge @list1_ab merged_all4.fits columns=@column_header_names.lis clobber=yes chatter=5
echo merged_all4.fits >> list1_ac
ftmerge @list1_ac merged_all3.fits columns=@column_header_names.lis clobber=yes chatter=5
echo merged_all3.fits >> list1_ad
ftmerge @list1_ad merged_all4.fits columns=@column_header_names.lis clobber=yes chatter=5
echo merged_all4.fits >> list1_ae
ftmerge @list1_ae merged_all_list1.fits columns=@column_header_names.lis clobber=yes chatter=5
#echo merged_all3.fits >> list1_af
#ftmerge @list1_af columns=@column_header_names.lis merged_all4.fits columns=@column_header_names.lis clobber=yes chatter=5
#rm list1_*
#
