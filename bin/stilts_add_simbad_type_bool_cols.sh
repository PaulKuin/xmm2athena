#!/bin/csh
if ($2 == "") then
   echo This script from Paul is to add Simbad columns 
   echo parameters: INFILE OUTFILE
endif
#alias stilts='java -jar /Users/kuin/bin/topcat-full.jar -stilts'
#echo $1
#echo $2
#echo $3
#
set IN=$1
set OUT=$2
#set SHORTOBSID=`echo $OBSID1 | bc`
#
#echo java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$IN out=$OUT cmd='select matches(OBSID,padWithZeros('${SHORTOBSID}',10))'
#echo java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$IN out=$OUT cmd='addcol AGN equals(main_type,"AGN")'
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$IN ifmt=fits out=$OUT ofmt=fits \
cmd='addcol AGN equals(main_type,\"AGN\")|equals(main_type,\"AGN_Candidate\")' \
cmd='addcol BLLac equals(main_type,\"BLLac\")' \
cmd='addcol Blazar equals(main_type,\"Blazar\")|equals(main_type,\"Blazar_Candidate\")' \
cmd='addcol QSO equals(main_type,\"QSO\")|equals(main_type,\"QSO_Candidate\")' \
cmd='addcol SN equals(main_type,\"Supernova\")|equals(main_type,\"Supernova_Candidate\")'\
cmd='addcol Galaxy equals(main_type,\"Galaxy\")|equals(main_type,\"GinCl\")|equals(main_type,\"Galaxy_Candidate\")|equals(main_type,\"GinPair\")|equals(main_type,\"BlueCompactG\")'\
cmd='addcol Seyfert equals(main_type,\"Seyfert\")|equals(main_type,\"Seyfert1\")|equals(main_type,\"Seyfert2\")|equals(main_type,\"StarBurstG\")' \
cmd='addcol LINER equals(main_type,\"LINER\")' \
cmd='addcol RadioG equals(main_type,\"RadioG\")' \
cmd='addcol AeBe equals(main_type,\"Ae*\")|equals(main_type,\"Ae*_Candidate\")|equals(main_type,\"Be*\")|equals(main_type,\"Be*_Candidate\")' \
cmd='addcol BlueStraggler equals(main_type,\"BlueStraggler\")|equals(main_type,\"BlueStraggler_Candidate\")' \
cmd='addcol BlueSG equals(main_type,\"BlueSG\")|equals(main_type,\"BlueSG_Candidate\")' \
cmd='addcol YellowSG equals(main_type,\"YellowSG\")|equals(main_type,\"YellowSG_Candidate\")' \
cmd='addcol RedSG equals(main_type,\"RedSG\")|equals(main_type,\"RedSG_Candidate\")' \
cmd='addcol CarbonStar equals(main_type,\"C*\")' \
cmd='addcol S_Star equals(main_type,\"S*\")|equals(main_type,\"S*_Candidate\")' \
cmd='addcol Pec equals(main_type,\"ChemPec*\")' \
cmd='addcol IrregularV equals(main_type,\"IrregularV*\")' \
cmd='addcol Star equals(main_type,\"Star\")' \
cmd='addcol WR equals(main_type,\"WolfRayet*\")|equals(main_type,\"WolfRayet*_Candidate\")' \
cmd='addcol Supergiant equals(main_type,\"Supergiant\")|equals(main_type,\"Supergiant_Candidate\")' \
cmd='addcol Binary equals(main_type,\"**\")' \
cmd='addcol SB equals(main_type,\"SB*\")|equals(main_type,\"SB*_Candidate\")' \
cmd='addcol EclBin equals(main_type,\"EclBin\")|equals(main_type,\"EclBin_Candidate\")' \
cmd='addcol XB equals(main_type,\"HighMassXBin\")|equals(main_type,\"LowMassXBin\")|equals(main_type,\"XrayBin\")|equals(main_type,\"HighMassXBin_Candidate\")|equals(main_type,\"LowMassXBin_Candidate\")|equals(main_type,\"XrayBin_Candidate\")' \
cmd='addcol EmLine equals(main_type,\"EmLine*\")' \
cmd='addcol CV equals(main_type,\"CataclyV*\")|equals(main_type,\"CataclyV*_Candidate\")' \
cmd='addcol Eruptive equals(main_type,\"Eruptive*\")' \
cmd='addcol Nova equals(main_type,\"Nova\")|equals(main_type,\"Nova_Candidate\")' \
cmd='addcol Symbiotic equals(main_type,\"Symbiotic*\")|equals(main_type,\"Symbiotic*_Candidate\")' \
cmd='addcol LowMassStar equals(main_type,\"Low-Mass*\")' \
cmd='addcol PulsV equals(main_type,\"PulsV*\")' \
cmd='addcol RRLyrae equals(main_type,\"RRLyrae\")|equals(main_type,\"RRLyrae_Candidate\")' \
cmd='addcol Mira equals(main_type,\"Mira\")|equals(main_type,\"Mira_Candidate\")'\
cmd='addcol Pulsar equals(main_type,\"Pulsar\")' \
cmd='addcol LPV equals(main_type,\"LongPeriodV*\")|equals(main_type,\"LongPeriodV*_Candidate\")' \
cmd='addcol Cepheid equals(main_type,\"Cepheid\")|equals(main_type,\"Cepheid_Candidate\")|equals(main_type,\"ClassicalCep\")' \
cmd='addcol BYDraV equals(main_type,\"BYDraV*\")' \
cmd='addcol RSCVnV equals(main_type,\"RSCVnV*\")' \
cmd='addcol RVTauV equals(main_type,\"RVTauV*\")' \
cmd='addcol RotV equals(main_type,\"RotV*\")' \
cmd='addcol alf2CVnV equals(main_type,\"alf2CVnV*\")' \
cmd='addcol Transient equals(main_type,\"Transient\")|equals(main_type,\"gammaBurst\")|equals(main_type,\"radioBurst\")' \
cmd='addcol bCepV equals(main_type,\"bCepV*\")' \
cmd='addcol delSctV equals(main_type,\"delSctV*\")' \
cmd='addcol gammaDorV equals(main_type,\"gammaDorV*\")' \
cmd='addcol HB equals(main_type,\"HorBranch*\")|equals(main_type,\"HorBranch*_Candidate\")' \
cmd='addcol AGB equals(main_type,\"AGB*\")|equals(main_type,\"AGB*_Candidate\")' \
cmd='addcol Assoc equals(main_type,\"Association\")' \
cmd='addcol OpenCl equals(main_type,\"OpenCluster\")' \
cmd='addcol StarCluster equals(main_type,\"Cluster*\")|equals(main_type,\"Cluster*_Candidate\")' \
cmd='addcol GlobCl equals(main_type,\"GlobCluster\")|equals(main_type,\"Glob_Cluster_Candidate\")' \
cmd='addcol postAGB equals(main_type,\"post-AGB*\")|equals(main_type,\"post-AGB*_Candidate\")' \
cmd='addcol PN equals(main_type,\"PlanetaryNeb\")|equals(main_type,\"PlanetaryNeb_Candidate\")' \
cmd='addcol Planet equals(main_type,\"Planet\")|equals(main_type,\"Planet_Candidate\")' \
cmd='addcol HighPM equals(main_type,\"HighPM*\")' \
cmd='addcol HighVel equals(main_type,\"HighVel*\")' \
cmd='addcol NS equals(main_type,\"Neutron*_Candidate\")' \
cmd='addcol ULX equals(main_type,\"ULX\")|equals(main_type,\"ULX_Candidate\")' \
cmd='addcol HotSubdwarf equals(main_type,\"HotSubdwarf\")|equals(main_type,\"HotSubdwarf_Candidate\")' \
#cmd='addcol DarkNeb equals(main_type,\"DarkNeb\")' \
cmd='addcol YSO equals(main_type,\"YSO\")|equals(main_type,\"YSO_Candidate\")' \
cmd='addcol TTauri equals(main_type,\"TTauri*\")|equals(main_type,\"TTauri*_Candidate\")' \
cmd='addcol OrionV equals(main_type,\"OrionV*\")' \
cmd='addcol VariableStar equals(main_type,\"Variable*\")|equals(main_type,\"Variable*_Candidate\")' \
cmd='addcol WD equals(main_type,\"WhiteDwarf\")|equals(main_type,\"WhiteDwarf_Candidate\")' 
# add various col as not any of the above
#old: !(AGN|BLLac|Blazar|QSO|SN|Galaxy|AeBe|Binary|SB|EclBin|XB|EmLine|CV|Eruptive|Nova|Symbiotic|LowMassStar|PulsV|RRLyrae|Mira|Pulsar|LPV|Cepheid|BYDraV|RSCVnV|RotV|alf2CVnV|Transient|bCepV|delSctV|gammaDorV|NS|ULX|HotSubdwarf|YSO|OrionV|Variable|WD)
#!(AGN|Assoc|BLLac|Blazar|QSO|SN|Galaxy|Seyfert|LINER|RadioG|AeBe|BlueStraggler|BlueSG|YellowSG|RedSG|CarbonStar|S_Star|Pec|IrregularV|Star|WR|Supergiant|StarCluster|GlobCl|Binary|SB|EclBin|XB|EmLine|CV|Eruptive|Nova|Symbiotic|LowMassStar|PulsV|RRLyrae|Mira|Pulsar|LPV|Cepheid|BYDraV|RSCVnV|RVTauV|RotV|alf2CVnV|Transient|bCepV|delSctV|gammaDorV|HB|OpenCl|NS|ULX|Star|HotSubdwarf|YSO|TTauri|OrionV|VariableStar|WD|postAGB)
# 
#set OTHER=\!\(AGN\|Assoc\|BLLac\|Blazar\|QSO\|SN\|Galaxy\|Seyfert\|LINER\|RadioG\|AeBe\|BlueStraggler\|BlueSG\|YellowSG\|
#RedSG\|CarbonStar\|S_Star\|Pec\|IrregularV\|Star\|WR\|SuperGiant\|StarCluster\|GlobCl\|
#Binary\|SB\|EclBin\|XB\|EmLine\|CV\|Eruptive\|Nova\|Symbiotic\|LowMassStar\|PulsV\|RRLyrae\|Mira\|Pulsar\|LPV\|Cepheid\|
#BYDraV\|RSCVnV\|RVTauV\|RotV\|alf2CVnV\|Transient\|bCepV\|delSctV\|gammaDorV\|HB\|OpenCl\|
#NS\|ULX\|Star\|HotSubdwarf\|YSO\|TTauri\|OrionV\|VariableStar\|WD\)

