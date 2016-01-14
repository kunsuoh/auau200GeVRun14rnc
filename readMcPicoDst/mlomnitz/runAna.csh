#!/bin/bash
starver dev
#root4star -b -q 'bfc.C('$2','$3',"pp2013a pxlRaw pxlDb pxlCluster pxlHit pxlQA mtd btof VFMinuit beamline BEmcChkStat Corr4 OSpaceZ2 OGridLeak3D -hitfilt","'$1'")'
echo $ROOTSYS
echo $LD_LIBRARY_PATH
cd /star/institutions/lbl/mlomnitz/D0_simu/Analysis
#root4star -b -q /star/institutions/ksu/lomnitz/SectorAlignmentProcedure_restructured/StRoot/StSectorAlign/sectorAlign.C\(\"/star/institutions/lbl_prod/hft/Run13/production/output/155/14155068/st_pxl_14155068_raw_*/*event.root\"\)>log
#root -b -q 'lMuDst.C' 'RecoSimGlobal.C+(100,"/star/institutions/lbl_prod/mlomnitz/hj_simu/Events/MuDst/hijing_'$1'.MuDst.root","testSIMU_'$1'")'
#Testin 
root -b -q 'lMuDst.C' 'RecoSimGlobal.C+(2000,"/star/institutions/lbl/mlomnitz/mlomnitz_prod/EffStudies/HaoFix/D0_MuDst/hijing_dzero_'$1'.MuDst.root","mcD0/testSIMU_'$1'")'
#tests
#root -b -q 'lMuDst.C' 'RecoGlobal.C+(2000,"/star/institutions/lbl_prod/mlomnitz/hj_simu/Events/*.MuDst.root","testSIMU")'>log
#from amilkar
#root4star -b -q /star/institutions/ksu/lomnitz/CosmicDataAlign_zeroField/StRoot/StSectorAlign/sectorAlign.C\(\"/star/u/aquinter/ksu/Simulation/CosmicRays/MisAlig/output/$1/*.event.root\"\)
#mv Alignment*\.align.root CosmicSimu_$1.root
#root4star -b -q /star/institutions/lbl_prod/hft/Run13/MatchedTree_Production/StRoot/StHftPool/HftMatchedTree/runHftTree.C\(\"st_pxl_14150049_raw_0780001.event.root\"\)
#echo $(basename *.root) > ./abc.dat
#rm -rf takealook.root
#abc=`echo $(basename *.root)`
#root4star -b -q '/star/institutions/lbl_prod/hft/Run13/test/test_shell.C("'$abc'")' 

#.x test_shell.C
