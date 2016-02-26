#!/bin/sh
echo "Kunsu: Start makeFZ + Reco + Pico"
starver SL15k
job=$1
#run=${job: -1}
run=`echo $job | cut -f2 -d"_"`

# ---- Make folder
mkdir ./Files_$job
mkdir ./Files_$job/fzd
mkdir ./Files_$job/hft_reco
mkdir ./Files_$job/pile_up
mkdir ./Files_$job/picodst

# ---- Pile up file
#at=`perl -e 'srand; print int(rand(99)+1)'`
#cp -p /star/data01/pwg/kunsu/pileup/pileupSet$at/pile**.root ./Files_$job/pile_up/pile_up$at.root


# ---- Producing sim file .fzd
root4star -b -l <<EOF
.L starsim.hijing.gamma.C
starsim(1,$run,$RANDOM)
.q
EOF

mv pi0Dalitz_* ./Files_$job/fzd/.

echo "Kunsu: HFT reco starting"
# ---- HFT reconstruction
start=0
end=19
#inFile=Files_$job/fzd/hijing_pi0real_$run.starsim.fzd
inFile=Files_$job/fzd/pi0Dalitz_$run.starsim.fzd
#inFile=/star/u/kunsu/pwg/Pi0Raw/fz/gamma_$run.fzd   # pi0 real
#inFile=/star/u/kunsu/pwg/gamma/fz/gamma_$run.fzd    # gamma flat pt
inPile=Files_$job/pile_up/pile_up$at.root
chain=y2014a,event,McEvent,MuDst,tpc,fzin,sim_T,gen_T,geantout,tpcrs,TpcHitMover,TpxClu,evout,-HitFilt,FieldOn,AgML,usexgeom,MakeEvent,ITTF,Sti,NoSsdIt,NoSvtIt,StiHftC,pxlFastSim,pxlCluster,pxlHit,istFastSim,Idst,BAna,l0,Tree,logger,genvtx,tpcDB,bbcSim,btofsim,tags,emcY2,EEfs,evout,-dstout,IdTruth,big,McEvout,MiniMcMk,StiPulls,ReadAll,clearmem,McAna
echo $chain

pwd 
#env
root4star -b -l <<EOF
.x bfc.C(-1,"$chain","$inFile");
StPxlSimMaker* pxl = chain->GetMaker("pxlSimMaker");
pxl->useIdealGeom(); // ideal geometry
//pxl->useDbGeom();  // survey geometry
pxl->useRandomSeed();

//pxl->addPileup();
//pxl->setPileupFile("$inPile");


StMcAnalysisMaker* mcAnalysisMaker = (StMcAnalysisMaker*)chain->GetMaker("StMcAnalysisMaker");


chain->Init();
chain->EventLoop($start,$end);
chain->Finish();

EOF

mv pi0Dalitz_*.root Files_$job/hft_reco/.
mv mcAnalysis.pxlSimQa.root Files_$job/hft_reco/.

# ---- PicoDst
#root4star -l -b -q makePicoDst.C\($run,\"Files_$job/hft_reco/hijing_pi0real_$run.MuDst.root\",\"Files_$job/hft_reco/hijing_pi0real_$run.McEvent.root\"\)
#root4star -l -b -q makePicoDst.C\($run,\"Files_$job/hft_reco/pi0Dalitz_$run.MuDst.root\",\"Files_$job/hft_reco/pi0Dalitz_$run.McEvent.root\"\)
#mv *.picoDst.root Pi0Raw_hijing_sim_production_v0_$job.picoDst.root
#mv *.picoDst.root sim_production_v0_$job.picoDst.root

#privilges
find Files_$job/ -type d -exec chgrp rhstar {} \;
find Files_$job/ -type d -exec chmod g+rw {} \;
find Files_$job/ -type f -exec chgrp rhstar {} \;
find Files_$job/ -type f -exec chmod g+rw {} \;

# ---- Done bring files back
tar -cvf Files_$job.tar Files_$job
mv Files_$job/hft_reco/mcAnalysis.pxlSimQa.root /star/u/kunsu/pwg/Pi0Dalitz/qaHist/qa_$job.root
