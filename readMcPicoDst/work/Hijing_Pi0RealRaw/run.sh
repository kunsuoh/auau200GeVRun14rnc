#!/bin/sh
echo "Kunsu: Start makeFZ + Reco + Pico"
starver SL15k
job=$1
#run=${job: -1}
run=`echo $job | cut -f2 -d"_"`

makeFolder=1
makeFZ=0
makeReco=1
makeRecoPileup=0
makePico=0
makeZip=0
makeQa=1
inputSource="gamma"

# ---- Make folder
if [ $makeFolder -eq 1 ]; then
echo "Kunsu: Make Folders"
mkdir ./Files_$job
mkdir ./Files_$job/fzd
mkdir ./Files_$job/hft_reco
mkdir ./Files_$job/pile_up
mkdir ./Files_$job/picodst
fi

# ---- Pile up file
if [ $makeRecoPileup -eq 1 ]; then
at=`perl -e 'srand; print int(rand(99)+1)'`
cp -p /star/data01/pwg/kunsu/pileup/pileupSet$at/pile**.root ./Files_$job/pile_up/pile_up$at.root
fi

# ---- Producing sim file .fzd
if [ $makeFZ -eq 1 ]; then
    echo "Kunsu: make sim file .fzd"
    root4star -b -l << EOF
    .L starsim.hijing.$inputSource.C
    starsim(9,$run,$RANDOM)
    .q
    EOF

    mv $inputSource_* ./Files_$job/fzd/.
else
    echo "Kunsu: Skip make sim file .fzd"
fi

if [ $makeReco -eq 1 ]; then
    echo "Kunsu: HFT reco starting"
    # ---- HFT reconstruction
    start=0
    end=19
    inFile=Files_$job/fzd/pi0Dalitz_$run.starsim.fzd
    if [ $makeFZ -eq 1 ]; then
        inFile=Files_$job/fzd/pi0Dalitz_$run.starsim.fzd
    else
        inFile=/star/u/kunsu/pwg/$inputSource/fz/$inputSource$run.fzd
    fi
    inPile=Files_$job/pile_up/pile_up$at.root
    if [ $makeQa -eq 1 ]; then
        chain=y2014a,event,McEvent,MuDst,tpc,fzin,sim_T,gen_T,geantout,tpcrs,TpcHitMover,TpxClu,evout,-HitFilt,FieldOn,AgML,usexgeom,MakeEvent,ITTF,Sti,NoSsdIt,NoSvtIt,StiHftC,pxlFastSim,pxlCluster,pxlHit,istFastSim,Idst,BAna,l0,Tree,logger,genvtx,tpcDB,bbcSim,btofsim,tags,emcY2,EEfs,evout,-dstout,IdTruth,big,McEvout,MiniMcMk,StiPulls,ReadAll,clearmem,McAna
    else
        chain=y2014a,event,McEvent,MuDst,tpc,fzin,sim_T,gen_T,geantout,tpcrs,TpcHitMover,TpxClu,evout,-HitFilt,FieldOn,AgML,usexgeom,MakeEvent,ITTF,Sti,NoSsdIt,NoSvtIt,StiHftC,pxlFastSim,pxlCluster,pxlHit,istFastSim,Idst,BAna,l0,Tree,logger,genvtx,tpcDB,bbcSim,btofsim,tags,emcY2,EEfs,evout,-dstout,IdTruth,big,McEvout,MiniMcMk,StiPulls,ReadAll,clearmem
    fi

    echo $chain
    pwd
    #env
    if [ -s .temprun.sh ]; then
        rm .temprun.sh
    fi
    echo "root4star -b -l <<EOF" > .temprun.sh
    echo ".x bfc.C(-1,"$chain","$inFile");" >> .temprun.sh
    echo "StPxlSimMaker* pxl = chain->GetMaker(\"pxlSimMaker\");" >> .temprun.sh
    echo "pxl->useIdealGeom(); // ideal geometry" >> .temprun.sh
    #echo "pxl->useDbGeom();  // survey geometry" >> .temprun.sh
    echo "pxl->useRandomSeed();" >> .temprun.sh
    if [ $makeRecoPileup -eq 1 ]; then
        echo "pxl->addPileup();" >> .temprun.sh
        echo "pxl->setPileupFile("$inPile");" >> .temprun.sh
    fi
    if [ $makeQa -eq 1 ]; then
        echo "StMcAnalysisMaker* mcAnalysisMaker = (StMcAnalysisMaker*)chain->GetMaker(\"StMcAnalysisMaker\");" >> .temprun.sh
    fi
    echo "chain->Init();" >> .temprun.sh
    echo "chain->EventLoop($start,$end);" >> .temprun.sh
    echo "chain->Finish();" >> .temprun.sh
    echo "EOF" >> .temprun.sh

    chmod +x .temprun.sh
    ./.temprun.sh
    #rm .temprun.sh

    mv $inputSource_*.root Files_$job/hft_reco/.
    if [ $makeQa -eq 1 ]; then
        mv mcAnalysis.pxlSimQa.root Files_$job/hft_reco/.
    fi
fi
if [ $makePico -eq 1 ]; then
    # ---- PicoDst
    root4star -l -b -q makePicoDst.C\($run,\"Files_$job/hft_reco/$inputSource_$run.MuDst.root\",\"Files_$job/hft_reco/$inputSource_$run.McEvent.root\"\)
    mv *.picoDst.root $inputSource_$job.picoDst.root
fi

#privilges
find Files_$job/ -type d -exec chgrp rhstar {} \;
find Files_$job/ -type d -exec chmod g+rw {} \;
find Files_$job/ -type f -exec chgrp rhstar {} \;
find Files_$job/ -type f -exec chmod g+rw {} \;

# ---- Done bring files back
if [ $makeZip -eq 1 ]; then
    tar -cvf Files_$job.tar Files_$job
fi

if [ $makeQa -eq 1 ]; then
    mv Files_$job/hft_reco/mcAnalysis.pxlSimQa.root /star/u/kunsu/pwg/$inputSource/qaHist/qa_$job.root
fi