#!/bin/sh
echo "Kunsu: Start makeFZ + Reco + Pico"
starver SL15k
job=$1
#run=${job: -1}
run=`echo $job | cut -f2 -d"_"`
inputSource=$2

makeFolder=1
makeReco=0
makeRecoPileup=1
makePico=1
makeZip=0
makeQa=0


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
if [ -s /star/u/kunsu/pwg/${inputSource}/fz/${inputSource}_${run}.starsim.fzd ]; then
echo "Kunsu: Skip make sim file .fzd"
else
echo "Kunsu: Make sim file .fzd"
root4star -b -l << EOF
.L starsim.hijing.$inputSource.C
starsim(19,$run,$RANDOM)
.q
EOF
mv ${inputSource}_* ./Files_$job/fzd/.
cp Files_$job/fzd/* /star/u/kunsu/pwg/$inputSource/fz/

fi

# ---- TPC reconstruction
if [ -s /star/u/kunsu/pwg/${inputSource}/tpc_reco/${inputSource}_${run}.event.root ]; then
echo "Kunsu: Skip make tpc_reco files"
else
echo "Kunsu: Make tpc_reco files"
start=0
end=19
if [ -s /star/u/kunsu/pwg/${inputSource}/fz/${inputSource}_${run}.starsim.fzd ]; then
inFile=/star/u/kunsu/pwg/${inputSource}/fz/${inputSource}_${run}.starsim.fzd
else
inFile=Files_$job/fzd/${inputSource}_$run.starsim.fzd
fi
chain=y2014a,event,tpc,fzin,geantout,tpcrs,TpcHitMover,TpxClu,evout,-HitFilt

root4star -b -l <<EOF
.x bfc.C($start,$end,"$chain","$inFile");
EOF
mv ${inputSource}_${run}* /star/u/kunsu/pwg/${inputSource}/tpc_reco/
fi

# ---- HFT reconstruction
inFile2=/star/u/kunsu/pwg/${inputSource}/tpc_reco/${inputSource}_${run}.event.root
inPile=Files_$job/pile_up/pile_up$at.root
if [ $makeQa -eq 1 ]; then
chain2=y2014a,event,McEvent,MuDst,in,sim_T,gen_T,geantout,evout,FieldOn,AgML,usexgeom,MakeEvent,ITTF,Sti,NoSsdIt,NoSvtIt,StiHftC,pxlFastSim,istFastSim,Idst,BAna,l0,Tree,logger,genvtx,tpcDB,bbcSim,btofsim,emcY2,EEfs,evout,-dstout,IdTruth,big,McEvout,MiniMcMk,McAna,ReadAll,clearmem,pxlCluster,pxlHit
else
chain2=y2014a,event,McEvent,MuDst,in,sim_T,gen_T,geantout,evout,FieldOn,AgML,usexgeom,MakeEvent,ITTF,Sti,NoSsdIt,NoSvtIt,StiHftC,pxlFastSim,istFastSim,Idst,BAna,l0,Tree,logger,genvtx,tpcDB,bbcSim,btofsim,emcY2,EEfs,evout,-dstout,IdTruth,big,McEvout,MiniMcMk,ReadAll,clearmem,pxlCluster,pxlHit
fi
echo $chain
pwd
#env
if [ -s .temprun.sh ]; then
rm .temprun.sh
fi
echo "starver SL15k" > .temprun.sh
echo "root4star -b -l <<EOF" > .temprun.sh
echo ".x bfc.C(-1,\"$chain2\",\"$inFile2\");" >> .temprun.sh
echo "StPxlSimMaker* pxl = chain->GetMaker(\"pxlSimMaker\");" >> .temprun.sh
#echo "pxl->useIdealGeom(); // ideal geometry" >> .temprun.sh
echo "pxl->useDbGeom();  // survey geometry" >> .temprun.sh
#echo "pxl->setFastSim();" >> .temprun.sh
echo "pxl->setFastSimRaw();" >> .temprun.sh
echo "pxl->setWrongRowRatio(0.5,0.5);" >> .temprun.sh # -1 : off, other : ratio
echo "pxl->useRandomSeed();" >> .temprun.sh
if [ $makeRecoPileup -eq 1 ]; then
echo "pxl->addPileup();" >> .temprun.sh
echo "pxl->setPileupFile(\"$inPile\");" >> .temprun.sh
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

mv ${inputSource}_*.root Files_$job/hft_reco/.
if [ $makeQa -eq 1 ]; then
mv mcAnalysis.pxlSimQa.root Files_$job/hft_reco/.
fi




if [ $makeReco -eq 1 ]; then
    echo "Kunsu: HFT reco starting"
    start=0
    end=19
    if [ -s /star/u/kunsu/pwg/${inputSource}/fz/${inputSource}_${run}.starsim.fzd ]; then
     inFile=/star/u/kunsu/pwg/${inputSource}/fz/${inputSource}_${run}.starsim.fzd
    else
        inFile=Files_$job/fzd/${inputSource}_$run.starsim.fzd
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
    echo "starver SL15k" > .temprun.sh
    echo "root4star -b -l <<EOF" > .temprun.sh
    echo ".x bfc.C(-1,\"$chain\",\"$inFile\");" >> .temprun.sh
    echo "StPxlSimMaker* pxl = chain->GetMaker(\"pxlSimMaker\");" >> .temprun.sh
    echo "pxl->useIdealGeom(); // ideal geometry" >> .temprun.sh
    #echo "pxl->useDbGeom();  // survey geometry" >> .temprun.sh
    #echo "pxl->setFastSim();" >> .temprun.sh
    echo "pxl->setFastSimRaw();" >> .temprun.sh
    echo "pxl->setWrongRowRatio(0.,0.);" >> .temprun.sh # -1 : off, other : ratio
    echo "pxl->useRandomSeed();" >> .temprun.sh
    if [ $makeRecoPileup -eq 1 ]; then
        echo "pxl->addPileup();" >> .temprun.sh
        echo "pxl->setPileupFile(\"$inPile\");" >> .temprun.sh
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

    mv ${inputSource}_*.root Files_$job/hft_reco/.
    if [ $makeQa -eq 1 ]; then
        mv mcAnalysis.pxlSimQa.root Files_$job/hft_reco/.
    fi
fi
if [ $makePico -eq 1 ]; then
    # ---- PicoDst
    root4star -l -b -q makePicoDst.C\($run,\"Files_$job/hft_reco/${inputSource}_$run.MuDst.root\",\"Files_$job/hft_reco/${inputSource}_$run.McEvent.root\",\"$inputSource\"\)
    mv *.picoDst.root ${inputSource}_$job.picoDst.root
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