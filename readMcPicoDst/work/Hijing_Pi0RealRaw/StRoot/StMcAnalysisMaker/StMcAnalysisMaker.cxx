#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TSystem.h"

#include "StParticleDefinition.hh"
#include "StTrack.h"
#include "StPrimaryTrack.h"
#include "StTrackGeometry.h"
#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcContainers.hh"
#include "SystemOfUnits.h"
#include "StEvent.h"
#include "StTrackNode.h"
#include "StGlobalTrack.h"
#include "StEvent/StTpcHit.h"
#include "StEvent/StTrackTopologyMap.h"

#include "StBFChain/StBFChain.h"
#include "StMcEvent.hh"
#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"
#include "StMcEvent/StMcTpcHit.hh"
#include "StMcAnalysisMaker.h"
#include "StMcTpcHitCollection.hh"

#include "StEvent/StEventTypes.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StEvent/StEventSummary.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTpcDedxPidAlgorithm.h"
#include "StMcEvent/StMcVertex.hh"
#include "StarClassLibrary/StParticleTypes.hh"

#include "StDetectorDbMaker/StDetectorDbTriggerID.h"

ClassImp(StMcAnalysisMaker);

StMcAnalysisMaker::StMcAnalysisMaker(const char *name, const char *title): StMaker(name), mFile(NULL), mTree(NULL)
{
    cout << "StMcAnalysisMaker::StMcAnalysisMaker - DONE" << endl;
}
//__________________________________
int StMcAnalysisMaker::Init()
{
    StBFChain *bfChain = (StBFChain *) StMaker::GetChain();
    
    /*if (!bfChain) return kStFatal;
     TString fileName( gSystem->BaseName(bfChain->GetFileOut().Data()) );
     fileName = fileName.ReplaceAll(".event.root","");
     */
    
    if(!mOutFileName.Length()) mOutFileName = "mcAnalysis";
    mOutFileName = mOutFileName.ReplaceAll(".root","");
    
    mFile = new TFile(Form("%s.pxlSimQa.root",mOutFileName.Data()), "recreate");
    assert(mFile && !mFile->IsZombie());
    
    
    mAssoc = (StAssociationMaker*)GetMaker("StAssociationMaker");
    
    if (!mAssoc)
    {
        cout << " empty StAssociateMaker, stop!!" << endl;
        exit(1);
    }
    
    // TTree
    mTree = new TTree("mTree","electron pair tree for QA");
    mTree->Branch("nPair",&nPair,"nPair/I");
    mTree->Branch("nMcPxl1Hits",&nMcPxl1Hits,"nMcPxl1Hits/I");
    mTree->Branch("nMcPxl2Hits",&nMcPxl2Hits,"nMcPxl2Hits/I");
    mTree->Branch("nMcIstHits", &nMcIstHits, "nMcIstHits/I");
    mTree->Branch("nRcPxl1Hits",&nRcPxl1Hits,"nRcPxl1Hits/I");
    mTree->Branch("nRcPxl2Hits",&nRcPxl2Hits,"nRcPxl2Hits/I");
    mTree->Branch("nRcIstHits", &nRcIstHits, "nRcIstHits/I");
    mTree->Branch("pairPt", &pairPt, "pairPt[nPair]/F");
    mTree->Branch("pairEta", &pairEta, "pairEta[nPair]/F");
    mTree->Branch("openangle", &openangle, "openangle[nPair]/F");
    mTree->Branch("mcopenangle", &mcopenangle, "mcopenangle[nPair]/F");
    mTree->Branch("mcDist_pxl1", &mcDist_pxl1, "mcDist_pxl1[nPair]/F");
    mTree->Branch("mcDist_pxl2", &mcDist_pxl2, "mcDist_pxl2[nPair]/F");
    mTree->Branch("mcDist_ist", &mcDist_ist, "mcDist_ist[nPair]/F");
    mTree->Branch("rcDist_pxl1", &rcDist_pxl1, "rcDist_pxl1[nPair]/F");
    mTree->Branch("rcDist_pxl2", &rcDist_pxl2, "rcDist_pxl2[nPair]/F");
    mTree->Branch("rcDist_ist", &rcDist_ist, "rcDist_ist[nPair]/F");
    mTree->Branch("convR", &convR, "convR[nPair]/F");
    mTree->Branch("parentGid", &parentGid, "parentGid[nPair]/I");
    mTree->Branch("mass", &mass, "mass[nPair]/F");
    mTree->Branch("pairDca", &pairDca, "pairDca[nPair]/F");
    mTree->Branch("pt1", &pt1, "pt1[nPair]/F");
    mTree->Branch("pt2", &pt2, "pt2[nPair]/F");
    mTree->Branch("eta1", &eta1, "eta1[nPair]/F");
    mTree->Branch("eta2", &eta2, "eta2[nPair]/F");
    mTree->Branch("clusterSize1_pxl1", &clusterSize1_pxl1, "clusterSize1_pxl1[nPair]/I");
    mTree->Branch("clusterSize1_pxl2", &clusterSize1_pxl2, "clusterSize1_pxl2[nPair]/I");
    mTree->Branch("clusterSize2_pxl1", &clusterSize2_pxl1, "clusterSize2_pxl1[nPair]/I");
    mTree->Branch("clusterSize2_pxl2", &clusterSize2_pxl2, "clusterSize2_pxl2[nPair]/I");
    mTree->Branch("idTruth", &idTruth, "idTruth[nPair]/I");
    mTree->Branch("rcHftHit1_pxl1", &rcHftHit1_pxl1, "rcHftHit1_pxl1[nPair]/I");
    mTree->Branch("rcHftHit1_pxl2", &rcHftHit1_pxl2, "rcHftHit1_pxl2[nPair]/I");
    mTree->Branch("rcHftHit1_ist",  &rcHftHit1_ist,  "rcHftHit1_ist[nPair]/I");
    mTree->Branch("truth1_pxl1", &truth1_pxl1, "truth1_pxl1[nPair]/I");
    mTree->Branch("truth1_pxl2", &truth1_pxl2, "truth1_pxl2[nPair]/I");
    mTree->Branch("truth1_ist",  &truth1_ist,  "truth1_ist[nPair]/I");
    mTree->Branch("nHits1_pxl1", &nHits1_pxl1, "nHits1_pxl1[nPair]/I");
    mTree->Branch("nHits1_pxl2", &nHits1_pxl2, "nHits1_pxl2[nPair]/I");
    mTree->Branch("nHits1_ist",  &nHits1_ist,  "nHits1_ist[nPair]/I");
    mTree->Branch("rcHftHit2_pxl1", &rcHftHit2_pxl1, "rcHftHit2_pxl1[nPair]/I");
    mTree->Branch("rcHftHit2_pxl2", &rcHftHit2_pxl2, "rcHftHit2_pxl2[nPair]/I");
    mTree->Branch("rcHftHit2_ist",  &rcHftHit2_ist,  "rcHftHit2_ist[nPair]/I");
    mTree->Branch("truth2_pxl1", &truth2_pxl1, "truth2_pxl1[nPair]/I");
    mTree->Branch("truth2_pxl2", &truth2_pxl2, "truth2_pxl2[nPair]/I");
    mTree->Branch("truth2_ist",  &truth2_ist,  "truth2_ist[nPair]/I");
    mTree->Branch("nHits2_pxl1", &nHits2_pxl1, "nHits2_pxl1[nPair]/I");
    mTree->Branch("nHits2_pxl2", &nHits2_pxl2, "nHits2_pxl2[nPair]/I");
    mTree->Branch("nHits2_ist",  &nHits2_ist,  "nHits2_ist[nPair]/I");

    cout << "StMcAnalysisMaker::Init - DONE" << endl;
    return StMaker::Init();
}

//__________________________________
int StMcAnalysisMaker::Make()
{
    cout << "StMcAnalysisMaker::Make() - call" << endl;
    StMcEvent* mcEvent = (StMcEvent*)GetDataSet("StMcEvent");
    
    if (!mcEvent)
    {
        LOG_WARN << "No StMcEvent" << endm;
        return kStWarn;
    }
    
    StEvent* event = (StEvent*)GetDataSet("StEvent");
    if (!event)
    {
        LOG_WARN << "No StEvent" << endm;
        return kStWarn;
    }
    
    cout << "StMcAnalysisMaker::Make() : event: " << event->id() << endl;
    return fillTracks(mcEvent,event);
}
//____________________________________
int StMcAnalysisMaker::fillTracks(StMcEvent* mcEvent,StEvent* event)
{
    
    
    StSPtrVecMcTrack& trks = mcEvent->tracks();
    cout << "Filling " << trks.size() << " mcTracks..." << endl;
    
    if (!mAssoc->rcTpcHitMap())
    {
        cout << "!!!There is no rcTpcHitMap in the association maker!!!!" << endl;
        return 1;
    }
    
    
    nPair=0;
    nMcPxl1Hits=0;
    nMcPxl2Hits=0;
    nMcIstHits=0;
    nRcPxl1Hits=0;
    nRcPxl2Hits=0;
    nRcIstHits=0;
    
    for (unsigned int i = 0;  i < trks.size(); i++){
        StMcTrack* mcTrack = trks[i];
        Int_t trackGid = mcTrack->geantId();
        if (trackGid==1 && mcTrack->stopVertex()) {
            StMcTrack * positron = 0;
            StMcTrack * electron = 0;
            StTrack * rcPositron = 0;
            StTrack * rcElectron = 0;
            for (int j = 0; j < mcTrack->stopVertex()->numberOfDaughters(); j++){
                StMcTrack * dauTrack = mcTrack->stopVertex()->daughter(j);
                Int_t dauTrackGid = dauTrack->geantId();
                
                if (dauTrackGid!=2 && dauTrackGid!=3 ) continue;
                
                int ncommonhits = 0;
                StTrack const* rcTrack = findPartner(dauTrack, ncommonhits);
                if(rcTrack)
                {
                    if (dauTrackGid==2) {
                        positron = dauTrack;
                        rcPositron = (StTrack *)rcTrack;
                    }
                    else if (dauTrackGid==3) {
                        electron = dauTrack;
                        rcElectron = (StTrack *)rcTrack;
                    }
                }
            }// end dauther loop
            if (rcElectron && rcPositron && electron && positron) {
                pairPt[nPair] = mcTrack->pt();
                pairEta[nPair] = mcTrack->pseudoRapidity();
                openangle[nPair] = 0;
                mcopenangle[nPair] = positron->momentum().angle(electron->momentum());
                mcDist_pxl1[nPair] = 0;
                mcDist_pxl2[nPair] = 0;
                mcDist_ist[nPair] = 0;
                rcDist_pxl1[nPair] = 0;
                rcDist_pxl2[nPair] = 0;
                rcDist_ist[nPair] = 0;
                convR[nPair] = TMath::Sqrt(electron->startVertex()->position().x()*electron->startVertex()->position().x()+electron->startVertex()->position().y()*electron->startVertex()->position().y());
                parentGid[nPair] = mcTrack->geantId();
                mass[nPair] = 0;
                pairDca[nPair] = 0;
                pt1[nPair] = positron->pt();
                pt2[nPair] = electron->pt();
                eta1[nPair] = positron->pseudoRapidity();
                eta2[nPair] = electron->pseudoRapidity();
                clusterSize1_pxl1[nPair] = 0;
                clusterSize1_pxl2[nPair] = 0;
                clusterSize2_pxl1[nPair] = 0;
                clusterSize2_pxl2[nPair] = 0;
                idTruth[nPair] = 1; // electron(2) or positron(1) idTruth
                rcHftHit1_pxl1[nPair]=0;
                rcHftHit1_pxl2[nPair]=0;
                rcHftHit1_ist[nPair]=0;
                truth1_pxl1[nPair]=0;
                truth1_pxl2[nPair]=0;
                truth1_ist[nPair]=0;
                nHits1_pxl1[nPair]=0;
                nHits1_pxl2[nPair]=0;
                nHits1_ist[nPair]=0;
                rcHftHit2_pxl1[nPair]=0;
                rcHftHit2_pxl2[nPair]=0;
                rcHftHit2_ist[nPair]=0;
                truth2_pxl1[nPair]=0;
                truth2_pxl2[nPair]=0;
                truth2_ist[nPair]=0;
                nHits2_pxl1[nPair]=0;
                nHits2_pxl2[nPair]=0;
                nHits2_ist[nPair]=0;
                
                
                if (1) nMcPxl1Hits++;
                if (1) nMcPxl2Hits++;
                if (1) nMcIstHits++;
                if (1) nRcPxl1Hits++;
                if (1) nRcPxl2Hits++;
                if (1) nRcIstHits++;
                nPair++;
            }
        }
    }
    mTree->Fill();
    
    return kStOk;
}
//________________________________________________
const StTrack* StMcAnalysisMaker::findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits)
{
    //..StMcTrack find partner from the StTracks
    pair<mcTrackMapIter, mcTrackMapIter> p = mAssoc->mcTrackMap()->equal_range(mcTrack);
    
    const StTrack* maxTrack = 0;
    maxCommonTpcHits = 0;
    for (mcTrackMapIter k = p.first; k != p.second; ++k)
    {
        int commonTpcHits = k->second->commonTpcHits();
        const StTrack* track = k->second->partnerTrack()->node()->track(global);//should be global
        if (track && commonTpcHits > maxCommonTpcHits)
        {
            maxTrack = track;
            maxCommonTpcHits = commonTpcHits;
        }
    }
    return maxTrack;
}


//________________________________________________
const StMcTrack* StMcAnalysisMaker::findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits)
{
    //.. StGlobalTracks find partner from StMcTracks.
    //.. See example from StRoot/StMcAnalysisMaker
    pair<rcTrackMapIter, rcTrackMapIter> p = mAssoc->rcTrackMap()->equal_range(rcTrack);
    
    const StMcTrack* maxTrack = 0;
    maxCommonTpcHits = 0;
    for (rcTrackMapIter k = p.first; k != p.second; ++k)
    {
        int commonTpcHits = k->second->commonTpcHits();
        
        const StMcTrack* track = k->second->partnerMcTrack();
        
        if (track && commonTpcHits > maxCommonTpcHits)
        {
            maxTrack = track;
            maxCommonTpcHits = commonTpcHits;
        }
    }
    return maxTrack;
}
//______________________________________________________
int StMcAnalysisMaker::Finish()
{
    mFile->cd();
    mTree->Write();
    mFile->Close();
    return kStOk;
}
