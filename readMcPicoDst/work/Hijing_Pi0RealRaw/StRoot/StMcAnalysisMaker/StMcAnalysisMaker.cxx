#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TSystem.h"
#include "TVector3.h"

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

#include "StEvent/StPxlHit.h"
#include "StEvent/StPxlHitCollection.h"
#include "StPxlRawHitMaker/StPxlRawHit.h"
#include "StPxlRawHitMaker/StPxlRawHitMaker.h"
#include "StPxlRawHitMaker/StPxlRawHitCollection.h"
#include "StLorentzVectorF.hh"


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
    mTree->Branch("nTrack",&nTrack,"nTrack/I");
    mTree->Branch("trksGeantId",&trksGeantId,"trksGeantId[nTrack]/I");
    mTree->Branch("trksParentGeantId",&trksParentGeantId,"trksParentGeantId[nTrack]/I");
    mTree->Branch("trksGeantProcess",&trksGeantProcess,"trksGeantProcess[nTrack]/I");
    mTree->Branch("trksGeantMedium",&trksGeantMedium,"trksGeantMedium[nTrack]/I");
    mTree->Branch("trksGeneratorProcess",&trksGeneratorProcess,"trksGeneratorProcess[nTrack]/I");
    mTree->Branch("trksNumberOfDaughters",&trksNumberOfDaughters,"trksNumberOfDaughters[nTrack]/I");
    mTree->Branch("trksPt",&trksPt,"trksPt[nTrack]/F");
    mTree->Branch("trksConvR",&trksConvR,"trksConvR[nTrack]/F");

    mTree->Branch("nPair",&nPair,"nPair/I");
    mTree->Branch("nMcPxl1Hits",&nMcPxl1Hits,"nMcPxl1Hits/I");
    mTree->Branch("nMcPxl2Hits",&nMcPxl2Hits,"nMcPxl2Hits/I");
    mTree->Branch("nMcIstHits", &nMcIstHits, "nMcIstHits/I");
    mTree->Branch("nRcPxl1Hits",&nRcPxl1Hits,"nRcPxl1Hits/I");
    mTree->Branch("nRcPxl2Hits",&nRcPxl2Hits,"nRcPxl2Hits/I");
    mTree->Branch("nRcIstHits", &nRcIstHits, "nRcIstHits/I");
    mTree->Branch("nRcPxl1HitsCheck",&nRcPxl1HitsCheck,"nRcPxl1HitsCheck/I");
    mTree->Branch("nRcPxl2HitsCheck",&nRcPxl2HitsCheck,"nRcPxl2HitsCheck/I");
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
    mTree->Branch("mcConvR", &mcConvR, "mcConvR[nPair]/F");
    mTree->Branch("rcConvR", &rcConvR, "rcConvR[nPair]/F");
    mTree->Branch("parentGid", &parentGid, "parentGid[nPair]/I");
    mTree->Branch("mass", &mass, "mass[nPair]/F");
    mTree->Branch("pairDca", &pairDca, "pairDca[nPair]/F");
    mTree->Branch("pt1", &pt1, "pt1[nPair]/F");
    mTree->Branch("pt2", &pt2, "pt2[nPair]/F");
    mTree->Branch("eta1", &eta1, "eta1[nPair]/F");
    mTree->Branch("eta2", &eta2, "eta2[nPair]/F");
    mTree->Branch("clusterSize1_pxl1", &clusterSize1_pxl1, "clusterSize1_pxl1[nPair]/I");
    mTree->Branch("clusterSize1_pxl2", &clusterSize1_pxl2, "clusterSize1_pxl2[nPair]/I");
    mTree->Branch("clusterSize1_pxl3", &clusterSize1_pxl3, "clusterSize1_pxl3[nPair]/I");
    mTree->Branch("clusterSize2_pxl1", &clusterSize2_pxl1, "clusterSize2_pxl1[nPair]/I");
    mTree->Branch("clusterSize2_pxl2", &clusterSize2_pxl2, "clusterSize2_pxl2[nPair]/I");
    mTree->Branch("clusterSize2_pxl3", &clusterSize2_pxl3, "clusterSize2_pxl3[nPair]/I");
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
    
    mTree->Branch("nHits",&nHits,"nHits/I");
    mTree->Branch("clusterSize_pxl1", &clusterSize_pxl1, "clusterSize_pxl1[nHits]/I");
    mTree->Branch("clusterSize_pxl2", &clusterSize_pxl2, "clusterSize_pxl2[nHits]/I");
    mTree->Branch("mcPt_pxl1", &mcPt_pxl1, "mcPt_pxl1[nHits]/F");
    mTree->Branch("mcPt1_pxl1", &mcPt1_pxl1, "mcPt1_pxl1[nHits]/F");
    mTree->Branch("mcPt2_pxl1", &mcPt2_pxl1, "mcPt2_pxl1[nHits]/F");
    mTree->Branch("hitGeantId", &hitGeantId, "hitGeantId[nHits]/I");

    
    hGeantId = new TH1F("hGeantId","hGeantId",100,0,100);
    
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
    
    float bField = event->runInfo()->magneticField();
    
    cout << "Init Tree..." << endl;
    initTree();
    
    StPxlHitCollection * pxlHitCol = event->pxlHitCollection();
    StMcPxlHitCollection * pxlMcHitCol = mcEvent->pxlHitCollection();
    
    // track loop
    cout << "Track Loop start" << endl;
    for (unsigned int i = 0;  i < trks.size(); i++){
        StMcTrack* mcTrack = trks[i];
        Int_t trackGid = mcTrack->geantId();
        
        trksGeantId[nTrack] = trackGid;
        if (mcTrack->parent()) trksParentGeantId[nTrack] = mcTrack->parent()->geantId();
        trksGeantProcess[nTrack] = mcTrack->startVertex()->geantProcess();
        trksGeantMedium[nTrack] = mcTrack->startVertex()->geantMedium();
        trksGeneratorProcess[nTrack] = mcTrack->startVertex()->generatorProcess();
        trksNumberOfDaughters[nTrack] = mcTrack->startVertex()->numberOfDaughters();
        trksPt[nTrack] = mcTrack->pt();
        trksConvR[nTrack] = TMath::Sqrt(mcTrack->startVertex()->position().x()*mcTrack->startVertex()->position().x()+mcTrack->startVertex()->position().y()*mcTrack->startVertex()->position().y());
        nTrack++;
        
        hGeantId->Fill(trackGid);
        
        if ((trackGid == 1 || trackGid==7 || trackGid==10007) &&
            mcTrack->stopVertex())
        {
            StMcTrack * positron = 0;
            StMcTrack * electron = 0;
            StTrack * rcPositron = 0;
            StTrack * rcElectron = 0;
            // find electron and positron from gamma or pi0dalitz
            for (unsigned int j = 0; j < mcTrack->stopVertex()->numberOfDaughters(); j++){
                StMcTrack * dauTrack = mcTrack->stopVertex()->daughter(j);
                Int_t dauTrackGid = 0;
                if (dauTrack) dauTrackGid = dauTrack->geantId();
                
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
            
            if (rcElectron && rcPositron && electron && positron)
            if (StGlobalTrack const* glRcPositron = dynamic_cast<StGlobalTrack const*>(rcPositron))
            if (StGlobalTrack const* glRcElectron = dynamic_cast<StGlobalTrack const*>(rcElectron)) {
                StPhysicalHelixD electronHelix = glRcElectron->dcaGeometry()->helix();
                StPhysicalHelixD partnerHelix = glRcPositron->dcaGeometry()->helix();
                
                // normal method
                pair<double,double> ss = electronHelix.pathLengths(partnerHelix);
                StThreeVectorD kAtDcaToPartner = electronHelix.at(ss.first);
                StThreeVectorD pAtDcaToElectron = partnerHelix.at(ss.second);
                
                // calculate DCA of partner to electron at their DCA
                StThreeVectorD VectorDca = kAtDcaToPartner - pAtDcaToElectron;
                
                // calculate Lorentz vector of electron-partner pair
                StThreeVectorF const electronMomAtDca = electronHelix.momentumAt(ss.first, bField * kilogauss);
                StThreeVectorF const partnerMomAtDca = partnerHelix.momentumAt(ss.second, bField * kilogauss);
                
                StLorentzVectorF const electronFourMom(electronMomAtDca, electronMomAtDca.massHypothesis(0.000511));
                StLorentzVectorF const partnerFourMom(partnerMomAtDca, partnerMomAtDca.massHypothesis(0.000511));
                StLorentzVectorF const epairFourMom = electronFourMom + partnerFourMom;
                
                StThreeVectorF const epairMomAtDca = epairFourMom.vect();
                StThreeVectorF const Position = (kAtDcaToPartner + pAtDcaToElectron)/2.0;
                
                float mPhiV, mOpenAngle;
                phiCalculation(partnerFourMom, electronFourMom, bField > 0 ? 1 : -1, mPhiV, mOpenAngle);
                
                // Pxl
                StPtrVecHit PartnerPxlHits1 = glRcPositron->detectorInfo()->hits(kPxlId);
                StPtrVecHit PartnerPxlHits2 = glRcElectron->detectorInfo()->hits(kPxlId);

                int nPartnerPxlHits1 = (int) PartnerPxlHits1.size();
                int nPartnerPxlHits2 = (int) PartnerPxlHits2.size();
                
                uint pxl1Truth1 = 1; // first Pxl positron
                uint pxl1Truth2 = 1; // first Pxl electron
                int pxl1Hits1 = 0;
                int pxl1Hits2 = 0;
                
                uint pxl2Truth1 = 1;  // second Pxl positron
                uint pxl2Truth2 = 1;  // second Pxl electron
                int pxl2Hits1 = 0;
                int pxl2Hits2 = 0;
                
                StThreeVectorF pxl1HitPosition1 = 0;
                StThreeVectorF pxl1HitPosition2 = 0;
                clusterSize1_pxl1[nPair] = 0;
                clusterSize1_pxl2[nPair] = 0;
                clusterSize1_pxl3[nPair] = 0;
                clusterSize2_pxl1[nPair] = 0;
                clusterSize2_pxl2[nPair] = 0;
                clusterSize2_pxl3[nPair] = 0;
                
                if (nPartnerPxlHits1==0 || nPartnerPxlHits2==0) continue;
                std::vector<int> idTruthList;
                
                for(int ipxlhit=0; ipxlhit<nPartnerPxlHits1; ipxlhit++) {
                    StPxlHit * pxlHit = dynamic_cast<StPxlHit *>(PartnerPxlHits1[ipxlhit]);
                    
                    StThreeVectorF pos = PartnerPxlHits1[ipxlhit]->position();
                    float const R = pow(pos.x(),2.0) + pow(pos.y(),2.0);
                    idTruthList.push_back(pxlHit->idTruth());
                    if(R < 3.5*3.5) {
                        pxl1HitPosition1 = pos;
                        clusterSize1_pxl1[nPair] = pxlHit->nRawHits();
                        nRcPxl1HitsCheck++;
                    }
                    else if (clusterSize1_pxl2[nPair]) {
                        clusterSize1_pxl3[nPair] = pxlHit->nRawHits();
                        nRcPxl2HitsCheck++;
                    }
                    else {
                        clusterSize1_pxl2[nPair] = pxlHit->nRawHits();
                        nRcPxl2HitsCheck++;
                    }
                    
                    // set truth
                    if(PartnerPxlHits1[ipxlhit]->idTruth() == positron->key()) continue;
                    if(R > 3.5*3.5){
                        pxl2Truth1 = 0;
                    }
                    else{
                        pxl1Truth1 = 0;
                    }
                }
                for(int ipxlhit=0; ipxlhit<nPartnerPxlHits2; ipxlhit++) {
                    StPxlHit * pxlHit = dynamic_cast<StPxlHit *>(PartnerPxlHits2[ipxlhit]);
                    StThreeVectorF pos = PartnerPxlHits2[ipxlhit]->position();
                    float const R = pow(pos.x(),2.0) + pow(pos.y(),2.0);
                    bool doubleCount = false;
                    for (int k=0; k<idTruthList.size(); k++) {
                        if (idTruthList.at(k)==pxlHit->idTruth()) {
                            doubleCount=true;
                        }
                    }
                    if (doubleCount == false) {
                        if(R < 3.5*3.5) {
                            pxl1HitPosition2 = pos;
                            clusterSize2_pxl1[nPair] = pxlHit->nRawHits();
                            nRcPxl1HitsCheck++;
                        }
                        else if (clusterSize2_pxl2[nPair]) {
                            clusterSize2_pxl3[nPair] = pxlHit->nRawHits();
                            nRcPxl2HitsCheck++;
                        }
                        else {
                            clusterSize2_pxl2[nPair] = pxlHit->nRawHits();
                            nRcPxl2HitsCheck++;
                        }
                        idTruth[nPair] = 2;
                    }
                    else idTruth[nPair] = 1;

                    if(PartnerPxlHits2[ipxlhit]->idTruth() == electron->key()) continue;
                    if(R > 3.5*3.5){
                        pxl2Truth2 = 0;
                    }
                    else{
                        pxl1Truth2 = 0;
                    }
                }
                
                const StPtrVecMcPxlHit mcPxlHits1 = positron->pxlHits();
                const StPtrVecMcPxlHit mcPxlHits2 = electron->pxlHits();
                StThreeVectorF mcPxl1HitPosition1 = 0;
                StThreeVectorF mcPxl1HitPosition2 = 0;
                //Loop over PXL hits to separate into layers
                for(int ipxlhit = 0; ipxlhit < (int)mcPxlHits1.size(); ipxlhit++){
                    if((int)mcPxlHits1.at(ipxlhit)->ladder() > 1){
                        pxl2Hits1++;
                    }
                    else{
                        pxl1Hits1++;
                        mcPxl1HitPosition1 = mcPxlHits1.at(ipxlhit)->position();
                    }
                }
                for(int ipxlhit = 0; ipxlhit < (int)mcPxlHits2.size(); ipxlhit++){
                    if((int)mcPxlHits2.at(ipxlhit)->ladder() > 1){
                        pxl2Hits2++;
                    }
                    else{
                        pxl1Hits2++;
                        mcPxl1HitPosition2= mcPxlHits2.at(ipxlhit)->position();
                        
                    }
                }
                
                Int_t hftHitMap1 = (Int_t)((UInt_t)(glRcPositron->topologyMap().data(0)) >> 1 & 0x7F);
                Int_t hftHitMap2 = (Int_t)((UInt_t)(glRcElectron->topologyMap().data(0)) >> 1 & 0x7F);
                
                pairPt[nPair] = mcTrack->pt();
                pairEta[nPair] = mcTrack->pseudoRapidity();
                openangle[nPair] = mOpenAngle;
                mcopenangle[nPair] = positron->momentum().angle(electron->momentum());
                mcDist_pxl1[nPair] = (mcPxl1HitPosition1-mcPxl1HitPosition2).mag();
                mcDist_pxl2[nPair] = 0;
                rcDist_pxl1[nPair] = (pxl1HitPosition1-pxl1HitPosition2).mag();
                rcDist_pxl2[nPair] = 0;
                mcConvR[nPair] = TMath::Sqrt(mcTrack->stopVertex()->position().x()*mcTrack->stopVertex()->position().x()+mcTrack->stopVertex()->position().y()*mcTrack->stopVertex()->position().y());
                rcConvR[nPair] = TMath::Sqrt(Position.x()*Position.x()+Position.y()+Position.y());
                parentGid[nPair] = mcTrack->geantId();
                mass[nPair] = epairFourMom.m();
                pairDca[nPair] = static_cast<float>(VectorDca.mag());
                pt1[nPair] = positron->pt();
                pt2[nPair] = electron->pt();
                eta1[nPair] = positron->pseudoRapidity();
                eta2[nPair] = electron->pseudoRapidity();
                //clusterSize1_pxl1[nPair] = 0;
                //clusterSize1_pxl2[nPair] = 0;
                //clusterSize1_pxl3[nPair] = 0;
                //clusterSize2_pxl1[nPair] = 0;
                //clusterSize2_pxl2[nPair] = 0;
                //clusterSize2_pxl3[nPair] = 0;
                //idTruth[nPair] = 1; // electron(2) or positron(1) idTruth
                rcHftHit1_pxl1[nPair]=hftHitMap1>>0 & 0x1;
                rcHftHit1_pxl2[nPair]=hftHitMap1>>1 & 0x3;
                rcHftHit1_ist[nPair]=hftHitMap1>>3 & 0x3;
                truth1_pxl1[nPair]=pxl1Truth1;
                truth1_pxl2[nPair]=pxl1Truth2;
                truth1_ist[nPair]=0;
                nHits1_pxl1[nPair]=pxl1Hits1;
                nHits1_pxl2[nPair]=pxl1Hits2;
                nHits1_ist[nPair]=0;
                rcHftHit2_pxl1[nPair]=hftHitMap2>>0 & 0x1;
                rcHftHit2_pxl2[nPair]=hftHitMap2>>1 & 0x3;
                rcHftHit2_ist[nPair]=hftHitMap2>>3 & 0x3;
                truth2_pxl1[nPair]=pxl2Truth1;
                truth2_pxl2[nPair]=pxl2Truth2;
                truth2_ist[nPair]=0;
                nHits2_pxl1[nPair]=pxl2Hits1;
                nHits2_pxl2[nPair]=pxl2Hits2;
                nHits2_ist[nPair]=0;
                nPair++;
                
            }
        }
    }
    
    // pxl hit loop
    cout << "Rc hit Loop start" << endl;
    nHits = 0;
    for (unsigned int iSec = 0; iSec<pxlHitCol->numberOfSectors(); iSec++){
        StPxlSectorHitCollection * pxlSecHitCol = pxlHitCol->sector(iSec);
        if (!pxlSecHitCol) continue;
        for (unsigned int iLad = 0; iLad < pxlSecHitCol->numberOfLadders(); iLad++) {
            StPxlLadderHitCollection * pxlLadHitCol = pxlSecHitCol->ladder(iLad);
            if (!pxlLadHitCol) continue;
            for (unsigned int iSen = 0; iSen<pxlLadHitCol->numberOfSensors(); iSen++) {
                StPxlSensorHitCollection * pxlSenHitCol = pxlLadHitCol->sensor(iSen);
                if (!pxlSenHitCol) continue;
                UInt_t nSenHits = pxlSenHitCol->hits().size();
                for (unsigned int iHit = 0; iHit < nSenHits; iHit++){
                    StPxlHit* pixHit = pxlSenHitCol->hits()[iHit];
                    if (!pixHit) continue;
                    StMcTrack * hitTrack = 0;
                    StMcTrack * hitTrack1 = 0;
                    StMcTrack * hitTrack2 = 0;
                    for (int i=0;i<trks.size();i++){
                        if (pixHit->idTruth()==trks[i]->key()) {
                            hitTrack = trks[i];
                            hitTrack1 = hitTrack->startVertex()->daughter(0)->geantId() == 2 ? hitTrack->startVertex()->daughter(0) : hitTrack->startVertex()->daughter(1);
                            hitTrack2 = hitTrack->startVertex()->daughter(0)->geantId() == 3 ? hitTrack->startVertex()->daughter(0) : hitTrack->startVertex()->daughter(1);
                            break;
                        }
                    }
                    hitGeantId[nHits] = hitTrack->geantId();
                    if (pixHit->ladder() == 1){
                        nRcPxl1Hits++;
                        clusterSize_pxl1[nHits] = pixHit->nRawHits();
                        mcPt_pxl1[nHits] = hitTrack->pt();
                        if (hitTrack1) mcPt1_pxl1[nHits] = hitTrack1->pt();
                        if (hitTrack2) mcPt2_pxl1[nHits] = hitTrack2->pt();
                    }
                    else {nRcPxl2Hits++;
                        clusterSize_pxl2[nHits] = pixHit->nRawHits();
                    }
                    nHits++;
                }
            }
        }
    }
    cout << "nRcPxl1Hits / nRcPxl2Hits / pxlHitCol->numberOfHits() : " << nRcPxl1Hits << " / " << nRcPxl2Hits << " / " << pxlHitCol->numberOfHits() << endl;
    
    cout << "Mc hit Loop start" << endl;
    for (unsigned int iSec = 0; iSec<pxlMcHitCol->numberOfSectors(); iSec++){
        StMcPxlSectorHitCollection * pxlSecHitCol = pxlMcHitCol->sector(iSec);
        if (!pxlSecHitCol) continue;
        for (unsigned int iLad = 0; iLad < pxlSecHitCol->numberOfLadders(); iLad++) {
            StMcPxlLadderHitCollection * pxlLadHitCol = pxlSecHitCol->ladder(iLad);
            if (!pxlLadHitCol) continue;
            for (unsigned int iSen = 0; iSen<pxlLadHitCol->numberOfSensors(); iSen++) {
                StMcPxlSensorHitCollection * pxlSenHitCol = pxlLadHitCol->sensor(iSen);
                if (!pxlSenHitCol) continue;
                UInt_t nSenHits = pxlSenHitCol->hits().size();
                for (unsigned int iHit = 0; iHit < nSenHits; iHit++){
                    StMcPxlHit* pixHit = pxlSenHitCol->hits()[iHit];
                    if (!pixHit) continue;
                    if (pixHit->ladder() == 1)nMcPxl1Hits++;
                    else nMcPxl2Hits++;
                }
            }
        }
    }
    cout << "nMcPxl1Hits / nMcPxl2Hits / pxlMcHitCol->numberOfHits() : " << nMcPxl1Hits << " / " << nMcPxl2Hits << " / " << pxlMcHitCol->numberOfHits() << endl;

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
    hGeantId->Write();
    
    mFile->Close();
    return kStOk;
}
//______________________________________________________
void StMcAnalysisMaker::phiCalculation(StLorentzVectorF const positron,StLorentzVectorF const electron, int mN, float &phiV, float &openangle)
{
    TVector3 ppp(positron.px(),positron.py(),positron.pz());
    TVector3 eee(electron.px(),electron.py(),electron.pz());
    TVector3 u=ppp+eee;
    TVector3 v=eee.Cross(ppp);
    TVector3 w=u.Cross(v);
    TVector3 nz(0.,0.,mN);
    TVector3 wc=u.Cross(nz);
    
    phiV =w.Angle(wc);
    openangle=ppp.Angle(eee);
    
}
void StMcAnalysisMaker::initTree(){
    
    nPair=0;
    nMcPxl1Hits=0;
    nMcPxl2Hits=0;
    nMcIstHits=0;
    nRcPxl1Hits=0;
    nRcPxl2Hits=0;
    nRcIstHits=0;
    nRcPxl1HitsCheck=0;
    nRcPxl2HitsCheck=0;
    for (int i=0; i<kMaxPair; i++) {
        pairPt[i]=-999;
        pairEta[i]=-999;
        openangle[i]=-999;
        mcopenangle[i]=-999;
        mcDist_pxl1[i]=-999;
        mcDist_pxl2[i]=-999;
        mcDist_ist[i]=-999;
        rcDist_pxl1[i]=-999;
        rcDist_pxl2[i]=-999;
        rcDist_ist[i]=-999;
        mcConvR[i]=-999;
        rcConvR[i]=-999;
        parentGid[i]=-999;
        mass[i]=-999;
        pairDca[i]=-999;
        pt1[i]=-999;
        pt2[i]=-999;
        eta1[i]=-999;
        eta2[i]=-999;
        clusterSize1_pxl1[i]=-999;
        clusterSize1_pxl2[i]=-999;
        clusterSize1_pxl3[i]=-999;
        clusterSize2_pxl1[i]=-999;
        clusterSize2_pxl2[i]=-999;
        clusterSize2_pxl3[i]=-999;
        idTruth[i]=-999;
        rcHftHit1_pxl1[i]=-999;
        rcHftHit1_pxl2[i]=-999;
        rcHftHit1_ist[i]=-999;
        truth1_pxl1[i]=-999;
        truth1_pxl2[i]=-999;
        truth1_ist[i]=-999;
        nHits1_pxl1[i]=-999;
        nHits1_pxl2[i]=-999;
        nHits1_ist[i]=-999;
        rcHftHit2_pxl1[i]=-999;
        rcHftHit2_pxl2[i]=-999;
        rcHftHit2_ist[i]=-999;
        truth2_pxl1[i]=-999;
        truth2_pxl2[i]=-999;
        truth2_ist[i]=-999;
        nHits2_pxl1[i]=-999;
        nHits2_pxl2[i]=-999;
        nHits2_ist[i]=-999;
    }
    nHits=0;
    for (int i=0; i<kMaxHits; i++) {
        clusterSize_pxl1[i]=-999;
        clusterSize_pxl2[i]=-999;
        mcPt_pxl1[i]=-999;
        mcPt1_pxl1[i]=-999;
        mcPt2_pxl1[i]=-999;
        hitGeantId[i]=-999;
        
    }
    nTrack=0;
    for (int i=0; i<kMaxTrack; i++) {
        trksGeantId[i]=-999;
        trksParentGeantId[i]=-999;
        trksGeantProcess[i]=-999;
        trksGeantMedium[i]=-999;
        trksGeneratorProcess[i]=-999;
        trksNumberOfDaughters[i]=-999;
        trksPt[i]=-999;
        trksConvR[i]=-999;
    }


}
