#include <vector>
#include <cmath>

#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "../StPicoDstMaker/StPicoDst.h"
#include "../StPicoDstMaker/StPicoDstMaker.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoDstMaker/StPicoTrack.h"
#include "../StPicoDstMaker/StPicoMcEvent.h"
#include "../StPicoDstMaker/StPicoMcTrack.h"

#include "StPicoMcNpeAnaMaker.h"
#include "StElectronPair.h"
#include "StCuts.h"


ClassImp(StPicoMcNpeAnaMaker)

//-----------------------------------------------------------------------------
StPicoMcNpeAnaMaker::StPicoMcNpeAnaMaker(char const* makerName, StPicoDstMaker* picoMaker, char const* fileBaseName)
: StMaker(makerName), mPicoDstMaker(picoMaker), mPicoEvent(NULL),
mOutputFile(NULL)
{
    
    TString baseName(fileBaseName);
    mOutputFile = new TFile(Form("%s.hist.root",fileBaseName), "RECREATE");
}

//-----------------------------------------------------------------------------
StPicoMcNpeAnaMaker::~StPicoMcNpeAnaMaker()
{
    /* mTree is owned by mOutputFile directory, it will be destructed once
     * the file is closed in ::Finish() */
}

//-----------------------------------------------------------------------------
Int_t StPicoMcNpeAnaMaker::Init()
{
    
    hEventVz = new TH1F("hEventVz","hEventVz",200,-10,10);
    hTrackParentGeantId = new TH1F("hTrackParentGeantId","hTrackParentGeantId",100,0,100);
    hTrackGeantId = new TH1F("hTrackGeantId","hTrackGeantId",100,0,100);
    hPairMass = new TH1F("hPairMass","hPairMass",100,-0.1,0.4);
    hPairDca = new TH1F("hPairDca","hPairDca",100,0,0.5);
    hPairPosition = new TH2F("hPairPosition","hPairPosition",200,-10,10,200,-10,10);
    
    tree = new TTree("T","Electron pair tree");
    tree->Branch("rc_x",&rc_x,"rc_x/F");
    tree->Branch("rc_y",&rc_y,"rc_y/F");
    tree->Branch("rc_z",&rc_z,"rc_z/F");
    
    tree->Branch("mc_x",&mc_x,"mc_x/F");
    tree->Branch("mc_y",&mc_y,"mc_y/F");
    tree->Branch("mc_z",&mc_z,"mc_z/F");
    
    tree->Branch("pt1",&pt1,"pt1/F");
    tree->Branch("pt2",&pt2,"pt2/F");
    tree->Branch("phiV",&phiV,"phiV/F");
    tree->Branch("openangle",&openangle,"openangle/F");
    tree->Branch("mcopenangle",&mcopenangle,"mcopenangle/F");
    tree->Branch("phi",&phi,"phi/F");
    tree->Branch("eta",&eta,"eta/F");
    tree->Branch("mass",&mass,"mass/F");
    tree->Branch("pairDca",&pairDca,"pairDca/F");
    tree->Branch("mcPairPt",&mcPairPt,"mcPairPt/F");
    tree->Branch("angle",&angle,"angle/F");
    tree->Branch("length",&length,"length/F");          //
    tree->Branch("parentGid",&parentGid,"parentGid/I");   //
    tree->Branch("refmult",&refmult,"refmult/I");   //
    tree->Branch("chi1",&chi1,"chi1/F");
    tree->Branch("chi2",&chi2,"chi2/F");
    tree->Branch("rcdca1",&rcdca1,"rcdca1/F");
    tree->Branch("rcdca2",&rcdca2,"rcdca2/F");
    tree->Branch("mcdca1",&mcdca1,"mcdca1/F");
    tree->Branch("mcdca2",&mcdca2,"mcdca2/F");
    
    tree->Branch("distHits",&distHits,"distHits/F");
    
    tree->Branch("nHits1_pxl1",&nHits1_pxl1,"nHits1_pxl1/b");   //
    tree->Branch("nHits1_pxl2",&nHits1_pxl2,"nHits1_pxl2/b");   //
    tree->Branch("nHits1_ist",&nHits1_ist,"nHits1_ist/b");   //
    tree->Branch("nHits1_ssd",&nHits1_ssd,"nHits1_ssd/b");   //
    
    tree->Branch("nHits2_pxl1",&nHits2_pxl1,"nHits2_pxl1/b");   //
    tree->Branch("nHits2_pxl2",&nHits2_pxl2,"nHits2_pxl2/b");   //
    tree->Branch("nHits2_ist",&nHits2_ist,"nHits2_ist/b");   //
    tree->Branch("nHits2_ssd",&nHits2_ssd,"nHits2_ssd/b");   //
    
    tree->Branch("truth1_pxl1",&truth1_pxl1,"truth1_pxl1/b");   //
    tree->Branch("truth1_pxl2",&truth1_pxl2,"truth1_pxl2/b");   //
    tree->Branch("truth1_ist",&truth1_ist,"truth1_ist/b");   //
    tree->Branch("truth1_ssd",&truth1_ssd,"truth1_ssd/b");   //
    
    tree->Branch("truth2_pxl1",&truth2_pxl1,"truth2_pxl1/b");   //
    tree->Branch("truth2_pxl2",&truth2_pxl2,"truth2_pxl2/b");   //
    tree->Branch("truth2_ist",&truth2_ist,"truth2_ist/b");   //
    tree->Branch("truth2_ssd",&truth2_ssd,"truth2_ssd/b");   //
    
    tree->Branch("rcHftHit1_pxl1",&rcHftHit1_pxl1,"rcHftHit1_pxl1/b");   //
    tree->Branch("rcHftHit1_pxl2",&rcHftHit1_pxl2,"rcHftHit1_pxl2/b");   //
    tree->Branch("rcHftHit1_ist",&rcHftHit1_ist,"rcHftHit1_ist/b");   //
    tree->Branch("rcHftHit1_ssd",&rcHftHit1_ssd,"rcHftHit1_ssd/b");   //
    
    tree->Branch("rcHftHit2_pxl1",&rcHftHit2_pxl1,"rcHftHit2_pxl1/b");   //
    tree->Branch("rcHftHit2_pxl2",&rcHftHit2_pxl2,"rcHftHit2_pxl2/b");   //
    tree->Branch("rcHftHit2_ist",&rcHftHit2_ist,"rcHftHit2_ist/b");   //
    tree->Branch("rcHftHit2_ssd",&rcHftHit2_ssd,"rcHftHit2_ssd/b");   //
    
    
    singleTree = new TTree("eT","Single Electron tree");
    singleTree->Branch("mc_x",&mc_x,"mc_x/F");
    singleTree->Branch("mc_y",&mc_y,"mc_y/F");
    singleTree->Branch("mc_z",&mc_z,"mc_z/F");
    singleTree->Branch("rcPt",&rcPt,"rcPt/F");
    singleTree->Branch("rcPhi",&rcPhi,"rcPhi/F");
    singleTree->Branch("rcEta",&rcEta,"rcEta/F");
    singleTree->Branch("mcPt",&mcPt,"mcPt/F");
    singleTree->Branch("mcPhi",&mcPhi,"mcPhi/F");
    singleTree->Branch("mcEta",&mcEta,"mcEta/F");
    singleTree->Branch("parentGid",&parentGid,"parentGid/I");   //
    singleTree->Branch("parentGid2",&parentGid2,"parentGid2/I");   //
    singleTree->Branch("chi",&chi,"chi/F");
    
    singleTree->Branch("nHits_pxl1",&nHits_pxl1,"nHits_pxl1/b");   //
    singleTree->Branch("nHits_pxl2",&nHits_pxl2,"nHits_pxl2/b");   //
    singleTree->Branch("nHits_ist",&nHits_ist,"nHits_ist/b");   //
    singleTree->Branch("nHits_ssd",&nHits_ssd,"nHits_ssd/b");   //

    singleTree->Branch("truth_pxl1",&truth_pxl1,"truth_pxl1/b");   //
    singleTree->Branch("truth_pxl2",&truth_pxl2,"truth_pxl2/b");   //
    singleTree->Branch("truth_ist",&truth_ist,"truth_ist/b");   //
    singleTree->Branch("truth_ssd",&truth_ssd,"truth_ssd/b");   //
    
    singleTree->Branch("rcHftHit_pxl1",&rcHftHit_pxl1,"rcHftHit_pxl1/b");   //
    singleTree->Branch("rcHftHit_pxl2",&rcHftHit_pxl2,"rcHftHit_pxl2/b");   //
    singleTree->Branch("rcHftHit_ist",&rcHftHit_ist,"rcHftHit_ist/b");   //
    singleTree->Branch("rcHftHit_ssd",&rcHftHit_ssd,"rcHftHit_ssd/b");   //
    
    singleTree->Branch("rcdca",&rcdca,"rcdca/F");
    singleTree->Branch("mcdca",&mcdca,"mcdca/F");
    
    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StPicoMcNpeAnaMaker::Finish()
{
    mOutputFile->cd();
    // write histograms
    hEventVz->Write();
    hTrackParentGeantId->Write();
    hTrackGeantId->Write();
    hPairMass->Write();
    hPairDca->Write();
    hPairPosition->Write();
    
    tree->Write();
    singleTree->Write();
    
    mOutputFile->Write();
    mOutputFile->Close();
    return kStOK;
}
//-----------------------------------------------------------------------------
void StPicoMcNpeAnaMaker::Clear(Option_t *opt)
{
}

//-----------------------------------------------------------------------------
Int_t StPicoMcNpeAnaMaker::Make()
{

    if (!mPicoDstMaker)
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    
    StPicoDst const * picoDst = mPicoDstMaker->picoDst();
    if (!picoDst)
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    
    mPicoEvent = picoDst->event();
    //cout << "before isGoodEvent()" << endl;
    if (isGoodEvent())
    {
        //cout << "in isGoodEvent()" << endl;
        std::vector<Int_t> idPicoDstRcElectrons;
        std::vector<Int_t> idPicoDstRcPositrons;
        std::vector<Int_t> idPicoDstMcElectrons;
        std::vector<Int_t> idPicoDstMcPositrons;

        StThreeVectorF pVtx = mPicoEvent->primaryVertex();
        float const bField = mPicoEvent->bField();
        int nMcTracks =  picoDst->numberOfMcTracks();
        refmult = mPicoEvent->refMult();
        
        //cout << "start mcTrack loop" << endl;
        for(int i_Mc=0; i_Mc<nMcTracks; i_Mc++){
            //cout << i_Mc << " " ;
            // get Mc Track
            StPicoMcTrack *mcTrk = (StPicoMcTrack*)picoDst->mctrack(i_Mc);
          //  if(mcTrk->Pxl1Truth()==0 || mcTrk->Pxl2Truth()==0) continue;

            
            // get Geant Id for track and parent
            parentGid= -999;
            if(mcTrk->parentId() != Pico::USHORTMAX) {
                StPicoMcTrack *mcParentTrk = (StPicoMcTrack*)picoDst->mctrack(mcTrk->parentId());
                parentGid=mcParentTrk->GePid();
                //if(mcParentTrk->parentId() != Pico::USHORTMAX) continue;
                delete mcParentTrk;
            }
            trackId=mcTrk->GePid();
            //cout << parentGid << " " << trackId << " " ;
            
            //if (parentGid!=1) continue;
            
            hTrackParentGeantId->Fill(parentGid);
            hTrackGeantId->Fill(trackId);
            
            // get Rc Trcak
            //if (parentGid != cuts::parentGid && cuts::parentGid != Pico::USHORTMAX) continue;
            if (mcTrk->parentId() != Pico::USHORTMAX) continue;
            if (trackId != cuts::dau1Gid && trackId != cuts::dau2Gid) continue;
            
            StPicoTrack *rcTrk=0;
            Int_t id=-999;
            isRcTrack(mcTrk,picoDst,id);
            //cout << id << " " ;
            if(id!=-999){
                
                rcTrk = (StPicoTrack*)picoDst->track(id);
                fillHistogram(rcTrk,mcTrk);
                
                nHits_pxl1 = (mcTrk->hitsPxl1());
                nHits_pxl2 = (mcTrk->hitsPxl2());
                nHits_ist = (mcTrk->hitsIst());
                nHits_ssd = (mcTrk->hitsSst());
                
                truth_pxl1 = (mcTrk->Pxl1Truth());
                truth_pxl2 = (mcTrk->Pxl2Truth());
                truth_ist = (mcTrk->IstTruth());
                truth_ssd = (mcTrk->SsdTruth());
                
                rcHftHit_pxl1 = rcTrk->nHitsMapHFT()>>0 & 0x1;
                rcHftHit_pxl2 = rcTrk->nHitsMapHFT()>>1 & 0x3;
                rcHftHit_ist = rcTrk->nHitsMapHFT()>>3 & 0x3;
                rcHftHit_ssd = 0;
                
                mc_x = mcTrk->Origin().x();
                mc_y = mcTrk->Origin().y();
                mc_z = mcTrk->Origin().z();
                
                rcPt = rcTrk->gPt();
                rcPhi = rcTrk->gMom(pVtx,bField).phi();
                rcEta = rcTrk->gMom(pVtx,bField).pseudoRapidity();
                mcPt = mcTrk->Mom().perp();
                mcPhi = mcTrk->Mom().phi();
                mcEta = mcTrk->pseudorapidity();
                chi = rcTrk->chi2();
                
                StPhysicalHelixD rcHelix = rcTrk->dcaGeometry().helix();
                rcdca = rcHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                
                StPhysicalHelixD mcHelix(mcTrk->Mom(), mcTrk->Origin(), bField, trackId == 2 ? 1 : -1);
                mcdca = mcHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                
                parentGid2 = ((StPicoMcTrack*)(picoDst->mctrack(mcTrk->parentId())))->GePid();
                
                singleTree->Fill();
                
                if (trackId==cuts::dau1Gid) {
                    idPicoDstRcPositrons.push_back(id);
                    idPicoDstMcPositrons.push_back(i_Mc);
                }
                else if (trackId==cuts::dau2Gid){
                    idPicoDstRcElectrons.push_back(id);
                    idPicoDstMcElectrons.push_back(i_Mc);
                }
                
            }
  //          cout << endl;
        }
        //cout << idPicoDstRcPositrons.size() << " " << idPicoDstRcElectrons.size() << endl;
        
        for (int i=0; i<idPicoDstMcPositrons.size(); i++) {
            StPicoMcTrack *mcPositron = (StPicoMcTrack*)picoDst->mctrack(idPicoDstMcPositrons[i]);
            for (int j=0; j<idPicoDstMcElectrons.size(); j++) {
                StPicoMcTrack *mcElectron = (StPicoMcTrack*)picoDst->mctrack(idPicoDstMcElectrons[j]);
                
                if (mcPositron->parentId() == Pico::USHORTMAX) continue;
                if (mcPositron->parentId() != mcElectron->parentId()) continue;
                //     if (mcElectron->hitsPxl1()==0) continue;
                //     if (mcElectron->hitsPxl2()==0) continue;
                //     if (mcPositron->hitsPxl1()==0) continue;
                //     if (mcPositron->hitsPxl2()==0) continue;
                StPicoTrack *rcPositron = picoDst->track(idPicoDstRcPositrons[i]);
                StPicoTrack *rcElectron = picoDst->track(idPicoDstRcElectrons[j]);
                StElectronPair * rcPair = new StElectronPair(rcPositron,rcElectron,i,j,bField,pVtx);
                
                if (isGoodTrack(rcPositron) && isGoodTrack(rcElectron)) {
                    fillHistogram(rcPair);
                    
                    length = rcPair->length();
                    angle = rcPair->angle();
                    mcPairPt = ((StPicoMcTrack*)(picoDst->mctrack(mcElectron->parentId())))->Mom().perp();
                    pairDca = rcPair->pairDca();
                    mass = rcPair->pairMass();
                    eta = rcPair->eta();
                    phi = rcPair->phi();
                    openangle = rcPair->openAngle();
                    TVector3 ppp(mcPositron->Mom().x(),mcPositron->Mom().y(),mcPositron->Mom().z());
                    TVector3 eee(mcElectron->Mom().x(),mcElectron->Mom().y(),mcElectron->Mom().z());
                    mcopenangle = ppp.Angle(eee);
                    phiV = rcPair->phiV();
                    pt1 = rcPositron->gPt();
                    pt2 = rcElectron->gPt() * -1;
                    
                    rc_x = rcPair->positionX();
                    rc_y = rcPair->positionY();
                    rc_z = rcPair->positionZ();
                    
                    mc_x = mcElectron->Origin().x();
                    mc_y = mcElectron->Origin().y();
                    mc_z = mcElectron->Origin().z();
                    
                    nHits1_pxl1 = (mcPositron->hitsPxl1());
                    nHits1_pxl2 = (mcPositron->hitsPxl2());
                    nHits1_ist = (mcPositron->hitsIst());
                    nHits1_ssd = (mcPositron->hitsSst());
                    
                    truth1_pxl1 = (mcPositron->Pxl1Truth());
                    truth1_pxl2 = (mcPositron->Pxl2Truth());
                    truth1_ist = (mcPositron->IstTruth());
                    truth1_ssd = (mcPositron->SsdTruth());
                    
                    nHits2_pxl1 = (mcElectron->hitsPxl1());
                    nHits2_pxl2 = (mcElectron->hitsPxl2());
                    nHits2_ist = (mcElectron->hitsIst());
                    nHits2_ssd = (mcElectron->hitsSst());
                    
                    truth2_pxl1 = (mcElectron->Pxl1Truth());
                    truth2_pxl2 = (mcElectron->Pxl2Truth());
                    truth2_ist = (mcElectron->IstTruth());
                    truth2_ssd = (mcElectron->SsdTruth());
                    
                    // get Geant Id for track and parent
                    parentGid= -999;
                    if(mcElectron->parentId() != Pico::USHORTMAX) {
                        StPicoMcTrack *mcParentTrk = (StPicoMcTrack*)picoDst->mctrack(mcElectron->parentId());
                        parentGid=mcParentTrk->GePid();
                        delete mcParentTrk;
                    }
                    
                    chi1 = rcPositron->chi2();
                    chi2 = rcElectron->chi2();
                    
                    rcHftHit1_pxl1 = rcPositron->nHitsMapHFT()>>0 & 0x1;
                    rcHftHit1_pxl2 = rcPositron->nHitsMapHFT()>>1 & 0x3;
                    rcHftHit1_ist = rcPositron->nHitsMapHFT()>>3 & 0x3;
                    rcHftHit1_ssd = 0;
                    
                    rcHftHit2_pxl1 = rcElectron->nHitsMapHFT()>>0 & 0x1;
                    rcHftHit2_pxl2 = rcElectron->nHitsMapHFT()>>1 & 0x3;
                    rcHftHit2_ist = rcElectron->nHitsMapHFT()>>3 & 0x3;
                    rcHftHit2_ssd = 0;
                    
                    StPhysicalHelixD rc1Helix = rcPositron->dcaGeometry().helix();
                    StPhysicalHelixD rc2Helix = rcElectron->dcaGeometry().helix();
                    rcdca1 = rc1Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                    rcdca2 = rc2Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                    
                    StPhysicalHelixD mc1Helix(mcPositron->Mom(), pVtx, bField, 1);
                    StPhysicalHelixD mc2Helix(mcElectron->Mom(), pVtx, bField, -1);
                    mcdca1 = mc1Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                    mcdca2 = mc2Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                    
                    distHits = sqrt(
                                    (rc1Helix.at(2.7).x() - rc2Helix.at(2.7).x()) * (rc1Helix.at(2.7).x() - rc2Helix.at(2.7).x()) +
                                    (rc1Helix.at(2.7).y() - rc2Helix.at(2.7).y()) * (rc1Helix.at(2.7).y() - rc2Helix.at(2.7).y())
                                    );
                    
                    tree->Fill();
                }
            }
        }
        
        
    } //.. end of good event
    
    return kStOK;
}

//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isGoodEvent()
{
    float vZ = mPicoEvent->primaryVertex().z();
    
    hEventVz->Fill(vZ);
    
    return fabs(vZ) < cuts::vz ;
    
}
//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
    return
    trk->gPt() > cuts::ptMin &&
    trk->gPt() < cuts::ptMax &&
    trk->nHitsFit() >= cuts::nHitsFit;// &&
//    trk->isHFTTrack();
//    isHftTrack(trk) ;
}
//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isHftTrack(StPicoTrack const * const trk) const
{
    return trk->nHitsMapHFT()>>1 & 0x3;
   // return !(trk->nHitsMapHFT()>>0 & 0x1) && (trk->nHitsMapHFT()>>1 & 0x3) && (trk->nHitsMapHFT()>>3 & 0x3);
}

//-----------------------------------------------------------------------------
void StPicoMcNpeAnaMaker::fillHistogram(StPicoTrack const * const rcTrk, StPicoMcTrack const * const mcTrk) const
{
//    hTrackPt->Fill(rcTrk->gPt());
//    hTrackNHitsFit->Fill(rcTrk->nHitsFit());
}
//-----------------------------------------------------------------------------
void StPicoMcNpeAnaMaker::fillHistogram(StElectronPair const * const pair) const
{
    
    hPairMass->Fill(pair->pairMass());
    hPairDca->Fill(pair->pairDca());
    hPairPosition->Fill(pair->positionX(),pair->positionY());

}
//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isRcTrack(StPicoMcTrack const * const PicoMcTrack, StPicoDst const * const  PicoDst,int &id)
{
    int nMcTracks =  PicoDst->numberOfMcTracks();
    if(PicoMcTrack->assoId() == Pico::USHORTMAX )
        return false;
    int temp = Pico::USHORTMAX ;
    for(int i_Rc =0; i_Rc<PicoDst->numberOfTracks(); ++i_Rc){
        StPicoTrack *Trk = (StPicoTrack*)PicoDst->track(i_Rc);
        if(PicoMcTrack->assoId() == Trk->id() ) {
            temp = i_Rc;
            break;
        }
    }
    if (temp == Pico::USHORTMAX) return false;
    id=temp;
    return true;
}
