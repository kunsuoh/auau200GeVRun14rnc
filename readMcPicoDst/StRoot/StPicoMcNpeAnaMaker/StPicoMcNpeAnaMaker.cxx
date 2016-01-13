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
    nt2 = new TNtuple("nt2","electron pair ntuple","pt1:pt2:phiV:openangle:v0x:v0y:v0z:phi:eta:mass:pairDca:mcv0x:mcv0y:mcv0z:mcPairPt");
    nt3 = new TNtuple("nt3","electron pair ntuple 3","pt1:pt2:v0x:v0y:v0z:phi:eta:mass:pairDca:mcv0x:mcv0y:mcv0z:mcPairPt:angle:length");
    
    tree = new TTree("T","Electron pair tree");
    tree->Branch("rc",&rc,"x:y:z");
    tree->Branch("mc",&mc,"x:y:z");
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
    tree->Branch("nHits1",&nHits1,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    tree->Branch("truth1",&truth1,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    tree->Branch("nHits2",&nHits2,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    tree->Branch("truth2",&truth2,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    tree->Branch("parentGid",&parentGid,"parentGid/s");   //
    tree->Branch("refmult",&refmult,"refmult/I");   //
    tree->Branch("chi1",&chi1,"chi1/F");
    tree->Branch("chi2",&chi2,"chi2/F");
    tree->Branch("rchfthit1",&rchfthit1,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    tree->Branch("rchfthit2",&rchfthit2,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    tree->Branch("rcdca1",&rcdca1,"rcdca1/F");
    tree->Branch("rcdca2",&rcdca2,"rcdca2/F");
    tree->Branch("mcdca1",&mcdca1,"mcdca1/F");
    tree->Branch("mcdca2",&mcdca2,"mcdca2/F");
    
    singleTree = new TTree("eT","Single Electron tree");
    singleTree->Branch("mc",&mc,"x:y:z");
    singleTree->Branch("rcPt",&rcPt,"rcPt/F");
    singleTree->Branch("rcPhi",&rcPhi,"rcPhi/F");
    singleTree->Branch("rcEta",&rcEta,"rcEta/F");
    singleTree->Branch("mcPt",&mcPt,"mcPt/F");
    singleTree->Branch("mcPhi",&mcPhi,"mcPhi/F");
    singleTree->Branch("mcEta",&mcEta,"mcEta/F");
    singleTree->Branch("parentGid",&parentGid,"parentGid/s");   //
    singleTree->Branch("parentGid2",&parentGid2,"parentGid2/s");   //
    singleTree->Branch("chi",&chi,"chi/F");
    singleTree->Branch("nHits",&nHits,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    singleTree->Branch("truth",&truth,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
    singleTree->Branch("rchfthit",&rchfthit,"pxl1/b:pxl2/b:ist/b:ssd/b");   //
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
    nt2->Write();
    nt3->Write();
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
    
    if (isGoodEvent())
    {
        std::vector<Int_t> idPicoDstRcElectrons;
        std::vector<Int_t> idPicoDstRcPositrons;
        std::vector<Int_t> idPicoDstMcElectrons;
        std::vector<Int_t> idPicoDstMcPositrons;

        StThreeVectorF pVtx = mPicoEvent->primaryVertex();
        float const bField = mPicoEvent->bField();
        int nMcTracks =  picoDst->numberOfMcTracks();
        refmult = mPicoEvent->refMult();

        for(int i_Mc=0; i_Mc<nMcTracks; i_Mc++){
            // get Mc Track
            StPicoMcTrack *mcTrk = (StPicoMcTrack*)picoDst->mctrack(i_Mc);
          //  if(mcTrk->Pxl1Truth()==0 || mcTrk->Pxl2Truth()==0) continue;

            
            // get Geant Id for track and parent
            Int_t parentGid= Pico::USHORTMAX;
            if(mcTrk->parentId() != Pico::USHORTMAX) {
                StPicoMcTrack *mcParentTrk = (StPicoMcTrack*)picoDst->mctrack(mcTrk->parentId());
                parentGid=mcParentTrk->GePid();
            //    if(mcParentTrk->parentId() != Pico::USHORTMAX) continue;
            }
            trackId=mcTrk->GePid();

            hTrackParentGeantId->Fill(parentGid);
            hTrackGeantId->Fill(trackId);
            
            // get Rc Trcak
            if ((parentGid == cuts::parentGid || cuts::parentGid == Pico::USHORTMAX) &&
                (trackId == cuts::dau1Gid || trackId == cuts::dau2Gid)) {
                StPicoTrack *rcTrk=0;
                Int_t id=-999;
                isRcTrack(mcTrk,picoDst,id);
                if(id!=-999){
                    rcTrk = (StPicoTrack*)picoDst->track(id);
                    fillHistogram(rcTrk,mcTrk);

                    nHits.pxl1 = (mcTrk->hitsPxl1());
                    nHits.pxl2 = (mcTrk->hitsPxl2());
                    nHits.ist = (mcTrk->hitsIst());
                    nHits.ssd = (mcTrk->hitsSst());
                    
                    truth.pxl1 = (mcTrk->Pxl1Truth());
                    truth.pxl2 = (mcTrk->Pxl2Truth());
                    truth.ist = (mcTrk->IstTruth());
                    truth.ssd = (mcTrk->SsdTruth());
                    
                    rchfthit.pxl1 = rcTrk->nHitsMapHFT()>>0 & 0x1;
                    rchfthit.pxl2 = rcTrk->nHitsMapHFT()>>1 & 0x3;
                    rchfthit.ist = rcTrk->nHitsMapHFT()>>3 & 0x3;
                    rchfthit.ssd = 0;

                    mc.x = mcTrk->Origin().x();
                    mc.y = mcTrk->Origin().y();
                    mc.z = mcTrk->Origin().z();

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
                    else {
                        idPicoDstRcElectrons.push_back(id);
                        idPicoDstMcElectrons.push_back(i_Mc);
                    }
                    
                }
            }
        }
        cout << idPicoDstRcPositrons.size() << " " << idPicoDstRcElectrons.size() << endl;
        
        for (int i=0; i<idPicoDstMcPositrons.size(); i++) {
            StPicoMcTrack *mcPositron = (StPicoMcTrack*)picoDst->mctrack(idPicoDstMcPositrons[i]);
            for (int j=0; j<idPicoDstMcElectrons.size(); j++) {
                StPicoMcTrack *mcElectron = (StPicoMcTrack*)picoDst->mctrack(idPicoDstMcElectrons[j]);
                if (mcPositron->parentId() != Pico::USHORTMAX && mcPositron->parentId() == mcElectron->parentId())
                {
               //     if (mcElectron->hitsPxl1()==0) continue;
               //     if (mcElectron->hitsPxl2()==0) continue;
               //     if (mcPositron->hitsPxl1()==0) continue;
               //     if (mcPositron->hitsPxl2()==0) continue;
                    StPicoTrack *rcPositron = picoDst->track(idPicoDstRcPositrons[i]);
                    StPicoTrack *rcElectron = picoDst->track(idPicoDstRcElectrons[j]);
                    StElectronPair * rcPair = new StElectronPair(rcPositron,rcElectron,i,j,bField,pVtx);
                    if (isGoodTrack(rcPositron) && isGoodTrack(rcElectron)) {
                        fillHistogram(rcPair);
                        nt2->Fill(rcPositron->gPt(),
                                  rcElectron->gPt()*-1,
                                  rcPair->phiV(),
                                  rcPair->openAngle(),
                                  rcPair->positionX(),rcPair->positionY(),rcPair->positionZ(),
                                  rcPair->phi(),rcPair->eta(),
                                  rcPair->pairMass(),rcPair->pairDca(),
                                  mcElectron->Origin().x(),mcElectron->Origin().y(),mcElectron->Origin().z(),
                                  ((StPicoMcTrack*)(picoDst->mctrack(mcElectron->parentId())))->Mom().perp()
                                  );
                        nt3->Fill(rcPositron->gPt(),
                                  rcElectron->gPt()*-1,
                                  rcPair->positionX(),rcPair->positionY(),rcPair->positionZ(),
                                  rcPair->phi(),rcPair->eta(),
                                  rcPair->pairMass(),rcPair->pairDca(),
                                  mcElectron->Origin().x(),mcElectron->Origin().y(),mcElectron->Origin().z(),
                                  ((StPicoMcTrack*)(picoDst->mctrack(mcElectron->parentId())))->Mom().perp(),
                                  rcPair->angle(),
                                  rcPair->length()
                                  );

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
                        
                        rc.x = rcPair->positionX();
                        rc.y = rcPair->positionY();
                        rc.z = rcPair->positionZ();
                        
                        mc.x = mcElectron->Origin().x();
                        mc.y = mcElectron->Origin().y();
                        mc.z = mcElectron->Origin().z();
                        
                        nHits1.pxl1 = (mcPositron->hitsPxl1());
                        nHits1.pxl2 = (mcPositron->hitsPxl2());
                        nHits1.ist = (mcPositron->hitsIst());
                        nHits1.ssd = (mcPositron->hitsSst());
                        
                        truth1.pxl1 = (mcPositron->Pxl1Truth());
                        truth1.pxl2 = (mcPositron->Pxl2Truth());
                        truth1.ist = (mcPositron->IstTruth());
                        truth1.ssd = (mcPositron->SsdTruth());

                        nHits2.pxl1 = (mcElectron->hitsPxl1());
                        nHits2.pxl2 = (mcElectron->hitsPxl2());
                        nHits2.ist = (mcElectron->hitsIst());
                        nHits2.ssd = (mcElectron->hitsSst());
                        
                        truth2.pxl1 = (mcElectron->Pxl1Truth());
                        truth2.pxl2 = (mcElectron->Pxl2Truth());
                        truth2.ist = (mcElectron->IstTruth());
                        truth2.ssd = (mcElectron->SsdTruth());
                        parentGid = ((StPicoMcTrack*)(picoDst->mctrack(mcElectron->parentId())))->GePid();
                        chi1 = rcPositron->chi2();
                        chi2 = rcElectron->chi2();
                        
                        rchfthit1.pxl1 = rcPositron->nHitsMapHFT()>>0 & 0x1;
                        rchfthit1.pxl2 = rcPositron->nHitsMapHFT()>>1 & 0x3;
                        rchfthit1.ist = rcPositron->nHitsMapHFT()>>3 & 0x3;
                        rchfthit1.ssd = 0;

                        rchfthit2.pxl1 = rcElectron->nHitsMapHFT()>>0 & 0x1;
                        rchfthit2.pxl2 = rcElectron->nHitsMapHFT()>>1 & 0x3;
                        rchfthit2.ist = rcElectron->nHitsMapHFT()>>3 & 0x3;
                        rchfthit2.ssd = 0;
                        
                        StPhysicalHelixD rc1Helix = rcPositron->dcaGeometry().helix();
                        StPhysicalHelixD mc1Helix(mcPositron->Mom(), mcPositron->Origin(), bField, 1);
                        StPhysicalHelixD rc2Helix = rcElectron->dcaGeometry().helix();
                        StPhysicalHelixD mc2Helix(mcElectron->Mom(), mcElectron->Origin(), bField, -1);
                        rcdca1 = rc1Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                        mcdca1 = mc1Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                        rcdca2 = rc2Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());
                        mcdca2 = mc2Helix.curvatureSignedDistance(pVtx.x(),pVtx.y());


                        tree->Fill();
                    }
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
