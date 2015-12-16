#include <vector>
#include <cmath>

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

        float const bField = mPicoEvent->bField();
        int nMcTracks =  picoDst->numberOfMcTracks();

        for(int i_Mc=0; i_Mc<nMcTracks; i_Mc++){
            // get Mc Track
            StPicoMcTrack *mcTrk = (StPicoMcTrack*)picoDst->mctrack(i_Mc);
            
            // get Geant Id for track and parent
            float parentGid= Pico::USHORTMAX;
            if(mcTrk->parentId() != Pico::USHORTMAX)
                parentGid=((StPicoMcTrack*)(picoDst->mctrack(mcTrk->parentId())))->GePid();
            float trackId=mcTrk->GePid();

            hTrackParentGeantId->Fill(picoDst->mctrack(mcTrk->parentId())->parentId()->GePid());
            if ( picoDst->mctrack(mcTrk->parentId())->parentId()!=Pico::USHORTMAX) continue;
            hTrackGeantId->Fill(trackId);
            
            // get Rc Trcak
            if (parentGid == cuts::parentGid && (trackId == cuts::dau1Gid || trackId == cuts::dau2Gid)) {
                StPicoTrack *rcTrk=0;
                Int_t id=-999;
                isRcTrack(mcTrk,picoDst,id);
                if(id!=-999){
                    rcTrk = (StPicoTrack*)picoDst->track(id);
                    fillHistogram(rcTrk,mcTrk);
                    
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
                    StPicoTrack *rcPositron = picoDst->track(idPicoDstRcPositrons[i]);
                    StPicoTrack *rcElectron = picoDst->track(idPicoDstRcElectrons[j]);
                    StElectronPair * rcPair = new StElectronPair(rcPositron,rcElectron,i,j,bField);
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
    trk->nHitsFit() >= cuts::nHitsFit &&
    trk->isHFTTrack() ;
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
