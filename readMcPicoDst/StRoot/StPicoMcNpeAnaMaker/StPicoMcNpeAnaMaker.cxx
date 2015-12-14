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
    
    cout << "check before isGoodEvent()" << endl;
    if (isGoodEvent())
    {
        
        float const bField = mPicoEvent->bField();
        int nMcTracks =  picoDst->numberOfMcTracks();
        for(int i_Mc=0; i_Mc<nMcTracks; i_Mc++){
            StPicoMcTrack *mcTrk = (StPicoMcTrack*)picoDst->mctrack(i_Mc);
            float parentId=((StPicoMcTrack*)(picoDst->mctrack(mcTrk->parentId())))->GePid();
            float gepid=mcTrk->GePid();
            hParentGeantId->Fill(parentId);
            hGeantId->Fill(gepid);
            

        }
        
        
    } //.. end of good event fill
    
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
    (float)trk->nHitsFit()/(float)trk->nHitsMax() > cuts::nHitsRatioMin &&
    (float)trk->nHitsFit()/(float)trk->nHitsMax() < cuts::nHitsRatioMax ;
}

//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isElectron(StPicoTrack const * const trk) const
{
    return
    isGoodTrack(trk) &&
     trk->nHitsDedx() >= cuts::nHitsDedx ;
}

//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isTaggedElectron(StPicoTrack const * const trk) const
{
    return
    isElectron(trk) ;
}

//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isPartnerElectron(StPicoTrack const * const trk) const
{
    return
    isGoodTrack(trk)  ;
}
//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isGoodElectronPair(StElectronPair const & epair, float pt) const
{
    return
    epair.pairMass() < cuts::pairMass &&
    epair.pairDca() < cuts::pairDca;
}
//-----------------------------------------------------------------------------
void StPicoMcNpeAnaMaker::fillHistogram(StPicoTrack const * const trk) const
{

}