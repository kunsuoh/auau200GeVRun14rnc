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
    
    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StPicoMcNpeAnaMaker::Finish()
{
    mOutputFile->cd();
    // write histograms
    
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
        
        UInt_t nTracks = picoDst->numberOfTracks();
        
        std::vector<unsigned short> idxPicoTaggedEs;
        std::vector<unsigned short> idxPicoPartnerEs;
        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
        {
            StPicoTrack const* trk = picoDst->track(iTrack);
            if (!trk) continue;
            if (isElectron(trk))
            {
                fillHistogram(trk);
                if (isTaggedElectron(trk)) idxPicoTaggedEs.push_back(iTrack);
            }
            
            if (isPartnerElectron(trk)) idxPicoPartnerEs.push_back(iTrack);
        } // .. end tracks loop
        
        
        float const bField = mPicoEvent->bField();
        /*
        for (unsigned short ik = 0; ik < idxPicoTaggedEs.size(); ++ik)
        {
            
            StPicoTrack const * electron = picoDst->track(idxPicoTaggedEs[ik]);
            
            // make electron pairs
            for (unsigned short ip = 0; ip < idxPicoPartnerEs.size(); ++ip)
            {
                
                if (idxPicoTaggedEs[ik] == idxPicoPartnerEs[ip]) continue;
                
                StPicoTrack const * partner = picoDst->track(idxPicoPartnerEs[ip]);
                
                StElectronPair electronPair(electron, partner, idxPicoTaggedEs[ik], idxPicoPartnerEs[ip], bField);
                
                if (!isGoodElectronPair(electronPair, electron->gPt())) continue;
                
            } // .. end make electron pairs
        } // .. end of tagged e loop
        */
        idxPicoTaggedEs.clear();
        idxPicoPartnerEs.clear();
        
    } //.. end of good event fill
    
    
    return kStOK;
}

//-----------------------------------------------------------------------------
bool StPicoMcNpeAnaMaker::isGoodEvent()
{
    if (fabs(mPicoEvent->primaryVertex().z()) < cuts::vz) {
        
        if (fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz) {
            return true;
        }
    }
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