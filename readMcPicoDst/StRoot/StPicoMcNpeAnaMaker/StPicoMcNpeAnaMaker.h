#ifndef StPicoMcNpeAnaMaker_h
#define StPicoMcNpeAnaMaker_h

/* **************************************************
 *  A Maker that reads StPicoEvents and creates
 *  StLowPtNpeAnas and stores them.
 *
 *  Authors:  **Kunsu OH        (kunsuoh@gmail.com)
 *
 *  **Code Maintainer
 * **************************************************
 */
#include "StMaker.h"
#include "TNtuple.h"
#include "TH2F.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StEvent/StDcaGeometry.h"

#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoMcTrack.h"

class TTree;
class TFile;
class TH1F;
class TH2F;
class TTree;

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;
class TString;
class StDcaGeometry;

class StElectronPair;

class StPicoMcNpeAnaMaker : public StMaker
{
public:
    StPicoMcNpeAnaMaker(char const* makerName, StPicoDstMaker* picoMaker, char const* fileBaseName);
    virtual ~StPicoMcNpeAnaMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    
private:
    bool  isGoodEvent();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isElectron(StPicoTrack const*) const;
    bool  isTaggedElectron(StPicoTrack const*) const;
    bool  isPartnerElectron(StPicoTrack const*) const;
    bool  isGoodElectronPair(StElectronPair const &, float) const;
    void  fillHistogram(StPicoTrack const*,StPicoMcTrack const*) const;
    void  fillHistogram(StElectronPair const*) const;
    bool  isRcTrack(StPicoMcTrack const * const trk ,StPicoDst const * const PicoDst, int &id);
    bool  isHftTrack(StPicoTrack const * const trk) const;

    StPicoDstMaker* mPicoDstMaker;
    StPicoEvent*    mPicoEvent;
    TTree *tree;
    TTree *singleTree;
    TFile* mOutputFile;
    TNtuple * nt2;
    TNtuple * nt3;
    
    TH1F * hEventVz;
    TH1F * hTrackParentGeantId;
    TH1F * hTrackGeantId;
    TH1F * hPairMass;
    TH1F * hPairDca;
    TH2F * hPairPosition;

    Float_t distHits;
    
    Float_t rc_x, rc_y, rc_z, mc_x, mc_y, mc_z;
    
    UChar_t nHits1_pxl1, nHits1_pxl2, nHits1_ist, nHits1_ssd;
    UChar_t nHits2_pxl1, nHits2_pxl2, nHits2_ist, nHits2_ssd;
    UChar_t nHits_pxl1, nHits_pxl2, nHits_ist, nHits_ssd;
    
    UChar_t truth1_pxl1, truth1_pxl2, truth1_ist, truth1_ssd;
    UChar_t truth2_pxl1, truth2_pxl2, truth2_ist, truth2_ssd;
    UChar_t truth_pxl1, truth_pxl2, truth_ist, truth_ssd;
    
    UChar_t rcHftHit1_pxl1, rcHftHit1_pxl2, rcHftHit1_ist, rcHftHit1_ssd;
    UChar_t rcHftHit2_pxl1, rcHftHit2_pxl2, rcHftHit2_ist, rcHftHit2_ssd;
    UChar_t rcHftHit_pxl1, rcHftHit_pxl2, rcHftHit_ist, rcHftHit_ssd;

    
    float length, angle, mcPairPt, pairDca, mass, eta, phi, openangle, mcopenangle, phiV, pt1, pt2, chi1, chi2;
    UShort_t parentGid, parentGid2, trackId;
    float rcdca1, rcdca2, rcdca, mcdca1, mcdca2, mcdca;
    
    int refmult;
    
    float rcPt, rcPhi, rcEta, mcPt, mcPhi, mcEta, chi;
    
    ClassDef(StPicoMcNpeAnaMaker, 0)
};

#endif
