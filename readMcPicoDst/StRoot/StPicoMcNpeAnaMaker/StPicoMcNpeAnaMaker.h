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
    void  fillHistogram(StPicoTrack const*) const;
    void  fillHistogram(StElectronPair const*) const;
    bool  isRcTrack(StPicoMcTrack const * const trk ,StPicoDst const * const PicoDst, int &id);
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoEvent*    mPicoEvent;
    
    TFile* mOutputFile;
    TTree* mTree;
    
    TH1F * hEventVz;
    TH1F * hTrackParentGeantId;
    TH1F * hTrackGeantId;
    TH1F * hPairMass;
    TH1F * hPairDca;
    TH2F * hPairPosition;

    
    ClassDef(StPicoMcNpeAnaMaker, 0)
};

#endif
