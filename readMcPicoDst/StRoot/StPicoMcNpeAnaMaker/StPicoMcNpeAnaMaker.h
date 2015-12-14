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
#include "TH2F.h"

class TTree;
class TFile;
class StPicoDstMaker;
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
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoEvent*    mPicoEvent;
    
    TFile* mOutputFile;
    TTree* mTree;
    
    TH1F * hEventVz;
    TH1F * hTrackParentGeantId;
    TH1F * hTrackGeantId;

    
    ClassDef(StPicoMcNpeAnaMaker, 0)
};

#endif
