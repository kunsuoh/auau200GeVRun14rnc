#ifndef StPicoNpeAnaMaker_h
#define StPicoNpeAnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoNpeEvent
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *          **Kunsu OH        (kunsuoh@gmail.com)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */



#include "TChain.h"
#include "StMaker.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TNtuple.h"
#include "StThreeVectorF.hh"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "StLorentzVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoNpeEvent;
class StElectronTrack;
class StElectronPair;
class StPicoDstMaker;
class StPicoTrack;
class StPicoDst;
class StPicoPrescales;
class TTree;

class StNpeCuts;

class StPicoNpeAnaMaker : public StMaker
{
  public:
    StPicoNpeAnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker);
    virtual ~StPicoNpeAnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setNpeCuts(StNpeCuts* cuts);
    void phiCalculation(StLorentzVectorF ,StLorentzVectorF , double , Float_t &, Float_t &);

    
  private:
    StPicoNpeAnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StElectronPair const*) const;
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoNpeEvent* mPicoNpeEvent;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StNpeCuts* mNpeCuts;

    
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    
    TH1I * h1dEvent;
    TH1F * h1dEventZDCx;
    TH1I * h1dEventRefMult;
    TH1I * h1dEventTrigger;
    TH1F * h1dEventZDCxCut;
    TH1I * h1dEventRefMultCut;
    TH1I * h1dEventTriggerCut;
    
    TH1I * h1dTrack;
    
    TH2F * h2dIncEDcaVsPt;
    TH2F * h2dIncENSigEVsPt;
    TH2F * h2dIncEDcaVsPtCut;
    TH2F * h2dIncENSigEVsPtCut;
    TH2F * h2dIncEDcaVsPtCut2;
    TH2F * h2dIncENSigEVsPtCut2;
    
    TH2F * h2dIncEBsmdNEtaPt;
    TH2F * h2dIncEBsmdNPhiPt;

    TH2F * h2dPhEDcaVsPt;
    TH2F * h2dPhENSigEVsPt;
    TH2F * h2dPhEConvRVsPt;
    
    TH2F * h2dPhELDcaVsPt;
    TH2F * h2dPhELNSigEVsPt;
    TH2F * h2dPhELConvRVsPt;
    
    TH1F * hQaPt;
    TH1F * hQaEta;
    TH1F * hQaDca;
    TH1F * hQaNHitFit;
    TH1F * hQaNHitDedx;
    
    TH1F * hQaPtCut;
    TH1F * hQaEtaCut;
    TH1F * hQaDcaCut;
    TH1F * hQaNHitFitCut;
    TH1F * hQaNHitDedxCut;
    
    TH2D * h2dPhENSigEVsZ;
    TH2D * h2dPhEConvRVsZ;
    TH3D * h2dPhEConvXYZ;
    TH2D * h2dPhEInvMassvsZ;
    TH2D * h2dPhELNSigEVsZ;
    TH2D * h2dPhELConvRVsZ;
    TH3D * h2dPhELConvXYZ;
    TH2D * h2dPhELInvMassvsZ;
    
    TH2D * h2dPhENSigEVsZ_HFT;
    TH2D * h2dPhEConvRVsZ_HFT;
    TH3D * h2dPhEConvXYZ_HFT;
    TH2D * h2dPhEInvMassvsZ_HFT;
    TH2D * h2dPhELNSigEVsZ_HFT;
    TH2D * h2dPhELConvRVsZ_HFT;
    TH3D * h2dPhELConvXYZ_HFT;
    TH2D * h2dPhELInvMassvsZ_HFT;

    TTree *tree;

    TH2D * h2dPhEMassVsPt;
    TH2D * h2dPhELMassVsPt;
    TH2D * h2dPhEPairDcaVsPt;
    TH2D * h2dPhELPairDcaVsPt;
    
    Float_t distHits;
    Float_t rc_x, rc_y, rc_z;
    UChar_t rcHftHit1_pxl1, rcHftHit1_pxl2, rcHftHit1_ist, rcHftHit1_ssd;
    UChar_t rcHftHit2_pxl1, rcHftHit2_pxl2, rcHftHit2_ist, rcHftHit2_ssd;
    UChar_t rcHftHit_pxl1, rcHftHit_pxl2, rcHftHit_ist, rcHftHit_ssd;
    Float_t length, angle, pairDca, mass, eta, phi, openangle, phiV, pt1, pt2, chi1, chi2, nsige1, nsige2, rcdca1, rcdca2;
    Int_t refmult;
    
    ClassDef(StPicoNpeAnaMaker, 0)
};

inline int StPicoNpeAnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoNpeAnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}
inline void StPicoNpeAnaMaker::setNpeCuts(StNpeCuts* cuts)
{
    mNpeCuts = cuts;
}


#endif
