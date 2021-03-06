#ifndef ST_MCANALYSISMAKER_H
#define ST_MCANALYSISMAKER_H

#include "TString.h"
#include "StLorentzVectorF.hh"


class TFile;
class TH3F;
class TH2F;
class TH1F;
class TTree;

class StMcTrack;
class StTrack;
class StGlobalTrack;
class StAssociationMaker;
class StMcEvent;
class StEvent;

#include "StMaker.h"


const Int_t kMaxPair= 500;
const Int_t kMaxHits= 5000;
const Int_t kMaxTrack= 5000;
class StMcAnalysisMaker : public StMaker
{
private:
    TFile* mFile;
    TTree * mTree;
    TH1F * hGeantId;
    Int_t nPair;
    Int_t nMcPxl1Hits;
    Int_t nMcPxl2Hits;
    Int_t nMcIstHits;
    Int_t nRcPxl1Hits;
    Int_t nRcPxl2Hits;
    Int_t nRcIstHits;
    Int_t nRcPxl1HitsCheck;
    Int_t nRcPxl2HitsCheck;
    Float_t pairPt[kMaxPair];
    Float_t pairEta[kMaxPair];
    Float_t openangle[kMaxPair];
    Float_t mcopenangle[kMaxPair];
    Float_t mcDist_pxl1[kMaxPair];
    Float_t mcDist_pxl2[kMaxPair];
    Float_t mcDist_ist[kMaxPair];
    Float_t rcDist_pxl1[kMaxPair];
    Float_t rcDist_pxl2[kMaxPair];
    Float_t rcDist_ist[kMaxPair];
    Float_t mcConvR[kMaxPair];
    Float_t rcConvR[kMaxPair];
    Int_t parentGid[kMaxPair];
    Float_t mass[kMaxPair];
    Float_t pairDca[kMaxPair];
    Float_t pt1[kMaxPair];
    Float_t pt2[kMaxPair];
    Float_t eta1[kMaxPair];
    Float_t eta2[kMaxPair];
    Int_t clusterSize1_pxl1[kMaxPair];
    Int_t clusterSize1_pxl2[kMaxPair];
    Int_t clusterSize1_pxl3[kMaxPair];
    Int_t clusterSize2_pxl1[kMaxPair];
    Int_t clusterSize2_pxl2[kMaxPair];
    Int_t clusterSize2_pxl3[kMaxPair];
    Int_t idTruth[kMaxPair];
    Int_t rcHftHit1_pxl1[kMaxPair];
    Int_t rcHftHit1_pxl2[kMaxPair];
    Int_t rcHftHit1_ist[kMaxPair];
    Int_t truth1_pxl1[kMaxPair];
    Int_t truth1_pxl2[kMaxPair];
    Int_t truth1_ist[kMaxPair];
    Int_t nHits1_pxl1[kMaxPair];
    Int_t nHits1_pxl2[kMaxPair];
    Int_t nHits1_ist[kMaxPair];
    Int_t rcHftHit2_pxl1[kMaxPair];
    Int_t rcHftHit2_pxl2[kMaxPair];
    Int_t rcHftHit2_ist[kMaxPair];
    Int_t truth2_pxl1[kMaxPair];
    Int_t truth2_pxl2[kMaxPair];
    Int_t truth2_ist[kMaxPair];
    Int_t nHits2_pxl1[kMaxPair];
    Int_t nHits2_pxl2[kMaxPair];
    Int_t nHits2_ist[kMaxPair];
    
    Int_t nHits;
    Int_t clusterSize_pxl1[kMaxHits];
    Int_t clusterSize_pxl2[kMaxHits];
    Float_t mcPt_pxl1[kMaxHits];
    Float_t mcPt1_pxl1[kMaxHits];
    Float_t mcPt2_pxl1[kMaxHits];
    Int_t hitGeantId[kMaxHits];
    
    Int_t nTrack;
    Int_t trksGeantId[kMaxTrack];
    Int_t trksParentGeantId[kMaxTrack];
    Int_t trksGeantProcess[kMaxTrack];
    Int_t trksGeantMedium[kMaxTrack];
    Int_t trksGeneratorProcess[kMaxTrack];
    Int_t trksNumberOfDaughters[kMaxTrack];
    Float_t trksPt[kMaxTrack];
    Float_t trksConvR[kMaxTrack];

    StAssociationMaker* mAssoc;
    const StTrack* findPartner(StMcTrack*, int&);
    const StMcTrack* findPartner(StGlobalTrack*, int&);
    void phiCalculation(StLorentzVectorF const,StLorentzVectorF const, int, float &, float &);
    int fillTracks(StMcEvent*,StEvent*);
    void initTree();
    
    TString mOutFileName;
public:
    StMcAnalysisMaker (const char *name="StMcAnalysisMaker", const char *title="event/StMcAnalysisMaker");
    
    int Init();
    int Make();
    int Finish();
    
    virtual void setOutFileName(TString in);
    
    ClassDef(StMcAnalysisMaker, 1)
};

inline void StMcAnalysisMaker::setOutFileName(TString in) { mOutFileName = in;}
#endif
