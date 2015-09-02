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
 *
 * **************************************************
 */
 

#include "TChain.h"
#include "StMaker.h"
#include "TH2F.h"
#include "StThreeVectorF.hh"
#include "THnSparse.h"

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

int const nnpt = 6;
int const nnpid = 8;
int const nntype = 8;//+8;
int const nnhisto = 17;
int const nnhistophe = 8;
int const nntrigger = 4;



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

    
  private:
    StPicoNpeAnaMaker() {}
    void readNextEvent();

    StPicoDstMaker* mPicoDstMaker;
    StPicoNpeEvent* mPicoNpeEvent;
    StPicoPrescales* mPrescales;


    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    bool isGoodEvent() const;
    bool isGoodTrack(StPicoTrack const*) const;
    bool isGoodTagged(StPicoTrack const*) const;
    bool isGoodPartner(StPicoTrack const*) const;
    bool isGoodTofTrack(StPicoTrack const*) const;
    bool isGoodEmcTrack(StPicoTrack const*) const;
    bool isGoodPair(StElectronPair const*) const;
    void initVariables();
    void setVariables(StPicoTrack *);
    void setVariables(StElectronPair *);
    void fillHistogram(TString);
    void fillHistogram(int);
    void fillHistogram(int,int,int);
    void fillHistogram(int,int,int,int);
    void setHistogram();
    int getPtBin(double);
    bool isBHTevent();
    bool isBemc();
    bool isBemc2();
    bool isBemc3();
    bool isBemc4();
    bool isBemc5();
    bool isBsmd();
    bool isTof();
    
    
    StPicoDst * picoDst;
    
    TTree * tInc;
    TTree * tPhE;
    
    float weight[nntrigger];
    float bField;
    StThreeVectorF pVtx;
    
    Float_t pairAngle3d;
    Float_t pairAnglePhi;
    Float_t pairAngleTheta;
    Float_t pairMass;
    Float_t pairMassP;
    char pairCharge;
    Float_t pairDca;
    Float_t pairPositionX;
    Float_t pairPositionY;
    Float_t pairPositionZ;
    Float_t pairConvRadious;
    
    
    Float_t dca;
    Float_t dcaCharge;
    Float_t pt;
    Float_t partner_pt;
    Float_t eta;
    Float_t nsige;
    Float_t nsigpion;
    Float_t partner_nsige;
    Float_t beta;
    Float_t tofmass;
    Float_t e;
    Float_t e0;
    Float_t e1;
    Float_t e2;
    Float_t e3;
    Float_t eoverp;
    Float_t adc0;
    
    unsigned char neta;
    unsigned char nphi;
    unsigned char nphieta;
    Float_t phiDist;
    Float_t zDist;
    Float_t etaTowDist;
    Float_t phiTowDist;
    
    UShort_t mZDCx;
    UShort_t mRefMult;
    unsigned char isHTEvents;

    
    TH1F * hEvent;
    TH1I * hTrigger;
    TH1I * hTriggerWt;
    TH1I * hCheckDoubleTrigger;

    TH1I * hTriggerCheck[nntrigger];
    TH1I * hTriggerCheckWt[nntrigger];
    TH1F * hZDCx[nntrigger];
    TH1F * hZDCxWt[nntrigger];
    TH1F * hRefMult[nntrigger];
    TH1F * hRefMultWt[nntrigger];
    TH1F * histo[nnpt][nnpid][nntype][nnhisto][nntrigger];
    TH1F * histoPureE[nnpt][nnhistophe];
    TH2F * histo2d[nnpid][nntrigger];
    TH2F * histo2dDcaPt[nnpid][nntrigger];
    int nptbin;
    int npid;
    int ntype;
    int nhisto;
    
    
    ClassDef(StPicoNpeAnaMaker, 2)
};


inline int StPicoNpeAnaMaker::getEntries() const
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoNpeAnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

#endif
