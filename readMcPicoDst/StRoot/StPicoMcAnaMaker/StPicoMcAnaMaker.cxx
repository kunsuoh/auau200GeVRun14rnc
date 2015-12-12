#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StPicoNpeEventMaker/StPicoNpeEvent.h"
#include "StPicoNpeEventMaker/StElectronPair.h"

#include "StLorentzVectorF.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"

#include "StPicoMcAnaMaker.h"
#include "StNpeCuts.h"


ClassImp(StPicoMcAnaMaker)

StPicoMcAnaMaker::StPicoMcAnaMaker(char const * name, char const * outName, StPicoDstMaker* picoDstMaker):
StMaker(name),mPicoDstMaker(picoDstMaker),mOutFileName(outName),
mOutputFile(NULL), mChain(NULL), mEventCounter(0), mNpeCuts(NULL)
{

}

Int_t StPicoMcAnaMaker::Init()
{

    mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
    mOutputFile->cd();
    
    if (!mNpeCuts)
        mNpeCuts = new StNpeCuts;
    mNpeCuts->init();
    
    
    // -------------- USER VARIABLES -------------------------
    
    return kStOK;
}
//-----------------------------------------------------------------------------
StPicoMcAnaMaker::~StPicoMcAnaMaker()
{
    /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoMcAnaMaker::Finish()
{
    LOG_INFO << " StPicoMcAnaMaker - writing data and closing output file " <<endm;
    mOutputFile->cd();
    // --------------- USER HISTOGRAM WRITE --------------------
    
    
    mOutputFile->Close();
    
    return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoMcAnaMaker::Make()
{
    readNextEvent();
    
    if (!mPicoDstMaker)
    {
        LOG_WARN << " StPicoMcAnaMaker - No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    
    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    if (!picoDst)
    {
        LOG_WARN << "StPicoMcAnaMaker - No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    
    // -------------- USER ANALYSIS -------------------------
    
    
    
    return kStOK;
}
//-----------------------------------------------------------------------------
void StPicoMcAnaMaker::phiCalculation(StLorentzVectorF positron,StLorentzVectorF electron, double mN, double &phiV, double &openangle)
{
    StThreeVector<double> ppp(positron.px(),positron.py(),positron.pz());
    StThreeVector<double> eee(electron.px(),electron.py(),electron.pz());
    StThreeVector<double> u=ppp+eee;
    StThreeVector<double> v=eee.cross(ppp);
    StThreeVector<double> w=u.cross(v);
    StThreeVector<double> nz(0.,0.,mN);
    StThreeVector<double> wc=u.cross(nz);
    
    phiV =w.angle(wc);
    openangle=ppp.angle(eee);
    
}