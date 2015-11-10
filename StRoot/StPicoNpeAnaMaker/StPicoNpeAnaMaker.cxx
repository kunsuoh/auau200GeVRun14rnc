#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoNpeEventMaker/StPicoNpeEvent.h"
#include "StPicoNpeEventMaker/StElectronPair.h"

#include "StBTofUtil/tofPathLength.hh"
#include "StLorentzVectorF.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"

#include "StPicoNpeAnaMaker.h"
#include "StNpeCuts.h"


ClassImp(StPicoNpeAnaMaker)

StPicoNpeAnaMaker::StPicoNpeAnaMaker(char const * name,char const * inputFilesList,
                                     char const * outName, StPicoDstMaker* picoDstMaker):
StMaker(name),mPicoDstMaker(picoDstMaker),mPicoNpeEvent(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
mOutputFile(NULL), mChain(NULL), mEventCounter(0), mNpeCuts(NULL)
{}

Int_t StPicoNpeAnaMaker::Init()
{
    mPicoNpeEvent = new StPicoNpeEvent();
    
    mChain = new TChain("T");
    std::ifstream listOfFiles(mInputFileList.Data());
    if (listOfFiles.is_open())
    {
        std::string file;
        while (getline(listOfFiles, file))
        {
            LOG_INFO << "StPicoNpeAnaMaker - Adding :" << file << endm;
            mChain->Add(file.c_str());
        }
    }
    else
    {
        LOG_ERROR << "StPicoNpeAnaMaker - Could not open list of files. ABORT!" << endm;
        return kStErr;
    }
    mChain->GetBranch("npeEvent")->SetAutoDelete(kFALSE);
    mChain->SetBranchAddress("npeEvent", &mPicoNpeEvent);
    
    mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
    mOutputFile->cd();
    
    if (!mNpeCuts)
        mNpeCuts = new StNpeCuts;
    mNpeCuts->init();
    

    // -------------- USER VARIABLES -------------------------
    h2dIncEDcaVsPt =    new TH2F("h2dIncEDcaVsPt",      "2D Dca vs pT Inclusive E",             200,0,20, 100,-0.1,0.1);
    h2dIncENSigEVsPt =  new TH2F("h2dIncENSigEVsPt",    "2D nSigE vs pT Inclusive E",           200,0,20, 289,-13,13);
 
    h2dIncEDcaVsPtCut =    new TH2F("h2dIncEDcaVsPtCut",      "2D Dca vs pT Inclusive E Cut",             200,0,20, 100,-0.1,0.1);
    h2dIncENSigEVsPtCut =  new TH2F("h2dIncENSigEVsPtCut",    "2D nSigE vs pT Inclusive E Cut",           200,0,20, 289,-13,13);

    h2dIncEDcaVsPtCut2 =    new TH2F("h2dIncEDcaVsPtCut2",      "2D Dca vs pT Inclusive E Cut2",             200,0,20, 100,-0.1,0.1);
    h2dIncENSigEVsPtCut2 =  new TH2F("h2dIncENSigEVsPtCut2",    "2D nSigE vs pT Inclusive E Cut2",           200,0,20, 289,-13,13);

    h2dPhEDcaVsPt =     new TH2F("h2dPhEDcaVsPt",       "2D Dca vs pT Photonic E",              200,0,20, 100,-0.1,0.1);
    h2dPhENSigEVsPt =   new TH2F("h2dPhENSigEVsPt",     "2D nSigE vs pT Photonic E",            200,0,20, 289,-13,13);
    h2dPhEConvRVsPt =   new TH2F("h2dPhEConvRVsPt",     "2D Conversion Radius vs pT Photonic E",200,0,20, 100,0.,20);
   
    h2dPhELDcaVsPt =     new TH2F("h2dPhELDcaVsPt",       "2D Dca vs pT Photonic E Like",              200,0,20, 100,-0.1,0.1);
    h2dPhELNSigEVsPt =   new TH2F("h2dPhELNSigEVsPt",     "2D nSigE vs pT Photonic E Like",            200,0,20, 289,-13,13);
    h2dPhELConvRVsPt =   new TH2F("h2dPhELConvRVsPt",     "2D Conversion Radius vs pT Photonic E Like",200,0,20, 100,0.,20);
    
    
    
    return kStOK;
}
//-----------------------------------------------------------------------------
StPicoNpeAnaMaker::~StPicoNpeAnaMaker()
{
    /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Finish()
{
    LOG_INFO << " StPicoNpeAnaMaker - writing data and closing output file " <<endm;
    mOutputFile->cd();
    // --------------- USER HISTOGRAM WRITE --------------------
    h2dIncEDcaVsPt->Write();
    h2dIncENSigEVsPt->Write();
    h2dIncEDcaVsPtCut->Write();
    h2dIncENSigEVsPtCut->Write();
    h2dIncEDcaVsPtCut2->Write();
    h2dIncENSigEVsPtCut2->Write();
    
    h2dPhEDcaVsPt->Write();
    h2dPhENSigEVsPt->Write();
    h2dPhEConvRVsPt->Write();
    
    h2dPhELDcaVsPt->Write();
    h2dPhELNSigEVsPt->Write();
    h2dPhELConvRVsPt->Write();
    

    mOutputFile->Close();
    
    return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Make()
{
    readNextEvent();
    
    if (!mPicoDstMaker)
    {
        LOG_WARN << " StPicoNpeAnaMaker - No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    
    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    
    if (!picoDst)
    {
        LOG_WARN << "StPicoNpeAnaMaker - No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    
    if(mPicoNpeEvent->runId() != picoDst->event()->runId() ||
       mPicoNpeEvent->eventId() != picoDst->event()->eventId())
    {
        LOG_ERROR <<" StPicoNpeAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoNpeEvent are not in sync."<<endm;
        exit(1);
    }
    
    // -------------- USER ANALYSIS -------------------------
    
    // check if good event (including bad run)
    if(!mNpeCuts->isGoodNpeEvent(const_cast<const StPicoDst*>(picoDst), NULL))
        return kStOk;
    StThreeVectorF pVtx = picoDst->event()->primaryVertex();

    // electron pair
    TClonesArray const * aElectronPair = mPicoNpeEvent->electronPairArray();
    for (int idx = 0; idx < aElectronPair->GetEntries(); ++idx)
    {
        // this is an example of how to get the ElectronPair pairs and their corresponsing tracks
        StElectronPair const* epair = (StElectronPair*)aElectronPair->At(idx);
        if(!mNpeCuts->isGoodElectronPair(epair)) continue;

        StPicoTrack const* electron = picoDst->track(epair->electronIdx());
        StPicoTrack const* partner = picoDst->track(epair->partnerIdx());
        
        StPhysicalHelixD eHelix = electron->dcaGeometry().helix();
        float dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
        float pt = electron->gPt();
        float nSigE = electron->nSigmaElectron();
        float pairPositionX = epair->positionX();
        float pairPositionY = epair->positionY();
        float convR = TMath::Sqrt((pairPositionX+0.2383) * (pairPositionX+0.2383) + (pairPositionY+0.1734) * (pairPositionY+0.1734));
        int pairCharge = electron->charge()+partner->charge();

        if (pairCharge==0) {
            h2dPhEDcaVsPt->Fill(pt, dca);
            h2dPhENSigEVsPt->Fill(pt, nSigE);
            h2dPhEConvRVsPt->Fill(pt, convR);
        }
        else {
            h2dPhELDcaVsPt->Fill(pt, dca);
            h2dPhELNSigEVsPt->Fill(pt, nSigE);
            h2dPhELConvRVsPt->Fill(pt, convR);
        }
        
        
    }
    
    
    // inclusive electron
    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (!track) continue;
        if (mNpeCuts->isGoodInclusiveElectron(track) && mNpeCuts->isBEMCElectron(track)) {
            StPhysicalHelixD eHelix = track->dcaGeometry().helix();
            float dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
            float pt = track->gPt();
            float nSigE = track->nSigmaElectron();
            
            h2dIncEDcaVsPt->Fill(pt, dca);
            h2dIncENSigEVsPt->Fill(pt, nSigE);

            if (mNpeCuts->isTPCElectron(track, 0, 3)){
                h2dIncEDcaVsPtCut->Fill(pt, dca);
                h2dIncENSigEVsPtCut->Fill(pt, nSigE);
            }
            if (mNpeCuts->isTPCElectron(track, 0, 3) && mNpeCuts->isBEMCElectron(track)){
                h2dIncEDcaVsPtCut2->Fill(pt, dca);
                h2dIncENSigEVsPtCut2->Fill(pt, nSigE);
            }
        }
    }

    return kStOK;
}
//-----------------------------------------------------------------------------
