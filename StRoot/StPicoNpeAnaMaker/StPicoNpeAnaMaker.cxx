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
    h1dEvent = new TH1I("h1dEvent", "Number of Events", 10, 0, 10);
    h1dEventZDCx = new TH1F("h1dEventZDCx", "ZDCx distribution", 1000, 0, 100000);
    h1dEventRefMult = new TH1I("h1dEventRefMult", "RefMult distribution", 1000, 0, 1000);
    h1dEventTrigger = new TH1I("h1dEventTrigger", "Trigger distribution", 30, 0, 30);
    h1dEventZDCxCut = new TH1F("h1dEventZDCxCut", "ZDCx distribution Cut", 1000, 0, 100000);
    h1dEventRefMultCut = new TH1I("h1dEventRefMultCut", "RefMult distribution Cut", 1000, 0, 1000);
    h1dEventTriggerCut = new TH1I("h1dEventTriggerCut", "Trigger distribution Cut", 30, 0, 30);
    
    
    h1dTrack = new TH1I("h1dTrack", "Number of Track", 10, 0, 10);
    
    h2dIncEDcaVsPt =    new TH2F("h2dIncEDcaVsPt",      "2D Dca vs pT Inclusive E",             200,0,20, 100,-0.1,0.1);
    h2dIncENSigEVsPt =  new TH2F("h2dIncENSigEVsPt",    "2D nSigE vs pT Inclusive E",           200,0,20, 289,-13,13);
    h2dIncEBsmdNEtaPt =  new TH2F("h2dIncEBsmdNEtaPt",    "2D nEta vs pT Inclusive E Cut2",           200,0,20, 10,-0.5,9.5);
    h2dIncEBsmdNPhiPt =  new TH2F("h2dIncEBsmdNPhiPt",    "2D nPhi vs pT Inclusive E Cut2",           200,0,20, 10,-0.5,9.5);
 
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
    
    hQaPt = new TH1F("hQaPt","hQaPt",200,0,20);
    hQaEta= new TH1F("hQaEta","hQaEta",100,-2,2);
    hQaDca= new TH1F("hQaDca","hQaDca",100,-5,5);
    hQaNHitFit = new TH1F("hQaNHitFit","hQaNHitFit",100,0,100);
    hQaNHitDedx = new TH1F("hQaNHitDedx","hQaNHitDedx",100,0,100);
    
    hQaPtCut = new TH1F("hQaPtCut","hQaPtCut",200,0,20);
    hQaEtaCut= new TH1F("hQaEtaCut","hQaEtaCut",100,-2,2);
    hQaDcaCut= new TH1F("hQaDcaCut","hQaDcaCut",100,-5,5);
    hQaNHitFitCut = new TH1F("hQaNHitFitCut","hQaNHitFitCut",100,0,100);
    hQaNHitDedxCut = new TH1F("hQaNHitDedxCut","hQaNHitDedxCut",100,0,100);
    
    
    
    
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
    h1dEvent->Write();
    h1dEventZDCx->Write();
    h1dEventRefMult->Write();
    h1dEventTrigger->Write();
    h1dEventZDCxCut->Write();
    h1dEventRefMultCut->Write();
    h1dEventTriggerCut->Write();

    h1dTrack->Write();
    h2dIncEDcaVsPt->Write();
    h2dIncENSigEVsPt->Write();
    h2dIncEDcaVsPtCut->Write();
    h2dIncENSigEVsPtCut->Write();
    h2dIncEDcaVsPtCut2->Write();
    h2dIncENSigEVsPtCut2->Write();
    
    h2dIncEBsmdNEtaPt->Write();
    h2dIncEBsmdNPhiPt->Write();
    
    h2dPhEDcaVsPt->Write();
    h2dPhENSigEVsPt->Write();
    h2dPhEConvRVsPt->Write();
    
    h2dPhELDcaVsPt->Write();
    h2dPhELNSigEVsPt->Write();
    h2dPhELConvRVsPt->Write();
    
    
    hQaPt->Write();
    hQaEta->Write();
    hQaDca->Write();
    hQaNHitFit->Write();
    hQaNHitDedx->Write();
    
    hQaPtCut->Write();
    hQaEtaCut->Write();
    hQaDcaCut->Write();
    hQaNHitFitCut->Write();
    hQaNHitDedxCut->Write();

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
    StThreeVectorF pVtx = picoDst->event()->primaryVertex();
    int mRefMult = picoDst->event()->refMult();
    float mZDCx = picoDst->event()->ZDCx();

    h1dEvent->Fill(0);
    h1dEventZDCx->Fill(mZDCx);
    h1dEventRefMult->Fill(mRefMult);
    for (int i=0; i<25; i++) if (picoDst->event()->triggerWord() >> i & 0x1) h1dEventTrigger->Fill(i);

    mNpeCuts->setPicoDst(const_cast<const StPicoDst*>(picoDst));
    if(!mNpeCuts->isGoodEvent(const_cast<const StPicoDst*>(picoDst), NULL))
        return kStOk;
    
    h1dEvent->Fill(1);
    h1dEventZDCxCut->Fill(mZDCx);
    h1dEventRefMultCut->Fill(mRefMult);
    for (int i=0; i<25; i++) if (picoDst->event()->triggerWord() >> i & 0x1) h1dEventTriggerCut->Fill(i);

    // electron pair
    TClonesArray const * aElectronPair = mPicoNpeEvent->electronPairArray();
    for (int idx = 0; idx < aElectronPair->GetEntries(); ++idx)
    {
        
        // this is an example of how to get the ElectronPair pairs and their corresponsing tracks
        StElectronPair const* epair = (StElectronPair*)aElectronPair->At(idx);
        if(!mNpeCuts->isGoodElectronPair(epair)) continue;

        StPicoTrack const* electron = picoDst->track(epair->electronIdx());
        StPicoTrack const* partner = picoDst->track(epair->partnerIdx());
        
        StPicoBTofPidTraits *tofPid = mNpeCuts->hasTofPid(electron);
        
        float beta;
        if (tofPid) {
            beta = tofPid->btofBeta();
            if (beta < 1e-4) {
                StThreeVectorF const btofHitPos = tofPid->btofHitPos();
                StPhysicalHelixD helix = electron->helix();
                float pathLength = tofPathLength(&pVtx, &btofHitPos, helix.curvature());
                float tof = tofPid->btof();
                beta = (tof > 0) ? pathLength / (tof * (C_C_LIGHT / 1.e9)) : std::numeric_limits<float>::quiet_NaN();
            }
        }
        else beta=999;
        
        if (fabs(1/beta -1) > 0.025) continue;
        StPicoBTofPidTraits *tofPid = mNpeCuts->hasTofPid(partner);
        
        float beta;
        if (tofPid) {
            beta = tofPid->btofBeta();
            if (beta < 1e-4) {
                StThreeVectorF const btofHitPos = tofPid->btofHitPos();
                StPhysicalHelixD helix = partner->helix();
                float pathLength = tofPathLength(&pVtx, &btofHitPos, helix.curvature());
                float tof = tofPid->btof();
                beta = (tof > 0) ? pathLength / (tof * (C_C_LIGHT / 1.e9)) : std::numeric_limits<float>::quiet_NaN();
            }
        }
        else beta=999;
        
        if (fabs(1/beta -1) > 0.025) continue;
        
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
        int jTrack=0;
        h1dTrack->Fill(jTrack);jTrack++;
        if (!track) continue;
        h1dTrack->Fill(jTrack);jTrack++;
        hQaPt->Fill(track->gPt());
        hQaEta->Fill(track->gMom(picoDst->event()->primaryVertex(), picoDst->event()->bField()).pseudoRapidity());
        hQaDca->Fill(track->dcaGeometry().helix().curvatureSignedDistance(picoDst->event()->primaryVertex().x(),picoDst->event()->primaryVertex().y()));
        hQaNHitFit->Fill(track->nHitsFit());
        hQaNHitDedx->Fill(track->nHitsDedx());
        if (mNpeCuts->isGoodInclusiveElectron(track)) {
            hQaPtCut->Fill(track->gPt());
            hQaEtaCut->Fill(track->gMom(picoDst->event()->primaryVertex(), picoDst->event()->bField()).pseudoRapidity());
            hQaDcaCut->Fill(track->dcaGeometry().helix().curvatureSignedDistance(picoDst->event()->primaryVertex().x(),picoDst->event()->primaryVertex().y()));
            hQaNHitFitCut->Fill(track->nHitsFit());
            hQaNHitDedxCut->Fill(track->nHitsDedx());
            
            StPhysicalHelixD eHelix = track->dcaGeometry().helix();
            float dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
            float pt = track->gPt();
            float nSigE = track->nSigmaElectron();

            h2dIncEDcaVsPt->Fill(pt, dca);
            h2dIncENSigEVsPt->Fill(pt, nSigE);

            h1dTrack->Fill(jTrack);jTrack++;
            if (mNpeCuts->isBEMCElectron(track)) {
                h1dTrack->Fill(jTrack);jTrack++;
                
                h2dIncEDcaVsPtCut->Fill(pt, dca);
                h2dIncENSigEVsPtCut->Fill(pt, nSigE);
                
                
                if (pt < 2 && mNpeCuts->isTPCElectron(track, 0, 1.5)){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                    h2dIncENSigEVsPtCut2->Fill(pt, nSigE);
                }
                else if (pt > 2 && pt < 4 && mNpeCuts->isTPCElectron(track, 0, 2.5)){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                    h2dIncENSigEVsPtCut2->Fill(pt, nSigE);
                }
                else if (pt > 4 && pt < 6 && mNpeCuts->isTPCElectron(track, 0.5, 3)){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                    h2dIncENSigEVsPtCut2->Fill(pt, nSigE);
                }
                else if (pt > 6 && pt < 8 && mNpeCuts->isTPCElectron(track, 1, 3)){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                    h2dIncENSigEVsPtCut2->Fill(pt, nSigE);
                }
                else if (pt > 8 && pt < 20 && mNpeCuts->isTPCElectron(track, 0, 3)){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                    h2dIncENSigEVsPtCut2->Fill(pt, nSigE);
                }
            }
        }
    }

    return kStOK;
}
//-----------------------------------------------------------------------------
