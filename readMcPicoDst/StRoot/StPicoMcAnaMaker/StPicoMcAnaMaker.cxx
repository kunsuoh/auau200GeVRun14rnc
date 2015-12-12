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

StPicoMcAnaMaker::StPicoMcAnaMaker(char const * name,char const * inputFilesList,
                                     char const * outName, StPicoDstMaker* picoDstMaker):
StMaker(name),mPicoDstMaker(picoDstMaker),mPicoNpeEvent(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
mOutputFile(NULL), mChain(NULL), mEventCounter(0), mNpeCuts(NULL)
{}

Int_t StPicoMcAnaMaker::Init()
{
    mPicoNpeEvent = new StPicoNpeEvent();
    
    mChain = new TChain("T");
    std::ifstream listOfFiles(mInputFileList.Data());
    if (listOfFiles.is_open())
    {
        std::string file;
        while (getline(listOfFiles, file))
        {
            LOG_INFO << "StPicoMcAnaMaker - Adding :" << file << endm;
            mChain->Add(file.c_str());
        }
    }
    else
    {
        LOG_ERROR << "StPicoMcAnaMaker - Could not open list of files. ABORT!" << endm;
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
    
    h2dIncEDcaVsPtCut =    new TH2F("h2dIncEDcaVsPtCut",      "2D Dca vs pT Inclusive E Cut",             200,0,20, 100,-0.1,0.1);
    
    h2dIncEDcaVsPtCut2 =    new TH2F("h2dIncEDcaVsPtCut2",      "2D Dca vs pT Inclusive E Cut2",             200,0,20, 100,-0.1,0.1);
    
    
    h2dPhEDcaVsPt =     new TH2F("h2dPhEDcaVsPt",       "2D Dca vs pT Photonic E",              200,0,20, 100,-0.1,0.1);
    h2dPhEConvRVsPt =   new TH2F("h2dPhEConvRVsPt",     "2D Conversion Radius vs pT Photonic E",200,0,20, 100,0.,20);
    
    h2dPhELDcaVsPt =     new TH2F("h2dPhELDcaVsPt",       "2D Dca vs pT Photonic E Like",              200,0,20, 100,-0.1,0.1);
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
    
    h2dPhEConvRVsZ = new TH2D("h2dPhEConvRVsZ","h2dPhEConvRVsZ",100,-20,20, 1000, 0, 100);
    h2dPhEConvXYZ = new TH3D("h2dPhEConvXYZ","h2dPhEConvXYZ",100,-50,50, 100,-50,50, 40,-20,20);
    h2dPhEInvMassvsZ = new TH2D("h2dPhEInvMassvsZ","h2dPhEInvMassvsZ",100,-20,20, 100, 0, 0.5);
    
    h2dPhELConvRVsZ = new TH2D("h2dPhELConvRVsZ","h2dPhELConvRVsZ",100,-20,20, 1000, 0, 100);
    h2dPhELConvXYZ = new TH3D("h2dPhELConvXYZ","h2dPhELConvXYZ",100,-50,50, 100,-50,50, 40,-20,20);
    h2dPhELInvMassvsZ = new TH2D("h2dPhELInvMassvsZ","h2dPhELInvMassvsZ",100,-20,20, 100, 0, 0.5);
    
    h2dPhEConvRVsZ_HFT = new TH2D("h2dPhEConvRVsZ_HFT","h2dPhEConvRVsZ_HFT",100,-20,20, 1000, 0, 100);
    h2dPhEConvXYZ_HFT = new TH3D("h2dPhEConvXYZ_HFT","h2dPhEConvXYZ_HFT",100,-50,50, 100,-50,50, 40,-20,20);
    h2dPhEInvMassvsZ_HFT = new TH2D("h2dPhEInvMassvsZ_HFT","h2dPhEInvMassvsZ_HFT",100,-20,20, 100, 0, 0.5);
    
    h2dPhELConvRVsZ_HFT = new TH2D("h2dPhELConvRVsZ_HFT","h2dPhELConvRVsZ_HFT",100,-20,20, 1000, 0, 100);
    h2dPhELConvXYZ_HFT = new TH3D("h2dPhELConvXYZ_HFT","h2dPhELConvXYZ_HFT",100,-50,50, 100,-50,50, 40,-20,20);
    h2dPhELInvMassvsZ_HFT = new TH2D("h2dPhELInvMassvsZ_HFT","h2dPhELInvMassvsZ_HFT",100,-20,20, 100, 0, 0.5);
    
    nt = new TNtuple("nt","electron pair ntuple","pt1:pt2:phiV:openangle:v0x:v0y:v0z:phi:eta:mass:pairDca");
    
    h2dPhEMassVsPt = new TH2D("h2dPhEMassVsPt","h2dPhEMassVsPt",100,0,0.5,100,0,20);
    h2dPhELMassVsPt = new TH2D("h2dPhELMassVsPt","h2dPhELMassVsPt",100,0,0.5,100,0,20);
    h2dPhEPairDcaVsPt = new TH2D("h2dPhEPairDcaVsPt","h2dPhEPairDcaVsPt",100,0,0.5,100,0,20);
    h2dPhELPairDcaVsPt = new TH2D("h2dPhELPairDcaVsPt","h2dPhELPairDcaVsPt",100,0,0.5,100,0,20);
    
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
    
    nt->Write();
    
    h2dPhEMassVsPt->Write();
    h2dPhELMassVsPt->Write();
    h2dPhEPairDcaVsPt->Write();
    h2dPhELPairDcaVsPt->Write();
    
    
    h2dPhEConvRVsZ->Write();
    h2dPhEConvXYZ->Write();
    h2dPhEInvMassvsZ->Write();
    
    h2dPhELConvRVsZ->Write();
    h2dPhELConvXYZ->Write();
    h2dPhELInvMassvsZ->Write();
    
    
    h2dPhEConvRVsZ_HFT->Write();
    h2dPhEConvXYZ_HFT->Write();
    h2dPhEInvMassvsZ_HFT->Write();

    h2dPhELConvRVsZ_HFT->Write();
    h2dPhELConvXYZ_HFT->Write();
    h2dPhELInvMassvsZ_HFT->Write();
    
    h1dEvent->Write();
    h1dEventZDCx->Write();
    h1dEventRefMult->Write();
    h1dEventTrigger->Write();
    h1dEventZDCxCut->Write();
    h1dEventRefMultCut->Write();
    h1dEventTriggerCut->Write();
    
    h1dTrack->Write();
    h2dIncEDcaVsPt->Write();
    h2dIncEDcaVsPtCut->Write();
    h2dIncEDcaVsPtCut2->Write();
    
    h2dPhEDcaVsPt->Write();
    h2dPhEConvRVsPt->Write();
    
    h2dPhELDcaVsPt->Write();
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
    
    if(mPicoNpeEvent->runId() != picoDst->event()->runId() ||
       mPicoNpeEvent->eventId() != picoDst->event()->eventId())
    {
        LOG_ERROR <<" StPicoMcAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
        LOG_ERROR <<" StPicoMcAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoNpeEvent are not in sync."<<endm;
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
        StPicoTrack const* electron = picoDst->track(epair->electronIdx());
        StPicoTrack const* partner = picoDst->track(epair->partnerIdx());
        
        if (!mNpeCuts->isGoodTaggedElectron(electron) || !mNpeCuts->isGoodPartnerElectron(partner)) continue;
        
        if (electron->charge()+partner->charge() == 0){
            h2dPhEMassVsPt->Fill(epair->pairMass(),electron->gPt());
            h2dPhEPairDcaVsPt->Fill(epair->pairDca(),electron->gPt());
        }
        else {
            h2dPhELMassVsPt->Fill(epair->pairMass(),electron->gPt());
            h2dPhELPairDcaVsPt->Fill(epair->pairDca(),electron->gPt());
        }
        
        if(!mNpeCuts->isGoodElectronPair(epair)) continue;
        
        
        StPhysicalHelixD electronHelix = electron->dcaGeometry().helix();
        StPhysicalHelixD partnerHelix = partner->dcaGeometry().helix();
        
        pair<double,double> ss = electronHelix.pathLengths(partnerHelix);
        StThreeVectorD kAtDcaToPartner = electronHelix.at(ss.first);
        StThreeVectorD pAtDcaToElectron = partnerHelix.at(ss.second);
        
        StThreeVectorF const electronMomAtDca = electronHelix.momentumAt(ss.first, picoDst->event()->bField() * kilogauss);
        StThreeVectorF const partnerMomAtDca = partnerHelix.momentumAt(ss.second, picoDst->event()->bField() * kilogauss);
        
        StLorentzVectorF const electronFourMom(electronMomAtDca, electronMomAtDca.massHypothesis(M_ELECTRON));
        StLorentzVectorF const partnerFourMom(partnerMomAtDca, partnerMomAtDca.massHypothesis(M_ELECTRON));
        StLorentzVectorF const epairFourMom = electronFourMom + partnerFourMom;
        
        double phiV = 0;
        double openangle = 0;
        if (electron->Charge() > 0 && partner->Charge() < 0) phiCalculation(electronFourMom,partnerFourMom,picoDst->event()->bField() > 0 ? 1. : -1.,phiV,openangle);
        else if (electron->Charge() < 0 && partner->Charge() > 0)  phiCalculation(partnerFourMom,electronFourMom,picoDst->event()->bField() > 0 ? 1. : -1.,phiV,openangle);
        else {}
        float pt1 = electron->gPt() * electron->charge();
        float pt2 = partner->gPt() * partner->charge();
        float eta = epairFourMom.pseudoRapidity();
        float phi = epairFourMom.phi();
        float v0x = epair->positionX();
        float v0y = epair->positionY();
        float v0z = epair->positionZ();
        float mass = epairFourMom.m();
        float pairDca = epair->pairDca();
        
        nt->Fill(pt1,pt2,phiV,openangle,v0x,v0y,v0z,phi,eta,mass,pairDca);
    }
    
    
    
    
    
    
    std::vector<unsigned short> idxPicoTaggedEs;
    std::vector<unsigned short> idxPicoPartnerEs;
    
    // inclusive electron
    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (mNpeCuts->isGoodTaggedElectron(track))  idxPicoTaggedEs.push_back(iTrack);
        if (mNpeCuts->isGoodPartnerElectron(track)) idxPicoPartnerEs.push_back(iTrack);
        
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
            
            h2dIncEDcaVsPt->Fill(pt, dca);
            
            h1dTrack->Fill(jTrack);jTrack++;
                h1dTrack->Fill(jTrack);jTrack++;
                
                h2dIncEDcaVsPtCut->Fill(pt, dca);
                
                
                if (pt < 2 ){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                }
                else if (pt > 2 && pt < 4 ){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                }
                else if (pt > 4 && pt < 6){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                }
                else if (pt > 6 && pt < 8){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                }
                else if (pt > 8 && pt < 20){
                    h1dTrack->Fill(jTrack);jTrack++;
                    h2dIncEDcaVsPtCut2->Fill(pt, dca);
                }
        }
    }
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