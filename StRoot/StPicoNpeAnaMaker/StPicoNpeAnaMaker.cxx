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
#include "StCuts.h"


ClassImp(StPicoNpeAnaMaker)

StPicoNpeAnaMaker::StPicoNpeAnaMaker(char const * name,char const * inputFilesList,
                                     char const * outName, StPicoDstMaker* picoDstMaker):
StMaker(name),mPicoDstMaker(picoDstMaker),mPicoNpeEvent(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
mOutputFile(NULL), mChain(NULL), mEventCounter(0)
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
    
    // -------------- USER VARIABLES -------------------------
    hEvent = new TH1F("hEvent","hEvent",10,0,10);
    hZDCx = new TH1F("hZDCx","hZDCx",100000,0,100000);
    hTrigger = new TH1I("hTrigger","hTrigger",30,0,30);
    
    tInc = new TTree("tInc","tree for electron");
    tPhE = new TTree("tPhE","tree for photonic electron");
    
    setTree(tInc,"T");
    setTree(tPhE,"P");
    
    const int nbin = 6;
    int ptbin[nbin+1] = {1.5, 1.7, 2.0, 2.5, 3.5, 5.5, 10.};
    for (int i=0; i<25; i++) {
        hRefMult[i] = new TH1F(Form("hRefMult_%d",i),Form("hRefMult_%d",i),1000,0,1000);
        for (int j=0; j<5; j++) { // 0: no cut, 1: ...
            hDcaByPt[i][j] = new TH2F(Form("hDcaByPt_%d_%d",i,j),Form("hDcaByPt_%d_%d",i,j),nbin,ptbin,100,-0.1,0.1);
            hNSigEByPt[i][j] = new TH2F(Form("hNSigEByPt_%d_%d",i,j),Form("hNSigEByPt_%d_%d",i,j),nbin,ptbni,289,-13,13);
            hEOverPByPt[i][j] = new TH2F(Form("hEOverPByPt_%d_%d",i,j),Form("hEOverPByPt_%d_%d",i,j),nbin,ptbin,100,0,4);
            hNEtaByPt[i][j] = new TH2F(Form("hNEtaByPt_%d_%d",i,j),Form("hNEtaByPt_%d_%d",i,j),nbin,ptbin,10,0,10);
            hNPhiByPt[i][j] = new TH2F(Form("hNPhiByPt_%d_%d",i,j),Form("hNPhiByPt_%d_%d",i,j),nbin,ptbin,10,0,10);
            hPairMassByPt[i][j] = new TH2F(Form("hPairMassByPt_%d_%d",i,j),Form("hPairMassByPt_%d_%d",i,j),nbin,ptbin,100,0,0.4);

        }
    }
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
    hEvent->Write();
    hZDCx->Write();
    hHFTInnerOuter->Write();
    hHFTInner->Write();
    hHFTOuter->Write();
    hTrigger->Write();
    
    tInc->Write();
    tPhE->Write();
    
    for (int i=0; i<25; i++) {
        hRefMult[i]->Write();
    }
    
    mOutputFile->Close();
    
    return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Make()
{

    readNextEvent();
    hEvent->Fill(0);

    if (!mPicoDstMaker)
    {
        LOG_WARN << " StPicoNpeAnaMaker - No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    hEvent->Fill(1);

    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    
    if (!picoDst)
    {
        LOG_WARN << "StPicoNpeAnaMaker - No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    hEvent->Fill(2);

    if(mPicoNpeEvent->runId() != picoDst->event()->runId() ||
       mPicoNpeEvent->eventId() != picoDst->event()->eventId())
    {
        LOG_ERROR <<" StPicoNpeAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoNpeEvent are not in sync."<<endm;
        exit(1);
    }
    hEvent->Fill(3);

  //  if (!isGoodEvent()) return kStOK;
    hEvent->Fill(4);

    // -------------- USER ANALYSIS -------------------------
    // Event informaiton
    mRefMult = std::numeric_limits<Int_t>::quiet_NaN();
    mZDCx = std::numeric_limits<Int_t>::quiet_NaN();
    
    mRefMult = picoDst->event()->refMult();
    mZDCx = picoDst->event()->ZDCx();

    bField = picoDst->event()->bField();
    pVtx = picoDst->event()->primaryVertex();

    hZDCx->Fill(mZDCx);
    hHFTInnerOuter->Fill(picoDst->event()->numberOfPxlInnerHits(),picoDst->event()->numberOfPxlOuterHits());
    hHFTInner->Fill(picoDst->event()->numberOfPxlInnerHits());
    hHFTOuter->Fill(picoDst->event()->numberOfPxlOuterHits());
    
    for (int i=0;i<25;i++)
    {
        if (picoDst->event()->triggerWord()>>i & 1)
        {
            hTrigger->Fill(i);
            hRefMult[i]->Fill(mRefMult);
        }
    }
    isHTEvents = 0;
    if (picoDst->event()->triggerWord()>>0 & 0x7F) isHTEvents += 1;
    if (picoDst->event()->triggerWord()>>7 & 0xFFF) isHTEvents += 2;
    if (picoDst->event()->triggerWord()>>19 & 0x3F) isHTEvents += 4;

    
    // hadrons & inclusive electron with StPicoTrack
    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (!track) continue;
        if (isGoodTrack(track)) {
            setVariables(track);
            tInc->Fill();
        }
    }

    
    // Photonic Electron
    TClonesArray const * aElectronPair = mPicoNpeEvent->electronPairArray();
    for (int idx = 0; idx < aElectronPair->GetEntries(); ++idx)
    {
        // this is an example of how to get the ElectronPair pairs and their corresponsing tracks
        StElectronPair * epair = (StElectronPair*)aElectronPair->At(idx);
        if (isGoodPair(epair))
        {
            setVariables(epair);
            tPhE->Fill();
        }
    }
    
    return kStOK;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodEvent() const
{
    return mPicoDstMaker->picoDst()->event()->triggerWord()>>cutsAna::trigger & cutsAna::triggerLength;
    
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodPair(StElectronPair const* const epair) const
{
    if(!epair) return false;
    
    StPicoTrack const* electron = mPicoDstMaker->picoDst()->track(epair->electronIdx());
    StPicoTrack const* partner = mPicoDstMaker->picoDst()->track(epair->partnerIdx());
    
    return
    isGoodTagged(electron) &&
    isGoodPartner(partner) &&
    epair->pairMass() < cutsAna::pairMass &&
    epair->pairDca() < cutsAna::pairDca
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
    return
    (!cutsAna::trackRequireHFT || trk->isHFTTrack()) &&
    trk->nHitsFit() >= cutsAna::trackNHitsFit &&
    trk->nHitsDedx() >= cutsAna::trackNhitsDedx &&
    fabs(trk->gMom(pVtx, bField).pseudoRapidity()) <= cutsAna::trackEta &&
    trk->gPt() >= cutsAna::trackPt
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodTagged(StPicoTrack const * const trk) const
{
    if (!isGoodTrack(trk)) return false;
    return fabs(trk->nSigmaElectron()) <= cutsAna::taggedNSigElectron;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodPartner(StPicoTrack const * const trk) const
{
    return
    trk->nHitsFit() >= cutsAna::partnerNHitsFit &&
    fabs(trk->gMom(pVtx, bField).pseudoRapidity()) <= cutsAna::partnerEta &&
    trk->gPt() >= cutsAna::partnerPt &&
    fabs(trk->nSigmaElectron()) <= cutsAna::partnerNSigElectron
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodTofTrack(StPicoTrack const * const trk) const
{
    if (trk->bTofPidTraitsIndex() < 0) return false;
    // StPicoBTofPidTraits * Tof = mPicoDstMaker->picoDst()->btofPidTraits(trk->bTofPidTraitsIndex());
    
    return
    // TMath::Abs(1-1/Tof->btofBeta()) < cutsAna::tofBeta
    true
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodEmcTrack(StPicoTrack const * const trk) const
{
    if (trk->emcPidTraitsIndex() < 0) return false;
    // StPicoEmcPidTraits * Emc =  mPicoDstMaker->picoDst()->emcPidTraits(trk->emcPidTraitsIndex());
    
    return
    /*   Emc->nPhi() > cutsAna::emcNPhi &&
     Emc->nEta() > cutsAna::emcNEta &&
     TMath::Abs(Emc->phiDist()) < cutsAna::emcPhiDist &&
     TMath::Abs(Emc->zDist()) < cutsAna::emcZDist &&
     TMath::Sqrt(Emc->phiTowDist()*Emc->phiTowDist() + Emc->etaTowDist()*Emc->etaTowDist()) < cutsAna::emcAssDist &&
     Emc->e0()/trk->dcaGeometry().momentum().mag() < cutsAna::emcEoverPHigh &&
     Emc->e0()/trk->dcaGeometry().momentum().mag() > cutsAna::emcEoverPLow
     */
    true
    ;
}
//-----------------------------------------------------------------------------
void StPicoNpeAnaMaker::setTree(TTree * tree, TString opt)
{
    tree->Branch("dca",&dca,"dca/F");
    tree->Branch("pt",&pt,"pt/F");
    tree->Branch("eta",&eta,"eta/F");
    tree->Branch("nsige",&nsige,"nsige/F");
    tree->Branch("nsigpion",&nsigpion,"nsigpion/F");
    tree->Branch("beta",&beta,"beta/F");
    tree->Branch("e",&e,"e/F");
    tree->Branch("e0",&e0,"e0/F");
    tree->Branch("e1",&e1,"e1/F");
    tree->Branch("e2",&e2,"e2/F");
    tree->Branch("e3",&e3,"e3/F");
    tree->Branch("neta",&neta,"neta/b");
    tree->Branch("nphi",&nphi,"nphi/b");
    tree->Branch("phiDist",&phiDist,"phiDist/F");
    tree->Branch("zDist",&zDist,"zDist/F");
    tree->Branch("etaTowDist",&etaTowDist,"etaTowDist/F");
    tree->Branch("phiTowDist",&phiTowDist,"phiTowDist/F");
    tree->Branch("mZDCx",&mZDCx,"mZDCx/s");
    tree->Branch("mRefMult",&mRefMult,"mRefMult/s");
    tree->Branch("isHTEvents",&isHTEvents,"isHTEvents/b");

    if (opt=="T") {}
    else if (opt=="P")
    {
        tree->Branch("parnter_pt",&partner_pt,"partner_pt/F");
        tree->Branch("partner_nsige",&partner_nsige,"partner_nsige/F");
        tree->Branch("pairAngle3d",&pairAngle3d,"pairAngle3d/F");
        tree->Branch("pairAnglePhi",&pairAnglePhi,"pairAnglePhi/F");
        tree->Branch("pairAngleTheta",&pairAngleTheta,"pairAngleTheta/F");
        tree->Branch("pairMass",&pairMass,"pairMass/F");
        tree->Branch("pairCharge",&pairCharge,"pairCharge/B");
        tree->Branch("pairDca",&pairDca,"pairDca/F");
        tree->Branch("pairPositionX",&pairPositionX,"pairPositionX/F");
        tree->Branch("pairPositionY",&pairPositionY,"pairPositionY/F");
    }
    else LOG_WARN << "Select tree options. (T is for track, P is for pair.)" << endm;
}
//-----------------------------------------------------------------------------
void StPicoNpeAnaMaker::initVariables()
{
    beta = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
    e0 = std::numeric_limits<float>::quiet_NaN();
    e1 = std::numeric_limits<float>::quiet_NaN();
    e2 = std::numeric_limits<float>::quiet_NaN();
    e3 = std::numeric_limits<float>::quiet_NaN();
    nphi = std::numeric_limits<unsigned char>::quiet_NaN();
    neta = std::numeric_limits<unsigned char>::quiet_NaN();
    phiDist = std::numeric_limits<float>::quiet_NaN();
    zDist = std::numeric_limits<float>::quiet_NaN();
    etaTowDist = std::numeric_limits<float>::quiet_NaN();
    phiTowDist = std::numeric_limits<float>::quiet_NaN();
}
//-----------------------------------------------------------------------------
void StPicoNpeAnaMaker::setVariables(StPicoTrack * track)
{
    initVariables();
    
    StPhysicalHelixD eHelix = track->dcaGeometry().helix();
    
    // Track
    dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
    pt = track->gPt();
    eta = track->gMom(pVtx, bField).pseudoRapidity();
    
    // PID
    nsige = track->nSigmaElectron();
    nsigpion = track->nSigmaPion();
    if (isGoodTofTrack(track))  {
        StPicoBTofPidTraits * Tof = picoDst->btofPidTraits(track->bTofPidTraitsIndex());
        
        // start global beta calculation
        double newBeta = Tof->btofBeta();
        if(newBeta<1e-4 && fabs(dca)<3.) {
            StThreeVectorF btofHitPos = Tof->btofHitPos();
            float L = tofPathLength(&pVtx, &btofHitPos, eHelix.curvature());
            float tof = Tof->btof();
            if (tof>0) newBeta = L/(tof*(C_C_LIGHT/1.e9));
        }
        beta = 1./newBeta;
        // end global beta calculation
    }
    if (isGoodEmcTrack(track)) {
        StPicoEmcPidTraits * Emc =  picoDst->emcPidTraits(track->emcPidTraitsIndex());
        e = Emc->e();
        e0 = Emc->e0();
        e1 = Emc->e1();
        e2 = Emc->e2();
        e3 = Emc->e3();
        nphi = Emc->nPhi();
        neta = Emc->nEta();
        phiDist = Emc->phiDist();
        zDist = Emc->zDist();
        etaTowDist = Emc->etaTowDist();
        phiTowDist = Emc->phiTowDist();
    }
    
}
//-----------------------------------------------------------------------------
void StPicoNpeAnaMaker::setVariables(StElectronPair * epair)
{
    pairMass = epair->pairMass();
    pairDca = epair->pairDca();
    
    StPicoTrack * electron = picoDst->track(epair->electronIdx());
    StPicoTrack * partner = picoDst->track(epair->partnerIdx());
    
    setVariables(electron);
    
    // calculate Lorentz vector of electron-partner pair
    StPhysicalHelixD electronHelix = electron->dcaGeometry().helix();
    StPhysicalHelixD partnerHelix = partner->dcaGeometry().helix();
    pair<double,double> ss = electronHelix.pathLengths(partnerHelix);
    StThreeVectorF const electronMomAtDca = electronHelix.momentumAt(ss.first, bField * kilogauss);
    StThreeVectorF const partnerMomAtDca = partnerHelix.momentumAt(ss.second, bField * kilogauss);
    StLorentzVectorF const electronFourMom(electronMomAtDca, electronMomAtDca.massHypothesis(M_ELECTRON));
    StLorentzVectorF const partnerFourMom(partnerMomAtDca, partnerMomAtDca.massHypothesis(M_ELECTRON));
    StLorentzVectorF const epairFourMom = electronFourMom + partnerFourMom;
    StPhysicalHelixD eHelix = electron->dcaGeometry().helix();
    
    pairMass = epairFourMom.m();
    pairAngle3d = electronMomAtDca.angle(partnerMomAtDca);
    pairAnglePhi = fabs(electronMomAtDca.phi() - partnerMomAtDca.phi());
    pairAngleTheta = fabs(electronMomAtDca.theta() - partnerMomAtDca.theta());
    pairCharge = electron->charge()+partner->charge();
    pairPositionX = epair->positionX();
    pairPositionY = epair->positionY();
    
    partner_pt = partner->gPt();
    partner_nsige = partner->nSigmaElectron();
}


