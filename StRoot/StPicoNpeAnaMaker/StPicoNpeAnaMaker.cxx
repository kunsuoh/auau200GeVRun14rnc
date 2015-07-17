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
#include "StPicoPrescales/StPicoPrescales.h"
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
mOutputFile(NULL), mChain(NULL), mEventCounter(0), mPrescales(NULL)
{}

Int_t StPicoNpeAnaMaker::Init()
{
    mPrescales = new StPicoPrescales(cutsAna::prescalesFilesDirectoryName);
    
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
    hCheckDoubleTrigger = new TH1I("hCheckDoubleTrigger","hCheckDoubleTrigger",20,0,20);
    hTrigger = new TH1I("hTrigger","hTrigger",30,0,30);
    hEvent = new TH1F("hEvent","hEvent",10,0,10);
    hZDCx = new TH1F("hZDCx","hZDCx",1000,0,100000);
    hZDCxWt = new TH1F("hZDCxWt","hZDCxWt",1000,0,100000);
    
    for (int i=0; i<2; i++) {
        hRefMult[i] = new TH1F(Form("hRefMult_%d",i),Form("hRefMult_%d",i),1000,0,1000);
        hRefMultWt[i] = new TH1F(Form("hRefMultWt_%d",i),Form("hRefMultWt_%d",i),1000,0,1000);
    }

    setHistogram(6,6,6,4);

    return kStOK;
}
//-----------------------------------------------------------------------------
StPicoNpeAnaMaker::~StPicoNpeAnaMaker()
{
    /*  */
    if (mPrescales)
        delete mPrescales;
    mPrescales = NULL;
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Finish()
{
    LOG_INFO << " StPicoNpeAnaMaker - writing data and closing output file " <<endm;
    mOutputFile->cd();

    // --------------- USER HISTOGRAM WRITE --------------------
    hEvent->Write();
    hZDCx->Write();
    hZDCxWt->Write();
    hTrigger->Write();
    hCheckDoubleTrigger->Write();
    
    hRefMult[0]->Write();
    hRefMult[1]->Write();
    hRefMultWt[0]->Write();
    hRefMultWt[1]->Write();
    
    for (int j=0;j<6;j++) // PID
        for (int i=1;i<6;i++) // PT
            for (int k=0;k<3;k++) // Particle:Type
                for (int l=0;l<2;l++) // Histograms
                    histo[i][j][k][l]->Write();

    for (int i=1;i<6;i++) histoTofMass[i]->Write(); // tofmass
    for (int i=1;i<6;i++) for (int j=0;j<3;j++) histoNSigE[i][j]->Write();
    for (int j=2;j<6;j++) for (int i=1;i<6;i++) histo[i][j][2][2]->Write();
    for (int j=0;j<6;j++) for (int i=1;i<6;i++) for (int k=0;k<2;k++) histo[i][j][k][3]->Write();
 
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
    


    // -------------- USER ANALYSIS -------------------------
    // Event informaiton
    mRefMult = std::numeric_limits<Int_t>::quiet_NaN();
    mZDCx = std::numeric_limits<Int_t>::quiet_NaN();
    
    mRefMult = picoDst->event()->refMult();
    mZDCx = picoDst->event()->ZDCx();

    bField = picoDst->event()->bField();
    pVtx = picoDst->event()->primaryVertex();


    isHTEvents = 0;
    if (picoDst->event()->triggerWord()>>0 & 0x7FF) isHTEvents += 1;
    if (picoDst->event()->triggerWord()>>19 & 0x3F) isHTEvents += 2;
    int checkDoubleTrigger = 0;
    int checkDoubleTrigger18 = 10;
    weight = 1;

    for (int i=0; i<30; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        //cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        //weight = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        checkDoubleTrigger++;
        if (i>18) checkDoubleTrigger18++;
    }
    hCheckDoubleTrigger->Fill(checkDoubleTrigger);
    hCheckDoubleTrigger->Fill(checkDoubleTrigger18);
    
    hEvent->Fill(4);
    hEvent->Fill(5,weight);

    hZDCx->Fill(mZDCx);
    hZDCxWt->Fill(mZDCx,weight);
    
    for (int i=0; i<2; i++) if (isHTEvents >> i & 0x1) {
        hRefMult[i]->Fill(mRefMult);
        hRefMultWt[i]->Fill(mRefMult,weight);
    }

    // hadrons & inclusive electron with StPicoTrack
    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (!track) continue;
        if (isGoodTrack(track)) {
            setVariables(track);
            fillHistogram(2); // electron
            fillHistogram(3); // hadron
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
            if (pairCharge == 0) fillHistogram(0); // US
            else fillHistogram(1);                 // LS
            
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
    StPhysicalHelixD eHelix = trk->dcaGeometry().helix();
    
    return
    (!cutsAna::trackRequireHFT || trk->isHFTTrack()) &&
    trk->nHitsFit() >= cutsAna::trackNHitsFit &&
    trk->nHitsDedx() >= cutsAna::trackNhitsDedx &&
    fabs(trk->gMom(pVtx, bField).pseudoRapidity()) <= cutsAna::trackEta &&
    trk->gPt() >= cutsAna::trackPt &&
    eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y()) < cutsAna::trackDca
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
void StPicoNpeAnaMaker::initVariables()
{
    beta = std::numeric_limits<float>::quiet_NaN();
    tofmass = std::numeric_limits<float>::quiet_NaN();
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
            if (tof>0) {
                newBeta = L/(tof*(C_C_LIGHT/1.e9));
            }
        }
        beta = 1./newBeta;
        tofmass = TMath::Sqrt(beta*beta-1)*pt*TMath::CosH(eta);
        // end global beta calculation
    }
    if (isGoodEmcTrack(track)) {
        StPicoEmcPidTraits * Emc =  picoDst->emcPidTraits(track->emcPidTraitsIndex());
        e = Emc->e();
        e0 = Emc->e0();
        e1 = Emc->e1();
        e2 = Emc->e2();
        e3 = Emc->e3();
        eoverp = e0/pt/TMath::CosH(eta);
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
//-----------------------------------------------------------------------------
void StPicoNpeAnaMaker::setHistogram(int nptbin,int npid,int ntype,int nhisto)
{
    
    double ptbin[10] = {0, 1.5, 1.8, 2.5, 4.0, 6.5, 10.};
    TString pid[10] = {"Tpc","TpcTof","TpcBemc","TpcBemcBsmd"};
    TString type[10] = {"PhEUS","PhELS","IncE","Pion","Kaon","Proton"};
    TString histoname[10] = {"nSigE","DCA","DCAafterPIDcut","pairDca"};
    
    int binHisto[10] = {289, 100, 100,100};
    double minHisto[10] = {-13, -0.1, -0.1,0};
    double maxHisto[10] = {13, 0.1, 0.1,0.5};
    
    for (int i=0;i<nptbin;i++){
        histoTofMass[i] = new TH1F(Form("histoTofMass_%d",i),Form("histoTofMass_%d",i),1000,-0.5,2.5);
        for (int j=0;j<npid;j++) histoNSigE[i][j] =  new TH1F(Form("histoNSigE_%d_%d",i,j),Form("histoNSigE_%d_%d",i,j),1301,-13,13);
        for (int j=0;j<npid;j++)
            for (int k=0;k<ntype;k++)
                for (int l=0;l<nhisto;l++)
                    histo[i][j][k][l] = new TH1F(
                                                 Form("histo_%d_%d_%d_%d", i,j,k,l),
                                                 Form("histo_pT%.1f_%.1f_%s_%s_%s",
                                                      ptbin[i],
                                                      ptbin[i+1],
                                                      pid[j].Data(),
                                                      type[k].Data(),
                                                      histoname[l].Data()
                                                      ),
                                                 binHisto[l],
                                                 minHisto[l],
                                                 maxHisto[l]
                                                 );

    }
    
    
}
//-------------------------------------------------------------------------------
int StPicoNpeAnaMaker::getPtBin(double pt) {
    int nPt = 0;
    if (pt < 1.5) nPt = 0;
    else if (pt < 1.8) nPt = 1;
    else if (pt < 2.5) nPt = 2;
    else if (pt < 4.0) nPt = 3;
    else if (pt < 6.5) nPt = 4;
    else if (pt < 10.) nPt = 5;
    else nPt = 0;
    
    return nPt;
}
//-------------------------------------------------------------------------------
void StPicoNpeAnaMaker::fillHistogram(int iPt, int iPid, int iType){
    histo[(const int)iPt][(const int)iPid][(const int)iType][0]->Fill(nsige,weight);
    if (iType < 2) {
        if (nsige > 0) {
            histo[(const int)iPt][(const int)iPid][(const int)iType][1]->Fill(dca,weight);
            histo[(const int)iPt][(const int)iPid][(const int)iType][3]->Fill(pairMass,weight);
        }
    }
    else {
        histo[(const int)iPt][(const int)iPid][(const int)iType][1]->Fill(dca,weight);
        histo[(const int)iPt][(const int)iPid][(const int)iType][3]->Fill(pairMass,weight);
    }
    if (iType==2) {
        float pidCutLw[2][6];
        float pidCutHi[2][6];
        pidCutLw[0]={0, -1.2, -1.2, -1.0, -1.0, 0};
        pidCutHi[0]={0, 1.8, 2.5, 3.0, 3.0, 0};
        pidCutLw[1]={0, -1.5, -1.4, -1.5, -1.1, 0};
        pidCutHi[1]={0, 1.8, 2.5, 3.0, 3.0, 0};
        if (iPid == 2 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) histo[iPt][2][2][2]->Fill(dca,weight);
        if (iPid == 3 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) histo[iPt][3][2][2]->Fill(dca,weight);
        if (iPid == 4 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) histo[iPt][4][2][2]->Fill(dca,weight);
        if (iPid == 5 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) histo[iPt][5][2][2]->Fill(dca,weight);

    }
}

    //-------------------------------------------------------------------------------
void StPicoNpeAnaMaker::fillHistogram(int iType){

    int iPt = getPtBin(pt);
    
    if (isHTEvents >> 0 & 0x1) {
        fillHistogram(iPt, 0, iType); // PID 1 : TPC
        if (abs(beta-1) < 0.025) fillHistogram(iPt, 1, iType); // PID 1 : TPC + TOF
        if (e0/pt/TMath::CosH(eta) > 0.8 && e0/pt/TMath::CosH(eta) < 2) { // PID 2 : TPC + BEMC
            fillHistogram(iPt, 2, iType);
            if (nphi > 1 && neta > 1) fillHistogram(iPt, 4, iType); // PID 4 : TPC + BEMC + BSMD for MB

        }
    }
    if (isHTEvents >> 1 & 0x1 && e0/pt/TMath::CosH(eta) > 0.8 && e0/pt/TMath::CosH(eta) < 2) { // PID 3 TPC + BEMC + BSMC
        if (nphi > 1 && neta > 1) fillHistogram(iPt, 3, iType); // PID 3 : TPC + BEMC + BSMD for HT
        fillHistogram(iPt, 5, iType); // PID 5 : TPC + BEMC for HT
    }
    if (iType==3) {
        histoTofMass[iPt]->Fill(tofmass,weight);
        if (tofmass < 1 && tofmass > 0.86) histoNSigE[iPt][0]->Fill(nsige,weight); // pion
        else if (tofmass < 0.55 && tofmass > 0.4) histoNSigE[iPt][1]->Fill(nsige,weight); // kona
        else if (tofmass < 0.15 && tofmass > 0.12) histoNSigE[iPt][2]->Fill(nsige,weight); // proton
    }
}

