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
#include "StPicoNpeEventMaker/StPicoNpeEventMaker.h"
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
    hTriggerWt = new TH1I("hTriggerWt","hTriggerWt",30,0,30);
    hEvent = new TH1F("hEvent","hEvent",10,0,10);
    hZDCx = new TH1F("hZDCx","hZDCx",1000,0,100000);
    hZDCxWt = new TH1F("hZDCxWt","hZDCxWt",1000,0,100000);
    hCheckDoubleTrigger->Sumw2();
    hTrigger->Sumw2();
    hTriggerWt->Sumw2();
    hEvent->Sumw2();
    hZDCx->Sumw2();
    hZDCxWt->Sumw2();
    
    for (int i=0; i<2; i++) {
        hRefMult[i] = new TH1F(Form("hRefMult_%d",i),Form("hRefMult_%d",i),1000,0,1000);
        hRefMultWt[i] = new TH1F(Form("hRefMultWt_%d",i),Form("hRefMultWt_%d",i),1000,0,1000);
        hRefMult[i]->Sumw2();
        hRefMultWt[i]->Sumw2();
    }

    setHistogram(6,8,4,12);  // nPt, nPid, nType, nHisto

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
    hTriggerWt->Write();
    hCheckDoubleTrigger->Write();
    
    hRefMult[0]->Write();
    hRefMult[1]->Write();
    hRefMultWt[0]->Write();
    hRefMultWt[1]->Write();
    
    for (int l=0;l<12;l++) // Histograms
    for (int j=0;j<8;j++) // PID
    for (int i=1;i<6;i++) // PT
    for (int k=0;k<4;k++) // Particle:Type
         histo[i][j][k][l]->Write();

    for (int j=0;j<8;j++) // PID
        histo2d[j]->Write();
    
    
//    for (int i=1;i<6;i++) histoTofMass[i]->Write(); // tofmass
//    for (int i=1;i<6;i++) for (int j=0;j<3;j++) histoNSigE[i][j]->Write();
//    for (int j=2;j<6;j++) for (int i=1;i<6;i++) histo[i][j][2][2]->Write();
//    for (int j=0;j<6;j++) for (int i=1;i<6;i++) for (int k=0;k<2;k++) histo[i][j][k][3]->Write();
 
    mOutputFile->Close();
    
    return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Make()
{

    hEvent->Fill(0);
    readNextEvent();
    hEvent->Fill(1);

    if (!mPicoDstMaker)
    {
        LOG_WARN << " StPicoNpeAnaMaker - No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    hEvent->Fill(2);

    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    
    if (!picoDst)
    {
        LOG_WARN << "StPicoNpeAnaMaker - No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    hEvent->Fill(3);

    if(mPicoNpeEvent->runId() != picoDst->event()->runId() ||
       mPicoNpeEvent->eventId() != picoDst->event()->eventId())
    {
        LOG_ERROR <<" StPicoNpeAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoNpeEvent are not in sync."<<endm;
        exit(1);
    }
    hEvent->Fill(4);

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
    int checkDoubleTrigger = 0;
    int checkDoubleTrigger18 = 10;
    weight = 1;

    for (int i=0; i<4; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        if ( mPrescales->prescale(mPicoNpeEvent->runId(), i) > 100)cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        weight = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        hTriggerWt->Fill(i,weight);
        checkDoubleTrigger++;
        isHTEvents=1;
    }
    for (int i=19; i<20; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        if ( mPrescales->prescale(mPicoNpeEvent->runId(), i) > 100)cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        weight = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        hTriggerWt->Fill(i,weight);
        checkDoubleTrigger18++;
        isHTEvents=2;
    }

    hCheckDoubleTrigger->Fill(checkDoubleTrigger);
    hCheckDoubleTrigger->Fill(checkDoubleTrigger18);
    
    hEvent->Fill(5);
    hEvent->Fill(6,weight);

    hZDCx->Fill(mZDCx);
    hZDCxWt->Fill(mZDCx,weight);
    
    for (int i=0; i<2; i++) if (isHTEvents >> i & 0x1) {
        hRefMult[i]->Fill(mRefMult);
        hRefMultWt[i]->Fill(mRefMult,weight);
    }

    // hadrons & inclusive electron with StPicoTrack
    mPicoNpeEvent = new StPicoNpeEvent();
    std::vector<unsigned short> idxPicoTaggedEs;
    std::vector<unsigned short> idxPicoPartnerEs;

    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (!track) continue;
        if (isGoodTrack(track)) {
            setVariables(track);
            fillHistogram(2); // electron
            fillHistogram(3); // hadron
            
            
            if (StPicoNpeEventMaker::isElectron(track)) {
                idxPicoTaggedEs.push_back(iTrack);
                StElectronTrack electronTrack((StPicoTrack const *)track, iTrack);
                mPicoNpeEvent->addElectron(&electronTrack);
            }
            if (isPartnerElectron(track)) idxPicoPartnerEs.push_back(iTrack);

        }
        
    }
    
    mPicoNpeEvent->nElectrons(idxPicoTaggedEs.size());
    mPicoNpeEvent->nPartners(idxPicoPartnerEs.size());

    for (unsigned short ik = 0; ik < idxPicoTaggedEs.size(); ++ik)
    {
        
        StPicoTrack const * electron = picoDst->track(idxPicoTaggedEs[ik]);
        
        // make electron pairs
        for (unsigned short ip = 0; ip < idxPicoPartnerEs.size(); ++ip)
        {
            
            if (idxPicoTaggedEs[ik] == idxPicoPartnerEs[ip]) continue;
            
            StPicoTrack const * partner = picoDst->track(idxPicoPartnerEs[ip]);
            
            StElectronPair electronPair(electron, partner, idxPicoTaggedEs[ik], idxPicoPartnerEs[ip], pVtx, bField);
            
            
            if (!isGoodElectronPair(electronPair, electron->gPt())) continue;
            
            mPicoNpeEvent->addElectronPair(&electronPair);
            }
        } // .. end make electron pairs
    } // .. end of tagged e loop
    
    
    // Photonic Electron
    TClonesArray const * aElectronPair = mPicoNpeEvent->electronPairArray();
    for (int idx = 0; idx < aElectronPair->GetEntries(); ++idx)
    {
        // this is an example of how to get the ElectronPair pairs and their corresponsing tracks
        StElectronPair * epair = (StElectronPair*)aElectronPair->At(idx);
        if (isGoodPair(epair))
        {
            setVariables(epair);
            if (pairCharge == 0) {                  // US
                fillHistogram(0);
            }
            else {                                  // LS
                fillHistogram(1);
            }
            
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
    trk->nSigmaElectron() <= cutsAna::partnerNSigElectronHigh &&
    trk->nSigmaElectron() >= cutsAna::partnerNSigElectronLow
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
        nphieta = nphi+neta;
        adc0 = Emc->adc0();
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
    TString pid[10] = {"TpcMB","TpcTofMB","TpcBemcMB","TpcBemcBsmdHT","TpcBemcBsmdMB","TpcBemcHT"};
    TString type[10] = {"PhEUS","PhELS","IncE","Pion","Kaon","Proton"};
    TString histoname[12] = {"nSigE","nSigEAfterCut","dca","pairDca","nEta","nPhi","e0/p","zDist","phiDist","etaTowDist","phiTowDist","nphieta"};
    
    int binHisto[12] = {    289,    289,    100,    100,    10, 10, 100,    100,    100,    100,    100,    20};
    double minHisto[12] = { -13,    -13,   -0.1,    0,      0,  0,  0,      -20,    -0.1,   -0.1,   -0.1,   0};
    double maxHisto[12] = { 13,      13,    0.1,    0.5,    10, 10, 3,      20,     0.1,    0.1,    0.1,    20};
    
    
    for (int j=0; j<npid; j++) {
        histo2d[j] = new TH2F(Form("histo2D_%d",j),Form("histo2D_%d",j),100,0,20,100,0,500);
        histo2d[j]->Sumw2();

    }
    
    for (int i=0;i<nptbin;i++){
        histoTofMass[i] = new TH1F(Form("histoTofMass_%d",i),Form("histoTofMass_%d",i),1000,-0.5,2.5);
        histoTofMass[i]->Sumw2();
        
        for (int j=0;j<npid;j++) {
            histoNSigE[i][j] =  new TH1F(Form("histoNSigE_%d_%d",i,j),Form("histoNSigE_%d_%d",i,j),1301,-13,13);
            histoNSigE[i][j]->Sumw2();
        }
        for (int j=0;j<npid;j++)
            for (int k=0;k<ntype;k++)
                for (int l=0;l<nhisto;l++) {
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
                    histo[i][j][k][l]->Sumw2();
                }
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
void StPicoNpeAnaMaker::fillHistogram(int iType){

    int iPt = getPtBin(pt);
    if (!isBHTevent()){
        fillHistogram(iPt, 0, iType);
        if (isTof())                fillHistogram(iPt, 1, iType);
        if (isBemc())               fillHistogram(iPt, 2, iType);
        if (isBemc() && isBsmd())   fillHistogram(iPt, 3, iType);
    }
    if (isBHTevent()){
        fillHistogram(iPt, 4, iType);
        if (isTof())               fillHistogram(iPt, 5, iType);
        if (isBemc())              fillHistogram(iPt, 6, iType);
        if (isBemc() && isBsmd())  fillHistogram(iPt, 7, iType);
    }
}
//-------------------------------------------------------------------------------
void StPicoNpeAnaMaker::fillHistogram(int iPt, int iPid, int iType){
    histo[(const int)iPt][(const int)iPid][(const int)iType][0]->Fill(nsige,weight);
    if (iType < 2 && nsige > 0) fillHistogram(iPt, iPid, iType, 0);
    else if (iType==2) {
        float pidCutLw[2][6];
        float pidCutHi[2][6];
        pidCutLw[0]={0, -1.2, -1.2, -1.0, -1.0, 0};
        pidCutHi[0]={0, 1.8, 2.5, 3.0, 3.0, 0};
        pidCutLw[1]={0, -1.5, -1.4, -1.5, -1.1, 0};
        pidCutHi[1]={0, 1.8, 2.5, 3.0, 3.0, 0};
//        if (iPid == 2 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) fillHistogram(iPt, iPid, iType, 0);
//        if (iPid == 3 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) fillHistogram(iPt, iPid, iType, 0);
//        if (iPid == 4 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) fillHistogram(iPt, iPid, iType, 0);
//        if (iPid == 5 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) fillHistogram(iPt, iPid, iType, 0);
        if (iPid%4 < 3 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) fillHistogram(iPt, iPid, iType, 0);
        if (iPid%4 ==3 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) fillHistogram(iPt, iPid, iType, 0);
        histo2d[iPid]->Fill(pt*TMath::CosH(eta),adc0,weight);
    }
    else if (iType==3 && fabs(nsigpion) < 2) fillHistogram(iPt, iPid, iType, 0);
}
//-------------------------------------------------------------------------------
void StPicoNpeAnaMaker::fillHistogram(int iPt, int iPid, int iType, int dummy){
    histo[(const int)iPt][(const int)iPid][(const int)iType][1]->Fill(nsige,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][2]->Fill(dca,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][3]->Fill(pairMass,weight);
    
    // PID QA
    histo[(const int)iPt][(const int)iPid][(const int)iType][4]->Fill(neta,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][5]->Fill(nphi,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][6]->Fill(e0/pt/TMath::CosH(eta),weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][7]->Fill(zDist,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][8]->Fill(phiDist,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][9]->Fill(etaTowDist,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][10]->Fill(phiTowDist,weight);
    histo[(const int)iPt][(const int)iPid][(const int)iType][11]->Fill(nphieta,weight);

}

//-----------------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isBemc(){
    if (e0/pt/TMath::CosH(eta) > 0.8 && e0/pt/TMath::CosH(eta) < 2) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBsmd(){
    if (nphi > 1 && neta > 1) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isTof(){
    if (fabs(beta-1) < 0.025) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBHTevent(){
    if (isHTEvents >> 1 & 0x1 ) return true;
    else return false;
}


