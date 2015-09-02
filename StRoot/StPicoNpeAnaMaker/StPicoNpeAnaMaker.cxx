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
    hCheckDoubleTrigger->Sumw2();
    hTrigger->Sumw2();
    hTriggerWt->Sumw2();
    hEvent->Sumw2();
    
    for (int i=0; i<nntrigger; i++) {
        hZDCx[i] = new TH1F(Form("hZDCx_%d",i),Form("hZDCx_%d",i),1000,0,100000);
        hZDCxWt[i] = new TH1F(Form("hZDCxWt_%d",i),Form("hZDCxWt_%d",i),1000,0,100000);
        hRefMult[i] = new TH1F(Form("hRefMult_%d",i),Form("hRefMult_%d",i),1000,0,1000);
        hRefMultWt[i] = new TH1F(Form("hRefMultWt_%d",i),Form("hRefMultWt_%d",i),1000,0,1000);
        hTriggerCheck[i] = new  TH1I(Form("hTriggerCheck_%d",i),Form("hTriggerCheck_%d",i),30,0,30);
        hTriggerCheckWt[i] = new  TH1I(Form("hTriggerCheckWt_%d",i),Form("hTriggerCheckWt_%d",i),30,0,30);
        hRefMult[i]->Sumw2();
        hRefMultWt[i]->Sumw2();
        hZDCx[i]->Sumw2();
        hZDCxWt[i]->Sumw2();
    }

    setHistogram();  // nPt, nPid, nType, nHisto, nTrigger

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
    hTrigger->Write();
    hTriggerWt->Write();
    hCheckDoubleTrigger->Write();
    
    for (int i=0;i<nntrigger;i++){
        hZDCx[i]->Write();
        hZDCxWt[i]->Write();
        hRefMult[i]->Write();
        hRefMultWt[i]->Write();
        hTriggerCheck[i]->Write();
        hTriggerCheckWt[i]->Write();
    }
    
    for (int l=0;l<nnhisto;l++) // Histograms
    for (int j=0;j<nnpid;j++) // PID
    for (int i=1;i<nnpt;i++) // PT
    for (int k=0;k<nntype;k++) // Particle:Type
    for (int m=0;m<nntrigger;m++) // Trigger
         histo[i][j][k][l][m]->Write();
    
    for (int j=0;j<nnpid;j++) {// PID
        for (int m=0;m<nntrigger;m++) { // Trigger
            histo2d[j][m]->Write();
      //      histo2dDcaPt[j][m]->Write();
        }
    }
    for (int i=1;i<nnpt;i++) // PT
    for (int j=0;j<6;j++)
        histoPureE[i][j]->Write();
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
        LOG_ERROR <<" StPicoNpeAnaMaker - mPicoNpeEvent->runId() = " << mPicoNpeEvent->runId() <<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - picoDst->event()->runId() = " << picoDst->event()->runId() <<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - mPicoNpeEvent->eventId() = " << mPicoNpeEvent->eventId() <<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - picoDst->event()->eventId() = " << picoDst->event()->eventId() <<endm;
        
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
    weight[0] = 1;
    weight[1] = 1;
    weight[2] = 1;
    weight[3] = 1;
    


    for (int i=0; i<4; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        if ( mPrescales->prescale(mPicoNpeEvent->runId(), i) > 100)cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        weight[0] = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        hTriggerWt->Fill(i,weight[0]);
        checkDoubleTrigger++;
        isHTEvents+=1;
        hEvent->Fill(6,weight[0]);
    }
    for (int i=19; i<21; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        if ( mPrescales->prescale(mPicoNpeEvent->runId(), i) > 100)cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        weight[1] = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        hTriggerWt->Fill(i,weight[1]);
        checkDoubleTrigger18++;
        isHTEvents+=2;
        hEvent->Fill(7,weight[1]);
    }
    for (int i=21; i<23; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        if ( mPrescales->prescale(mPicoNpeEvent->runId(), i) > 100)cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        weight[2] = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        hTriggerWt->Fill(i,weight[2]);
        checkDoubleTrigger18++;
        isHTEvents+=4;
        hEvent->Fill(8,weight[2]);
    }
    
    for (int i=23; i<25; i++) if (picoDst->event()->triggerWord() >> i & 0x1) {
        if ( mPrescales->prescale(mPicoNpeEvent->runId(), i) > 100)cout << "Prescale (" << mPicoNpeEvent->runId() << ", " << i << ", " << mPicoNpeEvent->eventId() << ") : " << mPrescales->prescale(mPicoNpeEvent->runId(), i) << endl;
        weight[3] = mPrescales->prescale(mPicoNpeEvent->runId(), i);
        hTrigger->Fill(i);
        hTriggerWt->Fill(i,weight[3]);
        checkDoubleTrigger18++;
        isHTEvents+=8;
        hEvent->Fill(9,weight[3]);
    }
    hEvent->Fill(5);

    hCheckDoubleTrigger->Fill(checkDoubleTrigger);
    hCheckDoubleTrigger->Fill(checkDoubleTrigger18);
    

    
    for (int i=0; i<nntrigger; i++) if (isHTEvents >> i & 0x1) {
        hZDCx[i]->Fill(mZDCx);
        hZDCxWt[i]->Fill(mZDCx,weight[i]);
        hRefMult[i]->Fill(mRefMult);
        hRefMultWt[i]->Fill(mRefMult,weight[i]);
        for (int j=0; j<30; j++) if (picoDst->event()->triggerWord() >> j & 0x1) {
            hTriggerCheck[i]->Fill(j);
            hTriggerCheckWt[i]->Fill(j,weight[i]);
        }
    }
    
    // hadrons & inclusive electron with StPicoTrack
    std::vector<unsigned short> idxPicoTaggedEs;
    std::vector<unsigned short> idxPicoPartnerEs;
    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (!track) continue;
        if (isGoodTrack(track)) {
            setVariables(track);
            if (track->isHFTTrack()){       // HFT
                fillHistogram(2);               // electron candidates
                fillHistogram(3);               // hadron
            }
  //          else {
  //              fillHistogram(8);               // electron candidates
  //              fillHistogram(9);               // hadron
  //          }
            
        }
  //      if (isGoodTagged(track))  idxPicoTaggedEs.push_back(iTrack);
  //      if (isGoodPartner(track)) idxPicoPartnerEs.push_back(iTrack);
        
        
    }
    /*
    for (unsigned short ik = 0; ik < idxPicoTaggedEs.size(); ++ik)
    {
        StPicoTrack const * electron = picoDst->track(idxPicoTaggedEs[ik]);
        // make electron pairs
        for (unsigned short ip = 0; ip < idxPicoPartnerEs.size(); ++ip)
        {
            if (idxPicoTaggedEs[ik] == idxPicoPartnerEs[ip]) continue;
            StPicoTrack const * partner = picoDst->track(idxPicoPartnerEs[ip]);
            StElectronPair * epair =  new StElectronPair(electron, partner, idxPicoTaggedEs[ik], idxPicoPartnerEs[ip], bField);
            if (isGoodPair(epair))
            {
                setVariables(epair);

                if (electron->isHFTTrack()){             // HFT
                    if (pairCharge == 0) fillHistogram(4); // US
                    else fillHistogram(5);                 // LS
                }
                else {                               // non HFT
                    if (pairCharge == 0) fillHistogram(6); // US
                    else fillHistogram(7);                 // LS
                    
                    // partner electron nsige variation
                    if (partner_nsige > -2){
                        if (pairCharge == 0) fillHistogram(10); // US
                        else fillHistogram(11);                 // LS
                        
                        if (partner_nsige > -1){
                            if (pairCharge == 0) fillHistogram(12); // US
                            else fillHistogram(13);                 // LS
                            
                            if (partner_nsige > 0){
                                if (pairCharge == 0) fillHistogram(14); // US
                                else fillHistogram(15);                 // LS
                                if (partner_nsige > 1){
                                    if (pairCharge == 0) fillHistogram(16); // US
                                    else fillHistogram(17);                 // LS
                                    
                                }

                            }
                        }
                    }
                }
                

                //cout << "1 " << pairMass << " " << epair->pairMass() << " " <<pairDca << " " << pt << " " << eta << " " << dca << " " << nsige << " " << pairPositionX << " " << pairPositionY << " " << pairPositionZ << " " << electron->nHitsFit() << " " << electron->nHitsDedx() <<  " " << partner->nHitsFit() << " " << partner->nHitsDedx() << " " << partner->gPt() << endl;
            }
            delete epair;
        } // .. end make electron pairs
    } // .. end of tagged e loop
    */
    // Photonic Electron
    TClonesArray const * aElectronPair = mPicoNpeEvent->electronPairArray();
    for (int idx = 0; idx < aElectronPair->GetEntries(); ++idx)
    {
        // this is an example of how to get the ElectronPair pairs and their corresponsing tracks
        StElectronPair * epair = (StElectronPair*)aElectronPair->At(idx);
        if (isGoodPair(epair))
        {
            setVariables(epair);
            if (pairCharge == 0) {
                fillHistogram(0); // US
                if (pairMass < 0.1) fillHistogram(4); // US
                if (pairMass < 0.01) fillHistogram(6); // US
                fillHistogram("PureE");
            }
            else {
                fillHistogram(1);                 // LS
                if (pairMass < 0.1) fillHistogram(5); // US
                if (pairMass < 0.01) fillHistogram(7); // US

            }
            //cout << "0 " << pairMass << " " << epair->pairMass() << " " << pairDca << " " << pt << " " << eta  << " " << dca << " " << nsige << " " << pairPositionX << " " << pairPositionY << " " << pairPositionZ << endl;
        }
    }
    
    idxPicoTaggedEs.clear();
    idxPicoPartnerEs.clear();
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
    epair->pairDca() < cutsAna::pairDca &&
    fabs(epair->positionX()) < 200. && fabs(epair->positionY()) < 200. && fabs(epair->positionZ()) < 200.
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
    StPhysicalHelixD eHelix = trk->dcaGeometry().helix();
    
    return
  //  (!cutsAna::trackRequireHFT || trk->isHFTTrack()) &&
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
    return true
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
    dcaCharge = eHelix.geometricSignedDistance(pVtx.x(),pVtx.y());
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
    cout << electronHelix.momentum(bField).perp() << " " << electronMomAtDca.perp() << " " << pt <<  " " << electron->pMom().perp() << endl;
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
    pairPositionZ = epair->positionZ();
    pairConvRadious = TMath::Sqrt((pairPositionX+0.2383) * (pairPositionX+0.2383) + (pairPositionY+0.1734) * (pairPositionY+0.1734));
    
    partner_pt = partner->gPt();
    partner_nsige = partner->nSigmaElectron();
}
//-----------------------------------------------------------------------------
void StPicoNpeAnaMaker::setHistogram()
{
    
    double ptbin[nnpt+1] = {0, 1.5, 1.8, 2.5, 4.0, 6.5, 10.};
    TString pid[nnpid+20] = {
        "Tpc", "TpcBemc3Bsmd", "TpcNoBemc3NoBsmd0", "TpcNoBemc3NoBsmd1", "TpcDalitz", "TpcConversion", "TpcBemc3BsmdDalitz", "TpcBemc3BsmdConversion"};
      //  "Tpc","TpcTof","TpcBemc","TpcBemc2","TpcBsmd","TpcBemcBsmd","TpcBemc2Bsmd","TpcBemc3Bsmd",
      //  "TpcBemc4","TpcBemc5","TpcBemc4Bsmd","TpcBemc5Bsmd"};
    TString type[nntype+20] = {
        "PhEUS","PhELS","IncE","Pion","PhEUS_Mass100","PhELS_Mass100","PhEUS_Mass10","PhELS_Mass10"};
        //"PhEUS","PhELS","IncE","Pion","RecoHFTPhEUS","RecoHFTPhELS","RecoNonHFTPhEUS","RecoNonHFTPhELS","IncENonHFT","PionNonHFT"};
    TString histoname[nnhisto] = {"nSigE","nSigEAfterCut","dca","pairMass","nEta","nPhi","e0/p","zDist","phiDist","etaTowDist","phiTowDist","nphieta","e/p","ConvRadious","pairDca","dcaCharge"};
    TString trigger[nntrigger] = {"MB", "BHT1", "BHT2", "BHT3"};
    
    int binHisto[nnhisto] = {    289,    289,    100,    100,    10, 10, 200,    100,    100,    100,    100,    20  ,200    ,500    ,200   ,100};
    double minHisto[nnhisto] = { -13,    -13,   -0.1,    0,      0,  0,  0,      -20,    -0.1,   -0.1,   -0.1,   0   ,0      ,0      ,0     ,-0.1};
    double maxHisto[nnhisto] = { 13,      13,    0.1,    0.2,    10, 10, 6,      20,     0.1,    0.1,    0.1,    20  ,6      ,50     ,2     ,0.1};
    
    
    for (int j=0; j<nnpid; j++)
        for (int m=0;m<nntrigger;m++) {
            histo2d[j][m] = new TH2F(Form("histo2D_%d_%d",j,m),Form("histo2D_%d_%d",j,m),100,0,10,200,0,1000);
       //     histo2dDcaPt[j][m] = new TH2F(Form("histo2DDcaPt_%d_%d",j,m),Form("histo2DDcaPt_%d_%d",j,m),100,0,10,100,-0.1,0.1);
            histo2d[j][m]->Sumw2();
       //     histo2dDcaPt[j][m]->Sumw2();
            
    }
    
    for (int i=0;i<nnpt;i++){
        histoPureE[i][0] = new TH1F(Form("histoPureE_%d_%d",i,0),Form("histoPureE_%d_%d",i,0),100,0,10);
        histoPureE[i][1] = new TH1F(Form("histoPureE_%d_%d",i,1),Form("histoPureE_%d_%d",i,1),100,0,0.05);
        histoPureE[i][2] = new TH1F(Form("histoPureE_%d_%d",i,2),Form("histoPureE_%d_%d",i,2),100,0,0.01);
        histoPureE[i][3] = new TH1F(Form("histoPureE_%d_%d",i,3),Form("histoPureE_%d_%d",i,3),100,0,1);
        histoPureE[i][4] = new TH1F(Form("histoPureE_%d_%d",i,4),Form("histoPureE_%d_%d",i,4),100,0,1);
        histoPureE[i][5] = new TH1F(Form("histoPureE_%d_%d",i,5),Form("histoPureE_%d_%d",i,5),100,-0.1,0.1);
        for (int j=0;j<nnpid;j++)
            for (int k=0;k<nntype;k++)
                for (int l=0;l<nnhisto;l++)
                    for (int m=0;m<nntrigger;m++) {
                        histo[i][j][k][l][m] = new TH1F(
                                                        Form("histo_%d_%d_%d_%d_%d", i,j,k,l,m),
                                                        Form("histo_pT%.1f_%.1f_%s_%s_%s_%s",
                                                             ptbin[i],
                                                             ptbin[i+1],
                                                             pid[j].Data(),
                                                             type[k].Data(),
                                                             histoname[l].Data(),
                                                             trigger[m].Data()
                                                             ),
                                                        binHisto[l],
                                                        minHisto[l],
                                                        maxHisto[l]
                                                        );
                        histo[i][j][k][l][m]->Sumw2();
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
    fillHistogram(iPt, 0, iType);
    if (isBemc3() && isBsmd())  fillHistogram(iPt, 1, iType);
    if (!isBemc3() && !isBsmd())fillHistogram(iPt, 2, iType);
    if (!isBemc3() || !isBsmd())fillHistogram(iPt, 3, iType);
    if (pairConvRadious < 2. )                         fillHistogram(iPt, 4, iType);
    if (pairConvRadious > 2. && pairConvRadious < 4.)  fillHistogram(iPt, 5, iType);
    if (isBemc3() && isBsmd() && pairConvRadious < 2. )                         fillHistogram(iPt, 6, iType);
    if (isBemc3() && isBsmd() && pairConvRadious > 2. && pairConvRadious < 4.)  fillHistogram(iPt, 7, iType);
    
    /*
     if (isTof())                fillHistogram(iPt, 1, iType);
     if (isBemc())               fillHistogram(iPt, 2, iType);
     if (isBemc2())              fillHistogram(iPt, 3, iType);
     if (isBsmd())               fillHistogram(iPt, 4, iType);
     if (isBemc() && isBsmd())   fillHistogram(iPt, 5, iType);
     if (isBemc2() && isBsmd())  fillHistogram(iPt, 6, iType);
     if (isBemc3() && isBsmd())  fillHistogram(iPt, 7, iType);
     if (isBemc4())              fillHistogram(iPt, 8, iType);
     if (isBemc5())              fillHistogram(iPt, 9, iType);
     if (isBemc4() && isBsmd())  fillHistogram(iPt, 10, iType);
     if (isBemc5() && isBsmd())  fillHistogram(iPt, 11, iType);
     */
}
//-------------------------------------------------------------------------------
void StPicoNpeAnaMaker::fillHistogram(int iPt, int iPid, int iType){
    for (int i=0; i<nntrigger; i++) {
        if (isHTEvents >> i & 0x1 ) {
            histo[(const int)iPt][(const int)iPid][(const int)iType][0][i]->Fill(nsige,weight[i]);
            if ((
                 iType == 0 || iType == 1 ||    // PhE from NPE tree
                 iType == 4 || iType == 5 ||    //  " and mass < 0.1 //PhE w/ HFT recon.
                 iType == 6 || iType == 7 ||    //  " and mass < 0.01 //PhE w/o HFT recon.
                 iType == 10 || iType == 11 ||  // PhE w/ HFT recon. nSigE_Part > -2
                 iType == 12 || iType == 13 ||  // PhE w/ HFT recon. nSigE_Part > -1
                 iType == 14 || iType == 15 ||  // PhE w/ HFT recon. nSigE_Part >  0
                 iType == 16 || iType == 17     // PhE w/ HFT recon. nSigE_Part >  1
                 ) &&
                nsige > -1)
            {
                fillHistogram(iPt, iPid, iType, i);

          //      if (iType==0) histo2dDcaPt[iPid][i]->Fill(pt,dca,weight[i]);
          //      if (iType==1) histo2dDcaPt[iPid][i]->Fill(pt,dca,-1*weight[i]);
            }
            if (iType==2 || iType==8) {
                float pidCutLw[2][6];
                float pidCutHi[2][6];
                pidCutLw[0]={0, -1.2, -1.2, -1.0, -1.0, -1.0};
                pidCutHi[0]={0, 1.8, 2.5, 3.0, 3.0, 3.0};
                pidCutLw[1]={0, -1.5, -1.4, -1.5, -1.1, 0};
                pidCutHi[1]={0, 1.8, 2.5, 3.0, 3.0, 0};
                //        if (iPid == 2 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) fillHistogram(iPt, iPid, iType, 0);
                //        if (iPid == 3 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) fillHistogram(iPt, iPid, iType, 0);
                //        if (iPid == 4 && nsige > pidCutLw[1][iPt] && nsige < pidCutHi[1][iPt]) fillHistogram(iPt, iPid, iType, 0);
                //        if (iPid == 5 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) fillHistogram(iPt, iPid, iType, 0);
                //        if (iPid%4 < 3 && nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) fillHistogram(iPt, iPid, iType, 0);
                if (nsige > pidCutLw[0][iPt] && nsige < pidCutHi[0][iPt]) {
                    fillHistogram(iPt, iPid, iType, i);
                    histo2d[iPid][i]->Fill(pt*TMath::CosH(eta),adc0,weight[i]);
                }
            }
            if ((iType==3 || iType==9) && fabs(nsigpion) < 2) fillHistogram(iPt, iPid, iType, i);
        }
    }
}
//-------------------------------------------------------------------------------
void StPicoNpeAnaMaker::fillHistogram(int iPt, int iPid, int iType, int iTrigger){
    histo[(const int)iPt][(const int)iPid][(const int)iType][1][iTrigger]->Fill(nsige,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][2][iTrigger]->Fill(dca,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][3][iTrigger]->Fill(pairMass,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][15][iTrigger]->Fill(dcaCharge,weight[iTrigger]);

    // PID QA
    histo[(const int)iPt][(const int)iPid][(const int)iType][4][iTrigger]->Fill(neta,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][5][iTrigger]->Fill(nphi,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][6][iTrigger]->Fill(eoverp,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][7][iTrigger]->Fill(zDist,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][8][iTrigger]->Fill(phiDist,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][9][iTrigger]->Fill(etaTowDist,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][10][iTrigger]->Fill(phiTowDist,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][11][iTrigger]->Fill(nphieta,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][12][iTrigger]->Fill(e/pt/TMath::CosH(eta),weight[iTrigger]);
    
    // pair QA
    histo[(const int)iPt][(const int)iPid][(const int)iType][13][iTrigger]->Fill(pairConvRadious,weight[iTrigger]);
    histo[(const int)iPt][(const int)iPid][(const int)iType][14][iTrigger]->Fill(pairDca,weight[iTrigger]);
    
}
void StPicoNpeAnaMaker::fillHistogram(TString check){
    histoPureE[getPtBin(pt)][1]->Fill(pairMass);
    if (check ==  "PureE" && pairMass < 0.01 && nsige > -1){
        histoPureE[getPtBin(pt)][0]->Fill(pairConvRadious);
        histoPureE[getPtBin(pt)][2]->Fill(pairMass);
        histoPureE[getPtBin(pt)][3]->Fill(pairAnglePhi);
        histoPureE[getPtBin(pt)][4]->Fill(pairDca);
        histoPureE[getPtBin(pt)][5]->Fill(dca);
    }
}


//-----------------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isBemc(){
    if (eoverp > cutsAna::emcEoverPLow && eoverp < cutsAna::emcEoverPHigh) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBemc2(){
    if (eoverp > cutsAna::emcEoverPLow2 && eoverp < cutsAna::emcEoverPHigh2) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBemc3(){
    if (isBemc2() && fabs(zDist) < cutsAna::emcZDist && fabs(phiDist) < cutsAna::emcPhiDist && sqrt(etaTowDist*etaTowDist + phiTowDist*phiTowDist) < cutsAna::emcAssDist ) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBemc4(){
    if (e/pt/TMath::CosH(eta) > cutsAna::emcEoverPLow4 && e/pt/TMath::CosH(eta) < cutsAna::emcEoverPHigh4) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBemc5(){
    if (isBemc4() && fabs(zDist) < cutsAna::emcZDist && fabs(phiDist) < cutsAna::emcPhiDist && sqrt(etaTowDist*etaTowDist + phiTowDist*phiTowDist) < cutsAna::emcAssDist ) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isBsmd(){
    if (nphi > cutsAna::emcNEta && neta > cutsAna::emcNPhi) return true;
    else return false;
}
bool StPicoNpeAnaMaker::isTof(){
    if (fabs(beta-1) < cutsAna::tofBeta) return true;
    else return false;
}


