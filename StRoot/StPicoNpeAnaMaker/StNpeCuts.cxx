
#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

#include "StNpeCuts.h"

#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StBTofUtil/tofPathLength.hh"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"


ClassImp(StNpeCuts)

// _________________________________________________________
StNpeCuts::StNpeCuts() : StPicoCutsBase("NpeCutsBase"), mPicoDst2(NULL),
mElectronPairDcaDaughtersMax(std::numeric_limits<float>::max()),
mElectronPairDecayLengthMin(std::numeric_limits<float>::min()), mElectronPairDecayLengthMax(std::numeric_limits<float>::max()),
mElectronPairCosThetaMin(std::numeric_limits<float>::min()),
mElectronPairMassMin(std::numeric_limits<float>::min()), mElectronPairMassMax(std::numeric_limits<float>::max()),
mElectronNHitsFitMax(std::numeric_limits<int>::max()),
mElectronNHitsdEdxMax(std::numeric_limits<int>::max()),
mElectronBsmdNEta(std::numeric_limits<int>::min()),
mElectronPtMin(std::numeric_limits<float>::min()),
mElectronBsmdNPhi(std::numeric_limits<int>::min()),
mElectronTofBeta(std::numeric_limits<float>::max()),
mElectronPtMax(std::numeric_limits<float>::max()),
mElectronEtaMin(std::numeric_limits<float>::min()),
mElectronRequireHFT(false),
mElectronEtaMax(std::numeric_limits<float>::max()),
mElectronDca(std::numeric_limits<float>::max()),
mElectronTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mElectronTPCNSigmaElectronMax(std::numeric_limits<float>::max()),
mElectronBemcEoverPMin(std::numeric_limits<float>::min()),
mElectronBemcEoverPMax(std::numeric_limits<float>::max()),
mElectronBemcPhiDistMax(std::numeric_limits<float>::max()),
mElectronBemcZDistMax(std::numeric_limits<float>::max()),
mElectronBemcAssDistMax(std::numeric_limits<float>::max()),
mPartnerElectronNHitsFitMax(std::numeric_limits<int>::min()),
mPartnerElectronNHitsdEdxMax(std::numeric_limits<int>::min()),
mPartnerElectronPtMin(std::numeric_limits<float>::min()),
mPartnerElectronPtMax(std::numeric_limits<float>::max()),
mPartnerElectronEtaMin(std::numeric_limits<float>::min()),
mPartnerElectronEtaMax(std::numeric_limits<float>::max()),
mPartnerElectronRequireHFT(false),
mPartnerTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mPartnerTPCNSigmaElectronMax(std::numeric_limits<float>::max()),
mElectronBemcPid(false),mElectronBsmdPid(false),mElectronTofPid(false){
    
    // -- default constructor
}

// _________________________________________________________
StNpeCuts::StNpeCuts(const Char_t *name) : StPicoCutsBase(name), mPicoDst2(NULL),
mElectronPairDcaDaughtersMax(std::numeric_limits<float>::max()),
mElectronPairDecayLengthMin(std::numeric_limits<float>::min()), mElectronPairDecayLengthMax(std::numeric_limits<float>::max()),
mElectronPairCosThetaMin(std::numeric_limits<float>::min()),
mElectronPairMassMin(std::numeric_limits<float>::min()), mElectronPairMassMax(std::numeric_limits<float>::max()),
mElectronNHitsFitMax(std::numeric_limits<int>::max()),
mElectronNHitsdEdxMax(std::numeric_limits<int>::max()),
mElectronBsmdNEta(std::numeric_limits<int>::min()),
mElectronPtMin(std::numeric_limits<float>::min()),
mElectronBsmdNPhi(std::numeric_limits<int>::min()),
mElectronTofBeta(std::numeric_limits<float>::max()),
mElectronPtMax(std::numeric_limits<float>::max()),
mElectronEtaMin(std::numeric_limits<float>::min()),
mElectronRequireHFT(false),
mElectronEtaMax(std::numeric_limits<float>::max()),
mElectronDca(std::numeric_limits<float>::max()),
mElectronTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mElectronTPCNSigmaElectronMax(std::numeric_limits<float>::max()),
mElectronBemcEoverPMin(std::numeric_limits<float>::min()),
mElectronBemcEoverPMax(std::numeric_limits<float>::max()),
mElectronBemcPhiDistMax(std::numeric_limits<float>::max()),
mElectronBemcZDistMax(std::numeric_limits<float>::max()),
mElectronBemcAssDistMax(std::numeric_limits<float>::max()),
mPartnerElectronNHitsFitMax(std::numeric_limits<int>::min()),
mPartnerElectronNHitsdEdxMax(std::numeric_limits<int>::min()),
mPartnerElectronPtMin(std::numeric_limits<float>::min()),
mPartnerElectronPtMax(std::numeric_limits<float>::max()),
mPartnerElectronEtaMin(std::numeric_limits<float>::min()),
mPartnerElectronEtaMax(std::numeric_limits<float>::max()),
mPartnerElectronRequireHFT(false),
mPartnerTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mPartnerTPCNSigmaElectronMax(std::numeric_limits<float>::max()),
mElectronBemcPid(false),mElectronBsmdPid(false),mElectronTofPid(false){

    
    // -- constructor
}

// _________________________________________________________
StNpeCuts::~StNpeCuts() {
    // destructor
    
}


// _________________________________________________________
bool StNpeCuts::isGoodElectronPair(StElectronPair const* epair) const {
    // -- check for good electron pairs
    StPicoTrack const* electron = mPicoDst2->track(epair->electronIdx());
    StPicoTrack const* partner = mPicoDst2->track(epair->partnerIdx());
    
    return
    isGoodTaggedElectron(electron) &&
    isGoodPartnerElectron(partner) &&
    epair->pairMass() > mElectronPairMassMin && epair->pairMass() < mElectronPairMassMax &&
    epair->pairDca() < mElectronPairDcaDaughtersMax ;
}
// _________________________________________________________
bool StNpeCuts::isGoodInclusiveElectron(StPicoTrack const *trk) const {
    // -- check for good tagged electron for electron pairs
    bool taggedElectronCut =
    trk->nHitsFit() >= mElectronNHitsFitMax &&
    trk->nHitsDedx() >= mElectronNHitsdEdxMax &&
    trk->gPt() >= mElectronPtMin && trk->gPt() < mElectronPtMax &&
    getEta(trk) > mElectronEtaMin && getEta(trk) < mElectronEtaMax &&
    fabs(getDca(trk)) < mElectronDca &&
    (!mElectronRequireHFT || trk->isHFTTrack());
    
    return taggedElectronCut
    && isTPCElectron(trk, -13, 13)
    ;
}
// _________________________________________________________
bool StNpeCuts::isGoodTaggedElectron(StPicoTrack const *trk) const {
    // -- check for good tagged electron for electron pairs
    return isGoodInclusiveElectron(trk)
    && isTPCElectron(trk, mElectronTPCNSigmaElectronMin, mElectronTPCNSigmaElectronMax)
    && isBEMCElectron(trk)
    && isBSMDElectron(trk)
    && isTOFElectron(trk)
    ;
}
// _________________________________________________________
bool StNpeCuts::isGoodPartnerElectron(StPicoTrack const *trk) const {
    // -- check for good partner electron for electron pairs

    bool partnerElectronCut =
    trk->nHitsFit() >= mPartnerElectronNHitsFitMax &&
    trk->gPt() >= mPartnerElectronPtMin && trk->gPt() < mPartnerElectronPtMax &&
    getEta(trk) > mPartnerElectronEtaMin && getEta(trk) < mPartnerElectronEtaMax &&
    (!mPartnerElectronRequireHFT || trk->isHFTTrack());
    
    return partnerElectronCut
    && isTPCElectron(trk, mPartnerTPCNSigmaElectronMin, mPartnerTPCNSigmaElectronMax)
    && isTOFElectron(trk)
;
}
// _________________________________________________________
bool StNpeCuts::isTPCElectron(StPicoTrack const *trk, float min, float max) const {
    // -- check for good TPC electrons
    float nSigma = trk->nSigmaElectron();
    
    return
    nSigma > min && nSigma < max;
}
// _________________________________________________________
bool StNpeCuts::isTOFElectron(StPicoTrack const *trk) const {
    // -- check for good TOF electrons
    if (!mElectronTofPid || trk->gPt() >= 2.) return true;
    StPicoBTofPidTraits *tofPid = hasTofPid(trk);
    
    float beta;
    if (tofPid) {
        beta = tofPid->btofBeta();
        if (beta < 1e-4) {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();
            StPhysicalHelixD helix = trk->helix();
            float pathLength = tofPathLength(&getpVtx(), &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            beta = (tof > 0) ? pathLength / (tof * (C_C_LIGHT / 1.e9)) : std::numeric_limits<float>::quiet_NaN();
        }
    }
    else beta=999;
    
    if (fabs(1/beta -1) > mElectronTofBeta) return true;
    else return false;
}

// _________________________________________________________
bool StNpeCuts::isBEMCElectron(StPicoTrack const *trk) const {
    // -- check for good BEMC electrons
    if (!mElectronBemcPid || trk->gPt() < 2.) return true;
    if (trk->emcPidTraitsIndex() < 0) return false;
    StPicoEmcPidTraits * Emc =  mPicoDst2->emcPidTraits(trk->emcPidTraitsIndex());
    float eoverp = Emc->e0()/trk->gPt()/TMath::CosH(getEta(trk));
    float phiDist = Emc->phiDist();
    float zDist = Emc->zDist();
    
    return eoverp > mElectronBemcEoverPMin && eoverp < mElectronBemcEoverPMax &&
    phiDist < mElectronBemcPhiDistMax && zDist < mElectronBemcZDistMax &&
    TMath::Sqrt(phiDist*phiDist + zDist*zDist) < mElectronBemcAssDistMax
    ;
}

bool StNpeCuts::isBSMDElectron(StPicoTrack const *trk) const {
    // -- check for good BSMD electrons
    if (!mElectronBsmdPid) return true;
    if (trk->emcPidTraitsIndex() < 0) return false;
    StPicoEmcPidTraits * Emc =  mPicoDst2->emcPidTraits(trk->emcPidTraitsIndex());
    int nphi = Emc->nPhi();
    int neta = Emc->nEta();
    
    return neta > mElectronBsmdNEta && nphi > mElectronBsmdNPhi ;
}
     
// _________________________________________________________
float StNpeCuts::getEta(StPicoTrack const *trk) const {
    return trk->gMom(getpVtx(), mPicoDst2->event()->bField()).pseudoRapidity();
}
     
// _________________________________________________________
float StNpeCuts::getDca(StPicoTrack const *trk) const {
    return trk->dcaGeometry().helix().curvatureSignedDistance(getpVtx().x(),getpVtx().y());
}
     
// _________________________________________________________
StThreeVectorF StNpeCuts::getpVtx() const {
    return mPicoDst2->event()->primaryVertex();
}


