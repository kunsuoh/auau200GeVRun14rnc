
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
mPartnerTPCNSigmaElectronMax(std::numeric_limits<float>::max()) {
    
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
mPartnerTPCNSigmaElectronMax(std::numeric_limits<float>::max()) {
    
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
bool StNpeCuts::isGoodTaggedElectron(StPicoTrack const *trk) const {
    // -- check for good tagged electron for electron pairs
    bool taggedElectronCut =
    trk->nHitsFit() >= mElectronNHitsFitMax &&
    trk->nHitsDedx() >= mElectronNHitsdEdxMax &&
    trk->gPt() >= mElectronPtMin && trk->gPt() < mElectronPtMax &&
    getEta(trk) > mElectronEtaMin && getEta(trk) < mElectronEtaMax &&
    getDca(trk) < mElectronDca &&
    (!mElectronRequireHFT || trk->isHFTTrack());
    
    return taggedElectronCut
//    && isTPCElectron(trk, mElectronTPCNSigmaElectronMin, mElectronTPCNSigmaElectronMax)
//    && isBEMCElectron(trk)
//    && isBSMDElectron(trk)
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
    
    return partnerElectronCut && isTPCElectron(trk, mPartnerTPCNSigmaElectronMin, mPartnerTPCNSigmaElectronMax)
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
bool StNpeCuts::isBEMCElectron(StPicoTrack const *trk) const {
    // -- check for good BEMC electrons
    StPicoEmcPidTraits * Emc =  mPicoDst2->emcPidTraits(trk->emcPidTraitsIndex());
    float eoverp = Emc->e0()/trk->gPt()/TMath::CosH(getEta(trk));
    float phiDist = Emc->phiDist();
    float zDist = Emc->zDist();
    
    return eoverp > mElectronBemcEoverPMin && eoverp < mElectronBemcEoverPMax &&
    phiDist < mElectronBemcPhiDistMax && zDist < mElectronBemcZDistMax &&
    TMath::Sqrt(phiDist*phiDist + zDist*zDist) < mElectronBemcAssDistMax
    ;
}
     
// _________________________________________________________
bool StNpeCuts::isBSMDElectron(StPicoTrack const *trk) const {
    // -- check for good BSMD electrons
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

// __________________________________________________________
bool StNpeCuts::isGoodNpeEvent(StPicoDst const * const picoDst, int *aEventCuts) {
    mPicoDst2 = picoDst;
    return isGoodEvent(picoDst,aEventCuts);
}

