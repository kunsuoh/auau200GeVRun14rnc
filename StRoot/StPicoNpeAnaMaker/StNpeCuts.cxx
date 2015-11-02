#include <limits>

#include "StNpeCuts.h"

ClassImp(StNpeCuts)

// _________________________________________________________
StNpeCuts::StNpeCuts() : StPicoCutsBase("NpeCutsBase"),
mElectronPairDcaDaughtersMax(std::numeric_limits<float>::max()),
mElectronPairDecayLengthMin(std::numeric_limits<float>::min()), mElectronPairDecayLengthMax(std::numeric_limits<float>::max()),
mElectronPairCosThetaMin(std::numeric_limits<float>::min()),
mElectronPairMassMin(std::numeric_limits<float>::min()), mElectronPairMassMax(std::numeric_limits<float>::max()),
mElectronNHitdEdxMax(std::numeric_limits<int>::max()),
mElectronBsmdNEta(std::numeric_limits<int>::min()),
mElectronBsmdNPhi(std::numeric_limits<int>::min()),
mElectronRequireHFT(false),
mElectronPtMin(std::numeric_limits<float>::min()),
mElectronPtMax(std::numeric_limits<float>::max()),
mElectronEtaMin(std::numeric_limits<float>::min()),
mElectronEtaMax(std::numeric_limits<float>::max()),
mElectronDca(std::numeric_limits<float>::max()),
mElectronTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mElectronTPCNSigmaElectronMax(std::numeric_limits<float>::max()),
mElectronBemcEoverPMin(std::numeric_limits<float>::min()),
mElectronBemcEoverPMax(std::numeric_limits<float>::max()),
mElectronBemcPhiDistMax(std::numeric_limits<float>::max()),
mElectronBemcZDistMax(std::numeric_limits<float>::max()),
mElectronBemcAssDistMax(std::numeric_limits<float>::max()),
mPartnerElectronNHitsdEdxMax(std::numeric_limits<int>::min()),
mPartnerElectronPtMin(std::numeric_limits<float>::min()),
mPartnerElectronPtMax(std::numeric_limits<float>::max()),
mPartnerElectronEtaMin(std::numeric_limits<float>::min()),
mPartnerElectronEtaMax(std::numeric_limits<float>::max()),
mPartnerTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mPartnerTPCNSigmaElectronMax(std::numeric_limits<float>::max())
{

    // -- default constructor
}

// _________________________________________________________
StNpeCuts::StNpeCuts(const Char_t *name) : StPicoCutsBase(name),
mElectronPairDcaDaughtersMax(std::numeric_limits<float>::max()),
mElectronPairDecayLengthMin(std::numeric_limits<float>::min()), mElectronPairDecayLengthMax(std::numeric_limits<float>::max()),
mElectronPairCosThetaMin(std::numeric_limits<float>::min()),
mElectronPairMassMin(std::numeric_limits<float>::min()), mElectronPairMassMax(std::numeric_limits<float>::max()),
mElectronNHitdEdxMax(std::numeric_limits<int>::max()),
mElectronBsmdNEta(std::numeric_limits<int>::min()),
mElectronBsmdNPhi(std::numeric_limits<int>::min()),
mElectronRequireHFT(false),
mElectronPtMin(std::numeric_limits<float>::min()),
mElectronPtMax(std::numeric_limits<float>::max()),
mElectronEtaMin(std::numeric_limits<float>::min()),
mElectronEtaMax(std::numeric_limits<float>::max()),
mElectronDca(std::numeric_limits<float>::max()),
mElectronTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mElectronTPCNSigmaElectronMax(std::numeric_limits<float>::max()),
mElectronBemcEoverPMin(std::numeric_limits<float>::min()),
mElectronBemcEoverPMax(std::numeric_limits<float>::max()),
mElectronBemcPhiDistMax(std::numeric_limits<float>::max()),
mElectronBemcZDistMax(std::numeric_limits<float>::max()),
mElectronBemcAssDistMax(std::numeric_limits<float>::max()),
mPartnerElectronNHitsdEdxMax(std::numeric_limits<int>::min()),
mPartnerElectronPtMin(std::numeric_limits<float>::min()),
mPartnerElectronPtMax(std::numeric_limits<float>::max()),
mPartnerElectronEtaMin(std::numeric_limits<float>::min()),
mPartnerElectronEtaMax(std::numeric_limits<float>::max()),
mPartnerTPCNSigmaElectronMin(std::numeric_limits<float>::min()),
mPartnerTPCNSigmaElectronMax(std::numeric_limits<float>::max())
{

    // -- constructor
}

// _________________________________________________________
StNpeCuts::~StNpeCuts() {
    // destructor
    
}


// _________________________________________________________
StNpeCuts::isGoodElectronPair(StElectronPair const* const epair) {
    // -- check for good electron pairs
    StPicoTrack const* electron = mPicoDstMaker->picoDst()->track(epair->electronIdx());
    StPicoTrack const* partner = mPicoDstMaker->picoDst()->track(epair->partnerIdx());

    return
    isGoodTaggedElectron(electron) && isGoodPartnerElectron(partner) &&
    epair->pairMass() > mElectronPairMassMin && epair->pairMass() < mElectronPairMassMax &&
    epair->pairDca() < mElectronPairDcaDaughtersMax
    ;
}// _________________________________________________________
StNpeCuts::isGoodTaggedElectron(StPicoTrack const *trk) {
    // -- check for good tagged electron for electron pairs

    bool taggedElectronCut =
    trk->nHitsFit() >= mNHitsFitMax &&
    trk->nHitsDedx() >= mElectronNHitdEdxMax &&
    trk->gPt() >= mElectronPtMin && trk->gPt() < mElectronPtMax &&
    getEta(trk) > mElectronEtaMin && getEta(trk) < mElectronEtaMax &&
    getDca(trk) < mElectronDca &&
    (!mElectronRequireHFT || trk->isHFTTrack());
    
    return taggedElectronCut && isTPCElectron(trk, mElectronTPCNSigmaElectronMin, mElectronTPCNSigmaElectronMax) && isBEMCElectron(trk) && isBSMDElectron(trk) ;
}
// _________________________________________________________
StNpeCuts::isGoodPartnerElectron(StPicoTrack const *trk) {
    // -- check for good partner electron for electron pairs
    
    bool partnerElectronCut =
    trk->nHitsFit() >= mNHitsFitMax &&
    trk->gPt() >= mPartnerElectronPtMin && trk->gPt() < mPartnerElectronPtMax &&
    getEta(trk) > mPartnerElectronEtaMin && getEta(trk) < mPartnerElectronEtaMax ;

    return partnerElectronCut && isTPCElectron(trk, mPartnerTPCNSigmaElectronMin, mPartnerTPCNSigmaElectronMax)
    ;
}
// _________________________________________________________
StNpeCuts::isTPCElectron(StPicoTrack const *trk, min, max) {
    // -- check for good TPC electrons
    float nSigma = trk->nSigmaElectron();

    return
    nSigma > min && nSigma < max;
}
// _________________________________________________________
StNpeCuts::isBEMCElectron(StPicoTrack const *trk) {
    // -- check for good BEMC electrons
    StPicoEmcPidTraits * Emc =  mPicoDst->emcPidTraits(trk->emcPidTraitsIndex());
    float eoverp = Emc->e0()/trk->gPt()/TMath::CosH(getEta(trk));
    float phiDist = Emc->phiDist();
    float zDist = Emc->zDist();
    
    return eoverp > mElectronBemcEoverPMin && eoverp < mElectronBemcEoverPMax &&
    phiDist < mElectronBemcPhiDistMax && zDist < mElectronBemcZDistMax &&
    TMath::Sqrt(phiDist*phiDist + zDist*zDist) < mElectronBemcAssDistMax
    ;
}
StNpeCuts::isBSMDElectron(StPicoTrack const *trk) {
    // -- check for good BSMD electrons
    StPicoEmcPidTraits * Emc =  mPicoDst->emcPidTraits(trk->emcPidTraitsIndex());
    int nphi = Emc->nPhi();
    int neta = Emc->nEta();

    return neta > mElectronBsmdNEta && nphi > mElectronBsmdNPhi ;
}
// _________________________________________________________
StNpeCuts::getEta(StPicoTrack const *trk) {
    return trk->gMom(getpVtx(), mPicoDst->event()->bField()).pseudoRapidity();
}
// _________________________________________________________
StNpeCuts::getDca(StPicoTrack const *trk) {
    return trk->dcaGeometry().helix().curvatureSignedDistance(getpVtx().x(),getpVtx().y());
}
