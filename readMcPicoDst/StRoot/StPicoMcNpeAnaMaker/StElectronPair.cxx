#include <limits>
#include <cmath>

#ifdef __ROOT__

#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "TVector3.h"

#include "StElectronPair.h"

ClassImp(StElectronPair)


StElectronPair::StElectronPair():
mElectronIdx(std::numeric_limits<unsigned short>::quiet_NaN()), 
mPartnerIdx(std::numeric_limits<unsigned short>::quiet_NaN()),
mMass(std::numeric_limits<unsigned short>::quiet_NaN()),
mPairDca(std::numeric_limits<float>::quiet_NaN()),
mPositionX(std::numeric_limits<float>::quiet_NaN()),
mPositionY(std::numeric_limits<float>::quiet_NaN()),
mPositionZ(std::numeric_limits<float>::quiet_NaN()),
mPhiV(std::numeric_limits<float>::quiet_NaN()),
mOpenAngle(std::numeric_limits<float>::quiet_NaN()),
mPhi(std::numeric_limits<float>::quiet_NaN()),
mEta(std::numeric_limits<float>::quiet_NaN()),
mAngle(std::numeric_limits<float>::quiet_NaN()),
mLength(std::numeric_limits<float>::quiet_NaN())
{
}
//------------------------------------
StElectronPair::StElectronPair(StElectronPair const * t) :
mElectronIdx(t->mElectronIdx), mPartnerIdx(t->mPartnerIdx),
mMass(t->mMass),
mPairDca(t->mPairDca),
mPositionX(t->mPositionX),
mPositionY(t->mPositionY),
mPositionZ(t->mPositionZ),
mPhiV(t->mPhiV),
mOpenAngle(t->mOpenAngle),
mPhi(t->mPhi),
mEta(t->mEta),
mAngle(t->mAngle),
mLength(t->mLength)
{
}
//------------------------------------
StElectronPair::StElectronPair(StPicoTrack const * const electron, StPicoTrack const * const partner,
                               unsigned short const electronIdx, unsigned short const partnerIdx, float const bField, StThreeVectorF pVtx) :
mElectronIdx(electronIdx), mPartnerIdx(partnerIdx),
mMass(std::numeric_limits<unsigned short>::quiet_NaN()),
mPairDca(std::numeric_limits<float>::quiet_NaN()),
mPositionX(std::numeric_limits<float>::quiet_NaN()),
mPositionY(std::numeric_limits<float>::quiet_NaN()),
mPositionZ(std::numeric_limits<float>::quiet_NaN()),
mPhiV(std::numeric_limits<float>::quiet_NaN()),
mOpenAngle(std::numeric_limits<float>::quiet_NaN()),
mPhi(std::numeric_limits<float>::quiet_NaN()),
mEta(std::numeric_limits<float>::quiet_NaN()),
mAngle(std::numeric_limits<float>::quiet_NaN()),
mLength(std::numeric_limits<float>::quiet_NaN())
{
    if ((!electron || !partner) || (electron->id() == partner->id()))
    {
        mElectronIdx = std::numeric_limits<unsigned short>::quiet_NaN();
        mPartnerIdx = std::numeric_limits<unsigned short>::quiet_NaN();
        return;
    }
    
    StPhysicalHelixD electronHelix = electron->dcaGeometry().helix();
    StPhysicalHelixD partnerHelix = partner->dcaGeometry().helix();
    
    // normal method
    pair<double,double> ss = electronHelix.pathLengths(partnerHelix);
    StThreeVectorD kAtDcaToPartner = electronHelix.at(ss.first);
    StThreeVectorD pAtDcaToElectron = partnerHelix.at(ss.second);
    
    // calculate DCA of partner to electron at their DCA
    StThreeVectorD VectorDca = kAtDcaToPartner - pAtDcaToElectron;
    mPairDca = static_cast<float>(VectorDca.mag());
    
    // calculate Lorentz vector of electron-partner pair
    StThreeVectorF const electronMomAtDca = electronHelix.momentumAt(ss.first, bField * kilogauss);
    StThreeVectorF const partnerMomAtDca = partnerHelix.momentumAt(ss.second, bField * kilogauss);
    
    StLorentzVectorF const electronFourMom(electronMomAtDca, electronMomAtDca.massHypothesis(M_ELECTRON));
    StLorentzVectorF const partnerFourMom(partnerMomAtDca, partnerMomAtDca.massHypothesis(M_ELECTRON));
    StLorentzVectorF const epairFourMom = electronFourMom + partnerFourMom;
    
    StThreeVectorF const epairMomAtDca = epairFourMom.vect();
    StThreeVectorF const Position = (kAtDcaToPartner + pAtDcaToElectron)/2.0;

    mAngle = epairMomAtDca.angle(Position-pVtx);
    mLength = TMath::Sin(mAngle)*Position.mag();
    
    mMass = epairFourMom.m();
    

    mPositionX = Position.x();

    mPositionY = Position.y();

    mPositionZ = Position.z();
    
    phiCalculation(partnerFourMom, electronFourMom, bField > 0 ? 1 : -1, mPhiV, mOpenAngle);
    
    mPhi = epairFourMom.phi();
    mEta = epairFourMom.pseudoRapidity();
    
}
void StElectronPair::phiCalculation(StLorentzVectorF const positron,StLorentzVectorF const electron, int mN, float &phiV, float &openangle)
{
    TVector3 ppp(positron.px(),positron.py(),positron.pz());
    TVector3 eee(electron.px(),electron.py(),electron.pz());
    TVector3 u=ppp+eee;
    TVector3 v=eee.Cross(ppp);
    TVector3 w=u.Cross(v);
    TVector3 nz(0.,0.,mN);
    TVector3 wc=u.Cross(nz);
    
    phiV =w.Angle(wc);
    openangle=ppp.Angle(eee);
    
}

#endif // __ROOT__
