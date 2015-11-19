#ifndef StNpeCuts_H
#define StNpeCuts_H

/* **************************************************
 *  Cut class for HF analysis
 *  - Based on PicoCuts class
 *
 *  Initial Authors:
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)
 *
 *  Contributing Authors
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Guannan Xie     (guannanxie@lbl.gov)
 *          **Kunsu OH        (kunsuoh@gmail.com)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

#include "StPicoNpeEventMaker/StElectronPair.h"
#include "StPicoCutsBase/StPicoCutsBase.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StThreeVectorF.hh"

class StPicoTrack;
class StPicoEvent;
class StPicoDst;

class StNpeCuts : public StPicoCutsBase
{
public:
    
    StNpeCuts();
    StNpeCuts(const Char_t *name);
    ~StNpeCuts();
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    virtual void init() { initBase(); }
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- SETTER for CUTS
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    void setCutElectronPair(float dcaDaughtersMax, float massMin, float massMax);
    void setCutElectronNHitsFitMax(int i);
    void setCutElectronNHitsdEdxMax(int i);
    void setCutPt(float fmin, float fmax);
    void setCutEta(float fmin, float fmax);
    void setCutDca(float f);
    void setCutElectronRequireHFT(bool b);
    
    void setCutPartnerElectronNHitsFitMax(int i);
    void setCutPartnerElectronNHitsdEdxMax(int i);
    void setCutPartnerPt(float fmin, float fmax);
    void setCutPartnerEta(float fmin, float fmax);
    void setCutPartnerElectronRequireHFT(bool b);

    void setCutTPCNSigmaElectron(float fmin, float fmax);
    void setCutPartnerTPCNSigmaElectron(float fmin, float fmax);
    void setCutBemcPid(bool pid, float epmin, float epmax, float phi, float z, float ass);
    void setCutBsmdPid(bool pid, int eta, int phi);
    void setCutTofPid(bool pid, float beta);
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- GETTER for single CUTS
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    StPicoEmcPidTraits * hasEmcPid(StPicoTrack const * const trk) const;
    
    void setPicoDst(StPicoDst const * picoDst);
    bool isTPCElectron(StPicoTrack const *trk, float min, float max) const;
    bool isTOFElectron(StPicoTrack const *trk) const;
    bool isBEMCElectron(StPicoTrack const *trk) const;
    bool isBSMDElectron(StPicoTrack const *trk) const;
    bool isGoodElectronPair(StElectronPair const* epair) const;

    bool isGoodInclusiveElectron(StPicoTrack const *trk) const;
    bool isGoodTaggedElectron(StPicoTrack const *trk) const;
    bool isGoodPartnerElectron(StPicoTrack const *trk) const;
    float getEta(StPicoTrack const *trk)         const;
    float getDca(StPicoTrack const *trk)         const;
    StThreeVectorF getpVtx()                     const;
    
    const float&    cutElectronPairDcaDaughtersMax()       const;
    const float&    cutElectronPairDecayLengthMin()        const;
    const float&    cutElectronPairDecayLengthMax()        const;
    const float&    cutElectronPairCosThetaMin()           const;
    const float&    cutElectronPairMassMin()               const;
    const float&    cutElectronPairMassMax()               const;

    
    
private:


    StNpeCuts(StNpeCuts const &);
    StNpeCuts& operator=(StNpeCuts const &);

    const StPicoDst*  mPicoDst2;   //! ptr to picoDst
    
    // ------------------------------------------
    // -- Pair cuts for electron pair
    // ------------------------------------------
    float mElectronPairDcaDaughtersMax;
    float mElectronPairDecayLengthMin;
    float mElectronPairDecayLengthMax;
    float mElectronPairCosThetaMin;
    float mElectronPairMassMin;
    float mElectronPairMassMax;
    
    // ------------------------------------------
    // -- Track cuts for tagged electron
    // ------------------------------------------
    int mElectronNHitsdEdxMax;
    int mElectronNHitsFitMax;
    float mElectronPtMin;
    float mElectronPtMax;
    float mElectronEtaMin;
    float mElectronEtaMax;
    float mElectronDca;
    int mElectronBsmdNEta;
    int mElectronBsmdNPhi;
    bool mElectronRequireHFT;
    float mElectronTPCNSigmaElectronMin;
    float mElectronTPCNSigmaElectronMax;
    float mElectronBemcEoverPMin;
    float mElectronBemcEoverPMax;
    float mElectronBemcPhiDistMax;
    float mElectronBemcZDistMax;
    float mElectronBemcAssDistMax;
    bool mElectronBemcPid;
    bool mElectronBsmdPid;
    bool mElectronTofPid;
    float mElectronTofBeta;
    
    // ------------------------------------------
    // -- Track cuts for partner electron
    // ------------------------------------------
    int mPartnerElectronNHitsdEdxMax;
    int mPartnerElectronNHitsFitMax;
    float mPartnerElectronPtMin;
    float mPartnerElectronPtMax;
    float mPartnerElectronEtaMin;
    float mPartnerElectronEtaMax;
    float mPartnerTPCNSigmaElectronMin;
    float mPartnerTPCNSigmaElectronMax;
    bool mPartnerElectronRequireHFT;

    
    ClassDef(StNpeCuts,1)
};

inline void StNpeCuts::setCutElectronPair(float dcaDaughtersMax, float massMin, float massMax)  {
    mElectronPairDcaDaughtersMax = dcaDaughtersMax;
    mElectronPairMassMin = massMin;
    mElectronPairMassMax = massMax;
}
inline void StNpeCuts::setCutElectronNHitsdEdxMax(int i)  {
    mElectronNHitsdEdxMax = i;
}
inline void StNpeCuts::setCutElectronNHitsFitMax(int i)  {
    mElectronNHitsFitMax = i;
}
inline void StNpeCuts::setCutPt(float fmin, float fmax)  {
    mElectronPtMin = fmin;
    mElectronPtMax = fmax;
}
inline void StNpeCuts::setCutEta(float fmin, float fmax)  {
    mElectronEtaMin = fmin;
    mElectronEtaMax = fmax;
}
inline void StNpeCuts::setCutDca(float f)  {
    mElectronDca = f;
}
inline void StNpeCuts::setCutElectronRequireHFT(bool b)  {
    mElectronRequireHFT = b;
}
inline void StNpeCuts::setCutPartnerElectronNHitsdEdxMax(int i)  {
    mPartnerElectronNHitsdEdxMax = i;
}
inline void StNpeCuts::setCutPartnerElectronNHitsFitMax(int i)  {
    mPartnerElectronNHitsFitMax = i;
}
inline void StNpeCuts::setCutPartnerPt(float fmin, float fmax)  {
    mPartnerElectronPtMin = fmin;
    mPartnerElectronPtMax = fmax;

}
inline void StNpeCuts::setCutPartnerEta(float fmin, float fmax)  {
    mPartnerElectronEtaMin = fmin;
    mPartnerElectronEtaMax = fmax;

}
inline void StNpeCuts::setCutPartnerElectronRequireHFT(bool b)  {
    mPartnerElectronRequireHFT = b;
}

inline void StNpeCuts::setCutTPCNSigmaElectron(float fmin, float fmax)  {
    mElectronTPCNSigmaElectronMin = fmin;
    mElectronTPCNSigmaElectronMax = fmax;
}
inline void StNpeCuts::setCutPartnerTPCNSigmaElectron(float fmin, float fmax)  {
    mPartnerTPCNSigmaElectronMin = fmin;
    mPartnerTPCNSigmaElectronMax = fmax;

}
inline void StNpeCuts::setCutBemcPid(bool pid, float epmin, float epmax, float phi, float z, float ass)  {
    mElectronBemcPid = pid;
    mElectronBemcEoverPMin = epmin;
    mElectronBemcEoverPMax = epmax;
    mElectronBemcPhiDistMax = phi;
    mElectronBemcZDistMax = z;
    mElectronBemcAssDistMax = ass;
    
}
inline void StNpeCuts::setCutBsmdPid(bool pid, int eta, int phi)  {
    mElectronBsmdPid = pid;
    mElectronBsmdNEta = eta;
    mElectronBsmdNPhi = phi;
}
inline void StNpeCuts::setCutTofPid(bool pid, float beta)  {
    mElectronTofPid = pid;
    mElectronTofBeta = beta;
}
inline void StNpeCuts::setPicoDst(StPicoDst const * picoDst)  {
    mPicoDst2=picoDst;
}


inline const float&    StNpeCuts::cutElectronPairDcaDaughtersMax()       const { return mElectronPairDcaDaughtersMax; }
inline const float&    StNpeCuts::cutElectronPairDecayLengthMin()        const { return mElectronPairDecayLengthMin; }
inline const float&    StNpeCuts::cutElectronPairDecayLengthMax()        const { return mElectronPairDecayLengthMax; }
inline const float&    StNpeCuts::cutElectronPairCosThetaMin()           const { return mElectronPairCosThetaMin; }
inline const float&    StNpeCuts::cutElectronPairMassMin()               const { return mElectronPairMassMin; }
inline const float&    StNpeCuts::cutElectronPairMassMax()               const { return mElectronPairMassMax; }
#endif
