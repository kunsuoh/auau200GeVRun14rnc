#ifndef StNpeCuts_H
#define StNpeCuts_H

/* **************************************************
 *  Cut class for HF analysis
 *  - Based on PicoCuts class
 *
 *  Initial Authors:
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *          **Jochen Thaeder  (jmthader@lbl.gov)
 *
 *  Contributing Authors
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Guannan Xie     (guannanxie@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

#include "StPicoCutsBase/StPicoCutsBase.h"

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
    void setCutSecondaryPair(float dcaDaughtersMax, float massMin, float massMax);

    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- GETTER for single CUTS
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    const float&    cutSecondaryPairDcaDaughtersMax()       const;
    const float&    cutSecondaryPairDecayLengthMin()        const;
    const float&    cutSecondaryPairDecayLengthMax()        const;
    const float&    cutSecondaryPairCosThetaMin()           const;
    const float&    cutSecondaryPairMassMin()               const;
    const float&    cutSecondaryPairMassMax()               const;
    
    
private:
    
    StNpeCuts(StNpeCuts const &);
    StNpeCuts& operator=(StNpeCuts const &);
    
    // ------------------------------------------
    // -- Pair cuts for secondary pair
    // ------------------------------------------
    float mSecondaryPairDcaDaughtersMax;
    float mSecondaryPairDecayLengthMin;
    float mSecondaryPairDecayLengthMax;
    float mSecondaryPairCosThetaMin;
    float mSecondaryPairMassMin;
    float mSecondaryPairMassMax;
    
    ClassDef(StNpeCuts,1)
};

inline void StNpeCuts::setCutSecondaryPair(float dcaDaughtersMax, float massMin, float massMax)  {
    mSecondaryPairDcaDaughtersMax = dcaDaughtersMax;
    mSecondaryPairMassMin = massMin; mSecondaryPairMassMax = massMax;
}


inline const float&    StNpeCuts::cutSecondaryPairDcaDaughtersMax()       const { return mSecondaryPairDcaDaughtersMax; }
inline const float&    StNpeCuts::cutSecondaryPairDecayLengthMin()        const { return mSecondaryPairDecayLengthMin; }
inline const float&    StNpeCuts::cutSecondaryPairDecayLengthMax()        const { return mSecondaryPairDecayLengthMax; }
inline const float&    StNpeCuts::cutSecondaryPairCosThetaMin()           const { return mSecondaryPairCosThetaMin; }
inline const float&    StNpeCuts::cutSecondaryPairMassMin()               const { return mSecondaryPairMassMin; }
inline const float&    StNpeCuts::cutSecondaryPairMassMax()               const { return mSecondaryPairMassMax; }

#endif
