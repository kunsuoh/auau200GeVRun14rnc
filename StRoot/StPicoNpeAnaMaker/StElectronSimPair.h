#ifndef StElectronSimPair_hh
#define StElectronSimPair_hh
#ifdef __ROOT__

/* **************************************************
 *  A specialized pair class for calculating electron pair
 *  lorentz vector and topological decay parameters
 *  and storing them.
 *
 *  Authors:  **Kunsu OH        (kunsuoh@gmail.com)
 *
 *  **Code Maintainer
 * **************************************************
 */

#include "TObject.h"
#include "StLorentzVectorF.hh"

class StPicoTrack;
class StPicoEvent;

class StElectronSimPair : public TObject
{
public:
    StElectronSimPair();
    StElectronSimPair(StElectronSimPair const *);
    StElectronSimPair(StPicoTrack const * Electron, StPicoTrack const * Partner,
                   unsigned short electronIdx,unsigned short partnerIdx, float bField);
    ~StElectronSimPair() {}// please keep this non-virtual and NEVER inherit from this class
    
    unsigned short   electronIdx() const;	// tagged electron idx
    unsigned short   partnerIdx() const;	// partner electron idx
    float pairMass()    const;                     // electron pair mass
    float pairDca() const;                  // DCA between tagged and partner
    float positionX() const;                // conversion position X
    float positionY() const;                // conversion position Y
    float positionZ() const;                // conversion position Z
    float phiV() const;
    float openAngle() const;
    float phi() const;
    float eta() const;
    float angle() const;
    float length() const;
    
    
private:
    // disable copy constructor and assignment operator by making them private (once C++11 is available in STAR you can use delete specifier instead)
    StElectronSimPair(StElectronSimPair const &);
    StElectronSimPair& operator=(StElectronSimPair const &);
    void phiCalculation(StLorentzVectorF const, StLorentzVectorF const, int, float &, float &);

    unsigned short mElectronIdx;    // index of electron track in StPicoDstEvent (2 Bytes)
    unsigned short mPartnerIdx;     // index of partner track in StPicoDstEvent (2 Bytes)
    float mMass;                    // mass (4 Bytes)
    float mPairDca;                 // pair dca (4 Bytes)
    float mPositionX;               // conversion position x (4 Bytes)
    float mPositionY;               // conversion position y (4 Bytes)
    float mPositionZ;               // conversion position z (4 Bytes)
    float mPhiV;
    float mOpenAngle;
    float mPhi;
    float mEta;
    float mAngle;
    float mLength;
    
    ClassDef(StElectronSimPair,2)
};
inline unsigned short   StElectronSimPair::electronIdx() const     { return mElectronIdx;                        }
inline unsigned short   StElectronSimPair::partnerIdx() const      { return mPartnerIdx;                         }
inline float StElectronSimPair::pairMass()    const                { return mMass;     }
inline float StElectronSimPair::pairDca() const                    { return mPairDca;                            }
inline float StElectronSimPair::positionX() const                  { return mPositionX;}
inline float StElectronSimPair::positionY() const                  { return mPositionY;}
inline float StElectronSimPair::positionZ() const                  { return mPositionZ;}
inline float StElectronSimPair::phiV() const                       { return mPhiV;}
inline float StElectronSimPair::openAngle() const                  { return mOpenAngle;}
inline float StElectronSimPair::phi() const                        { return mPhi;}
inline float StElectronSimPair::eta() const                        { return mEta ;}
inline float StElectronSimPair::angle() const                      { return mAngle;}
inline float StElectronSimPair::length() const                     { return mLength ;}

#endif
#endif

