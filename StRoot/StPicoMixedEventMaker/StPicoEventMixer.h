#ifndef StPicoEventMixer_hh
#define StPicoEventMixer_hh

/* **************************************************
 * Class stores event buffer used in event mixing. Mixing
 * is done automatically once buffer reaches defined maximum.
 * User should rpesonalize mixEvent() method to cosntruct 
 * desired background.
 *
 * **************************************************
 * 
 * Initial Authors:
 *          **Michael Lomnitz (mrlomnitz@lbl.gov)
 *          Musta Mustafa   (mmustafa@lbl.gov)
 *
 *  ** Code maintainer 
 *
 * **************************************************
 */

#include <vector>

#include "StThreeVectorF.hh"

class TTree;
class TH2F;
class StPicoEvent;
class StPicoTrack;
class StPicoDst;
class StMixerTrack;
class StMixerEvent;
class StMixerPair;
//class StThreeVector;

class StHFCuts;

class StPicoEventMixer {
 public: 
  StPicoEventMixer();
  ~StPicoEventMixer(){;};
  bool addPicoEvent(StPicoDst const* picoDst, StHFCuts const* mHFCuts);
  void setEventBuffer(int buffer);
  void mixEvents(StHFCuts *mHFCuts);
  bool isCloseTrack(StPicoTrack const& trk, StThreeVectorF const& pVtx);
  void finish();
 private:
  void fill(StMixerPair const* const);
  void fillFG(StMixerPair const* const);
  bool isMixerPion(StMixerTrack const&);
  bool isMixerKaon(StMixerTrack const&);

  TH2F* mVtx;
  TH2F* mFgVtx;
  TH2F* mForeground;
  TH2F* mBackground;
  //TTree * ntp_ME;
  std::vector <StMixerEvent*> mEvents; 

  unsigned short int mEventsBuffer; 
  unsigned short int filledBuffer;
  float dca1, dca2, dcaDaughters, theta_hs, decayL_hs;
  float pt_hs, mass_hs, eta_hs, phi_hs;
};

inline void StPicoEventMixer::setEventBuffer(int buffer){ mEventsBuffer = buffer;}
			    
    
#endif
