#ifndef StPicoMcAnaMaker_h
#define StPicoMcAnaMaker_h

#include "StMaker.h"
#include "TNtuple.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StEvent/StDcaGeometry.h"

#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoMcTrack.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;
class TString;
class TH1F;
class TH2F;
class TTree;
class StDcaGeometry;
class StPicoMcAnaHists;

class StPicoMcAnaMaker : public StMaker
{
 public: 
  StPicoMcAnaMaker( TString  name,   TString outname, StPicoDstMaker *picoMaker);
  virtual ~StPicoMcAnaMaker(){;};  
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  bool isGoodEvent(StPicoEvent *event);
  bool isGoodTrack(StPicoMcTrack const * const trk);
  bool isHftTrack(StPicoMcTrack const * const trk);
  bool isGoodTrack(StPicoTrack const * const trk, StPicoEvent const * const evt);
  bool isRcTrack(StPicoMcTrack const * const trk ,StPicoDst const * const PicoDst, int &id);

 private:
  StPicoDstMaker   *mPicoDstMaker;
  StPicoDst        *mPicoDst;

  TNtuple* nt;
  TFile* mFile;
  TString outfile;

  ClassDef(StPicoMcAnaMaker, 1)
};
inline bool StPicoMcAnaMaker::isGoodEvent(StPicoEvent *event)
{
  return( fabs(event->primaryVertex().z())<5.0 );   
}
inline bool StPicoMcAnaMaker::isGoodTrack(StPicoMcTrack const * const trk)
{
  return ( trk->hitsTpc()>15 && fabs(trk->Mom().pseudoRapidity())<1  &&
	   trk->Mom().perp() > 0.2  );
inline bool StPicoMcAnaMaker::isGoodTrack(StPicoTrack const * const trk, StPicoEvent const * const evt)
{
  StThreeVectorF pVtx = evt->primaryVertex();
  float B = evt->bField();
  StThreeVectorF mom = trk->gMom(pVtx,B);
  return( trk->nHitsMax()>15 && std::fabs(mom.pseudoRapidity()) &&
	  mom.perp()>0.2);
}
inline bool StPicoMcAnaMaker::isHftTrack(StPicoMcTrack const * const trk)
{
  if( trk->hitsPxl1() == 0 ||  trk->hitsPxl2() == 0 || trk->hitsIst() == 0 ) 
    return false;
  return true;
}
#endif

