#include "StPicoMcAnaMaker.h"
#include "StPicoMcAnaHists.h"
#include <assert.h>
#include <map>
#include <utility>
#include "Riostream.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TClassTable.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TStyle.h"
#include "SystemOfUnits.h"
#include "StarRoot/TPolynomial.h"
#include "StDcaGeometry.h"
#include "TRSymMatrix.h"
#include "THelixTrack.h"
#include "StBichsel/Bichsel.h"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoMcTrack.h"


ClassImp(StPicoMcAnaMaker)

//-----------------------------------------------------------------------------
StPicoMcAnaMaker::StPicoMcAnaMaker(const char* name, const char* outname, StPicoDstMaker *picoMaker)
: StMaker(name), mHists(NULL)
{
  mHists = new StPicoMcAnaHists(outname);
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
}
//-----------------------------------------------------------------------------
Int_t StPicoMcAnaMaker::Init(){
  return kStOk;
}
//-----------------------------------------------------------------------------
Int_t StPicoMcAnaMaker::Make(){
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  StThreeVectorF pVtxErr(0.,0.,0.);
  if(event) {
    pVtx = event->primaryVertex();
    pVtxErr = event->primaryVertexError();
  }

  if(!isGoodEvent(event)) return kStOk; 
  //mHists->addEvent(event);

  int nMcTracks =  mPicoDst->numberOfMcTracks();

  for(int i_Mc=0; i_Mc<nMcTracks; i_Mc++){
    StPicoMcTrack *mcTrk = (StPicoMcTrack*)mPicoDst->mctrack(i_Mc);
    if(!isGoodTrack(mcTrk) || !isPrimaryTrack(mcTrk)) continue;
    mHists->addMcTrack(mcTrk);

    if( mcTrk->assoId() == Pico::USHORTMAX )
      continue;
    if( !isHftTrack(mcTrk) )
      continue;
    int temp = Pico::USHORTMAX ;
    for(int i_Rc =0; i_Rc<mPicoDst->numberOfTracks(); ++i_Rc){
      StPicoTrack *Trk = (StPicoTrack*)mPicoDst->track(i_Rc);
      if(mcTrk->assoId() == Trk->id() ) {
	temp = i_Rc;
	break;
      }
    }
    if (temp == Pico::USHORTMAX) continue;
    StPicoTrack *Trk = (StPicoTrack*)mPicoDst->track(temp);
    if(!isGoodTrack(Trk,event)) continue;

    mHists->addMatchedTrack(mcTrk, Trk);
  }
  return kStOk;
}
Int_t StPicoMcAnaMaker::Finish()
{
  mHists->closeFile();
  return kStOk;
}
