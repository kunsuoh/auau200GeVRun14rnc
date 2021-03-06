#include "StPicoMcTrack.h"
#include "StPicoCut.h"
#include "StPicoConstants.h"
#include "StPicoDstMaker.h"
#include "TVector2.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
// -- MuMc
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMcEvent/StMcHit.hh"
#include "StMcEvent/StMcIstHitCollection.hh"
#include "StMcEvent/StMcIstLayerHitCollection.hh"
#include "StMcEvent/StMcIstHit.hh"
#include "StMcEvent/StMcPxlHitCollection.hh"
#include "StMcEvent/StMcPxlHit.hh"
#include "StMcEvent/StMcEvent.hh"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StRnDHit.h"
#include "StMcEvent/StMcTrack.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StPxlHit.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TDataSet.h"
#include "TObjectSet.h"
#include "StPxlDbMaker/StPxlDb.h"

ClassImp(StPicoMcTrack)
//----------------------------------------------------------------------------------
StPicoMcTrack::StPicoMcTrack() : mMcId(Pico::SHORTMAX),mGePid(Pico::SHORTMAX),mParentId(Pico::USHORTMAX), mAssoId(Pico::USHORTMAX),
nHitsPxl1(0), nHitsPxl2(0), nHitsIst(0),nHitsSst(0),
nHitsTpc(0),nHitsTpcCommon(0), mTruthBit(0),
mEta(std::numeric_limits<float>::quiet_NaN()), mMcMomentum(0.0,0.0,0.0), mOrigin(0.0,0.0,0.0)
{
}
//----------------------------------------------------------------------------------
StPicoMcTrack::StPicoMcTrack(StPicoMcTrack *t) :
mMcId(t->mMcId), mGePid(t->mGePid), mParentId(t->mParentId), mAssoId(t->mAssoId),
nHitsPxl1(t->nHitsPxl1), nHitsPxl2(t->nHitsPxl2), nHitsIst(t->nHitsIst),nHitsSst(t->nHitsSst),
nHitsTpc(t->nHitsTpc), nHitsTpcCommon(t->nHitsTpcCommon), mTruthBit(t->mTruthBit),
mEta(t->mEta), mMcMomentum(t->mMcMomentum), mOrigin(t->mOrigin)
{
}
//----------------------------------------------------------------------------------
StPicoMcTrack::StPicoMcTrack(StMcTrack *mcTrk, const StGlobalTrack *Trk = 0, const int nCommon = 0, const StPxlDb *mPxlDb = 0) :
mMcId(Pico::SHORTMAX), mGePid(Pico::SHORTMAX), mParentId(Pico::USHORTMAX), mAssoId(Pico::USHORTMAX),
nHitsPxl1(0), nHitsPxl2(0), nHitsIst(0),nHitsSst(0),
nHitsTpc(0),nHitsTpcCommon(nCommon), mTruthBit(0),
mEta(std::numeric_limits<float>::quiet_NaN()), mMcMomentum(0.0,0.0,0.0), mOrigin(0.0,0.0,0.0)
{
    mMcId = mcTrk->key();
    mGePid = mcTrk->geantId();
    
    if(mcTrk->parent()){
        cout << mcTrk->parent()->geantId() << " " << mcTrk->geantId() << endl;
        mParentId = mcTrk->parent()->key();
    }
    
    
    if(Trk)
        mAssoId = Trk->key();
    mMcMomentum = mcTrk->momentum();
    if(mcTrk->startVertex()){
        mOrigin = mcTrk->startVertex()->position();
    }
    mEta = mcTrk->pseudoRapidity();
    //Looping over hits to check for IdTruth
    if(Trk){
        StPtrVecHit PartnerSsdHits      = Trk->detectorInfo()->hits(kSsdId);
        StPtrVecHit PartnerIstHits      = Trk->detectorInfo()->hits(kIstId);
        StPtrVecHit PartnerPxlHits      = Trk->detectorInfo()->hits(kPxlId);
        
        int nPartnerSsdHits             = (int) PartnerSsdHits.size();
        int nPartnerIstHits             = (int) PartnerIstHits.size();
        int nPartnerPxlHits             = (int) PartnerPxlHits.size();
        
        uint ssdTruth=1;
        
        if(nPartnerSsdHits > 0) {
            for(int issdhit=0; issdhit<nPartnerSsdHits; issdhit++) {
                if(PartnerSsdHits[issdhit]->idTruth() == mMcId)
                    continue;
                ssdTruth=0;
                break;
            }
        }
        uint istTruth = 1;
        if(nPartnerIstHits > 0) {
            for(int iisthit=0; iisthit<nPartnerIstHits; iisthit++) {
                if( PartnerIstHits[iisthit]->idTruth() == mMcId)
                    continue;
                istTruth=0;
                break;
            }
        }
        uint pxlTruth1 = 1;
        uint pxlTruth2 = 1;
        int pxl1Hits = 0;
        int pxl2Hits = 0;
        
        if(nPartnerPxlHits > 0) {
            for(int ipxlhit=0; ipxlhit<nPartnerPxlHits; ipxlhit++) {
                StThreeVectorF pos = PartnerPxlHits[ipxlhit]->position();
                if (mcTrk->parent() && (mcTrk->parent()->geantId() == 1 || mcTrk->parent()->geantId() == 10007) && (mcTrk->geantId()==2 || mcTrk->geantId()==3))	cout << "rc position " <<  pos.x() << " " << pos.y() << endl;
                if(PartnerPxlHits[ipxlhit]->idTruth() == mMcId)
                    continue;
                float const R = pow(pos.x(),2.0) + pow(pos.y(),2.0);
                if(R > 3.5*3.5){
                    pxlTruth2 = 0;
                }
                else{
                    pxlTruth1 = 0;
                }
            }
        }
        const StPtrVecMcPxlHit mcPxlHits = mcTrk->pxlHits();
        //Loop over PXL hits to separate into layers
        for(int ipxlhit = 0; ipxlhit < (int)mcPxlHits.size(); ipxlhit++){
            StThreeVectorF pos = mcPxlHits.at(ipxlhit)->position();
            if (mcTrk->parent() && (mcTrk->parent()->geantId() == 1 || mcTrk->parent()->geantId() == 10007) && (mcTrk->geantId()==2 || mcTrk->geantId()==3))      {
                
                int sector = mcPxlHits.at(ipxlhit)->sector();
                int ladder = mcPxlHits.at(ipxlhit)->ladder();
                int sensor = mcPxlHits.at(ipxlhit)->sensor();
                Double_t local[3] = {0,0,0};
                Double_t master[3] = {0,0,0};
                
                TGeoHMatrix *combP = (TGeoHMatrix *)mPxlDb->geoHMatrixSensorOnGlobal(sector, ladder,sensor);
                combP->LocalToMaster(local, master);
                
                
                
                cout << "mc position " << pos.x() << " " << pos.y() << endl;
            }
            if((int)mcPxlHits.at(ipxlhit)->ladder() > 1)
                pxl2Hits++;
            else
                pxl1Hits++;
        }
        nHitsPxl1 = pxl1Hits;
        nHitsPxl2 = pxl2Hits;
        
        nHitsIst = (UShort_t)mcTrk->istHits().size();
        nHitsSst = (UShort_t)mcTrk->ssdHits().size();
        nHitsTpc = (UShort_t)mcTrk->tpcHits().size();
        
        mTruthBit = 0;
        mTruthBit= mTruthBit | pxlTruth1;
        mTruthBit= mTruthBit | (pxlTruth2<<1);
        mTruthBit= mTruthBit | (istTruth<<2);
        mTruthBit= mTruthBit | (ssdTruth<<3);
    }
}

