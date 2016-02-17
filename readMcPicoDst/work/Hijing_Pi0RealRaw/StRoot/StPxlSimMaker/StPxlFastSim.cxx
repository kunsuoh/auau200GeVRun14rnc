/*
 *
 * Author: A. Rose, LBL, Y. Fisyak, BNL, M. Miller, MIT, M. Mustafa
 *
 *
 **********************************************************
 * $Log: StPxlFastSim.cxx,v $
 * Revision 1.12  2015/05/14 18:57:52  smirnovd
 * Squashed commit of the following:
 *
 * StPxlFastSim: Streamlined creation of PXL hits by making use of StPxlUtil/StPxlDigiHit
 *
 * StPxlHitMaker: Updated comments
 *
 * StPxlHitMaker: Streamlined creation of PXL hits by making use of StPxlUtil/StPxlDigiHit
 *
 * StPxlDigiHit: A helper to manipulate local hit position in StPxlHit
 *
 * StPxlConsts: Define constants in namespace
 *
 * For safety reasons, the intentions is to move the constants into the namespace
 * and get rid of those defined in the global space.
 *
 * Revision 1.11  2015/05/07 21:24:31  smirnovd
 * StPxlSimMaker: Switched to using consts from StPxlUtil/
 *
 * Revision 1.10  2015/03/13 18:45:01  perev
 * Roll back
 *
 * Revision 1.8  2015/01/27 19:11:49  mstftsm
 * Set idTruth of StPxlHit to -999 if parentTrack of mcHit does not exist (for protection).
 *
 * Revision 1.7  2015/01/27 01:31:09  smirnovd
 * Minor refactoring of StPxlFastSim::distortHit() to include a new warning for unphysical hit position
 *
 * Revision 1.6  2014/07/03 19:46:37  mstftsm
 * Revereted the changes made for the pileup adder. That does not belong to the master branch.
 *
 * Revision 1.4  2014/03/13 17:00:19  mstftsm
 * StPxlSimMaker has a method to switch on random seed for StRandom generatos in simulators. Default is not a random seed.
 *
 * Revision 1.2  2013/11/14 19:10:27  mstftsm
 * StMcPxlHit has been changed to be on local coordinates. We no longer transfor from global to local before smearing
 *
 * Revision 1.1  2013/05/12 21:43:32  jeromel
 * Initial revision, code peer review closed 2013/05/06
 *
 * Revision 1.5  2013/05/09 02:58:36  mstftsm
 * Fixed a bug which called for sensor local Z value.
 *
 *
 */

#include <stdio.h>

#include "StMessMgr.h"
#include "Stypes.h"
#include "Stiostream.h"
#include "StPxlFastSim.h"
#include "StEvent/StPxlHit.h"
#include "StEvent/StPxlHitCollection.h"
#include "StPxlRawHitMaker/StPxlRawHit.h"
#include "StPxlRawHitMaker/StPxlRawHitMaker.h"
#include "StPxlRawHitMaker/StPxlRawHitCollection.h"
#include "StMcEvent/StMcPxlHitCollection.hh"
#include "StMcEvent/StMcPxlHit.hh"
#include "tables/St_HitError_Table.h"
#include "StarClassLibrary/StRandom.hh"
#include "StThreeVectorF.hh"
#include "StPxlDbMaker/StPxlDb.h"
#include "StPxlUtil/StPxlConstants.h"
#include "StPxlUtil/StPxlDigiHit.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"

#include "TDataSet.h"
#include "TObjectSet.h"
#include "TF1.h"

StPxlFastSim::~StPxlFastSim()
{
    if (mRandom) delete mRandom;
    if (mPxlDb) delete mPxlDb;
}
//____________________________________________________________
Int_t StPxlFastSim::initRun(const TDataSet& calib_db, const TObjectSet* pxlDbDataSet, const Int_t run)
{
    // run is not used in the current implementation, but might be necessary in the future.
    LOG_INFO << "StPxlFastSim::init()" << endm;
    
    if(pxlDbDataSet != 0)
    {
        mPxlDb = (StPxlDb *)pxlDbDataSet->GetObject();
        if (!mPxlDb)
        {
            LOG_ERROR << "StPxlFastSim - E - mPxlDb is not available" << endm;
            return kStErr;
        }
        else
        {
            LOG_INFO << "StPxlFastSim - Using geometry from pxlDB" <<endm;
        }
    }
    else
    {
        LOG_INFO << "StPxlFastSim - Using ideal geometry" <<endm;
    }
    
    if (!mRandom) mRandom = new StRandom();
    if(mUseRandomSeed)
    {
        Int_t seed = time(NULL);
        mRandom->setSeed(seed);
        LOG_INFO << "StPxlFastSim - smearing random generator is using random seed = " << seed <<endm;
    }
    else
    {
        LOG_INFO << "StPxlFastSim - smearing random generator is using default seed" <<endm;
    }
    
    St_HitError *pxlTableSet = (St_HitError*)calib_db.Find("PixelHitError");
    
    if (!pxlTableSet)
    {
        LOG_ERROR << "StPxlFastSim - E - PixelHitError is not available" << endm;
        return kStErr;
    }
    
    HitError_st* pxlHitError = pxlTableSet->GetTable();
    
    if (!pxlHitError)
    {
        LOG_ERROR << "StPxlFastSim - E - pxl hit table is not available in PixelHitError" << endm;
        return kStErr;
    }
    
    // please note that what is called local Y in the PXL sensor design
    // is actually called Z in STAR coordinates convention
    if (pxlHitError->coeff[0] <= 0 || pxlHitError->coeff[3] <= 0)
    {
        LOG_ERROR << "StPxlFastSim - E - negative or corrupted PXL hits errors in DB" << endm;
        return kStErr;
    }
    
    mResXPix = sqrt(pxlHitError->coeff[0]); // local x
    mResZPix = sqrt(pxlHitError->coeff[3]); // local Y
    mResYPix = 0;//sqrt(pxlHitError->coeff[2]); // needs to be updated in the DB later
    
    return kStOk;
}

//____________________________________________________________
Int_t StPxlFastSim::addPxlRawHits(const StMcPxlHitCollection& mcPxlHitCol,
                                  StPxlRawHitCollection& pxlRawHitCol)
{
    Float_t smearedX = 0, smearedY = 0, smearedZ = 0;
    // Loop over sectors
    for (UInt_t iSec = 0; iSec < mcPxlHitCol.numberOfSectors(); iSec++)
    {
        const StMcPxlSectorHitCollection* mcPxlSectorHitCol = mcPxlHitCol.sector(iSec);
        if (!mcPxlSectorHitCol) continue;
        
        for (UInt_t iLad = 0; iLad < mcPxlSectorHitCol->numberOfLadders(); iLad++)
        {
            const StMcPxlLadderHitCollection* mcPxlLadderHitCol = mcPxlSectorHitCol->ladder(iLad);
            if (!mcPxlLadderHitCol) continue;
            
            for (UInt_t iSen = 0; iSen < mcPxlLadderHitCol->numberOfSensors(); iSen++)
            {
                const StMcPxlSensorHitCollection* mcPxlSensorHitCol = mcPxlLadderHitCol->sensor(iSen);
                if (!mcPxlSensorHitCol) continue;
                
                UInt_t nSenHits = mcPxlSensorHitCol->hits().size();
                LOG_DEBUG << "Sector/Ladder/Sensor = " << iSec + 1 << "/" << iLad + 1 << "/" << iSen + 1 << ". Number of sensor hits = " << nSenHits << endm;
                
                // Loop over hits in the sensor
                for (UInt_t iHit = 0; iHit < nSenHits; iHit++)
                {
                    StMcPxlHit* mcPix = mcPxlSensorHitCol->hits()[iHit];
                    
                    Double_t localPixHitPos[3] = {mcPix->position().x(), mcPix->position().y(), mcPix->position().z()};

                    LOG_DEBUG << "localPixHitPos = " << localPixHitPos[0] << " " << localPixHitPos[1] << " " << localPixHitPos[2] << endm;
                    // please note that what is called local Y in the PXL sensor design
                    // is actually called Z in STAR coordinates convention and vice-versa
                    smearedX = distortHit(localPixHitPos[0], 10*0.1*0.001, StPxlConsts::kPxlActiveLengthX / 2.0); //
                    smearedZ = distortHit(localPixHitPos[2], 10*0.1*0.001, StPxlConsts::kPxlActiveLengthY / 2.0); //
                    smearedY = localPixHitPos[1];
                    
                    localPixHitPos[0] = smearedX;
                    localPixHitPos[2] = smearedZ;
                    localPixHitPos[1] = smearedY;
                    LOG_DEBUG << "smearedlocal = " << localPixHitPos[0] << " " << localPixHitPos[1] << " " << localPixHitPos[2] << endm;
                    
                    unsigned short idTruth = mcPix->parentTrack() ? mcPix->parentTrack()->key() : -999;

                    pxlRawHitCol.addRawHit(*addRawHit(localPixHitPos[0],localPixHitPos[2], iSec + 1, iLad + 1, iSen + 1, idTruth, 0));

                    cout << "i " << iSec + 1 << " " << iLad + 1 << " " << iSen + 1 << endl;
                    cout << "m " << mcPix->sector() << " " << mcPix->ladder() << " " << mcPix->sensor() << endl;
                    /*
                     Int_t row2=row;
                    Int_t column2=column;
                    Float_t smallR = 99999;
                    for (int iRow=-1; iRow<2; iRow++) {
                        if (row+iRow < 0) continue;
                        if (row+iRow > StPxlConsts::kPxlNumRowsPerSensor - 1) continue;
                        for (int iColumn=-1; iColumn<2; iColumn++) {
                            if (iRow == 0 && iColumn==0) continue;
                            if (column+iColumn < 0) continue;
                            if (column+iColumn > StPxlConsts::kPxlNumColumnsPerSensor - 1) continue;
                            Float_t r = (getLocalX(row+iRow)-localPixHitPos[0])*(getLocalX(row+iRow)-localPixHitPos[0]) + (getLocalZ(column+iColumn)-localPixHitPos[2])*(getLocalZ(column+iColumn)-localPixHitPos[2]);
                            if (r < smallR) {
                                smallR = r;
                                row2=row+iRow;
                                column2=column+iColumn;
                                cout << " " << row2 << " " << column2 << endl;

                            }
                        }
                    }
                    if (row==row2 && column==column2) continue;
                    StPxlRawHit* tempHit2;
                    tempHit2->setSector(iSec+1);
                    tempHit2->setLadder(mcPix->ladder());
                    tempHit2->setSensor(mcPix->sensor());
                    tempHit2->setRow(row2);
                    tempHit2->setColumn(column2);
                    tempHit2->setIdTruth(idTruth);
                    pxlRawHitCol.addRawHit(*tempHit2);
                    
                    
                    Int_t row3=row;
                    Int_t column3=column;
                    Float_t smallR = 99999;
                    for (int iRow=-1; iRow<2; iRow++) {
                        if (row+iRow < 0) continue;
                        if (row+iRow > StPxlConsts::kPxlNumRowsPerSensor - 1) continue;
                        for (int iColumn=-1; iColumn<2; iColumn++) {
                            if (iRow == 0 && iColumn==0) continue;
                            if (iRow+row == row2 && iColumn+column == column2) continue;
                            if (column+iColumn < 0) continue;
                            if (column+iColumn > StPxlConsts::kPxlNumColumnsPerSensor - 1) continue;
                            Float_t r = (getLocalX(row+iRow)-localPixHitPos[0])*(getLocalX(row+iRow)-localPixHitPos[0]) + (getLocalZ(column+iColumn)-localPixHitPos[2])*(getLocalZ(column+iColumn)-localPixHitPos[2]);
                            if (r < smallR) {
                                smallR = r;
                                row3=row+iRow;
                                column3=column+iColumn;
                            }
                        }
                    }
                    if (row==row3 && column==column3) continue;
                    StPxlRawHit* tempHit3;
                    tempHit3->setSector(iSec+1);
                    tempHit3->setLadder(mcPix->ladder());
                    tempHit3->setSensor(mcPix->sensor());
                    tempHit3->setRow(row3);
                    tempHit3->setColumn(column3);
                    tempHit3->setIdTruth(idTruth);
                    pxlRawHitCol.addRawHit(*tempHit3);
                    
                    */
                    
                }
            }
        }
    }
    
    return kStOK;
    
    
}
//____________________________________________________________
Int_t StPxlFastSim::getRow(float local){
    const double mFirstPixelX =  (StPxlConsts::kPxlNumRowsPerSensor - 1) * StPxlConsts::kPixelSize / 2;
    return -(local - mFirstPixelX) / StPxlConsts::kPixelSize;
}
Int_t StPxlFastSim::getColumn(float local){
    const double mFirstPixelZ = -(StPxlConsts::kPxlNumColumnsPerSensor - 1) * StPxlConsts::kPixelSize / 2;
    return (local - mFirstPixelZ) / StPxlConsts::kPixelSize;
}
Float_t StPxlFastSim::getLocalX(Int_t value){
    const double mFirstPixelX =  (StPxlConsts::kPxlNumRowsPerSensor - 1) * StPxlConsts::kPixelSize / 2;
    return mFirstPixelX - StPxlConsts::kPixelSize * value;
}
Float_t StPxlFastSim::getLocalZ(Int_t value){
    const double mFirstPixelZ = -(StPxlConsts::kPxlNumColumnsPerSensor - 1) * StPxlConsts::kPixelSize / 2;
    return mFirstPixelZ + StPxlConsts::kPixelSize * value;
}
StPxlRawHit * StPxlFastSim::addRawHit(float localX, float localZ, int iSec, int iLad, int iSen, int idTruth, int nHits=0){
    Int_t row = getRow(localX);
    Int_t column = getColumn(localZ);
    
    tempHit->setSector(iSec);
    tempHit->setLadder(iLad);
    tempHit->setSensor(iSen);
    tempHit->setRow(row);
    tempHit->setColumn(column);
    tempHit->setIdTruth(idTruth);
    return tempHit;

}

//____________________________________________________________
Int_t StPxlFastSim::addPxlHits(const StMcPxlHitCollection& mcPxlHitCol,
                               StPxlHitCollection& pxlHitCol)
{
    
    Float_t smearedX = 0, smearedY = 0, smearedZ = 0;
    
    
    // Loop over sectors
    for (UInt_t iSec = 0; iSec < mcPxlHitCol.numberOfSectors(); iSec++)
    {
        const StMcPxlSectorHitCollection* mcPxlSectorHitCol = mcPxlHitCol.sector(iSec);
        if (!mcPxlSectorHitCol) continue;
        
        for (UInt_t iLad = 0; iLad < mcPxlSectorHitCol->numberOfLadders(); iLad++)
        {
            const StMcPxlLadderHitCollection* mcPxlLadderHitCol = mcPxlSectorHitCol->ladder(iLad);
            if (!mcPxlLadderHitCol) continue;
            
            for (UInt_t iSen = 0; iSen < mcPxlLadderHitCol->numberOfSensors(); iSen++)
            {
                const StMcPxlSensorHitCollection* mcPxlSensorHitCol = mcPxlLadderHitCol->sensor(iSen);
                if (!mcPxlSensorHitCol) continue;
                
                UInt_t nSenHits = mcPxlSensorHitCol->hits().size();
                LOG_DEBUG << "Sector/Ladder/Sensor = " << iSec + 1 << "/" << iLad + 1 << "/" << iSen + 1 << ". Number of sensor hits = " << nSenHits << endm;
                
                // Loop over hits in the sensor
                for (UInt_t iHit = 0; iHit < nSenHits; iHit++)
                {
                    StMcPxlHit* mcPix = mcPxlSensorHitCol->hits()[iHit];
                    
                    //Long_t volId = mcPix->volumeId();
                    Int_t sector = mcPix->sector();
                    Int_t ladder = mcPix->ladder();
                    //Int_t sensor = mcPix->sensor();
                    
                    Double_t localPixHitPos[3] = {mcPix->position().x(), mcPix->position().y(), mcPix->position().z()};

                    LOG_DEBUG << "localPixHitPos = " << localPixHitPos[0] << " " << localPixHitPos[1] << " " << localPixHitPos[2] << endm;
                    // please note that what is called local Y in the PXL sensor design
                    // is actually called Z in STAR coordinates convention and vice-versa
                    smearedX = distortHit(localPixHitPos[0], mResXPix, StPxlConsts::kPxlActiveLengthX / 2.0);
                    smearedZ = distortHit(localPixHitPos[2], mResZPix, StPxlConsts::kPxlActiveLengthY / 2.0);
                    if (mResYPix) smearedY = distortHit(localPixHitPos[1], mResYPix, 0.0020); // Not properly constrained yet
                    else smearedY = localPixHitPos[1];
                    localPixHitPos[0] = smearedX;
                    localPixHitPos[2] = smearedZ;
                    localPixHitPos[1] = smearedY;
                    LOG_DEBUG << "smearedlocal = " << localPixHitPos[0] << " " << localPixHitPos[1] << " " << localPixHitPos[2] << endm;
                    Double_t smearedGlobalPixHitPos[3] = {0, 0, 0};
                    localToMatser(localPixHitPos,smearedGlobalPixHitPos,iSec+1,iLad+1,iSen+1);
                    
                    StThreeVectorF gpixpos(smearedGlobalPixHitPos);
                    StThreeVectorF mRndHitError(0., 0., 0.);
                    
                    UInt_t hw = sector * 10 + ladder; // needs to be updated later after clustering alogrithms are finalized
                    
                    unsigned short idTruth = mcPix->parentTrack() ? mcPix->parentTrack()->key() : -999;
                    unsigned short quality = mcPix->parentTrack() ? 100 : 0;
                    
                    StPxlDigiHit* tempHit = new StPxlDigiHit(localPixHitPos, iSec+1, mcPix->ladder(), mcPix->sensor(),
                                                             gpixpos, mRndHitError, hw, mcPix->dE(), 0, idTruth, quality, mcPix->key());
                    
                    LOG_DEBUG << "key() : " << mcPix->key() - 1 << " idTruth: " << mcPix->parentTrack()->key() << endm;
                    LOG_DEBUG << "from StMcPxlHit : x= " << mcPix->position().x() << ";  y= " << mcPix->position().y() << ";  z= " << mcPix->position().z() << endm;
                    LOG_DEBUG << "pxlHit location x= " << tempHit->position().x() << "; y= " << tempHit->position().y() << "; z= " << tempHit->position().z() << endm;
                    
                    pxlHitCol.addHit(tempHit);
                }
            }
        }
    }
    
    return kStOK;
}


/**
 * Calculates and returns new value for the local coordinate x by smearing it
 * acccording to a normal distribution N(mean, sigma) = N(x, res). The returned
 * value is constrained to be within the characteristic dimension detLength
 * provided by the user.
 */
double StPxlFastSim::distortHit(const double x, const double res, const double detLength) const
{
    // Do not smear x when it is outside the physical limits. Issue a warning instead
    if (fabs(x) > detLength) {
        LOG_WARN << "distortHit() - Generated hit is outside detector sensor plane" << endm;
        return x;
    }
    
    double smeared_x;
    
    do {
        smeared_x = mRandom->gauss(x, res);
    } while ( fabs(smeared_x) > detLength);
    
    return smeared_x;
}


//____________________________________________________________
void StPxlFastSim::localToMatser(Double_t* local,Double_t* master,Int_t sector,Int_t ladder,Int_t sensor)
{
    if(mPxlDb)
    {
        TGeoHMatrix *combP = (TGeoHMatrix *)mPxlDb->geoHMatrixSensorOnGlobal(sector, ladder,sensor);
        combP->LocalToMaster(local, master);
    }
    else
    {
        TString Path("");
        LOG_DEBUG << endm;
        Path = Form("/HALL_1/CAVE_1/TpcRefSys_1/IDSM_1/PXMO_1/PXLA_%i/LADR_%i/PXSI_%i/PLAC_1", sector, ladder, sensor);
        
        gGeoManager->RestoreMasterVolume();
        gGeoManager->CdTop();
        gGeoManager->cd(Path);
        
        gGeoManager->GetCurrentMatrix()->LocalToMaster(local, master);
    }
}
/*
 *
 * Author: M. Mustafa
 *
 *
 **********************************************************
 * $Log: StPxlFastSim.cxx,v $
 * Revision 1.12  2015/05/14 18:57:52  smirnovd
 * Squashed commit of the following:
 *
 * StPxlFastSim: Streamlined creation of PXL hits by making use of StPxlUtil/StPxlDigiHit
 *
 * StPxlHitMaker: Updated comments
 *
 * StPxlHitMaker: Streamlined creation of PXL hits by making use of StPxlUtil/StPxlDigiHit
 *
 * StPxlDigiHit: A helper to manipulate local hit position in StPxlHit
 *
 * StPxlConsts: Define constants in namespace
 *
 * For safety reasons, the intentions is to move the constants into the namespace
 * and get rid of those defined in the global space.
 *
 * Revision 1.11  2015/05/07 21:24:31  smirnovd
 * StPxlSimMaker: Switched to using consts from StPxlUtil/
 *
 * Revision 1.10  2015/03/13 18:45:01  perev
 * Roll back
 *
 * Revision 1.8  2015/01/27 19:11:49  mstftsm
 * Set idTruth of StPxlHit to -999 if parentTrack of mcHit does not exist (for protection).
 *
 * Revision 1.7  2015/01/27 01:31:09  smirnovd
 * Minor refactoring of StPxlFastSim::distortHit() to include a new warning for unphysical hit position
 *
 * Revision 1.6  2014/07/03 19:46:37  mstftsm
 * Revereted the changes made for the pileup adder. That does not belong to the master branch.
 *
 * Revision 1.4  2014/03/13 17:00:19  mstftsm
 * StPxlSimMaker has a method to switch on random seed for StRandom generatos in simulators. Default is not a random seed.
 *
 * Revision 1.2  2013/11/14 19:10:27  mstftsm
 * StMcPxlHit has been changed to be on local coordinates. We no longer transfor from global to local before smearing
 *
 * Revision 1.1  2013/05/12 21:43:32  jeromel
 * Initial revision, code peer review closed 2013/05/06
 *
 * Revision 1.5  2013/05/09 02:58:36  mstftsm
 * Fixed a bug which called for sensor local Z value.
 *
 */

