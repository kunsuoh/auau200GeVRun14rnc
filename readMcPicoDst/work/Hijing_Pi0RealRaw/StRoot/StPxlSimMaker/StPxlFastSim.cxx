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
#include "StPxlRawHit.h"
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
#include "TH1F.h"

StPxlFastSim::~StPxlFastSim()
{
    if (mRandom) delete mRandom;
    if (mPxlDb) delete mPxlDb;
    if (dataH2) delete dataH2;
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
    
    // Cluster size from 20 to 30 degrees in cosmic ray
    dataH2 = new TH1F("dataH2","Cluster multiplicity from 20 to 30 degrees",14,1,15);
    dataH2->SetBinContent(1,0.09685086);
    dataH2->SetBinContent(2,0.1921569);
    dataH2->SetBinContent(3,0.1858586);
    dataH2->SetBinContent(4,0.332145);
    dataH2->SetBinContent(5,0.07213309);
    dataH2->SetBinContent(6,0.05692216);
    dataH2->SetBinContent(7,0.02507427);
    dataH2->SetBinContent(8,0.01283422);
    dataH2->SetBinContent(9,0.008437314);
    dataH2->SetBinContent(10,0.004634581);
    dataH2->SetBinContent(11,0.004872252);
    dataH2->SetBinContent(12,0.004040404);
    dataH2->SetBinContent(13,0.002614379);
    dataH2->SetBinContent(14,0.001426025);
    dataH2->SetBinContent(15,0.002257873);
    dataH2->SetBinError(1,0.003392538);
    dataH2->SetBinError(2,0.004778602);
    dataH2->SetBinError(3,0.004699636);
    dataH2->SetBinError(4,0.006282562);
    dataH2->SetBinError(5,0.002927792);
    dataH2->SetBinError(6,0.00260084);
    dataH2->SetBinError(7,0.001726184);
    dataH2->SetBinError(8,0.001234974);
    dataH2->SetBinError(9,0.001001325);
    dataH2->SetBinError(10,0.0007421269);
    dataH2->SetBinError(11,0.0007609179);
    dataH2->SetBinError(12,0.0006929236);
    dataH2->SetBinError(13,0.0005573875);
    dataH2->SetBinError(14,0.0004116579);
    dataH2->SetBinError(15,0.0005179916);
    dataH2->SetEntries(8434);

histRandomHitDistribution = new TH2F("histRandomHitDistribution","",100,0,10000,100,0,2);
histRandomHitDistribution->SetBinContent(103,5);
histRandomHitDistribution->SetBinContent(104,135);
histRandomHitDistribution->SetBinContent(105,290);
histRandomHitDistribution->SetBinContent(106,409);
histRandomHitDistribution->SetBinContent(107,366);
histRandomHitDistribution->SetBinContent(108,406);
histRandomHitDistribution->SetBinContent(109,286);
histRandomHitDistribution->SetBinContent(110,200);
histRandomHitDistribution->SetBinContent(111,99);
histRandomHitDistribution->SetBinContent(112,61);
histRandomHitDistribution->SetBinContent(113,25);
histRandomHitDistribution->SetBinContent(114,7);
histRandomHitDistribution->SetBinContent(115,4);
histRandomHitDistribution->SetBinContent(117,1);
histRandomHitDistribution->SetBinContent(122,1);
histRandomHitDistribution->SetBinContent(126,1);
histRandomHitDistribution->SetBinContent(133,1);
histRandomHitDistribution->SetBinContent(135,1);
histRandomHitDistribution->SetBinContent(136,2);
histRandomHitDistribution->SetBinContent(139,1);
histRandomHitDistribution->SetBinContent(140,1);
histRandomHitDistribution->SetBinContent(205,3);
histRandomHitDistribution->SetBinContent(206,16);
histRandomHitDistribution->SetBinContent(207,149);
histRandomHitDistribution->SetBinContent(208,262);
histRandomHitDistribution->SetBinContent(209,289);
histRandomHitDistribution->SetBinContent(210,212);
histRandomHitDistribution->SetBinContent(211,152);
histRandomHitDistribution->SetBinContent(212,109);
histRandomHitDistribution->SetBinContent(213,68);
histRandomHitDistribution->SetBinContent(214,32);
histRandomHitDistribution->SetBinContent(215,10);
histRandomHitDistribution->SetBinContent(216,8);
histRandomHitDistribution->SetBinContent(217,1);
histRandomHitDistribution->SetBinContent(218,2);
histRandomHitDistribution->SetBinContent(244,1);
histRandomHitDistribution->SetBinContent(308,34);
histRandomHitDistribution->SetBinContent(309,120);
histRandomHitDistribution->SetBinContent(310,256);
histRandomHitDistribution->SetBinContent(311,275);
histRandomHitDistribution->SetBinContent(312,186);
histRandomHitDistribution->SetBinContent(313,129);
histRandomHitDistribution->SetBinContent(314,85);
histRandomHitDistribution->SetBinContent(315,49);
histRandomHitDistribution->SetBinContent(316,18);
histRandomHitDistribution->SetBinContent(317,11);
histRandomHitDistribution->SetBinContent(318,6);
histRandomHitDistribution->SetBinContent(319,3);
histRandomHitDistribution->SetBinContent(337,1);
histRandomHitDistribution->SetBinContent(339,1);
histRandomHitDistribution->SetBinContent(343,1);
histRandomHitDistribution->SetBinContent(346,1);
histRandomHitDistribution->SetBinContent(409,1);
histRandomHitDistribution->SetBinContent(410,25);
histRandomHitDistribution->SetBinContent(411,140);
histRandomHitDistribution->SetBinContent(412,220);
histRandomHitDistribution->SetBinContent(413,234);
histRandomHitDistribution->SetBinContent(414,141);
histRandomHitDistribution->SetBinContent(415,82);
histRandomHitDistribution->SetBinContent(416,76);
histRandomHitDistribution->SetBinContent(417,41);
histRandomHitDistribution->SetBinContent(418,20);
histRandomHitDistribution->SetBinContent(419,6);
histRandomHitDistribution->SetBinContent(420,3);
histRandomHitDistribution->SetBinContent(421,4);
histRandomHitDistribution->SetBinContent(422,2);
histRandomHitDistribution->SetBinContent(439,1);
histRandomHitDistribution->SetBinContent(442,1);
histRandomHitDistribution->SetBinContent(451,1);
histRandomHitDistribution->SetBinContent(512,25);
histRandomHitDistribution->SetBinContent(513,103);
histRandomHitDistribution->SetBinContent(514,209);
histRandomHitDistribution->SetBinContent(515,177);
histRandomHitDistribution->SetBinContent(516,111);
histRandomHitDistribution->SetBinContent(517,96);
histRandomHitDistribution->SetBinContent(518,75);
histRandomHitDistribution->SetBinContent(519,30);
histRandomHitDistribution->SetBinContent(520,21);
histRandomHitDistribution->SetBinContent(521,17);
histRandomHitDistribution->SetBinContent(522,6);
histRandomHitDistribution->SetBinContent(523,5);
histRandomHitDistribution->SetBinContent(524,1);
histRandomHitDistribution->SetBinContent(536,1);
histRandomHitDistribution->SetBinContent(539,1);
histRandomHitDistribution->SetBinContent(540,1);
histRandomHitDistribution->SetBinContent(541,1);
histRandomHitDistribution->SetBinContent(544,1);
histRandomHitDistribution->SetBinContent(545,1);
histRandomHitDistribution->SetBinContent(548,1);
histRandomHitDistribution->SetBinContent(550,1);
histRandomHitDistribution->SetBinContent(613,1);
histRandomHitDistribution->SetBinContent(614,32);
histRandomHitDistribution->SetBinContent(615,87);
histRandomHitDistribution->SetBinContent(616,177);
histRandomHitDistribution->SetBinContent(617,208);
histRandomHitDistribution->SetBinContent(618,145);
histRandomHitDistribution->SetBinContent(619,84);
histRandomHitDistribution->SetBinContent(620,75);
histRandomHitDistribution->SetBinContent(621,48);
histRandomHitDistribution->SetBinContent(622,17);
histRandomHitDistribution->SetBinContent(623,10);
histRandomHitDistribution->SetBinContent(624,3);
histRandomHitDistribution->SetBinContent(625,1);
histRandomHitDistribution->SetBinContent(626,3);
histRandomHitDistribution->SetBinContent(627,3);
histRandomHitDistribution->SetBinContent(628,1);
histRandomHitDistribution->SetBinContent(636,1);
histRandomHitDistribution->SetBinContent(640,2);
histRandomHitDistribution->SetBinContent(643,1);
histRandomHitDistribution->SetBinContent(652,1);
histRandomHitDistribution->SetBinContent(656,2);
histRandomHitDistribution->SetBinContent(716,24);
histRandomHitDistribution->SetBinContent(717,80);
histRandomHitDistribution->SetBinContent(718,152);
histRandomHitDistribution->SetBinContent(719,199);
histRandomHitDistribution->SetBinContent(720,141);
histRandomHitDistribution->SetBinContent(721,84);
histRandomHitDistribution->SetBinContent(722,63);
histRandomHitDistribution->SetBinContent(723,40);
histRandomHitDistribution->SetBinContent(724,16);
histRandomHitDistribution->SetBinContent(725,10);
histRandomHitDistribution->SetBinContent(726,6);
histRandomHitDistribution->SetBinContent(727,2);
histRandomHitDistribution->SetBinContent(728,2);
histRandomHitDistribution->SetBinContent(729,1);
histRandomHitDistribution->SetBinContent(730,1);
histRandomHitDistribution->SetBinContent(740,1);
histRandomHitDistribution->SetBinContent(742,1);
histRandomHitDistribution->SetBinContent(748,1);
histRandomHitDistribution->SetBinContent(750,2);
histRandomHitDistribution->SetBinContent(755,2);
histRandomHitDistribution->SetBinContent(758,1);
histRandomHitDistribution->SetBinContent(762,1);
histRandomHitDistribution->SetBinContent(818,23);
histRandomHitDistribution->SetBinContent(819,76);
histRandomHitDistribution->SetBinContent(820,136);
histRandomHitDistribution->SetBinContent(821,204);
histRandomHitDistribution->SetBinContent(822,171);
histRandomHitDistribution->SetBinContent(823,137);
histRandomHitDistribution->SetBinContent(824,65);
histRandomHitDistribution->SetBinContent(825,32);
histRandomHitDistribution->SetBinContent(826,27);
histRandomHitDistribution->SetBinContent(827,12);
histRandomHitDistribution->SetBinContent(828,10);
histRandomHitDistribution->SetBinContent(829,3);
histRandomHitDistribution->SetBinContent(835,1);
histRandomHitDistribution->SetBinContent(837,1);
histRandomHitDistribution->SetBinContent(838,1);
histRandomHitDistribution->SetBinContent(839,1);
histRandomHitDistribution->SetBinContent(843,1);
histRandomHitDistribution->SetBinContent(844,1);
histRandomHitDistribution->SetBinContent(845,2);
histRandomHitDistribution->SetBinContent(846,1);
histRandomHitDistribution->SetBinContent(847,1);
histRandomHitDistribution->SetBinContent(848,2);
histRandomHitDistribution->SetBinContent(849,1);
histRandomHitDistribution->SetBinContent(850,2);
histRandomHitDistribution->SetBinContent(851,2);
histRandomHitDistribution->SetBinContent(852,1);
histRandomHitDistribution->SetBinContent(854,1);
histRandomHitDistribution->SetBinContent(855,1);
histRandomHitDistribution->SetBinContent(857,2);
histRandomHitDistribution->SetBinContent(860,2);
histRandomHitDistribution->SetBinContent(861,2);
histRandomHitDistribution->SetBinContent(863,2);
histRandomHitDistribution->SetBinContent(865,1);
histRandomHitDistribution->SetBinContent(866,1);
histRandomHitDistribution->SetBinContent(919,1);
histRandomHitDistribution->SetBinContent(920,16);
histRandomHitDistribution->SetBinContent(921,86);
histRandomHitDistribution->SetBinContent(922,152);
histRandomHitDistribution->SetBinContent(923,200);
histRandomHitDistribution->SetBinContent(924,149);
histRandomHitDistribution->SetBinContent(925,102);
histRandomHitDistribution->SetBinContent(926,88);
histRandomHitDistribution->SetBinContent(927,56);
histRandomHitDistribution->SetBinContent(928,27);
histRandomHitDistribution->SetBinContent(929,12);
histRandomHitDistribution->SetBinContent(930,5);
histRandomHitDistribution->SetBinContent(931,8);
histRandomHitDistribution->SetBinContent(932,4);
histRandomHitDistribution->SetBinContent(935,1);
histRandomHitDistribution->SetBinContent(936,1);
histRandomHitDistribution->SetBinContent(943,1);
histRandomHitDistribution->SetBinContent(944,1);
histRandomHitDistribution->SetBinContent(947,1);
histRandomHitDistribution->SetBinContent(948,2);
histRandomHitDistribution->SetBinContent(950,1);
histRandomHitDistribution->SetBinContent(951,3);
histRandomHitDistribution->SetBinContent(952,1);
histRandomHitDistribution->SetBinContent(957,1);
histRandomHitDistribution->SetBinContent(958,1);
histRandomHitDistribution->SetBinContent(961,1);
histRandomHitDistribution->SetBinContent(965,1);
histRandomHitDistribution->SetBinContent(1021,2);
histRandomHitDistribution->SetBinContent(1022,12);
histRandomHitDistribution->SetBinContent(1023,61);
histRandomHitDistribution->SetBinContent(1024,155);
histRandomHitDistribution->SetBinContent(1025,194);
histRandomHitDistribution->SetBinContent(1026,171);
histRandomHitDistribution->SetBinContent(1027,129);
histRandomHitDistribution->SetBinContent(1028,94);
histRandomHitDistribution->SetBinContent(1029,48);
histRandomHitDistribution->SetBinContent(1030,24);
histRandomHitDistribution->SetBinContent(1031,20);
histRandomHitDistribution->SetBinContent(1032,6);
histRandomHitDistribution->SetBinContent(1033,6);
histRandomHitDistribution->SetBinContent(1034,1);
histRandomHitDistribution->SetBinContent(1040,1);
histRandomHitDistribution->SetBinContent(1044,1);
histRandomHitDistribution->SetBinContent(1050,3);
histRandomHitDistribution->SetBinContent(1051,1);
histRandomHitDistribution->SetBinContent(1053,2);
histRandomHitDistribution->SetBinContent(1055,3);
histRandomHitDistribution->SetBinContent(1056,3);
histRandomHitDistribution->SetBinContent(1057,2);
histRandomHitDistribution->SetBinContent(1058,3);
histRandomHitDistribution->SetBinContent(1061,1);
histRandomHitDistribution->SetBinContent(1065,1);
histRandomHitDistribution->SetBinContent(1066,2);
histRandomHitDistribution->SetBinContent(1068,3);
histRandomHitDistribution->SetBinContent(1123,1);
histRandomHitDistribution->SetBinContent(1124,21);
histRandomHitDistribution->SetBinContent(1125,85);
histRandomHitDistribution->SetBinContent(1126,171);
histRandomHitDistribution->SetBinContent(1127,197);
histRandomHitDistribution->SetBinContent(1128,191);
histRandomHitDistribution->SetBinContent(1129,137);
histRandomHitDistribution->SetBinContent(1130,83);
histRandomHitDistribution->SetBinContent(1131,51);
histRandomHitDistribution->SetBinContent(1132,38);
histRandomHitDistribution->SetBinContent(1133,22);
histRandomHitDistribution->SetBinContent(1134,11);
histRandomHitDistribution->SetBinContent(1136,3);
histRandomHitDistribution->SetBinContent(1137,1);
histRandomHitDistribution->SetBinContent(1138,1);
histRandomHitDistribution->SetBinContent(1142,1);
histRandomHitDistribution->SetBinContent(1143,1);
histRandomHitDistribution->SetBinContent(1148,4);
histRandomHitDistribution->SetBinContent(1150,1);
histRandomHitDistribution->SetBinContent(1151,1);
histRandomHitDistribution->SetBinContent(1153,2);
histRandomHitDistribution->SetBinContent(1155,3);
histRandomHitDistribution->SetBinContent(1156,1);
histRandomHitDistribution->SetBinContent(1157,2);
histRandomHitDistribution->SetBinContent(1158,2);
histRandomHitDistribution->SetBinContent(1159,2);
histRandomHitDistribution->SetBinContent(1160,1);
histRandomHitDistribution->SetBinContent(1161,2);
histRandomHitDistribution->SetBinContent(1164,2);
histRandomHitDistribution->SetBinContent(1168,1);
histRandomHitDistribution->SetBinContent(1171,2);
histRandomHitDistribution->SetBinContent(1172,4);
histRandomHitDistribution->SetBinContent(1175,1);
histRandomHitDistribution->SetBinContent(1179,1);
histRandomHitDistribution->SetBinContent(1181,1);
histRandomHitDistribution->SetBinContent(1226,17);
histRandomHitDistribution->SetBinContent(1227,88);
histRandomHitDistribution->SetBinContent(1228,141);
histRandomHitDistribution->SetBinContent(1229,190);
histRandomHitDistribution->SetBinContent(1230,186);
histRandomHitDistribution->SetBinContent(1231,175);
histRandomHitDistribution->SetBinContent(1232,91);
histRandomHitDistribution->SetBinContent(1233,55);
histRandomHitDistribution->SetBinContent(1234,38);
histRandomHitDistribution->SetBinContent(1235,24);
histRandomHitDistribution->SetBinContent(1236,16);
histRandomHitDistribution->SetBinContent(1237,12);
histRandomHitDistribution->SetBinContent(1238,2);
histRandomHitDistribution->SetBinContent(1239,3);
histRandomHitDistribution->SetBinContent(1245,2);
histRandomHitDistribution->SetBinContent(1249,2);
histRandomHitDistribution->SetBinContent(1250,1);
histRandomHitDistribution->SetBinContent(1251,3);
histRandomHitDistribution->SetBinContent(1252,1);
histRandomHitDistribution->SetBinContent(1253,3);
histRandomHitDistribution->SetBinContent(1254,2);
histRandomHitDistribution->SetBinContent(1255,5);
histRandomHitDistribution->SetBinContent(1256,1);
histRandomHitDistribution->SetBinContent(1257,1);
histRandomHitDistribution->SetBinContent(1258,2);
histRandomHitDistribution->SetBinContent(1259,1);
histRandomHitDistribution->SetBinContent(1260,3);
histRandomHitDistribution->SetBinContent(1261,1);
histRandomHitDistribution->SetBinContent(1262,5);
histRandomHitDistribution->SetBinContent(1263,3);
histRandomHitDistribution->SetBinContent(1264,1);
histRandomHitDistribution->SetBinContent(1266,4);
histRandomHitDistribution->SetBinContent(1268,2);
histRandomHitDistribution->SetBinContent(1269,1);
histRandomHitDistribution->SetBinContent(1270,1);
histRandomHitDistribution->SetBinContent(1273,1);
histRandomHitDistribution->SetBinContent(1274,1);
histRandomHitDistribution->SetBinContent(1277,1);
histRandomHitDistribution->SetBinContent(1279,1);
histRandomHitDistribution->SetBinContent(1281,1);
histRandomHitDistribution->SetBinContent(1285,1);
histRandomHitDistribution->SetBinContent(1327,1);
histRandomHitDistribution->SetBinContent(1328,13);
histRandomHitDistribution->SetBinContent(1329,89);
histRandomHitDistribution->SetBinContent(1330,168);
histRandomHitDistribution->SetBinContent(1331,187);
histRandomHitDistribution->SetBinContent(1332,210);
histRandomHitDistribution->SetBinContent(1333,176);
histRandomHitDistribution->SetBinContent(1334,115);
histRandomHitDistribution->SetBinContent(1335,79);
histRandomHitDistribution->SetBinContent(1336,41);
histRandomHitDistribution->SetBinContent(1337,20);
histRandomHitDistribution->SetBinContent(1338,7);
histRandomHitDistribution->SetBinContent(1339,5);
histRandomHitDistribution->SetBinContent(1340,2);
histRandomHitDistribution->SetBinContent(1341,1);
histRandomHitDistribution->SetBinContent(1343,1);
histRandomHitDistribution->SetBinContent(1346,1);
histRandomHitDistribution->SetBinContent(1347,1);
histRandomHitDistribution->SetBinContent(1348,2);
histRandomHitDistribution->SetBinContent(1349,1);
histRandomHitDistribution->SetBinContent(1350,1);
histRandomHitDistribution->SetBinContent(1351,4);
histRandomHitDistribution->SetBinContent(1352,1);
histRandomHitDistribution->SetBinContent(1353,1);
histRandomHitDistribution->SetBinContent(1354,1);
histRandomHitDistribution->SetBinContent(1355,3);
histRandomHitDistribution->SetBinContent(1356,4);
histRandomHitDistribution->SetBinContent(1357,3);
histRandomHitDistribution->SetBinContent(1358,3);
histRandomHitDistribution->SetBinContent(1359,3);
histRandomHitDistribution->SetBinContent(1360,6);
histRandomHitDistribution->SetBinContent(1361,2);
histRandomHitDistribution->SetBinContent(1362,6);
histRandomHitDistribution->SetBinContent(1363,1);
histRandomHitDistribution->SetBinContent(1364,1);
histRandomHitDistribution->SetBinContent(1365,1);
histRandomHitDistribution->SetBinContent(1366,2);
histRandomHitDistribution->SetBinContent(1367,1);
histRandomHitDistribution->SetBinContent(1368,2);
histRandomHitDistribution->SetBinContent(1370,1);
histRandomHitDistribution->SetBinContent(1371,2);
histRandomHitDistribution->SetBinContent(1372,1);
histRandomHitDistribution->SetBinContent(1374,1);
histRandomHitDistribution->SetBinContent(1376,2);
histRandomHitDistribution->SetBinContent(1378,1);
histRandomHitDistribution->SetBinContent(1383,1);
histRandomHitDistribution->SetBinContent(1430,20);
histRandomHitDistribution->SetBinContent(1431,76);
histRandomHitDistribution->SetBinContent(1432,151);
histRandomHitDistribution->SetBinContent(1433,206);
histRandomHitDistribution->SetBinContent(1434,202);
histRandomHitDistribution->SetBinContent(1435,175);
histRandomHitDistribution->SetBinContent(1436,134);
histRandomHitDistribution->SetBinContent(1437,87);
histRandomHitDistribution->SetBinContent(1438,51);
histRandomHitDistribution->SetBinContent(1439,26);
histRandomHitDistribution->SetBinContent(1440,12);
histRandomHitDistribution->SetBinContent(1441,9);
histRandomHitDistribution->SetBinContent(1442,2);
histRandomHitDistribution->SetBinContent(1443,2);
histRandomHitDistribution->SetBinContent(1444,2);
histRandomHitDistribution->SetBinContent(1446,1);
histRandomHitDistribution->SetBinContent(1447,1);
histRandomHitDistribution->SetBinContent(1449,1);
histRandomHitDistribution->SetBinContent(1450,3);
histRandomHitDistribution->SetBinContent(1451,1);
histRandomHitDistribution->SetBinContent(1452,1);
histRandomHitDistribution->SetBinContent(1453,1);
histRandomHitDistribution->SetBinContent(1454,2);
histRandomHitDistribution->SetBinContent(1455,3);
histRandomHitDistribution->SetBinContent(1456,2);
histRandomHitDistribution->SetBinContent(1457,3);
histRandomHitDistribution->SetBinContent(1458,5);
histRandomHitDistribution->SetBinContent(1459,5);
histRandomHitDistribution->SetBinContent(1460,5);
histRandomHitDistribution->SetBinContent(1461,4);
histRandomHitDistribution->SetBinContent(1462,2);
histRandomHitDistribution->SetBinContent(1463,4);
histRandomHitDistribution->SetBinContent(1464,7);
histRandomHitDistribution->SetBinContent(1465,9);
histRandomHitDistribution->SetBinContent(1466,1);
histRandomHitDistribution->SetBinContent(1467,4);
histRandomHitDistribution->SetBinContent(1468,4);
histRandomHitDistribution->SetBinContent(1469,3);
histRandomHitDistribution->SetBinContent(1470,2);
histRandomHitDistribution->SetBinContent(1471,2);
histRandomHitDistribution->SetBinContent(1472,2);
histRandomHitDistribution->SetBinContent(1473,2);
histRandomHitDistribution->SetBinContent(1474,2);
histRandomHitDistribution->SetBinContent(1475,1);
histRandomHitDistribution->SetBinContent(1477,2);
histRandomHitDistribution->SetBinContent(1479,2);
histRandomHitDistribution->SetBinContent(1481,2);
histRandomHitDistribution->SetBinContent(1483,1);
histRandomHitDistribution->SetBinContent(1484,2);
histRandomHitDistribution->SetBinContent(1487,2);
histRandomHitDistribution->SetBinContent(1498,1);
histRandomHitDistribution->SetBinContent(1531,2);
histRandomHitDistribution->SetBinContent(1532,22);
histRandomHitDistribution->SetBinContent(1533,86);
histRandomHitDistribution->SetBinContent(1534,169);
histRandomHitDistribution->SetBinContent(1535,205);
histRandomHitDistribution->SetBinContent(1536,278);
histRandomHitDistribution->SetBinContent(1537,193);
histRandomHitDistribution->SetBinContent(1538,137);
histRandomHitDistribution->SetBinContent(1539,81);
histRandomHitDistribution->SetBinContent(1540,63);
histRandomHitDistribution->SetBinContent(1541,36);
histRandomHitDistribution->SetBinContent(1542,13);
histRandomHitDistribution->SetBinContent(1543,6);
histRandomHitDistribution->SetBinContent(1544,2);
histRandomHitDistribution->SetBinContent(1545,5);
histRandomHitDistribution->SetBinContent(1549,2);
histRandomHitDistribution->SetBinContent(1552,2);
histRandomHitDistribution->SetBinContent(1553,4);
histRandomHitDistribution->SetBinContent(1554,5);
histRandomHitDistribution->SetBinContent(1555,5);
histRandomHitDistribution->SetBinContent(1556,4);
histRandomHitDistribution->SetBinContent(1557,5);
histRandomHitDistribution->SetBinContent(1558,2);
histRandomHitDistribution->SetBinContent(1559,3);
histRandomHitDistribution->SetBinContent(1560,8);
histRandomHitDistribution->SetBinContent(1561,8);
histRandomHitDistribution->SetBinContent(1562,8);
histRandomHitDistribution->SetBinContent(1563,6);
histRandomHitDistribution->SetBinContent(1564,7);
histRandomHitDistribution->SetBinContent(1565,9);
histRandomHitDistribution->SetBinContent(1566,3);
histRandomHitDistribution->SetBinContent(1567,7);
histRandomHitDistribution->SetBinContent(1568,2);
histRandomHitDistribution->SetBinContent(1569,4);
histRandomHitDistribution->SetBinContent(1570,5);
histRandomHitDistribution->SetBinContent(1571,2);
histRandomHitDistribution->SetBinContent(1572,2);
histRandomHitDistribution->SetBinContent(1573,3);
histRandomHitDistribution->SetBinContent(1574,3);
histRandomHitDistribution->SetBinContent(1576,2);
histRandomHitDistribution->SetBinContent(1577,3);
histRandomHitDistribution->SetBinContent(1578,4);
histRandomHitDistribution->SetBinContent(1579,1);
histRandomHitDistribution->SetBinContent(1581,1);
histRandomHitDistribution->SetBinContent(1582,2);
histRandomHitDistribution->SetBinContent(1583,1);
histRandomHitDistribution->SetBinContent(1584,2);
histRandomHitDistribution->SetBinContent(1585,1);
histRandomHitDistribution->SetBinContent(1586,1);
histRandomHitDistribution->SetBinContent(1587,1);
histRandomHitDistribution->SetBinContent(1588,2);
histRandomHitDistribution->SetBinContent(1590,2);
histRandomHitDistribution->SetBinContent(1592,1);
histRandomHitDistribution->SetBinContent(1634,20);
histRandomHitDistribution->SetBinContent(1635,103);
histRandomHitDistribution->SetBinContent(1636,174);
histRandomHitDistribution->SetBinContent(1637,201);
histRandomHitDistribution->SetBinContent(1638,260);
histRandomHitDistribution->SetBinContent(1639,210);
histRandomHitDistribution->SetBinContent(1640,145);
histRandomHitDistribution->SetBinContent(1641,96);
histRandomHitDistribution->SetBinContent(1642,76);
histRandomHitDistribution->SetBinContent(1643,47);
histRandomHitDistribution->SetBinContent(1644,24);
histRandomHitDistribution->SetBinContent(1645,8);
histRandomHitDistribution->SetBinContent(1646,4);
histRandomHitDistribution->SetBinContent(1647,3);
histRandomHitDistribution->SetBinContent(1648,1);
histRandomHitDistribution->SetBinContent(1649,3);
histRandomHitDistribution->SetBinContent(1650,3);
histRandomHitDistribution->SetBinContent(1652,2);
histRandomHitDistribution->SetBinContent(1653,5);
histRandomHitDistribution->SetBinContent(1654,3);
histRandomHitDistribution->SetBinContent(1655,4);
histRandomHitDistribution->SetBinContent(1656,5);
histRandomHitDistribution->SetBinContent(1657,5);
histRandomHitDistribution->SetBinContent(1658,8);
histRandomHitDistribution->SetBinContent(1659,8);
histRandomHitDistribution->SetBinContent(1660,7);
histRandomHitDistribution->SetBinContent(1661,9);
histRandomHitDistribution->SetBinContent(1662,6);
histRandomHitDistribution->SetBinContent(1663,7);
histRandomHitDistribution->SetBinContent(1664,7);
histRandomHitDistribution->SetBinContent(1665,9);
histRandomHitDistribution->SetBinContent(1666,7);
histRandomHitDistribution->SetBinContent(1667,3);
histRandomHitDistribution->SetBinContent(1668,11);
histRandomHitDistribution->SetBinContent(1669,6);
histRandomHitDistribution->SetBinContent(1670,7);
histRandomHitDistribution->SetBinContent(1671,5);
histRandomHitDistribution->SetBinContent(1672,4);
histRandomHitDistribution->SetBinContent(1673,2);
histRandomHitDistribution->SetBinContent(1674,2);
histRandomHitDistribution->SetBinContent(1675,2);
histRandomHitDistribution->SetBinContent(1676,1);
histRandomHitDistribution->SetBinContent(1677,1);
histRandomHitDistribution->SetBinContent(1678,5);
histRandomHitDistribution->SetBinContent(1679,2);
histRandomHitDistribution->SetBinContent(1680,1);
histRandomHitDistribution->SetBinContent(1682,1);
histRandomHitDistribution->SetBinContent(1683,2);
histRandomHitDistribution->SetBinContent(1684,2);
histRandomHitDistribution->SetBinContent(1685,5);
histRandomHitDistribution->SetBinContent(1686,1);
histRandomHitDistribution->SetBinContent(1687,1);
histRandomHitDistribution->SetBinContent(1689,1);
histRandomHitDistribution->SetBinContent(1690,5);
histRandomHitDistribution->SetBinContent(1691,1);
histRandomHitDistribution->SetBinContent(1694,1);
histRandomHitDistribution->SetBinContent(1736,14);
histRandomHitDistribution->SetBinContent(1737,90);
histRandomHitDistribution->SetBinContent(1738,147);
histRandomHitDistribution->SetBinContent(1739,226);
histRandomHitDistribution->SetBinContent(1740,292);
histRandomHitDistribution->SetBinContent(1741,280);
histRandomHitDistribution->SetBinContent(1742,163);
histRandomHitDistribution->SetBinContent(1743,95);
histRandomHitDistribution->SetBinContent(1744,61);
histRandomHitDistribution->SetBinContent(1745,36);
histRandomHitDistribution->SetBinContent(1746,23);
histRandomHitDistribution->SetBinContent(1747,8);
histRandomHitDistribution->SetBinContent(1748,8);
histRandomHitDistribution->SetBinContent(1749,3);
histRandomHitDistribution->SetBinContent(1750,2);
histRandomHitDistribution->SetBinContent(1751,2);
histRandomHitDistribution->SetBinContent(1752,4);
histRandomHitDistribution->SetBinContent(1753,1);
histRandomHitDistribution->SetBinContent(1755,3);
histRandomHitDistribution->SetBinContent(1756,4);
histRandomHitDistribution->SetBinContent(1757,3);
histRandomHitDistribution->SetBinContent(1758,8);
histRandomHitDistribution->SetBinContent(1759,6);
histRandomHitDistribution->SetBinContent(1760,12);
histRandomHitDistribution->SetBinContent(1761,7);
histRandomHitDistribution->SetBinContent(1762,6);
histRandomHitDistribution->SetBinContent(1763,11);
histRandomHitDistribution->SetBinContent(1764,8);
histRandomHitDistribution->SetBinContent(1765,6);
histRandomHitDistribution->SetBinContent(1766,2);
histRandomHitDistribution->SetBinContent(1767,18);
histRandomHitDistribution->SetBinContent(1768,13);
histRandomHitDistribution->SetBinContent(1769,12);
histRandomHitDistribution->SetBinContent(1770,4);
histRandomHitDistribution->SetBinContent(1771,10);
histRandomHitDistribution->SetBinContent(1772,7);
histRandomHitDistribution->SetBinContent(1773,6);
histRandomHitDistribution->SetBinContent(1774,4);
histRandomHitDistribution->SetBinContent(1775,2);
histRandomHitDistribution->SetBinContent(1776,5);
histRandomHitDistribution->SetBinContent(1777,8);
histRandomHitDistribution->SetBinContent(1778,5);
histRandomHitDistribution->SetBinContent(1779,4);
histRandomHitDistribution->SetBinContent(1780,4);
histRandomHitDistribution->SetBinContent(1781,4);
histRandomHitDistribution->SetBinContent(1782,2);
histRandomHitDistribution->SetBinContent(1783,7);
histRandomHitDistribution->SetBinContent(1784,2);
histRandomHitDistribution->SetBinContent(1785,3);
histRandomHitDistribution->SetBinContent(1786,1);
histRandomHitDistribution->SetBinContent(1787,3);
histRandomHitDistribution->SetBinContent(1788,2);
histRandomHitDistribution->SetBinContent(1789,1);
histRandomHitDistribution->SetBinContent(1793,1);
histRandomHitDistribution->SetBinContent(1794,1);
histRandomHitDistribution->SetBinContent(1795,1);
histRandomHitDistribution->SetBinContent(1796,1);
histRandomHitDistribution->SetBinContent(1801,2);
histRandomHitDistribution->SetBinContent(1802,1);
histRandomHitDistribution->SetBinContent(1838,17);
histRandomHitDistribution->SetBinContent(1839,110);
histRandomHitDistribution->SetBinContent(1840,196);
histRandomHitDistribution->SetBinContent(1841,240);
histRandomHitDistribution->SetBinContent(1842,297);
histRandomHitDistribution->SetBinContent(1843,251);
histRandomHitDistribution->SetBinContent(1844,211);
histRandomHitDistribution->SetBinContent(1845,98);
histRandomHitDistribution->SetBinContent(1846,63);
histRandomHitDistribution->SetBinContent(1847,42);
histRandomHitDistribution->SetBinContent(1848,31);
histRandomHitDistribution->SetBinContent(1849,21);
histRandomHitDistribution->SetBinContent(1850,8);
histRandomHitDistribution->SetBinContent(1851,1);
histRandomHitDistribution->SetBinContent(1852,1);
histRandomHitDistribution->SetBinContent(1853,2);
histRandomHitDistribution->SetBinContent(1854,1);
histRandomHitDistribution->SetBinContent(1855,2);
histRandomHitDistribution->SetBinContent(1857,9);
histRandomHitDistribution->SetBinContent(1858,4);
histRandomHitDistribution->SetBinContent(1859,7);
histRandomHitDistribution->SetBinContent(1860,8);
histRandomHitDistribution->SetBinContent(1861,13);
histRandomHitDistribution->SetBinContent(1862,10);
histRandomHitDistribution->SetBinContent(1863,19);
histRandomHitDistribution->SetBinContent(1864,10);
histRandomHitDistribution->SetBinContent(1865,7);
histRandomHitDistribution->SetBinContent(1866,9);
histRandomHitDistribution->SetBinContent(1867,9);
histRandomHitDistribution->SetBinContent(1868,10);
histRandomHitDistribution->SetBinContent(1869,12);
histRandomHitDistribution->SetBinContent(1870,10);
histRandomHitDistribution->SetBinContent(1871,13);
histRandomHitDistribution->SetBinContent(1872,7);
histRandomHitDistribution->SetBinContent(1873,10);
histRandomHitDistribution->SetBinContent(1874,15);
histRandomHitDistribution->SetBinContent(1875,8);
histRandomHitDistribution->SetBinContent(1876,11);
histRandomHitDistribution->SetBinContent(1877,6);
histRandomHitDistribution->SetBinContent(1878,9);
histRandomHitDistribution->SetBinContent(1879,8);
histRandomHitDistribution->SetBinContent(1880,4);
histRandomHitDistribution->SetBinContent(1881,5);
histRandomHitDistribution->SetBinContent(1882,5);
histRandomHitDistribution->SetBinContent(1883,8);
histRandomHitDistribution->SetBinContent(1884,4);
histRandomHitDistribution->SetBinContent(1885,4);
histRandomHitDistribution->SetBinContent(1886,3);
histRandomHitDistribution->SetBinContent(1887,2);
histRandomHitDistribution->SetBinContent(1888,4);
histRandomHitDistribution->SetBinContent(1889,3);
histRandomHitDistribution->SetBinContent(1890,3);
histRandomHitDistribution->SetBinContent(1891,2);
histRandomHitDistribution->SetBinContent(1892,4);
histRandomHitDistribution->SetBinContent(1894,2);
histRandomHitDistribution->SetBinContent(1895,3);
histRandomHitDistribution->SetBinContent(1896,1);
histRandomHitDistribution->SetBinContent(1898,1);
histRandomHitDistribution->SetBinContent(1899,1);
histRandomHitDistribution->SetBinContent(1904,1);
histRandomHitDistribution->SetBinContent(1940,16);
histRandomHitDistribution->SetBinContent(1941,99);
histRandomHitDistribution->SetBinContent(1942,191);
histRandomHitDistribution->SetBinContent(1943,216);
histRandomHitDistribution->SetBinContent(1944,267);
histRandomHitDistribution->SetBinContent(1945,254);
histRandomHitDistribution->SetBinContent(1946,185);
histRandomHitDistribution->SetBinContent(1947,141);
histRandomHitDistribution->SetBinContent(1948,75);
histRandomHitDistribution->SetBinContent(1949,32);
histRandomHitDistribution->SetBinContent(1950,28);
histRandomHitDistribution->SetBinContent(1951,17);
histRandomHitDistribution->SetBinContent(1952,6);
histRandomHitDistribution->SetBinContent(1953,1);
histRandomHitDistribution->SetBinContent(1955,1);
histRandomHitDistribution->SetBinContent(1956,2);
histRandomHitDistribution->SetBinContent(1957,1);
histRandomHitDistribution->SetBinContent(1958,2);
histRandomHitDistribution->SetBinContent(1959,3);
histRandomHitDistribution->SetBinContent(1960,6);
histRandomHitDistribution->SetBinContent(1961,5);
histRandomHitDistribution->SetBinContent(1962,8);
histRandomHitDistribution->SetBinContent(1963,12);
histRandomHitDistribution->SetBinContent(1964,8);
histRandomHitDistribution->SetBinContent(1965,14);
histRandomHitDistribution->SetBinContent(1966,9);
histRandomHitDistribution->SetBinContent(1967,12);
histRandomHitDistribution->SetBinContent(1968,14);
histRandomHitDistribution->SetBinContent(1969,12);
histRandomHitDistribution->SetBinContent(1970,11);
histRandomHitDistribution->SetBinContent(1971,15);
histRandomHitDistribution->SetBinContent(1972,18);
histRandomHitDistribution->SetBinContent(1973,14);
histRandomHitDistribution->SetBinContent(1974,12);
histRandomHitDistribution->SetBinContent(1975,9);
histRandomHitDistribution->SetBinContent(1976,16);
histRandomHitDistribution->SetBinContent(1977,10);
histRandomHitDistribution->SetBinContent(1978,5);
histRandomHitDistribution->SetBinContent(1979,11);
histRandomHitDistribution->SetBinContent(1980,9);
histRandomHitDistribution->SetBinContent(1981,15);
histRandomHitDistribution->SetBinContent(1982,12);
histRandomHitDistribution->SetBinContent(1983,9);
histRandomHitDistribution->SetBinContent(1984,15);
histRandomHitDistribution->SetBinContent(1985,9);
histRandomHitDistribution->SetBinContent(1986,4);
histRandomHitDistribution->SetBinContent(1987,6);
histRandomHitDistribution->SetBinContent(1988,10);
histRandomHitDistribution->SetBinContent(1989,7);
histRandomHitDistribution->SetBinContent(1990,2);
histRandomHitDistribution->SetBinContent(1991,4);
histRandomHitDistribution->SetBinContent(1992,3);
histRandomHitDistribution->SetBinContent(1993,3);
histRandomHitDistribution->SetBinContent(1994,4);
histRandomHitDistribution->SetBinContent(1996,6);
histRandomHitDistribution->SetBinContent(1997,2);
histRandomHitDistribution->SetBinContent(1998,1);
histRandomHitDistribution->SetBinContent(1999,1);
histRandomHitDistribution->SetBinContent(2000,2);
histRandomHitDistribution->SetBinContent(2001,2);
histRandomHitDistribution->SetBinContent(2004,1);
histRandomHitDistribution->SetBinContent(2006,1);
histRandomHitDistribution->SetBinContent(2013,1);
histRandomHitDistribution->SetBinContent(2042,14);
histRandomHitDistribution->SetBinContent(2043,90);
histRandomHitDistribution->SetBinContent(2044,161);
histRandomHitDistribution->SetBinContent(2045,228);
histRandomHitDistribution->SetBinContent(2046,255);
histRandomHitDistribution->SetBinContent(2047,213);
histRandomHitDistribution->SetBinContent(2048,201);
histRandomHitDistribution->SetBinContent(2049,110);
histRandomHitDistribution->SetBinContent(2050,91);
histRandomHitDistribution->SetBinContent(2051,37);
histRandomHitDistribution->SetBinContent(2052,24);
histRandomHitDistribution->SetBinContent(2053,26);
histRandomHitDistribution->SetBinContent(2054,10);
histRandomHitDistribution->SetBinContent(2055,5);
histRandomHitDistribution->SetBinContent(2056,1);
histRandomHitDistribution->SetBinContent(2057,1);
histRandomHitDistribution->SetBinContent(2058,3);
histRandomHitDistribution->SetBinContent(2060,4);
histRandomHitDistribution->SetBinContent(2061,5);
histRandomHitDistribution->SetBinContent(2062,5);
histRandomHitDistribution->SetBinContent(2063,6);
histRandomHitDistribution->SetBinContent(2064,5);
histRandomHitDistribution->SetBinContent(2065,9);
histRandomHitDistribution->SetBinContent(2066,14);
histRandomHitDistribution->SetBinContent(2067,12);
histRandomHitDistribution->SetBinContent(2068,9);
histRandomHitDistribution->SetBinContent(2069,16);
histRandomHitDistribution->SetBinContent(2070,20);
histRandomHitDistribution->SetBinContent(2071,13);
histRandomHitDistribution->SetBinContent(2072,24);
histRandomHitDistribution->SetBinContent(2073,19);
histRandomHitDistribution->SetBinContent(2074,15);
histRandomHitDistribution->SetBinContent(2075,17);
histRandomHitDistribution->SetBinContent(2076,12);
histRandomHitDistribution->SetBinContent(2077,21);
histRandomHitDistribution->SetBinContent(2078,16);
histRandomHitDistribution->SetBinContent(2079,21);
histRandomHitDistribution->SetBinContent(2080,23);
histRandomHitDistribution->SetBinContent(2081,20);
histRandomHitDistribution->SetBinContent(2082,13);
histRandomHitDistribution->SetBinContent(2083,14);
histRandomHitDistribution->SetBinContent(2084,13);
histRandomHitDistribution->SetBinContent(2085,13);
histRandomHitDistribution->SetBinContent(2086,13);
histRandomHitDistribution->SetBinContent(2087,10);
histRandomHitDistribution->SetBinContent(2088,13);
histRandomHitDistribution->SetBinContent(2089,10);
histRandomHitDistribution->SetBinContent(2090,7);
histRandomHitDistribution->SetBinContent(2091,9);
histRandomHitDistribution->SetBinContent(2092,10);
histRandomHitDistribution->SetBinContent(2093,5);
histRandomHitDistribution->SetBinContent(2094,5);
histRandomHitDistribution->SetBinContent(2095,8);
histRandomHitDistribution->SetBinContent(2096,3);
histRandomHitDistribution->SetBinContent(2097,4);
histRandomHitDistribution->SetBinContent(2098,3);
histRandomHitDistribution->SetBinContent(2099,4);
histRandomHitDistribution->SetBinContent(2101,1);
histRandomHitDistribution->SetBinContent(2102,1);
histRandomHitDistribution->SetBinContent(2103,2);
histRandomHitDistribution->SetBinContent(2104,1);
histRandomHitDistribution->SetBinContent(2105,3);
histRandomHitDistribution->SetBinContent(2107,1);
histRandomHitDistribution->SetBinContent(2114,1);
histRandomHitDistribution->SetBinContent(2144,14);
histRandomHitDistribution->SetBinContent(2145,97);
histRandomHitDistribution->SetBinContent(2146,153);
histRandomHitDistribution->SetBinContent(2147,262);
histRandomHitDistribution->SetBinContent(2148,232);
histRandomHitDistribution->SetBinContent(2149,243);
histRandomHitDistribution->SetBinContent(2150,227);
histRandomHitDistribution->SetBinContent(2151,175);
histRandomHitDistribution->SetBinContent(2152,142);
histRandomHitDistribution->SetBinContent(2153,62);
histRandomHitDistribution->SetBinContent(2154,24);
histRandomHitDistribution->SetBinContent(2155,17);
histRandomHitDistribution->SetBinContent(2156,11);
histRandomHitDistribution->SetBinContent(2157,3);
histRandomHitDistribution->SetBinContent(2158,8);
histRandomHitDistribution->SetBinContent(2159,3);
histRandomHitDistribution->SetBinContent(2160,1);
histRandomHitDistribution->SetBinContent(2163,8);
histRandomHitDistribution->SetBinContent(2164,6);
histRandomHitDistribution->SetBinContent(2165,14);
histRandomHitDistribution->SetBinContent(2166,6);
histRandomHitDistribution->SetBinContent(2167,12);
histRandomHitDistribution->SetBinContent(2168,15);
histRandomHitDistribution->SetBinContent(2169,25);
histRandomHitDistribution->SetBinContent(2170,11);
histRandomHitDistribution->SetBinContent(2171,20);
histRandomHitDistribution->SetBinContent(2172,23);
histRandomHitDistribution->SetBinContent(2173,16);
histRandomHitDistribution->SetBinContent(2174,20);
histRandomHitDistribution->SetBinContent(2175,19);
histRandomHitDistribution->SetBinContent(2176,33);
histRandomHitDistribution->SetBinContent(2177,20);
histRandomHitDistribution->SetBinContent(2178,24);
histRandomHitDistribution->SetBinContent(2179,16);
histRandomHitDistribution->SetBinContent(2180,32);
histRandomHitDistribution->SetBinContent(2181,17);
histRandomHitDistribution->SetBinContent(2182,31);
histRandomHitDistribution->SetBinContent(2183,14);
histRandomHitDistribution->SetBinContent(2184,24);
histRandomHitDistribution->SetBinContent(2185,21);
histRandomHitDistribution->SetBinContent(2186,17);
histRandomHitDistribution->SetBinContent(2187,19);
histRandomHitDistribution->SetBinContent(2188,15);
histRandomHitDistribution->SetBinContent(2189,16);
histRandomHitDistribution->SetBinContent(2190,14);
histRandomHitDistribution->SetBinContent(2191,15);
histRandomHitDistribution->SetBinContent(2192,9);
histRandomHitDistribution->SetBinContent(2193,8);
histRandomHitDistribution->SetBinContent(2194,9);
histRandomHitDistribution->SetBinContent(2195,11);
histRandomHitDistribution->SetBinContent(2196,4);
histRandomHitDistribution->SetBinContent(2197,9);
histRandomHitDistribution->SetBinContent(2198,5);
histRandomHitDistribution->SetBinContent(2199,3);
histRandomHitDistribution->SetBinContent(2200,3);
histRandomHitDistribution->SetBinContent(2201,2);
histRandomHitDistribution->SetBinContent(2202,3);
histRandomHitDistribution->SetBinContent(2203,1);
histRandomHitDistribution->SetBinContent(2204,1);
histRandomHitDistribution->SetBinContent(2205,2);
histRandomHitDistribution->SetBinContent(2206,1);
histRandomHitDistribution->SetBinContent(2207,3);
histRandomHitDistribution->SetBinContent(2208,1);
histRandomHitDistribution->SetBinContent(2210,1);
histRandomHitDistribution->SetBinContent(2217,2);
histRandomHitDistribution->SetBinContent(2246,17);
histRandomHitDistribution->SetBinContent(2247,77);
histRandomHitDistribution->SetBinContent(2248,206);
histRandomHitDistribution->SetBinContent(2249,246);
histRandomHitDistribution->SetBinContent(2250,238);
histRandomHitDistribution->SetBinContent(2251,255);
histRandomHitDistribution->SetBinContent(2252,194);
histRandomHitDistribution->SetBinContent(2253,204);
histRandomHitDistribution->SetBinContent(2254,123);
histRandomHitDistribution->SetBinContent(2255,58);
histRandomHitDistribution->SetBinContent(2256,40);
histRandomHitDistribution->SetBinContent(2257,25);
histRandomHitDistribution->SetBinContent(2258,9);
histRandomHitDistribution->SetBinContent(2259,6);
histRandomHitDistribution->SetBinContent(2260,2);
histRandomHitDistribution->SetBinContent(2261,1);
histRandomHitDistribution->SetBinContent(2263,2);
histRandomHitDistribution->SetBinContent(2264,2);
histRandomHitDistribution->SetBinContent(2265,9);
histRandomHitDistribution->SetBinContent(2266,9);
histRandomHitDistribution->SetBinContent(2267,9);
histRandomHitDistribution->SetBinContent(2268,9);
histRandomHitDistribution->SetBinContent(2269,16);
histRandomHitDistribution->SetBinContent(2270,24);
histRandomHitDistribution->SetBinContent(2271,12);
histRandomHitDistribution->SetBinContent(2272,31);
histRandomHitDistribution->SetBinContent(2273,33);
histRandomHitDistribution->SetBinContent(2274,20);
histRandomHitDistribution->SetBinContent(2275,32);
histRandomHitDistribution->SetBinContent(2276,27);
histRandomHitDistribution->SetBinContent(2277,25);
histRandomHitDistribution->SetBinContent(2278,18);
histRandomHitDistribution->SetBinContent(2279,35);
histRandomHitDistribution->SetBinContent(2280,32);
histRandomHitDistribution->SetBinContent(2281,34);
histRandomHitDistribution->SetBinContent(2282,32);
histRandomHitDistribution->SetBinContent(2283,26);
histRandomHitDistribution->SetBinContent(2284,25);
histRandomHitDistribution->SetBinContent(2285,43);
histRandomHitDistribution->SetBinContent(2286,35);
histRandomHitDistribution->SetBinContent(2287,30);
histRandomHitDistribution->SetBinContent(2288,22);
histRandomHitDistribution->SetBinContent(2289,24);
histRandomHitDistribution->SetBinContent(2290,15);
histRandomHitDistribution->SetBinContent(2291,22);
histRandomHitDistribution->SetBinContent(2292,18);
histRandomHitDistribution->SetBinContent(2293,12);
histRandomHitDistribution->SetBinContent(2294,9);
histRandomHitDistribution->SetBinContent(2295,11);
histRandomHitDistribution->SetBinContent(2296,17);
histRandomHitDistribution->SetBinContent(2297,16);
histRandomHitDistribution->SetBinContent(2298,7);
histRandomHitDistribution->SetBinContent(2299,9);
histRandomHitDistribution->SetBinContent(2300,11);
histRandomHitDistribution->SetBinContent(2301,9);
histRandomHitDistribution->SetBinContent(2302,5);
histRandomHitDistribution->SetBinContent(2303,6);
histRandomHitDistribution->SetBinContent(2304,4);
histRandomHitDistribution->SetBinContent(2305,4);
histRandomHitDistribution->SetBinContent(2306,2);
histRandomHitDistribution->SetBinContent(2307,5);
histRandomHitDistribution->SetBinContent(2308,1);
histRandomHitDistribution->SetBinContent(2310,3);
histRandomHitDistribution->SetBinContent(2314,1);
histRandomHitDistribution->SetBinContent(2317,1);
histRandomHitDistribution->SetBinContent(2348,14);
histRandomHitDistribution->SetBinContent(2349,73);
histRandomHitDistribution->SetBinContent(2350,194);
histRandomHitDistribution->SetBinContent(2351,259);
histRandomHitDistribution->SetBinContent(2352,261);
histRandomHitDistribution->SetBinContent(2353,240);
histRandomHitDistribution->SetBinContent(2354,209);
histRandomHitDistribution->SetBinContent(2355,162);
histRandomHitDistribution->SetBinContent(2356,118);
histRandomHitDistribution->SetBinContent(2357,64);
histRandomHitDistribution->SetBinContent(2358,47);
histRandomHitDistribution->SetBinContent(2359,26);
histRandomHitDistribution->SetBinContent(2360,12);
histRandomHitDistribution->SetBinContent(2361,6);
histRandomHitDistribution->SetBinContent(2362,7);
histRandomHitDistribution->SetBinContent(2365,4);
histRandomHitDistribution->SetBinContent(2366,5);
histRandomHitDistribution->SetBinContent(2367,8);
histRandomHitDistribution->SetBinContent(2368,12);
histRandomHitDistribution->SetBinContent(2369,12);
histRandomHitDistribution->SetBinContent(2370,16);
histRandomHitDistribution->SetBinContent(2371,19);
histRandomHitDistribution->SetBinContent(2372,15);
histRandomHitDistribution->SetBinContent(2373,18);
histRandomHitDistribution->SetBinContent(2374,22);
histRandomHitDistribution->SetBinContent(2375,17);
histRandomHitDistribution->SetBinContent(2376,26);
histRandomHitDistribution->SetBinContent(2377,33);
histRandomHitDistribution->SetBinContent(2378,29);
histRandomHitDistribution->SetBinContent(2379,23);
histRandomHitDistribution->SetBinContent(2380,27);
histRandomHitDistribution->SetBinContent(2381,24);
histRandomHitDistribution->SetBinContent(2382,36);
histRandomHitDistribution->SetBinContent(2383,41);
histRandomHitDistribution->SetBinContent(2384,30);
histRandomHitDistribution->SetBinContent(2385,28);
histRandomHitDistribution->SetBinContent(2386,33);
histRandomHitDistribution->SetBinContent(2387,21);
histRandomHitDistribution->SetBinContent(2388,36);
histRandomHitDistribution->SetBinContent(2389,34);
histRandomHitDistribution->SetBinContent(2390,27);
histRandomHitDistribution->SetBinContent(2391,25);
histRandomHitDistribution->SetBinContent(2392,28);
histRandomHitDistribution->SetBinContent(2393,22);
histRandomHitDistribution->SetBinContent(2394,20);
histRandomHitDistribution->SetBinContent(2395,18);
histRandomHitDistribution->SetBinContent(2396,18);
histRandomHitDistribution->SetBinContent(2397,20);
histRandomHitDistribution->SetBinContent(2398,20);
histRandomHitDistribution->SetBinContent(2399,14);
histRandomHitDistribution->SetBinContent(2400,11);
histRandomHitDistribution->SetBinContent(2401,14);
histRandomHitDistribution->SetBinContent(2402,7);
histRandomHitDistribution->SetBinContent(2403,10);
histRandomHitDistribution->SetBinContent(2404,5);
histRandomHitDistribution->SetBinContent(2405,7);
histRandomHitDistribution->SetBinContent(2406,6);
histRandomHitDistribution->SetBinContent(2407,4);
histRandomHitDistribution->SetBinContent(2408,4);
histRandomHitDistribution->SetBinContent(2409,6);
histRandomHitDistribution->SetBinContent(2410,3);
histRandomHitDistribution->SetBinContent(2411,3);
histRandomHitDistribution->SetBinContent(2412,3);
histRandomHitDistribution->SetBinContent(2413,2);
histRandomHitDistribution->SetBinContent(2414,3);
histRandomHitDistribution->SetBinContent(2415,2);
histRandomHitDistribution->SetBinContent(2416,2);
histRandomHitDistribution->SetBinContent(2417,2);
histRandomHitDistribution->SetBinContent(2418,1);
histRandomHitDistribution->SetBinContent(2419,1);
histRandomHitDistribution->SetBinContent(2421,1);
histRandomHitDistribution->SetBinContent(2428,1);
histRandomHitDistribution->SetBinContent(2450,18);
histRandomHitDistribution->SetBinContent(2451,95);
histRandomHitDistribution->SetBinContent(2452,225);
histRandomHitDistribution->SetBinContent(2453,242);
histRandomHitDistribution->SetBinContent(2454,284);
histRandomHitDistribution->SetBinContent(2455,232);
histRandomHitDistribution->SetBinContent(2456,206);
histRandomHitDistribution->SetBinContent(2457,156);
histRandomHitDistribution->SetBinContent(2458,135);
histRandomHitDistribution->SetBinContent(2459,64);
histRandomHitDistribution->SetBinContent(2460,53);
histRandomHitDistribution->SetBinContent(2461,29);
histRandomHitDistribution->SetBinContent(2462,15);
histRandomHitDistribution->SetBinContent(2463,13);
histRandomHitDistribution->SetBinContent(2464,4);
histRandomHitDistribution->SetBinContent(2465,1);
histRandomHitDistribution->SetBinContent(2466,2);
histRandomHitDistribution->SetBinContent(2467,3);
histRandomHitDistribution->SetBinContent(2468,5);
histRandomHitDistribution->SetBinContent(2469,6);
histRandomHitDistribution->SetBinContent(2470,9);
histRandomHitDistribution->SetBinContent(2471,4);
histRandomHitDistribution->SetBinContent(2472,14);
histRandomHitDistribution->SetBinContent(2473,17);
histRandomHitDistribution->SetBinContent(2474,17);
histRandomHitDistribution->SetBinContent(2475,20);
histRandomHitDistribution->SetBinContent(2476,28);
histRandomHitDistribution->SetBinContent(2477,22);
histRandomHitDistribution->SetBinContent(2478,31);
histRandomHitDistribution->SetBinContent(2479,38);
histRandomHitDistribution->SetBinContent(2480,24);
histRandomHitDistribution->SetBinContent(2481,24);
histRandomHitDistribution->SetBinContent(2482,23);
histRandomHitDistribution->SetBinContent(2483,34);
histRandomHitDistribution->SetBinContent(2484,38);
histRandomHitDistribution->SetBinContent(2485,38);
histRandomHitDistribution->SetBinContent(2486,37);
histRandomHitDistribution->SetBinContent(2487,32);
histRandomHitDistribution->SetBinContent(2488,39);
histRandomHitDistribution->SetBinContent(2489,35);
histRandomHitDistribution->SetBinContent(2490,36);
histRandomHitDistribution->SetBinContent(2491,34);
histRandomHitDistribution->SetBinContent(2492,35);
histRandomHitDistribution->SetBinContent(2493,25);
histRandomHitDistribution->SetBinContent(2494,26);
histRandomHitDistribution->SetBinContent(2495,27);
histRandomHitDistribution->SetBinContent(2496,25);
histRandomHitDistribution->SetBinContent(2497,20);
histRandomHitDistribution->SetBinContent(2498,22);
histRandomHitDistribution->SetBinContent(2499,27);
histRandomHitDistribution->SetBinContent(2500,23);
histRandomHitDistribution->SetBinContent(2501,17);
histRandomHitDistribution->SetBinContent(2502,14);
histRandomHitDistribution->SetBinContent(2503,15);
histRandomHitDistribution->SetBinContent(2504,17);
histRandomHitDistribution->SetBinContent(2505,12);
histRandomHitDistribution->SetBinContent(2506,11);
histRandomHitDistribution->SetBinContent(2507,11);
histRandomHitDistribution->SetBinContent(2508,10);
histRandomHitDistribution->SetBinContent(2509,8);
histRandomHitDistribution->SetBinContent(2510,8);
histRandomHitDistribution->SetBinContent(2511,6);
histRandomHitDistribution->SetBinContent(2512,8);
histRandomHitDistribution->SetBinContent(2513,3);
histRandomHitDistribution->SetBinContent(2514,7);
histRandomHitDistribution->SetBinContent(2515,4);
histRandomHitDistribution->SetBinContent(2516,2);
histRandomHitDistribution->SetBinContent(2517,2);
histRandomHitDistribution->SetBinContent(2518,1);
histRandomHitDistribution->SetBinContent(2519,2);
histRandomHitDistribution->SetBinContent(2520,1);
histRandomHitDistribution->SetBinContent(2521,1);
histRandomHitDistribution->SetBinContent(2522,1);
histRandomHitDistribution->SetBinContent(2523,2);
histRandomHitDistribution->SetBinContent(2524,1);
histRandomHitDistribution->SetBinContent(2525,1);
histRandomHitDistribution->SetBinContent(2552,7);
histRandomHitDistribution->SetBinContent(2553,83);
histRandomHitDistribution->SetBinContent(2554,201);
histRandomHitDistribution->SetBinContent(2555,288);
histRandomHitDistribution->SetBinContent(2556,301);
histRandomHitDistribution->SetBinContent(2557,278);
histRandomHitDistribution->SetBinContent(2558,228);
histRandomHitDistribution->SetBinContent(2559,185);
histRandomHitDistribution->SetBinContent(2560,144);
histRandomHitDistribution->SetBinContent(2561,79);
histRandomHitDistribution->SetBinContent(2562,64);
histRandomHitDistribution->SetBinContent(2563,34);
histRandomHitDistribution->SetBinContent(2564,22);
histRandomHitDistribution->SetBinContent(2565,7);
histRandomHitDistribution->SetBinContent(2566,6);
histRandomHitDistribution->SetBinContent(2567,8);
histRandomHitDistribution->SetBinContent(2568,4);
histRandomHitDistribution->SetBinContent(2569,4);
histRandomHitDistribution->SetBinContent(2570,7);
histRandomHitDistribution->SetBinContent(2571,5);
histRandomHitDistribution->SetBinContent(2572,8);
histRandomHitDistribution->SetBinContent(2573,21);
histRandomHitDistribution->SetBinContent(2574,17);
histRandomHitDistribution->SetBinContent(2575,16);
histRandomHitDistribution->SetBinContent(2576,19);
histRandomHitDistribution->SetBinContent(2577,20);
histRandomHitDistribution->SetBinContent(2578,26);
histRandomHitDistribution->SetBinContent(2579,28);
histRandomHitDistribution->SetBinContent(2580,29);
histRandomHitDistribution->SetBinContent(2581,16);
histRandomHitDistribution->SetBinContent(2582,35);
histRandomHitDistribution->SetBinContent(2583,29);
histRandomHitDistribution->SetBinContent(2584,25);
histRandomHitDistribution->SetBinContent(2585,33);
histRandomHitDistribution->SetBinContent(2586,37);
histRandomHitDistribution->SetBinContent(2587,36);
histRandomHitDistribution->SetBinContent(2588,38);
histRandomHitDistribution->SetBinContent(2589,51);
histRandomHitDistribution->SetBinContent(2590,35);
histRandomHitDistribution->SetBinContent(2591,34);
histRandomHitDistribution->SetBinContent(2592,28);
histRandomHitDistribution->SetBinContent(2593,39);
histRandomHitDistribution->SetBinContent(2594,31);
histRandomHitDistribution->SetBinContent(2595,31);
histRandomHitDistribution->SetBinContent(2596,38);
histRandomHitDistribution->SetBinContent(2597,31);
histRandomHitDistribution->SetBinContent(2598,30);
histRandomHitDistribution->SetBinContent(2599,27);
histRandomHitDistribution->SetBinContent(2600,19);
histRandomHitDistribution->SetBinContent(2601,16);
histRandomHitDistribution->SetBinContent(2602,18);
histRandomHitDistribution->SetBinContent(2603,22);
histRandomHitDistribution->SetBinContent(2604,17);
histRandomHitDistribution->SetBinContent(2605,12);
histRandomHitDistribution->SetBinContent(2606,13);
histRandomHitDistribution->SetBinContent(2607,12);
histRandomHitDistribution->SetBinContent(2608,9);
histRandomHitDistribution->SetBinContent(2609,12);
histRandomHitDistribution->SetBinContent(2610,11);
histRandomHitDistribution->SetBinContent(2611,9);
histRandomHitDistribution->SetBinContent(2612,14);
histRandomHitDistribution->SetBinContent(2613,7);
histRandomHitDistribution->SetBinContent(2614,9);
histRandomHitDistribution->SetBinContent(2615,4);
histRandomHitDistribution->SetBinContent(2616,6);
histRandomHitDistribution->SetBinContent(2617,1);
histRandomHitDistribution->SetBinContent(2618,6);
histRandomHitDistribution->SetBinContent(2619,6);
histRandomHitDistribution->SetBinContent(2620,2);
histRandomHitDistribution->SetBinContent(2621,2);
histRandomHitDistribution->SetBinContent(2623,4);
histRandomHitDistribution->SetBinContent(2624,2);
histRandomHitDistribution->SetBinContent(2625,3);
histRandomHitDistribution->SetBinContent(2628,2);
histRandomHitDistribution->SetBinContent(2630,1);
histRandomHitDistribution->SetBinContent(2654,16);
histRandomHitDistribution->SetBinContent(2655,89);
histRandomHitDistribution->SetBinContent(2656,199);
histRandomHitDistribution->SetBinContent(2657,260);
histRandomHitDistribution->SetBinContent(2658,338);
histRandomHitDistribution->SetBinContent(2659,286);
histRandomHitDistribution->SetBinContent(2660,298);
histRandomHitDistribution->SetBinContent(2661,211);
histRandomHitDistribution->SetBinContent(2662,143);
histRandomHitDistribution->SetBinContent(2663,115);
histRandomHitDistribution->SetBinContent(2664,61);
histRandomHitDistribution->SetBinContent(2665,38);
histRandomHitDistribution->SetBinContent(2666,23);
histRandomHitDistribution->SetBinContent(2667,8);
histRandomHitDistribution->SetBinContent(2668,4);
histRandomHitDistribution->SetBinContent(2669,4);
histRandomHitDistribution->SetBinContent(2670,3);
histRandomHitDistribution->SetBinContent(2671,3);
histRandomHitDistribution->SetBinContent(2672,8);
histRandomHitDistribution->SetBinContent(2673,5);
histRandomHitDistribution->SetBinContent(2674,12);
histRandomHitDistribution->SetBinContent(2675,15);
histRandomHitDistribution->SetBinContent(2676,17);
histRandomHitDistribution->SetBinContent(2677,25);
histRandomHitDistribution->SetBinContent(2678,20);
histRandomHitDistribution->SetBinContent(2679,20);
histRandomHitDistribution->SetBinContent(2680,25);
histRandomHitDistribution->SetBinContent(2681,30);
histRandomHitDistribution->SetBinContent(2682,24);
histRandomHitDistribution->SetBinContent(2683,33);
histRandomHitDistribution->SetBinContent(2684,37);
histRandomHitDistribution->SetBinContent(2685,45);
histRandomHitDistribution->SetBinContent(2686,53);
histRandomHitDistribution->SetBinContent(2687,37);
histRandomHitDistribution->SetBinContent(2688,35);
histRandomHitDistribution->SetBinContent(2689,53);
histRandomHitDistribution->SetBinContent(2690,42);
histRandomHitDistribution->SetBinContent(2691,51);
histRandomHitDistribution->SetBinContent(2692,37);
histRandomHitDistribution->SetBinContent(2693,55);
histRandomHitDistribution->SetBinContent(2694,45);
histRandomHitDistribution->SetBinContent(2695,39);
histRandomHitDistribution->SetBinContent(2696,43);
histRandomHitDistribution->SetBinContent(2697,40);
histRandomHitDistribution->SetBinContent(2698,44);
histRandomHitDistribution->SetBinContent(2699,29);
histRandomHitDistribution->SetBinContent(2700,38);
histRandomHitDistribution->SetBinContent(2701,27);
histRandomHitDistribution->SetBinContent(2702,16);
histRandomHitDistribution->SetBinContent(2703,27);
histRandomHitDistribution->SetBinContent(2704,24);
histRandomHitDistribution->SetBinContent(2705,20);
histRandomHitDistribution->SetBinContent(2706,18);
histRandomHitDistribution->SetBinContent(2707,18);
histRandomHitDistribution->SetBinContent(2708,20);
histRandomHitDistribution->SetBinContent(2709,14);
histRandomHitDistribution->SetBinContent(2710,10);
histRandomHitDistribution->SetBinContent(2711,4);
histRandomHitDistribution->SetBinContent(2712,14);
histRandomHitDistribution->SetBinContent(2713,11);
histRandomHitDistribution->SetBinContent(2714,14);
histRandomHitDistribution->SetBinContent(2715,11);
histRandomHitDistribution->SetBinContent(2716,9);
histRandomHitDistribution->SetBinContent(2717,5);
histRandomHitDistribution->SetBinContent(2718,12);
histRandomHitDistribution->SetBinContent(2719,2);
histRandomHitDistribution->SetBinContent(2720,8);
histRandomHitDistribution->SetBinContent(2721,1);
histRandomHitDistribution->SetBinContent(2722,3);
histRandomHitDistribution->SetBinContent(2724,2);
histRandomHitDistribution->SetBinContent(2725,1);
histRandomHitDistribution->SetBinContent(2726,1);
histRandomHitDistribution->SetBinContent(2727,2);
histRandomHitDistribution->SetBinContent(2728,1);
histRandomHitDistribution->SetBinContent(2729,2);
histRandomHitDistribution->SetBinContent(2730,2);
histRandomHitDistribution->SetBinContent(2732,2);
histRandomHitDistribution->SetBinContent(2735,1);
histRandomHitDistribution->SetBinContent(2738,1);
histRandomHitDistribution->SetBinContent(2756,17);
histRandomHitDistribution->SetBinContent(2757,95);
histRandomHitDistribution->SetBinContent(2758,216);
histRandomHitDistribution->SetBinContent(2759,284);
histRandomHitDistribution->SetBinContent(2760,316);
histRandomHitDistribution->SetBinContent(2761,312);
histRandomHitDistribution->SetBinContent(2762,287);
histRandomHitDistribution->SetBinContent(2763,224);
histRandomHitDistribution->SetBinContent(2764,148);
histRandomHitDistribution->SetBinContent(2765,92);
histRandomHitDistribution->SetBinContent(2766,55);
histRandomHitDistribution->SetBinContent(2767,36);
histRandomHitDistribution->SetBinContent(2768,36);
histRandomHitDistribution->SetBinContent(2769,20);
histRandomHitDistribution->SetBinContent(2770,6);
histRandomHitDistribution->SetBinContent(2772,5);
histRandomHitDistribution->SetBinContent(2773,2);
histRandomHitDistribution->SetBinContent(2774,8);
histRandomHitDistribution->SetBinContent(2775,8);
histRandomHitDistribution->SetBinContent(2776,8);
histRandomHitDistribution->SetBinContent(2777,16);
histRandomHitDistribution->SetBinContent(2778,15);
histRandomHitDistribution->SetBinContent(2779,26);
histRandomHitDistribution->SetBinContent(2780,17);
histRandomHitDistribution->SetBinContent(2781,29);
histRandomHitDistribution->SetBinContent(2782,29);
histRandomHitDistribution->SetBinContent(2783,24);
histRandomHitDistribution->SetBinContent(2784,34);
histRandomHitDistribution->SetBinContent(2785,42);
histRandomHitDistribution->SetBinContent(2786,32);
histRandomHitDistribution->SetBinContent(2787,31);
histRandomHitDistribution->SetBinContent(2788,53);
histRandomHitDistribution->SetBinContent(2789,58);
histRandomHitDistribution->SetBinContent(2790,46);
histRandomHitDistribution->SetBinContent(2791,37);
histRandomHitDistribution->SetBinContent(2792,41);
histRandomHitDistribution->SetBinContent(2793,42);
histRandomHitDistribution->SetBinContent(2794,55);
histRandomHitDistribution->SetBinContent(2795,39);
histRandomHitDistribution->SetBinContent(2796,55);
histRandomHitDistribution->SetBinContent(2797,45);
histRandomHitDistribution->SetBinContent(2798,60);
histRandomHitDistribution->SetBinContent(2799,39);
histRandomHitDistribution->SetBinContent(2800,52);
histRandomHitDistribution->SetBinContent(2801,51);
histRandomHitDistribution->SetBinContent(2802,36);
histRandomHitDistribution->SetBinContent(2803,31);
histRandomHitDistribution->SetBinContent(2804,29);
histRandomHitDistribution->SetBinContent(2805,36);
histRandomHitDistribution->SetBinContent(2806,26);
histRandomHitDistribution->SetBinContent(2807,20);
histRandomHitDistribution->SetBinContent(2808,24);
histRandomHitDistribution->SetBinContent(2809,31);
histRandomHitDistribution->SetBinContent(2810,20);
histRandomHitDistribution->SetBinContent(2811,20);
histRandomHitDistribution->SetBinContent(2812,13);
histRandomHitDistribution->SetBinContent(2813,12);
histRandomHitDistribution->SetBinContent(2814,12);
histRandomHitDistribution->SetBinContent(2815,8);
histRandomHitDistribution->SetBinContent(2816,9);
histRandomHitDistribution->SetBinContent(2817,14);
histRandomHitDistribution->SetBinContent(2818,8);
histRandomHitDistribution->SetBinContent(2819,5);
histRandomHitDistribution->SetBinContent(2820,6);
histRandomHitDistribution->SetBinContent(2821,4);
histRandomHitDistribution->SetBinContent(2822,3);
histRandomHitDistribution->SetBinContent(2823,4);
histRandomHitDistribution->SetBinContent(2824,5);
histRandomHitDistribution->SetBinContent(2825,1);
histRandomHitDistribution->SetBinContent(2826,2);
histRandomHitDistribution->SetBinContent(2828,1);
histRandomHitDistribution->SetBinContent(2829,2);
histRandomHitDistribution->SetBinContent(2830,1);
histRandomHitDistribution->SetBinContent(2831,2);
histRandomHitDistribution->SetBinContent(2833,3);
histRandomHitDistribution->SetBinContent(2846,1);
histRandomHitDistribution->SetBinContent(2858,9);
histRandomHitDistribution->SetBinContent(2859,67);
histRandomHitDistribution->SetBinContent(2860,211);
histRandomHitDistribution->SetBinContent(2861,295);
histRandomHitDistribution->SetBinContent(2862,329);
histRandomHitDistribution->SetBinContent(2863,281);
histRandomHitDistribution->SetBinContent(2864,299);
histRandomHitDistribution->SetBinContent(2865,196);
histRandomHitDistribution->SetBinContent(2866,163);
histRandomHitDistribution->SetBinContent(2867,88);
histRandomHitDistribution->SetBinContent(2868,51);
histRandomHitDistribution->SetBinContent(2869,26);
histRandomHitDistribution->SetBinContent(2870,21);
histRandomHitDistribution->SetBinContent(2871,9);
histRandomHitDistribution->SetBinContent(2872,9);
histRandomHitDistribution->SetBinContent(2873,8);
histRandomHitDistribution->SetBinContent(2874,5);
histRandomHitDistribution->SetBinContent(2875,2);
histRandomHitDistribution->SetBinContent(2876,4);
histRandomHitDistribution->SetBinContent(2877,7);
histRandomHitDistribution->SetBinContent(2878,17);
histRandomHitDistribution->SetBinContent(2879,17);
histRandomHitDistribution->SetBinContent(2880,14);
histRandomHitDistribution->SetBinContent(2881,14);
histRandomHitDistribution->SetBinContent(2882,30);
histRandomHitDistribution->SetBinContent(2883,39);
histRandomHitDistribution->SetBinContent(2884,30);
histRandomHitDistribution->SetBinContent(2885,33);
histRandomHitDistribution->SetBinContent(2886,37);
histRandomHitDistribution->SetBinContent(2887,49);
histRandomHitDistribution->SetBinContent(2888,41);
histRandomHitDistribution->SetBinContent(2889,55);
histRandomHitDistribution->SetBinContent(2890,49);
histRandomHitDistribution->SetBinContent(2891,55);
histRandomHitDistribution->SetBinContent(2892,42);
histRandomHitDistribution->SetBinContent(2893,57);
histRandomHitDistribution->SetBinContent(2894,56);
histRandomHitDistribution->SetBinContent(2895,63);
histRandomHitDistribution->SetBinContent(2896,60);
histRandomHitDistribution->SetBinContent(2897,66);
histRandomHitDistribution->SetBinContent(2898,56);
histRandomHitDistribution->SetBinContent(2899,57);
histRandomHitDistribution->SetBinContent(2900,55);
histRandomHitDistribution->SetBinContent(2901,60);
histRandomHitDistribution->SetBinContent(2902,59);
histRandomHitDistribution->SetBinContent(2903,41);
histRandomHitDistribution->SetBinContent(2904,37);
histRandomHitDistribution->SetBinContent(2905,51);
histRandomHitDistribution->SetBinContent(2906,45);
histRandomHitDistribution->SetBinContent(2907,44);
histRandomHitDistribution->SetBinContent(2908,36);
histRandomHitDistribution->SetBinContent(2909,37);
histRandomHitDistribution->SetBinContent(2910,41);
histRandomHitDistribution->SetBinContent(2911,34);
histRandomHitDistribution->SetBinContent(2912,29);
histRandomHitDistribution->SetBinContent(2913,22);
histRandomHitDistribution->SetBinContent(2914,29);
histRandomHitDistribution->SetBinContent(2915,23);
histRandomHitDistribution->SetBinContent(2916,19);
histRandomHitDistribution->SetBinContent(2917,13);
histRandomHitDistribution->SetBinContent(2918,12);
histRandomHitDistribution->SetBinContent(2919,19);
histRandomHitDistribution->SetBinContent(2920,12);
histRandomHitDistribution->SetBinContent(2921,5);
histRandomHitDistribution->SetBinContent(2922,8);
histRandomHitDistribution->SetBinContent(2923,8);
histRandomHitDistribution->SetBinContent(2924,3);
histRandomHitDistribution->SetBinContent(2925,5);
histRandomHitDistribution->SetBinContent(2926,4);
histRandomHitDistribution->SetBinContent(2927,4);
histRandomHitDistribution->SetBinContent(2928,5);
histRandomHitDistribution->SetBinContent(2929,5);
histRandomHitDistribution->SetBinContent(2930,3);
histRandomHitDistribution->SetBinContent(2931,2);
histRandomHitDistribution->SetBinContent(2932,1);
histRandomHitDistribution->SetBinContent(2933,4);
histRandomHitDistribution->SetBinContent(2935,2);
histRandomHitDistribution->SetBinContent(2936,2);
histRandomHitDistribution->SetBinContent(2944,1);
histRandomHitDistribution->SetBinContent(2957,1);
histRandomHitDistribution->SetBinContent(2960,16);
histRandomHitDistribution->SetBinContent(2961,107);
histRandomHitDistribution->SetBinContent(2962,223);
histRandomHitDistribution->SetBinContent(2963,275);
histRandomHitDistribution->SetBinContent(2964,328);
histRandomHitDistribution->SetBinContent(2965,285);
histRandomHitDistribution->SetBinContent(2966,248);
histRandomHitDistribution->SetBinContent(2967,177);
histRandomHitDistribution->SetBinContent(2968,107);
histRandomHitDistribution->SetBinContent(2969,87);
histRandomHitDistribution->SetBinContent(2970,68);
histRandomHitDistribution->SetBinContent(2971,32);
histRandomHitDistribution->SetBinContent(2972,17);
histRandomHitDistribution->SetBinContent(2973,12);
histRandomHitDistribution->SetBinContent(2974,7);
histRandomHitDistribution->SetBinContent(2975,3);
histRandomHitDistribution->SetBinContent(2976,8);
histRandomHitDistribution->SetBinContent(2977,7);
histRandomHitDistribution->SetBinContent(2978,7);
histRandomHitDistribution->SetBinContent(2979,11);
histRandomHitDistribution->SetBinContent(2980,18);
histRandomHitDistribution->SetBinContent(2981,16);
histRandomHitDistribution->SetBinContent(2982,23);
histRandomHitDistribution->SetBinContent(2983,22);
histRandomHitDistribution->SetBinContent(2984,31);
histRandomHitDistribution->SetBinContent(2985,33);
histRandomHitDistribution->SetBinContent(2986,38);
histRandomHitDistribution->SetBinContent(2987,28);
histRandomHitDistribution->SetBinContent(2988,48);
histRandomHitDistribution->SetBinContent(2989,49);
histRandomHitDistribution->SetBinContent(2990,51);
histRandomHitDistribution->SetBinContent(2991,52);
histRandomHitDistribution->SetBinContent(2992,49);
histRandomHitDistribution->SetBinContent(2993,55);
histRandomHitDistribution->SetBinContent(2994,51);
histRandomHitDistribution->SetBinContent(2995,46);
histRandomHitDistribution->SetBinContent(2996,67);
histRandomHitDistribution->SetBinContent(2997,75);
histRandomHitDistribution->SetBinContent(2998,67);
histRandomHitDistribution->SetBinContent(2999,66);
histRandomHitDistribution->SetBinContent(3000,53);
histRandomHitDistribution->SetBinContent(3001,53);
histRandomHitDistribution->SetBinContent(3002,59);
histRandomHitDistribution->SetBinContent(3003,49);
histRandomHitDistribution->SetBinContent(3004,70);
histRandomHitDistribution->SetBinContent(3005,55);
histRandomHitDistribution->SetBinContent(3006,50);
histRandomHitDistribution->SetBinContent(3007,48);
histRandomHitDistribution->SetBinContent(3008,52);
histRandomHitDistribution->SetBinContent(3009,45);
histRandomHitDistribution->SetBinContent(3010,43);
histRandomHitDistribution->SetBinContent(3011,42);
histRandomHitDistribution->SetBinContent(3012,46);
histRandomHitDistribution->SetBinContent(3013,28);
histRandomHitDistribution->SetBinContent(3014,23);
histRandomHitDistribution->SetBinContent(3015,37);
histRandomHitDistribution->SetBinContent(3016,27);
histRandomHitDistribution->SetBinContent(3017,25);
histRandomHitDistribution->SetBinContent(3018,23);
histRandomHitDistribution->SetBinContent(3019,23);
histRandomHitDistribution->SetBinContent(3020,16);
histRandomHitDistribution->SetBinContent(3021,18);
histRandomHitDistribution->SetBinContent(3022,11);
histRandomHitDistribution->SetBinContent(3023,10);
histRandomHitDistribution->SetBinContent(3024,5);
histRandomHitDistribution->SetBinContent(3025,12);
histRandomHitDistribution->SetBinContent(3026,8);
histRandomHitDistribution->SetBinContent(3027,7);
histRandomHitDistribution->SetBinContent(3028,5);
histRandomHitDistribution->SetBinContent(3029,2);
histRandomHitDistribution->SetBinContent(3030,3);
histRandomHitDistribution->SetBinContent(3031,4);
histRandomHitDistribution->SetBinContent(3032,3);
histRandomHitDistribution->SetBinContent(3034,2);
histRandomHitDistribution->SetBinContent(3035,1);
histRandomHitDistribution->SetBinContent(3036,1);
histRandomHitDistribution->SetBinContent(3037,1);
histRandomHitDistribution->SetBinContent(3038,1);
histRandomHitDistribution->SetBinContent(3044,1);
histRandomHitDistribution->SetBinContent(3062,11);
histRandomHitDistribution->SetBinContent(3063,89);
histRandomHitDistribution->SetBinContent(3064,228);
histRandomHitDistribution->SetBinContent(3065,319);
histRandomHitDistribution->SetBinContent(3066,291);
histRandomHitDistribution->SetBinContent(3067,333);
histRandomHitDistribution->SetBinContent(3068,255);
histRandomHitDistribution->SetBinContent(3069,175);
histRandomHitDistribution->SetBinContent(3070,151);
histRandomHitDistribution->SetBinContent(3071,100);
histRandomHitDistribution->SetBinContent(3072,56);
histRandomHitDistribution->SetBinContent(3073,32);
histRandomHitDistribution->SetBinContent(3074,22);
histRandomHitDistribution->SetBinContent(3075,8);
histRandomHitDistribution->SetBinContent(3076,4);
histRandomHitDistribution->SetBinContent(3077,5);
histRandomHitDistribution->SetBinContent(3078,2);
histRandomHitDistribution->SetBinContent(3079,7);
histRandomHitDistribution->SetBinContent(3080,9);
histRandomHitDistribution->SetBinContent(3081,7);
histRandomHitDistribution->SetBinContent(3082,15);
histRandomHitDistribution->SetBinContent(3083,23);
histRandomHitDistribution->SetBinContent(3084,23);
histRandomHitDistribution->SetBinContent(3085,32);
histRandomHitDistribution->SetBinContent(3086,26);
histRandomHitDistribution->SetBinContent(3087,35);
histRandomHitDistribution->SetBinContent(3088,39);
histRandomHitDistribution->SetBinContent(3089,40);
histRandomHitDistribution->SetBinContent(3090,50);
histRandomHitDistribution->SetBinContent(3091,49);
histRandomHitDistribution->SetBinContent(3092,58);
histRandomHitDistribution->SetBinContent(3093,47);
histRandomHitDistribution->SetBinContent(3094,74);
histRandomHitDistribution->SetBinContent(3095,75);
histRandomHitDistribution->SetBinContent(3096,74);
histRandomHitDistribution->SetBinContent(3097,54);
histRandomHitDistribution->SetBinContent(3098,66);
histRandomHitDistribution->SetBinContent(3099,58);
histRandomHitDistribution->SetBinContent(3100,76);
histRandomHitDistribution->SetBinContent(3101,63);
histRandomHitDistribution->SetBinContent(3102,57);
histRandomHitDistribution->SetBinContent(3103,52);
histRandomHitDistribution->SetBinContent(3104,46);
histRandomHitDistribution->SetBinContent(3105,64);
histRandomHitDistribution->SetBinContent(3106,57);
histRandomHitDistribution->SetBinContent(3107,72);
histRandomHitDistribution->SetBinContent(3108,56);
histRandomHitDistribution->SetBinContent(3109,51);
histRandomHitDistribution->SetBinContent(3110,42);
histRandomHitDistribution->SetBinContent(3111,45);
histRandomHitDistribution->SetBinContent(3112,35);
histRandomHitDistribution->SetBinContent(3113,42);
histRandomHitDistribution->SetBinContent(3114,36);
histRandomHitDistribution->SetBinContent(3115,29);
histRandomHitDistribution->SetBinContent(3116,27);
histRandomHitDistribution->SetBinContent(3117,25);
histRandomHitDistribution->SetBinContent(3118,24);
histRandomHitDistribution->SetBinContent(3119,27);
histRandomHitDistribution->SetBinContent(3120,23);
histRandomHitDistribution->SetBinContent(3121,12);
histRandomHitDistribution->SetBinContent(3122,13);
histRandomHitDistribution->SetBinContent(3123,10);
histRandomHitDistribution->SetBinContent(3124,9);
histRandomHitDistribution->SetBinContent(3125,11);
histRandomHitDistribution->SetBinContent(3126,9);
histRandomHitDistribution->SetBinContent(3127,7);
histRandomHitDistribution->SetBinContent(3128,9);
histRandomHitDistribution->SetBinContent(3129,6);
histRandomHitDistribution->SetBinContent(3130,7);
histRandomHitDistribution->SetBinContent(3131,6);
histRandomHitDistribution->SetBinContent(3132,3);
histRandomHitDistribution->SetBinContent(3133,2);
histRandomHitDistribution->SetBinContent(3134,1);
histRandomHitDistribution->SetBinContent(3135,5);
histRandomHitDistribution->SetBinContent(3136,1);
histRandomHitDistribution->SetBinContent(3137,2);
histRandomHitDistribution->SetBinContent(3138,2);
histRandomHitDistribution->SetBinContent(3140,1);
histRandomHitDistribution->SetBinContent(3141,1);
histRandomHitDistribution->SetBinContent(3164,14);
histRandomHitDistribution->SetBinContent(3165,69);
histRandomHitDistribution->SetBinContent(3166,217);
histRandomHitDistribution->SetBinContent(3167,312);
histRandomHitDistribution->SetBinContent(3168,342);
histRandomHitDistribution->SetBinContent(3169,326);
histRandomHitDistribution->SetBinContent(3170,239);
histRandomHitDistribution->SetBinContent(3171,209);
histRandomHitDistribution->SetBinContent(3172,149);
histRandomHitDistribution->SetBinContent(3173,78);
histRandomHitDistribution->SetBinContent(3174,48);
histRandomHitDistribution->SetBinContent(3175,26);
histRandomHitDistribution->SetBinContent(3176,14);
histRandomHitDistribution->SetBinContent(3177,10);
histRandomHitDistribution->SetBinContent(3178,5);
histRandomHitDistribution->SetBinContent(3179,3);
histRandomHitDistribution->SetBinContent(3180,5);
histRandomHitDistribution->SetBinContent(3181,3);
histRandomHitDistribution->SetBinContent(3182,14);
histRandomHitDistribution->SetBinContent(3183,12);
histRandomHitDistribution->SetBinContent(3184,12);
histRandomHitDistribution->SetBinContent(3185,15);
histRandomHitDistribution->SetBinContent(3186,19);
histRandomHitDistribution->SetBinContent(3187,23);
histRandomHitDistribution->SetBinContent(3188,42);
histRandomHitDistribution->SetBinContent(3189,32);
histRandomHitDistribution->SetBinContent(3190,41);
histRandomHitDistribution->SetBinContent(3191,44);
histRandomHitDistribution->SetBinContent(3192,43);
histRandomHitDistribution->SetBinContent(3193,38);
histRandomHitDistribution->SetBinContent(3194,44);
histRandomHitDistribution->SetBinContent(3195,44);
histRandomHitDistribution->SetBinContent(3196,52);
histRandomHitDistribution->SetBinContent(3197,37);
histRandomHitDistribution->SetBinContent(3198,56);
histRandomHitDistribution->SetBinContent(3199,69);
histRandomHitDistribution->SetBinContent(3200,48);
histRandomHitDistribution->SetBinContent(3201,71);
histRandomHitDistribution->SetBinContent(3202,68);
histRandomHitDistribution->SetBinContent(3203,59);
histRandomHitDistribution->SetBinContent(3204,67);
histRandomHitDistribution->SetBinContent(3205,69);
histRandomHitDistribution->SetBinContent(3206,50);
histRandomHitDistribution->SetBinContent(3207,63);
histRandomHitDistribution->SetBinContent(3208,67);
histRandomHitDistribution->SetBinContent(3209,45);
histRandomHitDistribution->SetBinContent(3210,55);
histRandomHitDistribution->SetBinContent(3211,60);
histRandomHitDistribution->SetBinContent(3212,58);
histRandomHitDistribution->SetBinContent(3213,42);
histRandomHitDistribution->SetBinContent(3214,52);
histRandomHitDistribution->SetBinContent(3215,37);
histRandomHitDistribution->SetBinContent(3216,33);
histRandomHitDistribution->SetBinContent(3217,33);
histRandomHitDistribution->SetBinContent(3218,31);
histRandomHitDistribution->SetBinContent(3219,31);
histRandomHitDistribution->SetBinContent(3220,25);
histRandomHitDistribution->SetBinContent(3221,21);
histRandomHitDistribution->SetBinContent(3222,23);
histRandomHitDistribution->SetBinContent(3223,17);
histRandomHitDistribution->SetBinContent(3224,16);
histRandomHitDistribution->SetBinContent(3225,16);
histRandomHitDistribution->SetBinContent(3226,11);
histRandomHitDistribution->SetBinContent(3227,14);
histRandomHitDistribution->SetBinContent(3228,12);
histRandomHitDistribution->SetBinContent(3229,16);
histRandomHitDistribution->SetBinContent(3230,10);
histRandomHitDistribution->SetBinContent(3231,7);
histRandomHitDistribution->SetBinContent(3232,9);
histRandomHitDistribution->SetBinContent(3233,7);
histRandomHitDistribution->SetBinContent(3234,1);
histRandomHitDistribution->SetBinContent(3235,2);
histRandomHitDistribution->SetBinContent(3236,3);
histRandomHitDistribution->SetBinContent(3237,4);
histRandomHitDistribution->SetBinContent(3238,3);
histRandomHitDistribution->SetBinContent(3239,4);
histRandomHitDistribution->SetBinContent(3241,1);
histRandomHitDistribution->SetBinContent(3242,1);
histRandomHitDistribution->SetBinContent(3243,2);
histRandomHitDistribution->SetBinContent(3244,2);
histRandomHitDistribution->SetBinContent(3245,1);
histRandomHitDistribution->SetBinContent(3248,1);
histRandomHitDistribution->SetBinContent(3251,3);
histRandomHitDistribution->SetBinContent(3254,1);
histRandomHitDistribution->SetBinContent(3266,13);
histRandomHitDistribution->SetBinContent(3267,93);
histRandomHitDistribution->SetBinContent(3268,189);
histRandomHitDistribution->SetBinContent(3269,316);
histRandomHitDistribution->SetBinContent(3270,327);
histRandomHitDistribution->SetBinContent(3271,320);
histRandomHitDistribution->SetBinContent(3272,260);
histRandomHitDistribution->SetBinContent(3273,158);
histRandomHitDistribution->SetBinContent(3274,132);
histRandomHitDistribution->SetBinContent(3275,77);
histRandomHitDistribution->SetBinContent(3276,52);
histRandomHitDistribution->SetBinContent(3277,26);
histRandomHitDistribution->SetBinContent(3278,19);
histRandomHitDistribution->SetBinContent(3279,13);
histRandomHitDistribution->SetBinContent(3280,7);
histRandomHitDistribution->SetBinContent(3281,4);
histRandomHitDistribution->SetBinContent(3282,5);
histRandomHitDistribution->SetBinContent(3283,1);
histRandomHitDistribution->SetBinContent(3284,7);
histRandomHitDistribution->SetBinContent(3285,9);
histRandomHitDistribution->SetBinContent(3286,7);
histRandomHitDistribution->SetBinContent(3287,12);
histRandomHitDistribution->SetBinContent(3288,21);
histRandomHitDistribution->SetBinContent(3289,26);
histRandomHitDistribution->SetBinContent(3290,24);
histRandomHitDistribution->SetBinContent(3291,28);
histRandomHitDistribution->SetBinContent(3292,34);
histRandomHitDistribution->SetBinContent(3293,33);
histRandomHitDistribution->SetBinContent(3294,43);
histRandomHitDistribution->SetBinContent(3295,35);
histRandomHitDistribution->SetBinContent(3296,43);
histRandomHitDistribution->SetBinContent(3297,39);
histRandomHitDistribution->SetBinContent(3298,41);
histRandomHitDistribution->SetBinContent(3299,41);
histRandomHitDistribution->SetBinContent(3300,53);
histRandomHitDistribution->SetBinContent(3301,46);
histRandomHitDistribution->SetBinContent(3302,36);
histRandomHitDistribution->SetBinContent(3303,47);
histRandomHitDistribution->SetBinContent(3304,55);
histRandomHitDistribution->SetBinContent(3305,40);
histRandomHitDistribution->SetBinContent(3306,53);
histRandomHitDistribution->SetBinContent(3307,58);
histRandomHitDistribution->SetBinContent(3308,46);
histRandomHitDistribution->SetBinContent(3309,40);
histRandomHitDistribution->SetBinContent(3310,43);
histRandomHitDistribution->SetBinContent(3311,40);
histRandomHitDistribution->SetBinContent(3312,52);
histRandomHitDistribution->SetBinContent(3313,33);
histRandomHitDistribution->SetBinContent(3314,36);
histRandomHitDistribution->SetBinContent(3315,34);
histRandomHitDistribution->SetBinContent(3316,36);
histRandomHitDistribution->SetBinContent(3317,43);
histRandomHitDistribution->SetBinContent(3318,47);
histRandomHitDistribution->SetBinContent(3319,34);
histRandomHitDistribution->SetBinContent(3320,26);
histRandomHitDistribution->SetBinContent(3321,19);
histRandomHitDistribution->SetBinContent(3322,22);
histRandomHitDistribution->SetBinContent(3323,26);
histRandomHitDistribution->SetBinContent(3324,20);
histRandomHitDistribution->SetBinContent(3325,13);
histRandomHitDistribution->SetBinContent(3326,12);
histRandomHitDistribution->SetBinContent(3327,18);
histRandomHitDistribution->SetBinContent(3328,11);
histRandomHitDistribution->SetBinContent(3329,18);
histRandomHitDistribution->SetBinContent(3330,8);
histRandomHitDistribution->SetBinContent(3331,7);
histRandomHitDistribution->SetBinContent(3332,12);
histRandomHitDistribution->SetBinContent(3333,5);
histRandomHitDistribution->SetBinContent(3334,8);
histRandomHitDistribution->SetBinContent(3335,2);
histRandomHitDistribution->SetBinContent(3336,4);
histRandomHitDistribution->SetBinContent(3337,2);
histRandomHitDistribution->SetBinContent(3338,1);
histRandomHitDistribution->SetBinContent(3339,1);
histRandomHitDistribution->SetBinContent(3340,1);
histRandomHitDistribution->SetBinContent(3341,1);
histRandomHitDistribution->SetBinContent(3343,2);
histRandomHitDistribution->SetBinContent(3344,1);
histRandomHitDistribution->SetBinContent(3345,1);
histRandomHitDistribution->SetBinContent(3349,1);
histRandomHitDistribution->SetBinContent(3352,1);
histRandomHitDistribution->SetBinContent(3368,9);
histRandomHitDistribution->SetBinContent(3369,85);
histRandomHitDistribution->SetBinContent(3370,207);
histRandomHitDistribution->SetBinContent(3371,316);
histRandomHitDistribution->SetBinContent(3372,370);
histRandomHitDistribution->SetBinContent(3373,308);
histRandomHitDistribution->SetBinContent(3374,258);
histRandomHitDistribution->SetBinContent(3375,169);
histRandomHitDistribution->SetBinContent(3376,101);
histRandomHitDistribution->SetBinContent(3377,71);
histRandomHitDistribution->SetBinContent(3378,41);
histRandomHitDistribution->SetBinContent(3379,30);
histRandomHitDistribution->SetBinContent(3380,24);
histRandomHitDistribution->SetBinContent(3381,8);
histRandomHitDistribution->SetBinContent(3382,4);
histRandomHitDistribution->SetBinContent(3383,2);
histRandomHitDistribution->SetBinContent(3384,6);
histRandomHitDistribution->SetBinContent(3385,2);
histRandomHitDistribution->SetBinContent(3386,6);
histRandomHitDistribution->SetBinContent(3387,5);
histRandomHitDistribution->SetBinContent(3388,10);
histRandomHitDistribution->SetBinContent(3389,13);
histRandomHitDistribution->SetBinContent(3390,13);
histRandomHitDistribution->SetBinContent(3391,23);
histRandomHitDistribution->SetBinContent(3392,17);
histRandomHitDistribution->SetBinContent(3393,25);
histRandomHitDistribution->SetBinContent(3394,23);
histRandomHitDistribution->SetBinContent(3395,30);
histRandomHitDistribution->SetBinContent(3396,41);
histRandomHitDistribution->SetBinContent(3397,39);
histRandomHitDistribution->SetBinContent(3398,34);
histRandomHitDistribution->SetBinContent(3399,41);
histRandomHitDistribution->SetBinContent(3400,42);
histRandomHitDistribution->SetBinContent(3401,51);
histRandomHitDistribution->SetBinContent(3402,44);
histRandomHitDistribution->SetBinContent(3403,43);
histRandomHitDistribution->SetBinContent(3404,39);
histRandomHitDistribution->SetBinContent(3405,52);
histRandomHitDistribution->SetBinContent(3406,44);
histRandomHitDistribution->SetBinContent(3407,46);
histRandomHitDistribution->SetBinContent(3408,29);
histRandomHitDistribution->SetBinContent(3409,49);
histRandomHitDistribution->SetBinContent(3410,50);
histRandomHitDistribution->SetBinContent(3411,53);
histRandomHitDistribution->SetBinContent(3412,50);
histRandomHitDistribution->SetBinContent(3413,30);
histRandomHitDistribution->SetBinContent(3414,44);
histRandomHitDistribution->SetBinContent(3415,48);
histRandomHitDistribution->SetBinContent(3416,39);
histRandomHitDistribution->SetBinContent(3417,33);
histRandomHitDistribution->SetBinContent(3418,38);
histRandomHitDistribution->SetBinContent(3419,36);
histRandomHitDistribution->SetBinContent(3420,36);
histRandomHitDistribution->SetBinContent(3421,30);
histRandomHitDistribution->SetBinContent(3422,26);
histRandomHitDistribution->SetBinContent(3423,27);
histRandomHitDistribution->SetBinContent(3424,19);
histRandomHitDistribution->SetBinContent(3425,23);
histRandomHitDistribution->SetBinContent(3426,19);
histRandomHitDistribution->SetBinContent(3427,13);
histRandomHitDistribution->SetBinContent(3428,11);
histRandomHitDistribution->SetBinContent(3429,16);
histRandomHitDistribution->SetBinContent(3430,6);
histRandomHitDistribution->SetBinContent(3431,9);
histRandomHitDistribution->SetBinContent(3432,7);
histRandomHitDistribution->SetBinContent(3433,5);
histRandomHitDistribution->SetBinContent(3434,8);
histRandomHitDistribution->SetBinContent(3435,3);
histRandomHitDistribution->SetBinContent(3436,5);
histRandomHitDistribution->SetBinContent(3437,7);
histRandomHitDistribution->SetBinContent(3438,1);
histRandomHitDistribution->SetBinContent(3439,3);
histRandomHitDistribution->SetBinContent(3440,5);
histRandomHitDistribution->SetBinContent(3441,2);
histRandomHitDistribution->SetBinContent(3442,1);
histRandomHitDistribution->SetBinContent(3444,2);
histRandomHitDistribution->SetBinContent(3445,1);
histRandomHitDistribution->SetBinContent(3447,2);
histRandomHitDistribution->SetBinContent(3448,2);
histRandomHitDistribution->SetBinContent(3452,1);
histRandomHitDistribution->SetBinContent(3453,2);
histRandomHitDistribution->SetBinContent(3470,13);
histRandomHitDistribution->SetBinContent(3471,80);
histRandomHitDistribution->SetBinContent(3472,199);
histRandomHitDistribution->SetBinContent(3473,300);
histRandomHitDistribution->SetBinContent(3474,334);
histRandomHitDistribution->SetBinContent(3475,329);
histRandomHitDistribution->SetBinContent(3476,264);
histRandomHitDistribution->SetBinContent(3477,164);
histRandomHitDistribution->SetBinContent(3478,95);
histRandomHitDistribution->SetBinContent(3479,43);
histRandomHitDistribution->SetBinContent(3480,39);
histRandomHitDistribution->SetBinContent(3481,31);
histRandomHitDistribution->SetBinContent(3482,18);
histRandomHitDistribution->SetBinContent(3483,9);
histRandomHitDistribution->SetBinContent(3485,3);
histRandomHitDistribution->SetBinContent(3486,2);
histRandomHitDistribution->SetBinContent(3487,5);
histRandomHitDistribution->SetBinContent(3488,4);
histRandomHitDistribution->SetBinContent(3489,5);
histRandomHitDistribution->SetBinContent(3490,10);
histRandomHitDistribution->SetBinContent(3491,20);
histRandomHitDistribution->SetBinContent(3492,18);
histRandomHitDistribution->SetBinContent(3493,22);
histRandomHitDistribution->SetBinContent(3494,17);
histRandomHitDistribution->SetBinContent(3495,29);
histRandomHitDistribution->SetBinContent(3496,20);
histRandomHitDistribution->SetBinContent(3497,26);
histRandomHitDistribution->SetBinContent(3498,30);
histRandomHitDistribution->SetBinContent(3499,33);
histRandomHitDistribution->SetBinContent(3500,38);
histRandomHitDistribution->SetBinContent(3501,27);
histRandomHitDistribution->SetBinContent(3502,37);
histRandomHitDistribution->SetBinContent(3503,29);
histRandomHitDistribution->SetBinContent(3504,45);
histRandomHitDistribution->SetBinContent(3505,32);
histRandomHitDistribution->SetBinContent(3506,49);
histRandomHitDistribution->SetBinContent(3507,47);
histRandomHitDistribution->SetBinContent(3508,50);
histRandomHitDistribution->SetBinContent(3509,56);
histRandomHitDistribution->SetBinContent(3510,41);
histRandomHitDistribution->SetBinContent(3511,41);
histRandomHitDistribution->SetBinContent(3512,41);
histRandomHitDistribution->SetBinContent(3513,47);
histRandomHitDistribution->SetBinContent(3514,45);
histRandomHitDistribution->SetBinContent(3515,34);
histRandomHitDistribution->SetBinContent(3516,34);
histRandomHitDistribution->SetBinContent(3517,30);
histRandomHitDistribution->SetBinContent(3518,29);
histRandomHitDistribution->SetBinContent(3519,30);
histRandomHitDistribution->SetBinContent(3520,28);
histRandomHitDistribution->SetBinContent(3521,34);
histRandomHitDistribution->SetBinContent(3522,22);
histRandomHitDistribution->SetBinContent(3523,20);
histRandomHitDistribution->SetBinContent(3524,14);
histRandomHitDistribution->SetBinContent(3525,24);
histRandomHitDistribution->SetBinContent(3526,16);
histRandomHitDistribution->SetBinContent(3527,22);
histRandomHitDistribution->SetBinContent(3528,12);
histRandomHitDistribution->SetBinContent(3529,11);
histRandomHitDistribution->SetBinContent(3530,11);
histRandomHitDistribution->SetBinContent(3531,15);
histRandomHitDistribution->SetBinContent(3532,8);
histRandomHitDistribution->SetBinContent(3533,9);
histRandomHitDistribution->SetBinContent(3534,10);
histRandomHitDistribution->SetBinContent(3535,9);
histRandomHitDistribution->SetBinContent(3536,10);
histRandomHitDistribution->SetBinContent(3537,7);
histRandomHitDistribution->SetBinContent(3538,6);
histRandomHitDistribution->SetBinContent(3539,3);
histRandomHitDistribution->SetBinContent(3540,7);
histRandomHitDistribution->SetBinContent(3541,1);
histRandomHitDistribution->SetBinContent(3542,3);
histRandomHitDistribution->SetBinContent(3543,3);
histRandomHitDistribution->SetBinContent(3544,3);
histRandomHitDistribution->SetBinContent(3546,1);
histRandomHitDistribution->SetBinContent(3547,1);
histRandomHitDistribution->SetBinContent(3548,1);
histRandomHitDistribution->SetBinContent(3549,1);
histRandomHitDistribution->SetBinContent(3559,1);
histRandomHitDistribution->SetBinContent(3560,1);
histRandomHitDistribution->SetBinContent(3572,12);
histRandomHitDistribution->SetBinContent(3573,103);
histRandomHitDistribution->SetBinContent(3574,182);
histRandomHitDistribution->SetBinContent(3575,291);
histRandomHitDistribution->SetBinContent(3576,306);
histRandomHitDistribution->SetBinContent(3577,312);
histRandomHitDistribution->SetBinContent(3578,221);
histRandomHitDistribution->SetBinContent(3579,172);
histRandomHitDistribution->SetBinContent(3580,115);
histRandomHitDistribution->SetBinContent(3581,76);
histRandomHitDistribution->SetBinContent(3582,47);
histRandomHitDistribution->SetBinContent(3583,48);
histRandomHitDistribution->SetBinContent(3584,20);
histRandomHitDistribution->SetBinContent(3585,7);
histRandomHitDistribution->SetBinContent(3586,8);
histRandomHitDistribution->SetBinContent(3587,4);
histRandomHitDistribution->SetBinContent(3588,3);
histRandomHitDistribution->SetBinContent(3589,4);
histRandomHitDistribution->SetBinContent(3590,2);
histRandomHitDistribution->SetBinContent(3591,11);
histRandomHitDistribution->SetBinContent(3592,10);
histRandomHitDistribution->SetBinContent(3593,9);
histRandomHitDistribution->SetBinContent(3594,14);
histRandomHitDistribution->SetBinContent(3595,14);
histRandomHitDistribution->SetBinContent(3596,13);
histRandomHitDistribution->SetBinContent(3597,15);
histRandomHitDistribution->SetBinContent(3598,21);
histRandomHitDistribution->SetBinContent(3599,16);
histRandomHitDistribution->SetBinContent(3600,20);
histRandomHitDistribution->SetBinContent(3601,32);
histRandomHitDistribution->SetBinContent(3602,26);
histRandomHitDistribution->SetBinContent(3603,18);
histRandomHitDistribution->SetBinContent(3604,29);
histRandomHitDistribution->SetBinContent(3605,32);
histRandomHitDistribution->SetBinContent(3606,24);
histRandomHitDistribution->SetBinContent(3607,31);
histRandomHitDistribution->SetBinContent(3608,27);
histRandomHitDistribution->SetBinContent(3609,28);
histRandomHitDistribution->SetBinContent(3610,34);
histRandomHitDistribution->SetBinContent(3611,38);
histRandomHitDistribution->SetBinContent(3612,29);
histRandomHitDistribution->SetBinContent(3613,36);
histRandomHitDistribution->SetBinContent(3614,30);
histRandomHitDistribution->SetBinContent(3615,17);
histRandomHitDistribution->SetBinContent(3616,34);
histRandomHitDistribution->SetBinContent(3617,23);
histRandomHitDistribution->SetBinContent(3618,30);
histRandomHitDistribution->SetBinContent(3619,26);
histRandomHitDistribution->SetBinContent(3620,24);
histRandomHitDistribution->SetBinContent(3621,27);
histRandomHitDistribution->SetBinContent(3622,23);
histRandomHitDistribution->SetBinContent(3623,21);
histRandomHitDistribution->SetBinContent(3624,18);
histRandomHitDistribution->SetBinContent(3625,15);
histRandomHitDistribution->SetBinContent(3626,13);
histRandomHitDistribution->SetBinContent(3627,10);
histRandomHitDistribution->SetBinContent(3628,15);
histRandomHitDistribution->SetBinContent(3629,14);
histRandomHitDistribution->SetBinContent(3630,10);
histRandomHitDistribution->SetBinContent(3631,13);
histRandomHitDistribution->SetBinContent(3632,13);
histRandomHitDistribution->SetBinContent(3633,8);
histRandomHitDistribution->SetBinContent(3634,1);
histRandomHitDistribution->SetBinContent(3635,9);
histRandomHitDistribution->SetBinContent(3636,5);
histRandomHitDistribution->SetBinContent(3637,4);
histRandomHitDistribution->SetBinContent(3638,6);
histRandomHitDistribution->SetBinContent(3639,4);
histRandomHitDistribution->SetBinContent(3641,4);
histRandomHitDistribution->SetBinContent(3642,1);
histRandomHitDistribution->SetBinContent(3643,3);
histRandomHitDistribution->SetBinContent(3644,1);
histRandomHitDistribution->SetBinContent(3645,4);
histRandomHitDistribution->SetBinContent(3646,2);
histRandomHitDistribution->SetBinContent(3647,2);
histRandomHitDistribution->SetBinContent(3649,1);
histRandomHitDistribution->SetBinContent(3652,1);
histRandomHitDistribution->SetBinContent(3653,2);
histRandomHitDistribution->SetBinContent(3654,1);
histRandomHitDistribution->SetBinContent(3656,2);
histRandomHitDistribution->SetBinContent(3657,3);
histRandomHitDistribution->SetBinContent(3664,1);
histRandomHitDistribution->SetBinContent(3667,1);
histRandomHitDistribution->SetBinContent(3671,1);
histRandomHitDistribution->SetBinContent(3674,21);
histRandomHitDistribution->SetBinContent(3675,87);
histRandomHitDistribution->SetBinContent(3676,186);
histRandomHitDistribution->SetBinContent(3677,307);
histRandomHitDistribution->SetBinContent(3678,323);
histRandomHitDistribution->SetBinContent(3679,281);
histRandomHitDistribution->SetBinContent(3680,274);
histRandomHitDistribution->SetBinContent(3681,162);
histRandomHitDistribution->SetBinContent(3682,127);
histRandomHitDistribution->SetBinContent(3683,80);
histRandomHitDistribution->SetBinContent(3684,52);
histRandomHitDistribution->SetBinContent(3685,38);
histRandomHitDistribution->SetBinContent(3686,20);
histRandomHitDistribution->SetBinContent(3687,13);
histRandomHitDistribution->SetBinContent(3688,5);
histRandomHitDistribution->SetBinContent(3689,3);
histRandomHitDistribution->SetBinContent(3690,7);
histRandomHitDistribution->SetBinContent(3691,4);
histRandomHitDistribution->SetBinContent(3692,6);
histRandomHitDistribution->SetBinContent(3693,6);
histRandomHitDistribution->SetBinContent(3694,12);
histRandomHitDistribution->SetBinContent(3695,10);
histRandomHitDistribution->SetBinContent(3696,9);
histRandomHitDistribution->SetBinContent(3697,12);
histRandomHitDistribution->SetBinContent(3698,11);
histRandomHitDistribution->SetBinContent(3699,11);
histRandomHitDistribution->SetBinContent(3700,20);
histRandomHitDistribution->SetBinContent(3701,10);
histRandomHitDistribution->SetBinContent(3702,28);
histRandomHitDistribution->SetBinContent(3703,28);
histRandomHitDistribution->SetBinContent(3704,26);
histRandomHitDistribution->SetBinContent(3705,17);
histRandomHitDistribution->SetBinContent(3706,18);
histRandomHitDistribution->SetBinContent(3707,30);
histRandomHitDistribution->SetBinContent(3708,27);
histRandomHitDistribution->SetBinContent(3709,30);
histRandomHitDistribution->SetBinContent(3710,29);
histRandomHitDistribution->SetBinContent(3711,28);
histRandomHitDistribution->SetBinContent(3712,28);
histRandomHitDistribution->SetBinContent(3713,30);
histRandomHitDistribution->SetBinContent(3714,19);
histRandomHitDistribution->SetBinContent(3715,25);
histRandomHitDistribution->SetBinContent(3716,29);
histRandomHitDistribution->SetBinContent(3717,20);
histRandomHitDistribution->SetBinContent(3718,26);
histRandomHitDistribution->SetBinContent(3719,24);
histRandomHitDistribution->SetBinContent(3720,31);
histRandomHitDistribution->SetBinContent(3721,23);
histRandomHitDistribution->SetBinContent(3722,19);
histRandomHitDistribution->SetBinContent(3723,19);
histRandomHitDistribution->SetBinContent(3724,11);
histRandomHitDistribution->SetBinContent(3725,25);
histRandomHitDistribution->SetBinContent(3726,18);
histRandomHitDistribution->SetBinContent(3727,14);
histRandomHitDistribution->SetBinContent(3728,10);
histRandomHitDistribution->SetBinContent(3729,12);
histRandomHitDistribution->SetBinContent(3730,5);
histRandomHitDistribution->SetBinContent(3731,10);
histRandomHitDistribution->SetBinContent(3732,5);
histRandomHitDistribution->SetBinContent(3733,8);
histRandomHitDistribution->SetBinContent(3734,5);
histRandomHitDistribution->SetBinContent(3735,4);
histRandomHitDistribution->SetBinContent(3736,7);
histRandomHitDistribution->SetBinContent(3737,2);
histRandomHitDistribution->SetBinContent(3738,4);
histRandomHitDistribution->SetBinContent(3739,3);
histRandomHitDistribution->SetBinContent(3740,4);
histRandomHitDistribution->SetBinContent(3741,4);
histRandomHitDistribution->SetBinContent(3742,2);
histRandomHitDistribution->SetBinContent(3743,2);
histRandomHitDistribution->SetBinContent(3745,1);
histRandomHitDistribution->SetBinContent(3746,2);
histRandomHitDistribution->SetBinContent(3747,2);
histRandomHitDistribution->SetBinContent(3749,1);
histRandomHitDistribution->SetBinContent(3751,1);
histRandomHitDistribution->SetBinContent(3752,1);
histRandomHitDistribution->SetBinContent(3753,3);
histRandomHitDistribution->SetBinContent(3756,1);
histRandomHitDistribution->SetBinContent(3758,1);
histRandomHitDistribution->SetBinContent(3759,1);
histRandomHitDistribution->SetBinContent(3776,20);
histRandomHitDistribution->SetBinContent(3777,92);
histRandomHitDistribution->SetBinContent(3778,206);
histRandomHitDistribution->SetBinContent(3779,237);
histRandomHitDistribution->SetBinContent(3780,283);
histRandomHitDistribution->SetBinContent(3781,265);
histRandomHitDistribution->SetBinContent(3782,234);
histRandomHitDistribution->SetBinContent(3783,174);
histRandomHitDistribution->SetBinContent(3784,83);
histRandomHitDistribution->SetBinContent(3785,66);
histRandomHitDistribution->SetBinContent(3786,49);
histRandomHitDistribution->SetBinContent(3787,24);
histRandomHitDistribution->SetBinContent(3788,8);
histRandomHitDistribution->SetBinContent(3789,12);
histRandomHitDistribution->SetBinContent(3790,5);
histRandomHitDistribution->SetBinContent(3791,4);
histRandomHitDistribution->SetBinContent(3792,1);
histRandomHitDistribution->SetBinContent(3793,3);
histRandomHitDistribution->SetBinContent(3794,1);
histRandomHitDistribution->SetBinContent(3795,4);
histRandomHitDistribution->SetBinContent(3796,1);
histRandomHitDistribution->SetBinContent(3797,5);
histRandomHitDistribution->SetBinContent(3798,5);
histRandomHitDistribution->SetBinContent(3799,5);
histRandomHitDistribution->SetBinContent(3800,6);
histRandomHitDistribution->SetBinContent(3801,11);
histRandomHitDistribution->SetBinContent(3802,16);
histRandomHitDistribution->SetBinContent(3803,11);
histRandomHitDistribution->SetBinContent(3804,17);
histRandomHitDistribution->SetBinContent(3805,17);
histRandomHitDistribution->SetBinContent(3806,17);
histRandomHitDistribution->SetBinContent(3807,19);
histRandomHitDistribution->SetBinContent(3808,14);
histRandomHitDistribution->SetBinContent(3809,25);
histRandomHitDistribution->SetBinContent(3810,28);
histRandomHitDistribution->SetBinContent(3811,22);
histRandomHitDistribution->SetBinContent(3812,19);
histRandomHitDistribution->SetBinContent(3813,24);
histRandomHitDistribution->SetBinContent(3814,20);
histRandomHitDistribution->SetBinContent(3815,10);
histRandomHitDistribution->SetBinContent(3816,21);
histRandomHitDistribution->SetBinContent(3817,19);
histRandomHitDistribution->SetBinContent(3818,13);
histRandomHitDistribution->SetBinContent(3819,9);
histRandomHitDistribution->SetBinContent(3820,14);
histRandomHitDistribution->SetBinContent(3821,19);
histRandomHitDistribution->SetBinContent(3822,12);
histRandomHitDistribution->SetBinContent(3823,11);
histRandomHitDistribution->SetBinContent(3824,15);
histRandomHitDistribution->SetBinContent(3825,7);
histRandomHitDistribution->SetBinContent(3826,7);
histRandomHitDistribution->SetBinContent(3827,8);
histRandomHitDistribution->SetBinContent(3828,8);
histRandomHitDistribution->SetBinContent(3829,8);
histRandomHitDistribution->SetBinContent(3830,4);
histRandomHitDistribution->SetBinContent(3831,7);
histRandomHitDistribution->SetBinContent(3832,7);
histRandomHitDistribution->SetBinContent(3833,6);
histRandomHitDistribution->SetBinContent(3834,5);
histRandomHitDistribution->SetBinContent(3835,5);
histRandomHitDistribution->SetBinContent(3836,4);
histRandomHitDistribution->SetBinContent(3837,4);
histRandomHitDistribution->SetBinContent(3838,5);
histRandomHitDistribution->SetBinContent(3839,4);
histRandomHitDistribution->SetBinContent(3840,6);
histRandomHitDistribution->SetBinContent(3841,3);
histRandomHitDistribution->SetBinContent(3842,5);
histRandomHitDistribution->SetBinContent(3843,2);
histRandomHitDistribution->SetBinContent(3844,1);
histRandomHitDistribution->SetBinContent(3847,2);
histRandomHitDistribution->SetBinContent(3848,2);
histRandomHitDistribution->SetBinContent(3850,1);
histRandomHitDistribution->SetBinContent(3853,1);
histRandomHitDistribution->SetBinContent(3878,23);
histRandomHitDistribution->SetBinContent(3879,108);
histRandomHitDistribution->SetBinContent(3880,217);
histRandomHitDistribution->SetBinContent(3881,252);
histRandomHitDistribution->SetBinContent(3882,255);
histRandomHitDistribution->SetBinContent(3883,290);
histRandomHitDistribution->SetBinContent(3884,231);
histRandomHitDistribution->SetBinContent(3885,166);
histRandomHitDistribution->SetBinContent(3886,117);
histRandomHitDistribution->SetBinContent(3887,57);
histRandomHitDistribution->SetBinContent(3888,31);
histRandomHitDistribution->SetBinContent(3889,22);
histRandomHitDistribution->SetBinContent(3890,12);
histRandomHitDistribution->SetBinContent(3891,3);
histRandomHitDistribution->SetBinContent(3892,1);
histRandomHitDistribution->SetBinContent(3893,4);
histRandomHitDistribution->SetBinContent(3894,2);
histRandomHitDistribution->SetBinContent(3895,3);
histRandomHitDistribution->SetBinContent(3896,1);
histRandomHitDistribution->SetBinContent(3897,3);
histRandomHitDistribution->SetBinContent(3898,5);
histRandomHitDistribution->SetBinContent(3899,4);
histRandomHitDistribution->SetBinContent(3900,4);
histRandomHitDistribution->SetBinContent(3901,8);
histRandomHitDistribution->SetBinContent(3902,8);
histRandomHitDistribution->SetBinContent(3903,3);
histRandomHitDistribution->SetBinContent(3904,10);
histRandomHitDistribution->SetBinContent(3905,6);
histRandomHitDistribution->SetBinContent(3906,3);
histRandomHitDistribution->SetBinContent(3907,13);
histRandomHitDistribution->SetBinContent(3908,13);
histRandomHitDistribution->SetBinContent(3909,4);
histRandomHitDistribution->SetBinContent(3910,10);
histRandomHitDistribution->SetBinContent(3911,15);
histRandomHitDistribution->SetBinContent(3912,8);
histRandomHitDistribution->SetBinContent(3913,14);
histRandomHitDistribution->SetBinContent(3914,9);
histRandomHitDistribution->SetBinContent(3915,10);
histRandomHitDistribution->SetBinContent(3916,11);
histRandomHitDistribution->SetBinContent(3917,14);
histRandomHitDistribution->SetBinContent(3918,12);
histRandomHitDistribution->SetBinContent(3919,8);
histRandomHitDistribution->SetBinContent(3920,11);
histRandomHitDistribution->SetBinContent(3921,16);
histRandomHitDistribution->SetBinContent(3922,12);
histRandomHitDistribution->SetBinContent(3923,10);
histRandomHitDistribution->SetBinContent(3924,5);
histRandomHitDistribution->SetBinContent(3925,8);
histRandomHitDistribution->SetBinContent(3926,11);
histRandomHitDistribution->SetBinContent(3927,4);
histRandomHitDistribution->SetBinContent(3928,8);
histRandomHitDistribution->SetBinContent(3929,5);
histRandomHitDistribution->SetBinContent(3930,11);
histRandomHitDistribution->SetBinContent(3931,5);
histRandomHitDistribution->SetBinContent(3932,10);
histRandomHitDistribution->SetBinContent(3933,1);
histRandomHitDistribution->SetBinContent(3934,4);
histRandomHitDistribution->SetBinContent(3935,4);
histRandomHitDistribution->SetBinContent(3936,2);
histRandomHitDistribution->SetBinContent(3937,3);
histRandomHitDistribution->SetBinContent(3938,3);
histRandomHitDistribution->SetBinContent(3939,3);
histRandomHitDistribution->SetBinContent(3940,4);
histRandomHitDistribution->SetBinContent(3941,3);
histRandomHitDistribution->SetBinContent(3942,1);
histRandomHitDistribution->SetBinContent(3943,2);
histRandomHitDistribution->SetBinContent(3945,2);
histRandomHitDistribution->SetBinContent(3946,1);
histRandomHitDistribution->SetBinContent(3947,1);
histRandomHitDistribution->SetBinContent(3955,2);
histRandomHitDistribution->SetBinContent(3979,1);
histRandomHitDistribution->SetBinContent(3980,8);
histRandomHitDistribution->SetBinContent(3981,88);
histRandomHitDistribution->SetBinContent(3982,189);
histRandomHitDistribution->SetBinContent(3983,263);
histRandomHitDistribution->SetBinContent(3984,258);
histRandomHitDistribution->SetBinContent(3985,276);
histRandomHitDistribution->SetBinContent(3986,219);
histRandomHitDistribution->SetBinContent(3987,148);
histRandomHitDistribution->SetBinContent(3988,90);
histRandomHitDistribution->SetBinContent(3989,57);
histRandomHitDistribution->SetBinContent(3990,34);
histRandomHitDistribution->SetBinContent(3991,14);
histRandomHitDistribution->SetBinContent(3992,10);
histRandomHitDistribution->SetBinContent(3993,3);
histRandomHitDistribution->SetBinContent(3994,4);
histRandomHitDistribution->SetBinContent(3998,2);
histRandomHitDistribution->SetBinContent(3999,1);
histRandomHitDistribution->SetBinContent(4000,6);
histRandomHitDistribution->SetBinContent(4001,2);
histRandomHitDistribution->SetBinContent(4002,4);
histRandomHitDistribution->SetBinContent(4003,2);
histRandomHitDistribution->SetBinContent(4004,5);
histRandomHitDistribution->SetBinContent(4005,9);
histRandomHitDistribution->SetBinContent(4006,10);
histRandomHitDistribution->SetBinContent(4007,11);
histRandomHitDistribution->SetBinContent(4008,6);
histRandomHitDistribution->SetBinContent(4009,5);
histRandomHitDistribution->SetBinContent(4010,9);
histRandomHitDistribution->SetBinContent(4011,10);
histRandomHitDistribution->SetBinContent(4012,6);
histRandomHitDistribution->SetBinContent(4013,8);
histRandomHitDistribution->SetBinContent(4014,12);
histRandomHitDistribution->SetBinContent(4015,6);
histRandomHitDistribution->SetBinContent(4016,4);
histRandomHitDistribution->SetBinContent(4017,8);
histRandomHitDistribution->SetBinContent(4018,7);
histRandomHitDistribution->SetBinContent(4019,10);
histRandomHitDistribution->SetBinContent(4020,8);
histRandomHitDistribution->SetBinContent(4021,9);
histRandomHitDistribution->SetBinContent(4022,8);
histRandomHitDistribution->SetBinContent(4023,5);
histRandomHitDistribution->SetBinContent(4024,3);
histRandomHitDistribution->SetBinContent(4025,7);
histRandomHitDistribution->SetBinContent(4026,7);
histRandomHitDistribution->SetBinContent(4027,4);
histRandomHitDistribution->SetBinContent(4028,1);
histRandomHitDistribution->SetBinContent(4029,7);
histRandomHitDistribution->SetBinContent(4030,2);
histRandomHitDistribution->SetBinContent(4031,3);
histRandomHitDistribution->SetBinContent(4032,4);
histRandomHitDistribution->SetBinContent(4033,2);
histRandomHitDistribution->SetBinContent(4034,3);
histRandomHitDistribution->SetBinContent(4035,6);
histRandomHitDistribution->SetBinContent(4036,2);
histRandomHitDistribution->SetBinContent(4037,1);
histRandomHitDistribution->SetBinContent(4038,2);
histRandomHitDistribution->SetBinContent(4039,1);
histRandomHitDistribution->SetBinContent(4040,1);
histRandomHitDistribution->SetBinContent(4041,1);
histRandomHitDistribution->SetBinContent(4042,1);
histRandomHitDistribution->SetBinContent(4046,1);
histRandomHitDistribution->SetBinContent(4082,9);
histRandomHitDistribution->SetBinContent(4083,87);
histRandomHitDistribution->SetBinContent(4084,154);
histRandomHitDistribution->SetBinContent(4085,213);
histRandomHitDistribution->SetBinContent(4086,302);
histRandomHitDistribution->SetBinContent(4087,273);
histRandomHitDistribution->SetBinContent(4088,176);
histRandomHitDistribution->SetBinContent(4089,122);
histRandomHitDistribution->SetBinContent(4090,74);
histRandomHitDistribution->SetBinContent(4091,52);
histRandomHitDistribution->SetBinContent(4092,41);
histRandomHitDistribution->SetBinContent(4093,18);
histRandomHitDistribution->SetBinContent(4094,14);
histRandomHitDistribution->SetBinContent(4095,2);
histRandomHitDistribution->SetBinContent(4096,1);
histRandomHitDistribution->SetBinContent(4097,3);
histRandomHitDistribution->SetBinContent(4098,1);
histRandomHitDistribution->SetBinContent(4100,1);
histRandomHitDistribution->SetBinContent(4102,1);
histRandomHitDistribution->SetBinContent(4103,1);
histRandomHitDistribution->SetBinContent(4104,5);
histRandomHitDistribution->SetBinContent(4105,4);
histRandomHitDistribution->SetBinContent(4106,4);
histRandomHitDistribution->SetBinContent(4107,3);
histRandomHitDistribution->SetBinContent(4108,5);
histRandomHitDistribution->SetBinContent(4109,7);
histRandomHitDistribution->SetBinContent(4110,6);
histRandomHitDistribution->SetBinContent(4111,5);
histRandomHitDistribution->SetBinContent(4112,4);
histRandomHitDistribution->SetBinContent(4113,2);
histRandomHitDistribution->SetBinContent(4114,6);
histRandomHitDistribution->SetBinContent(4115,6);
histRandomHitDistribution->SetBinContent(4116,4);
histRandomHitDistribution->SetBinContent(4117,6);
histRandomHitDistribution->SetBinContent(4118,3);
histRandomHitDistribution->SetBinContent(4119,3);
histRandomHitDistribution->SetBinContent(4120,4);
histRandomHitDistribution->SetBinContent(4121,7);
histRandomHitDistribution->SetBinContent(4122,5);
histRandomHitDistribution->SetBinContent(4123,1);
histRandomHitDistribution->SetBinContent(4124,2);
histRandomHitDistribution->SetBinContent(4125,2);
histRandomHitDistribution->SetBinContent(4126,6);
histRandomHitDistribution->SetBinContent(4128,5);
histRandomHitDistribution->SetBinContent(4129,4);
histRandomHitDistribution->SetBinContent(4130,5);
histRandomHitDistribution->SetBinContent(4131,4);
histRandomHitDistribution->SetBinContent(4132,1);
histRandomHitDistribution->SetBinContent(4133,2);
histRandomHitDistribution->SetBinContent(4134,3);
histRandomHitDistribution->SetBinContent(4135,5);
histRandomHitDistribution->SetBinContent(4136,1);
histRandomHitDistribution->SetBinContent(4137,2);
histRandomHitDistribution->SetBinContent(4138,1);
histRandomHitDistribution->SetBinContent(4139,1);
histRandomHitDistribution->SetBinContent(4140,1);
histRandomHitDistribution->SetBinContent(4162,1);
histRandomHitDistribution->SetBinContent(4183,1);
histRandomHitDistribution->SetBinContent(4184,15);
histRandomHitDistribution->SetBinContent(4185,77);
histRandomHitDistribution->SetBinContent(4186,160);
histRandomHitDistribution->SetBinContent(4187,220);
histRandomHitDistribution->SetBinContent(4188,284);
histRandomHitDistribution->SetBinContent(4189,217);
histRandomHitDistribution->SetBinContent(4190,115);
histRandomHitDistribution->SetBinContent(4191,64);
histRandomHitDistribution->SetBinContent(4192,68);
histRandomHitDistribution->SetBinContent(4193,45);
histRandomHitDistribution->SetBinContent(4194,29);
histRandomHitDistribution->SetBinContent(4195,11);
histRandomHitDistribution->SetBinContent(4196,8);
histRandomHitDistribution->SetBinContent(4197,5);
histRandomHitDistribution->SetBinContent(4198,2);
histRandomHitDistribution->SetBinContent(4199,2);
histRandomHitDistribution->SetBinContent(4200,1);
histRandomHitDistribution->SetBinContent(4202,1);
histRandomHitDistribution->SetBinContent(4203,3);
histRandomHitDistribution->SetBinContent(4205,1);
histRandomHitDistribution->SetBinContent(4206,2);
histRandomHitDistribution->SetBinContent(4207,4);
histRandomHitDistribution->SetBinContent(4208,6);
histRandomHitDistribution->SetBinContent(4209,2);
histRandomHitDistribution->SetBinContent(4210,5);
histRandomHitDistribution->SetBinContent(4211,3);
histRandomHitDistribution->SetBinContent(4212,5);
histRandomHitDistribution->SetBinContent(4213,5);
histRandomHitDistribution->SetBinContent(4214,2);
histRandomHitDistribution->SetBinContent(4215,6);
histRandomHitDistribution->SetBinContent(4216,2);
histRandomHitDistribution->SetBinContent(4217,6);
histRandomHitDistribution->SetBinContent(4218,5);
histRandomHitDistribution->SetBinContent(4219,1);
histRandomHitDistribution->SetBinContent(4220,2);
histRandomHitDistribution->SetBinContent(4221,3);
histRandomHitDistribution->SetBinContent(4222,2);
histRandomHitDistribution->SetBinContent(4223,2);
histRandomHitDistribution->SetBinContent(4224,3);
histRandomHitDistribution->SetBinContent(4225,3);
histRandomHitDistribution->SetBinContent(4226,1);
histRandomHitDistribution->SetBinContent(4227,1);
histRandomHitDistribution->SetBinContent(4228,2);
histRandomHitDistribution->SetBinContent(4229,2);
histRandomHitDistribution->SetBinContent(4230,2);
histRandomHitDistribution->SetBinContent(4232,2);
histRandomHitDistribution->SetBinContent(4233,1);
histRandomHitDistribution->SetBinContent(4234,1);
histRandomHitDistribution->SetBinContent(4235,1);
histRandomHitDistribution->SetBinContent(4236,1);
histRandomHitDistribution->SetBinContent(4237,1);
histRandomHitDistribution->SetBinContent(4238,1);
histRandomHitDistribution->SetBinContent(4242,1);
histRandomHitDistribution->SetBinContent(4245,1);
histRandomHitDistribution->SetBinContent(4251,1);
histRandomHitDistribution->SetBinContent(4256,1);
histRandomHitDistribution->SetBinContent(4286,6);
histRandomHitDistribution->SetBinContent(4287,63);
histRandomHitDistribution->SetBinContent(4288,168);
histRandomHitDistribution->SetBinContent(4289,216);
histRandomHitDistribution->SetBinContent(4290,227);
histRandomHitDistribution->SetBinContent(4291,161);
histRandomHitDistribution->SetBinContent(4292,126);
histRandomHitDistribution->SetBinContent(4293,86);
histRandomHitDistribution->SetBinContent(4294,45);
histRandomHitDistribution->SetBinContent(4295,31);
histRandomHitDistribution->SetBinContent(4296,16);
histRandomHitDistribution->SetBinContent(4297,9);
histRandomHitDistribution->SetBinContent(4298,1);
histRandomHitDistribution->SetBinContent(4300,2);
histRandomHitDistribution->SetBinContent(4301,2);
histRandomHitDistribution->SetBinContent(4302,2);
histRandomHitDistribution->SetBinContent(4306,1);
histRandomHitDistribution->SetBinContent(4307,1);
histRandomHitDistribution->SetBinContent(4308,2);
histRandomHitDistribution->SetBinContent(4310,2);
histRandomHitDistribution->SetBinContent(4313,1);
histRandomHitDistribution->SetBinContent(4314,2);
histRandomHitDistribution->SetBinContent(4315,2);
histRandomHitDistribution->SetBinContent(4316,2);
histRandomHitDistribution->SetBinContent(4317,3);
histRandomHitDistribution->SetBinContent(4318,3);
histRandomHitDistribution->SetBinContent(4319,1);
histRandomHitDistribution->SetBinContent(4320,2);
histRandomHitDistribution->SetBinContent(4321,2);
histRandomHitDistribution->SetBinContent(4324,6);
histRandomHitDistribution->SetBinContent(4325,1);
histRandomHitDistribution->SetBinContent(4326,3);
histRandomHitDistribution->SetBinContent(4327,4);
histRandomHitDistribution->SetBinContent(4329,4);
histRandomHitDistribution->SetBinContent(4331,1);
histRandomHitDistribution->SetBinContent(4332,2);
histRandomHitDistribution->SetBinContent(4337,1);
histRandomHitDistribution->SetBinContent(4338,1);
histRandomHitDistribution->SetBinContent(4339,2);
histRandomHitDistribution->SetBinContent(4388,17);
histRandomHitDistribution->SetBinContent(4389,56);
histRandomHitDistribution->SetBinContent(4390,167);
histRandomHitDistribution->SetBinContent(4391,180);
histRandomHitDistribution->SetBinContent(4392,196);
histRandomHitDistribution->SetBinContent(4393,172);
histRandomHitDistribution->SetBinContent(4394,101);
histRandomHitDistribution->SetBinContent(4395,94);
histRandomHitDistribution->SetBinContent(4396,65);
histRandomHitDistribution->SetBinContent(4397,23);
histRandomHitDistribution->SetBinContent(4398,6);
histRandomHitDistribution->SetBinContent(4399,9);
histRandomHitDistribution->SetBinContent(4400,4);
histRandomHitDistribution->SetBinContent(4401,5);
histRandomHitDistribution->SetBinContent(4402,2);
histRandomHitDistribution->SetBinContent(4404,1);
histRandomHitDistribution->SetBinContent(4405,1);
histRandomHitDistribution->SetBinContent(4409,1);
histRandomHitDistribution->SetBinContent(4410,3);
histRandomHitDistribution->SetBinContent(4411,2);
histRandomHitDistribution->SetBinContent(4412,2);
histRandomHitDistribution->SetBinContent(4413,1);
histRandomHitDistribution->SetBinContent(4414,1);
histRandomHitDistribution->SetBinContent(4415,3);
histRandomHitDistribution->SetBinContent(4416,1);
histRandomHitDistribution->SetBinContent(4417,1);
histRandomHitDistribution->SetBinContent(4418,3);
histRandomHitDistribution->SetBinContent(4419,7);
histRandomHitDistribution->SetBinContent(4420,2);
histRandomHitDistribution->SetBinContent(4421,2);
histRandomHitDistribution->SetBinContent(4423,3);
histRandomHitDistribution->SetBinContent(4424,3);
histRandomHitDistribution->SetBinContent(4425,2);
histRandomHitDistribution->SetBinContent(4426,3);
histRandomHitDistribution->SetBinContent(4429,2);
histRandomHitDistribution->SetBinContent(4430,1);
histRandomHitDistribution->SetBinContent(4431,1);
histRandomHitDistribution->SetBinContent(4432,2);
histRandomHitDistribution->SetBinContent(4490,10);
histRandomHitDistribution->SetBinContent(4491,53);
histRandomHitDistribution->SetBinContent(4492,149);
histRandomHitDistribution->SetBinContent(4493,188);
histRandomHitDistribution->SetBinContent(4494,231);
histRandomHitDistribution->SetBinContent(4495,128);
histRandomHitDistribution->SetBinContent(4496,93);
histRandomHitDistribution->SetBinContent(4497,62);
histRandomHitDistribution->SetBinContent(4498,36);
histRandomHitDistribution->SetBinContent(4499,23);
histRandomHitDistribution->SetBinContent(4500,11);
histRandomHitDistribution->SetBinContent(4501,7);
histRandomHitDistribution->SetBinContent(4502,5);
histRandomHitDistribution->SetBinContent(4508,1);
histRandomHitDistribution->SetBinContent(4510,1);
histRandomHitDistribution->SetBinContent(4511,2);
histRandomHitDistribution->SetBinContent(4512,1);
histRandomHitDistribution->SetBinContent(4515,2);
histRandomHitDistribution->SetBinContent(4516,1);
histRandomHitDistribution->SetBinContent(4517,1);
histRandomHitDistribution->SetBinContent(4518,3);
histRandomHitDistribution->SetBinContent(4519,1);
histRandomHitDistribution->SetBinContent(4520,1);
histRandomHitDistribution->SetBinContent(4521,1);
histRandomHitDistribution->SetBinContent(4522,1);
histRandomHitDistribution->SetBinContent(4524,1);
histRandomHitDistribution->SetBinContent(4525,3);
histRandomHitDistribution->SetBinContent(4527,1);
histRandomHitDistribution->SetBinContent(4529,1);
histRandomHitDistribution->SetBinContent(4592,18);
histRandomHitDistribution->SetBinContent(4593,69);
histRandomHitDistribution->SetBinContent(4594,129);
histRandomHitDistribution->SetBinContent(4595,188);
histRandomHitDistribution->SetBinContent(4596,167);
histRandomHitDistribution->SetBinContent(4597,128);
histRandomHitDistribution->SetBinContent(4598,93);
histRandomHitDistribution->SetBinContent(4599,37);
histRandomHitDistribution->SetBinContent(4600,20);
histRandomHitDistribution->SetBinContent(4601,18);
histRandomHitDistribution->SetBinContent(4602,11);
histRandomHitDistribution->SetBinContent(4603,7);
histRandomHitDistribution->SetBinContent(4604,1);
histRandomHitDistribution->SetBinContent(4607,1);
histRandomHitDistribution->SetBinContent(4620,1);
histRandomHitDistribution->SetBinContent(4627,1);
histRandomHitDistribution->SetBinContent(4628,1);
histRandomHitDistribution->SetBinContent(4631,1);
histRandomHitDistribution->SetBinContent(4632,1);
histRandomHitDistribution->SetBinContent(4634,1);
histRandomHitDistribution->SetBinContent(4636,1);
histRandomHitDistribution->SetBinContent(4694,21);
histRandomHitDistribution->SetBinContent(4695,65);
histRandomHitDistribution->SetBinContent(4696,120);
histRandomHitDistribution->SetBinContent(4697,153);
histRandomHitDistribution->SetBinContent(4698,90);
histRandomHitDistribution->SetBinContent(4699,60);
histRandomHitDistribution->SetBinContent(4700,40);
histRandomHitDistribution->SetBinContent(4701,23);
histRandomHitDistribution->SetBinContent(4702,15);
histRandomHitDistribution->SetBinContent(4703,10);
histRandomHitDistribution->SetBinContent(4704,1);
histRandomHitDistribution->SetBinContent(4705,2);
histRandomHitDistribution->SetBinContent(4706,1);
histRandomHitDistribution->SetBinContent(4716,1);
histRandomHitDistribution->SetBinContent(4726,1);
histRandomHitDistribution->SetBinContent(4796,18);
histRandomHitDistribution->SetBinContent(4797,74);
histRandomHitDistribution->SetBinContent(4798,126);
histRandomHitDistribution->SetBinContent(4799,120);
histRandomHitDistribution->SetBinContent(4800,96);
histRandomHitDistribution->SetBinContent(4801,53);
histRandomHitDistribution->SetBinContent(4802,37);
histRandomHitDistribution->SetBinContent(4803,25);
histRandomHitDistribution->SetBinContent(4804,16);
histRandomHitDistribution->SetBinContent(4805,6);
histRandomHitDistribution->SetBinContent(4806,5);
histRandomHitDistribution->SetBinContent(4807,3);
histRandomHitDistribution->SetBinContent(4808,2);
histRandomHitDistribution->SetBinContent(4809,1);
histRandomHitDistribution->SetBinContent(4811,1);
histRandomHitDistribution->SetBinContent(4897,1);
histRandomHitDistribution->SetBinContent(4898,39);
histRandomHitDistribution->SetBinContent(4899,59);
histRandomHitDistribution->SetBinContent(4900,93);
histRandomHitDistribution->SetBinContent(4901,154);
histRandomHitDistribution->SetBinContent(4902,93);
histRandomHitDistribution->SetBinContent(4903,65);
histRandomHitDistribution->SetBinContent(4904,33);
histRandomHitDistribution->SetBinContent(4905,16);
histRandomHitDistribution->SetBinContent(4906,12);
histRandomHitDistribution->SetBinContent(4907,3);
histRandomHitDistribution->SetBinContent(4908,5);
histRandomHitDistribution->SetBinContent(4909,3);
histRandomHitDistribution->SetBinContent(4910,3);
histRandomHitDistribution->SetBinContent(4911,2);
histRandomHitDistribution->SetBinContent(4999,1);
histRandomHitDistribution->SetBinContent(5000,23);
histRandomHitDistribution->SetBinContent(5001,41);
histRandomHitDistribution->SetBinContent(5002,86);
histRandomHitDistribution->SetBinContent(5003,94);
histRandomHitDistribution->SetBinContent(5004,50);
histRandomHitDistribution->SetBinContent(5005,37);
histRandomHitDistribution->SetBinContent(5006,25);
histRandomHitDistribution->SetBinContent(5007,12);
histRandomHitDistribution->SetBinContent(5008,10);
histRandomHitDistribution->SetBinContent(5009,1);
histRandomHitDistribution->SetBinContent(5102,10);
histRandomHitDistribution->SetBinContent(5103,42);
histRandomHitDistribution->SetBinContent(5104,55);
histRandomHitDistribution->SetBinContent(5105,55);
histRandomHitDistribution->SetBinContent(5106,47);
histRandomHitDistribution->SetBinContent(5107,18);
histRandomHitDistribution->SetBinContent(5108,7);
histRandomHitDistribution->SetBinContent(5203,16);
histRandomHitDistribution->SetBinContent(5204,38);
histRandomHitDistribution->SetBinContent(5205,37);
histRandomHitDistribution->SetBinContent(5206,47);
histRandomHitDistribution->SetBinContent(5207,34);
histRandomHitDistribution->SetBinContent(5208,20);
histRandomHitDistribution->SetBinContent(5209,3);
histRandomHitDistribution->SetBinContent(5210,1);

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
            Int_t nHitsOnLadder = mcPxlLadderHitCol->numberOfHits()*4/100+1;
            Float_t randomRatio = histRandomHitDistribution->ProjectionX("_px",nHitsOnLadder)->GetRandom();
            cout << nHitsOnLadder << " " << randomRatio << endl;
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

                    //Int_t clusterSize = 4;
                    //while (clusterSize > 9) clusterSize = (Int_t)dataH2->GetRandom();
                    Int_t clusterSize = (Int_t)dataH2->GetRandom();
                    //cout << "cluster size: " << clusterSize << endl;
                    // init rowColumn[25][2];
                    for (int i=0;i<25;i++) for (int j=0;j<2;j++) rowColumn[i][j] = 0;
                    for (int i = 0 ; i < clusterSize ; i++) pxlRawHitCol.addRawHit(makeRawHit(localPixHitPos[0],localPixHitPos[2], iSec + 1, iLad + 1, iSen + 1, idTruth, i));

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
    return mFirstPixelX - StPxlConsts::kPixelSize * value - StPxlConsts::kPixelSize/2;
}
Float_t StPxlFastSim::getLocalZ(Int_t value){
    const double mFirstPixelZ = -(StPxlConsts::kPxlNumColumnsPerSensor - 1) * StPxlConsts::kPixelSize / 2;
    return mFirstPixelZ + StPxlConsts::kPixelSize * value + StPxlConsts::kPixelSize/2;
}
StPxlRawHit StPxlFastSim::makeRawHit(float localX, float localZ, int iSec, int iLad, int iSen, int idTruth, int nHits){
    Int_t row = getRow(localX);
    Int_t column = getColumn(localZ);
    Float_t smallR = 999999;
    for (int iRow=-2; iRow<3; iRow++) {
        if (row+iRow < 0) continue;
        if (row+iRow > StPxlConsts::kPxlNumRowsPerSensor - 1) continue;
        for (int iColumn=-2; iColumn<3; iColumn++) {
            if (column+iColumn < 0) continue;
            if (column+iColumn > StPxlConsts::kPxlNumColumnsPerSensor - 1) continue;
            
            Int_t beep = 0;
            //for (int i=0; i<nHits; i++) if (iRow == rowColumn[0][i] && iColumn==rowColumn[0][i]) beep=1;
            for (int i=0; i<nHits; i++) if (iRow == rowColumn[i][0] && iColumn==rowColumn[i][1]) beep=1;
            if (beep==1) {
                //cout << "skip!!" << endl;
                continue;
            }

            Float_t r = (getLocalX(row+iRow)-localX)*(getLocalX(row+iRow)-localX) + (getLocalZ(column+iColumn)-localZ)*(getLocalZ(column+iColumn)-localZ);
            if (r < smallR) {
                smallR = r;
                rowColumn[nHits][0]=iRow;
                rowColumn[nHits][1]=iColumn;
                //cout << iRow << " " << iColumn << " " << r << " " << smallR << " <-" << endl;
            }
            //else cout << iRow << " " << iColumn << " " << r << " " << smallR << endl;

        }
    }
    //cout << nHits << " " << rowColumn[nHits][0] << " " << rowColumn[nHits][1] << endl;
    tempHit.setSector(iSec);
    tempHit.setLadder(iLad);
    tempHit.setSensor(iSen);
    if (mPxlWrongRow) {
        //cout << "StPxlFastSim::mPxlWrongRow is on!" << endl;
        if (mRandom->flat(0.,1.) < (iLad == 1 ? mPxlWrongRowRatio1 : mPxlWrongRowRatio2)) tempHit.setRow(mRandom->flat(0,928));
        else tempHit.setRow(row+rowColumn[nHits][0]);
    }
    else tempHit.setRow(row+rowColumn[nHits][0]);
    tempHit.setColumn(column+rowColumn[nHits][1]);
    tempHit.setIdTruth(idTruth);
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

