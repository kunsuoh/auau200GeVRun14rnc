/* **************************************************
 *  A macro to run StPicoNpeAnaMaker
 *
 *  Authors:  **Kunsu OH (kunsu OH)
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;

StChain * npeChain;

void runPicoNpeAnaMaker(TString npeList, TString outFileName, TString badRunListFileName = "picoList_bad_MB.list")
{
    //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
    string SL_version = "SL15c";
    string env_SL = getenv("STAR");
    if (env_SL.find(SL_version) == string::npos)
    {
        cout << "Environment Star Library does not match the requested library in runPicoNpeEventMaker.C. Exiting..." << endl;
        exit(1);
    }
    
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    
    gSystem->Load("StBTofUtil");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StPicoPrescales");
    gSystem->Load("StPicoCutsBase");
    gSystem->Load("StPicoNpeEventMaker");
    gSystem->Load("StPicoNpeAnaMaker");

    npeChain = new StChain();
    
    // create list of picoDst files
    TString command = "sed 's/hft\\\/npeTree/picodsts/g' " + npeList + " >correspondingPico.list";
    gSystem->Exec(command.Data());
    command = "sed -i 's/picoNpe/picoDst/g' correspondingPico.list";
    gSystem->Exec(command.Data());
    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, "correspondingPico.list", "picoDstMaker");
    StPicoNpeAnaMaker*  picoNpeAnaMaker = new StPicoNpeAnaMaker("picoNpeAnaMaker", npeList, outFileName.Data(), picoDstMaker);
    
    StNpeCuts* npeCuts = new StNpeCuts("npeCuts");
    picoNpeAnaMaker->setNpeCuts(npeCuts);

    // -------------- USER variables -------------------------
    
    // -- File name of bad run list
    npeCuts->setBadRunListFileName(badRunListFileName);

    // Add your cuts here.
    
    // Event cuts
    npeCuts->setCutVzMax(6);
    npeCuts->setCutVzVpdVzMax(3);
    npeCuts->setCutTriggerWord(0xFFFFFFF);

    // Tagged electron cuts
    npeCuts->setCutElectronNHitsFitMax(15);
    npeCuts->setCutElectronNHitsdEdxMax(0);
    npeCuts->setCutPt(0.2, 20);
    npeCuts->setCutEta(-1., 1.);
    npeCuts->setCutDca(100);
    npeCuts->setCutElectronRequireHFT(true);
    npeCuts->setCutTPCNSigmaElectron(-2.0, 2.0);
    
    // Partner electron cuts
    npeCuts->setCutPartnerElectronNHitsFitMax(15);
    npeCuts->setCutPartnerPt(0.2, 20);
    npeCuts->setCutPartnerEta(-1., 1.);
    npeCuts->setCutPartnerTPCNSigmaElectron(-3.0, 3.0);
    npeCuts->setCutPartnerElectronRequireHFT(true);
    
    // Electron pair cuts
    float dcaDaughtersMax = 0.05;  // maximum
    float minMass         = 0;
    float maxMass         = 0.05;
    npeCuts->setCutElectronPair(dcaDaughtersMax, minMass, maxMass);

    
    // BEMC PID
    bool const bemc = true;
    float const minEoverP   = 0.8;
    float const maxEoverP   = 2.0;
    float const phiDist     = 0.015;
    float const zDist       = 3;
    float const assDist     = 0.06;
    npeCuts->setCutBemcPid(bemc, minEoverP, maxEoverP, phiDist, zDist, assDist);

    // BSMD PID
    bool const bsmd = false;
    int const nEta     = 1;
    int const nPhi     = 1;
    npeCuts->setCutBsmdPid(bsmd, nEta, nPhi);
    
    // TOF PID
    bool const tof = true;
    float const tofBeta     = 0.025;
    npeCuts->setCutTofPid(tof, tofBeta);
    
    npeChain->Init();
    int nEntries = picoNpeAnaMaker->getEntries();
    cout << " Total entries = " << nEntries << endl;
    for (int iEvent = 0; iEvent < nEntries; ++iEvent)
    {
        if (iEvent%100==0) cout << iEvent << endl;
        npeChain->Clear();
        int iret = npeChain->Make();
        if (iret)
        {
            cout << "Bad return code!" << iret << endl;
            break;
        }
    }
    
    npeChain->Finish();
    delete npeChain;
    
    // delete list of picos
    command = "rm -f correspondingPico.list";
    gSystem->Exec(command.Data());
    
}
