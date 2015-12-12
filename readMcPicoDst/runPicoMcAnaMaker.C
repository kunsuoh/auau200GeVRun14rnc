/* **************************************************
 *  A macro to run StPicoMcAnaMaker
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

void runPicoMcAnaMaker(TString mcPicoList="test.list", TString outFileName="test")
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
    
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StPicoPrescales");
    gSystem->Load("StPicoCutsBase");
    gSystem->Load("StPicoNpeEventMaker");
    gSystem->Load("StPicoMcAnaMaker");

    npeChain = new StChain();
    
    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, mcPicoList, "picoDstMaker");
    StPicoMcAnaMaker*  picoMcAnaMaker = new StPicoMcAnaMaker("picoMcAnaMaker", outFileName.Data(), picoDstMaker);
    
    StNpeCuts* npeCuts = new StNpeCuts("npeCuts");
    picoMcAnaMaker->setNpeCuts(npeCuts);

    // -------------- USER variables -------------------------

    // Event cuts
    npeCuts->setCutVzMax(6);
    npeCuts->setCutVzVpdVzMax(3);
    npeCuts->setCutTriggerWord(0xFFFFFFF);

    // Tagged electron cuts
    npeCuts->setCutElectronNHitsFitMax(15);
    npeCuts->setCutElectronNHitsdEdxMax(0);
    npeCuts->setCutPt(0.6, 20);
    npeCuts->setCutEta(-1., 1.);
    npeCuts->setCutDca(100);
    npeCuts->setCutElectronRequireHFT(true);
    
    // Partner electron cuts
    npeCuts->setCutPartnerElectronNHitsFitMax(15);
    npeCuts->setCutPartnerPt(0.6, 20);
    npeCuts->setCutPartnerEta(-1., 1.);
    npeCuts->setCutPartnerElectronRequireHFT(true);
    
    // Electron pair cuts
    float dcaDaughtersMax = 1.;  // maximum
    float minMass         = 0;
    float maxMass         = 0.4;
    npeCuts->setCutElectronPair(dcaDaughtersMax, minMass, maxMass);

    
    npeChain->Init();
    int nEntries = picoMcAnaMaker->getEntries();
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
