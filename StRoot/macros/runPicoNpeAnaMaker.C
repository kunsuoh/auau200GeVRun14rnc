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

    // add your cuts here.

    // tracking
    npeCuts->setCutNHitsFitMax(20);
    
    // Electron pair cuts
    float dcaDaughtersMax = 0.5;  // maximum
    float minMass         = 0;
    float maxMass         = 0.1;
    npeCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);


    npeChain->Init();
    int nEntries = picoNpeAnaMaker->getEntries();
    cout << " Total entries = " << nEntries << endl;
    for (int iEvent = 0; iEvent < nEntries; ++iEvent)
    {
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
