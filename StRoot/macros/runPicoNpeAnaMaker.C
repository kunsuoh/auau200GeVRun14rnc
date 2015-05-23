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


void runPicoNpeAnaMaker(TString Npelist, TString outFileName, TString badRunListFileName = "picoList_bad_MB.list")
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
    gSystem->Load("StPicoNpeEventMaker");
    gSystem->Load("StPicoNpeAnaMaker");
    
    chain = new StChain();
    
    // create list of picoDst files
    TString command = "sed 's/hft\\\/npeTree/picodsts/g' " + Npelist + " >correspondingPico.list";
    gSystem->Exec(command.Data());
    command = "sed -i 's/picoNpe/picoDst/g' correspondingPico.list";
    gSystem->Exec(command.Data());
    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, "correspondingPico.list", "picoDstMaker");
    StPicoNpeAnaMaker*  picoNpeAnaMaker = new StPicoNpeAnaMaker("picoNpeAnaMaker", Npelist, outFileName.Data(), picoDstMaker);
    
    cout << "DEBUG!" << endl
    // -------------- USER variables -------------------------
    
    // add your cuts here.
    
    chain->Init();
    cout << "DEBUG!" << endl
    int nEntries = picoNpeAnaMaker->getEntries();
    cout << "DEBUG!" << endl
    for (int iEvent = 0; iEvent < nEntries; ++iEvent)
    {
        chain->Clear();
        int iret = chain->Make();
        if (iret)
        {
            cout << "Bad return code!" << iret << endl;
            break;
        }
    }
    cout << "DEBUG!" << endl
    
    chain->Finish();
    delete chain;
    
    // delete list of picos
    command = "rm -f correspondingPico.list";
    gSystem->Exec(command.Data());
    
}