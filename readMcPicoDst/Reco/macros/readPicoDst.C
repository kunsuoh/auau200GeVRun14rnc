
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
void readPicoDst(const Char_t *inputFile="st_dezero_999999.picoDst.root", const Char_t *outputFile="test.root")
{
TStopwatch*   stopWatch = new TStopwatch();
 stopWatch->Start();

        Int_t nEvents = 10000000;
//	Int_t nEvents = 10;	
//Load all the System libraries
	
        gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StPicoDstMaker");
        gSystem->Load("StPicoMcAnaMaker");

	chain = new StChain();

	StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
	StPicoMcAnaMaker *picoMcAnaMaker = new StPicoMcAnaMaker(inputFile,"test",picoMaker);
        //StMyAnalysisMaker *anaMaker = new StMyAnalysisMaker("ana",picoMaker,outputFile);
	picoMaker->setMcMode(true);
	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
        if(nEvents>total) nEvents = total;
	for (Int_t i=0; i<nEvents; i++){
	  if(i%100==0)
		cout << "Working on eventNumber " << i << endl;
		
	  chain->Clear();
	  int iret = chain->Make(i);
		
	  if (iret) { cout << "Bad return code!" << iret << endl; break;}

	  total++;
		
	}
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;
	
	   stopWatch->Stop();   
	     stopWatch->Print();
	
}
