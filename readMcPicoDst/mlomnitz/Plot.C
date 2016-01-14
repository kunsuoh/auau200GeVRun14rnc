#include <string>
#include "Riostream.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
#include "TDatime.h"
#include "TMath.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
//#include "AQHistKFPart.h"
//#include "TSystemDirectory.h"
//#include "TList.h"
//#include "TDirIter.h"
//#include "TSystemFile.h"

//--->START MAIN PROGRAM
//________________________________________________________________________________
void Plot(const Char_t *dirname= "out/", const Char_t *ext=".root") {

  //gROOT->cd();
  TDatime now;                                          //Set time in Root
  now.Print();

  //TFile *fout = new TFile("pout.root","recreate") ;
  //TNtuple *vertex;
    //////////////////////////////////////////////////////->Define the variables
  //DefineHisto();
  ///////////////////////////////////////////////////////<-
  Int_t size =30;
  Int_t bin = size*10;
  TH2D *his1 = new TH2D("his1","MCout",bin,-size,size,bin,-size,size);
  TH2D *hist = new TH2D("hist","Helix",bin,-size,size,bin,-size,size);
  TH2D *his2 = new TH2D("his2","KFPar",bin,-size,size,bin,-size,size);
  //TH1D *back = new TH1D("back","Reconstructed back",100,0.45,0.55);
  TCanvas *c1 = new TCanvas("c1","Decay lenght",0,0,1200,450);
  c1->cd();
  c1->Divide(3,1);
  // ----------------------------------------------
  Int_t count=0;
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
	//cout << fname.Data() << endl;
	//gStyle->SetPalette(1);
	Char_t eachfile[50];
	sprintf(eachfile,"%s%s",dirname,fname.Data());
	//cout << eachfile << endl;
	TFile f(eachfile);
	//vertex = (TNtuple*)f->Get("primaryvtx");
	TTree *verte1 = (TTree*)f.Get("MCTuple");
	TTree *vertex = (TTree*)f.Get("secHE");
	TTree *verte2 = (TTree*)f.Get("secKF");
	TCut cuth = "q==0 && proj>0.9 && Mass<0.1";
	TCut cutk = "q==0 && proj>0.9 && Mass<0.1";
	if (!verte1) {cout<< "No Vertex" << endl; verte1=0;}
	his1->SetDirectory(gDirectory); // "attach" it to the current root file
	hist->SetDirectory(gDirectory); // "attach" it to the current root file
	his2->SetDirectory(gDirectory); // "attach" it to the current root file
	c1->cd(1);
	if (verte1) verte1->Draw("MCY:MCX>>+his1","abs(MCY)<30 && abs(MCX)<30 && sqrt(MCY*MCY+MCX*MCX)>1","goff");
	c1->cd(1)->SetLogz();
	his1->Draw("colz");
	c1->cd(2);
	vertex->Draw("HEY:HEX>>+hist","abs(HEY)<30 && abs(HEX)<30 && sqrt(HEY*HEY+HEX*HEX)>1" && cuth,"goff");
	c1->cd(2)->SetLogz();
	hist->Draw("colz");
	c1->cd(3);
	verte2->Draw("KFY:KFX>>+his2","abs(KFY)<30 && abs(KFX)<30 && sqrt(KFY*KFY+KFX*KFX)>1" && cutk,"goff");
	c1->cd(3)->SetLogz();
	his2->Draw("colz");
	c1->Draw();
	his1->SetDirectory(gROOT); // "detach" it from the current root file
	hist->SetDirectory(gROOT); // "detach" it from the current root file
	his2->SetDirectory(gROOT); // "detach" it from the current root file
	//hist->Draw();
	//back->SetLineColor(3);
	//vertex->Draw("Mass>>+back","q!=0 && dist>1","same");
	count++;
	if (count%50==0){cout << "Files: " << count << endl; gPad->Update();}
      }
    }
  }
  
  //////<-
  //fout->Close();
  TDatime now1;
  now1.Print();
}
  
