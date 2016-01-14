#ifndef AQHistograms_hh
#define AQHistograms_hh

TH2D *vertexall;   
TH1D *vertexZall;
TH1D *multiplicityP;
TH1D *refmultP;
  ////////
TH1D *pthist;
TH2D *etaphihist;
TH2D *dEdxphist;
TH2D *m2tofhist;
  /////////////////
TH1D *istposshist;
TH1D *istfithist;
TH1D *pxlposshist;
TH1D *pxlfithist;
  /////////////////
TH1D *istposscuthist;
TH1D *istfitcuthist;
TH1D *pxlposscuthist;
TH1D *pxlfitcuthist;
  /////////////////
TH2D *dcaxypthist;
TH2D *dcazpthist;
  //TH1D *dcaxypthist;
  //TH1D *dcazpthist;
 /////////////////
TH1D *disthist;
TH1D *projhist;
TH1D *masshist;
TH1D *masRhist;
TH1D *lifehist;

void DefineHisto(){
  
  vertexall = new TH2D("vertexall","All the vertices in the run",200,-1,1,200,-1,1);
  vertexall->SetTitle("All vertices");
  vertexall->GetXaxis()->SetTitle("X (cm)");
  vertexall->GetYaxis()->SetTitle("Y (cm)");
  vertexZall = new TH1D("vertexZall","All the vertices in the run",100,-10,10);
  vertexZall->SetTitle("All vertices");
  vertexZall->GetXaxis()->SetTitle("Z (cm)");
  multiplicityP = new TH1D("multiplicityP","Multiplicity of Primary vertices",50,0,500);
  refmultP = new TH1D("refmultP","Reference Multiplicity of Primary vertices",30,0,300);   //For Simulations

  ////////
  pthist = new TH1D("pthist","Pt all tracks",500,0,5);
  pthist->GetXaxis()->SetTitle("pt");
  etaphihist =  new TH2D("etaphihist","Eta vs phi",1000,-3.2,3.2,1000,-2,2);
  etaphihist->GetYaxis()->SetTitle("eta");
  etaphihist->GetXaxis()->SetTitle("phi");
  dEdxphist =  new TH2D("dEdxphist","Pid",1000,0,5,1000,0,20);
  dEdxphist->GetXaxis()->SetTitle("p");
  dEdxphist->GetYaxis()->SetTitle("dEdx (keV/cm)");
  m2tofhist = new TH2D("m2tofhist","Mass^2 from TOF",600,0,3,400,-0.5,1.5);
  m2tofhist->GetXaxis()->SetTitle("p");
  m2tofhist->GetYaxis()->SetTitle("Mass^{2}(GeV/c^{2})");
 
  /////////////////
  istposshist = new TH1D("istposshist","Possible Ist hits",7,0,7);
  istfithist = new TH1D("istfithist","Fitted Ist hits",5,0,5);
  pxlposshist = new TH1D("pxlposshist","Possible Pxl hits",7,0,7);
  pxlfithist = new TH1D("pxlfithist","Fitted Pxl hits",5,0,5); 

  /////////////////
  istposscuthist = new TH1D("istposscuthist","Possible Ist hits {pt>0.2 && abs(eta)<1}",7,0,7);
  istfitcuthist = new TH1D("istfitcuthist","Fitted Ist hits  {pt>0.2 && abs(eta)<1}",5,0,5);
  pxlposscuthist = new TH1D("pxlposscuthist","Possible Pxl hits  {pt>0.2 && abs(eta)<1}",7,0,7);
  pxlfitcuthist = new TH1D("pxlfitcuthist","Fitted Pxl hits  {pt>0.2 && abs(eta)<1}",5,0,5); 

  /////////////////
  dcaxypthist = new TH2D("dcaxypthist","DCAxy vs 1/p",100,0,5,500,-0.2,0.2);
  dcazpthist  = new TH2D("dcazpthist", "DCAz vs 1/p",100,0,5,500,-0.2,0.2);
  //dcaxypthist = new TH1D("pxlposssechist","Possible Secondary Pxl hits",7,0,7);
  //dcazpthist = new TH1D("pxlfitsechist","Fitted Secondary Pxl hits",5,0,5);

 /////////////////
  disthist = new TH1D("disthist","Distance from the Primary vertex",30,0,30);
  projhist = new TH1D("projhist","Cosine of the angle with the Primary vertex",201,-1,1);
  masshist = new TH1D("masshist","Mass of the KFParticle",200,0.45,0.55);
  masRhist = new TH1D("masRhist","Rotated mass of the KFParticle",200,0.45,0.55);
  lifehist = new TH1D("lifehist","Life time of the KFParticle",201,0,1);
}

void PlotHisto(){
  //DefineHisto();
  gStyle->SetPalette(1);
  
  TCanvas *c1 = new TCanvas("c1","Vertex",0,0,800,600);
  c1->cd();
  gStyle->SetOptStat(11);
  c1->Divide(2,2);
  c1->cd(1);
  //c1->SetLogz();
  vertexall->Draw("colz");
  c1->cd(2);
  vertexZall->Draw();
  c1->cd(3);
  multiplicityP->Draw();
  c1->cd(4);
  refmultP->Draw();
  c1->cd();
  c1->Draw();
  //c1->Print("PrimVertices.png");

  TCanvas *c2 = new TCanvas("c2","Tracks",0,0,800,600);
  c2->cd();
  gStyle->SetOptStat(11);
  c2->Divide(2,2);
  c2->cd(1);
  pthist->Draw();
  c2->cd(2);
  c2->cd(2)->SetLogz();
  etaphihist->Draw("colz");
  c2->cd(3);
  c2->cd(3)->SetLogz();	  
  dEdxphist->Draw("colz");
  c2->cd(4);
  //c2->cd(4)->SetLogz();
  m2tofhist->Draw("colz");
  c2->Draw();
  //c2->Print("TracksPrimary.png");

  TCanvas *c3 = new TCanvas("c3","HFT Primary hits",0,0,800,600);
  c3->cd();
  gStyle->SetOptStat(11);
  c3->Divide(2,2);
  c3->cd(1);
  c3->cd(1)->SetLogy();
  istposshist->Draw();
  c3->cd(2);
  c3->cd(2)->SetLogy();
  pxlposshist->Draw();
  c3->cd(3);
  c3->cd(3)->SetLogy();	  
  istfithist->Draw();
  c3->cd(4);
  c3->cd(4)->SetLogy();
  pxlfithist->Draw();
  c3->Draw();
  //c3->Print("HftHits.png");

  TCanvas *c4 = new TCanvas("c4","HFT Primary hits cut",0,0,800,600);
  c4->cd();
  gStyle->SetOptStat(11);
  c4->Divide(2,2);
  c4->cd(1);
  c4->cd(1)->SetLogy();
  istposscuthist->Draw();
  c4->cd(2);
  c4->cd(2)->SetLogy();
  pxlposscuthist->Draw();
  c4->cd(3);
  c4->cd(3)->SetLogy();	  
  istfitcuthist->Draw();
  c4->cd(4);
  c4->cd(4)->SetLogy();
  pxlfitcuthist->Draw();
  c4->Draw();
  //c4->Print("HftHitscut.png");

  TCanvas *c5 = new TCanvas("c5","DCA Primary",0,0,800,600);
  c5->cd();
  gStyle->SetOptStat(11);
  c5->Divide(2,2);
  c5->cd(1);
  dcaxypthist->Draw("colz");
  c5->cd(2);
  //c5->cd(2)->SetLogz();
  dcazpthist->Draw("colz");
  c5->cd(3);
  //c5->cd(3)->SetLogz();
  dcaxypthist->FitSlicesY();
  TH1D *dcaxypthist_2 = (TH1D*)gDirectory->Get("dcaxypthist_2");
  dcaxypthist_2->SetMarkerStyle(21);
  dcaxypthist_2->SetMarkerSize(0.8);
  dcaxypthist_2->SetXTitle("1/p");
  dcaxypthist_2->SetMinimum(0);
  dcaxypthist_2->SetMaximum(0.1);
  dcaxypthist_2->SetStats(0);
  dcaxypthist_2->Draw();
  c5->cd(4);
  //c3->cd(4)->SetLogz();
  dcazpthist->FitSlicesY();
  TH1D *dcazpthist_2 = (TH1D*)gDirectory->Get("dcazpthist_2");
  dcazpthist_2->SetMarkerStyle(21);
  dcazpthist_2->SetMarkerSize(0.8);
  dcazpthist_2->SetXTitle("1/p");
  dcazpthist_2->SetMinimum(0);
  dcazpthist_2->SetMaximum(0.1);
  dcazpthist_2->SetStats(0);
  dcazpthist_2->Draw();
  c5->Draw();
  //c5->Print("DCAsprim.png");

  TCanvas *c6 = new TCanvas("c6","Secondaries Vertices",0,0,800,600);
  c6->cd();
  gStyle->SetOptStat(11);
  c6->Divide(2,2);
  c6->cd(1);
  disthist->Draw();
  c6->cd(2);
  projhist->Draw();
  c6->cd(3);
  masshist->Draw();
  masRhist->SetLineColor(2);
  masRhist->Draw("same");
  c6->cd(4);
  lifehist->Draw();
  c6->Draw();
  
  //////<-


}
#endif
