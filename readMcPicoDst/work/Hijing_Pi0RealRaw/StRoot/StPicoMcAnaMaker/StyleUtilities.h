#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
//
#include "TBox.h"
#include <fstream>
//
#include "TStyle.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TList.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
//
#include "TLine.h"
//Raa
//float const length = 0.027;
//float const width = 0.07;
//v2
float const length = 0.005;
float const width = 0.07;
void setGraphicsStyle()
{
    // **************************** Set graphic style ***************************************
    gStyle->SetEndErrorSize(5);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    //gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    //gStyle->SetPadRightMargin(0.115);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetTitleSize(0.056,"X");
    gStyle->SetTitleSize(0.056,"Y");
    gStyle->SetTitleOffset(2.0,"X");
    gStyle->SetTitleOffset(2.0,"Y");
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleFont(42, "X");
    gStyle->SetTitleFont(42, "Y");
    gStyle->SetTitleFont(42, "Z");
    gStyle->SetTitleOffset(1.0, "Y");
    //added
    //gStyle->SetLineWidth(2);
    gStyle->SetEndErrorSize(0.01);
    gStyle->SetFrameLineWidth(5);
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->UseCurrentStyle();
    // **************************************************************************************
}

void setStyle(TH1* h,int marker,Color_t color,std::string title="")
{
  h->GetYaxis()->SetTitle(title.c_str());
  h->SetMarkerStyle(marker);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(3);
  //v2
  h->SetMarkerSize(2);
  //Raa
  //h->SetMarkerSize(2.4);
  h->SetTitleFont(42, "X");
  h->SetTitleFont(42, "Y");
  h->SetLabelFont(42, "X");
  h->SetLabelFont(42, "Y");
  h->SetTitleSize(.056, "X");
  h->SetTitleSize(.056, "Y");
  h->SetLabelSize(.045, "X");
  h->SetLabelSize(.045, "Y");
  h->SetTitleOffset(1.1, "X");
  h->SetTitleOffset(1.0, "Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

void setStyle(TH2* h,int marker=0,Color_t color=0)
{
  if(marker) 
  {
    h->SetMarkerStyle(marker);
    //v2
    h->SetMarkerSize(2);
    //Raa
    //h->SetMarkerSize(2.4);
  }

  if(color)
  {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
  }

  h->SetTitleFont(42, "X");
  h->SetTitleFont(42, "Y");
  h->SetLabelFont(42, "X");
  h->SetLabelFont(42, "Y");
  h->SetTitleSize(.053, "X");
  h->SetTitleSize(.053, "Y");
  h->SetLabelSize(.053, "X");
  h->SetLabelSize(.053, "Y");
  h->SetTitleOffset(1., "X");
  h->SetTitleOffset(1.0, "Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

void setStyle(TGraph* gr,int marker,Color_t color)
{
  gr->SetMarkerStyle(marker);
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
  gr->SetLineWidth(1.4);
  //v2
  gr->SetMarkerSize(2);
  //RAA
  //gr->SetMarkerSize(2.4);
}

TLegend* legend(std::string header="",float y1Ndc=0.55,float y2Ndc=0.88,float x1Ndc=0.3,float x2Ndc=0.8)
{
  static int nLegends = 0;
  ++nLegends;
  TLegend* leg = new TLegend(x1Ndc,y1Ndc,x2Ndc,y2Ndc);
  leg->SetName(Form("LEGEND_%i",nLegends));
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  if(!header.empty()) leg->SetHeader(header.c_str());

  return leg;
}

TPaveText* paveText(float x1=0.2,float y1=0.6,float x2=0.3,float y2=0.85)
{
  TPaveText* pv = new TPaveText(x1,y1,x2,y2,"NDC");
  pv->SetFillColor(0);
  pv->SetTextColor(1);
  pv->SetTextFont(42);
  pv->SetTextSize(0.045);
  pv->SetBorderSize(0);

  return pv;
}

TCanvas* canvas(bool square=false)
{
  static int nCvs = 0;
  ++nCvs;
  TCanvas* cv = NULL;
  //if(!square) cv = new TCanvas(Form("cv%i",nCvs),Form("cv%i",nCvs),2*600,2*338);
  if(!square) cv = new TCanvas(Form("cv%i",nCvs),Form("cv%i",nCvs),2*400,2*300);
  //origin
  else cv = new TCanvas(Form("cv%i",nCvs),Form("cv%i",nCvs),1000,1000);


  return cv;
}

TPaveStats* statBox(TH1* h)
{
  gPad->Update(); 
  TPaveStats* st0 = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");

  st0->SetName(Form("st_%s",h->GetName()));
  st0->SetTextColor(h->GetMarkerColor());
  st0->SetTextSize(0.035);
  st0->SetX1NDC(0.55);
  st0->SetX2NDC(0.98);

  return st0;
}

void alignStatBox(TPaveStats* const st0,TPaveStats* const st1)
{
  st1->SetX1NDC(st0->GetX1NDC());
  st1->SetX2NDC(st0->GetX2NDC());
  st1->SetY2NDC(st0->GetY1NDC());
  st1->SetY1NDC(st0->GetY1NDC()-(st0->GetY2NDC()-st0->GetY1NDC()));
}
void setTitles(TGraph* gr, string title="", string xAxis = "", string yAxis = "")
{
  gr->SetTitle(title.c_str());
  gr->GetXaxis()->SetTitle(xAxis.c_str());
  gr->GetYaxis()->SetTitle(yAxis.c_str());
}
void setTitles(TH1* gr, string title="", string xAxis = "", string yAxis = "")
{
  gr->SetTitle(title.c_str());
  gr->GetXaxis()->SetTitle(xAxis.c_str());
  gr->GetYaxis()->SetTitle(yAxis.c_str());
}
void drawErrors(TGraphErrors* gr)
{
  Color_t col = gr->GetMarkerColor();

  for(int ii = 0; ii< gr->GetN(); ++ii){
    double x1 = gr->GetX()[ii]-width;
    double x2 = gr->GetX()[ii]+width;
    double y1 = gr->GetY()[ii]-gr->GetEY()[ii];
    double y2 = gr->GetY()[ii]+gr->GetEY()[ii];
    
    TLine *la = new TLine(x1, y1, x1, y1+length);
    la->SetLineColor(col);
    //la->Draw("same");
    if(y1>-0.03) la->Draw("same");
    
    TLine *lb = new TLine(x2, y1, x2, y1+length);
    lb->SetLineColor(col);
    //lb->Draw("same");
    if(y1>-0.03) lb->Draw("same");
    
    TLine *lc = new TLine(x1, y2, x1, y2-length);
    lc->SetLineColor(col);
    lc->Draw("same");

    TLine *ld = new TLine(x2, y2, x2, y2-length);
    ld->SetLineColor(col);
    ld->Draw("same");
    
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->SetLineColor(col);
    //le->Draw("same");
    if(y1>-0.03) le->Draw("same");

    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineColor(col);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }   
}
void drawErrors(TGraphAsymmErrors* gr)
{
  Color_t col = gr->GetMarkerColor();

  for(int ii = 0; ii< gr->GetN(); ++ii){
    double x1 = gr->GetX()[ii]-width;
    double x2 = gr->GetX()[ii]+width;
    double y1 = gr->GetY()[ii]-gr->GetEYlow()[ii];
    double y2 = gr->GetY()[ii]+gr->GetEYhigh()[ii];
    
    TLine *la = new TLine(x1, y1, x1, y1+length);
    la->SetLineColor(col);
    la->Draw("same");
    //if(y1>-0.03) la->Draw("same");
    
    TLine *lb = new TLine(x2, y1, x2, y1+length);
    lb->SetLineColor(col);
    lb->Draw("same");
    //if(y1>-0.03) lb->Draw("same");
    
    TLine *lc = new TLine(x1, y2, x1, y2-length);
    lc->SetLineColor(col);
    lc->Draw("same");

    TLine *ld = new TLine(x2, y2, x2, y2-length);
    ld->SetLineColor(col);
    ld->Draw("same");
    
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->SetLineColor(col);
    le->Draw("same");

    //    if(y1>-0.03) le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineColor(col);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }
}
