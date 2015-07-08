void fitHistogram(int i = 2){
    
    TFile * infile = new TFile("outfile_2_1.root");
    TH1D * histoUS = (TH1D*)infile->Get(Form("cdpt%dUS/h%d3",i,i*2+6));
    TH1D * histoLS = (TH1D*)infile->Get(Form("cdpt%dLS/h%d3",i,i*2+7));
    TH1D * histo = histoUS->Clone();
    histoUS->Sumw2();
    histoLS->Sumw2();
    histo->Add(histoUS,histoLS,1,-1.11);
    double par[9];
    TF1 * fit1 = new TF1("fit1","gaus",-13,-4);
    TF1 * fit2 = new TF1("fit2","gaus",-1,2);
    TF1 * fit3 = new TF1("fit3","gaus",3,7);
    TF1 * fitfun = new TF1("fitfun","gaus(0)+gaus(3)+gaus(6)",-13,13);
    
    
 //   fit1->SetParameters(30000, -5, 1);
    fit2->SetParameters(10000, 0, 1);
 //   fit3->SetParameters(1000, 5, 1);
    
 //   histo->Fit(fit1,"LNOR+");
    histo->Fit(fit2,"LNOR+");
 //   histo->Fit(fit3,"LNOR+");
    
 //   fit1->GetParameters(&par[0]);
    fit2->GetParameters(&par[3]);
 //   fit3->GetParameters(&par[6]);
    
//    fitfun->SetParameters(par);
//    histo->Fit(fitfun,"LNOR+");
//    fitfun->GetParameters(&par[0]);

    
    TCanvas * cc = new TCanvas("cc","cc",500,500);
    cc->SetLogy();
    
    histoUS->SetLineColor(4);
    histoLS->SetLineColor(2);
    histo->SetLineColor(1);
    
    histoUS->GetXaxis()->SetRangeUser(-5,5);
    histoUS->SetMinimum(0.7);
    histoUS->Draw();
    histoLS->Draw("same");
    histo->SetLineWidth(2);
    histo->Draw("same");
    fit1->SetParameters(par[0],par[1],par[2]);
    fit2->SetParameters(par[3],par[4],par[5]);
    fit3->SetParameters(par[6],par[7],par[8]);
    
    fit1->SetLineColor(4);
    fit1->SetLineStyle(1);
    fit1->SetLineWidth(0.5);
    fit1->SetRange(-13,13);
    fit1->Draw("same");
    
    fit2->SetLineColor(4);
    fit2->SetLineStyle(1);
    fit2->SetLineWidth(0.5);
    fit2->SetRange(-13,13);
    fit2->Draw("same");
    
    fit3->SetLineColor(4);
    fit3->SetLineStyle(1);
    fit3->SetLineWidth(0.5);
    fit3->SetRange(-13,13);
    fit3->Draw("same");
    
    fitfun->SetLineColor(2);
    fitfun->SetLineStyle(2);
    fitfun->SetRange(-13,13);
    fitfun->Draw("same");
    
    cout << par[3] << " " << par[4] << " " << par[5] << endl;
    

}
