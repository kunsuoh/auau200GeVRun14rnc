
void drawHistogram2(){

    TFile * infile = new TFile("out_18.root");
    TCanvas * cc = new TCanvas("cc","cc",500,500);
    TCanvas * cc2 = new TCanvas("cc2","cc2",500,1000);
    TCanvas * cc3 = new TCanvas("cc3","cc3",500,500);
    TCanvas * cc4 = new TCanvas("cc4","cc4",500,500);
    TCanvas * cc5 = new TCanvas("cc5","cc5",500,500);
    TCanvas * cc6 = new TCanvas("cc6","cc6",500,1000);
    
    cc2->Divide(1,2);

    cc->SetLogy();
    
    TString pid[6] = {"TPC","TPC+TOF","TPC+BEMC","TPC+BEMC+BSMD (BHT)","TPC+BEMC+BSMD (MB)","TPC+BEMC (BHT)"};
    double pt[6] = {1.5, 1.8, 2.5, 4.0, 6.5, 10};
    double hadronNSigEShape[7][6];
    // PID #2
    //float pidCutLw[6] = {0, -1.2, -1.2, -1.0, -1.0,
    //float pidCutHi[6] = {0, 1.8, 2.5, 3.0, 3.0, 0};
    // PID # 3
    float pidCutLw[6] = {0, -1.5, -1.4, -1.5, -1.1, 0};
    float pidCutHi[6] = {0, 1.8, 2.5, 3.0, 3.0, 0};

    const int nPid = 6;

    TH1F * hYield[nPid];
    TH1F * hRatio[nPid];
    TH1F * hMean[nPid];
    TH1F * hSigma[nPid];
    TH1F * hRawYield[nPid];
    TH1F * hUS;
    TH1F * hLS;
    TH1F * hSignal;

    TH1F * hIncE;
    
    TH1F * hNSigEPion;
    TH1F * hNSigEKaon;
    TH1F * hNSigEProton;
    TF1 * constant = new TF1("constant" ,"pol0", -0.1, 0.1);
    constant->SetParameter(0,1);

    cc6->cd();
    cc6->Divide(1,2);
    for (int iPid=2; iPid<6; iPid++) for (int iPt=1; iPt<5; iPt++) {
        cc6->cd(1)->SetLogy();
        TH1F * dum = (TH1F*)infile->Get(Form("histo_%d_%d_2_1",iPt, iPid));
        TH1F * his = (TH1F*)infile->Get(Form("histo_%d_%d_2_2",iPt, iPid));
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_1",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_1",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,-0.1,0.1);
        dum->Sumw2();
        his->Sumw2();
        hisLS->Sumw2();
        hisUS->Sumw2();
        hisPE->Sumw2();
        
        hisPE->Add(hisUS,hisLS,1,-1);
        double glo = 0.1;
     //   hisPE->Divide(constant,hisPE->GetMaximum()*glo);
        dum->Divide(constant,dum->GetMaximum()*glo);
        his->Divide(constant,his->GetMaximum());
        
        his->SetMarkerStyle(20);
        his->SetMarkerColor(1);
        his->SetMarkerSize(0.5);
        his->SetMaximum(2);
        his->SetMinimum(1e-4);
        
        dum->SetMarkerStyle(20);
        dum->SetMarkerColor(2);
        dum->SetMarkerSize(0.5);
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(0.5);

        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(0.5);
        hisUS->SetMinimum(0.5);

        
        hisPE->SetMarkerStyle(20);
        hisPE->SetMarkerColor(1);
        hisPE->SetMarkerSize(0.5);


     //   his->Draw("p");
     //   dum->Draw("psame");
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("psame");

        cc6->cd(2)->SetLogy();
   //     infile->Get(Form("histo_%d_%d_1_3",iPt, iPid))->Draw();
   //     infile->Get(Form("histo_%d_%d_0_3",iPt, iPid))->Draw("same");
        cc6->SaveAs(Form("~/Desktop/DcaAterPid_Pid%d_Pt%d.pdf",iPid,iPt));
 
        TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
        mc->Add(hisPE);
        mc->Add(dum);
        TF1 * fHF = new TF1("fHF",funHF,-0.1,0.1,3);
        fHF->SetParameters(1,20,0);
        TH1F * hHF = new TH1F("hHF","hHF",100,-0.1,0.1);
        hHF->FillRandom("fHF",100000);
        hHF->SetMarkerStyle(20);
        hHF->SetMarkerColor(5);
        hHF->SetMarkerSize(0.5);
     //   hHF->Divide(constant,hHF->GetMaximum());

     //   mc->Add(hHF);

        TFractionFitter* fit = new TFractionFitter(his, mc); // initialise
    //    fit->Constrain(1,0.0,10.0);               // constrain fraction 1 to be between 0 and 1
    //    fit->Constrain(2,0.0,10.0);               // constrain fraction 1 to be between 0 and 1
        fit->SetRangeX(1,100);                    // use only the first 15 bins in the fit
        Int_t status = fit->Fit();               // perform the fit
        std::cout << "fit status: " << status << std::endl;
        if (status == 0) {                       // check on fit status
            TH1F* result = (TH1F*) fit->GetPlot();
            his->Draw("p");
            result->SetMarkerStyle(24);
            result->SetMarkerColor(1);
            result->SetMarkerSize(1);
            double aPE, aHD, aHF, ePE, eHD, eHF;
            fit->GetResult(0,aPE,ePE);
            fit->GetResult(1,aHD,eHD);
          //  fit->GetResult(2,aHF,eHF);
            
            hisPE->Multiply(constant,aPE*glo);
            dum->Multiply(constant,aHD*glo);
            hHF->Multiply(constant,aHF);

            result->Draw("samep");
            hHF->Draw("psame");
            dum->Draw("psame");
            hisPE->Draw("psame");
            cc6->SetLogy();
            cc6->SaveAs(Form("~/Desktop/DcaAterPid_Fit_Pid%d_Pt%d.pdf",iPid,iPt));
        }

    }
    
    
    cout << "=========>START iPt loop ! " << endl;
    for (int iPt=1; iPt<6; iPt++){
        cout << "=========>=========>iPt : " << iPt << endl;
        hNSigEPion = (TH1F*)infile->Get(Form("histoNSigE_%d_0",iPt));
        hNSigEKaon = (TH1F*)infile->Get(Form("histoNSigE_%d_1",iPt));
        hNSigEProton = (TH1F*)infile->Get(Form("histoNSigE_%d_2",iPt));
        TF1 * fitfunHadron = new TF1("fitfunHadron","gaus",-13,13);
        fitfunHadron->SetParameters(1000,-5,1);
        
        cc4->cd()->SetLogy();
        // pion
        hNSigEPion->Draw();
        hNSigEPion->Fit(fitfunHadron,"NOR+");
        fitfunHadron->SetLineColor(2);
        fitfunHadron->SetLineStyle(2);
        fitfunHadron->Draw("same");
        cc4->SaveAs(Form("~/Desktop/hadron_PionNSigE_pt%d.pdf",iPt));
        hadronNSigEShape[iPt][0] = fitfunHadron->GetParameter(1);
        hadronNSigEShape[iPt][1] = fitfunHadron->GetParameter(2);
        
        // kaon
        hNSigEKaon->Draw();
        hNSigEKaon->Fit(fitfunHadron,"NOR+");
        fitfunHadron->SetLineColor(2);
        fitfunHadron->SetLineStyle(2);
        fitfunHadron->Draw("same");
        cc4->SaveAs(Form("~/Desktop/hadron_KaonNSigE_pt%d.pdf",iPt));
        hadronNSigEShape[iPt][2] = fitfunHadron->GetParameter(1);
        hadronNSigEShape[iPt][3] = fitfunHadron->GetParameter(2);
        
        // proton
        hNSigEProton->Draw();
        hNSigEProton->Fit(fitfunHadron,"NOR+");
        fitfunHadron->SetLineColor(2);
        fitfunHadron->SetLineStyle(2);
        fitfunHadron->Draw("same");
        cc4->SaveAs(Form("~/Desktop/hadron_ProtonNSigE_pt%d.pdf",iPt));
        hadronNSigEShape[iPt][4] = fitfunHadron->GetParameter(1);
        hadronNSigEShape[iPt][5] = fitfunHadron->GetParameter(2);
        
    } // end iPt loop
    cout << "=========>END iPt loop !" << endl;

    cout << "=========>START iPt loop ! " << endl;
    for (int iPid=0; iPid<6; iPid++){
        cout << "=========>=========>iPid : " << iPid << endl;
        hRatio[iPid] = new TH1F(Form("hRatio_%d",iPid),Form("hRatio_%d",iPid),5,pt);
        hYield[iPid] = new TH1F(Form("hYield_%d",iPid),Form("hYield_%d",iPid),5,pt);
        hMean[iPid] = new TH1F(Form("hMean_%d",iPid),Form("hMean_%d",iPid),5,pt);
        hSigma[iPid] = new TH1F(Form("hSigma_%d",iPid),Form("hSigma_%d",iPid),5,pt);
        hRawYield[iPid] = new TH1F(Form("hRawYield_%d",iPid),Form("hRawYield_%d",iPid),5,pt);

        cout << "START iPt loop ! " << endl;
        for (int iPt=1; iPt<6; iPt++){
            cout << "=========>=========>=========>iPt : " << iPt << endl;
            double dpt = pt[iPt]-pt[iPt-1];
            
            hUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_0",iPt,iPid));hUS->Sumw2();
            hLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_0",iPt,iPid));hLS->Sumw2();
            hIncE = (TH1F*)infile->Get(Form("histo_%d_%d_2_0",iPt,iPid));hIncE->Sumw2();
            
            hSignal = (TH1F*)hUS->Clone();
            hSignal->Add(hLS,-1);
            
            TF1 * fitfun = new TF1("fitfun","gaus",-2,3);
            fitfun->SetParameters(hSignal->GetMaximum(),0,1);
            hSignal->Fit(fitfun,"NOR+");
            
            hUS->SetTitle(Form("%s, %.1f < pT < %.1f",pid[iPid].Data(),pt[iPt-1],pt[iPt]));
            hUS->SetMinimum(0.7);
            hUS->GetXaxis()->SetRangeUser(-5,5);
            hUS->SetLineColor(4);
            hLS->SetLineColor(2);
            hSignal->SetLineColor(1);
            hSignal->SetLineWidth(2);
            fitfun->SetRange(-13,13);
            fitfun->SetLineColor(2);
            fitfun->SetLineStyle(2);
            fitfun->SetLineWidth(2);
            
            cc->cd();
            hUS->Draw();
            hLS->Draw("same");
            hSignal->Draw("same");
            fitfun->Draw("same");
            cc->SaveAs(Form("~/Desktop/PhE_Pid%d_Pt%d.pdf",iPid, iPt));
            
            hUS->Delete();
            hLS->Delete();
            hSignal->Delete();
            
//            cc->cd();
            hYield[iPid]->SetBinContent(iPt,fitfun->Integral(-13,13)/dpt);
            hYield[iPid]->SetBinError(iPt,fitfun->IntegralError(-13,13)/dpt);
            
            hMean[iPid]->SetBinContent(iPt,fitfun->GetParameter(1));
            hMean[iPid]->SetBinError(iPt,fitfun->GetParError(1));
            hSigma[iPid]->SetBinContent(iPt,fitfun->GetParameter(2));
            hSigma[iPid]->SetBinError(iPt,fitfun->GetParError(2));
            
            cc3->cd()->SetLogy();
       //     hIncE->SetTitle(Form("%s, %.1f < pT < %.1f",pid[iPid].Data(),pt[iPt-1],pt[iPt]));
            hIncE->Draw();
            TF1 * fitfunPion = new TF1("fitfunPion","gaus",-10,-6);
            TF1 * fitfunKaon = new TF1("fitfunKaon","gaus",-6,-3);
            TF1 * fitfunMerged= new TF1("fitfunMerged","gaus",4.5,8);
            TF1 * fitfunE = new TF1("fitfunE","gaus",-1,2);
            TF1 * fitfunAll= new TF1("fitfunAll","gaus(0)+gaus(3)+gaus(6)+gaus(9)",-13,13);
            
            fitfunPion->SetParameters(20000,-7,1);
            fitfunKaon->SetParameters(20000,-3,1);
            fitfunMerged->SetParameters(10,5,1);
            fitfunE->SetParameters(1000,-0.297203,0.963029);

       //     fitfunE->FixParameter(1,fitfun->GetParameter(1));
       //     fitfunE->FixParameter(2,fitfun->GetParameter(2));
            
            hIncE->Fit(fitfunPion,"NOR+");
            hIncE->Fit(fitfunKaon,"NOR+");
            hIncE->Fit(fitfunMerged,"NOR+");
            hIncE->Fit(fitfunE,"NOR+");
        
            double par[12];
            fitfunPion->GetParameters(&par[0]);
            fitfunKaon->GetParameters(&par[3]);
            fitfunMerged->GetParameters(&par[6]);
            fitfunE->GetParameters(&par[9]);
        
            fitfunAll->SetParameters(par);
            
        //   fitfunAll->FixParameter(7,4.5);
        //    fitfunAll->FixParameter(8,1);
     //       fitfunAll->FixParameter(10,fitfunE->GetParameter(1));
     //       fitfunAll->FixParameter(11,fitfunE->GetParameter(2));
    //        fitfunAll->FixParameter(10, -0.303934);
    //        fitfunAll->FixParameter(11, 0.952401);
            
            if (iPt==1) {
                fitfunAll->FixParameter(1, hadronNSigEShape[iPt][0]);
                fitfunAll->FixParameter(2, hadronNSigEShape[iPt][1]);
                
                fitfunAll->FixParameter(4, hadronNSigEShape[iPt][2]);
                fitfunAll->FixParameter(5, hadronNSigEShape[iPt][3]);
                
                fitfunAll->FixParameter(7, hadronNSigEShape[iPt][4]);
                fitfunAll->FixParameter(8, hadronNSigEShape[iPt][5]);
            }
            
            
           if (iPt==4) {
                fitfunAll->FixParameter(6, 0);
                fitfunAll->SetParLimits(4,-3.5.,-2.5);
                fitfunAll->SetParLimits(1,-6,-5);
            }
            
            fitfunAll->SetParLimits(0,1,hIncE->GetMaximum());
            fitfunAll->SetParLimits(3,1,hIncE->GetMaximum());
       //     fitfunAll->SetParLimits(6,1,1e10);
            fitfunAll->SetParLimits(9,1,hIncE->GetMaximum());
            
            fitfunAll->SetParLimits(8,0.5,1.5);
            fitfunAll->SetParLimits(7,3,7);

            fitfunAll->SetParLimits(11,0.5,1.5);
            fitfunAll->SetParLimits(10,-0.5,0.5);

            hIncE->Fit(fitfunAll,"NOR+");
            
            fitfunAll->GetParameters(&par[0]);
            
            fitfunPion->SetParameters(&par[0]);
            fitfunKaon->SetParameters(&par[3]);
            fitfunMerged->SetParameters(&par[6]);
            fitfunE->SetParameters(&par[9]);
            
            fitfunPion->SetRange(-13,13);
            fitfunKaon->SetRange(-13,13);
            fitfunMerged->SetRange(-13,13);
            fitfunE->SetRange(-13,13);
            
            fitfunPion->SetLineColor(2);
            fitfunKaon->SetLineColor(2);
            fitfunMerged->SetLineColor(2);
            fitfunE->SetLineColor(4);
            fitfunAll->SetLineStyle(2);
            fitfunAll->SetLineColor(6);
            
            fitfunPion->Draw("same");
            fitfunKaon->Draw("same");
            fitfunMerged->Draw("same");
            fitfunE->Draw("same");
            
            fitfunAll->Draw("same");
            
            cc3->SaveAs(Form("~/Desktop/IncE_%d_%d.pdf",iPid,iPt));
            cout << iPid << " " << iPt  << " Purity (" << pidCutLw[iPt] << ", " << pidCutHi[iPt] << ") : " << fitfunE->Integral(pidCutLw[iPt],pidCutHi[iPt])/fitfunAll->Integral(pidCutLw[iPt],pidCutHi[iPt]) * 100 << "%" << endl;
            
            if (iPt < 5) hRawYield[iPid]->SetBinContent(iPt,fitfunE->Integral(-13,13)/dpt);
          //  hRawYield[iPid]->SetBinError(iPt,fitfunE->IntegralError(-13,13)/dpt);

        } // end iPt loop
        cout << "=========>=========>END iPt loop ! " << endl;
        cc2->cd(1)->SetLogy();
        hYield[iPid]->SetTitle(pid[iPid]);
        hYield[iPid]->SetMarkerStyle(20);
        hYield[iPid]->SetMarkerColor(2);
        hYield[iPid]->SetMinimum(hYield[iPid]->GetBinContent(1)*10e-7);
        hYield[iPid]->SetMaximum(hYield[iPid]->GetBinContent(1)*10);
        hYield[iPid]->Draw("p");
        
        
        cc2->cd(2);
        hMean[iPid]->SetMinimum(-2);
        hMean[iPid]->SetMaximum(3);
        hMean[iPid]->SetTitle("");
        hMean[iPid]->SetMarkerStyle(20);
        hSigma[iPid]->SetMarkerStyle(24);
        hMean[iPid]->Draw("p");
        hSigma[iPid]->Draw("psame");
        TF1 * fitfunMean = new TF1("fitfunMean","pol0",1.5,10);
        TF1 * fitfunSigma = new TF1("fitfunSigma","pol0",1.5,10);
        
        hMean[iPid]->Fit(fitfunMean,"NOR+");
        hSigma[iPid]->Fit(fitfunSigma,"NOR+");
        fitfunMean->SetLineStyle(2);
        fitfunSigma->SetLineStyle(2);
        fitfunMean->Draw("same");
        fitfunSigma->Draw("same");
        
        
        cc2->SaveAs(Form("~/Desktop/Yield_Pid%d.pdf",iPid));
        
        // Raw Yield for Inclusive Electrons
        cc5->cd();
        cc5->SetLogy();
        hRawYield[iPid]->Sumw2();
        hRawYield[iPid]->SetTitle(pid[iPid]);
        hRawYield[iPid]->SetMarkerStyle(20);
        hRawYield[iPid]->SetMinimum(hYield[iPid]->GetBinContent(1)*10e-7);
        hRawYield[iPid]->SetMaximum(hYield[iPid]->GetBinContent(1)*10e2);
        hRawYield[iPid]->Draw("p");
        hYield[iPid]->Draw("psame");

        cc5->SaveAs(Form("~/Desktop/RawYield_Pid%d.pdf",iPid));
        
    } // end iPid loop
    cout << "=========>END iPid loop ! " << endl;
    
    TCanvas * cc4 = new TCanvas("cc4","cc4",500,500);
    //cc4->SetLogy();
    int i = 0;
    hRatio[i]->Divide(hYield[3],hYield[2]);
    hRatio[i]->SetTitle("(BHT)TPC+BEMC+BSMD / (MB)TPC+BEMC");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    
    int i = 1;
    hRatio[i]->Divide(hYield[2],hYield[0]);
    hRatio[i]->SetTitle("TPC+BEMC / TPC only");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 2;
    hRatio[i]->Divide(hYield[1],hYield[0]);
    hRatio[i]->SetTitle("TPC+TOF / TPC only");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 3;
    hRatio[i]->Divide(hYield[3],hYield[5]);
    hRatio[i]->SetTitle("(BHT)TPC+BEMC+BSMD / (BHT)TPC+BEMC");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 4;
    hRatio[i]->Divide(hYield[4],hYield[2]);
    hRatio[i]->SetTitle("(MB)TPC+BEMC+BSMD / (MB)TPC+BEMC");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    
}

double funHF(double * x, double * par){
    if (x[0]>0) return TMath::Exp(-1*x[0]*par[1]+par[2])*par[0];
    else return TMath::Exp(x[0]*par[1]+par[2])*par[0];
}
