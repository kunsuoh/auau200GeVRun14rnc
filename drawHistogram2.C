
void drawHistogram2(int in = 33){
    
    TFile * infile = new TFile(Form("out_%d.root", in));
    TFile * outfile = new TFile(Form("outfile_%d.root",in),"RECREATE");
    TCanvas * cc = new TCanvas("cc","cc",500,500);
    TCanvas * cc2 = new TCanvas("cc2","cc2",500,1000);
    TCanvas * cc3 = new TCanvas("cc3","cc3",500,500);
    TCanvas * cc4 = new TCanvas("cc4","cc4",500,500);
    TCanvas * cc5 = new TCanvas("cc5","cc5",500,1000);
    TCanvas * cc6 = new TCanvas("cc6","cc6",500,1000);
    TCanvas * cc7 = new TCanvas("cc7","cc7",500,1000);
    TCanvas * cc8 = new TCanvas("cc8","cc8",1000,1000);
    TCanvas * cc9 = new TCanvas("cc9","cc9",1000,500);
    TCanvas * cc10 = new TCanvas("cc10","cc10",500,500);
    TCanvas * cc11 = new TCanvas("cc11","cc11",500,500);
    
    cc2->Divide(1,2);
    cc5->Divide(1,2);
    cc6->Divide(1,2);
    cc7->Divide(1,2);
    cc8->Divide(4,4);
    cc9->Divide(4,2);

    cc->SetLogy();
    
    const int nPid = 8; // 0-7
    const int nPt = 6;  // 1-5

    TString pid[nPid] = {"TPC (MB)","TPC+TOF (MB)","TPC+BEMC (MB)","TPC+BEMC+BSMD (MB)","TPC (BHT)","TPC+TOF (BHT)","TPC+BEMC (BHT)","TPC+BEMC+BSMD (BHT)"};
    double pt[nPt] = {1.5, 1.8, 2.5, 4.0, 6.5, 10};
    double hadronNSigEShape[7][6];
    // PID #2
    //float pidCutLw[6] = {0, -1.2, -1.2, -1.0, -1.0,
    //float pidCutHi[6] = {0, 1.8, 2.5, 3.0, 3.0, 0};
    // PID # 3
    float pidCutLw[nPt] = {0, -1.5, -1.4, -1.5, -1.1, 0};
    float pidCutHi[nPt] = {0, 1.8, 2.5, 3.0, 3.0, 0};


    TH1F * hYield[nPid];
    TH1F * hRatio[nPid];
    TH1F * hMean[nPid];
    TH1F * hSigma[nPid];
    TH1F * hRawYield[nPid];
    TH1F * hRatio[nPid];
    TH1F * hUS;
    TH1F * hLS;
    TH1F * hSignal;

    TH1F * hIncE;
    
    TH1F * hNSigEPion;
    TH1F * hNSigEKaon;
    TH1F * hNSigEProton;
    TH1F * hEffBemc[nPid];
    TH1F * hEffBsmd[nPid][2];
    
    TF1 * fitfunHadron = new TF1("fitfunHadron","gaus",-13,13);
    
    for (int iPid=0; iPid<nPid; iPid++){
        hEffBemc[iPid] = new TH1F(Form("hEffBemc_%d",iPid),Form("hEffBemc_%d",iPid),5, pt);
        hEffBsmd[iPid][0] = new TH1F(Form("hEffBsmd_%d_0",iPid),Form("hEffBsmd_%d_0",iPid),5, pt);
        hEffBsmd[iPid][1] = new TH1F(Form("hEffBsmd_%d_1",iPid),Form("hEffBsmd_%d_1",iPid),5, pt);
    }
    TF1 * constant = new TF1("constant" ,"pol0", -0.1, 0.1);
    constant->SetParameter(0,1);


    outfile->cd();
    infile->Get("hTrigger")->Write();
    
    cout << "=========>START iPt loop ! " << endl;
    for (int iPid=0; iPid<nPid; iPid++){
        cout << "=========>=========>iPid : " << iPid << endl;
        hRatio[iPid] = new TH1F(Form("hRatio_%d",iPid),Form("hRatio_%d",iPid),5,pt);
        hYield[iPid] = new TH1F(Form("hYield_%d",iPid),Form("hYield_%d",iPid),5,pt);
        hMean[iPid] = new TH1F(Form("hMean_%d",iPid),Form("hMean_%d",iPid),5,pt);
        hSigma[iPid] = new TH1F(Form("hSigma_%d",iPid),Form("hSigma_%d",iPid),5,pt);
        hRawYield[iPid] = new TH1F(Form("hRawYield_%d",iPid),Form("hRawYield_%d",iPid),5,pt);
        hRatio[iPid] = new TH1F(Form("hRatio_%d",iPid),Form("hRatio_%d",iPid),5,pt);
        cout << "START iPt loop ! " << endl;
        for (int iPt=1; iPt<nPt; iPt++){
            cout << "=========>=========>=========>iPt : " << iPt << endl;
            double dpt = pt[iPt]-pt[iPt-1];
            
            hUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_0",iPt,iPid));
            hLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_0",iPt,iPid));
            hIncE = (TH1F*)infile->Get(Form("histo_%d_%d_2_0",iPt,iPid));
            
            hSignal = (TH1F*)hUS->Clone();
            hSignal->Add(hLS,-1);
            
            TF1 * fitfun = new TF1("fitfun","gaus",-3,3);
            fitfun->SetParameters(hSignal->GetMaximum(),0,1);
            fitfun->SetParLimits(1,-0.8,0.8);
            fitfun->SetParLimits(2,0.5 ,1.5);

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
         //   hSignal->Delete();
            
//            cc->cd();
       //     hYield[iPid]->SetBinContent(iPt,fitfun->Integral(-13,13)/dpt);
       //     hYield[iPid]->SetBinError(iPt,fitfun->IntegralError(-13,13)/dpt);
            hYield[iPid]->SetBinContent(iPt,hSignal->Integral()/dpt);
       //     hYield[iPid]->SetBinError(iPt,fitfun->IntegralError(-13,13)/dpt);
            
            hMean[iPid]->SetBinContent(iPt,fitfun->GetParameter(1));
            hMean[iPid]->SetBinError(iPt,fitfun->GetParError(1));
            hSigma[iPid]->SetBinContent(iPt,fitfun->GetParameter(2));
            hSigma[iPid]->SetBinError(iPt,fitfun->GetParError(2));
            
            cc3->cd()->SetLogy();
       //     hIncE->SetTitle(Form("%s, %.1f < pT < %.1f",pid[iPid].Data(),pt[iPt-1],pt[iPt]));
            hIncE->Draw();
            TF1 * fitfunPion = new TF1("fitfunPion","gaus",-6,-2.5);
            TF1 * fitfunKaon = new TF1("fitfunKaon","gaus",-10,-6);
            TF1 * fitfunMerged= new TF1("fitfunMerged","gaus",3,8);
            TF1 * fitfunE = new TF1("fitfunE","gaus",-1,1);
            TF1 * fitfunAll= new TF1("fitfunAll","gaus(0)+gaus(3)+gaus(6)+gaus(9)",-13,13);

            fitfunPion->SetParameters(hIncE->GetMaximum()-1,-3,1);
            fitfunKaon->SetParameters(hIncE->GetMaximum()-1,-7,1);
            fitfunMerged->SetParameters(hIncE->GetMaximum()/10-1,5,1);
            fitfunE->SetParameters(hIncE->GetMaximum()-1,-0.297203,0.963029);

            fitfunE->FixParameter(1,fitfun->GetParameter(1));
            fitfunE->FixParameter(2,fitfun->GetParameter(2));
            
     //       hIncE->Fit(fitfunPion,"NOR+");
     //       hIncE->Fit(fitfunKaon,"NOR+");
     //       hIncE->Fit(fitfunMerged,"NOR+");
     //       hIncE->Fit(fitfunE,"NOR+");
        
            double par[12];
            fitfunPion->GetParameters(&par[0]);
            fitfunKaon->GetParameters(&par[3]);
            fitfunMerged->GetParameters(&par[6]);
            fitfunE->GetParameters(&par[9]);


            fitfunAll->SetParameters(par);
            fitfunAll->SetParLimits(0,1,hIncE->GetMaximum());
            fitfunAll->SetParLimits(3,1,hIncE->GetMaximum());
            fitfunAll->SetParLimits(6,1,hIncE->GetMaximum()/10);
            fitfunAll->SetParLimits(9,1,hIncE->GetMaximum());

            fitfunAll->SetParLimits(1,-5,-2);
            fitfunAll->SetParLimits(2,0.5,1.5);
            
            fitfunAll->SetParLimits(4,-10,-4);
            fitfunAll->SetParLimits(5,0.5,1.5);
            
            fitfunAll->SetParLimits(7,2.,10);
            fitfunAll->SetParLimits(8,0.5,1.5);
        
            fitfunAll->SetParLimits(10,-0.5,0.5);
            fitfunAll->SetParLimits(11,0.5 ,1.5);
            //fitfunAll->FixParameter(10,-0.297203);
            //fitfunAll->FixParameter(11,0.963029);
            

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
            
            fitfunAll->SetLineWidth(4);
            fitfunAll->SetLineStyle(2);
            fitfunAll->SetLineColor(6);
            
            fitfunPion->Draw("same");
            fitfunKaon->Draw("same");
            fitfunMerged->Draw("same");
            fitfunE->Draw("same");
            
            fitfunAll->Draw("same");
            
            cc3->SaveAs(Form("~/Desktop/IncE_%d_%d.pdf",iPid,iPt));
            cout << iPid << " " << iPt  << " Purity (" << pidCutLw[iPt] << ", " << pidCutHi[iPt] << ") : " << fitfunE->Integral(pidCutLw[iPt],pidCutHi[iPt])/fitfunAll->Integral(pidCutLw[iPt],pidCutHi[iPt]) * 100 << "%" << endl;
            
            hRawYield[iPid]->SetBinContent(iPt,fitfunE->Integral(-13,13)/dpt);
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
        cc5->cd(1)->SetLogy();
        hRawYield[iPid]->SetTitle(pid[iPid]);
        hRawYield[iPid]->SetMarkerStyle(20);
        hRawYield[iPid]->SetMinimum(hYield[iPid]->GetBinContent(1)*10e-7);
        hRawYield[iPid]->SetMaximum(hYield[iPid]->GetBinContent(1)*10e2);
        hRawYield[iPid]->Draw("p");
        hYield[iPid]->Draw("psame");
        
        cc5->cd(2);
        
        hRatio[iPid]->Divide(hYield[iPid],hRawYield[iPid]);
        hRatio[iPid]->SetMaximum(0.1);
        hRatio[iPid]->SetMinimum(0.);
        hRatio[iPid]->SetTitle("");
        hRatio[iPid]->SetMarkerStyle(20);
        hRatio[iPid]->Draw("p");
        cc5->SaveAs(Form("~/Desktop/RawYield_Pid%d.pdf",iPid));
        outfile->cd();
        hRawYield[iPid]->Write();
        hYield[iPid]->Write();
    } // end iPid loop
    cout << "=========>END iPid loop ! " << endl;
    outfile->Close();
    
    //cc4->SetLogy();
    cc4->cd();
    int i = 0;
    hRatio[i]->Divide(hYield[1],hYield[0],1,1,"B");
    hRatio[i]->SetTitle("TPC+TOF / TPC (MB)");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 1;
    hRatio[i]->Divide(hYield[2],hYield[0],1,1,"B");
    hRatio[i]->SetTitle("TPC+BEMC / TPC (MB)");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 2;
    hRatio[i]->Divide(hYield[3],hYield[2],1,1,"B");
    hRatio[i]->SetTitle("TPC+BEMC+BSMD / TPC+BEMC (MB)");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    
    int i = 3;
    hRatio[i]->Divide(hYield[5],hYield[4],1,1,"B");
    hRatio[i]->SetTitle("TPC+TOF / TPC (BHT)");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 4;
    hRatio[i]->Divide(hYield[6],hYield[4],1,1,"B");
    hRatio[i]->SetTitle("TPC+BEMC / TPC (BHT)");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    
    int i = 5;
    hRatio[i]->Divide(hYield[7],hYield[6],1,1,"B");
    hRatio[i]->SetTitle("TPC+BEMC+BSMD / TPC+BEMC (BHT)");
    hRatio[i]->SetMaximum(2);
    hRatio[i]->SetMinimum(0);
    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->Draw("p");
    cc4->SaveAs(Form("~/Desktop/Ratio_%d.pdf",i));
    

    for (int iPid=0; iPid<nPid; iPid++) for (int iPt=1; iPt<nPt; iPt++) {
        
        cc8->cd(1)->SetLogy(); // neta
        TH1I * hisQaLS = (TH1I*)infile->Get(Form("histo_%d_%d_1_4",iPt, iPid));
        TH1I * hisQaUS = (TH1I*)infile->Get(Form("histo_%d_%d_0_4",iPt, iPid));
        TH1I * hisQaPE = new TH1I("hisQaPE","hisQaPE",10,0,10);
        
        hisQaPE->Add(hisQaUS,hisQaLS,1,-1);
        
        
        hisQaLS->SetMarkerStyle(20);
        hisQaLS->SetMarkerColor(2);
        hisQaLS->SetMarkerSize(0.5);
        
        
        hisQaUS->SetMarkerStyle(20);
        hisQaUS->SetMarkerColor(4);
        hisQaUS->SetMarkerSize(0.5);
        hisQaUS->SetMinimum(0.5);
        
        
        hisQaPE->SetFillStyle(3001);
        hisQaPE->SetFillColor(1);
        hisQaPE->SetLineColor(1);
        
        hisQaUS->Draw("p");
        hisQaLS->Draw("psame");
        hisQaPE->Draw("BARsame");
        
        
        double err1, err2;
        if (iPid%4 < 2) {
            hEffBsmd[iPid][0]->SetBinContent(iPt, hisQaPE->IntegralAndError(3,10,err1)/hisQaPE->IntegralAndError(2,10,err2));
            hEffBsmd[iPid][0]->SetBinError(iPt, TMath::Sqrt(err1*err1/hisQaPE->Integral(3,10)/hisQaPE->Integral(3,10) + err2*err2/hisQaPE->Integral(2,10)/hisQaPE->Integral(2,10)));
        }

        
        
        cc8->cd(2)->SetLogy(); // nphi
        TH1I * hisQaLS = (TH1I*)infile->Get(Form("histo_%d_%d_1_5",iPt, iPid));
        TH1I * hisQaUS = (TH1I*)infile->Get(Form("histo_%d_%d_0_5",iPt, iPid));
        TH1I * hisQaPE = new TH1I("hisQaPE","hisQaPE",10,0,10);
        
        hisQaPE->Add(hisQaUS,hisQaLS,1,-1);
        
        
        hisQaLS->SetMarkerStyle(20);
        hisQaLS->SetMarkerColor(2);
        hisQaLS->SetMarkerSize(0.5);
        
        
        hisQaUS->SetMarkerStyle(20);
        hisQaUS->SetMarkerColor(4);
        hisQaUS->SetMarkerSize(0.5);
        hisQaUS->SetMinimum(0.5);
        
        
        hisQaPE->SetFillStyle(3001);
        hisQaPE->SetFillColor(1);
        hisQaPE->SetLineColor(1);
        
        hisQaUS->Draw("p");
        hisQaLS->Draw("psame");
        hisQaPE->Draw("BARsame");
        if (iPid < 2) {
            hEffBsmd[iPid][1]->SetBinContent(iPt, hisQaPE->IntegralAndError(3,10,err1)/hisQaPE->IntegralAndError(2,10,err2));
            hEffBsmd[iPid][1]->SetBinError(iPt, TMath::Sqrt(err1*err1/hisQaPE->Integral(3,10)/hisQaPE->Integral(3,10) + err2*err2/hisQaPE->Integral(2,10)/hisQaPE->Integral(2,10)));
        }
        
        
        cc8->cd(3)->SetLogy(); // nphieta
        TH1I * hisQaLS = (TH1I*)infile->Get(Form("histo_%d_%d_1_11",iPt, iPid));
        TH1I * hisQaUS = (TH1I*)infile->Get(Form("histo_%d_%d_0_11",iPt, iPid));
        TH1I * hisQaPE = new TH1I("hisQaPE","hisQaPE",20,0,20);
        
        hisQaPE->Add(hisQaUS,hisQaLS,1,-1);
        
        
        hisQaLS->SetMarkerStyle(20);
        hisQaLS->SetMarkerColor(2);
        hisQaLS->SetMarkerSize(0.5);
        
        
        hisQaUS->SetMarkerStyle(20);
        hisQaUS->SetMarkerColor(4);
        hisQaUS->SetMarkerSize(0.5);
        hisQaUS->SetMinimum(0.5);
        
        
        hisQaPE->SetFillStyle(3001);
        hisQaPE->SetFillColor(1);
        hisQaPE->SetLineColor(1);
        
        hisQaUS->Draw("p");
        hisQaLS->Draw("psame");
        hisQaPE->Draw("BARsame");
        
        
        
        cc8->cd(4)->SetLogy(); // e0/p
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_6",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_6",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,0,3);
        
        hisPE->Add(hisUS,hisLS,1,-1);
        
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(0.5);
        
        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(0.5);
        hisUS->SetMinimum(0.5);
        
        
        hisPE->SetFillStyle(3001);
        hisPE->SetFillColor(1);
        hisPE->SetLineColor(1);
        
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("BARsame");
        
        if(iPid%4 < 2){
            fitfunHadron->SetParameters(hisPE->GetMaximum(),0.8,1);
            fitfunHadron->SetRange(0.4,1.5);
            hisPE->Fit(fitfunHadron,"NOR+");
            
            fitfunHadron->SetRange(0.,3);
            fitfunHadron->SetLineColor(2);
            fitfunHadron->SetLineStyle(2);
            fitfunHadron->Draw("same");
        }
        
        
        
        cc8->cd(5)->SetLogy(); // zDist
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_7",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_7",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,-20,20);
        
        hisPE->Add(hisUS,hisLS,1,-1);
        
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(0.5);
        
        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(0.5);
        hisUS->SetMinimum(0.5);
        
        
        hisPE->SetFillStyle(3001);
        hisPE->SetFillColor(1);
        hisPE->SetLineColor(1);
        
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("BARsame");
        

        
        cc8->cd(6)->SetLogy(); // phiDist
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_8",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_8",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,-0.1,0.1);
        
        hisPE->Add(hisUS,hisLS,1,-1);
        
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(0.5);
        
        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(0.5);
        hisUS->SetMinimum(0.5);
        
        
        hisPE->SetFillStyle(3001);
        hisPE->SetFillColor(1);
        hisPE->SetLineColor(1);
        
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("BARsame");
        
        
        
        cc8->cd(7)->SetLogy(); // etaTowDist
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_9",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_9",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,-0.1,0.1);
        
        hisPE->Add(hisUS,hisLS,1,-1);
        
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(0.5);
        
        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(0.5);
        hisUS->SetMinimum(0.5);
        
        
        hisPE->SetFillStyle(3001);
        hisPE->SetFillColor(1);
        hisPE->SetLineColor(1);
        
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("BARsame");
        

        cc8->cd(8)->SetLogy(); // phiTowDist
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_10",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_10",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,-0.1,0.1);
        
        hisPE->Add(hisUS,hisLS,1,-1);
        
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(0.5);
        
        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(0.5);
        hisUS->SetMinimum(0.5);
        
        
        hisPE->SetFillStyle(3001);
        hisPE->SetFillColor(1);
        hisPE->SetLineColor(1);
        
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("BARsame");
        

        
        
        
        
        
        
   //     cc8->SaveAs(Form("~/Desktop/QaBemcBsmd_Pid%d_Pt%d.pdf",iPid,iPt));
        
        
        
        // hadron
        int jj;
        for (int j=0; j<8; j++) {
            cc8->cd(j+1+8)->SetLogy(); // neta
            jj=j;
            if (j==2) jj = 7;
            else if (j > 2) jj=j-1;
            TH1I * hisQaPE = (TH1I*)infile->Get(Form("histo_%d_%d_3_%d",iPt, iPid, jj+4));
            
            hisQaPE->SetFillStyle(3001);
            hisQaPE->SetFillColor(2);
            hisQaPE->SetLineColor(1);
            
            hisQaPE->Draw("BAR");

        }
        
        
        cc8->SaveAs(Form("~/Desktop/QaHadronBemcBsmd_Pid%d_Pt%d.pdf",iPid,iPt));
        if (iPid%4 < 2) {
            
            cout << fitfunHadron->Integral(0,3)<< " " << fitfunHadron->Integral(0.8,2) /fitfunHadron->Integral(0.,3) << endl;
            cout << fitfunHadron->GetParameter(0) << " " << fitfunHadron->GetParameter(1) << " " << fitfunHadron->GetParameter(2) << endl;
            // end hadron
            hEffBemc[iPid]->SetBinContent(iPt, fitfunHadron->Integral(0.8,2) /fitfunHadron->Integral(0.,3));
            hEffBemc[iPid]->SetBinError(iPt, TMath::Sqrt(fitfunHadron->IntegralError(0.8,2)*fitfunHadron->IntegralError(0.8,2)/fitfunHadron->Integral(0.8,2)/fitfunHadron->Integral(0.8,2) + fitfunHadron->IntegralError(0.,3)*fitfunHadron->IntegralError(0.,3)/fitfunHadron->Integral(0.,3)/fitfunHadron->Integral(0.,3)));
        }
        if (iPid == 1 && iPt == 5)    {
            
            
            
            
            cc11->cd();
            
            hEffBemc[0]->SetMaximum(1.1);
            hEffBemc[0]->SetMinimum(0);
            hEffBemc[0]->SetLineColor(2);
            hEffBemc[0]->SetTitle("BEMC PID Efficiency");
            hEffBemc[0]->Draw("");
            hEffBemc[1]->Draw("same");
            hEffBemc[4]->Draw("same");
            hEffBemc[5]->Draw("same");
            cc11->SaveAs("~/Desktop/EffBemc.png");
            
            hEffBsmd[0][0]->SetMaximum(1.5);
            hEffBsmd[0][0]->SetMinimum(0);
            hEffBsmd[0][0]->SetLineColor(2);
            hEffBsmd[0][0]->SetTitle("BSMD PID Efficiency");
            
            hEffBsmd[0][1]->SetLineColor(7);
            hEffBsmd[1][1]->SetLineColor(3);

            hEffBsmd[0][0]->Draw("");
            hEffBsmd[1][0]->Draw("same");
            hEffBsmd[4][0]->Draw("same");
            hEffBsmd[5][0]->Draw("same");
            
            hEffBsmd[0][1]->Draw("same");
            hEffBsmd[1][1]->Draw("same");
            hEffBsmd[4][1]->Draw("same");
            hEffBsmd[5][1]->Draw("same");
            cc11->SaveAs("~/Desktop/EffBsmd.png");
        }
    }

    cc6->cd();
    
    for (int iPid=0; iPid<nPid; iPid++) for (int iPt=1; iPt<nPt; iPt++) {
        cc6->cd(1)->SetLogy();
        TH1F * dum = (TH1F*)infile->Get(Form("histo_%d_%d_2_1",iPt, iPid));
        TH1F * his = (TH1F*)infile->Get(Form("histo_%d_%d_2_2",iPt, iPid));
        TH1F * hisLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_2",iPt, iPid));
        TH1F * hisUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_2",iPt, iPid));
        TH1F * hisPE = new TH1F("hisPE","hisPE",100,-0.1,0.1);
        
        hisPE->Add(hisUS,hisLS,1,-1);
        
        
        hisLS->SetMarkerStyle(20);
        hisLS->SetMarkerColor(2);
        hisLS->SetMarkerSize(1);
        
        
        hisUS->SetMarkerStyle(20);
        hisUS->SetMarkerColor(4);
        hisUS->SetMarkerSize(1);
        hisUS->SetMinimum(0.5);
        
        
        hisPE->SetFillStyle(3004);
        hisPE->SetFillColor(1);
        hisPE->SetLineColor(1);
        
        hisUS->Draw("p");
        hisLS->Draw("psame");
        hisPE->Draw("BARsame");
        
        cc6->cd(2);
        TH1F * hPairMassUS = (TH1F*)infile->Get(Form("histo_%d_%d_0_3",iPt, iPid));
        TH1F * hPairMassLS = (TH1F*)infile->Get(Form("histo_%d_%d_1_3",iPt, iPid));
        TH1F * hdum = new TH1F("hdum","hdum",100,0,0.2);


        
        hPairMassUS->GetXaxis()->SetRangeUser(0,0.1);
        hPairMassUS->SetMinimum(-1*hPairMassUS->GetMaximum()*0.1);
        hPairMassUS->SetMarkerStyle(20);
        hPairMassUS->SetMarkerColor(4);
        hPairMassUS->SetMarkerSize(1);
        
        hPairMassLS->SetMarkerStyle(20);
        hPairMassLS->SetMarkerColor(2);
        hPairMassLS->SetMarkerSize(1);
        
        hdum->GetXaxis()->SetRangeUser(0,0.1);
        hdum->Add(hPairMassUS,hPairMassLS,1,-1);
        hdum->SetFillStyle(3004);
        hdum->SetFillColor(1);
        hdum->SetLineColor(1);
        
        
        
        
        hPairMassUS->Draw("p");
        hPairMassLS->Draw("psame");
        hdum->Draw("BARsame");
        TLine * line = new TLine(0,0,0.1,0);
        line->SetLineColor(1);
        line->SetLineStyle(2);
        line->SetLineWidth(4);
        line->Draw("same");
        cc6->SaveAs(Form("~/Desktop/DcaAterPid_Pid%d_Pt%d.pdf",iPid,iPt));
        
        
        
        cc7->cd(1)->SetLogy();
        double glo = 0.1;
        hisPE->Divide(constant,hisPE->Integral(1,100)*10);
        dum->Divide(constant,dum->Integral(1,100)*90);
        his->Divide(constant,his->Integral(1,100)*10);
        
        his->SetMarkerStyle(20);
        his->SetMarkerColor(1);
        his->SetMarkerSize(0.5);
        his->SetMaximum(0.1);
        his->SetMinimum(1e-6);
        
        dum->SetMarkerStyle(20);
        dum->SetMarkerColor(2);
        dum->SetMarkerSize(0.5);
        
        
        his->DrawClone("p");
        dum->DrawClone("psame");
        //  hisPE->DrawClone("psame");
        
        
        
        cc7->cd(2)->SetLogy();
        TH1F * hisdum = new TH1F("hisdum","hisdum",100,-0.1,0.1);
        hisdum->Add(his,dum,1,-1);
        hisdum->SetMaximum(0.1);
        hisdum->SetMinimum(1e-6);
        hisdum->Draw("p");
        hisPE->DrawClone("psame");
        
        cc7->SaveAs(Form("~/Desktop/DcaAterPid_Fit_Pid%d_Pt%d.pdf",iPid,iPt));
        
        delete hdum;
        delete hisdum;
        
        delete line;
 
    }
    //return 0 ;
    

    cc10->cd();
    infile->Get("hTrigger")->Draw();
    cc10->SaveAs("~/Desktop/Trigger.png");
    cc10->cd();
    infile->Get("hTriggerWt")->Draw();
    cc10->SaveAs("~/Desktop/TriggerWt.png");
}

double funHF(double * x, double * par){
    if (x[0]>0) return TMath::Exp(-1*x[0]*par[1]+par[2])*par[0];
    else return TMath::Exp(x[0]*par[1]+par[2])*par[0];
}
