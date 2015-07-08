void drawHistogram(int event = 1, int source = 1) // event 0:MB, 1:CC, 2:BHT1, 3:BHT2, 4:BHT3
                                                  // source 0:inclusive, 1:photonic
{
    const Int_t nbin = 14;
    const Int_t bins[nbin] = {85, 100, 100, 289, 200, 10, 10, 100, 100, 100, 100, 100, 100, 3};
    const Double_t xmin[nbin] = {1.5, -0.7, -0.1, -13, 0, 0, 0, -20, -0.1, -0.1, -0.1, 0, 0 , -1.5};
    const Double_t xmax[nbin] = {10, 0.7, 0.1, 13, 4, 10, 10, 20, 0.1, 0.1, 0.1, 0.4, 1, 1.5};
    
    int coord[nbin];
    int knbin[nbin];
    double kxmin[nbin];
    double kxmax[nbin];
    double val[nbin];
    const int nHisto = 22;
    TH1D * h[nHisto][nbin];
    for (int i=0; i<nbin; i++) for (int j=0; j<nHisto; j++) {
        h[j][i] = new TH1D(Form("h%d%d",j,i),Form("h%d%d",j,i),bins[i],xmin[i],xmax[i]);
    }
    for (int i=0; i<708; i++)
    {
        if (i%10==0) cout << i << endl;
        if (i==433 || i==639 || i==640) continue;
        TFile * infile = new TFile(Form("production/356C0EB17CA7336CE6C00B6590B1CB12_%d.root",i));
        THnSparseF * hs = (THnSparseF*)infile->Get(Form("sparse_%d_%d", event, source));
        int ntrack = hs->GetNbins();
        for (int j=0; j<ntrack; j++)
        {
            hs->GetBinContent(j,coord);
            for (int k = 0; k<nbin; k++)
            {
                knbin[k] = hs->GetAxis(k)->GetNbins();
                kxmin[k] = hs->GetAxis(k)->GetXmin();
                kxmax[k] = hs->GetAxis(k)->GetXmax();
                if (j==0 && i==0) cout << knbin[k] << " " << kxmin[k] << " " << kxmax[k] << endl;
                val[k] = coord[k]*(kxmax[k]-kxmin[k])/knbin[k]+kxmin[k]-(kxmax[k]-kxmin[k])/knbin[k]*0.5;
                
                h[1][k]->Fill(val[k]);
                
                if (source == 1) {
                    if (val[11]<0.05 && val[13]==0) h[2][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0) h[3][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1) {
                        h[4][k]->Fill(val[k]);
                        h[(const int)getPt(val[0])*2+6][k]->Fill(val[k]);
                    }
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1) {
                        h[5][k]->Fill(val[k]);
                        h[(const int)getPt(val[0])*2+7][k]->Fill(val[k]);

                    }
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[3] > 0) h[6][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[3] > 0) h[7][k]->Fill(val[k]);
                    /*
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] < 1.7) h[8][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] < 1.7) h[9][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] > 1.7 && val[0] < 2.0) h[10][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] > 1.7 && val[0] < 2.0) h[11][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] > 2.0 && val[0] < 2.5) h[12][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] > 2.0 && val[0] < 2.5) h[13][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] > 2.5 && val[0] < 3.0) h[14][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] > 2.5 && val[0] < 3.0) h[15][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] > 3.0 && val[0] < 4.0) h[16][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] > 3.0 && val[0] < 4.0) h[17][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] > 4.0 && val[0] < 6.0) h[18][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] > 4.0 && val[0] < 6.0) h[19][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]==0 && val[12] < 0.1 && val[0] > 6.0 && val[0] < 10.) h[20][k]->Fill(val[k]);
                    if (val[11]<0.05 && val[13]!=0 && val[12] < 0.1 && val[0] > 6.0 && val[0] < 10.) h[21][k]->Fill(val[k]);
                     */
                }
                if (source == 0){
                    if (val[4] > 0.8) h[2][k]->Fill(val[k]);
                    if (val[4] > 0.8 && val[5] > 1 ) h[3][k]->Fill(val[k]);
                    if (val[4] > 0.8 && val[6] > 1 ) h[4][k]->Fill(val[k]);
                    if (val[4] > 0.8 && val[5] > 1 && val[6] > 1 ) h[5][k]->Fill(val[k]);
                    if (val[4] > 0.8 && val[5] > 1 && val[6] > 1 && abs(val[7]) < 2 ) h[6][k]->Fill(val[k]);
                    if (val[4] > 0.8 && val[5] > 1 && val[6] > 1 && abs(val[8]) < 0.1 ) h[7][k]->Fill(val[k]);
                    if (val[4] > 0.8 && val[5] > 1 && val[6] > 1 && abs(val[7]) < 2 && abs(val[8]) < 0.1 ) h[8][k]->Fill(val[k]);
                    if (val[0] < 1.7) h[9][k]->Fill(val[k]);
                    if (val[0] > 1.7 && val[0] < 2.0 ) h[10][k]->Fill(val[k]);
                    if (val[0] > 2.0 && val[0] < 2.5 ) h[11][k]->Fill(val[k]);
                    if (val[0] > 2.5 && val[0] < 3.0 ) h[12][k]->Fill(val[k]);
                    if (val[0] > 3.0 && val[0] < 4.0 ) h[13][k]->Fill(val[k]);
                    if (val[0] > 4.0 && val[0] < 6.0 ) h[14][k]->Fill(val[k]);
                    if (val[0] > 6.0 && val[0] < 10. ) h[15][k]->Fill(val[k]);
                }
            }
            
        }
        
        infile->Close();
        hs->Delete();
    }
    
    TFile * f = new TFile(Form("outfile_%d_%d.root",event,source),"recreate");
    f->cd();
    
    int nMax;
    if (source==1) nMax=22;
    if (source==0) nMax=16;

    TDirectory * cd[nHisto];
    if (source==1) {
        cd[1] = f->mkdir("noCut");
        cd[2] = f->mkdir("massCutUS");
        cd[3] = f->mkdir("massCutLS");
        cd[4] = f->mkdir("massCutpairDcaCutUS");
        cd[5] = f->mkdir("massCutpairDcaCutLS");
        cd[6] = f->mkdir("massCutpairDcaCutnSigECutUS");
        cd[7] = f->mkdir("massCutpairDcaCutnSigECutLS");
        cd[8] = f->mkdir("cdpt1US");
        cd[9] = f->mkdir("cdpt1LS");
        cd[10] = f->mkdir("cdpt2US");
        cd[11] = f->mkdir("cdpt2LS");
        cd[12] = f->mkdir("cdpt3US");
        cd[13] = f->mkdir("cdpt3LS");
        cd[14] = f->mkdir("cdpt4US");
        cd[15] = f->mkdir("cdpt4LS");
        cd[16] = f->mkdir("cdpt5US");
        cd[17] = f->mkdir("cdpt5LS");
        cd[18] = f->mkdir("cdpt6US");
        cd[19] = f->mkdir("cdpt6LS");
        cd[20] = f->mkdir("cdpt7US");
        cd[21] = f->mkdir("cdpt7LS");
    }
    if (source==0) {
        cd[1] = f->mkdir("noCut"); // TPC only
        cd[2] = f->mkdir("eoverpCut"); // TPC+BEMC
        cd[3] = f->mkdir("eoverpnEta");
        cd[4] = f->mkdir("eoverpnPhi");
        cd[5] = f->mkdir("eoverpnEtanPhi"); // TPC+BEMC+BSMD
        cd[6] = f->mkdir("eoverpnEtanPhizDist");
        cd[7] = f->mkdir("eoverpnEtanPhiphiDist");
        cd[8] = f->mkdir("eoverpnEtanPhizDistphiDist"); // TPC+BEMC+BSMD+Enhancement
        cd[9] = f->mkdir("cdpt1");
        cd[10] = f->mkdir("cdpt2");
        cd[11] = f->mkdir("cdpt3");
        cd[12] = f->mkdir("cdpt4");
        cd[13] = f->mkdir("cdpt5");
        cd[14] = f->mkdir("cdpt6");
        cd[15] = f->mkdir("cdpt7");
    }
    for (int j=1; j<nMax; j++){
        cd[j]->cd();
        for (int i=0; i<nbin; i++) h[j][i]->Write();
    }
    f->Close();
    
}

int getPt(double pt)
{
    int nPt = 0;
    if (pt < 1.7) nPt = 1;
    else if (pt < 2.0) nPt = 2;
    else if (pt < 2.5) nPt = 3;
    else if (pt < 3.0) nPt = 4;
    else if (pt < 4.0) nPt = 5;
    else if (pt < 6.0) nPt = 6;
    else nPt = 7;

    return nPt;
}



