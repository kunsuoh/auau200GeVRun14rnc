void purity(){
    double par[9] = {  4.04863e+01,
        -3.24209e+00,
        1.78641e+00,
        9.76785e+00,
        3.05655e-01,
        6.92029e-01};
    TF1 * fitfun = new TF1("fitfun","gaus(0)+gaus(3)+gaus(6)",-13,13);
    TF1 * fit2 = new TF1("fit2","gaus",-13,13);
    fitfun->SetParameters(par);
    fit2->SetParameters(par[3],par[4],par[5]);
    
    
    TCanvas * cc = new TCanvas("cc","cc",600,600);
    cc->SetLogy();
    TH1D * dum = new TH1D("dum","dum",1000,-13,13);
    dum->SetMinimum(0.7);
    dum->SetMaximum(1000);
    dum->Draw();
    fit2->SetLineColor(4);
    fit2->SetLineStyle(1);
    fit2->SetLineWidth(0.5);
    fit2->SetRange(-13,13);
    fit2->Draw("same");
    
    fitfun->SetLineColor(2);
    fitfun->SetLineStyle(2);
    fitfun->SetRange(-13,13);
    fitfun->Draw("same");
    
    cout << fit2->Integral(0,2) / fitfun->Integral(0,2) << endl;

}
