// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C
//
// By Y. Zhang  07/29/2014
// Modified from Jason 's macro
// Added real distributions for pT, y


class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarKinematics;
StarKinematics *kinematics = 0;

TF1 *ptDist  = 0;
TF1 *yDist = 0;

//Initialize the settings:
Float_t vx     = 0.;
Float_t vy     = 0.;
Float_t vz     = 0.;
Float_t vx_sig = 0.01;
Float_t vy_sig = 0.01;
Float_t vz_sig = 2.0;
//Float_t minVz  = -5.0;
//Float_t maxVz  = +5.0;
Float_t minPt  = +0.1;
Float_t maxPt  = +10.0;
Float_t minY   = -1.0;
Float_t maxY   = +1.0;

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag;
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=0 )
{
  for ( Int_t i=0; i<n+1; i++ ) {
    chain->Clear();
    // Throw 5 D0 flat in -1 to 1 w/ momentum btwn 0.1 and 5 GeV
    //if(kinematics) kinematics->Kine( 5, "D0", minPt, maxPt, minY, maxY );

    // Generate 5 D0 according to a PT and Y distribution
    //if(kinematics) kinematics->Dist(500, "pi0", ptDist, yDist );
    if(kinematics) kinematics->Kine(100, "gamma", 0.0, 20.0, -1.0, +1.0 );;

    chain->Make();
    //    command("gprint kine");
  }
}
// ----------------------------------------------------------------------------
void myKine()
{ 
  kinematics = new StarKinematics("gamma");
  _primary -> AddGenerator(kinematics);
}
// ----------------------------------------------------------------------------
void Hijing()
{
  StarHijing *hijing = new StarHijing("hijing");
  hijing->SetTitle("Hijing 1.383");

  // Setup collision frame, energy and beam species
  hijing->SetFrame("CMS",200.0);
  hijing->SetBlue("Au");
  hijing->SetYell("Au");  
  //hijing->SetImpact(0.0, 30.0);       // Impact parameter min/max (fm)    0.   30.
  hijing->SetImpact(0.0, 15.0);       // Minimum Bias
  hijing->hiparnt().ihpr2(4) = 0;     // Jet quenching (1=yes/0=no)       0
  hijing->hiparnt().ihpr2(3) = 0;     // Hard scattering (1=yes/0=no)
  hijing->hiparnt().hipr1(10) = 2.0;  //    pT jet
  hijing->hiparnt().ihpr2(8)  = 10;   // Max number of jets / nucleon
  hijing->hiparnt().ihpr2(11) = 1;    // Set baryon production
  hijing->hiparnt().ihpr2(12) = 1;    // Turn on/off decay of particles [1=recommended]
  hijing->hiparnt().ihpr2(18) = 1;    // Turn on/off B production
  hijing->hiparnt().hipr1(7) = 5.35;  // Set B production ???? Not really used... Really ????

  // For more configuration options, see the HIJING manual
  // http://ntc0.lbl.gov/~xnwang/hijing/doc.html

  _primary -> AddGenerator(hijing);
  _primary -> SetCuts( 1.0E-6 , -1., -2.6, 2.6 );
  
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=1, Int_t Index = 0, Int_t rngSeed=4321 )
{ 

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "y2014 geant gstar usexgeom agml ";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  gSystem->Load( "libMathMore.so"   );  
  gSystem->Load( "libHijing1_383.so");
  gSystem->Load( "libKinematics.so");
  gSystem->Load( "xgeometry.so"     );

  // force gstar load/call
  gSystem->Load( "gstar.so" );
  command("call gstar");

  // Setup RNG seed and map all ROOT TRandom here
  StarRandom::seed( rngSeed );
  StarRandom::capture();

  char rootname[100],fzname[100];
  sprintf(rootname,"gamma_%d.starsim.root",Index);
  sprintf(fzname,"gfile o gamma_%d.starsim.fzd",Index);

  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName(rootname);
    chain -> AddBefore( "geant", _primary );
  }

  //
  // These should be adjusted to your best vertex estimates
  //
  _primary -> SetVertex( vx,vy,vz );
  _primary -> SetSigma( vx_sig,vy_sig,vz_sig );




  // Setup an event generator
  //
  // Hijing();
  //
  // Setup single particle
  //  
  myKine(); 
 

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();

  //
  // Setup geometry and set starsim to use agusread for input
  //
  //geometry("y2014");
  command("gkine -4 0");
  command(fzname);

  //Double_t pt0 = 3.0;
  //ptDist = new TF1("ptDist","(x/[0])/(1+(x/[0])^2)^6",0.0,10.0);
  //ptDist->SetParameter(0, pt0);
  //ptDist->Draw();
//  ptDist = new TF1("ptDist","[0]*x*TMath::Exp(-x/[1])",minPt,maxPt); //dN/pT/dpT is exp 
//  ptDist->SetParameters(1.,1.);//slope = 1.;
//  ptDist = new TF1("ptDist", funModiHagedorn_pt,1.0,10,6);
//  ptDist->SetParameters(5.655567e+02,5.072933e-01,1.215052e-01,7.709202e-01,8.340029e+00,1.395702e-01);
  //yDist = new TF1("yDist","-TMath::Erf(x+2.6)*TMath::Erf(x-2.6)",minY,maxY);
//  yDist = new TF1("yDist","pol0",minY,maxY);
//  yDist->SetParameter(0,1.);
  
  //phi, default 0 ~ TMath::TwoPi() flat

  //
  // Trigger on nevents
  //
  trig( nevents );

  _primary->event()->Print();

  //  command("gprint kine");

  command("call agexit");  // Make sure that STARSIM exits properly

}
// ----------------------------------------------------------------------------
double funModiHagedorn_pt(double * x, double * par) {
    double scale = par[0];
    double a = par[1];
    double b = par[2];
    double p0= par[3];
    double n = par[4];
    double pt = TMath::Sqrt(x[0]*x[0] + par[5]*par[5] - 0.13957018*0.13957018);
    double dndy = scale * TMath::Power(TMath::Exp( (-1) * (a*pt + b*TMath::Power(pt,2)) ) + pt/p0, n*(-1));
    
    return dndy*pt;
}
