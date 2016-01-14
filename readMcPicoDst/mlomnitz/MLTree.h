#ifndef KFPartTree_hh
#define KFPartTree_hh

TTree *primaryvtx;
struct primVtxPoint_t{ 
  Int_t eventp,l;
  Float_t primX,primY,primZ;
  Float_t MCpvX,MCpvY,MCpvZ;
  Int_t MultP,refMultP;
};
primVtxPoint_t primVtx;

TTree *primtrack;
struct TrackPoint_t {
  Int_t id, GPid;
  Float_t dEdx;
  Float_t pMag;
  Float_t pt, MCpt; 
  Float_t eta;
  Float_t phi;
  Float_t dcaxy;
  Float_t dcaz;
  Int_t charge;
  Int_t tpc;
  Int_t istpos, istfit, isttop;
  Int_t pxlpos, pxlfit, pxltop, MCpxl;
  Float_t SigPi; Float_t SigKa;
  Float_t SigEl; Float_t SigPr;
  //FOR Siimu use only
  Int_t isD0, isK0;
  };  
TrackPoint_t Track;

TTree *mcInfo;
struct MCPoint_t {
  Int_t         GePid;
  Float_t       dist; 
  Float_t       proj; 
  Int_t         daughterId[5]={-999};
  Float_t       phi,eta;
  Float_t       pX,pY,pZ;
  //  Float_t       Mass;
  //Float_t      dcaxy, dcaz, dcamag, Prob;    
};
BPoint_t mcD0;
TTree *secKF;
struct BPoint_t {
  Int_t         q;
  Float_t       dist, KFX, KFY, KFZ; 
  Float_t       proj; 
  Float_t       dcadaughter;
  Float_t       phi;
  Float_t       lambda;
  Float_t       mom;
  Float_t       dca2vtx;//ist,pxl;
  Float_t       lifet, Mass, MassRot;
  Float_t      dcaxy, dcaz, dcamag, Prob;    
};
BPoint_t B,D0_B;

TTree *secHE;
struct B2Point_t {
  Int_t         q;
  Float_t       dist, HEX, HEY, HEZ; 
  Float_t       proj; 
  Float_t       phi;
  Float_t       trk1dca, trk2dca;
  Float_t       lambda;
  Float_t       mom;
  Float_t       dca2vtx;//ist,pxl;
  Float_t       lifet, Mass, MassRot;
  Float_t      dcaxy, dcaz, dcamag;    
};
B2Point_t B2,D0_B2;

void DefineTree(){
 
  primaryvtx = new TTree("primaryvtx","The Primary Vertices");
  primaryvtx->SetAutoSave(1000000);
  primaryvtx->Branch("eventp",&primVtx.eventp,"eventp/I");
  primaryvtx->Branch("l",&primVtx.l,"l/I");
  primaryvtx->Branch("primX",&primVtx.primX,"primX/F");
  primaryvtx->Branch("primY",&primVtx.primY,"primY/F");
  primaryvtx->Branch("primZ",&primVtx.primZ,"primZ/F");
  primaryvtx->Branch("MCpvX",&primVtx.MCpvX,"MCpvX/F");
  primaryvtx->Branch("MCpvY",&primVtx.MCpvY,"MCpvY/F");
  primaryvtx->Branch("MCpvZ",&primVtx.MCpvZ,"MCpvZ/F");
  primaryvtx->Branch("MultP",&primVtx.MultP,"MultP/I");
  primaryvtx->Branch("refMultP",&primVtx.refMultP,"refMultP/I");
  
  primtrack = new TTree("primtrack","The Primary Tracks");
  primtrack->SetAutoSave(1000000);
  primtrack->Branch("id",&Track.id,"id/I");
  primtrack->Branch("GPid",&Track.GPid,"GPid/I");
  primtrack->Branch("dEdx",&Track.dEdx,"dEdx/F");
  primtrack->Branch("pMag",&Track.pMag,"pMag/F");
  primtrack->Branch("pt",&Track.pt,"pt/F");
  primtrack->Branch("MCpt",&Track.MCpt,"MCpt/F");
  primtrack->Branch("eta",&Track.eta,"eta/F");
  primtrack->Branch("phi",&Track.phi,"phi/F");
  primtrack->Branch("dcaxy",&Track.dcaxy,"dcaxy/F");
  primtrack->Branch("dcaz",&Track.dcaz,"dcaz/F");
  primtrack->Branch("charge",&Track.charge,"charge/I");
  primtrack->Branch("tpc",&Track.tpc,"tpc/I");
  primtrack->Branch("istpos",&Track.istpos,"istpos/I");
  primtrack->Branch("istfit",&Track.istfit,"istfit/I");
  primtrack->Branch("isttop",&Track.isttop,"isttop/I");
  primtrack->Branch("pxlpos",&Track.pxlpos,"pxlpos/I");
  primtrack->Branch("MCpxl",&Track.MCpxl,"MCpxl/I");
  primtrack->Branch("pxlfit",&Track.pxlfit,"pxlfit/I");
  primtrack->Branch("pxltop",&Track.pxltop,"pxltop/I");
  primtrack->Branch("SigPi",&Track.SigPi,"SigPi/F");
  primtrack->Branch("SigKa",&Track.SigKa,"SigKa/F");
  primtrack->Branch("SigPr",&Track.SigPr,"SigPr/F");
  primtrack->Branch("SigEl",&Track.SigEl,"SigEl/F");
  primtrack->Branch("isD0",&Track.isD0,"isD0/I");
  primtrack->Branch("isK0",&Track.isK0,"isK0/I");

  secKF = new TTree("secKF","The Secondary Vertices via KFParticles");
  secKF->SetAutoSave(1000000);
  //secKF->Branch("eventt",&B.eventt,"eventt/I");
  secKF->Branch("q",&B.q,"q/I");
  secKF->Branch("dist",&B.dist,"dist/F");
  secKF->Branch("KFX",&B.KFX,"KFX/F");
  secKF->Branch("KFY",&B.KFY,"KFY/F");
  secKF->Branch("KFZ",&B.KFZ,"KFZ/F");
  secKF->Branch("proj",&B.proj,"proj/F");
  secKF->Branch("phi",&B.phi,"phi/F");
  secKF->Branch("lambda",&B.lambda,"lambda/F");
  secKF->Branch("dcadaughter",&B.dcadaughter,"dcadaughter/F");
  secKF->Branch("mom",&B.mom,"mom/F");
  secKF->Branch("dca2vtx",&B.dca2vtx,"dca2vtx/F");
  secKF->Branch("lifet",&B.lifet,"lifet/F");
  secKF->Branch("Mass",&B.Mass,"Mass/F");
  secKF->Branch("MassRot",&B.MassRot,"MassRot/F");
  secKF->Branch("dcaxy",&B.dcaxy,"dcaxy/F");
  secKF->Branch("dcaz",&B.dcaz,"dcaz/F");
  secKF->Branch("dcamag",&B.dcamag,"dcamag/F");
  secKF->Branch("Prob",&B.Prob,"Prob/F");

  secHE = new TTree("secHE","The Secondary Vertices via HEParticles");
  secHE->SetAutoSave(1000000);
  //secHE->Branch("eventt",&B.eventt,"eventt/I");
  secHE->Branch("q",&B2.q,"q/I");
  secHE->Branch("dist",&B2.dist,"dist/F");
  secHE->Branch("HEX",&B2.HEX,"HEX/F");
  secHE->Branch("HEY",&B2.HEY,"HEY/F");
  secHE->Branch("HEZ",&B2.HEZ,"HEZ/F");
  secHE->Branch("proj",&B2.proj,"proj/F");
  secHE->Branch("phi",&B2.phi,"phi/F");
  secHE->Branch("lambda",&B2.lambda,"lambda/F");
  secHE->Branch("mom",&B2.mom,"mom/F");
  secHE->Branch("dca2vtx",&B2.dca2vtx,"dca2vtx/F");
  secHE->Branch("trk1dca",&B2.trk1dca,"trk1dca/F");
  secHE->Branch("trk2dca",&B2.trk2dca,"trk2dca/F");
  secHE->Branch("lifet",&B2.lifet,"lifet/F");
  secHE->Branch("Mass",&B2.Mass,"Mass/F");
  secHE->Branch("MassRot",&B2.MassRot,"MassRot/F");
  secHE->Branch("dcaxy",&B2.dcaxy,"dcaxy/F");
  secHE->Branch("dcaz",&B2.dcaz,"dcaz/F");
  secHE->Branch("dcamag",&B2.dcamag,"dcamag/F");

  mcInfo= new TTree("mcInfo","MC information for selected particles");
  mcInfo->SetAutoSave(1000000);
  mcInfo->Branch("GePid",&mcD0.GePid,"GePid/I");
  mcInfo->Branch("nDaughter",&mcD0.nDaughter,"nDaughter/I");
  mcInfo->Branch("daughterId",&mcD0.daughterId,"daughterId/I[5]");
  mcInfo->Branch("phi",&mcD0.phi,"phi/F");
  mcInfo->Branch("eta",&mcD0.eta,"eta/F");
  mcInfo->Branch("pX",&mcD0.pX,"pX/F");
  mcInfo->Branch("pY",&mcD0.pY,"pY/F");
  mcInfo->Branch("pZ",&mcD0.pZ,"pZ/F");    
}

#endif
