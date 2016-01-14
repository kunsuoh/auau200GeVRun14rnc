// Read the MuDst of a simulation
//Amilkar Quintero 
// September 2014
// The file AQHistosgrams.h define the used histograms
//   to run > root.exe lMuDst.C 'RecoSim.C+(9999999,"/path/to/files/*MuDst.root")'

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <string>
#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
//#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
//#include "TObjString.h"
//#include "TArrayF.h"
//#include "TArrayD.h"
#include "TVector3.h"
#include "TRVector.h"
#include "TLorentzVector.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
class StMuDstMaker;

#define ClassStMessMgr
#define StMessMgr Int_t
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#undef  StMessMgr
#undef ClassStMessMgr
//StMuDstMaker* maker = 0;
#include "TDatime.h"
#include "TBranch.h"
#include "TMath.h"
#include "TStyle.h"
//#include "THelixTrack.h"
#include "StBTofHeader.h"

//#include "KFVertex.h"
#include "KFParticle.h"
#include "MVertex.h"
#include "MTrack.h"
#endif

StMuDstMaker* maker = 0;
#include "AQTree1.h"
#include "AppKFPart.h"
#include "AppHelix.h"
//#include "AppTCFit.h"
#include "/star/institutions/ksu/aquinter/V0Finder/Ask.h"
#include "AQHist.h"

//________________________________________________________________________________
// Conditions to accept the vertices
Bool_t AcceptVtx(const StMuPrimaryVertex *Vtx = 0) {
  if (! Vtx) return kFALSE;
  if (abs(Vtx->position().z()) > 5) return kFALSE;
  if (Vtx->refMult()<30) return kFALSE;                              //REF multiplicity CUT
  return kTRUE;
}
//________________________________________________________________________________
Bool_t AcceptTrk(const StMuTrack *gTrack = 0) {
  if (! gTrack)            return kFALSE;
  if ( gTrack->nHitsFit(kTpcId) < 15)                            return kFALSE;
  if ( TMath::Abs(gTrack->eta())> 2)          return kFALSE;
  if ( (1.0*gTrack->nHitsFit(kTpcId))/gTrack->nHitsPoss(kTpcId) < 0.51) return kFALSE;
  //if(gTrack->nHitsFit(kPxlId)==0) return kFALSE;        //One pixel hit in the track
  if (  gTrack->flag() < 100 ||  gTrack->flag()%100 == 11) return kFALSE; // bad fit or short track pointing to EEMC
  if (  gTrack->flag() > 1000) return kFALSE;  // pile up track in TPC
  if (  gTrack->nHitsFit() < 10) return kFALSE;
  if (TMath::Abs(gTrack->charge())!=1) return kFALSE;
  if ( gTrack->p().mag()< 0.1)                  return kFALSE;
  // if ( gTrack->p().mag()> 2)                  return kFALSE;
  return kTRUE;
}
//________________________________________________________________________________
Bool_t IsElectron(const StMuMcTrack *mcTrack = 0){
  if (mcTrack->GePid()>3) return kFALSE;
  else if (mcTrack->GePid()==1) return kFALSE;     
  return kTRUE;
}
//________________________________________________________________________________
//--->START MAIN PROGRAM
//________________________________________________________________________________
void RecoSim(Long64_t nevent = 999999,const char* file="./*.MuDst.root",const  char* outFile="test1") {
  
  gROOT->cd();
  
  TString OutFile(outFile);
  OutFile += ".root";
  TFile fOut(OutFile,"recreate");         //Create the file to save the data
  //fOut->cd();
  
  //////////////////////////////////////////////////////->Define the variables 
  DefineTree();      //From KFPartTree.h
  DefineHisto();     //From the AQHistKFPart.h. It define the used histograms
  ///////////////////////////////////////////////////////<-
  
  
  // ----------------------------------------------
  StMuDebug::setLevel(0);  
  maker = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1e9);   // set up maker in read mode
  //                       0,0                        this mean read mode
  //                           dir                    read all files in this directory
  //                               file               bla.lis real all file in this list, if (file!="") dir is ignored
  //                                    filter        apply filter to filenames, multiple filters are separated by ':'
  //                                          10      maximum number of file to read
  maker->SetStatus("*",0);
  const Char_t *ActiveBranches[] = {
    "MuEvent",
    "PrimaryVertices",
    "PrimaryTracks",
    //"CovPrimTrack",
    "GlobalTracks",
    "CovGlobTrack",
    "StStMuMcVertex",
    "StStMuMcTrack",
  }; 

  Int_t Nb = sizeof(ActiveBranches)/sizeof(Char_t *);
  for (Int_t i = 0; i < Nb; i++) maker->SetStatus(ActiveBranches[i],1); // Set Active braches
  StMuDebug::setLevel(0);  
  TChain *tree = maker->chain();
  Long64_t nentries = tree->GetEntries();
  nevent = TMath::Min(nevent,nentries);
  cout << nentries << " events in chain " << nevent << " will be read." << endl;
  //  if (nentries < 100) return;
  tree->SetCacheSize(-1);        //by setting the read cache to -1 we set it to the AutoFlush value when writing
  tree->SetCacheLearnEntries(1); //one entry is sufficient to learn
  tree->SetCacheEntryRange(0,nevent);

  TDatime now;                                          //Set time in Root
  now.Print();
  Int_t count=1;
  for (Long64_t ev = 0; ev < nevent; ev++) {
    if (maker->Make()) break;
    StMuDst* mu = maker->muDst();   // get a pointer to the StMuDst class, the class that points to all the data
    StMuEvent* muEvent = mu->event(); // get a pointer to the class holding event-wise information
    if (ev%100==0){ cout << "Event Number: " << ev << endl; TDatime now2; now2.Print();}
    if (_debugAsk) cout << "Read event #" << ev << "\tRun\t" << muEvent->runId() << "\tId: " << muEvent->eventId() << endl;
    TClonesArray *PrimaryVertices   = mu->primaryVertices(); 
    Int_t NoPrimaryVertices = PrimaryVertices->GetEntriesFast();  // cout << "\tPrimaryVertices " << NoPrimaryVertices;
    TClonesArray *PrimaryTracks    = mu->array(muPrimary);  
    Int_t NoPrimaryTracks = PrimaryTracks->GetEntriesFast();  // cout << "\tPrimaryTracks " << NoPrimaryTracks;
    //TClonesArray *CovPrimTrack     = mu->covPrimTrack(); // cout << "\tCovPrimTrack " << CovPrimTrack->GetEntriesFast();
    TClonesArray *GlobalTracks    = mu->array(muGlobal);  
    Int_t NoGlobalTracks = GlobalTracks->GetEntriesFast();  // cout << "\tPrimaryTracks " << NoPrimaryTracks;
    TClonesArray *CovGlobTrack     = mu->covGlobTrack();
    TClonesArray *MuMcVertices   = mu->mcArray(0); 
    Int_t NoMuMcVertices = MuMcVertices->GetEntriesFast();
    TClonesArray *MuMcTracks     = mu->mcArray(1); 
    Int_t NoMuMcTracks = MuMcTracks->GetEntriesFast(); // cout << "\t" << StMuArrays::mcArrayTypes[1] << " " << NoMuMcTracks;
    //if (_debug) cout << endl;
    if (! NoMuMcVertices || ! NoMuMcTracks) {
      cout << "Ev. " << ev << " has no MC information ==> skip it" << endl;
      continue;
      }
	
    ////////////////////////--->Look for the maximum multiplicity vertex
    // This is the primary vertex but sometimes is bad ranked 
    Int_t MaxMult = 0;
    Int_t nummult = 0;
    for (Int_t ll = 0; ll < NoPrimaryVertices; ll++) {
      StMuPrimaryVertex *Vtx = (StMuPrimaryVertex *) PrimaryVertices->UncheckedAt(ll);
      Float_t Multiplicity = Vtx->nTracksUsed();
      //Float_t Multiplicity = Vtx->noTracks();
      if(MaxMult < Multiplicity) {            //Amilkar: check if the multiplicity is higher than previous
	MaxMult = Multiplicity;               //Amilkar: asing the new maximum value
	nummult = ll;                                   
      }
    }
    ////////////////////////<--- END Look for the maximum multiplicity vertex 


    ///////////////////////---> Fill the primary vertices
    StMuMcVertex *mcVertex1 = (StMuMcVertex *) MuMcVertices->UncheckedAt(0); //The MC info for the primary vertex
    StMuPrimaryVertex *Vtx1 = (StMuPrimaryVertex *) PrimaryVertices->UncheckedAt(nummult);
    if (! Vtx1) continue;
    //if (!AcceptVtx(Vtx1)) continue;
    //if (NoPrimaryTracks<10) continue;    //Only events with more than 10 tracks
    
    TVector3 xyzP(Vtx1->position().xyz());
    //TVector3 *MCxyzP = new TVector3(mcVertex1->XyzV());
    primVtx.eventp    = muEvent->eventId();
    primVtx.l         = nummult;
    primVtx.MultP = Vtx1->noTracks();
    primVtx.refMultP = Vtx1->refMult();
    primVtx.primX = xyzP.X();
    primVtx.primY = xyzP.Y();
    primVtx.primZ = xyzP.Z();
    primVtx.MCpvX = mcVertex1->XyzV().x();
    primVtx.MCpvY = mcVertex1->XyzV().y();
    primVtx.MCpvZ = mcVertex1->XyzV().z();
    primaryvtx->Fill();
    //if (TMath::Abs(primVtx.primZ-primVtx.MCpvZ)>3) continue;
    ////////////->Fill  primary vertices histos
    vertexall->Fill(xyzP.X(),xyzP.Y());
    vertexZall->Fill(xyzP.Z());
    multiplicityP->Fill(Vtx1->noTracks());
    refmultP->Fill(Vtx1->refMult());
    /////////////////////////////////////

    ////////////--->Primary vertex tracks
    for (Int_t k = 0; k < NoPrimaryTracks; k++) {
      StMuTrack *Trk = (StMuTrack *) PrimaryTracks->UncheckedAt(k);  
      if (Trk->vertexIndex() != nummult) continue;
      if (! AcceptTrk(Trk)) continue;
      if (Trk->idTruth()>= NoMuMcTracks) continue;    //This will skip the tracks with no MC partner
      StMuMcTrack *mcTrack = (StMuMcTrack*) MuMcTracks->UncheckedAt(Trk->idTruth()-1);
      //if (!mcTrack) continue;
      //if(!IsElectron(mcTrack)) continue;
      if ( mcTrack->Id() != Trk->idTruth()) {
	cout << "Mismatched idTruth " << Trk->idTruth() << " and mcTrack Id " 
	     <<  mcTrack->Id() << " The track is ignored" <<  endl;
      }
     
      ///////->Calculate DCA
      Int_t kg = Trk->index2Global();    //Like Spiros looping over PrimarTracks 
      if (kg < 0 || kg > NoGlobalTracks) continue; 
      StMuTrack *gTrk = (StMuTrack *) GlobalTracks->UncheckedAt(kg);  // NoPut a -1 shift 
      Int_t kgc = gTrk->index2Cov();         //Like Yuri   looping over GlobalTracks *****Yuri
      if (kgc<=0) continue;
      StDcaGeometry *dcaG = (StDcaGeometry *) CovGlobTrack->UncheckedAt(kgc);   //No Put a -1 shift
      if (! dcaG) continue;
      THelixTrack thelix =  dcaG->thelix();
      Double_t ermx[3];
      Double_t pars[2];
      Double_t vtx[3] = {Vtx1->position().x(), Vtx1->position().y(), Vtx1->position().z()};
      thelix.Move(thelix.Path(vtx));
      thelix.Dca(vtx,pars[0],pars[1],ermx,2);
      Track.dcaxy       = pars[0];
      Track.dcaz        = pars[1];
      ////////////////<- END Calculate DCA
 
      //Track info
      Track.id    = Trk->id();
      Track.GPid   =  mcTrack->GePid();
      Track.MCpxl  =  mcTrack->No_pix_hit();
      Track.dEdx  = Trk->dEdx();
      //Track.beta  = Trk->btofPidTraits().beta();
      Track.pMag  = Trk->momentum().mag(); //AQ
      Track.pt    = Trk->pt(); //AQ
      Track.MCpt  = mcTrack->pT();
      Track.eta   = Trk->eta();
      Track.phi   = Trk->phi();//TMath::ATan(p.y()/p.x());
      Track.charge = Trk->charge();
      Track.tpc   = Trk->nHitsFit(kTpcId);
      Track.istfit   = Trk->nHitsFit(kIstId);
      Track.istpos  = Trk->nHitsPoss(kIstId);
      Track.isttop   = Trk->topologyMap().hasHitInIstLayer(1);//Trk->nHitsFit(kIstId);
      Track.pxlfit   = Trk->nHitsFit(kPxlId);
      Track.pxlpos  = Trk->nHitsPoss(kPxlId);
      Track.pxltop   = Trk->topologyMap().hasHitInPxlLayer(1) + 10*(Trk->topologyMap().hasHitInPxlLayer(2));      

      Track.SigPi  = fabs(Trk->nSigmaPion())>100. ? 100. : Trk->nSigmaPion();
      Track.SigKa  = fabs(Trk->nSigmaKaon())>100. ? 100. : Trk->nSigmaKaon();
      Track.SigPr  = fabs(Trk->nSigmaProton())>100. ? 100. : Trk->nSigmaProton();
      Track.SigEl  = fabs(Trk->nSigmaElectron())>100. ? 100. : Trk->nSigmaElectron();
      
      //primtrack->Fill();
      
      ////Fill Tracks histograms
      pthist->Fill(Track.pt);
      etaphihist->Fill(Track.phi,Track.eta);
      dEdxphist->Fill(Track.pMag,Track.dEdx*1000000);
      //if (Track.beta>-999){
      //Float_t m2tof = Track.pMag*Track.pMag*(1-(Track.beta*Track.beta))/(Track.beta*Track.beta);
      //	m2tofhist->Fill(Track.pMag,m2tof);
      // }
      istfithist->Fill(gTrk->nHitsFit(kIstId));
      istposshist->Fill(gTrk->nHitsPoss(kIstId));
      pxlfithist->Fill(gTrk->nHitsFit(kPxlId));
      pxlposshist->Fill(gTrk->nHitsPoss(kPxlId));
      if(Track.pt>0.2 && TMath::Abs(gTrk->eta())<1){
	istfitcuthist->Fill(gTrk->nHitsFit(kIstId));
	istposscuthist->Fill(gTrk->nHitsPoss(kIstId));
	pxlfitcuthist->Fill(gTrk->nHitsFit(kPxlId));
	pxlposscuthist->Fill(gTrk->nHitsPoss(kPxlId));
      }
      
      // if(Track.pt>0.1 && TMath::Abs(Track.eta)<0.5 && Track.isttop==1 && Track.pxltop==11){
      //dcaxypthist->Fill(1./Track.pMag,Track.dcaxy);
      //dcazpthist->Fill(1./Track.pMag,Track.dcaz);
      // }
      /////////////END FillTracks histograms
  
      if (k==NoPrimaryTracks-1) continue;                     //Does not count the last track
      
      /////////////////////////////////////////
      
      for (Int_t kk = k+1; kk < NoPrimaryTracks; kk++) { //Second track
	StMuTrack *Trkp2 = (StMuTrack *) PrimaryTracks->UncheckedAt(kk); 
	if (Trkp2->vertexIndex() != nummult) continue;
	if (! AcceptTrk(Trkp2)) continue;
	if (Trkp2->idTruth()>= NoMuMcTracks) continue;    //This will skip the tracks with no MC partner
	StMuMcTrack *mcTrack2 = (StMuMcTrack*) MuMcTracks->UncheckedAt(Trkp2->idTruth()-1);
	//if(!IsElectron(mcTrack2)) continue;
	if ( mcTrack2->Id() != Trkp2->idTruth()) {
	  cout << "Mismatched idTruth " << Trkp2->idTruth() << " and mcTrack Id " <<  mcTrack2->Id() 
	       << " The track is ignored" <<  endl;
	}
	
	///////->Calculate DCA
	Int_t kg2 = Trkp2->index2Global();    //Like Spiros looping over PrimarTracks 
	if (kg2 < 0 || kg2 > NoGlobalTracks) continue; 
	StMuTrack *gTrk2 = (StMuTrack *) GlobalTracks->UncheckedAt(kg2);
	Int_t kgc2 = gTrk2->index2Cov();         //Like Yuri   looping over GlobalTracks *****Yuri
	if (kgc2<=0) continue;
	StDcaGeometry *dcaG2 = (StDcaGeometry *) CovGlobTrack->UncheckedAt(kgc2);
	if (! dcaG2) continue;
	////////////////<- END Calculate DCA 
	
	//Apply KFParticles
	AppKFPart(muEvent->magneticField(),Vtx1,dcaG,gTrk,dcaG2,gTrk2,secKF);
	//END Apply KFParticles
	//Apply Helix swing
	AppHelix(muEvent->magneticField(),Vtx1,gTrk,gTrk2,secHE);
	//END Apply Helix swing

	//Fill Histograms
	disthist->Fill(B2.dist);
	projhist->Fill(B2.proj);
	if(B2.q==0 && B2.dist>2 && B2.dist<15 && B2.proj>0.99 && B2.dca2vtx<0.5 && B2.dcamag<0.5 && B2.lambda<0.08) masshist->Fill(B2.Mass);
	if(B2.q!=0 && B2.dist>2 && B2.dist<15 && B2.proj>0.99 && B2.dca2vtx<0.5 && B2.dcamag<0.5 && B2.lambda<0.08) masRhist->Fill(B2.Mass);//masRhist->Fill(B2.MassRot);
	lifehist->Fill(B2.lifet);
	//END Fill Histograms
	 
      }  ///////END TRACK 2
    }  ///////END TRACK 1
    
    
    // }//END VERTEX
    if(_debugAsk) {cout << Form("[%i]",nummult) <<Form(" %8.3f  %8.3f  %8.3f",Vtx1->position().x(),Vtx1->position().y(),Vtx1->position().z()) << "   MCz: " << primVtx.MCpvZ << Form("  Rank:%1.0f",Vtx1->ranking()) << "    Multiplicity: " << Vtx1->noTracks();;
      if(nummult!=0) cout << "   Bad Ranked   " << count++ << endl;
      else cout << endl;}
    
    if (! gROOT->IsBatch()) {
      if (Ask()) return;
    } else {_debugAsk = 0;}
    
  }     //END EVENTS
  
  fOut.Write();
  
  TDatime now1;
  now1.Print();
  
  //fOut->Close();      //Need to be used for condor
}


