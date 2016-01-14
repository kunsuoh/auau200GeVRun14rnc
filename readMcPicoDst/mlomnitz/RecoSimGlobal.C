// Read the MuDst of a simulation
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
#include "TNtuple.h"
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
    if (abs(Vtx->position().z()) > 6) return kFALSE;
    if (Vtx->refMult()<30) return kFALSE;                              //REF multiplicity CUT
    return kTRUE;
}
//________________________________________________________________________________
Bool_t AcceptTrk(const StMuTrack *gTrack = 0) {
    if (! gTrack)                                                           return kFALSE;
    if ( gTrack->nHitsFit(kTpcId) < 15)                                     return kFALSE;
    if ( TMath::Abs(gTrack->eta())> 1)                                      return kFALSE;
    if ( (1.0*gTrack->nHitsFit(kTpcId))/gTrack->nHitsPoss(kTpcId) < 0.51)   return kFALSE;
    //if(gTrack->nHitsFit(kPxlId)==0) return kFALSE;        //One pixel hit in the track
    if (  gTrack->flag() < 100 ||  gTrack->flag()%100 == 11)                return kFALSE; // bad fit or short track pointing to EEMC
    if (  gTrack->flag() > 1000)                                            return kFALSE;  // pile up track in TPC
    if (  gTrack->nHitsFit() < 10)                                          return kFALSE;
    if (TMath::Abs(gTrack->charge())!=1)                                    return kFALSE;
    if ( gTrack->p().mag()< 0.2)                                            return kFALSE;
    // if ( gTrack->p().mag()> 2)                  return kFALSE;
    return kTRUE;
}
//________________________________________________________________________________
Bool_t IsElectron(const StMuMcTrack *mcTrack = 0){
    if (mcTrack->GePid()>3) return kFALSE;
    else if (mcTrack->GePid()==1) return kFALSE;
    return kTRUE;
}
Bool_t IsPion(const StMuMcTrack *mcTrack = 0){
    if (mcTrack->GePid()>9 || mcTrack->GePid()<7) return kFALSE;
    else return kTRUE;
}
Bool_t IsKaon(const StMuMcTrack *mcTrack = 0){
    if (mcTrack->GePid()>12 || mcTrack->GePid()<11) return kFALSE;
    else return kTRUE;
}
//________________________________________________________________________________
//--->START MAIN PROGRAM
//________________________________________________________________________________
void RecoSimGlobal(Long64_t nevent = 999999,const char* file="./*.MuDst.root",const  char* outFile="test1") {
    
    gROOT->cd();
    
    TString OutFile(outFile);
    OutFile += ".root";
    TFile fOut(OutFile,"recreate");         //Create the file to save the data
    //fOut->cd();
    //MICHAEL L
    Char_t *tempChar=Form("ks_%s.root",outFile);
    //TFile *myoutfile=new TFile(tempChar,"RECREATE");
    //TNtuple *tTree=new TNtuple("tTree","pid:mass","pid:mass");
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
        //cout<<"Here MC informatio. Vertices: "<<NoMuMcVertices<<" tracks "<< NoMuMcTracks<<endl;
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
        if (TMath::Abs(primVtx.primZ-primVtx.MCpvZ)>3) continue;
        //MICHAEL FOR TEST purposes
        std::vector<std::pair<int,int> > decayId;
        for(Int_t k = 0; k < NoMuMcTracks; k++) {
            //      StMuTrack *pTrk = (StMuTrack *) PrimaryTracks->UncheckedAt(k);
            //if (pTrk->idTruth()>= NoMuMcTracks) continue;    //This will skip the tracks with no MC partner
            StMuMcTrack *pmcTrack=(StMuMcTrack*) MuMcTracks->UncheckedAt(k);
            //if(pmcTrack->GePid()!=16) continue;//K0short
            ////////////  Geant PID info //////////
            //
            //  -p+ 8 p- 9
            //  -16 For K0s
            //  -37 Fror D0 https://drupal.star.bnl.gov/STAR/blog/yfzhang/instruction-hft-simulation
            //
            /////////////////////
            //cout<<"Lomnitz: Curent track PID"<<pmcTrack->GePid()<<endl;
            if(pmcTrack->GePid()==1 || pmcTrack->GePid()==10007){ //D0 according to http://www.star.bnl.gov/public/comp/simu/gstar/kumacs/NewParticle.html
                //cout<<"Lomnitz: Adding to list"<<endl;
                decayId.push_back(std::make_pair(pmcTrack->IdVxEnd(),pmcTrack->GePid()));
            }
        }
        
        
        
        /////////////////////////////////////////
        vector<int> piplus,pivertex,kminus,kvertex;
        ////////////--->Primary vertex tracks
        for (Int_t k = 0; k < NoGlobalTracks; k++) {
            StMuTrack *Trk = (StMuTrack *) GlobalTracks->UncheckedAt(k);
            StTrack *stTrk = (StTrack *) mu->createStTrack(Trk);
            
            if (Trk->vertexIndex() != nummult) continue;
            if (! AcceptTrk(Trk)) continue;
            if (Trk->idTruth()>= NoMuMcTracks) continue;    //This will skip the tracks with no MC partner
            StMuMcTrack *mcTrack = (StMuMcTrack*) MuMcTracks->UncheckedAt(Trk->idTruth()-1);
            //if (!mcTrack) continue;
            //if(!IsElectron(mcTrack)) continue;
            //Bool_t isPi=IsPion(mcTrack);
            //MICHAEL L
            if ( mcTrack->Id() != Trk->idTruth()) {
                cout << "Mismatched idTruth " << Trk->idTruth() << " and mcTrack Id "
                <<  mcTrack->Id() << " The track is ignored" <<  endl;
            }
            //if(mcTrack->GePid()!=16) continue; //Lomnitz Only keeping kshort for my tests
            //cout<<"Lomnitz: Gots Kaon"<<endl;
            //kshort decay vertex
            //Int_t decayVtxId=mcTrack->IdVxEnd();
            //Now will fin particles created at the decay vertex
            
            ///////->Calculate DCA
            //Int_t kg = Trk->index2Global();    //Like Spiros looping over PrimarTracks
            //if (kg < 0 || kg > NoGlobalTracks) continue;
            //StMuTrack *gTrk = (StMuTrack *) GlobalTracks->UncheckedAt(kg);  // NoPut a -1 shift
            Int_t kgc = Trk->index2Cov();         //Like Yuri   looping over GlobalTracks *****Yuri
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
            //Finding track parent for pc distribution test
            Track.isGamma=0;
            Track.isPi0Dalitz=0;
            Int_t prodvtx_ID=mcTrack->IdVx();
            if(!IsElectron(mcTrack)) continue;
            cout<<"Checking against "<<decayId.size()<<" decay vertices for daugter id "<<mcTrack->GePid()<<endl;
            for(Int_t jjj=0; jjj<decayId.size(); jjj++){
                pair<int, int> temp=decayId.at(jjj);
                if(temp.first!=mcTrack->IdVx()) continue;
                cout<<"Found parent vertex, parent particle is "<<temp.second<<endl;
                if(temp.second==1){
                    Track.isGamma=1;
                    jjj=decayId.size();
                    continue;
                }
                if(temp.second==10007){
                    Track.isPi0Dalitz=1;
                    jjj=decayId.size();
                    continue;
                }
            }
            primtrack->Fill();
            
            StThreeVectorF pos = Trk->firstPoint();
            //StThreeVectorF pos2 = mcTrack->firstPoint();
            cout << pos.x() << "RC: " << pos.y() << " " << pos.z() << "   idTruth(): " << stTrk->idTruth() << endl;

        }
        
    }     //END EVENTS
    
    fOut.Write();
    //myoutfile->cd();
    //tTree->Write();
    //myoutfile->Close();
    TDatime now1;
    now1.Print();
    
    //fOut->Close();      //Need to be used for condor
}


