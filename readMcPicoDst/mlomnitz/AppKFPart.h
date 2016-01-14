#ifndef ApplyKFPart_hh
#define ApplyKFPart_hh

void AppKFPart(Double_t Mag, StMuPrimaryVertex* Vertex, StDcaGeometry *dca1, StMuTrack *Trk1, StDcaGeometry *dca2,StMuTrack *Trk2, TTree* secondaryvtx){
  //TEst on daighter dca
  StPhysicalHelixD helix_pos(Trk1->helix().momentum(Mag*kilogauss), Trk1->helix().origin(), Mag*kilogauss, Trk1->charge());
  StPhysicalHelixD helix_neg(Trk2->helix().momentum(Mag*kilogauss), Trk2->helix().origin(), Mag*kilogauss, Trk2->charge());
  pair<double,double> s = helix_pos.pathLengths(helix_neg);
  StThreeVectorF dcaP_pos = helix_pos.at(s.first);
  StThreeVectorF dcaP_neg = helix_neg.at(s.second);
  Float_t dcaDaughters = (dcaP_pos - dcaP_neg).mag();
  //if(dcaDaughters>V0Finder::mV0DcaDaughtersMax) continue;
  //if(helix_pos.distance(Vertex->position())<1.5 || helix_neg.distance(Vertex->position())<1.5) return;
  if(dcaDaughters>1.0) return;
  B.dcadaughter=dcaDaughters;
  KFParticle::SetField(Mag);    //Need Magnetic field for KFPArticles
  ////Set vertex
  MVertex vert;
  vert.SetXYZ(Vertex->position().x(),Vertex->position().y(),Vertex->position().z());
  Float_t xvrtxres      = Vertex->posError().x()*Vertex->posError().x();
  Float_t yvrtxres      = Vertex->posError().y()*Vertex->posError().y();
  Float_t zvrtxres      = Vertex->posError().z()*Vertex->posError().z();
  Double_t  CovVert[6] = {              
    xvrtxres,                       //Why the covariance matrix is this one???????
    0, yvrtxres,                   // The matrix is 3x3 but it takes only the lower part
    0, 0, zvrtxres};         //Should try switch 0 -> xvrtxrex*yvrtxres,xvrtxrex*zvrtxres,yvrtxrex*zvrtxres    
  if (CovVert[0] < 1e-8) CovVert[0] = 1e-8;
  if (CovVert[2] < 1e-8) CovVert[2] = 1e-8;
  if (CovVert[5] < 1e-8) CovVert[5] = 1e-8;
  vert.SetCovarianceMatrix(CovVert);
  vert.SetNContributors(2);// 2 contributors (=2 tracks)
  vert.SetChi2(1.01);
  //cout << " vertex is X : " << vert.GetX() << " X : "<<  vert.GetY()<<"  Z : " << vert.GetZ() << endl;
  KFParticle KFVertex;
  KFVertex = KFParticle(vert);
  ////END set vertex

  Int_t pid = 211;   //PID is required 

  //Track 1
  MTrack     Mtrack1;
  TRVector posmom1(6);   //position and momentum
  TRSymMatrix C1(21);
  dca1->GetXYZ(posmom1.GetArray(),C1.GetArray());
  Mtrack1.SetParameters(posmom1.GetArray());
  Mtrack1.SetCovarianceMatrix(C1.GetArray());
  Mtrack1.SetNDF(1);
  Mtrack1.SetChi2(1.5);
  Mtrack1.SetCharge(Trk1->charge());
  if (abs(Trk1->charge())!=1) return;//continue;
  //Track one is pi+
  pid=211;
  //  if (Trk1->charge()==-1) pid=-211;
  //else pid=211;
  KFParticle KFtrack1;
  KFtrack1 = KFParticle(Mtrack1,pid);     
  //END Track 1
  
  Int_t pid2 = 211;   //PID is required
  //Track 2
  MTrack     Mtrack2;
  TRVector posmom2(6);    //position and momentum
  TRSymMatrix C2(21);
  dca2->GetXYZ(posmom2.GetArray(),C2.GetArray());
  Mtrack2.SetParameters(posmom2.GetArray());
  Mtrack2.SetCovarianceMatrix(C2.GetArray());
  Mtrack2.SetNDF(1);
  Mtrack2.SetChi2(1.5);
  Mtrack2.SetCharge(Trk2->charge());
  if (abs(Trk2->charge())!=1) return;//continue;
  //if (Trk2->charge()==-1) pid2=-211;
  //else pid2=211;
  pid2=321;
  KFParticle KFtrack2;
  KFtrack2 = KFParticle(Mtrack2,pid2);    //Default Pid is pion
  //END Track 2
  
  //START KFParticle
  TVector3 p1(KFtrack1.GetPx(),KFtrack1.GetPy(),KFtrack1.GetPz());
  TVector3 p2(KFtrack2.GetPx(),KFtrack2.GetPy(),KFtrack2.GetPz());
  TVector3 p2rot((-1)*KFtrack2.GetPx(),(-1)*KFtrack2.GetPy(),(-1)*KFtrack2.GetPz());
  //TVector3 p1U =p1.Unit();
  //TVector3 p2U = p2.Unit();
  //cout << p1U.X() << "  " << p1U.Y() << "  " << p1U.Z() << endl;
  //cout << p2U.X() << "  " << p2U.Y() << "  " << p2U.Z() << endl; 
  TVector3 p3 = p1+p2;
  TLorentzVector Lp1,Lp2,Lp2rot;
  //Float_t masselectron = 0.00051;
  Float_t masspion = 0.13957;
  Float_t masskaon=0.49367;
  Lp1.SetVectMag(p1,masspion);
  Lp2.SetVectMag(p2,masskaon);
  Lp2rot.SetVectMag(p2rot,masskaon);
/// This is the example of the KFParticleBase::Construct function usage.
/// parameter 1 - array of the daughter particles
/// parameter 2 - number of the daughter particles
/// parameter 3 - vertex (it should be the object of the KFParticle class)
/// parameter 4 - the value we force the particle mass to be equial to.
  const KFParticle pVertex = KFVertex;           //vertex (it should be the object of the KFParticle class)
  const KFParticle *vDaughters[2] = {&KFtrack1,&KFtrack2};
  KFParticle RecoKFPart;
  RecoKFPart.Construct(vDaughters,2,&pVertex,-1,0);
  RecoKFPart.SetProductionVertex(pVertex);
  //END KFParticle
  //cout << RecoKFPart.GetX() << " " <<  RecoKFPart.GetY() << " " << RecoKFPart.GetZ() << endl;
//Fill secondary vertex tree
  // B.eventt = 0;//muEvent->eventId();
  B.q      = Trk1->charge() + Trk2->charge();
  B.dist   = RecoKFPart.GetDecayLength();
  B.KFX   = RecoKFPart.GetX();
  B.KFY   = RecoKFPart.GetY();
  B.KFZ   = RecoKFPart.GetZ();
  TVector3 dist2vtx(RecoKFPart.GetX()-KFVertex.GetX(),RecoKFPart.GetY()-KFVertex.GetY(),RecoKFPart.GetZ()-KFVertex.GetZ());
  TVector3 recopartmom(RecoKFPart.GetPx(),RecoKFPart.GetPy(),RecoKFPart.GetPz());    //Momentum of the KFParticle
  TVector3 xyzD(RecoKFPart.GetX(),RecoKFPart.GetY(),RecoKFPart.GetZ());
  //cout << dist2vtx.Mag() << "  " << << endl;
  TVector3 dist2vtxU = dist2vtx.Unit();
  TVector3 recopartmomU = p3.Unit();//recopartmom.Unit();
  TVector3 D(xyzD.X()-Vertex->position().x(),xyzD.Y()-Vertex->position().y(),xyzD.Z()-Vertex->position().z());
  //cout << dist2vtxU.X() << "  " << dist2vtxU.Y() << "   " << dist2vtxU.Z() << endl;
   //cout << recopartmomU.X() << "  " << recopartmomU.Y() << "   " << recopartmomU.Z() << endl;
  TVector3 U = D.Unit();
  //cout << U.X() << "  " << U.Y() << "   " << U.Z() << endl;
  //cout << recopartmomU.X() << "  " << recopartmomU.Y() << "   " << recopartmomU.Z() << endl;
  //B.lambda   = fabs(Trk1->p().mag() - Trk2->p().mag());//dist2vtx.Angle(recopartmom);
  B.lambda=dist2vtx.Angle(recopartmom);
  B.phi    = RecoKFPart.GetPhi();
  B.proj = U.Dot(recopartmomU);//dist2vtxU.Dot(recopartmomU);//TMath::Cos(B.lambda);
  //cout << B.proj << endl;
  B.mom    = RecoKFPart.GetP();
  B.dca2vtx = dist2vtx.Mag()*TMath::Sin(dist2vtx.Angle(recopartmom));
  B.lifet  = RecoKFPart.GetLifeTime();
  B.Mass   = RecoKFPart.GetMass();
  B.MassRot = (Lp1+Lp2rot).M();
  //cout << "Mass: " << B.Mass << endl;
  /////////////////////////////////////
  cout<<"KF particle"<<endl;
  cout<<"Passed! Will continue all calc."<<endl;
  cout<<"trk1 dca 2 trk2: ("<<dcaP_pos.x()<<" ,"<<dcaP_pos.y()<<" ,"<<dcaP_pos.z()<<")"<<endl;
  cout<<"trk2 dca 2 trk1: ("<<dcaP_neg.x()<<" ,"<<dcaP_neg.y()<<" ,"<<dcaP_neg.z()<<")"<<endl;
  cout << "RC vertex"<<RecoKFPart.GetX() << " " << RecoKFPart.GetY() << " " << RecoKFPart.GetZ() << endl;  
  cout<<"Original vertex"<<Vertex->position().x()<<" "<<Vertex->position().y()<<" "<<Vertex->position().z()<<endl;
  cout<<"Decay length"<<RecoKFPart.GetDecayLength()<<endl;
  /////////////////////////////////////
  THelixTrack   tHelix1, tHelix2; 
  Double_t d1=0,d2=0;
  Double_t x1[3]={0,0,0},x2[3]={0,0,0};
  tHelix1  =  dca1->thelix();
  tHelix2  =  dca2->thelix();
  d1  = tHelix1.PathX(tHelix2,&d2);
  tHelix1.Eval(d1,x1);
  tHelix2.Eval(d2,x2);
  TVector3 dcavpi(x1[0]-x2[0],x1[1]-x2[1],x1[2]-x2[2]);
  B.dcaxy = dcavpi.X();
  B.dcaz  = dcavpi.Z();
  B.dcamag= dcavpi.Mag();
  B.Prob  = TMath::Prob(RecoKFPart.Chi2(),RecoKFPart.GetNDF());
  //if (dcavpi.Mag()>0.85) return;
  //if(TMath::Prob(RecoKFPart.Chi2(),RecoKFPart.GetNDF())<0.1) return;
  if (B.Mass<5 && B.Mass>0.05) secondaryvtx->Fill();
  //END Fill secondary vertex tree
  //StPhysicalHelix toRet(RecoKFPart.GetP(),recopartmom,Mag*kilogauss,B.q);
  return;
}
#endif
