#ifndef ApplyHelix_hh
#define ApplyHelix_hh

void AppHelix(Double_t Bfield, StMuPrimaryVertex* Vertex, StMuTrack *t_pos, StMuTrack *t_neg,TTree *secondaryvtx){
  //t_pos is the first track
  //t_neg is the second track not positive of negative
  //if(t_pos->momentum().perp()<1 && t_pos->dca().mag()<0.75) return;
  //if(t_neg->momentum().perp()<1 && t_neg->dca().mag()<0.75) return;
  //Original
  StPhysicalHelixD helix_pos(t_pos->helix().momentum(Bfield*kilogauss), t_pos->helix().origin(), Bfield*kilogauss, t_pos->charge());
  StPhysicalHelixD helix_neg(t_neg->helix().momentum(Bfield*kilogauss), t_neg->helix().origin(), Bfield*kilogauss, t_neg->charge());
  //cout<<"momentum first = "<<t_pos->helix().momentum(Bfield*kilogauss).mag()<<" momentum last = "<<t_pos->outerHelix().momentum(Bfield*kilogauss).mag()<<endl;
  Double_t res = 3.0;  // allow resolution
  Double_t xc_pos = helix_pos.xcenter();
  Double_t yc_pos = helix_pos.ycenter();
  Double_t xc_neg = helix_neg.xcenter();
  Double_t yc_neg = helix_neg.ycenter();
  Double_t dd = TMath::Sqrt((xc_pos-xc_neg)*(xc_pos-xc_neg)+(yc_pos-yc_neg)*(yc_pos-yc_neg));
  Double_t r_pos = 1./helix_pos.curvature();
  Double_t r_neg = 1./helix_neg.curvature();
  //Switched following && for  ||
  //if( fabs(r_pos-r_neg)-res > dd || dd > r_pos+r_neg+res ) return;
  B2.trk1dca=helix_pos.distance(Vertex->position());
  B2.trk2dca=helix_neg.distance(Vertex->position());
  //TPair *s = new TPair();  
  pair<double,double> s = helix_pos.pathLengths(helix_neg);
  //Original
  StThreeVectorF dcaP_pos = helix_pos.at(s.first);
  StThreeVectorF dcaP_neg = helix_neg.at(s.second);
  Float_t dcaDaughters = (dcaP_pos - dcaP_neg).mag();
  //if(dcaDaughters>V0Finder::mV0DcaDaughtersMax) continue;
  if(dcaDaughters>1.0) return;
  StThreeVectorF primaryVertex = Vertex->position().xyz();
  StThreeVectorF v0 = (dcaP_pos+dcaP_neg)*0.5;
  /////////////////////////////////////
  cout<<"Helix swimming"<<endl;
  cout<<"Passed! Will continue all calc."<<endl;
  cout<<"Track 1 origin"<<t_pos->dca().x()<<" "<<t_pos->dca().y()<<" "<<t_pos->dca().z()<<endl;
  cout<<"Track 2 origin"<<t_neg->dca().x()<<" "<<t_neg->dca().y()<<" "<<t_neg->dca().z()<<endl;
  cout<<"trk1 dca 2 trk2: ("<<dcaP_pos.x()<<" ,"<<dcaP_pos.y()<<" ,"<<dcaP_pos.z()<<")"<<endl;
  cout<<"trk2 dca 2 trk1: ("<<dcaP_neg.x()<<" ,"<<dcaP_neg.y()<<" ,"<<dcaP_neg.z()<<")"<<endl;
  cout<<"RC vertex"<<v0.x()<<" "<<v0.y()<<" "<<v0.z()<<endl;
  cout<<"Original vertex"<<Vertex->position().x()<<" "<<Vertex->position().y()<<" "<<Vertex->position().z()<<endl;
  cout<<"Decay length"<<(v0-primaryVertex).mag()<<endl;
  ////////////////////////////////////
  StThreeVectorF mom_pos = helix_pos.momentumAt(s.first, Bfield*kilogauss);
  StThreeVectorF mom_neg = helix_neg.momentumAt(s.second, Bfield*kilogauss);
  StThreeVectorF mom_v0 = mom_pos + mom_neg;
  Float_t angle = (v0-primaryVertex).angle(mom_v0);
  Float_t dca2vtx = (v0-primaryVertex).mag()*TMath::Sin(angle);
  Float_t decaylength = (v0-primaryVertex).mag();
  
  //TVector3 p1(KFtrack1.GetPx(),KFtrack1.GetPy(),KFtrack1.GetPz());
  //TVector3 p2(KFtrack2.GetPx(),KFtrack2.GetPy(),KFtrack2.GetPz());
  TVector3 p2rot((-1)*mom_neg.x(),(-1)*mom_neg.y(),(-1)*mom_neg.z());
  //TVector3 p2rot((-1)*mom_neg.x(),(-1)*mom_neg.y(),mom_neg.z());  //This one does not work
  TLorentzVector Lp1,Lp2,Lp2rot;
  //Float_t masselectron = 0.00051;
  Float_t masspion = 0.13957;
  Float_t masskaon=0.49367;
  Lp1.SetXYZM(mom_pos.x(),mom_pos.y(),mom_pos.z(),masspion);
  Lp2.SetXYZM(mom_neg.x(),mom_neg.y(),mom_neg.z(),masskaon);
  Lp2rot.SetVectM(p2rot,masskaon);
  //Lp1.SetXYZM(mom_pos.x(),mom_pos.y(),mom_pos.z(),TMath::Sqrt(mom_pos.mag2()+masspion*masspion));   //Like the StHelixV0.cxx line 113
  //Lp2.SetXYZM(mom_neg.x(),mom_neg.y(),mom_neg.z(),TMath::Sqrt(mom_neg.mag2()+masspion*masspion));
  //Lp2rot.SetVectM(p2rot,sqrt(p2rot.Mag2()+masspion*masspion));
  //Fill secondary vertex tree
  // B2.eventt = muEvent->eventId();
  B2.q      = t_pos->charge() + t_neg->charge();
  B2.dist   = decaylength;
  B2.HEX    = v0.x();
  B2.HEY    = v0.y();
  B2.HEZ    = v0.z();
  B2.lambda   = angle;
  B2.phi    = mom_v0.phi();
  B2.proj = TMath::Cos(angle);
  //cout << B2.proj << endl;
  B2.mom    = mom_v0.mag();
  B2.dca2vtx = dca2vtx;
  B2.lifet  = dcaDaughters;
  B2.Mass   = (Lp1+Lp2).M();
  B2.MassRot = (Lp1+Lp2rot).M();
  
  B2.dcaxy = (dcaP_pos - dcaP_neg).x();
  B2.dcaz  = (dcaP_pos - dcaP_neg).z();
  B2.dcamag = dcaDaughters;
  //cout << B2.dcaxy << "   " << B2.dcaz << endl;
  if (B2.Mass<5.0 && B2.Mass>0.5) secondaryvtx->Fill();
  /*StThreeVectorF origin(v0.x(),v0.y(),v0.z());
  StPhysicalHelix toRet(mom_v0,origin, Bfield*kilogauss,0);
  return(toRet);*/
  return;
}
#endif
