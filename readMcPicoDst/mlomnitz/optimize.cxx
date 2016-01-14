Double_t *Gradient(Double_t p[6],Double_t step[6]){
  Double_t D[6]={0};
  for(Int_t entry=0; entry<6; entry++){
    Double_t pplus[6],pminus[6];
    for(int ii=0; ii<6;ii++){
      if(ii==entry){
	pplus[ii]=p[ii]+step[ii];
	pminus[ii]=p[ii]-step[ii];
      }
      else{
	pplus[ii]=p[ii];
	pminus[ii]=p[ii];
      }
    }
    D[entry]=(significance(pplus)-significance(pminus))/(2.0*step[entry]);
  }
  Double_t *toRet=D;
  return(toRet);
}
void optimize(){

  return;
}
