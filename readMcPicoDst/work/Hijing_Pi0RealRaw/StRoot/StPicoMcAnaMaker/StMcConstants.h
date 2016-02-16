#ifndef StMcConstants_H
#define StMcConstants_H
#include "Rtypes.h"
#include <string>

#include "StMcConstants.h"

namespace mcConst
{
  const Int_t phibins=300;
  const Int_t pTbins=200;
  const int nPtBins=18;
  const int dcaBins=100;
  float ptBins[nPtBins+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25,
			     0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 
			     0.7, 0.8, 0.9, 1., 1.2, 1.4, 2.};
}
#endif
