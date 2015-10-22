#include <limits>

#include "StNpeCuts.h"

ClassImp(StNpeCuts)

// _________________________________________________________
StNpeCuts::StNpeCuts() : StPicoCutsBase("NpeCutsBase"),
mSecondaryPairDcaDaughtersMax(std::numeric_limits<float>::max()),
mSecondaryPairDecayLengthMin(std::numeric_limits<float>::min()), mSecondaryPairDecayLengthMax(std::numeric_limits<float>::max()),
mSecondaryPairCosThetaMin(std::numeric_limits<float>::min()),
mSecondaryPairMassMin(std::numeric_limits<float>::min()), mSecondaryPairMassMax(std::numeric_limits<float>::max())  {

    // -- default constructor
}

// _________________________________________________________
StNpeCuts::StNpeCuts(const Char_t *name) : StPicoCutsBase(name),
mSecondaryPairDcaDaughtersMax(std::numeric_limits<float>::max()),
mSecondaryPairDecayLengthMin(std::numeric_limits<float>::min()), mSecondaryPairDecayLengthMax(std::numeric_limits<float>::max()),
mSecondaryPairCosThetaMin(std::numeric_limits<float>::min()),
mSecondaryPairMassMin(std::numeric_limits<float>::min()), mSecondaryPairMassMax(std::numeric_limits<float>::max()) {

    // -- constructor
}

// _________________________________________________________
StNpeCuts::~StNpeCuts() {
    // destructor
    
}
