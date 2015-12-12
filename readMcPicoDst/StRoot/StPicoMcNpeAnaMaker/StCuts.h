#ifndef NPE_EVENT_CUTS_H
#define NPE_EVENT_CUTS_H
/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  **Kunsu OH        (kunsuoh@gmail.com)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>

namespace cuts
{

    //event
    float const vz = 30.0;
    float const vzVpdVz = 3.0;
    
    // track
    int const nHitsFit = 20;
    float const nHitsRatioMin = 0.52; float const nHitsRatioMax = 1.2;
    float const ptMin = .2;
    float const ptMax = 5.;

    
    // electrons
    int const nHitsDedx = 15;
    float const globalDca = 1.;
    float const etaTagged = 0.5;
    
    // partner
    float const etaPartner = 0.7;
    
    // electron pair cuts
    float const pairMass = 0.5;
    float const pairDca = 3.;
    
}
#endif
