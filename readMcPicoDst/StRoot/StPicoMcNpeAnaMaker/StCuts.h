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
    float const vz = 6.0;
    float const vzVpdVz = 3.0;
    
    // track
    int const nHitsFit = 15;
    float const ptMin = .6;
    float const ptMax = 20.;
    float const eta = 1.;
    
}
#endif
