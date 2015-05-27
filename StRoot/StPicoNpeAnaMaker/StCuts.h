#ifndef CUTS_H
#define CUTS_H
/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  **Kunsu OH        (kunsuoh@gmail.com)
 *              Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>

namespace cutsAna
{
    // path to lists of triggers prescales
    // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/

    // event
    float const vz = 6.0;// cm.
    float const vzVpdVz = 3.0; // 3 cm.
    unsigned char const trigger = 0; // 19: BHT1, 21: BHT2, 23: BHT3
    unsigned char const triggerLength = 0xFFFF;
    
    // hadron selectron
    float const nSigPion = 2;
    bool const pionRequireHFT = true;
    int const pionNHitsFit = 20;
    int const pionNhitsDedx = 15;
    float const pionEta = 0.7;
    float const pionPt = 1.5;
    
    // electron
    bool const electronRequireHFT = true;
    int const electronNHitsFit = 20;
    int const electronNhitsDedx = 15;
    float const electronEta = 0.7;
    float const electronPt = 1.5;
    
    
    // partner
    int const partnerNHitsFit = 15;
    float const partnerEta = 1.;
    float const partnerPt = 0.2;
    
    
    // electron + partner pair cuts
    float const pairMass = 0.4;
    float const pairDca = 3.;
    
    
    // pid
    float const nSigElectron = 3;
    float const partnerNSigElectron = 3;
    
    int const emcNEta = 0;//1;
    int const emcNPhi = 0;//1;
    float const emcEoverPLow = 0;//0.5;
    float const emcEoverPHigh = 10;//1.7;
    float const emcPhiDist = 10;//0.013;
    float const emcZDist = 100;//2;
    float const emcAssDist = 10;//0.05;
    float const tofBeta = 1;//0.1;
}
#endif
