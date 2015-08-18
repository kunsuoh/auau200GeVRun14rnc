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
    std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";
    
    // event
    float const vz = 6.0;// cm.
    float const vzVpdVz = 3.0; // 3 cm.
    unsigned char const trigger = 0; // 19: BHT1, 21: BHT2, 23: BHT3
    unsigned int const triggerLength = 0xFFFF;
    
    // track
    bool const trackRequireHFT = false;
    int const trackNHitsFit = 20;
    int const trackNhitsDedx = 15;
    float const trackEta = 0.7;
    float const trackPt = 1.5;
    float const trackDca = 0.1;
    
    
    
    // partner
    int const partnerNHitsFit = 15;
    float const partnerEta = 1.;
    float const partnerPt = 0.2;
    
    
    // electron + partner pair cuts
    bool const isRecoPhE = true;
    float const pairMass = 0.1;
    float const pairDca = 1.5; // 0.5 for HFT on
    
    
    // pid
    float const taggedNSigElectron = 3;
    float const partnerNSigElectronLow = 0;
    float const partnerNSigElectronHigh = 3;
    
    int const emcNEta = 1;
    int const emcNPhi = 1;
    float const emcEoverPLow = 0.8;
    float const emcEoverPHigh = 2.0;
    float const emcEoverPLow2 = 0.333;
    float const emcEoverPHigh2 = 3;
    float const emcPhiDist = 0.015;//0.013;
    float const emcZDist = 3;//2;
    float const emcAssDist = 0.06;//0.05;
    float const tofBeta = 0.025;
}
#endif
