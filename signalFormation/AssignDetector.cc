#include "AssignDetector.hh"
    /* 
    The idea is to have an int ID, nDet, which will range from 0 to the number of Detectors (e.g. [0,17]).
    - If nothing is found, then returns the initialisation value of nDet: -1. 
    - The TCutG is assumed to be even!!!    < size/2 --> Front.     > size/2 --> Back
    */


int AssignDetector(float &x, float &y, float &z, int &nDetSim, std::vector<TCutG> *Borders)
{
    if (!Borders) return nDetSim;       // If no TCutG are passed, simply return what the sim is already telling us

    int nDet = -1;                      // If everything goes wrong, just return -1
    if ( z > 0 )                         // z > 0 --> Front detectors
    {
        for( int iDet = 0; iDet < (int) Borders->size() / 2.; iDet ++) if ( Borders->at(iDet).IsInside(x, y))                           { nDet = iDet; break; }  // nDet will range from 0 to nDetector in the front, e.g. [0,8]
    }

    else if (z < 0 )                     // z < 0 --> Back detectors
    {
        for( int iDet = 0; iDet < (int) Borders->size() / 2.; iDet ++) if ( Borders->at(iDet + int(Borders->size()/2.)).IsInside(x, y)) { nDet = iDet + int(Borders->size()/2.); break; }    // nDet will range from 0 to nDetector in the back + nDetector in the front, e.g. [9,17]
    }


    return nDet;              
}
