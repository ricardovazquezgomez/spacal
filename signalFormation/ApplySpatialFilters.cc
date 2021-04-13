#include "ApplySpatialFilters.hh"
#include <iostream>
bool ApplySpatialFilters(const float& UnifSpatEff, TGraph2D *SEGraph, float &x, float &y, float &z, TRandom *randGen)
{



    if (SEGraph)                                                                    // If the efficiency graph exists
    {                       
        if ( randGen -> Rndm() <  SEGraph -> Interpolate(x,y)) return true;         // If the random number is below the efficiency value, return true.
        else return false;                                                          // else, return false
    }



    // ****************************************************************
    // TEMPORARY SOLUTION
    if ( randGen -> Rndm() < UnifSpatEff) return true;
    else return false;
    // ****************************************************************


    // if ( randGen -> Rndm() <  0.6) return true;         // If the random number is below the efficiency value, return true.
    // else return false;                                                          // else, return false
    // else return true;                                                              // else, return true.
}
