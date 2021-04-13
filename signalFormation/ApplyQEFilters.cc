#include "ApplyQEFilters.hh"

bool ApplyQEFilters(const TGraph &QEGraph, float &wave, TRandom *randGen)
{

    if ( randGen -> Rndm() <  QEGraph.Eval(wave)) return true;          // If the random number is below the efficiency value, return true.
    else return false;                                                  // else, return false
    
}
