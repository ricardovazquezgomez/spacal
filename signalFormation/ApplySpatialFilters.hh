#ifndef APPLYSPATIALFILTERS_HH
#define APPLYSPATIALFILTERS_HH
#include "TGraph2D.h"
#include "TRandom3.h"

bool ApplySpatialFilters(const float& UnifSpatEff, TGraph2D *SEGraph, float &x, float &y, float &z, TRandom *randGen);

#endif