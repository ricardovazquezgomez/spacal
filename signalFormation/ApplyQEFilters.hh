#ifndef APPLYQEFILTERS_HH
#define APPLYQEFILTERS_HH
#include "TGraph.h"
#include "TRandom3.h"

bool ApplyQEFilters(const TGraph &QEGraph, float &wave, TRandom *randGen);

#endif