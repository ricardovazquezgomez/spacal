#ifndef ASSIGNDETECTOR_HH
#define ASSIGNDETECTOR_HH

#include "TCutG.h"
#include <vector>
/* -1 --> Out of every detector. */
int AssignDetector(float &x, float &y, float &z, int &nDetSim, std::vector<TCutG> *Borders);

#endif