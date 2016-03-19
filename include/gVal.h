#ifndef GVAL_HPP_
#define GVAL_HPP_

#include "app.h"
#include "arch.h"
#include "param.h"


extern struct Param param;

extern vector<pProcess> gProcess;
extern pLink gLink;
extern pTile gTile;

extern int gProcNum;
extern int gTileNum;
extern int gLinkNum;

extern int g_edge_size;

// normalization values, used by the combined objective algorithm;
extern float COST_ENERGY_NORMALIZATION;
extern float COST_RELIABILITY_NORMALIZATION;

#endif
