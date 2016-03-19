#ifndef SA_HPP_
#define SA_HPP_

#include "gVal.h"
#include "arch.h"
#include "nocmap.h"
#include <cmath>

#define A   147453245
#define C   226908347
#define M   1073741824   
#define TOLERANCE 1
#define TEMPS 5
#define MINACCEPT 0.001

double uniform_rv();
void init_rand (int seed);
void init_rand (int seed);
long int uniform_int_rv(long int imin, long int imax);

int accept(double deltac, double temperature);
double anneal_at_temp(double t, double& currentcost, double *acceptratio);
void anneal();

void makeRandomSwap(int & tile1, int & tile2);

void rollBack(int tile1, int tile2);

void swapProcesses(int tile1, int tile2);

float calc_total_cost();
float calc_total_cost_with_reliability();

float fixed_routing_overload_calc();
float adaptive_routing_overload_calc();
bool sa_route_traffic(int src_tile, int dst_tile, int BW);

bool sa_program_routers();

#endif
