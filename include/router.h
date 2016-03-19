#ifndef ROUTER_HPP_
#define ROUTER_HPP_

#include <iostream>
#include "app.h"
#include "arch.h"
#include "metrics.h"


typedef class Router{
    ARCH_TYPE arch_type;
    ROUTING_TYPE routing_type;
    Tile * host_tile;               
    int **routingTable;      //routingTable[i][j] shows how to which link to send the 
    //from tile i to tile j. If it's -2, means not reachable.
    //if it's -1, means the destination is the current tile.

    bool generate_xy_routing_table();
    bool generate_dally_xy_routing_table();
public:
    Router();
    bool initialize(ARCH_TYPE arch_type, routing_type type, pTile host);
    int route_to_link(int src_tile, int dst_tile) const;
    int set_routing_entry(int src_tile, int dst_tile, int link_id);
    ~Router();
    int choose_egress(int src_tile, int dst_tile);
} *pRouter, Router;


#endif

