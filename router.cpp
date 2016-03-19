#include "arch.h"
#include "gVal.h"
#include "nocmap.h"

Tile::Router::Router() {
    routing_table = new int*[gTileNum];
    for(int i=0; i<gTileNum; i++) 
        routing_table[i] = new int[gTileNum];
    for(int i=0; i<gTileNum; i++) 
        for(int j=0; j<gTileNum; j++) 
            routing_table[i][j] = -2;
}

Tile::Router::~Router() {
    if(routing_table) {
        for(int i=0; i<gTileNum; i++) 
            delete []routing_table[i];
        delete []routing_table;
    }
}

bool Tile::Router::initialize(ROUTING_TYPE r_type, pTile host) 
{
    host_tile = host;
    routing_type = r_type;
    if (routing_type == routing_xy) 
        return generate_xy_routing_table();
    else {
        cerr<<"This architecture and routing combination is not supported."<<endl;
        return false;
    }
}

bool Tile::Router::generate_xy_routing_table() 
{
	//This method is used to generate the fixed, non-adaptive, minimal routing 
	//table. This method is only applicable to Mesh or Torus

    for(int i=0; i<totalTiles; i++) 
        for(int j=0; j<totalTiles; j++) 
            routing_table[i][j] = -2;
  
    for(int dstTile=0; dstTile<gTileNum; dstTile++) {
        if(dstTile == host_tile->id) {                     //deliver to me
            routing_table[0][dstTile] = -1;
            continue;
        }

        //check out the dst Tile's position first
        Position dstPos;
        dstPos.row = dstTile / g_edge_size;
        dstPos.col = dstTile % g_edge_size;

        Position nextStep = host_tile->pos;

        if(dstPos.col != host_tile->pos.col) {            //We should go horizontally
            if(host_tile->pos.col > dstPos.col) 
                nextStep.col--;
            else 
                nextStep.col++;
        }
        else {            //We should go vertically
#ifdef DEBUG
            //Sanity checking
            if(dstPos.row == host_tile->pos.row) {
                cerr << "Routing table creation error. " << endl;
                return false;
            }
#endif //DEBUG
            if(host_tile->pos.row > dstPos.row) 
                nextStep.row--;
            else 
                nextStep.row++;
        }

        int i=0;
        for(i=0; i<host_tile->goLinkNum; i++) {
            pLink pL = host_tile->goLink[i];
            if(pL->ToTilePos() == nextStep) {
                routing_table[0][dstTile] = pL->GetId();
                break;
            }
        }
#ifdef DEBUG
        //Sanity checking
        if(i==host_tile->goLinkNum) {
            cerr << "Routing table creation error." << endl;
            return false;
        }
#endif //DEBUG
    }

    //Duplicate this routing row to the other routing rows.
    for(int i=1; i<gTileNum; i++) 
        for(int j=0; j<gTileNum; j++) 
            routing_table[i][j] = routing_table[0][j];

    return true;
}

int Tile::Router::route_to_link(int src_tile, int dst_tile) const {
    return routing_table[src_tile][dst_tile];
}

int Tile::Router::set_routing_entry(int src_tile, int dst_tile, int link_id) {
    routing_table[src_tile][dst_tile] = link_id;
    return link_id;
}

int Tile::Router::choose_egress(int src_tile, int dst_tile) {
    return routing_table[src_tile][dst_tile];
}
