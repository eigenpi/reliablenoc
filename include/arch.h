#ifndef ARCH_HPP_
#define ARCH_HPP_

#include <iostream>
#include <vector>
#include "param.h"


////////////////////////////////////////////////////////////////////////////////
//
// Position
//
////////////////////////////////////////////////////////////////////////////////

// aid == operation for float numbers
#define DELTA 0.0001

#define NORTH 0
#define SOUTH 1
#define EAST 2
#define WEST 3
#define SOURCE 4 // Local processing element data source
#define SINK 5 // Local processing element data sink

typedef struct Position {
	int row;
	int col;
} *pPosition, Position;

typedef struct Position3D {
	Position pos;
	int dir;           // North, south, east ...
} *pPosition3D, Position3D;

////////////////////////////////////////////////////////////////////////////////
//
// Link
//
////////////////////////////////////////////////////////////////////////////////

typedef class Link {
    friend ostream & operator<<(ostream &os, const Link & link);
    static int cnt;
    float powerCoeff; // Average energy consumption of sending one bit data
                      // over that link
    int id;
    int fromTileId;
    int toTileId;
 public:
    Position fromTile;
    Position toTile;
    int used_BW;
 public:
    int bandwidth;
    float Cost() {return powerCoeff;}
    int GetId() {return id;}
    Link();
    int const GetCnt();
    const Position & FromTilePos();
    const Position & ToTilePos();
    int FromTile() { return fromTileId;}
    int ToTile() {return toTileId;}
    int direction();
    friend bool verify_BW_requirement();
    friend int locate_link(Position & pos, int direction);
} *pLink, Link;

////////////////////////////////////////////////////////////////////////////////
//
// Tile
//
////////////////////////////////////////////////////////////////////////////////

typedef class Tile {
    static int cnt;
    static int totalTiles;
    //    static int rowNum;
    int id;
    int switch_delay;        //The delay of the switch in the tile (in cycles unit)
    int arbitration_delay;   //The delay of the routing. In wormhole routing, only
                             //processing the header flit will incur this delay
    float powerCoeff;        //The average energy consumption of routing one bit data
    Position pos;
    int goLinkNum;
    int comeLinkNum;
    vector<pLink> goLink;
    vector<pLink> comeLink;
    int procId;

    class Router {
        ROUTING_TYPE routing_type;
        Tile * host_tile;               
        int **routing_table;      //routing_table[i][j] shows how to which link to send the 
        //from tile i to tile j. If it's -2, means not reachable.
        //if it's -1, means the destination is the current tile.
    
        bool generate_xy_routing_table();
    public:
        Router();
        bool initialize(ROUTING_TYPE type, Tile * host);
        int route_to_link(int src_tile, int dst_tile) const;
        int set_routing_entry(int src_tile, int dst_tile, int link_id);
        ~Router();
        int choose_egress(int src_tile, int dst_tile);
    } router;

 public:
    Tile();
    float Cost() {return powerCoeff;}
    int mapToProc() { return procId;}
    int mapToProc(int pId) { procId = pId; return pId;}
    int GetId() const;
    int GetGoLinkNum() const;
    int GetComeLinkNum() const;

    bool initialize_router(ROUTING_TYPE routing_type);
    int RouteToLink(int srcId, int dstId) const;
    int RouteToLink(const Tile& srcTile, const Tile& dstTile) const;
    int RouteToLink(int srcId, int dstId, int linkId);
    int set_routing_entry(int src_tile, int dst_tile, int link_id);
    Position GetPosition() const;
    pLink GoLink(int i);
    pLink ComeLink(int i);
    void AttachLink(pLink gLink);
    ~Tile();
    friend ostream & operator<<(ostream & os, const Tile & tile);
} *pTile, Tile;

/************************************************************
findLink is used to find out the id of the link that connects
the source tile with the destination tile. If no link exists, 
return -1
************************************************************/
int findLink(int srcTile, int dstTile);

#endif

