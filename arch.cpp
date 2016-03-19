#include "arch.h"
#include "nocmap.h"
#include "gVal.h"
#include <cstdlib>

int Tile::GetId() const {
    return id;
}

int Tile::GetGoLinkNum() const {
    return goLinkNum;
}

int Tile::GetComeLinkNum() const {
    return comeLinkNum;
}

pLink Tile::GoLink(int i) {
    return goLink[i];
}

pLink Tile::ComeLink(int i) {
    return comeLink[i];
}

Tile::Tile() {
    totalTiles = gTileNum;

    id = cnt++;

    pos.row = id / g_edge_size;
    pos.col = id % g_edge_size;

    powerCoeff = param.switch_ebit;
}

int Tile::set_routing_entry(int src_tile, int dst_tile, int link_id) {
    return router.set_routing_entry(src_tile, dst_tile, link_id);
}

Tile::~Tile() {
    goLink.clear();
    comeLink.clear();
}

void Tile::AttachLink(pLink gLink)
{
    goLinkNum = comeLinkNum = 0;
    for (int i=0; i<gLink[0].GetCnt(); i++) {
        if (gLink[i].FromTilePos() == pos) {
            goLink.push_back(&gLink[i]);
            goLinkNum++;
        }
        if (gLink[i].ToTilePos() == pos) {
            comeLink.push_back(&gLink[i]);
            comeLinkNum++;
        }
    }
}

bool Tile::initialize_router(ROUTING_TYPE routing_type) {
    return router.initialize(routing_type, this);
}

Position Tile::GetPosition() const {
    return pos;
}

int Tile::RouteToLink(int srcId, int dstId) const {
    return router.route_to_link(srcId, dstId);
}

int Tile::RouteToLink(const Tile& srcTile, const Tile& dstTile) const {
    return router.route_to_link(srcTile.id, dstTile.id);
}

//This method is used to program the routing table
int Tile::RouteToLink(int srcId, int dstId, int linkId) {
    return router.set_routing_entry(srcId, dstId, linkId);
}

ostream & operator << (ostream &os, const Tile & tile) {
    os << "Tile " << tile.id << ":" << endl;
    os << "Leaving Links:" << endl;
    for (int i=0; i<tile.goLinkNum; i++) 
        os << (*tile.goLink[i]);
    os << "Entering Links:" << endl;
    for (int i=0; i<tile.comeLinkNum; i++) 
        os << (*tile.comeLink[i]);
    return os;
}


// definition for class Link.
Link::Link() 
{ 
    id = cnt++;
  
    // There are totally 2*(g_edge_size-1)*g_edge_size*2 links. The first half links 
	// are horizontal the second half links are veritcal links. 
    if (id < 2*(g_edge_size-1)*g_edge_size) {
        fromTile.row = id/(2*(g_edge_size-1));
        toTile.row = id/(2*(g_edge_size-1));
        int localId = id%(2*(g_edge_size-1));
        if (localId < (g_edge_size-1)) {
            //from west to east
            fromTile.col = localId;
            toTile.col = localId + 1;
        }
        else {
            //from east to west
            localId = localId - (g_edge_size-1);
            fromTile.col = localId + 1;
            toTile.col = localId;
        }
    }
    else {
        int localId = id - 2*(g_edge_size-1)*g_edge_size;
        fromTile.col = localId/(2*(g_edge_size-1));
        toTile.col = localId/(2*(g_edge_size-1));
        localId = localId%(2*(g_edge_size-1));
        if (localId < (g_edge_size-1)) {
            //from south to north
            fromTile.row = localId;
            toTile.row = localId + 1;
        }
        else {
            //from north to south
            localId = localId - (g_edge_size-1);
            fromTile.row = localId + 1;
            toTile.row = localId;
        }
    }
    //For mesh, all the links have the same power coefficency
    powerCoeff = param.link_ebit;

    fromTileId = fromTile.row*g_edge_size + fromTile.col;
    toTileId = toTile.row*g_edge_size + toTile.col;

    bandwidth = param.link_bandwidth;
}

int const Link::GetCnt() {
    return cnt;
}

const Position & Link::FromTilePos() {
    return fromTile;
}

const Position & Link::ToTilePos() {
    return toTile;
}

// return the direction of the link
int Link::direction() {
    Position from, to;
    from = gTile[fromTileId].GetPosition();
    to = gTile[toTileId].GetPosition();
    if (from.row - to.row == 1)
        return SOUTH;
    else if (from.row - to.row == -1)
        return NORTH;
    else if (from.col - to.col == 1)
        return WEST;
    else if (from.col - to.col == -1)
        return EAST;

    //should never reach here. 
    assert(0);
    return -1;
}

ostream & operator<<(ostream &os, const Link & link) {
    os << "Link " << link.id << ": " << "from tile("
       << link.fromTile.row << "," << link.fromTile.col << ") to tile("
       << link.toTile.row << "," << link.toTile.col << ")" << endl;
    return os;
}

int findLink(int srcTile, int dstTile) {
    for (int i=0; i<gLinkNum; i++) {
        if (gLink[i].FromTile()==srcTile && gLink[i].ToTile()==dstTile)
            return i;
    }
    return -1;
}

