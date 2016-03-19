#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include "arch.h"
#include "app.h"
#include "gVal.h"
#include "sa.h"


void branchAndBoundMapping();

////////////////////////////////////////////////////////////////////////////////
//
// MAPPING_NODE
//
////////////////////////////////////////////////////////////////////////////////

class MappingNode
{
public:
    static int cnt;
private:
    bool illegal;    // It is an illegal node if it violates the spec
                     // constructor will init this
    int stage;       // How many processes have been mapped
    int *mappingSequency;
    bool *tileOccupancyTable;
    float cost;
    float lowerBound;
    float upperBound;
    bool occupancyTableReady;
    bool illegal_child_mapping;
    int *link_BW_usage;   //The array record how much BW have been used 
                          //for each link, used for xy-routing
    int ***R_syn_link_BW_usage;  //this array is used when we also
                                 //need to synthesize the routing table
                                 //This can be indexed by 
                                 //[src_row][src_col][direction]
    int ***R_syn_link_BW_usage_temp; //working space
    int ****routing_table;           //[row][col][src_tile][dst_tile]

    /**********************************************************************
     * The following three member are useful only in routing_synthesis    *
     * mode and the routing_effort is HARD                                *
     **********************************************************************/
    int *routing_bit_array;        //0: route in X; 1: route in Y
    int routing_int;               //The routing_bit_array in integer form
    int *best_routing_bit_array; 
    bool first_routing_path;
    int MAX_ROUTING_INT;

    class MappingNode* next;
    friend class PQueue;
    friend void branchAndBoundMapping();
    void GreedyMapping();
    float LowerBound();
    float UpperBound();

    // In the fixed routing case, this is used to check the BW usage
    bool fixed_verify_BW_usage();

    float lowestUnitCost(int tileId);
    float lowestUnmappedUnitCost();
	float lowestUnitCost_multiobjective(int tileId, int comm_vol_ij);
	float lowestUnmappedUnitCost_reliability();
    void MapNode(int procId, float goodRow, float goodCol);

    // route the traffic related to the processes specifed between
    // begin_stage and end_stage, syn_routing controls whether it
    // needs to write the routing table
    bool route_traffics(int begin_stage, int end_stage, 
                        bool commit=true, bool update_routing_table=false);  
    bool route_traffic_easy(int src_tile, int dst_tile, int BW, 
                            bool commit, bool update_routing_table);
    int calc_adaptivity(int src_tile, int dst_tile, int BW) const;


    /**********************************************************************
     * For the routing table generation in hard mode                      *
     **********************************************************************/
    bool route_traffic_hard(int src_tile, int dst_tile, int BW, 
                            bool commit, bool update_routing_table);
    bool init_routing_path_generator(int x_hop, int y_hop);
    bool next_routing_path(int x_hop, int y_hop);
    int path_BW_usage(Position & src, Position & dst, int ***BW_usage, int BW);
    bool one_bits(int r, int onebits);

    void create_BW_temp_memory();
    void release_BW_temp_memory();
    void generate_routing_table();
    void remove_routing_table();

public:
    MappingNode(const MappingNode& parent, int tileId, bool calcBound=true);
    MappingNode(int tileId);
    MappingNode(const MappingNode& origin);
    ~MappingNode();

    int mapToTile(int i) {return mappingSequency[i];}
    int bestUpperBoundCandidate();
    int bestCostCandidate();
	int bestCostCandidate_multiobjective();
    bool Expandable(int tileId);
    int Stage() {return stage;}
    bool program_routers();
    bool is_illegal(void) const { return illegal; }
	double reliability_cost( int dx, int dy);
};
typedef class MappingNode MappingNode, *pMappingNode;

////////////////////////////////////////////////////////////////////////////////
//
// PQueue
//
////////////////////////////////////////////////////////////////////////////////

// A class for priority queue. 
// TODO: use heapsort data structure for speeding up
typedef class PQueue {
private:
    int length;
    pMappingNode head;
    int minCost;
    int minLowerBound;
    int minUpperBound;
public:
    PQueue();
    ~PQueue();
    int Length() {return length;}
    bool empty();
    void insert(pMappingNode node);
    pMappingNode next();
} PQueue, *pPQueue;

////////////////////////////////////////////////////////////////////////////////
//
// methods; these should be part of MAPPING class;
//
////////////////////////////////////////////////////////////////////////////////

void optimizeMapping(); // contains the main show of branch-and-bound;

void BBMInit();
void BBMSortProcesses();
void BBMBuildProcessMatrix();
void BBMBuildArchitectureMatrix();
void BBMMap(pMappingNode bestMapping);
void BBMClear();
// functions related to reliability;
void BBM_Compute_normalization_factors();
void reliability_costs_investigation();
double reliability_cost_clone( int dx, int dy);

#endif
