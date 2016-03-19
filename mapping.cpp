#include "mapping.h"
#include <sys/times.h>
#include <string.h>
#include <vector>

// again globals? they should be sketch variables inside MAPPING class;
int ** procMatrix = NULL;
float ** archMatrix = NULL;

// proc_map_array[i] represents the actual process that the i-th mapped process
// corresponding to
int * proc_map_array = NULL;

const float MAX_VALUE = INT_MAX - 100;

// each newly mapped communication transaction should be less than
// this value. Useful for non-regular region mapping
float MAX_PER_TRAN_COST = MAX_VALUE; 

// normalization values;
float COST_ENERGY_NORMALIZATION = 1.0;
float COST_RELIABILITY_NORMALIZATION = 1.0;

////////////////////////////////////////////////////////////////////////////////
//
// MAPPING_NODE
//
////////////////////////////////////////////////////////////////////////////////

MappingNode::MappingNode(const MappingNode& parent, int tileId, bool calcBound) 
{
    illegal = false;

    tileOccupancyTable = NULL;
    mappingSequency = NULL;
    link_BW_usage = NULL;
    R_syn_link_BW_usage = NULL;
    R_syn_link_BW_usage_temp = NULL;
    routing_table = NULL;

    routing_bit_array = NULL;  
    best_routing_bit_array = NULL; 

    occupancyTableReady = false;
    lowerBound = -1;
    cnt ++;
    mappingSequency = new int[gProcNum];
    for (int i=0; i<gProcNum; i++) 
        mappingSequency[i] = -1;

    stage = parent.stage;

    int proc1 = proc_map_array[stage];
    if (gProcess[proc1]->is_locked() && gProcess[proc1]->lock_to != tileId) {
        illegal = true;
        return;
    }
    lowerBound = parent.lowerBound;
    upperBound = parent.upperBound;


    // (1) Copy the parent's partial mapping
    memcpy(mappingSequency, parent.mappingSequency, gProcNum*sizeof(int));

    if (param.routing_table_synthesis) {
        // Copy the parent's link bandwidth usage
        R_syn_link_BW_usage = new int**[g_edge_size];
        for (int i=0; i<g_edge_size; i++) {
            R_syn_link_BW_usage[i] = new int*[g_edge_size];
            for (int j=0; j<g_edge_size; j++) {
                R_syn_link_BW_usage[i][j] = new int[4];
                memcpy(R_syn_link_BW_usage[i][j], parent.R_syn_link_BW_usage[i][j],
                       sizeof(int)*4);
            }
        }
    }
    else {
        link_BW_usage = new int[gLinkNum];
        memcpy(link_BW_usage, parent.link_BW_usage, sizeof(int)*gLinkNum);
    }

    // (2) Map the next process to tile tileId;
	// compute cost of mapping to tileId;
    mappingSequency[stage] = tileId;
    next = NULL;
    cost = parent.cost;

	float cost_e = 0;
	float cost_r = 0;
	int temp_dist = 0, delta_x = 0, delta_y = 0;
	
    for (int i=0; i<stage; i++) {
        int tile1 = tileId;
        int tile2 = mappingSequency[i]; // processes which have been mapped before;
        float this_tran_cost = procMatrix[i][stage];
        this_tran_cost = this_tran_cost * archMatrix[tile1][tile2];
        cost_e += this_tran_cost;
        if (this_tran_cost > MAX_PER_TRAN_COST) {
            illegal = true;
            return;
        }
		if ( param.do_reliability) { // reliability;
			if ( procMatrix[i][stage] > 0) {
				Position src = gTile[tile1].GetPosition();
				Position dst = gTile[tile2].GetPosition();
				delta_x = abs(dst.col - src.col);
				delta_y = abs(dst.row - src.row);
				temp_dist = delta_x + delta_y;
				cost_r += reliability_cost(delta_x,delta_y);
			}	
		}
    }

	if ( !param.do_reliability) { // OLD
		cost += cost_e;
	} else { // NEW
		//assert(cost_e < COST_ENERGY_NORMALIZATION);
		//assert(cost_r < COST_RELIABILITY_NORMALIZATION);
		cost_e /= COST_ENERGY_NORMALIZATION;
		cost_r /= COST_RELIABILITY_NORMALIZATION;
		cost += ( (1 - param.alpha) * cost_e + param.alpha * cost_r);
	}

    if (param.routing_table_synthesis) {
        if (!route_traffics(stage, stage)) {
            cost = MAX_VALUE+1;
            illegal = true;
            return;
        }
    } else {
        for (int i=0; i<stage; i++) {
            int tile1 = tileId;
            int tile2 = mappingSequency[i];
            int proc1 = proc_map_array[stage];
            int proc2 = proc_map_array[i];
            if (gProcess[proc1]->to_BW_requirement[proc2]) {
                for (unsigned int i=0; i<param.link_usage_list[tile1][tile2].size(); i++) {
                    int link_id = param.link_usage_list[tile1][tile2][i];
                    link_BW_usage[link_id] += gProcess[proc1]->to_BW_requirement[proc2];
                    if (link_BW_usage[link_id] > gLink[link_id].bandwidth) {
                        cost = MAX_VALUE+1;
                        illegal = true;
                        return;
                    }
                }
            }
            if (gProcess[proc1]->from_BW_requirement[proc2]) {
                for (unsigned int i=0; i<param.link_usage_list[tile2][tile1].size(); i++) {
                    int link_id = param.link_usage_list[tile2][tile1][i];
                    link_BW_usage[link_id] += gProcess[proc1]->from_BW_requirement[proc2];
                    if (link_BW_usage[link_id] > gLink[link_id].bandwidth) {
                        cost = MAX_VALUE+1;
                        illegal = true;
                        return;
                    }
                }
            }
        }
	}

	// this process becomes a mapped one and therefore the stage is now of one
	// more process;
    stage ++;

    if (!calcBound)
        return;

    tileOccupancyTable = new bool[gTileNum];
    for (int i=0; i<gTileNum; i++)
        tileOccupancyTable[i] = false;

    lowerBound = LowerBound();
    upperBound = UpperBound();
}

MappingNode::MappingNode(int tileId)
{
    illegal = false;
    tileOccupancyTable = NULL;
    mappingSequency = NULL;
    link_BW_usage = NULL;
    R_syn_link_BW_usage = NULL;
    R_syn_link_BW_usage_temp = NULL;
    routing_table = NULL;
    occupancyTableReady = false;
    lowerBound = -1;

    routing_bit_array = NULL;  
    best_routing_bit_array = NULL; 

    cnt ++; // this is a static variable!
    mappingSequency = new int[gProcNum];
    for (int i=0; i<gProcNum; i++) 
        mappingSequency[i] = -1;

    stage = 1;
    mappingSequency[0] = tileId;
    next = NULL;
    cost = 0;

    int proc1 = proc_map_array[0];
    if (gProcess[proc1]->is_locked() && gProcess[proc1]->lock_to != tileId) {
        illegal = true;
        return;
    }


    tileOccupancyTable = new bool[gTileNum];
    for (int i=0; i<gTileNum; i++)
        tileOccupancyTable[i] = false;

    if (param.routing_table_synthesis) {
        R_syn_link_BW_usage = new int**[g_edge_size];
        for (int i=0; i<g_edge_size; i++) {
            R_syn_link_BW_usage[i] = new int*[g_edge_size];
            for (int j=0; j<g_edge_size; j++) {
                R_syn_link_BW_usage[i][j] = new int[4];
                for (int k=0; k<4; k++) 
                    R_syn_link_BW_usage[i][j][k] = 0;
            }
        }
    }
    else{ 
        link_BW_usage = new int[gLinkNum];
        for (int i=0; i<gLinkNum; i++) 
            link_BW_usage[i] = 0;
    }

    lowerBound = LowerBound();
    upperBound = UpperBound();
}

MappingNode::MappingNode(const MappingNode & origin) 
{
	// essentially, this is to generate a mapping node which is a copy of
	// the node origin

    tileOccupancyTable = NULL;
    mappingSequency = NULL;
    link_BW_usage = NULL;
    R_syn_link_BW_usage = NULL;
    R_syn_link_BW_usage_temp = NULL;
    routing_table = NULL;

    routing_bit_array = NULL;  
    best_routing_bit_array = NULL; 

    occupancyTableReady = false;
    lowerBound = -1;
    cnt ++;
    mappingSequency = new int[gProcNum];
    for (int i=0; i<gProcNum; i++)
        mappingSequency[i] = -1;
    stage = origin.stage;
    illegal = origin.illegal;

    //Copy the parent's partial mapping
    for (int i=0; i<gProcNum; i++)
        mappingSequency[i] = origin.mappingSequency[i];

    if (param.routing_table_synthesis) {
        //Copy the parent's link bandwidth usage
        R_syn_link_BW_usage = new int**[g_edge_size];
        for (int i=0; i<g_edge_size; i++) {
            R_syn_link_BW_usage[i] = new int*[g_edge_size];
            for (int j=0; j<g_edge_size; j++) {
                R_syn_link_BW_usage[i][j] = new int[4];
                memcpy(R_syn_link_BW_usage[i][j], origin.R_syn_link_BW_usage[i][j],
                       sizeof(int)*4);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// lower and upper bound costs calculations
//
////////////////////////////////////////////////////////////////////////////////


float MappingNode::LowerBound() 
{
	//This calculate the lower bound cost of the unmapped process nodes
	//in the current mapping

    for (int i=0; i<gTileNum; i++) 
        tileOccupancyTable[i] = false;
    for (int i=0; i<stage; i++) 
        tileOccupancyTable[mappingSequency[i]] = true;

    occupancyTableReady = true;
    lowerBound = cost;

    // (1) The first part of the cost is the communication between those
    // nodes that have been mapped and those that have not yet. 
    // We assume that the unmapped node can occupy the unoccupied tile 
    // which has the lowest cost to the occupied node
    for (int i=0; i<stage; i++) {
        for (int j=stage; j<gProcNum; j++) {
            if (procMatrix[i][j]==0) {
				// there is no communication between mapped process i and 
				// unmapped process j;
                continue;
			} else {
				if ( !param.do_reliability) { // OLD
					// mappingSequency[i] tells me to what tile id this process it
					// has been mapped; procMatrix[i][j] gives me the total comm volume
					// between processes ranked i and j (remember - the hightest ranked
					// processes have been mapped first);
					lowerBound += procMatrix[i][j] * lowestUnitCost(mappingSequency[i]);
				} else { // NEW
					lowerBound += lowestUnitCost_multiobjective(mappingSequency[i], procMatrix[i][j]);
				}
			}
        }
    }

    // (2) Now add the cost of the communication among all the un-mapped nodes
	float cost_e = 0; // energy;
	float cost_r = 0; // reliability;
    int vol = 0;
	int count_comms = 0;
    for (int i=stage; i<gProcNum; i++) {
        for (int j=i+1; j<gProcNum; j++) {
            vol += procMatrix[i][j]; // aggregate comm volumes;
			if ( procMatrix[i][j] > 0) {
				count_comms ++;
			}
		}
    }
	if ( !param.do_reliability) { // OLD
		lowerBound += vol * // total comm volume, example: 1214278;
			lowestUnmappedUnitCost(); // energy, example: 1.017000;
	} else { // NEW
		cost_e = vol * lowestUnmappedUnitCost();
		cost_r = count_comms * lowestUnmappedUnitCost_reliability(); // example: 190;
		//assert(cost_e < COST_ENERGY_NORMALIZATION);
		//assert(cost_r < COST_RELIABILITY_NORMALIZATION);
		cost_e /= COST_ENERGY_NORMALIZATION;
		cost_r /= COST_RELIABILITY_NORMALIZATION;
		lowerBound += ( (1 - param.alpha) * cost_e + param.alpha * cost_r );
	}

    return ( lowerBound);
}


float MappingNode::UpperBound() 
{
	//This calculate the upper bound cost of the this partial mapping
	//in the current mapping

    if (!occupancyTableReady) {
        for (int i=0; i<gTileNum; i++) 
            tileOccupancyTable[i] = false;
        for (int i=0; i<stage; i++) 
            tileOccupancyTable[mappingSequency[i]] = true;
    }

	// (1) do "temporary" best mapping then compute its cost;
    GreedyMapping();
    upperBound = cost;

    illegal_child_mapping = false;

	// (2) special cases only;
    if (param.routing_table_synthesis) {
        create_BW_temp_memory();
        if (!route_traffics(stage, gProcNum-1, false)) {
            illegal_child_mapping = true;
            upperBound = MAX_VALUE;
            release_BW_temp_memory();
            return upperBound;
        }
        release_BW_temp_memory();
    }
    else if (!fixed_verify_BW_usage()) {
        illegal_child_mapping = true;
        upperBound = MAX_VALUE;
        return upperBound;
    }

	// (3) actual upper bound aggregation;
	float cost_e = 0.0;
	float cost_r = 0.0;
	// energy;
    for (int i=0; i<stage; i++) {
        int tile1 = mappingSequency[i]; 
        for (int j=stage; j<gProcNum; j++) {
            int tile2 = mappingSequency[j];
            cost_e += procMatrix[i][j] * archMatrix[tile1][tile2];
        }
    }
    for (int i=stage; i<gProcNum; i++) {
        int tile1 = mappingSequency[i];
        for (int j=i+1; j<gProcNum; j++) {
            int tile2 = mappingSequency[j];
            cost_e += procMatrix[i][j] * archMatrix[tile1][tile2];
        }
    }
	
	if ( !param.do_reliability) { // OLD
		upperBound += cost_e; // example: 2946321.70;
	} else { // NEW
		// reliability;
		if ( param.alpha > 0) {
			// otherwise it is not worth doing calculations; Note: also I do not why
			// upperBound would even with alpha = 0 be making the trajectory of the
			// algo go nuts and not find a solution? there is a bug somewhere...
			int temp_dist = 0, delta_x = 0, delta_y = 0;
			for (int i=0; i<stage; i++) { // between mapped and "un-mapped";
				int tile1 = mappingSequency[i]; 
				for (int j=stage; j<gProcNum; j++) {
					int tile2 = mappingSequency[j];
					if ( procMatrix[i][j] > 0) {
						Position src = gTile[tile1].GetPosition();
						Position dst = gTile[tile2].GetPosition();
						delta_x = abs(dst.col - src.col);
						delta_y = abs(dst.row - src.row);
						temp_dist = delta_x + delta_y;
						cost_r += reliability_cost(delta_x,delta_y);
					}
				}
			}
			for (int i=stage; i<gProcNum; i++) { // between "un-mapped";
				int tile1 = mappingSequency[i];
				for (int j=i+1; j<gProcNum; j++) {
					int tile2 = mappingSequency[j];
					if ( procMatrix[i][j] > 0) {
						Position src = gTile[tile1].GetPosition();
						Position dst = gTile[tile2].GetPosition();
						delta_x = abs(dst.col - src.col);
						delta_y = abs(dst.row - src.row);
						temp_dist = delta_x + delta_y;
						cost_r += reliability_cost(delta_x,delta_y);
					}
				}
			}
		}
		//assert(cost_e < COST_ENERGY_NORMALIZATION);
		//assert(cost_r < COST_RELIABILITY_NORMALIZATION);
		cost_e /= COST_ENERGY_NORMALIZATION;
		cost_r /= COST_RELIABILITY_NORMALIZATION;
		upperBound += ( (1 - param.alpha) * cost_e + param.alpha * cost_r );
		//printf(" (%.2f %.2f)",cost_e, cost_r);
	}

    return ( upperBound);
}


bool MappingNode::fixed_verify_BW_usage() 
{
    int *link_BW_usage_temp = new int[gLinkNum];
    memcpy(link_BW_usage_temp, link_BW_usage, sizeof(int)*gLinkNum);

    for (int i=0; i<stage; i++) {
        int tile1 = mappingSequency[i];
        int proc1 = proc_map_array[i];
        for (int j=stage; j<gProcNum; j++) {
            int tile2 = mappingSequency[j];
            int proc2 = proc_map_array[j];
            if (gProcess[proc1]->to_BW_requirement[proc2]) {
                for (unsigned int i=0; i<param.link_usage_list[tile1][tile2].size(); i++) {
                    int link_id = param.link_usage_list[tile1][tile2][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->to_BW_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > gLink[link_id].bandwidth) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }

            if (gProcess[proc1]->from_BW_requirement[proc2]) {
                for (unsigned int i=0; i<param.link_usage_list[tile2][tile1].size(); i++) {
                    int link_id = param.link_usage_list[tile2][tile1][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->from_BW_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > gLink[link_id].bandwidth) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }
        }
    }
    for (int i=stage; i<gProcNum; i++) {
        int tile1 = mappingSequency[i];
        int proc1 = proc_map_array[i];
        for (int j=i+1; j<gProcNum; j++) {
            int tile2 = mappingSequency[j];
            int proc2 = proc_map_array[j];
            if (gProcess[proc1]->to_BW_requirement[proc2]) {
                for (unsigned int i=0; i<param.link_usage_list[tile1][tile2].size(); i++) {
                    int link_id = param.link_usage_list[tile1][tile2][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->to_BW_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > gLink[link_id].bandwidth) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }

            if (gProcess[proc1]->from_BW_requirement[proc2]) {
                for (unsigned int i=0; i<param.link_usage_list[tile2][tile1].size(); i++) {
                    int link_id = param.link_usage_list[tile2][tile1][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->from_BW_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > gLink[link_id].bandwidth) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }
        }
    }
    delete []link_BW_usage_temp;
    return true;
}


//This method returns the tile to be mapped for the next node which will
//lead to the smallest partial mapping cost
int MappingNode::bestCostCandidate() 
{
    float minimal = MAX_VALUE;
    for (int i=0; i<gTileNum; i++) 
        tileOccupancyTable[i] = false;
    for (int i=0; i<stage; i++) 
        tileOccupancyTable[mappingSequency[i]] = true;

    int index = -1;
    for (int tileId=0; tileId<gTileNum; tileId++) {
        if (tileOccupancyTable[tileId])
            continue;
        float additionalCost = 0;
        for (int i=0; i<stage; i++) {
            int tile1 = tileId;
            int tile2 = mappingSequency[i];
            additionalCost += procMatrix[i][stage] * archMatrix[tile1][tile2];
            if (additionalCost >= minimal)
                break;
        }
        if (additionalCost < minimal) {
            minimal = additionalCost;
			index = tileId;
		}
    }
    return index;
}
int MappingNode::bestCostCandidate_multiobjective() 
{
    float minimal = INT_MAX;
    for (int i=0; i<gTileNum; i++) 
        tileOccupancyTable[i] = false;
    for (int i=0; i<stage; i++) 
        tileOccupancyTable[mappingSequency[i]] = true;

    int index = -1;
    for (int tileId=0; tileId<gTileNum; tileId++) {
        if (tileOccupancyTable[tileId])
            continue;
		float cost_e = 0;
		float cost_r = 0;
		int temp_dist = 0, delta_x = 0, delta_y = 0;
        float additionalCost = 0;
        for (int i=0; i<stage; i++) {
            int tile1 = tileId;
            int tile2 = mappingSequency[i];
            cost_e += procMatrix[i][stage] * archMatrix[tile1][tile2];
			if ( procMatrix[i][stage] > 0) {
				Position src = gTile[tile1].GetPosition();
				Position dst = gTile[tile2].GetPosition();
				delta_x = abs(dst.col - src.col);
				delta_y = abs(dst.row - src.row);
				temp_dist = delta_x + delta_y;
				cost_r += reliability_cost(delta_x,delta_y);
			}
        }
		//assert(cost_e < COST_ENERGY_NORMALIZATION);
		//assert(cost_r < COST_RELIABILITY_NORMALIZATION);
		cost_e /= COST_ENERGY_NORMALIZATION;
		cost_r /= COST_RELIABILITY_NORMALIZATION;
		additionalCost = ( (1 - param.alpha) * cost_e + param.alpha * cost_r);

        if (additionalCost < minimal) {
            minimal = additionalCost;
			index = tileId;
		}
    }
    return index;
}

//This method returns the tile to be mapped for the next node with the
//criteria of the greedy mapping of the current one. 
int MappingNode::bestUpperBoundCandidate() 
{
    return mappingSequency[stage];
}


float MappingNode::lowestUnitCost(int tileId) 
{
	// used in lower bound calculations;
	// This function returns the lowest cost of the tileId to any unoccupied
	// tile.
    float min = INT_MAX;
    for (int i=0; i<gTileNum; i++) {
        if (i==tileId)
            continue;
        if (tileOccupancyTable[i])
            continue;
        if (archMatrix[tileId][i] < min)
            min = archMatrix[tileId][i]; // example:  53571 < 71428;
    }
    return min;
}
float MappingNode::lowestUnitCost_multiobjective(int tileId, int comm_vol_ij) 
{
	// combined objectives;
    float min_cost = INT_MAX;
	float temp_cost = 0;
	float cost_e = 0;
	float cost_r = 0;
	int temp_dist = 0, delta_x = 0, delta_y = 0;
    for (int i=0; i<gTileNum; i++) {
        if (i==tileId)
            continue;
        if (tileOccupancyTable[i])
            continue;

		if ( !param.do_reliability || param.alpha == 0) { // OLD
			temp_cost = comm_vol_ij * archMatrix[tileId][i];
		} else { // NEW
			cost_e = comm_vol_ij * archMatrix[tileId][i];
			Position src = gTile[tileId].GetPosition();
			Position dst = gTile[i].GetPosition();
			delta_x = abs(dst.col - src.col);
			delta_y = abs(dst.row - src.row);		
			temp_dist = delta_x + delta_y;
			cost_r = reliability_cost(delta_x,delta_y); // attention: not aggregation!
			//assert(cost_e < COST_ENERGY_NORMALIZATION);
			//assert(cost_r < COST_RELIABILITY_NORMALIZATION);
			cost_e /= COST_ENERGY_NORMALIZATION;
			cost_r /= COST_RELIABILITY_NORMALIZATION;
			// combined cost:
			temp_cost = (1 - param.alpha) * cost_e + param.alpha * cost_r;
		}
		if ( temp_cost < min_cost)
			min_cost = temp_cost;
    }
    return min_cost;
}

float MappingNode::lowestUnmappedUnitCost() 
{
	// This function returns the lowest cost between any two unoccupied tiles
    float min = INT_MAX;
    for (int i=0; i<gTileNum; i++) {
        if (tileOccupancyTable[i])
            continue;
        for (int j=i+1; j<gTileNum; j++) {
            if (tileOccupancyTable[j])
                continue;
            if (archMatrix[i][j]<min)
                min = archMatrix[i][j];
        }
    }
    return min;
}
float MappingNode::lowestUnmappedUnitCost_reliability() 
{
	// this function returns the lowest reliability cost between any two 
	// unoccupied tiles; for now this is the smallest Manhattan distance?
    float min_cost = INT_MAX;
	int temp_dist = 0, delta_x = 0, delta_y = 0;
	float temp_cost = 0.0;
    for (int i=0; i<gTileNum; i++) {
        if (tileOccupancyTable[i])
            continue;
        for (int j=i+1; j<gTileNum; j++) {
            if (tileOccupancyTable[j])
                continue;
			Position src = gTile[i].GetPosition();
			Position dst = gTile[j].GetPosition();
			delta_x = abs(dst.col - src.col);
			delta_y = abs(dst.row - src.row);		
			temp_dist = delta_x + delta_y;
			temp_cost = reliability_cost(delta_x,delta_y);
            if ( temp_cost < min_cost)
                min_cost = temp_cost;
        }
    }
    return float( min_cost);
}

void MappingNode::GreedyMapping() 
{
	// Map the other unmapped process node using greedy mapping
	// Note: this is only a "temporary" mapping for the sake of computing the
	// UBC (upper bound cost); this mapping is stored in mappingSequency[];

    for (int i=stage; i<gProcNum; i++) {
        int sumRow = 0;
        int sumCol = 0;
        int vol = 0;
        for (int j=0; j<i; j++) {
            if (procMatrix[i][j]==0)
                continue;
			// tile where mapped process j is located and with which process 
			// unmapped i communicates;
            int tileId = mappingSequency[j];
            int row = tileId / g_edge_size;
            int col = tileId % g_edge_size;
            sumRow += procMatrix[i][j]*row;
            sumCol += procMatrix[i][j]*col;
            vol += procMatrix[i][j];
        }
        // this is somehow the ideal position (it is kind of a weighted
		// center of gravity)
        float myRow, myCol;
        if (vol==0) {
            myRow = -1;
            myCol = -1;
        } 
        else{
            myRow = ((float) sumRow)/vol;
            myCol = ((float) sumCol)/vol;
        }
        MapNode(i, myRow, myCol);
    }
}

//Try to map the node to an unoccupied tile which is closest to 
//the tile(goodRow, goodCol)
void MappingNode::MapNode(int procId, float goodRow, float goodCol) 
{
    float minDist = 10000;
    int bestId = -1;
    for (int i=0; i<gTileNum; i++) {
        if (tileOccupancyTable[i])
            continue;
        if (goodRow<0) {
            bestId = i;
            break;
        }
        int row = i / g_edge_size;
        int col = i % g_edge_size;
        float dist = fabs(goodRow-row) + fabs(goodCol-col);
        if (dist<minDist) {
            minDist = dist;
            bestId = i;
        }
    }
    mappingSequency[procId] = bestId;
    tileOccupancyTable[bestId] = true;
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

bool MappingNode::Expandable(int tileId) 
{
    //If it's an illegal mapping, then just return false
    for (int i=0; i<stage; i++) {
		// Note: remember - mappingSequency[0] is initialized to the tile id to 
		// which the process is assigned to during the node ctor;
        if (mappingSequency[i] == tileId) //the tile has already been occupied
            return false;
    }
    return true;
}

void MappingNode::create_BW_temp_memory() 
{
    //Copy the bandwidth usage status to R_syn_link_BW_usage_temp
    R_syn_link_BW_usage_temp = new int**[g_edge_size];
    for (int i=0; i<g_edge_size; i++) {
        R_syn_link_BW_usage_temp[i] = new int*[g_edge_size];
        for (int j=0; j<g_edge_size; j++) {
            R_syn_link_BW_usage_temp[i][j] = new int[4];
            memcpy(R_syn_link_BW_usage_temp[i][j], 
                   R_syn_link_BW_usage[i][j], sizeof(int)*4);
        }
    }
}

void MappingNode::release_BW_temp_memory() {
    if (R_syn_link_BW_usage_temp) {
        for (int i=0; i<g_edge_size; i++) {
            for (int j=0; j<g_edge_size; j++) 
                delete []R_syn_link_BW_usage_temp[i][j];
            delete []R_syn_link_BW_usage_temp[i];
        }
        delete []R_syn_link_BW_usage_temp;
    }
}

MappingNode::~MappingNode(){
    delete []mappingSequency;

    if (link_BW_usage) {
        delete []link_BW_usage;
    }

    if (tileOccupancyTable) 
        delete []tileOccupancyTable;

    if (R_syn_link_BW_usage) {
        for (int i=0; i<g_edge_size; i++) {
            for (int j=0; j<g_edge_size; j++) 
                delete []R_syn_link_BW_usage[i][j];
            delete []R_syn_link_BW_usage[i];
        }
        delete []R_syn_link_BW_usage;
    }

    if (routing_bit_array) 
        delete []routing_bit_array;
    if (best_routing_bit_array)
        delete []best_routing_bit_array;
}

//fixing the routing tables of the tiles which are occupied by the processes
//from begin_stage to end_stage
bool MappingNode::route_traffics(int begin_stage, int end_stage, bool commit,
	bool update_routing_table)
{
    vector<Proc_Comm> Q;
    Q.clear();
    for (int cur_stage=begin_stage; cur_stage<=end_stage; cur_stage++) {
        int new_proc = proc_map_array[cur_stage];
        //Sort the request in the queue according to the BW
        //However, if the src and the dst are in the same row or in the same column, 
        //then we should insert it at the head of the queue.
        for (int i=0; i<cur_stage; i++) {
            int old_proc = proc_map_array[i];
            if (gProcess[new_proc]->to_BW_requirement[old_proc]) {
                Proc_Comm proc_comm;
                proc_comm.src_proc = cur_stage;   //we put virtual proc id
                proc_comm.dst_proc = i;
                proc_comm.BW = gProcess[new_proc]->to_BW_requirement[old_proc];
                proc_comm.adaptivity = calc_adaptivity(mappingSequency[proc_comm.src_proc], 
                                                       mappingSequency[proc_comm.dst_proc],
                                                       proc_comm.BW);
                vector<Proc_Comm>::iterator iter = Q.begin();
                while (iter != Q.end()) {
                    if ((proc_comm.adaptivity<iter->adaptivity) ||
                       (proc_comm.adaptivity==iter->adaptivity && proc_comm.BW>iter->BW))
                        break;
                    iter ++;
                }
                Q.insert(iter, proc_comm);
            }
            if (gProcess[new_proc]->from_BW_requirement[old_proc]) {
                Proc_Comm proc_comm;
                proc_comm.src_proc = i;
                proc_comm.dst_proc = cur_stage;
                proc_comm.BW = gProcess[new_proc]->from_BW_requirement[old_proc];
                proc_comm.adaptivity = calc_adaptivity(mappingSequency[proc_comm.src_proc],
                                                       mappingSequency[proc_comm.dst_proc],
                                                       proc_comm.BW);
                vector<Proc_Comm>::iterator iter = Q.begin();
                while (iter != Q.end()) {
                    if ((proc_comm.adaptivity<iter->adaptivity) ||
                       (proc_comm.adaptivity==iter->adaptivity && proc_comm.BW>iter->BW))
                        break;
                    iter ++;
                }
                Q.insert(iter, proc_comm);
            }
        }
    }
    //now route the traffic
    for (int i=0; ((unsigned) i)<Q.size(); i++) {
        int src_proc = Q[i].src_proc;
        int dst_proc = Q[i].dst_proc;
        int src_tile = mappingSequency[src_proc];
        int dst_tile = mappingSequency[dst_proc];
        int BW = Q[i].BW;
        if (param.routing_effort == easy) {
            if (!route_traffic_easy(src_tile, dst_tile, BW, commit, update_routing_table)) 
                return false;
        }
        else {
            if (!route_traffic_hard(src_tile, dst_tile, BW, commit, update_routing_table))
                return false;
        }
    }
    return true;
}

//currently there is only two levels of adaptivity. 0 for no adaptivity, 1 for (maybe) 
//some adaptivity
int MappingNode::calc_adaptivity(int src_tile, int dst_tile, int BW) const{
    Position src = gTile[src_tile].GetPosition();
    Position dst = gTile[dst_tile].GetPosition();
    int adaptivity;
    if (src.row==dst.row || src.col==dst.col) {
        adaptivity = 0;
        return adaptivity;
    }

    int ***BW_usage = R_syn_link_BW_usage;

    adaptivity = 1;
    int row = src.row;
    int col = src.col;
    int direction = -2;
    while (row!=dst.row || col!=dst.col) {
        //For west-first routing
        if (param.legal_turn_set == west_first) {
            if (col > dst.col) //step west
                return 0;
            else if (col == dst.col) 
                return 0;
            else if (row == dst.row) 
                return 0;
            //Here comes the flexibility. We can choose whether to go vertical or horizontal
            else { 
                int direction1 = (row<dst.row)?NORTH:SOUTH;
                if (BW_usage[row][col][direction1]+BW < param.link_bandwidth &&
                    BW_usage[row][col][EAST]+BW < param.link_bandwidth)
                    return 1;
                direction = (BW_usage[row][col][direction1]<BW_usage[row][col][EAST])? direction1:EAST;
            }
        }
        //For odd-even routing
        else if (param.legal_turn_set == odd_even) {
            int e0 = dst.col - col;
            int e1 = dst.row - row;
            if (e0==0) //currently the same column as destination
                direction = (e1>0)?NORTH:SOUTH;
            else {
                if (e0>0) {        //eastbound messages
                    if (e1==0)
                        direction = EAST;
                    else {
                        int direction1=-1, direction2=-1;
                        if (col%2==1 || col==src.col) 
                            direction1 = (e1>0)?NORTH:SOUTH;
                        if (dst.col%2==1 || e0!=1) 
                            direction2 = EAST;
                        if (direction1==-1)
                            direction = direction2;
                        else if (direction2==-1) 
                            direction = direction1;
                        else {//we have two choices
                            if (BW_usage[row][col][direction1]+BW < param.link_bandwidth &&
                                BW_usage[row][col][direction2]+BW < param.link_bandwidth)
                                return 1;
                            direction = 
                                (BW_usage[row][col][direction1]<BW_usage[row][col][direction2])?
                                direction1:direction2;
                        }
                    }
                }
                else {   //westbound messages
                    if (col%2!=0 || e1==0)
                        direction = WEST;
                    else {
                        int direction1 = (e1>0)?NORTH:SOUTH;
                        if (BW_usage[row][col][direction1]+BW < param.link_bandwidth &&
                            BW_usage[row][col][WEST]+BW < param.link_bandwidth)
                            return 1;
                        direction = 
                            (BW_usage[row][col][WEST]<BW_usage[row][col][direction1])?WEST:direction1;
                    }
                }
            }
        }
        switch(direction) {
        case SOUTH:
            row --;
            break;
        case NORTH:
            row ++;
            break;
        case EAST:
            col ++;
            break;
        case WEST:
            col --;
            break;
        default:
            cerr<<"Error"<<endl;
            break;
        }
    }
    return 0;
}

//Route the traffic from src_tile to dst_tile using BW
//Two routing methods are provided: west-first and odd-even routing
bool MappingNode::route_traffic_easy(int src_tile, int dst_tile, int BW,
	bool commit, bool update_routing_table) 
{
    Position src = gTile[src_tile].GetPosition();
    Position dst = gTile[dst_tile].GetPosition();
    int row = src.row;
    int col = src.col;

    int ***BW_usage = commit ? R_syn_link_BW_usage:R_syn_link_BW_usage_temp;
    int direction = -2;
    while (row!=dst.row || col!=dst.col) {
        //For west-first routing
        if (param.legal_turn_set == west_first) {
            if (col > dst.col) //step west
                direction = WEST;
            else if (col == dst.col) 
                direction = (row<dst.row)? NORTH:SOUTH;
            else if (row == dst.row) 
                direction = EAST;
            // Here comes the flexibility. We can choose whether to go
            // vertical or horizontal
 
			//else { 
			//int direction1 = (row<dst.row)?NORTH:SOUTH;
			//if (BW_usage[row][col][direction1] < BW_usage[row][col][EAST])
			//direction = direction1;
			//else if (BW_usage[row][col][direction1] > BW_usage[row][col][EAST])
			//direction = EAST;
			//else { 
			//In this case, we select the direction which has the longest 
			//distance to the destination
			//if ((dst.col-col)*(dst.col-col) <= (dst.row-row)*(dst.row-row)) 
			//direction = direction1;
			//else //Horizontal move
			//direction = EAST;
			//}
			//}
 
            else { 
                direction = EAST;
                if (BW_usage[row][col][direction]+BW > param.link_bandwidth)
                    direction = (row<dst.row)?NORTH:SOUTH;
            }
        }
        //For odd-even routing
        else if (param.legal_turn_set == odd_even) {
            int e0 = dst.col - col;
            int e1 = dst.row - row;
            if (e0==0) //currently the same column as destination
                direction = (e1>0)?NORTH:SOUTH;
            else {
                if (e0>0) {        //eastbound messages
                    if (e1==0)
                        direction = EAST;
                    else {
                        int direction1=-1, direction2=-1;
                        if (col%2==1 || col==src.col) 
                            direction1 = (e1>0)?NORTH:SOUTH;
                        if (dst.col%2==1 || e0!=1) 
                            direction2 = EAST;
                        if (direction1==-1&&direction2==-1) {
                            cerr<<"Error, debug me."<<endl;
                            exit(1);
                        }
                        if (direction1==-1)
                            direction = direction2;
                        else if (direction2==-1) 
                            direction = direction1;
                        else //we have two choices
                            direction = 
                                (BW_usage[row][col][direction1]<BW_usage[row][col][direction2])?
                                direction1:direction2;
                    }
                }
                else {   //westbound messages
                    if (col%2!=0 || e1==0)
                        direction = WEST;
                    else {
                        int direction1 = (e1>0)?NORTH:SOUTH;
                        direction = 
                            (BW_usage[row][col][WEST]<BW_usage[row][col][direction1])?WEST:direction1;
                    }
                }
            }
        }

        BW_usage[row][col][direction] += BW;
        if (BW_usage[row][col][direction] > param.link_bandwidth && (!update_routing_table))
            return false;
        if (update_routing_table)
            routing_table[row][col][src_tile][dst_tile] = direction;

        switch(direction) {
        case SOUTH:
            row --;
            break;
        case NORTH:
            row ++;
            break;
        case EAST:
            col ++;
            break;
        case WEST:
            col --;
            break;
        default:
            cerr<<"Error"<<endl;
            break;
        }
    }
    return true;
}

/**********************************************************************
 * Route the traffic from src_tile to dst_tile using BW in complex    *
 * model, by which it means select the path from all its candidate    *
 * pathes which has the minimal maximal BW usage :)                   *
 **********************************************************************/
bool MappingNode::route_traffic_hard(int src_tile, int dst_tile, int BW, 
				     bool commit, bool update_routing_table) {
    Position src = gTile[src_tile].GetPosition();
    Position dst = gTile[dst_tile].GetPosition();
    int row = src.row;
    int col = src.col;

    int ***BW_usage = commit ? R_syn_link_BW_usage:R_syn_link_BW_usage_temp;

    //We can arrive at any destination with 2*g_edge_size hops
    if (!routing_bit_array)
        routing_bit_array = new int[2*g_edge_size];
    if (!best_routing_bit_array)
        best_routing_bit_array = new int[2*g_edge_size];

    //In the following, we find the routing path which has the minimal maximal
    //link BW usage and store that routing path to best_routing_bit_array
    int min_path_BW = INT_MAX;
    int x_hop = src.col - dst.col;
    x_hop = (x_hop>=0) ? x_hop:(0-x_hop);
    int y_hop = src.row - dst.row;
    y_hop = (y_hop>=0) ? y_hop:(0-y_hop);

    init_routing_path_generator(x_hop, y_hop);

    while (next_routing_path(x_hop, y_hop)) {    //For each path
        int usage = path_BW_usage(src, dst, BW_usage, BW);
        if (usage<min_path_BW) {
            min_path_BW = usage;
            memcpy(best_routing_bit_array, routing_bit_array, (x_hop+y_hop)*sizeof(int));
        }
    }

    if (min_path_BW == INT_MAX)
        return false;


    int direction = -2;

    int hop_id = 0;
    while (row!=dst.row || col!=dst.col) {

        if (best_routing_bit_array[hop_id++]) 
            direction = (row<dst.row) ? NORTH:SOUTH;
        else
            direction = (col<dst.col) ? EAST:WEST;

        BW_usage[row][col][direction] += BW;

        if ((BW_usage[row][col][direction] > param.link_bandwidth) && (!update_routing_table))
            return false;
        if (update_routing_table)
            routing_table[row][col][src_tile][dst_tile] = direction;

        switch(direction) {
        case SOUTH:
            row --;
            break;
        case NORTH:
            row ++;
            break;
        case EAST:
            col ++;
            break;
        case WEST:
            col --;
            break;
        default:
            cerr<<"Error"<<endl;
            break;
        }
    }
    return true;
}

/**********************************************************************
 * This function needs to fulfill two tasks. First, check to see      *
 * whether it's a valid path according to the selected routing alg.   *
 * Second, it checks to see if any BW requirement violates. If either *
 * of the two conditions is not met, return INT_MAX.                  *
 **********************************************************************/
int MappingNode::path_BW_usage(Position & src, Position & dst, 
			       int ***BW_usage, int BW) {
    int row = src.row;
    int col = src.col;

    int max_BW = 0;
    int hop_id = 0;

    while (row!=dst.row || col!=dst.col) {
        int direction = -2;
        // For west-first routing
        if (param.legal_turn_set == west_first) {
            if (col > dst.col) {   // step west
                direction = WEST;
                if (routing_bit_array[hop_id])
                    return INT_MAX;
            }
            else if (col == dst.col) {
                direction = (row<dst.row) ? NORTH:SOUTH;
                if (!routing_bit_array[hop_id])
                    return INT_MAX;
            }
            else if (row == dst.row) {
                direction = EAST;
                if (routing_bit_array[hop_id])
                    return INT_MAX;
            }
            // Here comes the flexibility. We can choose whether to go vertical or horizontal
            else { 
                int direction1 = (row<dst.row) ? NORTH : SOUTH;
                int direction2 = EAST;
                direction = (routing_bit_array[hop_id]) ? direction1 : direction2;
            }
        }
        // For odd-even routing
        else if (param.legal_turn_set == odd_even) {
            int e0 = dst.col - col;
            int e1 = dst.row - row;
            if (e0==0) {          // currently the same column as destination
                direction = (e1>0) ? NORTH : SOUTH;
                if (!routing_bit_array[hop_id])
                    return INT_MAX;
            }
            else {
                if (e0>0) {       // eastbound messages
                    if (e1==0) {
                        direction = EAST;
                        if (routing_bit_array[hop_id])
                            return INT_MAX;
                    }
                    else {
                        int direction1=-1, direction2=-1;
                        if (col%2==1 || col==src.col) 
                            direction1 = (e1>0) ? NORTH : SOUTH;
                        if (dst.col%2==1 || e0!=1) 
                            direction2 = EAST;
                        assert(!(direction1==-1&&direction2==-1));
                        if (direction1==-1) {
                            direction = direction2;
                            if (routing_bit_array[hop_id])
                                return INT_MAX;
                        }
                        else if (direction2==-1) {
                            direction = direction1;
                            if (!routing_bit_array[hop_id])
                                return INT_MAX;
                        }
                        else // we have two choices
                            direction = (routing_bit_array[hop_id]) ? direction1:direction2;
                    }
                }
                else {   // westbound messages
                    if (col%2!=0 || e1==0) {
                        direction = WEST;
                        if (routing_bit_array[hop_id])
                            return INT_MAX;
                    }
                    else {
                        int direction1 = (e1>0)?NORTH:SOUTH;
                        int direction2 = WEST;
                        direction = (routing_bit_array[hop_id]) ? direction1:direction2;
                    }
                }
            }
        }

        if (BW_usage[row][col][direction] > max_BW)
            max_BW = BW_usage[row][col][direction];

        if (BW_usage[row][col][direction]+BW > param.link_bandwidth)
            return INT_MAX;

        switch(direction) {
        case SOUTH:
            row --;
            break;
        case NORTH:
            row ++;
            break;
        case EAST:
            col ++;
            break;
        case WEST:
            col --;
            break;
        default:
            cerr<<"Error"<<endl;
            break;
        }
        hop_id++;

#ifdef DEBUG
        int x_hop = src.col - dst.col;
        x_hop = (x_hop>0) ? x_hop:(-x_hop);
        int y_hop = src.row - dst.row;
        y_hop = (y_hop>0) ? y_hop:(-y_hop);
        if (hop_id > (x_hop+y_hop)) {
            cerr<<"Error in routing. quit"<<endl;
            exit(1);
        }
#endif //DEBUG

    }
    return true;
}

bool MappingNode::init_routing_path_generator(int x_hop, int y_hop) {
    first_routing_path = true;
    MAX_ROUTING_INT = 0;
    for (int index=0; index<y_hop; index++) 
        MAX_ROUTING_INT = (MAX_ROUTING_INT<<1) + 1;
    MAX_ROUTING_INT = MAX_ROUTING_INT<<x_hop;
    return true;
}

bool MappingNode::next_routing_path(int x_hop, int y_hop) {
    if (first_routing_path) {
        first_routing_path = false;
        int index = 0;
        routing_int = 0;
        for (index=0; index<y_hop; index++) {
            routing_bit_array[index] = 1;
            routing_int = (routing_int<<1) + 1;
        }
        for (int x_index=0 ; x_index<x_hop; x_index++) 
            routing_bit_array[index+x_index] = 0;
        return true;
    }

    /**********************************************************************
     * find the next routing path based on the current routing_bit_array  *
     * the next one is the one which is the minimal array which is larger *
     * than the current routing_bit_array but with the same number of 1s  *
     **********************************************************************/
    while (routing_int<=MAX_ROUTING_INT) {
        if (routing_int%2==0)  //For an even number
            routing_int += 2;
        else
            routing_int ++;
        if (one_bits(routing_int, y_hop))
            break;
    }
    if (routing_int<=MAX_ROUTING_INT)
        return true;
    else
        return false;
}

//Returns true if the binary representation of r contains y_hop number of
//1s. It also assigns the bit form to routing_bit_array
bool MappingNode::one_bits(int r, int onebits) {
    memset(routing_bit_array, 0, 2*g_edge_size*sizeof(int));
    int index = 0;
    int cur_one_bits = 0;
    while (r) {
        routing_bit_array[index] = r&1;
        if (routing_bit_array[index]) 
            cur_one_bits ++;
        if (cur_one_bits > onebits)
            return false;
        index ++;
        r = r>>1;
    }
    if (cur_one_bits == onebits)
        return true;
    else
        return false;
}

bool MappingNode::program_routers()
{
    generate_routing_table();
    //clean all the old routing table
    for (int tile_id=0; tile_id<gTileNum; tile_id++) {
        for (int src_tile=0; src_tile<gTileNum; src_tile++) {
            for (int dst_tile=0; dst_tile<gTileNum; dst_tile++) {
                if (tile_id == dst_tile) 
                    gTile[tile_id].set_routing_entry(src_tile, dst_tile, -1);
                else
                    gTile[tile_id].set_routing_entry(src_tile, dst_tile, -2);
            }
        }
    }

    for (int row=0; row<g_edge_size; row++) {
        for (int col=0; col<g_edge_size; col++) {
            int tile_id = row*g_edge_size + col;
            for (int src_tile=0; src_tile<gTileNum; src_tile++) {
                for (int dst_tile=0; dst_tile<gTileNum; dst_tile++) {
                    Position pos;
                    pos.row = row; pos.col = col;
                    int link_id = locate_link(pos, routing_table[row][col][src_tile][dst_tile]);
                    if (link_id!=-1)
                        gTile[tile_id].set_routing_entry(src_tile, dst_tile, link_id);
                }
            }
        }
    }
    remove_routing_table();
    return true;
}

void MappingNode::generate_routing_table() 
{
    //reset all the BW_usage 
    for (int i=0; i<g_edge_size; i++) 
        for (int j=0; j<g_edge_size; j++) 
            for (int k=0; k<4; k++) 
                R_syn_link_BW_usage[i][j][k] = 0;

    routing_table = new int***[g_edge_size];
    for (int i=0; i<g_edge_size; i++) {
        routing_table[i] = new int**[g_edge_size];
        for (int j=0; j<g_edge_size; j++) {
            routing_table[i][j] = new int*[gTileNum];
            for (int k=0; k<gTileNum; k++) {
                routing_table[i][j][k] = new int[gTileNum];
                for (int m=0; m<gTileNum; m++)
                    routing_table[i][j][k][m] = -2;
            }
        }
    }

    //if it's a real child mappingnode.
    if (stage==gProcNum) 
        route_traffics(0, gProcNum-1, true, true);
    //if it's the node which generate min upperBound
    else {
        for (int i=0; i<stage; i++) 
            route_traffics(i, i, true, true);
        route_traffics(stage, gProcNum-1, true, true);
    }
}

void MappingNode::remove_routing_table() {
    if (routing_table) {
        for (int i=0; i<g_edge_size; i++) {
            for (int j=0; j<g_edge_size; j++) {
                for (int k=0; k<gTileNum; k++) 
                    delete []routing_table[i][j][k];
                delete []routing_table[i][j];
            }
            delete []routing_table[i];
        }
        delete []routing_table;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// PQueue
//
////////////////////////////////////////////////////////////////////////////////

PQueue::PQueue() {
    length = 0;
    head = NULL;
}

PQueue::~PQueue() {
}

bool PQueue::empty() {
    if (length==0) 
        return true;
    else 
        return false;
}

void PQueue::insert(pMappingNode node) 
{
    // here we should insert the node at the position which
    // is decided by the cost of the node
    if (length == 0) {
        head = node;
        node->next = NULL;
        length++;
        return;
    }
    pMappingNode parentNode = NULL;
    pMappingNode curNode = head;
    int i=0;
    for (i=0; i<length; i++) {
        if (curNode->cost > node->cost)
            //if (curNode->upperBound > node->upperBound)
            break;
        parentNode = curNode;
        curNode = curNode->next;
    }
    if (parentNode == NULL) {
        pMappingNode oldHead = head;
        head = node;
        node->next = oldHead;
        length++;
        return;
    }
    pMappingNode pNode = parentNode->next;
    parentNode->next = node;
    node->next = pNode;
    length++;
}

pMappingNode PQueue::next() 
{
    if (length==0)
        return NULL;
    pMappingNode oldHead = head;
    head = oldHead->next;
    length --;
    return oldHead;
}

////////////////////////////////////////////////////////////////////////////////
//
// methods of MAPPING; now these methods are globals but they should be made
// part of class MAPPING;
//
////////////////////////////////////////////////////////////////////////////////

//Initialization function for Branch-and-bound mapping
void BBMInit() 
{
    cout<<"Initialize for branch-and-bound"<<endl;
    if (proc_map_array)
        delete []proc_map_array;
    proc_map_array = new int[gProcNum];
    BBMSortProcesses();
    BBMBuildProcessMatrix();
    BBMBuildArchitectureMatrix();
	if ( param.do_reliability) {
		BBM_Compute_normalization_factors();
	}

    if (exist_non_regular_regions()) {
        // let's calculate the maximum ebit of sending a bit 
        // from one tile to its neighboring tile
        float max_e = -1.;
        for (int i=0; i<gLinkNum; i++) {
            if (gLink[i].Cost() > max_e)
                max_e = gLink[i].Cost();
        }
        float eb = max_e;
        max_e = -1;
        for (int i=0; i<gTileNum; i++) {
            if (gTile[i].Cost() > max_e)
                max_e = gTile[i].Cost();
        }
        eb += max_e*2;
        MAX_PER_TRAN_COST = eb * DUMMY_VOL * 1.3;   // let's put some overhead
    }
}

void BBMSortProcesses() 
{
	// sort the processes so that the branch-and-bound mapping can be 
	// accelerated
	// Note: this should be improved by using faster sorting;

    //sort them according to the sum of each process's ingress and egress 
    //communication volume
    for (unsigned int i=0; i<gProcess.size(); i++) {
        gProcess[i]->totalCommVol = 0;
        for (unsigned int k=0; k<gProcess.size(); k++) {
            gProcess[i]->totalCommVol += gProcess[i]->toComm[k];
            gProcess[i]->totalCommVol += gProcess[i]->fromComm[k];
        }
    }
    // Now rank them
    int cur_rank = 0;
    // locked PEs have the highest priority
    for (unsigned int i=0; i<gProcess.size(); i++) {
        if (gProcess[i]->is_locked()) {
            proc_map_array[cur_rank] = i;
            gProcess[i]->rank = cur_rank ++;
        }
        else
            gProcess[i]->rank = -1;
    }
    // the remaining PEs are sorted based on their comm volume
    for (unsigned int i=cur_rank; i<gProcess.size(); i++) {
        int max = -1;
        int maxid = -1;
        for (int k=0; k<gProcNum; k++) {
            if (gProcess[k]->rank != -1)
                continue;
            if (gProcess[k]->totalCommVol > max) {
                max = gProcess[k]->totalCommVol;
                maxid = k;
            }
        }
        gProcess[maxid]->rank = i;
        proc_map_array[i] = maxid;
    }
}

void BBMBuildProcessMatrix() 
{
	// procMatrix (defined at the top of this file) is a gProcNum x gProcNum (e.g., 16x16)
	// matrix that has entries (r,c) with the comm volume between two processes
	// (i,j) ranked with ranks r,c;
    procMatrix = new int*[gProcNum];
    for (int i=0; i<gProcNum; i++) 
        procMatrix[i] = new int[gProcNum];
    //fill it with corresponding value
    for (int i=0; i<gProcNum; i++) {
        int row = gProcess[i]->rank;
        for (int j=0; j<gProcNum; j++) {
            int col = gProcess[j]->rank;
            procMatrix[row][col] = gProcess[i]->fromComm[j] + gProcess[i]->toComm[j];
        }
    }
    //Sanity checking
#ifdef DEBUG
    for (int i=0; i<gProcNum; i++) {
        for (int j=0; j<gProcNum; j++) {
            if (procMatrix[i][j]<0) {
                cerr<<"Error for <0"<<endl;
                exit(1);
            }
            if (procMatrix[i][j]!=procMatrix[j][i]) {
                cerr<<"Error. The process matrix is not symetric."<<endl;
                exit(1);
            }
        }
    }
#endif //DEBUG
}

void BBMBuildArchitectureMatrix()
{
	// archMatrix (defined at the top of this file) is a gTileNum x gTileNum (e.g., 16x16)
	// matrix that has entries (src,dst) with the energy required to spend for
	// transfering a bit form tile src to tile dst;
    archMatrix = new float*[gTileNum];
    for (int i=0; i<gTileNum; i++) 
        archMatrix[i] = new float[gTileNum];
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            float energy = 0;
            pTile currentTile = &gTile[src];
            energy += currentTile->Cost();
            while (currentTile->GetId() != dst) {
                int linkId = currentTile->RouteToLink(src, dst);
                pLink pL = &gLink[linkId];
                energy += pL->Cost();
                currentTile = &gTile[pL->ToTile()];
                energy += currentTile->Cost();
            }
            archMatrix[src][dst] = energy;
        }
    }
    //Sanity checking
#ifdef DEBUG
    for (int i=0; i<gTileNum; i++) {
        for (int j=0; j<gTileNum; j++) {
            if (archMatrix[i][j]!=archMatrix[j][i]) {
                cerr<<"Error. The architecture matrix is not symetric."<<endl;
                exit(1);
            }
        }
    }
#endif //DEBUG
}

void BBMClear()
{
    cout<<"Clear for Branch-and-Bound mapping"<<endl;
    if (archMatrix) {
        for (int i=0; i<gTileNum; i++) 
            delete []archMatrix[i];
        delete []archMatrix;
    }
    if (procMatrix) {
        for (int i=0; i<gProcNum; i++) 
            delete []procMatrix[i];
        delete []procMatrix;
    }
    if (proc_map_array)
        delete []proc_map_array;
}

void BBM_Compute_normalization_factors()
{
	// Note: this has to be called only one time;
	// compute maximum possible energy and reliability costs; will be used
	// throughout calculations of bound costs and costs;
	// Note: COST_ENERGY_NORMALIZATION and COST_RELIABILITY_NORMALIZATION are 
	// initialized as 1;
	int count_comms = 0;
    for (int i=0; i<gProcNum; i++) {
        for (int j=i; j<gProcNum; j++) {
			COST_ENERGY_NORMALIZATION += procMatrix[i][j] * archMatrix[0][gProcNum-1];
			if ( procMatrix[i][j] > 0) {
				count_comms ++; // should we count ij and ji?
			}
        }
    }
	int delta_x = g_edge_size - 1, delta_y = g_edge_size - 1;
	COST_RELIABILITY_NORMALIZATION +=
		count_comms * reliability_cost_clone(delta_x,delta_y);
	printf("Compute normalization factors: %.2f %.2f \n",COST_ENERGY_NORMALIZATION,
		COST_RELIABILITY_NORMALIZATION);
}


void branchAndBoundMapping() 
{
	// this is the core engine of the branch-and-bound mapping;

    BBMInit();
	//reliability_costs_investigation(); exit(1); // debug;

    float minCost = MAX_VALUE;
    float minUpperBound = MAX_VALUE;
    pPQueue Q = new PQueue();

	// (1)
    if (exist_locked_pe()) {
        // this ruins the symmetric structure of the system completely. 
        // although for some corner cases, symmetry still exists, we don't 
        // consider it here. 
        for (int i=0; i<g_edge_size; i++) {
            for (int j=0; j<g_edge_size; j++) {
                pMappingNode pNode = new MappingNode(i*g_edge_size+j);
                if (pNode->is_illegal())
                    delete pNode;
                else
                    Q->insert(pNode);
            }
        }
    }
    else {
        // To exploit the symmetric structure of the system, we only need
        // to map the first processes to one corner of the chip, as shown
        // in the following code.
        /*****************************************************************
         * And if we need to synthesize the routing table, then there is
         * not much symmetry property to be exploited
         *****************************************************************/
        if (!param.routing_table_synthesis) {
            int size = (g_edge_size+1)/2;
            for (int i=0; i<size; i++) {
                for (int j=0; j<=i; j++) {
                    pMappingNode pNode = new MappingNode(i*g_edge_size+j);
                    if (pNode->is_illegal())
                        delete pNode;
                    else
                        Q->insert(pNode);
                }
            }
        }
        else {
            // for west-first or odd-even, we only need to consider the bottom half;
			// Note: these first half of nodes for the time being are only the nodes saying
			// that half of the application processes are each mapped on one of the 
			// first tiles; this is the "level one" of nodes in the search tree just
			// below the root node (e.g., process i mapped to tile j; "xxx..i..xxx");
			// each of these nodes will have children from the level-two of the search
			// tree (meaning two processes already mapped to two tiles); and so on;
            int size = (g_edge_size+1)/2;
            for (int i=0; i<size; i++) {
                for (int j=0; j<g_edge_size; j++) {
                    pMappingNode pNode = new MappingNode(i*g_edge_size+j);
                    if (pNode->is_illegal())
                        delete pNode;
                    else
                        Q->insert(pNode);
                }
            }
        }
    }

	// (2)
    pMappingNode bestMapping = NULL;
    int min_upperbound_hit_cnt = 0;
	
    while (!Q->empty()) {
		// get next node with smallest cost from the queue;
        pMappingNode pNode = Q->next();
        if (pNode->cost > minCost || pNode->lowerBound > minUpperBound) {
            delete pNode;
            continue;
        }

        bool insertAllFlag = false;
        int prev_insert = 0;
		// (a) regular functioning, insert better child nodes as they
		// are discovered;
        /**********************************************************************
         *Change this to adjust the tradeoff between the solution quality     *
         *and the run time                                                    *
         **********************************************************************/
        if (Q->Length() < param.pq_size) { 
        insertAll:
            if (pNode->upperBound == minUpperBound && minUpperBound < MAX_VALUE 
                && min_upperbound_hit_cnt <= param.min_hit_threshold) 
                insertAllFlag = true;
            for (int i=prev_insert; i<gTileNum; i++) {
                if (pNode->Expandable(i)) {
					// example: initially at level-one first half of the tiles
					// have a process locked to (or assigned) to them; here the
					// first expandable tile would be the first one from the 
					// remaining half; 
                    pMappingNode child = new MappingNode(*pNode, i);
                    if ( child->lowerBound>minUpperBound || child->cost>minCost
                        || (child->cost==minCost && bestMapping)
						|| child->is_illegal())
                        delete child;
                    else { // this is a child node worth pursuing further;
                        if (child->upperBound < minUpperBound) {
                            minUpperBound = child->upperBound;
                            cout << "Current minimum cost upper bound is (dia) " 
                                 << minUpperBound << " at stage " << child->stage << endl;
                            min_upperbound_hit_cnt = 0;

                            //some new stuff here: we keep the mapping with min upperBound
                            if (param.routing_table_synthesis) {
                                if (child->upperBound<minCost) {
                                    if (bestMapping)
                                        delete bestMapping;
                                    bestMapping = new MappingNode(*child);
                                    minCost = child->upperBound;
                                }
                                else if (child->upperBound<minUpperBound && !bestMapping) 
                                    bestMapping = new MappingNode(*child);
                            }
                        }
                        if (child->stage==gProcNum) {
							// have all processes been mapped yet?
                            minCost = child->cost;
                            if (minCost < minUpperBound)
                                minUpperBound = minCost;
                            cout<<"Current minimum cost is (dia) "<<minCost<<endl;
                            if (bestMapping)
                                delete bestMapping;
                            bestMapping = child;
                        }
                        else { // processes still exist to be mapped; add this child to queue;
                            Q->insert(child);
                            if (Q->Length()>=param.pq_size && !insertAllFlag) {
                                prev_insert = i;
                                goto selective_insert;
                            }
                        }
                    }
                }
            }
            continue;
        }
		// (b) when pqueue grows beyond some size;
        else { 
        selective_insert:
            if ((abs(pNode->upperBound-minUpperBound)<=0.01) && minUpperBound<MAX_VALUE 
                && min_upperbound_hit_cnt <= param.min_hit_threshold) {
                min_upperbound_hit_cnt++;
                goto insertAll;
            }
            // In this case, we only select one child which has the
            // smallest partial cost. However, if the node is currently
            // the one with the minUpperBound, then its child which
            // is generated by the corresponding minUpperBound is
            // also generated
            int index = -1;
			if ( !param.do_reliability) { // OLD
				index = pNode->bestCostCandidate();
			} else { // NEW
				index = pNode->bestCostCandidate_multiobjective();
			}
            pMappingNode child = new MappingNode(*pNode, index);
            if (child->lowerBound>minUpperBound || child->cost>minCost 
                || (child->cost==minCost && bestMapping)
                || child->is_illegal())
                delete child;
            else {
                if (child->upperBound < minUpperBound - 0.0) { // -0.01;
                    // In this case, we should also insert other children
                    delete child;
                    insertAllFlag = true;
                    goto insertAll;
                }
                if (child->stage==gProcNum || child->lowerBound==child->upperBound) {       
                    minCost = child->cost;
                    if (child->stage<gProcNum)
                        minCost = child->upperBound;
                    if (minCost<minUpperBound)
                        minUpperBound = minCost;
                    cout<<"Current minimum cost is <dsi1> "<<minCost<<endl;
                    if (bestMapping)
                        delete bestMapping;
                    bestMapping = child;
                } else {
                    Q->insert(child);
				}
            }

            if (pNode->upperBound > minUpperBound || pNode->upperBound==MAX_VALUE) {
                delete pNode;
                continue;
            }

            if (index == pNode->bestUpperBoundCandidate()) {
                delete pNode;
                continue;
            }

            index = pNode->bestUpperBoundCandidate();
#ifdef DEBUG
            if (!pNode->Expandable(index)) {
                cerr<<"Error in expanding at stage "<<pNode->Stage()<<endl;
                cerr<<"index = "<<index<<endl;
                exit(-1);
            }
#endif //DEBUG
            child = new MappingNode(*pNode, index);
            if (child->lowerBound>minUpperBound || child->cost>minCost)
                delete child;
            else {
                if (child->upperBound < minUpperBound) {
                    minUpperBound = child->upperBound;
                    cout<<"Current minimum cost upper bound is <dsi> "<<minUpperBound<<endl;
                    min_upperbound_hit_cnt = 0;
                }
                if (child->stage==gProcNum || child->lowerBound==child->upperBound) {
                    if (minCost==child->cost && bestMapping) 
                        delete child;
                    else {
                        minCost = child->cost;
                        if (child->stage<gProcNum)
                            minCost = child->upperBound;
                        if (minCost<minUpperBound)
                            minUpperBound = minCost;
                        cout<<"Current minimum cost is <dsi2> "<<minCost<<endl;
                        if (bestMapping)
                            delete bestMapping;
                        bestMapping = child;
                    }
                } else {
                    Q->insert(child);
				}
            }
        }
        delete pNode;
    }
    cout<<"Totally "<<MappingNode::cnt<<" nodes have been generated"<<endl;
    if ( bestMapping) {
        BBMMap( bestMapping);
        delete bestMapping;
    } else {
        cout<<"Can not find a suitable solution."<<endl;
	}
	
    BBMClear();
}

void BBMMap(pMappingNode bestMapping) {
    for (int i=0; i<gProcNum; i++) 
        gProcess[i]->mapToTile(-1);
    for (int i=0; i<gTileNum; i++) 
        gTile[i].mapToProc(-1);
    for (int i=0; i<gProcNum; i++) {
        int procId = proc_map_array[i];
        gProcess[procId]->mapToTile(bestMapping->mapToTile(i));
        gTile[bestMapping->mapToTile(i)].mapToProc(procId);
    }
    if (param.routing_table_synthesis)
        bestMapping->program_routers();
}

void optimizeMapping()
{
	// this is the main function which does the mapping; it creates the tree
	// for branch-and-bound;

    clock_t start, end;
    struct tms tmsstart, tmsend;
    start = times(&tmsstart);
    cout<<"Start mapping..."<<endl;

    assert(gProcNum == ((int) gProcess.size()));


    if ( param.optimization_method == opt_annealing) {
        anneal();
    } else {
        branchAndBoundMapping();
	}

	
    end = times(&tmsend);
    cout<<"Mapping process finished successfully."<<endl;
    //printTimes(end-start, &tmsstart, &tmsend);
}

////////////////////////////////////////////////////////////////////////////////
//
// more globals;
//
////////////////////////////////////////////////////////////////////////////////

void reliability_costs_investigation() 
{
	// debug purposes only;
	printf("\nReliability costs spread:");
	int edge_size = 4;
	float RC = 1000;
	float B = 1;
	int d = 0;
	int nx = edge_size-1;
	int w = 0, h = 0;
	double c = 0.0;
	for ( int dist = 1; dist <= nx + nx; dist ++) {
		d = dist / 2;
		B = RC / (d * (dist-1) + 1);
		printf("\ndist = %d", dist);
		if ( dist < 2) {
			w = 1; h = 0;
			c = 3*RC/4;
			printf("\n\t dx = %d, dy = %d  c = %.2f", w, h, c);
		} else {
			for ( w = 0; w <= dist/2; w ++) {
				h = dist - w;
				c = RC - B * w * h;
				printf("\n\t dx = %d, dy = %d  c = %.2f", w, h, c);
			}
		}
	}
}

double reliability_cost_clone( int dx, int dy)
{
	// this is an exact clone of MappingNode::reliability_cost(); I need it
	// also as a global 'coz I use it in BBM_Compute_normalization_factors();
	double rc = 0.0;
	int dist = dx + dy;
	double RC = 1000;
	if ( dist < 2) {
		return double(3*RC/4);
	}
	int d = dist / 2;
	double B = RC / (d * (dist-d) + 1);
	rc = RC - B * dx * dy;
	if ( dist > 2) {
		rc = rc + RC*(dist - 2);
	}
	return rc;
}

////////////////////////////////////////////////////////////////////////////////
//
// MAPPING_NODE
//
////////////////////////////////////////////////////////////////////////////////

double MappingNode::reliability_cost( int dx, int dy)
{
	// reliability cost for two processes mapped to tiles i,j which form
	// boundig box of width = dx and height = dy;
	double rc = 0.0;
	int dist = dx + dy;
	double RC = 1000;
	if ( dist < 2) {
		// this is to discourage somewhat the use of one hop connections;
		// I make the cost of dist=1 higher than of dist=2;
		return double(3*RC/4);
	}
	int d = dist / 2;
	double B = RC / (d * (dist-d) + 1);
	rc = RC - B * dx * dy;
	if ( dist > 2) {
		rc = rc + RC*(dist - 2);
	}
	return rc;
}
