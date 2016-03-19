#include "nocmap.h"
#include "gVal.h"
#include "mapping.h"
#include "param.h"
#include "msg.h"
#include <cstdio>
#include <iomanip>
#include <sys/times.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// utilities; used throughout the source files...
//
////////////////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream & os, const Position & pos) {
    os << "(" << pos.row << "," << pos.col << ")";
    return os;
}

ostream & operator<<(ostream & os, const Position3D & pos) {
    os << "(" << pos.pos.row << "," << pos.pos.col << ",";
    switch(pos.dir) {
    case NORTH:  os << "north)"; break;
    case EAST:   os << "east)"; break;
    case WEST:   os << "west)"; break;
    case SOUTH:  os << "south)"; break;
    case SOURCE: os << "source)"; break;
    case SINK:   os << "sink)"; break;
    default:     assert(0);
    }
    return os;
}

int toAbsAddr(const Position & pos) {
    return pos.row * g_edge_size + pos.col;
}

int toAbsAddr(const Position3D & pos) {
    return (toAbsAddr(pos.pos)*6 + pos.dir);
}

bool exist_locked_pe(void) {
    for (int i=0; i<gProcNum; i++) {
        if (gProcess[i]->is_locked())
            return true;
    }
    return false;
}

bool exist_non_regular_regions(void) {
    for (unsigned int i=0; i<gProcess.size(); i++) {
        if (gProcess[i]->get_pe_size() > 1)
            return true;
    }
    return false;
}

static bool spawn_dummy_process(void) 
{
	// spawn dummy processes for those processes whose PE has a size larger than 1
	// hmm, it never returns false;
    unsigned int real_proc_num = gProcess.size();
    for (unsigned int i=0; i<real_proc_num; i++) {
        if (gProcess[i]->get_pe_size() == 1) 
            continue;
        assert(gProcess[i]->get_pe_size() > 1);
        gProcess[i]->spawn_dummy_process();
    }
    return true;
}


// perform sanity check on the system. Return false if there is anything 
// wrong with the system
static bool check_system(void) {
    bool success = true;
    for (unsigned int i=0; i<gProcess.size(); i++) {
        if (!gProcess[i]->check_myself()) {
            cerr << "Process " << i << " fails sanity check." << endl;
            success = false;
        }
    }
    return success;
}

bool initialize(void) 
{
	// initialize the system, return false if the initialization encounts any error

    // check whether all the PEs can be fit in the chip
    int pe_sum_size = 0;
    for (unsigned int i=0; i<gProcess.size(); i++) 
        pe_sum_size += gProcess[i]->get_pe_size();
    if (pe_sum_size > gTileNum) {
        cerr << "Sorry, noc arch is too small " << gTileNum
			 << " to fit all the PEs " << pe_sum_size << endl;
        return false;
    }

    if (!spawn_dummy_process()) {
        cerr << "Error in spawning dummy process." << endl;
        return false;
    }

	// gLinkNum is *2 because links are in both direction for each channel;
    gLinkNum = 2 * (g_edge_size-1) * g_edge_size * 2;
    gLink = new Link[gLinkNum]();
    gTile = new Tile[gTileNum]();
    for (int i=0; i<gTileNum; i++) {
        gTile[i].AttachLink(gLink);
        gTile[i].initialize_router(param.routing_type);
    }

    if (param.verbose)
        cout << "Building link usage matrix." << endl;
    build_link_usage_matrix();

    return check_system();
}

void build_link_usage_matrix()
{
    if (param.link_usage_matrix) 
        delete_link_usage_matrix();

    // (1) Allocate the space for the link usage table

    param.link_usage_matrix = new int**[gTileNum];
    for (int i=0; i<gTileNum; i++) {
        param.link_usage_matrix[i] = new int*[gTileNum];
        for (int j=0; j<gTileNum; j++) 
            param.link_usage_matrix[i][j] = new int[gLinkNum];
    }

    for (int i=0; i<gTileNum; i++) 
        for (int j=0; j<gTileNum; j++) 
            for (int k=0; k<gTileNum; k++) 
                param.link_usage_matrix[i][j][k] = 0;

    // Setting up the link usage matrix
    for (int src_id=0; src_id<gTileNum; src_id++) {
        for (int dst_id=0; dst_id<gTileNum; dst_id++) {
            if (src_id == dst_id)
                continue;
            pTile current_tile = &gTile[src_id];
            while (current_tile->GetId() != dst_id) {
                int link_id = current_tile->RouteToLink(src_id, dst_id);
                pLink pL = &gLink[link_id];
                param.link_usage_matrix[src_id][dst_id][link_id] = 1;
                current_tile = &gTile[pL->ToTile()];
            }
        }
    }
  
    if (param.link_usage_list) 
        delete_link_usage_list();

    // (2) Now build the g_link_usage_list

    param.link_usage_list = new vector<int>*[gTileNum];
    for (int i=0; i<gTileNum; i++) 
        param.link_usage_list[i] = new vector<int>[gTileNum];
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            param.link_usage_list[src][dst].clear();
            if (src==dst) 
                continue;
            for (int link_id=0; link_id<gLinkNum; link_id++) {
                if (param.link_usage_matrix[src][dst][link_id]) 
                    param.link_usage_list[src][dst].push_back(link_id);
            }
        }
    }
}

void delete_link_usage_list() {
    if (param.link_usage_list) {
        for (int i=0; i<gTileNum; i++) 
            delete []param.link_usage_list[i];
        delete []param.link_usage_list;
    }
}

void delete_link_usage_matrix() {
    //free g_link_usage_matrix
    if (param.link_usage_matrix) {
        for (int i=0; i<gTileNum; i++) 
            for (int j=0; j<gTileNum; j++) 
                delete []param.link_usage_matrix[i][j];

        for (int i=0; i<gTileNum; i++) 
            delete []param.link_usage_matrix[i];
        delete []param.link_usage_matrix;
    }
}

void clean() {
    if (gProcess.size()) {
        for (vector<pProcess>::iterator iter = gProcess.begin(); iter < gProcess.end(); iter++) {
            delete *iter;
        }
        gProcess.clear();
    }
    if (!gLink) delete []gLink;
    if (!gTile) delete []gTile;

    delete_link_usage_matrix();
    delete_link_usage_list();
}

bool operator==(const Position & pos1, const Position & pos2) {
    return (pos1.row==pos2.row && pos1.col==pos2.col);
}

bool operator!=(const Position & pos1, const Position & pos2) {
    return !(pos1 == pos2);
}

bool operator==(const Position3D & pos1, const Position3D & pos2) {
    if (pos1.pos != pos2.pos)
        return false;
    if (pos1.dir != pos2.dir)
        return false;
    return true;
}

bool operator!=(const Position3D & pos1, const Position3D & pos2) {
    return !(pos1 == pos2);
}

bool testRoutability(Tile& srcTile, Tile& dstTile) {
    if (param.verbose) {
        cout << "Testing routablity from Tile" << srcTile.GetId() << " ("
             << srcTile.GetPosition().row << "," << srcTile.GetPosition().col
             << ") to Tile" << dstTile.GetId() << " ("
             << dstTile.GetPosition().row << "," << dstTile.GetPosition().col
             << "):" << endl;
    }
    if (param.verbose)
        cout << "Tile" << srcTile.GetId() << "(" << srcTile.GetPosition().row
             << "," << srcTile.GetPosition().col << ")";
    pTile currentTile = &srcTile;
    int srcProc = srcTile.mapToProc();
    int dstProc = dstTile.mapToProc();
    if (!gProcess[srcProc]->ToComm(dstProc)) {
        if (param.verbose)
            cout << "no data flow ... Done!"<<endl;
        return true;
    } 
    int hopNum = 0;
    while (currentTile->GetId() != dstTile.GetId()) {
        int linkId = currentTile->RouteToLink(srcTile.GetId(), dstTile.GetId());
        pLink pL = &gLink[linkId];
        if (param.verbose) 
            cout << " using link " << pL->GetId() << " --> ";
        currentTile = &gTile[pL->ToTile()];
        if (param.verbose)
            cout << "Tile" << currentTile->GetId() << "("<<currentTile->GetPosition().row
                 << "," << currentTile->GetPosition().col << ")";
        hopNum ++;
        if (hopNum == gTileNum)
            break;
    }

    if (hopNum == gTileNum) {
        if (param.verbose) 
            cout << "....Fail!"<<endl;
        return false;
    }
    else {
        if (param.verbose)
            cout << "....Done!"<<endl;
        return true;
    }
}

bool verify() {
    cout << "Total tile number: " << gTileNum << endl;
    cout << "Number of tiles in one row/column: " << g_edge_size << endl;
    cout << "Total link number: "<<gLinkNum<<endl;
    for (int i=0; i<gTileNum; i++) 
        for (int j=0; j<gTileNum; j++) 
            if (!testRoutability(gTile[i], gTile[j])) {
                cout << "Not routable from tile " << i << "to tile " << j << endl;
                return false;
            }
    return true;
}

bool verify_BW_requirement() {
    if (!param.link_usage_matrix) 
        build_link_usage_matrix();

    for (int i=0; i<gLinkNum; i++) 
        gLink[i].used_BW = 0;

    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            if (src == dst)
                continue;
            int src_proc = gTile[src].mapToProc();
            int dst_proc = gTile[dst].mapToProc();
            int comm_load = gProcess[src_proc]->to_BW_requirement[dst_proc];
            if (comm_load == 0)
                continue;
            pTile current_tile = &gTile[src];
            while (current_tile->GetId() != dst) {
                int link_id = current_tile->RouteToLink(src, dst);
                pLink pL = &gLink[link_id];
                current_tile = &gTile[pL->ToTile()];
                gLink[link_id].used_BW += comm_load;
            }
        }
    }
    //check for the overloaded links
    int violations = 0;
    for (int i=0; i<gLinkNum; i++) {
        if (gLink[i].used_BW > gLink[i].bandwidth) {
            cout << "Link " << i << " is overloaded: " << gLink[i].used_BW << " > "
                 << gLink[i].bandwidth << endl;
            violations ++;
        }
    }
    if (violations)
        return false;
    return true;
}


bool mapProcToTile(int type)
{
    if (gProcNum != gTileNum || gProcNum != ((int) gProcess.size())) {
        cerr << "Process number and tile number do not match." << endl;
        cerr << "Mapping failed" << endl;
        return false;
    }

    //First clear the old mapping.
    for (int i=0; i<gTileNum; i++) {
        gProcess[i]->mapToTile(-1);
        gTile[i].mapToProc(-1);
    }
    cout << "Mapping processes to the tile ";
    cout << "(type: ";
  
    switch(type) {
    case SEQUENTIAL_MAPPING:
        cout << "SEQUENTIAL) ..."<<endl; 
        for (int i=0; i<gTileNum; i++) {
            gProcess[i]->mapToTile(i);
            gTile[i].mapToProc(i);
        }
        break;
    case RANDOM_MAPPING:
        cout << "RANDOM) ..."<<endl;
        for (int i=0; i<gTileNum; i++) {
            int k = rand()%gTileNum;
            while (gTile[k].mapToProc() != -1) 
                k = rand()%gTileNum;
            gProcess[i]->mapToTile(k);
            gTile[k].mapToProc(i);
        }
        break;
    default:
        cout << "UNSUPPORTED!) ..."<<endl;
        cout << "quit ..."<<endl;
        return false;
    }
    return true;
}


bool parse_options(int argn, char ** argv)
{
	// return false if error 

    if (argn == 1) {
        print_help();
        param.quit_flag = true;
        return true;
    }

    //prcessing command line arguments
    for (int i=1; i<argn; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            print_help();
            param.quit_flag = true;
            return true;
        }

        if (strcmp(argv[i], "-seed") == 0) {
            param.seed = atoi(argv[++i]);
            srand(param.seed);
            continue;
        }

        if (strcmp(argv[i], "-PQ_size") == 0) {
            sscanf(argv[++i], "%d", &param.pq_size);
            cout << "Maximum priority queue length set up to " << param.pq_size << endl;
            continue;
        }

        if (strcmp(argv[i], "-bandwidth") == 0) {
            sscanf(argv[++i], "%d", &param.link_bandwidth);
            cout << "Link bandwidth is set to " << param.link_bandwidth << endl;
            continue;
        }

        if (strcmp(argv[i], "-link_ebit") == 0) {
            sscanf(argv[++i], "%f", &param.link_ebit);
            cout << "Link power coefficient set to " << param.link_ebit << endl;
            continue;
        }

        if (strcmp(argv[i], "-switch_ebit") == 0) {
            sscanf(argv[++i], "%f", &param.switch_ebit);
            cout << "Switch power coefficient set to " << param.switch_ebit << endl;
            continue;
        }

        if (strcmp(argv[i], "-buffer_read_ebit") == 0) {
            sscanf(argv[++i], "%f", &param.buf_read_ebit);
            continue;
        }

        if (strcmp(argv[i], "-buffer_write_ebit") == 0) {
            sscanf(argv[++i], "%f", &param.buf_write_ebit);
            continue;
        }

        if (strcmp(argv[i], "-tilenum") == 0) {
            gTileNum = atoi(argv[++i]);
            g_edge_size = (int) sqrt((double) gTileNum);
            if (g_edge_size * g_edge_size != gTileNum) {
                cerr << ERRO_NETWORK_SIZE;
                param.quit_flag = true;
                return false;
            }
            continue;
        }

        if (strcmp(argv[i], "-routing_xy") == 0) {
            param.routing_type = routing_xy;
            continue;
        }

        if (strcmp(argv[i], "-routing_oddeven") == 0) {
            param.routing_type = routing_oddeven;
            continue;
        }

        if (strcmp(argv[i], "-routing_effort")  ==  0) {
            i++;
            if (strcmp(argv[i], "hard")  ==  0) 
                param.routing_effort = hard;
            else if (strcmp(argv[i], "easy")  ==  0)
                param.routing_effort = easy;
            else {
                param.quit_flag = true;
                param.parse_error = true;
                return false;
            }
            continue;
        }

        if (strcmp(argv[i], "-me") == 0) {
            param.do_mapping = true;
			// if only energy is the mapping objective, by default no routing
			// adaptivity is assumed; this way the old tool behaviour is maintained;
			param._tables_update = 0; // false;
            continue;
        }
        if (strcmp(argv[i], "-mr") == 0) {
            param.do_reliability = true;
			// if mapping is done also for reliability optimization, then 
			// minimal adaptivity (via adaptive routing) is assumed to exist 
			// or be provided by the NoC architecture; in this case the reliability
			// engine is selected via tables_update;
			param._tables_update = 1; // true;
            continue;
        }
        if (strcmp(argv[i], "-alpha") == 0) {
            sscanf(argv[++i], "%lf", &param.alpha);
			// if alpha is set to zero, it means that user wants to do only
			// energy oriented mapping but with the cost component normalized
			// as in the case of energy and reliability oriented mapping; in this
			// case the reliability cost component will still be computed but
			// basically not used;
			if ( param.alpha <= 0.0) {
				param.alpha = 1e-9;
				param._tables_update = 0; // false;
			}
			cout << "Alpha weight parameter set to " << param.alpha << endl;
			continue;
        }
        if (strcmp(argv[i], "-q") == 0) {
            sscanf(argv[++i], "%lf", &param._q);
            continue;
        }
        if (strcmp(argv[i], "-M") == 0) {
            sscanf(argv[++i], "%d", &param._M);
            continue;
        }
		// tables_update is used to directly control what reliability engine
		// is used irrespective of the type of mapping objective(s);
        if (strcmp(argv[i], "-tables_update") == 0) {
            sscanf(argv[++i], "%d", &param._tables_update);
            continue;
        }

        if (strcmp(argv[i], "-synthesis_routing") == 0) {
            param.routing_table_synthesis = true;
            i++;
            if (strcmp(argv[i], "odd_even") == 0) 
                param.legal_turn_set = odd_even;
            else if (strcmp(argv[i], "west_first") == 0) 
                param.legal_turn_set = west_first;
            else {
                param.parse_error = true;
                return -1;
            }
            continue;
        }

        if (strcmp(argv[i], "-annealing") == 0) {
            param.optimization_method = opt_annealing;
            continue;
        }

        if (strcmp(argv[i], "-dump_mapping") == 0) {
            param.dump_mapping = true;
            param.dump_mapping_table_file = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "-read_mapping") == 0) {
            param.read_mapping = true;
            param.read_mapping_table_file = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "-dump_apcg") == 0) {
            param.dump_apcg = true;
            param.dump_apcg_file = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "-dump_routing") == 0) {
            param.dump_routing = true;
            param.routing_table_file = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "-apcg") == 0) {
            if (gProcess.empty()) {
                gProcNum = gTileNum;
                for (int i=0; i<gProcNum; i++) {
                    pProcess p = new Process();
                    gProcess.push_back(p);
                }
            }
            if (!parse_apcg(argv[++i])) {
                cerr << ERRO_APCG_FILE << argv[i] << endl;
                return false;
            }
            continue;
        }
		// () read-in traffic file;
        if (strcmp(argv[i], "-traffic_config") == 0) {
            if (gProcess.empty()) {
                gProcNum = gTileNum;
                for (int i=0; i<gProcNum; i++) {
                    pProcess p = new Process();
                    gProcess.push_back(p);
                }
            }
            if (!parse_traffic_config(argv[++i])) {
                cerr << "Error in parsing input APCG file " << argv[i] << endl;
                return -1;
            }
            continue;
        }

        if (strcmp(argv[i], "-v") == 0) {
            param.verbose = true;
            continue;
        }

        if (strcmp(argv[i], "-V") == 0) {
            cout << "Version " << 1.2 << endl;
            param.quit_flag = true;
            return false;
        }

        cerr << ERRO_UNRECOGNIZED_OPTION << argv[i] << endl;
        param.parse_error = true;
        param.quit_flag = true;
        return false;
    }

    return true;
}


void print_help() {
    cout << "Usage: " << endl;
    cout << "-h:" << endl;
    cout << "\tprint out this help." << endl;
    cout << "-V:" << endl;
    cout << "\tprint out the version." << endl;
    cout << "-v:" << endl;
    cout << "\trun in verbose mode." << endl;
    cout << "-tilenum tile_number:" << endl;
    cout << "\tspecify the number of tiles in the architecture." << endl;
    cout << "-seed seed:" << endl;
    cout << "\tinitialize the random number generator using seed." << endl;
    cout << "-packet_size size:" << endl;
    cout << "\tset the packet size in bit. (default = 64bit)" << endl;
    cout << "-bandwidth BW:" << endl;
    cout << "\tset the bandwidth of links to BW. (default = 1000000)" << endl;
    cout << "-router_buffer router_buffer:" << endl;
    cout << "\tset the buffer size (# of pkts) in each router's input port. (default = 2 pkts)" << endl;
    cout << "-link_ebit link_ebit:" << endl;
    cout << "\tspecify the ebit of the link with unit length. (default = 1)" << endl;
    cout << "-switch_ebit switch_ebit:" << endl;
    cout << "\tspecify the ebit of the switch. (default = 1)" << endl;
    cout << "-buffer_read_ebit buffer_read_ebit:" << endl;
    cout << "\tspecify the buffer read energy consumption per bit." << endl;
    cout << "-buffer_write_ebit buffer_write_ebit:" << endl;
    cout << "\tspecify the buffering energy consumption per bit." << endl;
    cout << "-routing_xy:" << endl;
    cout << "\tuse XY routing." << endl;
    cout << "-routing_oddeven:" << endl;
    cout << "\tuse odd-even routing." << endl;
    cout << "-annealing:" << endl;
    cout << "\twhen used with -me option, it uses simulated annealing to find the mapping." << endl;
    cout << "-synthesis_routing odd_even|west_first:" << endl;
    cout << "\twhen combined with -me option, it also generates the routing table in" << endl
         << "\tstead of limiting itself to xy routing. " << endl
         << "\tChoose between odd_even or west_first for corresponding legal turn set" << endl;
    cout << "-PQ_size size:" << endl;
    cout << "\tthe maximum allowed priority queue length (default = 2000)" << endl;
    cout << "-dump_mapping file_name:" << endl;
    cout << "\tsave the generated mapping to the file specified by file_name." << endl;
    cout << "-read_mapping file_name:" << endl;
    cout << "\tread the mapping from the file specified by file_name." << endl;
    cout << "-dump_apcg file_name:" << endl;
    cout << "\tsave the internal APCG to the file specified by file_name." << endl;
    cout << "-dump_routing file_name:" << endl;
    cout << "\tsave the routing table to the file specified by file_name." << endl;
    cout << "-apcg file_name:" << endl;
    cout << "\tuse apcg format as system input." << endl;
    cout << "-traffic_config file_name:" << endl;
    cout << "\tuse traffic config format (which is used in worm_sim) as system input." << endl;
    cout << "\tbandwidth constraints will be set up proportional to the rate." << endl;
    cout << "-me:" << endl;
    cout << "\tenergy oriented mapping to map processes to tiles" << endl;
    cout << "-mr:" << endl;
    cout << "\tinclude reliability also into the cost function" << endl;
    cout << "-alpha:" << endl;
    cout << "\tweight of reliability cost component; used with -mr only; default is 0.6" << endl;
    cout << "-q:" << endl;
    cout << "\tprobability link is in 'down' state; default is 0.01" << endl;
    cout << "-M:" << endl;
    cout << "\tnumber of runs inside the Monte Carlo reliability engine; default is 10k" << endl;

}

float calc_comm_energy() {
    float sw_energy = calc_switch_energy();
    float link_energy = calc_link_energy();
    float buffer_energy = calc_buffer_energy();
    return (sw_energy + link_energy + buffer_energy);
}

float calc_switch_energy() {
    float energy = 0.;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue;  // no go for vitual volume of dummy PE
            if (!commVol)
                continue;
            pTile currentTile = &gTile[src];
            energy += currentTile->Cost() * commVol;
            while (currentTile->GetId() != dst) {
                int linkId = currentTile->RouteToLink(src, dst);
                pLink pL = &gLink[linkId];
                currentTile = &gTile[pL->ToTile()];
                energy += currentTile->Cost() * commVol;
            }
        }
    }
    return energy;
}

float calc_link_energy() {
    float energy = 0.;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue;  // no go for vitual volume of dummy PE
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (!commVol)
                continue;
            pTile currentTile = &gTile[src];
            while (currentTile->GetId() != dst) {
                int linkId = currentTile->RouteToLink(src, dst);
                pLink pL = &gLink[linkId];
                energy += pL->Cost() * commVol;
                currentTile = &gTile[pL->ToTile()];
            }
        }
    }
    return energy;
}

float calc_buffer_energy() {
    float energy = 0.;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue;  // no go for vitual volume of dummy PE
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (!commVol)
                continue;
            pTile currentTile = &gTile[src];
            while (currentTile->GetId() != dst) {
                int linkId = currentTile->RouteToLink(src, dst);
                pLink pL = &gLink[linkId];
                energy += (param.buf_read_ebit + param.buf_write_ebit) * commVol;
                currentTile = &gTile[pL->ToTile()];
            }
            energy += param.buf_write_ebit * commVol; // 2.831 * commVol;
        }
    }
    return energy;
}


int calc_total_wirelength()
{
	// computation goes hop-by-hop; this will be useful when the routing
	// will not be minimal but could be adaptive as well;
    int wl = 0;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue; // no go for vitual volume of dummy PE
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (!commVol)
                continue;
            pTile currentTile = &gTile[src];
            while (currentTile->GetId() != dst) {
                int linkId = currentTile->RouteToLink(src, dst);
                pLink pL = &gLink[linkId];
                wl ++;
                currentTile = &gTile[pL->ToTile()];
            }
        }
    }
    return wl;
}
int calc_total_area()
{
	// computation of area of all bboxes;
    int area = 0;
	int delta_x = 0, delta_y = 0;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue; // no go for vitual volume of dummy PE
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (!commVol)
                continue;
			delta_x = abs( gTile[src].GetPosition().col -  gTile[dst].GetPosition().col);
			delta_y = abs( gTile[src].GetPosition().row -  gTile[dst].GetPosition().row);
			area += (delta_x * delta_y);
        }
    }
    return area;
}
double calc_total_reliability_cost()
{
	// utilized only by the Simulated Annealing algo;
    double cost_r = 0.0;
	int delta_x = 0, delta_y = 0;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue; // no go for vitual volume of dummy PE
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (!commVol)
                continue;
			delta_x = abs( gTile[src].GetPosition().col -  gTile[dst].GetPosition().col);
			delta_y = abs( gTile[src].GetPosition().row -  gTile[dst].GetPosition().row);
			cost_r += reliability_cost_clone( delta_x, delta_y);
			// I do not like the above global function; TO DO: re-architect the whole code;
        }
    }
    return cost_r;
}

void Anneal_Compute_normalization_factors()
{
	// Note: this has to be called only one time; used only when Simulated
	// Anealing algorithm is utilized;
	// Note: COST_ENERGY_NORMALIZATION and COST_RELIABILITY_NORMALIZATION are 
	// initialized as 1;
	int delta_x = g_edge_size - 1, delta_y = g_edge_size - 1;
	int dist = delta_x + delta_y, dist_plus_1 = dist + 1;
	double cost_e = 0.0;
	double cost_r = 0.0;
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            int srcProc = gTile[src].mapToProc();
            int dstProc = gTile[dst].mapToProc();
            int commVol = gProcess[srcProc]->ToComm(dstProc);
            if (gProcess[srcProc]->is_dummy() || gProcess[dstProc]->is_dummy())
                continue; // no go for vitual volume of dummy PE
            if (!commVol)
                continue;
			// assume the worst case: as if the two cores would be mapped in
			// opposite corners of the NoC architecture;
			pTile currentTile = &gTile[src];
			int linkId = currentTile->RouteToLink(src, dst);
			pLink pL = &gLink[linkId];
			// energy;
			cost_e += dist_plus_1 * currentTile->Cost() * commVol; // switch;
			cost_e += dist * pL->Cost() * commVol; // link;
			cost_e += dist_plus_1 * (param.buf_read_ebit + param.buf_write_ebit) * commVol; // buffer;
			// reliability;
			//printf(" (%d, dist %d)", commVol, dist);
			cost_r += 1 * reliability_cost_clone(delta_x,delta_y);
        }
    }

	COST_ENERGY_NORMALIZATION += cost_e;
	COST_RELIABILITY_NORMALIZATION += cost_r;

	printf("Compute normalization factors: %.2f %.2f \n",COST_ENERGY_NORMALIZATION,
		   COST_RELIABILITY_NORMALIZATION);
}

void analyzeIt()
{
    cout << "Verify the communication load of each link...";
    if (verify_BW_requirement()) {
        cout << "Succeed." << endl;
    } else {
        cout << "Fail." << endl;
    }
	if (param.verbose) {
		cout << "Total: WL is " << calc_total_wirelength() << " bbox area is " << calc_total_area() << endl;
	}	
    cout << "Energy consumption estimation " << endl
         << "(note that this is not exact numbers, but serve as a relative energy indication) " << endl;
    cout << "Energy consumed in link is " << calc_link_energy() << endl;
    cout << "Energy consumed in switch is " << calc_switch_energy() << endl;
    cout << "Energy consumed in buffer is " << calc_buffer_energy() << endl;
    cout << "Total communication energy consumption is " << calc_comm_energy() << endl;
}

void dispMapping()
{
	// very important: this new way of reporting the mapping result is
	// according to how it is dumped into the out file; please see
	// dump_mapping_table() for more details on this; 
	cout << "Mapping:  process i  mapped to tile j" << endl;
    for (int i=0; i<gProcNum; i++) {
        int proc_id = i;
        if (gProcess[i]->is_dummy()) {
            pProcess master = gProcess[i]->get_master();
            proc_id = master->GetId();
        }
        int tid = gProcess[i]->mapToTile();
        int row = tid / g_edge_size;
        int col = tid % g_edge_size;
        tid = col * g_edge_size + row;
		cout << proc_id << "\t" << tid << endl;
    }
	// the original version:
    //for (int i=0; i<gProcNum; i++) {
    //    int tileId = gProcess[i]->mapToTile();
    //    Position pos = gTile[tileId].GetPosition();
    //    cout << "Process "<<i<<" mapped to tile ("<<pos.row<<","<<pos.col<<")"<<endl;
    //}
}

void printTimes(clock_t real, struct tms *tmsstart, struct tms *tmsend) {
    static long clktck = 0;
    if (clktck == 0) {
        if ((clktck = sysconf(_SC_CLK_TCK)) < 0) 
            cerr<<"sysconf error"<<endl;
    }
    cout << "Used time:" << endl;
    cout << "\tReal:\t" << real/((double) clktck) << endl;
    cout << "\tUser:\t" << (tmsend->tms_utime - tmsstart->tms_utime)/(double) clktck << endl;
    cout << "\tSys:\t" << (tmsend->tms_stime-tmsstart->tms_stime)/(double) clktck << endl;
    cout << "\tCPU:\t"
		 <<(tmsend->tms_utime-tmsstart->tms_utime + tmsend->tms_stime-tmsstart->tms_stime)/(double) clktck
		 <<endl;
}

bool dump_mapping_table(char * file_name) {
    remove(file_name);
    ofstream out_file(file_name);
    if (!out_file) {
        remove(file_name);
        cerr << "Can not open file " << file_name << " for writing!" << endl;
        return false;
    }
    cout << "Dumping the mapping talbe to file " << file_name << "..." << endl;

    /************************************************************
     * The format of the file mapping_table                    
     * conatins two columns, the first one is the IP id,  
     * The second column contains the tile id.                  
     ************************************************************/
    for (int i=0; i<gProcNum; i++) {
        int proc_id = i;
        if (gProcess[i]->is_dummy()) {
            pProcess master = gProcess[i]->get_master();
            proc_id = master->GetId();
        }
        // FIXME (2005-01-22) the id calculation in nocmap and worm_sim/metropolis noc
        // is not consistent. here id of a tile at (col, row) is calculated as 
        // row*g_edge_size + col, while in worm_sim the id of the router at (x,y)
        // is calculated by x*g_edge_size + y. 
        // To make this result useable for worm_sim and metropolis noc simulation, I temporarily
        // modified it. However, this is very ugly and can potentially cause lots of problems:
        // e.g. read_mapping_table function should also changed, etc. 
        // The final goal is to modify nocmap throughout so the address representation can be in
        // coherent w/ worm_sim
        //        out_file << proc_id << "\t" << gProcess[i]->mapToTile() << endl;
        int tid = gProcess[i]->mapToTile();
        int row = tid / g_edge_size;
        int col = tid % g_edge_size;
        tid = col * g_edge_size + row;
        out_file << proc_id << "\t" << tid << endl;
    }

    out_file.close();
    return true;
}

bool dump_apcg(char * file_name) {
    remove(file_name);
    ofstream outFile(file_name);
    if (!outFile) {
        remove(file_name);
        cerr << "Can not open file " << file_name << " for write" << endl;
        return false;
    }
    cout << "Dumping the APCG to file " << file_name << "..." << endl;

    outFile << "#Total IPs " << gProcNum << endl;
    outFile << "#Lines starting with # are comments." << endl;
    outFile << "#Format:" << endl;
    outFile << "#SrcIP\tSinkIP\tCommVol\tCommBW" << endl;

    for (int src=0; src<gProcNum; src++) {
        for (int dst=0; dst<gProcNum; dst++) {
            if (gProcess[src]->toComm[dst] == 0)
                continue;
            outFile << src << "\t" << dst << "\t" << gProcess[src]->toComm[dst]
                    << "\t" << gProcess[src]->to_BW_requirement[dst] << endl;
        }
    }

    outFile.close();
    return true;
}


bool read_mapping_table(char* file_name) {
    char input_line[MAX_LINE];
    int proc_id, tile_id;
    ifstream mapping_file(file_name);
    if (!mapping_file) {
        cerr << "Error in openning " << file_name << endl;
        return false;
    }
    //delete the original mapping, if any
    for (int i=0; i<gProcNum; i++) {
        gTile[i].mapToProc(-1);
        gProcess[i]->mapToTile(-1);
    }
    int line = 0;
    while (mapping_file.getline(input_line, MAX_LINE, '\n')) {
        if (!sscanf(input_line, "%d\t%d", &proc_id, &tile_id))
            return false;
        gProcess[proc_id]->mapToTile(tile_id);
        gTile[tile_id].mapToProc(proc_id);
        line ++;
    }
    if (line!=gProcNum) {
        cerr << "Error: Line number and the process number does not match" << endl;
        return false;
    }
    return true;
}

//find out the link_id. if the direction is not set, return -1
int locate_link(Position & pos, int direction) {
    Position next_hop = pos;
    switch(direction) {
    case NORTH:
        next_hop.row++;
        break;
    case SOUTH:
        next_hop.row--;
        break;
    case EAST:
        next_hop.col++;
        break;
    case WEST:
        next_hop.col--;
        break;
    default:
        return -1;
    }
    int link_id;
    for (link_id=0; link_id<gLinkNum; link_id++) {
        if (gLink[link_id].fromTile==pos && gLink[link_id].toTile==next_hop)
            break;
    }
    if (link_id == gLinkNum) {
        cerr<<"Error in locating link";
        exit(-1);
    }
    return link_id;
}


int * generate_random_array(int array_size, int min, int max) {
    if (min>=max||array_size>(max-min+1)) {
        cerr<<"Error: wrong parameters"<<endl;
        return NULL;
    }
    int *array = new int[array_size];
    for (int i=0; i<array_size; i++) {
        while (true) {
            float range = 1.0 + max - min;
            int r = min + (int) (range*rand()/(RAND_MAX+1.0));
            //check whether this number is already in the array
            int cnt;
            for (cnt=0; cnt<i; cnt++) {
                if (array[cnt] == r)
                    break;
            }
            if (cnt==i) {
                array[cnt] = r;
                break;
            }
        }
    }
    return array;
}


void print_comm_statistics() {
    cout << "Communication statistics:"<<endl;
    for (int pe=0; pe<gProcNum; pe++) {
        cout << "For PE "<<pe<<":"<<endl;
        cout << "To communication"<<endl;
        for (int i=0; i<gProcNum; i++) {
            cout << gProcess[pe]->toComm[i]<<"\t";
            if (gProcess[pe]->toComm[i] != gProcess[i]->fromComm[pe]) 
                cerr<<"Error in checking comm"<<endl;
        }
        cout << endl;
        cout << "To communication BW requirement"<<endl;
        for (int i=0; i<gProcNum; i++) {
            cout << gProcess[pe]->to_BW_requirement[i]<<"\t";
            if (gProcess[pe]->to_BW_requirement[i] != gProcess[i]->from_BW_requirement[pe]) 
                cerr<<"Error in checking comm"<<endl;
        }
        cout << endl;

        cout << "From communication"<<endl;
        for (int i=0; i<gProcNum; i++) {
            cout << gProcess[pe]->fromComm[i]<<"\t";
            if (gProcess[pe]->fromComm[i] != gProcess[i]->toComm[pe]) 
                cerr<<"Error in checking comm"<<endl;
        }
        cout << endl;
        cout << "From communication BW requirement"<<endl;
        for (int i=0; i<gProcNum; i++) {
            cout << gProcess[pe]->from_BW_requirement[i]<<"\t";
            if (gProcess[pe]->from_BW_requirement[i] != gProcess[i]->to_BW_requirement[pe]) 
                cerr<<"Error in checking comm"<<endl;
        }
        cout << endl;
    }
}

bool clear_routing_tables(void) {
    for (int tId=0; tId<gTileNum; tId++) 
        for (int srcId=0; srcId<gTileNum; srcId++) 
            for (int dstId=0; dstId<gTileNum; dstId++) 
                gTile[tId].RouteToLink(srcId, dstId, -2);
    return false;
}


bool dump_routing_tables(char * file_name) {
    remove(file_name);
    ofstream out_file(file_name);
    if (!out_file) {
        remove(file_name);
        cerr << "Can not open file " << file_name << " for write" << endl;
        return false;
    }
    cout << "Dumping the routing table to file " << file_name << "..." << endl;

    out_file << "#Total tiles " << gTileNum << endl;
    out_file << "#Lines starting with # are comments." << endl;
    out_file << "#Format:" << endl;
    out_file << "#Keyword Tile begins a new tile's routing table, followed by the tile id" << endl;
    out_file << "#Each entry is in the format of SrcID\tSinkID\tNextHopID" << endl;

    for (int tId=0; tId<gTileNum; tId++) {
        out_file << endl << endl;
        out_file << "Tile\t" << tId << "\t# Routing table for tile " << gTile[tId].GetPosition() << endl;
        for (int srcId=0; srcId<gTileNum; srcId++) 
            for (int dstId=0; dstId<gTileNum; dstId++) {
                int linkId = gTile[tId].RouteToLink(srcId, dstId);
                if (linkId == -2)
                    continue;
                if (dstId == srcId) 
                    continue;
                if (dstId == tId) 
                    continue;
                out_file << srcId << "\t" << dstId << "\t" << gLink[linkId].ToTile() << endl;
            }
    }

    out_file.close();
    return true;
}


////////////////////////////////////////////////////////////////////////////////
//
// NOCMAP
//
////////////////////////////////////////////////////////////////////////////////
