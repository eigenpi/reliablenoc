#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/times.h>
#include "app.h"
#include "nocmap.h"
#include "arch.h"
#include "mapping.h"
#include "param.h"
#include "msg.h"
#include "reliability.h"

using namespace std;

// the next ones are globals that are worked-on/accessed from a lot of files;
// TO DO: make them variables of NOCMAP class for clean coding; these
// globals e o sursa de probleme; ar trebui sa scap de ele insa e mult de
// munca, doar de restructurat code...
Param param;
pLink gLink = NULL;
pTile gTile = NULL;
vector<pProcess> gProcess;
int gLinkNum;
int gProcNum;
int gTileNum = 16;
int g_edge_size = 4; // how many tiles per edge

// I should get rid of these statics too;
int Process::proc_num = 0;
int Link::cnt = 0;
int Tile::cnt = 0;
int Tile::totalTiles = 0;
int MappingNode::cnt = 0;


////////////////////////////////////////////////////////////////////////////////
//
// launching point;
//
////////////////////////////////////////////////////////////////////////////////

int test_single_path();

int main( int argn, char** argv)
{
	// () cpu time;
    clock_t start, end;
    struct tms tmsstart, tmsend;
    start = times(&tmsstart);


    srand( param.seed); // srand( time(NULL));
    gProcNum = 16;
    gTileNum = 16;
    gProcess.clear();

    strcpy( param.exec_str, argv[0]);
    parse_options( argn, argv); // works on globals such as param... :(
	                            // create also gProcess;


    if ( param.parse_error) {
        cerr << ERRO_PARSING_OPTIONS;
        return -1;
    }

    if ( param.quit_flag)
        return 0;
	

    // () By calling initialize, classes are constructed, with some randomly
    // generated communication pattern and the default routing table. 
    // The processes parameters can be overloaded by calling parseTGFile()
    // later, and the routing table can be later be override by calling
    // parseRoutingTable()
    if ( !initialize()) {
        cerr << ERRO_SYSTEM_INIT;
        return -1;
    }

	// () assign randomly application's processes to tiles; part of NOCMAP methods;
	// Note: nocmap works with a number of processes same as the number of
	// tiles; if application has less processes, then we err out; this should
	// be corrected;
    mapProcToTile(RANDOM_MAPPING);

	// () optimizeMapping is member of MAPPING methods;
    if ( param.do_mapping) 
        optimizeMapping();

	// () print's and save's;
    if ( param.read_mapping) 
        read_mapping_table(param.read_mapping_table_file);
    if ( param.dump_mapping) 
        dump_mapping_table(param.dump_mapping_table_file);
    if ( param.verbose) {
		dispMapping(); // print the final mapping;
	}

    if ( param.dump_apcg)
        dump_apcg(param.dump_apcg_file);
    if ( param.dump_routing)
        dump_routing_tables(param.routing_table_file);

	// () sanity check;
    analyzeIt();


	// () cpu time - taken by mapping;
	if ( param.verbose) {
		end = times(&tmsend);
		printTimes(end-start, &tmsstart, &tmsend);
	}
	

	// () create light weight noc_architecture and build the reliability engine;
	//test_single_path(); exit(1); // debug;
	NOC_ARCHITECTURE noc_architecture;
	noc_architecture.create_nodes_and_arcs();
	//noc_architecture.print_noc_architecture();
    if ( param._tables_update == 1) {
		// this is the case advocated in the paper: assume minimal adaptivity
		// (supported via minimal adaptive routing) which does link failure 
		// detection and routing tables update to set new healthy routing paths 
		// via "monotonic XY-YX" routing;
		noc_architecture.run_NetlifeSpectrum(); // reliability engine 1;
	} else {
		// in this case, NoC does not have minimal adaptivity;
		noc_architecture.run_NetlifeSpectrum_no_tables_update(); // reliability engine 2;
	}


	// () clean up;
    //print_comm_statistics();
    clean();

	// () cpu time - includes also reliability estimation;
    end = times(&tmsend);
    printTimes(end-start, &tmsstart, &tmsend);

    return 1;
}

////////////////////////////////////////////////////////////////////////////////
//
// testbenches;
//
////////////////////////////////////////////////////////////////////////////////

int test_single_path()
{
	// debugging purposes only;

	NOC_ARCHITECTURE dummy_noc_architecture;
	dummy_noc_architecture.create_nodes_and_arcs_debug();
	dummy_noc_architecture.run_NetlifeSpectrum();

	return 1;
}
