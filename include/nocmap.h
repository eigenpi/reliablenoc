#ifndef _NOCMAP_H_
#define _NOCMAP_H_

#include <fstream>
#include <cassert>
#include "gVal.h"
#include "app.h"
#include "arch.h"


#define SEQUENTIAL_MAPPING 0
#define RANDOM_MAPPING 1


////////////////////////////////////////////////////////////////////////////////
//
// utilities; used throughout the source files...
//
////////////////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream & os, const Position & pos);
ostream & operator<<(ostream & os, const Position3D & pos);
int toAbsAddr(const Position & pos);
int toAbsAddr(const Position3D & pos);

bool parse_options(int argn, char ** argv);

bool initialize(void);
void clean();

//Free or build the link usage matrix
//The matrix is three dimensional. 
//If g_link_usage_matrix[tile1][tile2][linkId] = 1 means the 
//communication uses the corresponding link. Otherwise, the 
//link is not used for the communication.
void delete_link_usage_matrix();
void delete_link_usage_list();
void build_link_usage_matrix();

bool clearRoutingTables(void);

bool operator==(const Position & pos1, const Position &pos2);
bool operator!=(const Position & pos1, const Position &pos2);
bool operator==(const Position3D & pos1, const Position3D & pos2);
bool operator!=(const Position3D & pos1, const Position3D & pos2);

/***********************************************************
testRoutability is used to test whether packets from srcTile 
is routable to dstTile.
return TRUE: routable
       FALSE: Not routable.
************************************************************/
bool testRoutability(Tile& srcTile, Tile& dstTile);

/************************************************************
Function verify() is used to verify whether the chip's 
communication architecture can work correctly
************************************************************/
bool verify();

/************************************************************
This function is used to verify whether each link's 
communication load exceeds its bandwidth capability
************************************************************/
bool verify_BW_requirement();

/************************************************************
Function mapProcToFile maps processes to the tiles.
type: SEQUENTIAL_MAPPING
      RANDOM_MAPPING
return TRUE:  mapping succeed
       FALSE: mapping failed
************************************************************/
bool mapProcToTile(int type);

void padWriteFile(ofstream & oFile, int d1, int d2, int d3);
void padWriteFile(ofstream & oFile, int lId);

void print_help();

int peakComm();
int totalComm();

float calc_comm_energy();
float calc_switch_energy();
float calc_link_energy();
float calc_buffer_energy();
int calc_total_wirelength();
int calc_total_area();
double calc_total_reliability_cost();
void Anneal_Compute_normalization_factors();

void synthesizeRoutingTable();

void analyzeIt();

void dispMapping();

void printTimes(clock_t real, struct tms *tmsstart, struct tms *tmsend);

bool dump_mapping_table(char* file_name);
bool read_mapping_table(char* file_name);
bool dump_apcg(char * file_name);

int locate_link(Position & pos, int direction);

bool exist_locked_pe(void);
bool exist_non_regular_regions(void);

/**********************************************************************
This function generate array_size numbers between min and max, all the 
numbers are guaranteed to be different
**********************************************************************/ 
int * generate_random_array(int array_size, int min, int max);

void print_comm_statistics();

bool dump_routing_tables(char * file_name);


// Note: all the above should be made part of NOCMAP class; there should
// be no globals and no utility function called here and there thru the
// source files!

////////////////////////////////////////////////////////////////////////////////
//
// NOCMAP
//
////////////////////////////////////////////////////////////////////////////////

//class NOCMAP {
// private:
	//Param param;
	//vector<pProcess> gProcess;
	//pLink gLink = NULL;
	//pTile gTile = NULL;
	//int gProcNum;
	//int gLinkNum;
	//int gTileNum;
	//int g_edge_size; // How many tiles per edge; default 4;
	
// public:
//	NOCMAP() {
		//gProcess.clear();
		//gLink = NULL;
		//gTile = NULL;
		//gProcNum = 16;
		//gTileNum = 16;
		//gTileNum = 16;
		//g_edge_size = 4;
//	}
//	~NOCMAP() {};
//};

#endif
