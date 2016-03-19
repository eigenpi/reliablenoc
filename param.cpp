#include <cstdio>
#include "param.h"

Param::Param(void) {
    var_max_length = 64;
    precision = 4;

    do_mapping = false;
    dump_mapping = false;
    read_mapping = false;
    dump_routing = false;
    dump_apcg = false;

	// by default we do not do reliability optimization;
	do_reliability = false;
	alpha = 0.6;
	_q = 0.01;
	_M = 10000;
	// by default we assume normal NoC topology without any adaptivity;
	// this is to preserve the original tool behaviour;
	_tables_update = 0;
	
    // if true, then XY-routing is not necessarily used. Instead, the routing 
    // table will be synthesized to reflect the application requirement
    routing_table_synthesis = false;
    min_hit_threshold = 400; // 400;

    dump_mapping_table_file = NULL;
    read_mapping_table_file = NULL;
    routing_table_file = NULL;
    dump_apcg_file = NULL;
    verbose = false;

    link_bandwidth = 1000000;
    buf_read_ebit  = 1.056; // energy consumption per bit read
    buf_write_ebit = 2.831; // energy consumption per bit write
    link_ebit = 0.449; 
    switch_ebit = 0.284;
    legal_turn_set = odd_even;
    optimization_method = opt_bbm;
    routing_type = routing_xy;
    routing_effort = easy;

    pq_size = 2000; // 2000;
    link_usage_matrix = NULL;
    link_usage_list = NULL;

    seed = 5;

    parse_error = false;
    quit_flag = false;
}
