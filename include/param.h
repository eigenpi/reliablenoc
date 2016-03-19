#ifndef PARAM_HPP_
#define PARAM_HPP_

#include <vector>


using namespace std;


#define MAX_LINE 1024

typedef enum Optimization_method {opt_annealing, opt_bbm, opt_method_invalid} Optimization_method;
typedef enum Effort_type {easy, hard} Effort_type;
typedef enum ROUTING_TYPE {routing_xy, routing_oddeven} ROUTING_TYPE;
typedef enum Adaptive_routing_scheme {west_first, odd_even} Adaptive_routing_scheme;


////////////////////////////////////////////////////////////////////////////////
//
// PARAM
//
////////////////////////////////////////////////////////////////////////////////

typedef struct Param {

    // some program constants
    int var_max_length;
    int precision;

    // network system parameters
    int link_bandwidth;
    float buf_read_ebit;
    float buf_write_ebit;
    float link_ebit; 
    float switch_ebit;
    ROUTING_TYPE routing_type;
    Adaptive_routing_scheme legal_turn_set;
    Optimization_method optimization_method;
    Effort_type routing_effort;

    // mapping tool configruation
    int pq_size;   // priority queue size

    // mapping tool internal data structure
    int ***link_usage_matrix;
    vector<int> ** link_usage_list;

    int seed;

    // program control
    double alpha; 
	bool do_reliability; // if true, we include reliability in cost calculations;
	double _q; // F0(t) probability link is "down" (of link failure); default 0.05;
	int _M; // Monte Carlo runs; default is 10000;
	int _tables_update; // "1" for yes; "0" for no;
    bool routing_table_synthesis;
    bool do_mapping;
    bool dump_mapping;
    bool read_mapping;
    bool dump_routing;
    bool dump_apcg;
    bool verbose;
    int min_hit_threshold;

    // file names
    char exec_str[MAX_LINE];
    char* dump_mapping_table_file;
    char* read_mapping_table_file;
    char* routing_table_file;
    char* dump_apcg_file;

    // option parsing related 
    bool parse_error;
    bool quit_flag;

    Param(void);
};
typedef Param *pParam;

#endif
