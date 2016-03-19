#ifndef _RELIABILITY_H_
#define _RELIABILITY_H_

#include <deque>
#include "gVal.h"


// I re-define my own directions here to make easy the implementation
// of the reliability engines;
enum DIRECTION { NORTH_DIR = 0, EAST_DIR = 1, SOUTH_DIR = 2, WEST_DIR = 3 };


////////////////////////////////////////////////////////////////////////////////
//
// PAIR_TWO
//
////////////////////////////////////////////////////////////////////////////////

class PAIR_TWO
{
 private:
	int _id;
	int _value;
 public:
	PAIR_TWO() { _id = 0; _value = 0; }
	PAIR_TWO(int id, int value) : _id(id), _value(value) {}
	PAIR_TWO(const PAIR_TWO &pt) : 
		_id(pt._id), _value(pt._value) {}
	~PAIR_TWO() {}

	int id() const { return _id; }
	int value() const { return _value; }
	void set_id(int id) { _id = id; }
	void set_value(int value) { _value = value; }
	void operator=(const PAIR_TWO &pt) { _id = pt._id; _value = pt._value; }
	bool operator==(const PAIR_TWO &pt) const { return (_value == pt._value); }
	bool operator!=(const PAIR_TWO &pt) const { return (_value != pt._value); }
	bool operator<(const PAIR_TWO &pt)  const { return (_value  < pt._value); }
};

////////////////////////////////////////////////////////////////////////////////
//
// ROUTERS_NODE
//
////////////////////////////////////////////////////////////////////////////////

class ROUTERS_NODE {
 private:
	int _id; // id of itself;
	int _x, _y; // coordinates;
	vector<int> _adj_nodes; // id's of adjacent nodes;
	vector<int> _adj_arcs; // id's of links which connect to this router;
	// _comm_pairs store the id of the routers (or tiles) with wich this router
	// form a communication pair in the original application graph; this router
	// is the source and the destinations are recorded here; note that a 
	// communication pair is stored once only at the router which is source;
	vector<int> _comm_pairs;
	// id of core mapped to the tile of this router; if no core is mapped to it,
	// then this stays -1;
	int _assigned_core;
	// variables utilized in the maximal-weight path calculations;
	int _parent;
	int _parent_arc;
	int _path_maxweight;
		
 public:
	ROUTERS_NODE( int id, int x, int y) : _id(id), _x(x), _y(y) {
		_assigned_core = -1;
		_parent = -1; _parent_arc = -1; _path_maxweight = 0;
	}
	ROUTERS_NODE( int id, int x, int y, int assigned_core) : // used for bbox graphs only;
		_id(id), _x(x), _y(y), _assigned_core(assigned_core) {
		_parent = -1; _parent_arc = -1; _path_maxweight = 0;
	}
	~ROUTERS_NODE() {}
	
	int id() const { return _id; }
	int x() const { return _x; }
	int y() const { return _y; }
	int assigned_core() const { return _assigned_core; }
	void set_assigned_core( int id) { _assigned_core = id; }
	vector<int> &adj_nodes() { return _adj_nodes; }
	void add_adj_nodes( int id) { // id on node;
		_adj_nodes.push_back( id);
	}
	vector<int> &adj_arcs() { return _adj_arcs; }
	void add_adj_arcs( int id) { // id of arc;
		_adj_arcs.push_back( id);
	}
	void set_adj_arcs(int dir, int id) { // id of arc;
		_adj_arcs[dir] = id;
	}
	vector<int> &comm_pairs() { return _comm_pairs; }
	void add_comm_pairs( int id) {
		_comm_pairs.push_back( id);
	}
	// maximal-weight spanning tree;
	int parent() const { return _parent; }
	void set_parent( int p) { _parent = p; }
	int parent_arc() const { return _parent_arc; }
	void set_parent_arc( int p_arc) { _parent_arc = p_arc; }
	int path_maxweight() const { return _path_maxweight; }
	void set_path_maxweight( int w) { _path_maxweight = w; }
	void clear_parent() { _parent = -1; }
	void clear_parent_arc() { _parent_arc = -1; }
	void clear_path_maxweight() { _path_maxweight = 0; }
	bool parent_is_set() { return ( _parent >= 0 ? true : false ); }
};

////////////////////////////////////////////////////////////////////////////////
//
// ROUTERS_ARC
//
////////////////////////////////////////////////////////////////////////////////

class ROUTERS_ARC {
 private:
	int _id;
	int _src_id;
	int _des_id;
	int _weight;
	// _id_in_big_noc is the counterpart id of the arc in the big noc host;
	// used for graph_bbox purposes only;
	int _id_in_big_noc;
	int _marked;
		
 public:
	ROUTERS_ARC( int id, int src, int des) : _id(id), _src_id(src), _des_id(des) {
		_weight = -1;
		_id_in_big_noc = -1;
		_marked = -1;
	}
	ROUTERS_ARC( int id, int src, int des, int id_in_big_noc, int weight) :
		_id(id), _src_id(src), _des_id(des), _weight(weight),
		_id_in_big_noc(id_in_big_noc) {
		_marked = -1;
	}
	~ROUTERS_ARC() {}

	int id() const { return _id; }
	int id_in_big_noc() const { return _id_in_big_noc; }
	int src_id() const { return _src_id; }
	int des_id() const { return _des_id; }
	int weight() const { return _weight; }
	void set_weight( int weight) { _weight = weight; }
	int marked() const { return _marked; }
	void set_marked( int marked) { _marked = marked; }
};

////////////////////////////////////////////////////////////////////////////////
//
// NOC_ARCHITECTURE
//
////////////////////////////////////////////////////////////////////////////////

class NOC_ARCHITECTURE {
 private:
	// things related to the actual regular architecture;
	int _nx, _ny; // number of routers in each direction;
	int _nodes_count;
	int _arcs_count;
	vector<ROUTERS_NODE> _nodes;
	vector<ROUTERS_ARC> _arcs;
	int _comm_pairs_count;

	// things related to reliability estimation - engine 1;
	// this is Algorithm - NetlifeSpectrum NoC (see paper);
	int _M_runs; // Monte Carlo runs; default is 10000;
	double _q;  // F0(t) probability link is "down" (of link failure);
	double _p;  // 1-F0(t) probability link is "up";
	vector<int> _permutation_pi;
	vector<double> _Fr;
	// we look at all communication pairs as separate s-t communications
	// for which we compute an individual reliability; then we take the minimum 
	// reliability among all pairs as the weakest link; the computations share 
	// the NoC platform and so correlations are considered?
	vector<vector<int> > _all_N;
	vector<vector<double> > _all_fhat_r; // spectrum for each communication pair;
	vector<double> _all_RN;
	// final result;
	double _FN; // FN(t);
	double _RN; // RN(t) = 1 - FN(t);

 public:
	NOC_ARCHITECTURE() {
		_nodes_count = 0; _arcs_count = 0; _comm_pairs_count = 0;
		_M_runs = 10000; // default 10000
		_FN = 0.0; _RN = 0.0;
		_q = 0.01; // default 0.01;
		_p = 1 - _q;
	}
	~NOC_ARCHITECTURE() {}
	
	int nodes_count() { return _nodes_count; }
	int arcs_count() { return _arcs_count; }
	vector<ROUTERS_NODE> &nodes() { return _nodes; }
	ROUTERS_NODE *get_node( int id) {
		//assert(id >= 0 && id < _nodes_count);
		return &_nodes[ id]; 
	}
	ROUTERS_ARC *get_arc( int id) {
		return &_arcs[ id]; 
	}
	int node_id( int x_big, int y_big) {
		// location gives the id based on how we construct the router nodes;
		return ( y_big * _nx + x_big);
	}
	int arc_id( int &weight, int x_big, int y_big, DIRECTION dir) {
		int node_id = y_big * _nx + x_big;
		int arc_id = _nodes[node_id].adj_arcs()[dir];
		weight = _arcs[arc_id].weight();
		return ( arc_id);
	}
	void create_nodes_and_arcs();
	void create_nodes_and_arcs( NOC_ARCHITECTURE *big_noc,
		int src, int dst, bool walk_NE, bool walk_NW);
	void simulate_random_permutation();
	bool run_NetlifeSpectrum();
	bool run_NetlifeSpectrum_no_tables_update();
	void compute_Frs();
	int min_edge_of_max_monotonic_path( int src, int dst);
	int min_edge_of_max_spanning_tree( int src, int dst);
	int min_edge_of_this_unique_routing_path( int src, int dst);
	
	// debug;
	void create_nodes_and_arcs_debug();
	void print_noc_architecture();
	void print_permutation_pi();
	void print_all_fhat_r();
	void print_Fr();

	// utils;
	int dist_manhattan( int x1, int y1, int x2, int y2) {
		return ( abs(x1 - x2) + abs(y1 - y2));
	}
	double factorial( int c) {
		double res = 1;
		for ( int j = 2; j <= c; j++) {
			res *= j;
		} return res;
	}
	// Kruskal spanning tree;
	int compute_spanning_tree( int src, int dst);
	int num_marked_arcs( int node_id) {
		int count = 0;
		int arc_id = -1;
    	for ( int i=0; i<4; i++) {
			arc_id = _nodes[node_id].adj_arcs()[i];
			if ( arc_id >= 0) { // valid neighbor?
				if ( _arcs[ arc_id].marked() == 1) count++;
			}
		}
		return count;
	}
};

#endif
