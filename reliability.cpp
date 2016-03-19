#include <cassert>
#include <math.h>
#include "reliability.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// NOC_ARCHITECTURE
//
////////////////////////////////////////////////////////////////////////////////

void NOC_ARCHITECTURE::create_nodes_and_arcs()
{
	// this must be called only once;

	//      10    11
	//   6 --- 7 --- 8
	//  5|    7|     |9
	//   3 --- 4 --- 5
	//  0|  6 2|  8  |4  ...
	//   0 --- 1 --- 2   <--- node id's
	//      1     3      <--- arc id's
	//

	// (1) basics;
	_nx = g_edge_size;
	_ny = g_edge_size;
	_nodes_count = gTileNum;

	_q = param._q;
	_p = 1 - _q;
	_M_runs = param._M;	

	// (2) nodes;
	int x = 0, y = 0;
	int north_id = -1, east_id = -1, south_id = -1, west_id = -1;
	for ( int r_id = 0; r_id < _nodes_count; r_id++) {
		x = r_id % _nx; // col;
		y = r_id / _nx; // row;
		_nodes.push_back( ROUTERS_NODE( r_id, x, y));
		// set also adjacent neighbors at north, south, east, and west;
		north_id = -1, east_id = -1, south_id = -1, west_id = -1;
		if ( x > 0) west_id = r_id - 1;
		if ( x < _nx-1) east_id = r_id + 1;
		if ( y > 0) south_id = r_id - _nx;
		if ( y < _ny-1) north_id = r_id + _nx;
		// add neighbors in clockwise order;
		_nodes[ r_id].add_adj_nodes( north_id); // 0
		_nodes[ r_id].add_adj_nodes( east_id);  // 1
		_nodes[ r_id].add_adj_nodes( south_id); // 2
		_nodes[ r_id].add_adj_nodes( west_id);  // 3
		// add also dummy adj_arcs here; they will later be populated
		// with the correct info when we'll create the actual links;
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
	}

	// (3) set their assigned core id's; this is same as inside dispMapping();
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
		_nodes[ tid].set_assigned_core( proc_id);
    }
		
	// (4) arcs; go thru each node and take their North and East neighbors;
	int id1 = 0, id2 = 0;
	for ( int i = 0; i < _nodes_count; i++) {
		id1 = i;
		ROUTERS_NODE *this_node = get_node( i);
		id2 = this_node->adj_nodes()[ NORTH_DIR];
		if ( id2 >= 0) {
			//assert( id2 < _nodes_count);
			_arcs.push_back( ROUTERS_ARC( _arcs_count, id1, id2));
			_nodes[ id1].set_adj_arcs( NORTH_DIR, _arcs_count);
			_nodes[ id2].set_adj_arcs( SOUTH_DIR, _arcs_count);
			_arcs_count ++; // prepare for next one;
		}
		id2 = this_node->adj_nodes()[ EAST_DIR];		
		if ( id2 >= 0) {
			//assert( id2 < _nodes_count);
			_arcs.push_back( ROUTERS_ARC( _arcs_count, id1, id2));
			_nodes[ id1].set_adj_arcs( EAST_DIR, _arcs_count);
			_nodes[ id2].set_adj_arcs( WEST_DIR, _arcs_count);
			_arcs_count ++; // prepare for next one;
		}
	}

	// (5) create communication pairs; similar to BBMSortProcesses();
	double comm_vol = 0;
	for (unsigned int i=0; i<gProcess.size(); i++) {
        gProcess[i]->totalCommVol = 0;
        for (unsigned int k=0; k<gProcess.size(); k++) {
            comm_vol = gProcess[i]->toComm[k];
			if ( comm_vol > 0) {
				// source tile;
				int tid_src = gProcess[i]->mapToTile();
				int row = tid_src / g_edge_size;
				int col = tid_src % g_edge_size;
				tid_src = col * g_edge_size + row;
				// destination tile;
				int tid_dst = gProcess[k]->mapToTile();
				row = tid_dst / g_edge_size;
				col = tid_dst % g_edge_size;
				tid_dst = col * g_edge_size + row;
				_nodes[ tid_src].add_comm_pairs( tid_dst);
				_comm_pairs_count ++;
			}
        }
    }

	// (6) allocate mem and initialize sketch arrays for reliability engine;
	_permutation_pi.resize( _arcs_count);
	_Fr.resize( _arcs_count);
	for ( int i = 0; i < _arcs_count; i++) {
		_permutation_pi[ i] = -1;
		_Fr[ i] = 0.0;	
	}
	_all_N.resize( _comm_pairs_count);
	_all_fhat_r.resize( _comm_pairs_count);
	_all_RN.resize( _comm_pairs_count);
	for ( int i = 0; i < _comm_pairs_count; i++) {
		_all_N[i].resize( _arcs_count);
		_all_fhat_r[i].resize( _arcs_count);
		for ( int j = 0; j < _arcs_count; j++) {
			_all_N[i][j] = 0;
			_all_fhat_r[i][j] = 0.0;
		}
		_all_RN[i] = 0.0;
	}


	// (7) sanity checks;
	assert ( _arcs_count == _nx*(_ny-1) + _ny*(_nx-1));
	//print_noc_architecture(); exit(1); // debug;
}


void NOC_ARCHITECTURE::create_nodes_and_arcs( NOC_ARCHITECTURE *big_noc,
	int src, int dst, bool walk_NE, bool walk_NW)
{
	// Note: only walk_NE or walk_NW can be true;
	// Note: this function is used only for graph_bbox objects, which represent areas
	// of the big noc;

	// () coordinates within big noc;
	int x1 = -1, x2 = -1, y1 = -1, y2 = -1;
	int x_big = -1, y_big = -1;
	x1 = big_noc->nodes()[src].x();
	y1 = big_noc->nodes()[src].y();	
	x2 = big_noc->nodes()[dst].x();
	y2 = big_noc->nodes()[dst].y();

	// local ones;
	_nx = abs( x2 - x1) + 1;
	_ny = abs( y2 - y1) + 1;
	_nodes_count = _nx * _ny;

	// () nodes; node 0 here corresponds to node "s" from the big noc
	// and last node here corresponds to node "d" from big noc;
	int x = 0, y = 0;
	int north_id = -1, east_id = -1, south_id = -1, west_id = -1;
	for ( int r_id = 0; r_id < _nodes_count; r_id++) {
		x = r_id % _nx; // col within this graph (local);
		y = r_id / _nx; // row;

		x_big = (walk_NE == true) ? x1 + x : x1 - x;
		y_big = y1 + y;
		int node_id_in_big_noc = big_noc->node_id( x_big, y_big);
		// only for the scope of this type of objects (graph_bbox) I use 
		// "assigned_core" storage to store the id of node in the big noc 
		// (to keep it "simple" and save mem);
		_nodes.push_back( ROUTERS_NODE( r_id, x, y, node_id_in_big_noc));

		// set also adjacent neighbors at north, south, east, and west;
		north_id = -1, east_id = -1, south_id = -1, west_id = -1;
		if ( x > 0) west_id = r_id - 1;
		if ( x < _nx-1) east_id = r_id + 1;
		if ( y > 0) south_id = r_id - _nx;
		if ( y < _ny-1) north_id = r_id + _nx;
		// add neighbors in clockwise order;
		_nodes[ r_id].add_adj_nodes( north_id); // 0
		_nodes[ r_id].add_adj_nodes( east_id);  // 1
		_nodes[ r_id].add_adj_nodes( south_id); // 2
		_nodes[ r_id].add_adj_nodes( west_id);  // 3
		// add also dummy adj_arcs here; they will later be populated
		// with the correct info when we'll create the actual links;
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
	}
	
	// () arcs; go thru each node and take their North and East neighbors;
	int id1 = 0, id2 = 0;
	int arc_id_in_big_noc = -1;
	DIRECTION dir_in_big_noc = NORTH_DIR;
	int this_weight = 0;
	for ( int i = 0; i < _nodes_count; i++) {
		id1 = i;

		// x_big,y_big will be used to retreive arc id's in the big host noc;
		x = id1 % _nx; // col local;
		y = id1 / _nx; // row local;
		x_big = (walk_NE == true) ? x1 + x : x1 - x;
		y_big = y1 + y;

		ROUTERS_NODE *this_node = get_node( i);
		id2 = this_node->adj_nodes()[ NORTH_DIR];		
		if ( id2 >= 0) {
			dir_in_big_noc = NORTH_DIR;
			arc_id_in_big_noc = big_noc->arc_id( this_weight,x_big,y_big, dir_in_big_noc);
			_arcs.push_back( ROUTERS_ARC( _arcs_count,id1,id2, arc_id_in_big_noc,this_weight));
			_nodes[ id1].set_adj_arcs( NORTH_DIR, _arcs_count);
			_nodes[ id2].set_adj_arcs( SOUTH_DIR, _arcs_count);
			_arcs_count ++; // prepare for next one;
		}
		id2 = this_node->adj_nodes()[ EAST_DIR];		
		if ( id2 >= 0) {
			dir_in_big_noc = (walk_NE == true) ? EAST_DIR : WEST_DIR;
			arc_id_in_big_noc = big_noc->arc_id( this_weight,x_big, y_big, dir_in_big_noc);
			_arcs.push_back( ROUTERS_ARC( _arcs_count,id1,id2, arc_id_in_big_noc,this_weight));
			_nodes[ id1].set_adj_arcs( EAST_DIR, _arcs_count);
			_nodes[ id2].set_adj_arcs( WEST_DIR, _arcs_count);
			_arcs_count ++; // prepare for next one;
		}
	}

	// () sanity checks;
	assert ( _arcs_count == _nx*(_ny-1) + _ny*(_nx-1));
}

////////////////////////////////////////////////////////////////////////////////
//
// some functions only for creating testbenches - debugging purposes only;
//
////////////////////////////////////////////////////////////////////////////////

void NOC_ARCHITECTURE::create_nodes_and_arcs_debug()
{
	// Note: only used for debug purposes;

	_nx = 3;
	_ny = 5;
	_nodes_count = _nx * _ny;
	_q = param._q;
	_p = 1 - _q;
	_M_runs = param._M;	
	printf("\nR if s,t are on a straight line: %f\n", pow(_p, (_nx-1)+(_ny-1)) );
	
	// nodes;
	int x = 0, y = 0;
	int north_id = -1, east_id = -1, south_id = -1, west_id = -1;
	for ( int r_id = 0; r_id < _nodes_count; r_id++) {
		x = r_id % _nx; // col;
		y = r_id / _nx; // row;
		_nodes.push_back( ROUTERS_NODE( r_id, x, y));
		// set also adjacent neighbors at north, south, east, and west;
		north_id = -1, east_id = -1, south_id = -1, west_id = -1;
		if ( x > 0) west_id = r_id - 1;
		if ( x < _nx-1) east_id = r_id + 1;
		if ( y > 0) south_id = r_id - _nx;
		if ( y < _ny-1) north_id = r_id + _nx;
		// add neighbors in clockwise order;
		_nodes[ r_id].add_adj_nodes( north_id); // 0
		_nodes[ r_id].add_adj_nodes( east_id);  // 1
		_nodes[ r_id].add_adj_nodes( south_id); // 2
		_nodes[ r_id].add_adj_nodes( west_id);  // 3
		// add also dummy adj_arcs here; they will later be populated
		// with the correct info when we'll create the actual links;
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);
		_nodes[ r_id].add_adj_arcs( -1);

		// set also assigned core id as router id;
		_nodes[ r_id].set_assigned_core( r_id);
	}

	// arcs;
	int id1 = 0, id2 = 0;
	for ( int i = 0; i < _nodes_count; i++) {
		id1 = i;
		ROUTERS_NODE *this_node = get_node( i);
		id2 = this_node->adj_nodes()[ NORTH_DIR];
		if ( id2 >= 0) {
			//assert( id2 < _nodes_count);
			_arcs.push_back( ROUTERS_ARC( _arcs_count, id1, id2));
			_nodes[ id1].set_adj_arcs( NORTH_DIR, _arcs_count);
			_nodes[ id2].set_adj_arcs( SOUTH_DIR, _arcs_count);
			_arcs_count ++; // prepare for next one;
		}
		id2 = this_node->adj_nodes()[ EAST_DIR];		
		if ( id2 >= 0) {
			//assert( id2 < _nodes_count);
			_arcs.push_back( ROUTERS_ARC( _arcs_count, id1, id2));
			_nodes[ id1].set_adj_arcs( EAST_DIR, _arcs_count);
			_nodes[ id2].set_adj_arcs( WEST_DIR, _arcs_count);
			_arcs_count ++; // prepare for next one;
		}
	}

	// create one communication pair only;
	int tid_src = 0;
	int tid_dst = _nodes_count - 1;
	_nodes[ tid_src].add_comm_pairs( tid_dst);
	_comm_pairs_count ++;

	// allocate mem and initialize sketch arrays for reliability engine;
	_permutation_pi.resize( _arcs_count);
	_Fr.resize( _arcs_count);
	for ( int i = 0; i < _arcs_count; i++) {
		_permutation_pi[ i] = -1;
		_Fr[ i] = 0.0;	
	}
	_all_N.resize( _comm_pairs_count);
	_all_fhat_r.resize( _comm_pairs_count);
	_all_RN.resize( _comm_pairs_count);
	for ( int i = 0; i < _comm_pairs_count; i++) {
		_all_N[i].resize( _arcs_count);
		_all_fhat_r[i].resize( _arcs_count);
		for ( int j = 0; j < _arcs_count; j++) {
			_all_N[i][j] = 0;
			_all_fhat_r[i][j] = 0.0;
		}
		_all_RN[i] = 0.0;
	}

}

////////////////////////////////////////////////////////////////////////////////
//
// NOC_ARCHITECTURE
//
////////////////////////////////////////////////////////////////////////////////

void NOC_ARCHITECTURE::print_noc_architecture()
{
	printf("\n\nRouters graph");
	printf("\nCommunication pairs: %d", _comm_pairs_count);
	printf("\nNodes: %d", _nodes.size());
	for ( int i = 0; i < _nodes_count; i ++) {
		printf("\n Node id: %d (%d, %d) Assigned_core: %d",
			i, _nodes[i].x(), _nodes[i].y(), _nodes[i].assigned_core());
		printf("\n\tadj nodes:");
		int adj_count = _nodes[i].adj_nodes().size();
		for ( int k = 0; k < adj_count; k ++) {
			printf(" %d", _nodes[i].adj_nodes()[k]);
		}
		printf("\n\tadj arcs:");
		int adj_arcs_count = _nodes[i].adj_arcs().size();
		for ( int k = 0; k < adj_arcs_count; k ++) {
			printf(" %d", _nodes[i].adj_arcs()[k]);
		}
		printf("\n\tcomm pairs:");
		int comm_count = _nodes[i].comm_pairs().size();
		for ( int k = 0; k < comm_count; k ++) {
			printf(" %d", _nodes[i].comm_pairs()[k]);
		}
	}
	printf("\nArcs: %d", _arcs.size());
	for ( int j = 0; j < _arcs_count; j++) {
		printf("\nArc id:%d    %d --- %d    m:%d    \t%d", _arcs[j].id(),
			   _arcs[j].src_id(), _arcs[j].des_id(), 
				_arcs[j].marked(), _arcs[j].weight());
	}
	printf("\n");
}
void NOC_ARCHITECTURE::print_permutation_pi()
{
	printf("\npermutation_pi: \n");
	for ( int i = 0; i < _arcs_count; i++) {
		printf("%d \t", i);
	}
	printf("\n");
	for ( int i = 0; i < _arcs_count; i++) {
		printf("%d \t", _permutation_pi[ i]);
	}
}
void NOC_ARCHITECTURE::print_all_fhat_r()
{
	printf("\nSpectrums: \n");
	for ( int i = 0; i < _comm_pairs_count; i++) {
		for ( int j = 0; j < _arcs_count; j++) {
			printf("%.2f ", _all_fhat_r[i][j]);
		}
		printf("\n");
	}
}
void NOC_ARCHITECTURE::print_Fr()
{
	printf("\nOrder statistics CDF Fr's: \n");
	for ( int i = 0; i < _arcs_count; i ++) {
		cout << _Fr[i] << "  ";
	}
	cout << endl;
}
void NOC_ARCHITECTURE::simulate_random_permutation()
{
	// (1) clear stuff;
	for ( int i = 0; i < _arcs_count; i++) {
		_permutation_pi[ i] = -1;
	}
	// (2) generate random permutation;
	for ( int i = 0; i < _arcs_count; i++) {
		int k = rand() % _arcs_count;
		while ( _permutation_pi[ k] != -1) {
			k = rand() % _arcs_count;
		}
		_permutation_pi[ k] = i;
	}
}

////////////////////////////////////////////////////////////////////////////////
//
// network spectrum related;
//
////////////////////////////////////////////////////////////////////////////////

void NOC_ARCHITECTURE::compute_Frs()
{
	// see paper for the formula;
	int m = _arcs_count;
	double m_factorial = factorial( m);
	double coef = 0.0;
	double prod = 0.0;
	for ( int r = 1; r <= _arcs_count; r++) {
		double this_Fr = 0.0;
		for ( int j = r; j <= _arcs_count; j++) {
			coef = m_factorial / ( factorial(j) * factorial(m-j) );
			prod = pow( _q, j) * pow( _p, m-j);
			this_Fr += coef * prod;
		}
		_Fr[ r-1] = this_Fr;
	}
}

////////////////////////////////////////////////////////////////////////////////
//
// NetlifeSpectrum engine 1
//
////////////////////////////////////////////////////////////////////////////////

bool NOC_ARCHITECTURE::run_NetlifeSpectrum()
{
	// reliability engine 1;
	// the main algorithm for estimating reliability of a given mapping solution;
	// Note: this is under the assumption that the NoC architecture has
	// support for minimal adaptivity;

    cout << "Estimate reliability..." << endl;


	// (1) compute first all Fr's (CDF's of order statistics tau_r; see paper);
	compute_Frs();


	// (2) the main algorithm;
	int this_min_edge = -1;
	int src = 0, dst = 0;
	for ( int run_i = 0; run_i < _M_runs; run_i++) {
		// () generate random permutation;
		simulate_random_permutation();
		//print_permutation_pi(); // debug;

		// () assign arcs' weights;
		for ( int j = 0; j < _arcs_count; j++) {
			_arcs[ j].set_weight( _permutation_pi[j]);
		}

		// () find the maximal-weight "monotonic XY-YX" spanning tree 
		// for each communication pair; then record the minimal-weight
		// edge of each of these trees;
		int pair_counter = 0;
		for ( int i = 0; i < _nodes_count; i++) {
			src = i;
			int comm_count = _nodes[i].comm_pairs().size();
			for ( int k = 0; k < comm_count; k++) {
				dst = _nodes[i].comm_pairs()[k];
				// find the max-weight spanning tree and get its min-weight edge;
				// version 1: simpler, but under-estimates reliability:
				//this_min_edge = min_edge_of_max_monotonic_path( src, dst);
				// version 2: more accurate;
				this_min_edge = min_edge_of_max_spanning_tree( src, dst);

				// all individually computed;
				_all_N[pair_counter][this_min_edge] = _all_N[pair_counter][this_min_edge] + 1;

				pair_counter++;
			}
		}
	}


	// (3) compute fhat_r's;
	for ( int i = 0; i < _comm_pairs_count; i++) {
		for ( int j = 0; j < _arcs_count; j++) {
			_all_fhat_r[i][j] = double(_all_N[i][j]) / double(_M_runs);
		}
	}
	//print_all_fhat_r(); // debug;
	

	// (4) estimate _FN, _RN;
	double max_FN = 0.0;
	double this_comm_FN = 0.0;
	int weakest_id = -1;
	for ( int i = 0; i < _comm_pairs_count; i++) {
		this_comm_FN = 0.0;
		for ( int j = 0; j < _arcs_count; j++) {
			this_comm_FN += _all_fhat_r[i][j] * _Fr[j];
		}
		if ( this_comm_FN > max_FN) {
			max_FN = this_comm_FN;
			weakest_id = i;
		}
	}
	_FN = max_FN;
	_RN = 1 - _FN;
	//cout << "FN(t) = " << max_FN << endl;
	cout << "Reliability = " << 1.0 - max_FN << endl; // RN(t);

	return true;
}

int NOC_ARCHITECTURE::min_edge_of_max_monotonic_path( int src, int dst)
{
	// Note: currently not used;
	// Note: this is a simpler, but, under-estimates reliability?
	int min_edge = INT_MAX; // result;

	int s = src, t = dst;
	int x1 = -1, x2 = -1, y1 = -1, y2 = -1;
	x1 = _nodes[src].x();
	y1 = _nodes[src].y();	
	x2 = _nodes[dst].x();
	y2 = _nodes[dst].y();


	// (1) horizontal path;
	if ( y1 == y2) {
		assert( x1 != x2); // we do not have comm pairs with same src and dst;
		if ( x1 > x2) {
			int temp = x1; x1 = x2; x2 = temp;
			s = dst; t = src;
		}
		// this is now something like: s---o---o---o---t
		int arc_id = -1;
		while ( s != t) { // walk east-wards;
			arc_id = _nodes[s].adj_arcs()[EAST_DIR];
			if ( _arcs[ arc_id ].weight() < min_edge) {
				min_edge = _arcs[ arc_id ].weight();
			}
			s = s + 1;
		}
		return min_edge;
	}


	// (2) vertical path;
	else if ( x1 == x2) {
		assert( y1 != y2);
		if ( y1 > y2) {
			int temp = y1; y1 = y2; y2 = temp;
			s = dst; t = src;
		}
		int arc_id = -1;
		while ( s != t) { // walk north-wards;
			arc_id = _nodes[s].adj_arcs()[NORTH_DIR];
			if ( _arcs[ arc_id ].weight() < min_edge) {
				min_edge = _arcs[ arc_id ].weight();
			}
			s = s + _ny;
		}
		return min_edge;
	}


	// (3) general case;
	else {
		// set up things so that we walk either NE or NW only;
		bool walk_NW = false, walk_NE = false;
		if ( x1 < x2 && y1 > y2) {
			walk_NW = true;
			int temp = x1; x1 = x2; x2 = temp;
			temp = y1; y1 = y2; y2 = temp;
			s = dst; t = src;
		} else if ( x1 < x2 && y1 < y2) {
			walk_NE = true;
		} else if ( x1 > x2 && y1 > y2) {
			walk_NE = true;
			int temp = x1; x1 = x2; x2 = temp;
			temp = y1; y1 = y2; y2 = temp;
			s = dst; t = src;
		} else if ( x1 > x2 && y1 < y2) {
			walk_NW = true;
		}
		// clear parents and path_maxweigth; instead of cleaning all nodes, we
		// clean only those between "s" and "t";
		for ( int i = 0; i < _nodes_count; i++) {
			_nodes[i].clear_parent();
			_nodes[i].clear_parent_arc();
			_nodes[i].clear_path_maxweight();
		}


		// (a) walk_NE;
		int this_id = -1, next_id = -1, prev_id = -1;
		int arc_id = -1;
		int this_path_maxweight = -1;
		if ( walk_NE == true) {

			deque<long> ne_queue; // fifo;
			ne_queue.push_back( s);
			while ( !ne_queue.empty()) {
				// get node from queue and process its north and east neighbors, if any;
				this_id = ne_queue.front();
				ne_queue.pop_front(); // remove it also;
				// (a1) visit EAST_DIR neighbor if any;
				if ( _nodes[this_id].x() < x2) {
					next_id = _nodes[this_id].adj_nodes()[EAST_DIR];
					arc_id = _nodes[this_id].adj_arcs()[EAST_DIR];
					// set the neighbors parent as this_id, compute its path-length too
					// and add it to the queue;
					if ( _nodes[next_id].parent_is_set() == false) {
						// this is the first time when we arrived to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						_nodes[next_id].set_path_maxweight( this_path_maxweight);
						_nodes[next_id].set_parent( this_id);
						_nodes[next_id].set_parent_arc( arc_id);
					} else {
						// need to check if this path has bigger weight; only in that case
						// we overwrite the previously found max path to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						if ( this_path_maxweight > _nodes[next_id].path_maxweight()) {
							_nodes[next_id].set_path_maxweight( this_path_maxweight);
							_nodes[next_id].set_parent( this_id);
							_nodes[next_id].set_parent_arc( arc_id);
						}
					}
					ne_queue.push_back( next_id);
				}
				// (a2) visit NORTH_DIR neighbor if any;
				if ( _nodes[this_id].y() < y2) {
					next_id = _nodes[this_id].adj_nodes()[NORTH_DIR];
					arc_id = _nodes[this_id].adj_arcs()[NORTH_DIR];
					// set the neighbors parent as this_id, compute its path-length too
					// and add it to the queue;
					if ( _nodes[next_id].parent_is_set() == false) {
						// this is the first time when we arrived to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						_nodes[next_id].set_path_maxweight( this_path_maxweight);
						_nodes[next_id].set_parent( this_id);
						_nodes[next_id].set_parent_arc( arc_id);
					} else {
						// need to check if this path has bigger weight; only in that case
						// we overwrite the previously found max path to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						if ( this_path_maxweight > _nodes[next_id].path_maxweight()) {
							_nodes[next_id].set_path_maxweight( this_path_maxweight);
							_nodes[next_id].set_parent( this_id);
							_nodes[next_id].set_parent_arc( arc_id);
						}
					}
					ne_queue.push_back( next_id);
				}
			} //while ( !ne_queue.empty())

			// (a3) at this time all nodes within the bbox have set their parents and "t"
			// stores the maximum weight path to it from "s"; now, we backtrace and
			// record what is the minimum weight edge along this maximal-weight path;
			int parent_arc_weight = 0;
			prev_id = t;
			while ( prev_id != s) {
				this_id = prev_id;
				parent_arc_weight = _arcs[ _nodes[this_id].parent_arc() ].weight();
				if ( parent_arc_weight < min_edge) {
					min_edge = parent_arc_weight;
				}
				prev_id = _nodes[this_id].parent();
			}
		} //if ( walk_NE == true)
		
		
		// (b) walk_NW;
		if ( walk_NW == true) {

			deque<long> nw_queue; // fifo;
			nw_queue.push_back( s);
			while ( !nw_queue.empty()) {
				// get node from queue and process its north and west neighbors, if any;
				this_id = nw_queue.front();
				nw_queue.pop_front(); // remove it also;
				// (b1) visit WEST_DIR neighbor if any;
				if ( _nodes[this_id].x() > x2) {
					next_id = _nodes[this_id].adj_nodes()[WEST_DIR];
					arc_id = _nodes[this_id].adj_arcs()[WEST_DIR];
					// set the neighbors parent as this_id, compute its path-length too
					// and add it to the queue;
					if ( _nodes[next_id].parent_is_set() == false) {
						// this is the first time when we arrived to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						_nodes[next_id].set_path_maxweight( this_path_maxweight);
						_nodes[next_id].set_parent( this_id);
						_nodes[next_id].set_parent_arc( arc_id);
					} else {
						// need to check if this path has bigger weight; only in that case
						// we overwrite the previously found max path to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						if ( this_path_maxweight > _nodes[next_id].path_maxweight()) {
							_nodes[next_id].set_path_maxweight( this_path_maxweight);
							_nodes[next_id].set_parent( this_id);
							_nodes[next_id].set_parent_arc( arc_id);
						}
					}
					nw_queue.push_back( next_id);
				}
				// (b2) visit NORTH_DIR neighbor if any;
				if ( _nodes[this_id].y() < y2) {
					next_id = _nodes[this_id].adj_nodes()[NORTH_DIR];
					arc_id = _nodes[this_id].adj_arcs()[NORTH_DIR];
					// set the neighbors parent as this_id, compute its path-length too
					// and add it to the queue;
					if ( _nodes[next_id].parent_is_set() == false) {
						// this is the first time when we arrived to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						_nodes[next_id].set_path_maxweight( this_path_maxweight);
						_nodes[next_id].set_parent( this_id);
						_nodes[next_id].set_parent_arc( arc_id);
					} else {
						// need to check if this path has bigger weight; only in that case
						// we overwrite the previously found max path to the neighbor node;
						this_path_maxweight =
							_nodes[this_id].path_maxweight() + _arcs[arc_id].weight();
						if ( this_path_maxweight > _nodes[next_id].path_maxweight()) {
							_nodes[next_id].set_path_maxweight( this_path_maxweight);
							_nodes[next_id].set_parent( this_id);
							_nodes[next_id].set_parent_arc( arc_id);
						}
					}
					nw_queue.push_back( next_id);
				}
			} //while ( !nw_queue.empty())

			// (b3) at this time all nodes within the bbox have set their parents and "t"
			// stores the maximum weight path to it from "s"; now, we backtrace and
			// record what is the minimum weight edge along this maximal-weight path;
			int parent_arc_weight = 0;
			prev_id = t;
			while ( prev_id != s) {
				this_id = prev_id;
				parent_arc_weight = _arcs[ _nodes[this_id].parent_arc() ].weight();
				if ( parent_arc_weight < min_edge) {
					min_edge = parent_arc_weight;
				}
				prev_id = _nodes[this_id].parent();
			}
		} //if ( walk_NW == true)


	} //(3)

	return min_edge;
}

int NOC_ARCHITECTURE::min_edge_of_max_spanning_tree( int src, int dst)
{
	// Note: this is the function which is currently used;

	int min_edge = INT_MAX; // result;
	int s = src, t = dst;
	int x1 = -1, x2 = -1, y1 = -1, y2 = -1;
	x1 = _nodes[src].x();
	y1 = _nodes[src].y();	
	x2 = _nodes[dst].x();
	y2 = _nodes[dst].y();

	// (1) horizontal path;
	if ( y1 == y2) {
		assert( x1 != x2); // we do not have comm pairs with same src and dst;
		if ( x1 > x2) {
			int temp = x1; x1 = x2; x2 = temp;
			s = dst; t = src;
		}
		int arc_id = -1;
		while ( s != t) { // walk east-wards;
			arc_id = _nodes[s].adj_arcs()[EAST_DIR];
			if ( _arcs[ arc_id ].weight() < min_edge) {
				min_edge = _arcs[ arc_id ].weight();
			}
			s = s + 1;
		}
		return min_edge;
	}

	// (2) vertical path;
	else if ( x1 == x2) {
		assert( y1 != y2);
		if ( y1 > y2) {
			int temp = y1; y1 = y2; y2 = temp;
			s = dst; t = src;
		}
		int arc_id = -1;
		while ( s != t) { // walk north-wards;
			arc_id = _nodes[s].adj_arcs()[NORTH_DIR];
			if ( _arcs[ arc_id ].weight() < min_edge) {
				min_edge = _arcs[ arc_id ].weight();
			}
			s = s + _ny;
		}
		return min_edge;
	}

	// (3) general case;
	else {
		// set up things so that we walk either NE or NW only;
		bool walk_NW = false, walk_NE = false;
		if ( x1 < x2 && y1 > y2) {
			walk_NW = true;
			int temp = x1; x1 = x2; x2 = temp;
			temp = y1; y1 = y2; y2 = temp;
			s = dst; t = src;
		} else if ( x1 < x2 && y1 < y2) {
			walk_NE = true;
		} else if ( x1 > x2 && y1 > y2) {
			walk_NE = true;
			int temp = x1; x1 = x2; x2 = temp;
			temp = y1; y1 = y2; y2 = temp;
			s = dst; t = src;
		} else if ( x1 > x2 && y1 < y2) {
			walk_NW = true;
		}
		// clear parents and path_maxweigth; instead of cleaning all nodes, we
		// clean only those between "s" and "t";
		for ( int i = 0; i < _nodes_count; i++) {
			_nodes[i].clear_parent();
			_nodes[i].clear_parent_arc();
			_nodes[i].clear_path_maxweight();
		}

		// TO DO: maybe introduce a new class for this purpose?
		// build noc graph corresponding to routers within bbox only; this is 
		// just a light weight NOC_ARCHITECTURE object;
		NOC_ARCHITECTURE graph_bbox;

		// (a) walk_NE;
		if ( walk_NE == true) {
			// walk thru the routers within the bounding box and create
			// temporary graph (inside which we'll do the kruskal action);
			graph_bbox.create_nodes_and_arcs( this, s, t, true, false); // walk_NE,walk_NW;
		}

		// (b) walk_NW;
		if ( walk_NW == true) {
			graph_bbox.create_nodes_and_arcs( this, s, t, false, true); // walk_NE,walk_NW;
		}

		// (c) compute maximum weight spanning tree, then prune (see paper and book),
		// then get minimum weight edge;
		min_edge = graph_bbox.compute_spanning_tree( src,dst);

	} //(3)

	return min_edge;
}

int NOC_ARCHITECTURE::compute_spanning_tree( int src, int dst)
{
	// TO DO: speed this up; currently is written in a hurry;
	// compute maximum spanning tree;
	int min_edge = INT_MAX; // result;

	// ()
	vector<PAIR_TWO> sorted_arcs;
	for ( int j = 0; j < _arcs_count; j++) {
		sorted_arcs.push_back( PAIR_TWO( j, _arcs[j].weight()));
	}
	sort( sorted_arcs.begin(), sorted_arcs.end());

	int nodes_count_1 = _nodes_count - 1; // number of nodes minus 1;

	// () maximum weight spanning tree is constructed and recorded by
	// "marking" arcs and nodes (nodes are marked by seting _parent); 
	int id1 = -1, id2 = -1; // each arc has source node id1 and dest node id2;
	int arc_id = -1;
	// first add the arc with the largest weight;
	arc_id = sorted_arcs[_arcs_count - 1].id();
	_arcs[arc_id].set_marked( 1); // add arc to spanning tree;
	_nodes[ _arcs[arc_id].src_id() ].set_parent( 1);
	_nodes[ _arcs[arc_id].des_id() ].set_parent( 1);
	int discovered_arcs = 1;
	while ( discovered_arcs < nodes_count_1) {	
		for ( int j = _arcs_count - 2; j >= 0; j--) {
			arc_id = sorted_arcs[j].id(); // next arc with largest weight;
			if ( _arcs[ arc_id].marked() == -1) { // not a marked arc;
				id1 = _arcs[arc_id].src_id();
				id2 = _arcs[arc_id].des_id();		
				if ( (_nodes[id1].parent() != 1 && _nodes[id2].parent() != 1) ||
					 (_nodes[id1].parent() == 1 && _nodes[id2].parent() == 1)) {
					continue;
				} else {
					_arcs[arc_id].set_marked( 1); // add arc to spanning tree;
					_nodes[ id1 ].set_parent( 1);
					_nodes[ id2 ].set_parent( 1);
					discovered_arcs++;
					break;
				}
			}
		}
	}

	// () prune "hanging" arcs so that only "s"=0 and "t"=_nodes;
	bool pruned_at_least_one = true;

	while ( pruned_at_least_one) {
		pruned_at_least_one = false;
		for ( int j = 0; j < _arcs_count; j++) {
			if ( _arcs[j].marked() == 1) {
				// check if any of its terminals is hanging; if so prune this arc;
				id1 = _arcs[j].src_id();
				id2 = _arcs[j].des_id();
				// check if any node is a "dead-end";
				bool id1_dead_end = false;
				if ( num_marked_arcs( id1) <= 1) {
					id1_dead_end = true;
				}
				bool id2_dead_end = false;
				if ( num_marked_arcs( id2) <= 1) {
					id2_dead_end = true;
				}
				// to prune or not to prune? 
				if ( id1_dead_end == true && id1 != 0 && id1 != nodes_count_1) {
					_nodes[ id1 ].set_parent( -1);
					_arcs[j].set_marked( -1);
					pruned_at_least_one = true;
				}
				if ( id2_dead_end == true && id2 != 0 && id2 != nodes_count_1) {
					_nodes[ id2 ].set_parent( -1);
					_arcs[j].set_marked( -1);
					pruned_at_least_one = true;
				}
			}
		}
	}	

	// () find the minimum weight edge among all that remained;
	for ( int j = 0; j < _arcs_count; j++) {
		if ( _arcs[j].marked() == 1) {
			if ( _arcs[ j ].weight() < min_edge) {
				min_edge = _arcs[ j ].weight();
			}
		}
	}

	return min_edge;
}

////////////////////////////////////////////////////////////////////////////////
//
// NetlifeSpectrum engine 2
//
////////////////////////////////////////////////////////////////////////////////

bool NOC_ARCHITECTURE::run_NetlifeSpectrum_no_tables_update()
{
	// reliability engine 2 (no routing tables update);
	// Note: this is under the assumption that the NoC architecture has no
	// support for minimal adaptivity (routing tables are not updated
	// when link failures occur); that is the initial routing path
	// is the only one used; it is as if "each path is a straight path" - 
	// its spanning tree is always the only path itself;

    cout << "Estimate reliability (no minimal adaptivity)..." << endl;

	// (1) compute first all Fr's (CDF's of order statistics tau_r; see paper);
	compute_Frs();

	// (2) the main algorithm;
	int this_min_edge = -1;
	int src = 0, dst = 0;
	for ( int run_i = 0; run_i < _M_runs; run_i++) {
		// () generate random permutation;
		simulate_random_permutation();

		// () assign arcs' weights;
		for ( int j = 0; j < _arcs_count; j++) {
			_arcs[ j].set_weight( _permutation_pi[j]);
		}

		// () find the maximal-weight "monotonic XY-YX" spanning tree 
		// for each communication pair; then record the minimal-weight
		// edge of each of these trees;
		int pair_counter = 0;
		for ( int i = 0; i < _nodes_count; i++) {
			src = i;
			int comm_count = _nodes[i].comm_pairs().size();
			for ( int k = 0; k < comm_count; k++) {
				dst = _nodes[i].comm_pairs()[k];
				// find the min-weight edge along the unique routing path;
				this_min_edge = min_edge_of_this_unique_routing_path( src, dst);
				
				// all individually computed;
				_all_N[pair_counter][this_min_edge] = _all_N[pair_counter][this_min_edge] + 1;

				pair_counter ++;
			}
		}
	}


	// (3) compute fhat_r's;
	for ( int i = 0; i < _comm_pairs_count; i++) {
		for ( int j = 0; j < _arcs_count; j++) {
			_all_fhat_r[i][j] = double(_all_N[i][j]) / double(_M_runs);
		}
	}
	//print_all_fhat_r(); // debug;
	

	// (4) estimate _FN, _RN;
	double max_FN = 0.0;
	double this_comm_FN = 0.0;
	int weakest_id = -1;
	for ( int i = 0; i < _comm_pairs_count; i++) {
		this_comm_FN = 0.0;
		for ( int j = 0; j < _arcs_count; j++) {
			this_comm_FN += _all_fhat_r[i][j] * _Fr[j];
		}
		if ( this_comm_FN > max_FN) {
			max_FN = this_comm_FN;
			weakest_id = i;
		}
	}
	_FN = max_FN;
	_RN = 1 - _FN;
	//cout << "FN(t) = " << max_FN << endl;
	cout << "Reliability = " << 1.0 - max_FN << endl; // RN(t);

	return true;
}

int NOC_ARCHITECTURE::min_edge_of_this_unique_routing_path( int src, int dst)
{
	int min_edge = INT_MAX; // result;

	// Note: here I have to work with the globals again :( that is because 
	// I do not record the routing paths inside NOC_ARCHITECTURE object 
	// (maybe we should do it to avoid globals);
	// Note: another thing is that inside the nocmap main code tiles are
	// counted 0,1,2,... bottom-up from left-to-right; dumping of routing
	// however and my indexing too is counting 0,1,2,... left-right from
	// bottom-to-up; so, to use the globals I have to take care of this
	// conversion all the time :(

	// src and dst are in the 
	int s = src;
	int row = s / g_edge_size;
	int col = s % g_edge_size;
	s = col * g_edge_size + row;
	int t = dst;
	row = t / g_edge_size;
	col = t % g_edge_size;
	t = col * g_edge_size + row;

	int prev_tid_in_new_coord = src;
	int this_tid_in_new_coord = src;
	pTile currentTile = &gTile[s];
	while ( currentTile->GetId() != t) {
		int linkId = currentTile->RouteToLink(s, t);
		pLink pL = &gLink[linkId];

		int arc_id = -1;
		// get tile id in the new coordinate system;
		this_tid_in_new_coord = pL->ToTile();
		row = this_tid_in_new_coord / g_edge_size;
		col = this_tid_in_new_coord % g_edge_size;
		this_tid_in_new_coord = col * g_edge_size + row;
		// see where we go from prev to this and find the arc_id;
		if ( this_tid_in_new_coord == prev_tid_in_new_coord + 1) {
			// we walked from "prev" to east;
			arc_id = _nodes[ prev_tid_in_new_coord].adj_arcs()[1];
		} else if ( this_tid_in_new_coord == prev_tid_in_new_coord - 1) {
			// we walked from "prev" to West;
			arc_id = _nodes[ prev_tid_in_new_coord].adj_arcs()[3];
		} else if ( this_tid_in_new_coord > prev_tid_in_new_coord) {
			// we walked from "prev" to North;
			arc_id = _nodes[ prev_tid_in_new_coord].adj_arcs()[0];
		} else if ( this_tid_in_new_coord < prev_tid_in_new_coord) {
			// we walked from "prev" to South;
			arc_id = _nodes[ prev_tid_in_new_coord].adj_arcs()[2];
		}

		// is this the edge along the path with min weight?
		if ( _arcs[ arc_id ].weight() < min_edge) {
			min_edge = _arcs[ arc_id ].weight();
		}

		currentTile = &gTile[pL->ToTile()]; // continue to "walk";
		prev_tid_in_new_coord = this_tid_in_new_coord;
	}

	return min_edge;
}
