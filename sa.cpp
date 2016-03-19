#include "sa.h"
#include "mapping.h"

int last_x;
int attempts;
int zeroCostAcceptance;
bool needStop;

const float OVERLOAD_UNIT_COST = 1000000000.;

int * link_BW_usage = NULL;
int *** syn_link_BW_usage = NULL;
int **** sa_routing_table = NULL;

/*----------------------------------------------------*/
/* ways to gen Random Vars with specific distributions */
/*
 * Simple random number generator based on the linear-congruential method
 * using parameters from example D, p 40, Knuth Vol 2.
 *
 * uniform_rv() returns a real number uniformly distributed on [0,1]. This
 * version has the advantae that it should behave the same on different
 * machines, since the generator and starting point are explicitly
 * specified.
 */
double uniform_rv()
{
    /*
     * one small problem: the sequence we use can produce integers
     * larger than the word size used,  ie they can wrap around
     * negative.  We wimp out  on this matter and just make them
     * positive again.
     */
    last_x = ((A * last_x) + C) % M;
    if (last_x < 0)
        last_x = -last_x;
    return (((double)last_x) / ((double) M));
}

/* Initialize the random number stream */
void init_rand (int seed)
{
    last_x = seed;
}

/* produces a random INTEGER in [imin, imax],
**  ie, including the endpoints imin and imax
*/
long int uniform_int_rv(long int imin, long int imax)
{
    double  u;
    int     m;

    u = uniform_rv();
    m = (int) imin + ((int) floor((double)(imax + 1 - imin)*u));
    return m;
}

/* the usual metropolis accept criterion,
** returns 1 for accept, 0 reject
*/
int accept(double deltac, double temperature)
{
    double  pa;

    /* annealing accept criterion */
    if(deltac == 0) //accpet it, but record the number of zero cost accpetance
        zeroCostAcceptance++;

    if(deltac <= 0)
        return 1;
    else {
        pa = exp( (double)(-deltac)/temperature);
        if(uniform_rv() <= pa)
            return 1;
        else
            return 0;
    }
}

/* this does the actual evolution of the placement
** by annealing at a fixed tempt t passed in
** as an input.  startcost is also passed in
** so we can compute some statistics.
** acceptratio (fraction of moves tried that were accepted)
** is passed back since the caller wants to use it
** to help decide if we are frozen yet
*/
double anneal_at_temp(double t, double &currentcost, double *acceptratio)
{
    double newCost;
    int acceptcount = 0;
    double total_delta_cost = 0;

    static int zeroTempCnt = 0;
  
    int unit = attempts/10;

    //clear the zeroCostAcceptance
    zeroCostAcceptance = 0;

    /* this is the main loop doing moves.
    ** we do  'attempts' moves in all, then quit
    ** at this temperature
    */
    //cout<<"attempts = "<<attempts<<endl;
    for(int m=1; m<attempts; m++) {
        int tile1=0, tile2=0;
        makeRandomSwap(tile1, tile2);
		if ( param.do_reliability) { // NEW;
	        newCost = calc_total_cost_with_reliability();
		} else { // OLD;
	        newCost = calc_total_cost();
		}
        double delta_cost = newCost - currentcost;
        if(accept(delta_cost/currentcost*100, t)) {
            acceptcount++ ;
            total_delta_cost += delta_cost;
            currentcost = newCost;
        } /* if accept */
        else {  /* reject ! */
            rollBack(tile1, tile2);
        }
        if(m%unit==0) {
            //This is just to print out the process of the algorithm
            cerr<<"#";
            //cout<<"Current cost = "<<currentcost<<endl;
            //cout<<"Delta cost "<<delta_cost<<endl;
        }
    }
    cout<<endl;
    *acceptratio = ((double) acceptcount)/attempts;

    if(zeroCostAcceptance == acceptcount) 
        zeroTempCnt ++;
    else
        zeroTempCnt = 0;

    if(zeroTempCnt == 10)
        needStop = true;

    return total_delta_cost;
}

void anneal()
{
    double cost3, cost2, costcurrent;
    int   done;
    //double tol3, tol2, temp, new_temp;
    double tol3, tol2, temp;
    //double acceptratio, acceptratio_t, idea_acceptratio;
    double acceptratio;
    double delta_cost;

    if(!param.routing_table_synthesis)
        link_BW_usage = new int[gLinkNum];
    else{
        syn_link_BW_usage = new int**[g_edge_size];
        for(int i=0; i<g_edge_size; i++) {
            syn_link_BW_usage[i] = new int*[g_edge_size];
            for(int j=0; j<g_edge_size; j++) 
                syn_link_BW_usage[i][j] = new int[4];
        }
        sa_routing_table = new int***[g_edge_size];
        for(int i=0; i<g_edge_size; i++) {
            sa_routing_table[i] = new int**[g_edge_size];
            for(int j=0; j<g_edge_size; j++) {
                sa_routing_table[i][j] = new int*[gTileNum];
                for(int k=0; k<gTileNum; k++) {
                    sa_routing_table[i][j][k] = new int[gTileNum];
                    for(int m=0; m<gTileNum; m++) 
                        sa_routing_table[i][j][k][m] = -2;
                }
            }
        }
    }


	if ( param.do_reliability) {
		// I need this BB initialization only for the sake of computing the 
		// normalization factors; because of this I also have at the end
		// a call for BBMClear();
		Anneal_Compute_normalization_factors();
		
		//BBMInit();
	}

  
    /* set up the global control parameters for this
    ** annealing run
    */
    int tempcount = 0;
    cost3 = 999999999;
    cost2 = 999999999;
    costcurrent = cost2;

    attempts = gTileNum * gTileNum * 100;
    //attempts = gTileNum * 10;
  
    init_rand(1234567);
  
    //Determin initial temperature by accepting all moves and 
    //calculate variance.
    /* compute initial temperature 
       anneal_at_temp(10000.0, costcurrent, &acceptratio, 1);
       temp = 20.0 * VAR;
       init_anneal();*/

    temp = 100;
  
    /* here is the temperature cooling loop of the annealer */
    done = 0;
    do {
        needStop = false;
        cout<<"Round "<<tempcount<<":"<<endl;
        cout<<"Current Annealing temperature "<<temp<<endl;
		delta_cost = anneal_at_temp(temp, costcurrent, &acceptratio);
        cout<<"total delta cost "<<delta_cost<<endl;
        cout<<"Current cost "<<costcurrent<<endl;
        cout<<"Accept ratio "<<acceptratio<<endl;
    
        /* OK, if we got here the cost function is working
        ** fine. We can now look at whether we are
        ** frozen, or whether we should cool some more.
        ** We basically just look at the last 2
        ** temperatures, and see if the cost is not
        ** changing much (thats the TOLERANCE test)
        ** and if the we have done enough temperatures
        ** (thats the TEMPS test), and if the accept
        ** ratio fraction is small enough (that is the
        ** MINACCEPT test).  If all are satisfied,
        ** we quit.
        */
        tol3 = ((double)cost3 - (double)cost2)/ (double)cost3;
        if(tol3 < 0) tol3 = -tol3;
        tol2 = ((double)cost2 - (double)costcurrent)/ (double)cost2;
        if(tol2 < 0) tol2 = -tol2;

        if( tol3 < TOLERANCE
            && tol2 < TOLERANCE
            && tempcount > TEMPS
            && (acceptratio < MINACCEPT || needStop)) {
            done = 1;
        }
        else {  /* no, not frozen */
            /* save the relevant info to test for frozen
            ** after the NEXT temperature.
            */
            cost3 = cost2;
            cost2 = costcurrent;
		
            /* cooling schedule 
               acceptratio_t = get_ideal_ratio();
               new_temp = (1.0 - (acceptratio - acceptratio_t)/K)*temp;
               temp = new_temp;*/
		
            /*
              new_temp = temp * exp((-0.7)*temp/VAR);
              if (new_temp/temp<0.995)
              temp = 0.995 * temp;
              else
              temp = new_temp;
            *********/
            temp = 0.9*temp;
            tempcount ++;
        }
    }while(!done);   /* go back and do annealing at next cooler temp */
    if (param.routing_table_synthesis)
        sa_program_routers();
    if (link_BW_usage)
        delete []link_BW_usage;
    if (syn_link_BW_usage) {
        for (int i=0; i<g_edge_size; i++) {
            for (int j=0; j<g_edge_size; j++) 
                delete []syn_link_BW_usage[i][j];
            delete []syn_link_BW_usage[i];
        }
        delete []syn_link_BW_usage;
    }
    if (sa_routing_table) {
        for(int i=0; i<g_edge_size; i++) {
            for(int j=0; j<g_edge_size; j++) {
                for(int k=0; k<gTileNum; k++) 
                    delete []sa_routing_table[i][j][k];
                delete []sa_routing_table[i][j];
            }
            delete []sa_routing_table[i];
        }
        delete []sa_routing_table;
    }


	if ( param.do_reliability) {
		// Note: see the comments above, when I used BBMInit();
		//BBMClear();
	}
}

void makeRandomSwap(int & tile1, int & tile2) {
    tile1 = uniform_int_rv(0, gTileNum-1);

#ifdef DEBUG
    if(tile1>=gTileNum||tile1<0) {
        cerr<<"Error in random swap..."<<endl;
        exit(0);
    }
#endif //DEBUG
    while(1) {
        //select two tiles to swap
        tile2 = uniform_int_rv(0, gTileNum-1);
        if(tile1!=tile2 && (gTile[tile1].mapToProc()!=-1||gTile[tile2].mapToProc()!=-1))
            break;
    }

    //Swap the processes attached to these two tiles
    swapProcesses(tile1, tile2);
}

void swapProcesses(int tile1, int tile2) {
    int p1 = gTile[tile1].mapToProc();
    int p2 = gTile[tile2].mapToProc();
    gTile[tile1].mapToProc(p2);
    gTile[tile2].mapToProc(p1);
    if(p1!=-1) 
        gProcess[p1]->mapToTile(tile2);
    if(p2!=-1)
        gProcess[p2]->mapToTile(tile1);
}

void rollBack(int tile1, int tile2) {
    swapProcesses(tile1, tile2);
}

float calc_total_cost()
{
	// Calculate the total cost in terms of the sum of the energy consumption and
	// the penalty of the link overloading
    float energy_cost = calc_comm_energy(); // the communication energy part
    float overload_cost;
    // now calculate the overloaded BW cost
    if(!param.routing_table_synthesis) 
        overload_cost = fixed_routing_overload_calc();
    else
        overload_cost = adaptive_routing_overload_calc();
    return ( energy_cost + overload_cost);
}

float calc_total_cost_with_reliability()
{
	// this is derived from calc_total_cost(); we consider also the 
	// reliability cost too here;
    float energy_cost = calc_comm_energy(); // the communication energy part
    float overload_cost = 0.0;
    float cost_e = 0.0;
    float cost_r = 0.0;
    // now calculate the overloaded BW cost
    if ( !param.routing_table_synthesis) {
        overload_cost = fixed_routing_overload_calc();
    } else {
        overload_cost = adaptive_routing_overload_calc();
	}
	cost_e = energy_cost + overload_cost;
	// reliability cost also;
	cost_r = calc_total_reliability_cost();
	// normalization;
	cost_e = cost_e / COST_ENERGY_NORMALIZATION;
	cost_r = cost_r / COST_RELIABILITY_NORMALIZATION;
	//printf(" (%.2f %.2f)",cost_e, cost_r);
	float cost = ( (1 - param.alpha) * cost_e + param.alpha * cost_r);
    return ( cost);
}

float fixed_routing_overload_calc() {
    memset((void *) link_BW_usage, 0, sizeof(int)*gLinkNum);
    for (int proc1=0; proc1<gProcNum; proc1++) {
        for (int proc2=proc1+1; proc2<gProcNum; proc2++) {
            if (gProcess[proc1]->to_BW_requirement[proc2] > 0) {
                int tile1 = gProcess[proc1]->mapToTile();
                int tile2 = gProcess[proc2]->mapToTile();
                for (unsigned int i=0; i<param.link_usage_list[tile1][tile2].size(); i++) {
                    int link_id = param.link_usage_list[tile1][tile2][i];
                    link_BW_usage[link_id] += gProcess[proc1]->to_BW_requirement[proc2];
                }
            }
            if(gProcess[proc1]->from_BW_requirement[proc2] > 0) {
                int tile1 = gProcess[proc1]->mapToTile();
                int tile2 = gProcess[proc2]->mapToTile();
                for (unsigned int i=0; i<param.link_usage_list[tile2][tile1].size(); i++) {
                    int link_id = param.link_usage_list[tile2][tile1][i];
                    link_BW_usage[link_id] += gProcess[proc1]->from_BW_requirement[proc2];
                }
            }
        }
    }
    float overload_cost = 0;
    for (int i=0; i<gLinkNum; i++) {
        if (link_BW_usage[i] > gLink[i].bandwidth) 
            overload_cost = ((float) link_BW_usage[i])/gLink[i].bandwidth - 1.;
    }
    overload_cost *= OVERLOAD_UNIT_COST;
    return overload_cost;
}

float adaptive_routing_overload_calc() {
    float overload_cost = 0.0;

    //Clear the link usage
    for(int i=0; i<g_edge_size; i++) {
        for(int j=0; j<g_edge_size; j++) 
            memset((void *) syn_link_BW_usage[i][j], 0, sizeof(int)*4);
    }

    for(int src=0; src<gProcNum; src++) {
        for(int dst=0; dst<gProcNum; dst++) {
            int tile1 = gProcess[src]->mapToTile();
            int tile2 = gProcess[dst]->mapToTile();
            if(gProcess[src]->to_BW_requirement[dst]) 
                sa_route_traffic(tile1, tile2, gProcess[src]->to_BW_requirement[dst]);
        }
    }

    for(int i=0; i<g_edge_size; i++) {
        for(int j=0; j<g_edge_size; j++) {
            for(int k=0; k<4; k++) {
                if(syn_link_BW_usage[i][j][k]> gLink[0].bandwidth)
                    overload_cost += ((float) syn_link_BW_usage[i][j][k])/gLink[0].bandwidth - 1.;
            }
        }
    }
    overload_cost *= OVERLOAD_UNIT_COST;
    return overload_cost;
}

bool sa_route_traffic(int src_tile, int dst_tile, int BW) {
    bool commit = true;
    Position src = gTile[src_tile].GetPosition();
    Position dst = gTile[dst_tile].GetPosition();
    int row = src.row;
    int col = src.col;

    int ***BW_usage = syn_link_BW_usage;
    int direction = -2;
    while(row!=dst.row || col!=dst.col) {
        //For west-first routing
        if(param.legal_turn_set == west_first) {
            if(col > dst.col) //step west
                direction = WEST;
            else if(col == dst.col) 
                direction = (row<dst.row)? NORTH:SOUTH;
            else if(row == dst.row) 
                direction = EAST;
            //Here comes the flexibility. We can choose whether to go vertical or horizontal
            else { 
                int direction1 = (row<dst.row)?NORTH:SOUTH;
                if(BW_usage[row][col][direction1] < BW_usage[row][col][EAST])
                    direction = direction1;
                else if(BW_usage[row][col][direction1] > BW_usage[row][col][EAST])
                    direction = EAST;
                else { 
                    //In this case, we select the direction which has the longest 
                    //distance to the destination
                    if((dst.col-col)*(dst.col-col) <= (dst.row-row)*(dst.row-row)) 
                        direction = direction1;
                    else //Horizontal move
                        direction = EAST;
                }
            }
        }
        //For odd-even routing
        else if(param.legal_turn_set == odd_even) {
            int e0 = dst.col - col;
            int e1 = dst.row - row;
            if(e0==0) //currently the same column as destination
                direction = (e1>0)?NORTH:SOUTH;
            else {
                if(e0>0) {        //eastbound messages
                    if(e1==0)
                        direction = EAST;
                    else {
                        int direction1=-1, direction2=-1;
                        if(col%2==1 || col==src.col) 
                            direction1 = (e1>0)?NORTH:SOUTH;
                        if(dst.col%2==1 || e0!=1) 
                            direction2 = EAST;
                        if(direction1==-1&&direction2==-1) {
                            cerr<<"Error, debug me."<<endl;
                            exit(1);
                        }
                        if(direction1==-1)
                            direction = direction2;
                        else if(direction2==-1) 
                            direction = direction1;
                        else //we have two choices
                            direction = 
                                (BW_usage[row][col][direction1]<BW_usage[row][col][direction2])?direction1:direction2;
                    }
                }
                else {   //westbound messages
                    if(col%2!=0 || e1==0)
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

        if(commit)
            sa_routing_table[row][col][src_tile][dst_tile] = direction;
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

bool sa_program_routers() {
    //clean all the old routing table
    for(int tile_id=0; tile_id<gTileNum; tile_id++) {
        for(int src_tile=0; src_tile<gTileNum; src_tile++) {
            for(int dst_tile=0; dst_tile<gTileNum; dst_tile++) {
                if(tile_id == dst_tile) 
                    gTile[tile_id].set_routing_entry(src_tile, dst_tile, -1);
                else
                    gTile[tile_id].set_routing_entry(src_tile, dst_tile, -2);
            }
        }
    }

    for(int row=0; row<g_edge_size; row++) {
        for(int col=0; col<g_edge_size; col++) {
            int tile_id = row*g_edge_size + col;
            for(int src_tile=0; src_tile<gTileNum; src_tile++) {
                for(int dst_tile=0; dst_tile<gTileNum; dst_tile++) {
                    Position pos;
                    pos.row = row; pos.col = col;
                    int link_id = locate_link(pos, sa_routing_table[row][col][src_tile][dst_tile]);
                    if(link_id!=-1)
                        gTile[tile_id].set_routing_entry(src_tile, dst_tile, link_id);
                }
            }
        }
    }
    return true;
}

