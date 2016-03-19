#ifndef APP_HPP_
#define APP_HPP_

#include <iostream>
#include <vector>
#include "arch.h"

//#define MAX_COMM_PER_TASK 30
#define DUMMY_VOL INT_MAX/100

using namespace std;

typedef enum Direction{Send, Recv} Direction;

struct Proc_Comm{ 
    int src_proc;
    int dst_proc;
    int BW;
    int adaptivity;        //only useful in routing synthesis

    // only useful energy aware routing
    int volume;
    float rate;
};
typedef struct Proc_Comm *pProc_Comm, Proc_Comm;

////////////////////////////////////////////////////////////////////////////////
//
// Process
//
////////////////////////////////////////////////////////////////////////////////

/********************************************************************************
 * Note that we assume each process will occupy an exclusive PE
 ********************************************************************************/
class Process {
    static int proc_num;              // the number of processes in the application
    Position pos;
    int tileId;
    int BW_scale;
    bool dummy;                // whether it's a dummy process
    class Process * master;    // if dummy, who is my master

public:
    int id;                    // the unique id of that process, from 0 -> app_size-1
	// toComm, fromComm are a bit of an overkill by storing potential
	// communication to any other process; this would be useful only in 
	// complete graphs;
    int *toComm; 
    int *fromComm;
    int *to_BW_requirement;    // the bandwidth requirement of the out-going traffic
    int *from_BW_requirement;  // the bandwidth requirement of the incoming traffic

    int rank;
    int totalCommVol;          // TODO: update this for ...

    // shape of the host PE
    int pe_width;
    int pe_length;

    // which tile it should be locked to. -1 if no lock
    int lock_to;

    vector<class Process *> spawned;

    void set_BW_scale(int i) {BW_scale = i;}
    int mapToTile() {return tileId;}
    int mapToTile(int tid) { tileId = tid; return tid; }
    Process();
    int GetId() const { return id; }
    ~Process();
    int ToComm(int i) { return toComm[i]; }     
    int FromComm(int i) { return fromComm[i]; }
    int get_pe_size(void) const { return pe_width * pe_length; }
    bool is_locked(void) const { return (lock_to != -1); }
    void spawn_dummy_process(void); 
    int add_from_comm_vol(int from_proc, int vol);
    int add_to_comm_vol(int to_proc, int vol);
    bool is_dummy(void) const { return dummy; }
    bool is_my_dummy_process(class Process * proc) const;
    bool check_myself(void) const;
    class Process * get_master(void) { return master; }
    friend ostream & operator << (ostream & os, const Process & proc);
};
typedef class Process *pProcess, Process;

ostream & operator<<(ostream & os, const Process & proc);

bool parse_apcg(char * fileName);
bool parse_traffic_config(char * fileName);

#endif
