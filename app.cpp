#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "app.h"
#include "gVal.h"
#include "nocmap.h"
#include "msg.h"

#define MAX_LINE 1024

Process::Process() {
    id = proc_num++;
    toComm = NULL;
    fromComm = NULL;
    to_BW_requirement = NULL;
    from_BW_requirement = NULL;
    toComm = new int[gTileNum];
    fromComm = new int[gTileNum];
    to_BW_requirement = new int[gTileNum];
    from_BW_requirement = new int[gTileNum];
    for (int i=0; i<gTileNum; i++) {
        toComm[i] = 0;
        fromComm[i] = 0;
        to_BW_requirement[i] = 0;
        from_BW_requirement[i] = 0;
    }

    // default shape of the hosting PE
    pe_width = 1;
    pe_length = 1;

    // by default it is not locked to a particular tile
    lock_to = -1;
    dummy = false;
    master = NULL;   // no master

    spawned.clear(); // no associated dummy pe
}

Process::~Process() {
    if (toComm)
        delete []toComm;
    if (fromComm)
        delete []fromComm;
    if (to_BW_requirement)
        delete []to_BW_requirement;
    if (from_BW_requirement)
        delete []from_BW_requirement; 
}

int Process::add_from_comm_vol(int from_proc, int vol) {
    if (vol == INT_MAX) 
        fromComm[from_proc] = INT_MAX;
    else
        fromComm[from_proc] += vol;
    return fromComm[from_proc];
}

int Process::add_to_comm_vol(int to_proc, int vol) {
    if (vol == INT_MAX)
        toComm[to_proc] = INT_MAX;
    else
        toComm[to_proc] += vol;
    return toComm[to_proc];
}

ostream & operator<<(ostream & os, const Process & proc) {
    os<<"Process "<<proc.id<<":"<<endl;
    os<<"Out going traffic: "<<endl;
    for (int i=0; i<proc.proc_num; i++) 
        os<<proc.toComm[i]<<"\t";
    os<<endl;
    return os;
} 

void Process::spawn_dummy_process(void) 
{
    if (get_pe_size() == 1) 
        return;
    int new_proc_cnt = get_pe_size() - 1;
    int edge_len = (int) sqrt((double) get_pe_size());
    assert(edge_len * edge_len == get_pe_size());
    spawned.clear();
    for (int i=0; i<new_proc_cnt; i++) {
        pProcess proc = new Process;
        proc->dummy = true;
        proc->master = this;
        proc->pe_length = proc->pe_width = 0;
        proc->BW_scale = BW_scale;
        spawned.push_back(proc);
        gProcess.push_back(proc);
    }

    /* fake pseudo comm volume of those nodes, 
     * the bandwidth requirement remains zero
     * currently, we only support square shape of 2x2 and 3x3. This is done manually.
     * for larger cases, we can later develop an algorithm to automatically do that
     */
    if (get_pe_size() == 4) {
        add_from_comm_vol(spawned[0]->GetId(), DUMMY_VOL);
        spawned[0]->add_to_comm_vol(GetId(), DUMMY_VOL);
        add_from_comm_vol(spawned[1]->GetId(), DUMMY_VOL);
        spawned[1]->add_to_comm_vol(GetId(), DUMMY_VOL);

        spawned[2]->add_from_comm_vol(spawned[0]->GetId(), DUMMY_VOL);
        spawned[0]->add_to_comm_vol(spawned[2]->GetId(), DUMMY_VOL);
        spawned[2]->add_from_comm_vol(spawned[1]->GetId(), DUMMY_VOL);
        spawned[1]->add_to_comm_vol(spawned[2]->GetId(), DUMMY_VOL);
    }

    if (get_pe_size() == 9) {
        add_from_comm_vol(spawned[0]->GetId(), DUMMY_VOL);
        spawned[0]->add_to_comm_vol(GetId(), DUMMY_VOL);
        add_from_comm_vol(spawned[1]->GetId(), DUMMY_VOL);
        spawned[1]->add_to_comm_vol(GetId(), DUMMY_VOL);
        add_from_comm_vol(spawned[2]->GetId(), DUMMY_VOL);
        spawned[2]->add_to_comm_vol(GetId(), DUMMY_VOL);
        add_from_comm_vol(spawned[3]->GetId(), DUMMY_VOL);
        spawned[3]->add_to_comm_vol(GetId(), DUMMY_VOL);

        spawned[4]->add_from_comm_vol(spawned[0]->GetId(), DUMMY_VOL);
        spawned[0]->add_to_comm_vol(spawned[4]->GetId(), DUMMY_VOL);
        spawned[4]->add_from_comm_vol(spawned[1]->GetId(), DUMMY_VOL);
        spawned[1]->add_to_comm_vol(spawned[4]->GetId(), DUMMY_VOL);

        spawned[5]->add_from_comm_vol(spawned[0]->GetId(), DUMMY_VOL);
        spawned[0]->add_to_comm_vol(spawned[5]->GetId(), DUMMY_VOL);
        spawned[5]->add_from_comm_vol(spawned[2]->GetId(), DUMMY_VOL);
        spawned[2]->add_to_comm_vol(spawned[5]->GetId(), DUMMY_VOL);

        spawned[6]->add_from_comm_vol(spawned[1]->GetId(), DUMMY_VOL);
        spawned[1]->add_to_comm_vol(spawned[6]->GetId(), DUMMY_VOL);
        spawned[6]->add_from_comm_vol(spawned[3]->GetId(), DUMMY_VOL);
        spawned[3]->add_to_comm_vol(spawned[6]->GetId(), DUMMY_VOL);

        spawned[7]->add_from_comm_vol(spawned[2]->GetId(), DUMMY_VOL);
        spawned[2]->add_to_comm_vol(spawned[7]->GetId(), DUMMY_VOL);
        spawned[7]->add_from_comm_vol(spawned[3]->GetId(), DUMMY_VOL);
        spawned[3]->add_to_comm_vol(spawned[7]->GetId(), DUMMY_VOL);
    }
}


// return true if @proc@ is my dummy process
bool Process::is_my_dummy_process(pProcess proc) const {
    if (spawned.size() == 0)
        return false;
    for (vector<pProcess>::const_iterator iter=spawned.begin(); iter<spawned.end(); iter++) {
        if (proc->GetId() == (*iter)->GetId())
            return true;
    }
    return false;
}


// check the process and return false if detects anything wong 
bool Process::check_myself(void) const {
    if (id < 0 || (id >= (int) gProcess.size()))
        return false;
    if (dummy) {
        for (unsigned int i=0; i<gProcess.size(); i++) {
            if (to_BW_requirement[i] != 0 || from_BW_requirement[i] != 0)
                return false;
            if (gProcess[i] != master) {
                // check to see if process[i] is my sibling
                bool sibling = master->is_my_dummy_process(gProcess[i]);
                if (sibling) {
                    // if my sibling, the comm volume should be either 0 or DUMMY_VOL
                    if ((toComm[i] != DUMMY_VOL && toComm[i] != 0) ||
                        (fromComm[i] != DUMMY_VOL && fromComm[i] != 0))
                        return false;
                }
                else {
                    if ((toComm[i] != 0 || fromComm[i] != 0))
                        return false;
                }
            }
            else {
                if ((toComm[i] != DUMMY_VOL && toComm[i] != 0) ||
                    (fromComm[i] != DUMMY_VOL && fromComm[i] != 0))
                    return false;
            }
            if (toComm[i] < 0 || fromComm[i] < 0 ||
                to_BW_requirement[i] < 0 || from_BW_requirement[i] < 0)
                return false;
            if (toComm[i] != gProcess[i]->fromComm[id] || fromComm[i] != gProcess[i]->toComm[id])
                return false;
            if (to_BW_requirement[i] != gProcess[i]->from_BW_requirement[id] || 
                from_BW_requirement[i] != gProcess[i]->to_BW_requirement[id])
                return false;
        }
        if (pe_width != 0 || pe_length != 0)
            return false;
    }
    else {
        if (pe_width < 1 || pe_length < 1) 
            return false;
        if (((int) spawned.size() + 1) != pe_width * pe_length)
            return false;
        // check the comm to and from my spawned pes
        for (unsigned int i=0; i<spawned.size(); i++) {
            int tmp = spawned[i]->GetId();
            if (to_BW_requirement[tmp] != 0 || from_BW_requirement[tmp] != 0)
                return false;
            if ((toComm[tmp] != DUMMY_VOL && toComm[tmp] != 0) || 
                (fromComm[tmp] != DUMMY_VOL && fromComm[tmp] != 0))
                return false;
        }
        for (unsigned int i=0; i<gProcess.size(); i++) {
            if (toComm[i] < 0 || fromComm[i] < 0 || to_BW_requirement[i] < 0 || 
                from_BW_requirement[i] < 0)
                return false;
            if (toComm[i] != gProcess[i]->fromComm[id] || 
                fromComm[i] != gProcess[i]->toComm[id])
                return false;
            if (to_BW_requirement[i] != gProcess[i]->from_BW_requirement[id] || 
                from_BW_requirement[i] != gProcess[i]->to_BW_requirement[id])
                return false;
        }
    }
    return true;
}

// APCG format:
// source    destination    communication-volume    bandwidth
// all in integer
bool parse_apcg(char * fileName) {
    char inputLine[MAX_LINE];
    ifstream inputFile(fileName);
    int src, dst, commVol, BW;
    if (!inputFile) {
        cerr<<"Error in openning file "<<fileName<<endl;
        return false;
    }

    while (inputFile.getline(inputLine, MAX_LINE, '\n')) {
        if (strchr(inputLine, '#'))
            continue;
        if (strlen(inputLine) == 0)
            continue;
        if (inputLine[0] == '\r' || inputLine[0] == '\n')
            continue;
        if (!sscanf(inputLine, "%d\t%d\t%d\t%d", &src, &dst, &commVol, &BW)) 
            return false;
        gProcess[src]->toComm[dst] = commVol;
        gProcess[src]->to_BW_requirement[dst] = BW;
        gProcess[dst]->fromComm[src] = commVol;
        gProcess[dst]->from_BW_requirement[src] = BW;
    }
    inputFile.close();
    return true;
}

// traffic_config file format:
// @NODE    i
// packet_to_destination_rate     k     "rate-in-float"
// We use a very mechanism to convert it to communicaiton volume
// and bandwidth metric needed by mapping:
bool parse_traffic_config(char * file_name) {
    if (param.verbose) 
        cout << "parsing traffic configuration file " << file_name << "..." << endl;

    char input_line[MAX_LINE];
    char * sub_string;

    ifstream in_file(file_name);
    if (!in_file.is_open() || !in_file.good())
        return false;

    int id = -1;

    while (in_file.getline(input_line, MAX_LINE, '\n')) {

        // line starting with "#" are coomments
        if (input_line[0] == '#')
            continue;

        if (strstr(input_line, "@NODE")) {
            sub_string = strstr(input_line, "@NODE");
            sub_string += strlen("@NODE");

            if (!sscanf(sub_string, "%d", &id))
                return -1;

            if (param.verbose) 
                cout << INFO_CONFIG_NODE << id << endl;

            continue;
        }

        if (strstr(input_line, "packet_to_destination_rate")) {
            sub_string = strstr(input_line, "packet_to_destination_rate");
            sub_string += strlen("packet_to_destination_rate");
            unsigned int dst_id;
            double rate = 0;
            if (!sscanf(sub_string, "%d\t%lf", &dst_id, &rate)) 
                return false;
            if (rate > 1.) {
                cerr << ERRO_PACKET_RATE;
                return false;
            }
            gProcess[id]->toComm[dst_id] = (int) (rate * 1000000);
            gProcess[id]->to_BW_requirement[dst_id] = (int) (rate * 3 * param.link_bandwidth);
            gProcess[dst_id]->fromComm[id] = (int) (rate * 1000000);
            gProcess[dst_id]->from_BW_requirement[id] = (int) (rate * 3 * param.link_bandwidth);
            continue;
        }

        // un-recognized options
        if (strlen(input_line) > 0) {
            cerr << WARN_OPTION_CONFIG_FILE << input_line << endl;
            continue;
        }
    }

    in_file.close();
    return true;
}

