//
//  main.cpp
//  tmc
//
//  Created by Penghang Liu on 10/27/21.
//  Copyright © 2021 Penghang Liu. All rights reserved.
//

#include <iostream>
#include "tmc.hpp"

int main(int argc, char * argv[]) {
    const auto t1 = chrono::steady_clock::now();
    
    string tmp (argv[1]);
    string gname = tmp.substr (tmp.find_last_of("/") + 1);
    int Max_event = stoi(argv[2]);
    int Max_memory = stoi(argv[3]);
    string consecutive (argv[4]);
    string out_file = "out_" + gname.substr(0,gname.size()-4) + "_" + argv[2] + "_" + argv[3] + "_" + argv[4];

    FILE* fp = fopen (out_file.c_str(), "w");
    
    cout << "Max event: " << Max_event << endl;
    cout << "Memory length: " << Max_memory << endl;
    
//Read file and create a sorted list of temporal events
    vector<event> events;
    createEvents(tmp, events);
        
    cout << "# of events:" << events.size() << endl;
    const auto t2 = chrono::steady_clock::now();
    print_time (fp, "Read data time: ", t2 - t1);
    
    instancemap instances;
    set<vector<event>> keys;
    for (size_t i=0; i<events.size(); i++) {
        if (i % 1000 == 0) {
            cout << i << " of " << events.size() << endl;
        }
        countInstance(events[i], instances, keys, Max_event, Max_memory, consecutive);
    }
    
    map<string, int> motif_count;
    for (auto it=instances.begin(); it!=instances.end(); ++it) {
        string motif = encodeMotif(it->first);
        vector<event> motif_events = it->first;
        if (motif_events.size() < 2) continue;
        motif_count[motif] += it->second.first;
    }
    
    const auto t3 = chrono::steady_clock::now();
    print_time (fp, "Count motifs time: ", t3 - t2);
    print_time (fp, "End-to-end Time: ", t3 - t1);
    
    for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
        fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
    }
    
    fclose (fp);
    cout << "ALL DONE" << endl;
    return 0;
}
