//
//  main.cpp
//  mtm
//
//  Created by Penghang Liu on 10/27/21.
//  Copyright © 2021 Penghang Liu. All rights reserved.
//

#include <iostream>
#include "mtm.hpp"

int main(int argc, char * argv[]) {
    const auto t1 = chrono::steady_clock::now();
    
    string tmp (argv[1]);
    string gname = tmp.substr (tmp.find_last_of("/") + 1);
    int Max_event = stoi(argv[2]);
    int Max_memory = stoi(argv[3]);
    int K = stoi(argv[4]);
    string out_file = gname.substr(0,gname.find_last_of(".")) + "_" + argv[2] + "_" + argv[3] + "_" + argv[4] + ".txt";

    FILE* fp = fopen (out_file.c_str(), "w");
    
    cout << "Max event: " << Max_event << endl;
    cout << "Memory length: " << Max_memory << endl;
    cout << "K: " << K << endl;
    
//Read file and create a sorted list of temporal events
    vector<event> events;
    set<vertex> V;
    createEvents(tmp, events, V);
        
    cout << "# of events:" << events.size() << endl;
    const auto t2 = chrono::steady_clock::now();
    print_time (fp, "Read data time: ", t2 - t1);
    transition TR;
    set<motif> prefixes;
    vector<event> initial_E;
    motif_count MC;
    cout << "calculating transition" << endl;
    for (size_t i=0; i<events.size(); i++) {
        if (i % 1000 == 0) {
            cout << i << " of " << events.size() << endl;
        }
        countMTP(events[i], MC, prefixes, Max_event, Max_memory, initial_E);
    }
    cout << events.size() << " of " << events.size() << endl;
    double new_ratio;
    countTransiton(MC, TR, initial_E.size(), new_ratio);
    int graph_size = countSize(events);
    const auto t3 = chrono::steady_clock::now();
    print_time (fp, "Learn transition time: ", t3 - t2);
    
    cout << "generating graph" << endl;
    vector<event> out_graph;
    string out_time = out_file + "_time";
    vector<event> RIE = RandomizeIE(initial_E, V);
    generateGraph(TR, RIE, V, Max_event, out_graph, K, out_time, graph_size, new_ratio);
    
    const auto t4 = chrono::steady_clock::now();
    print_time (fp, "Generate graph time: ", t4 - t3);
    print_time (fp, "End-to-end Time: ", t4 - t1);
    
    for (int i=0; i<out_graph.size(); i++) {
        fprintf (fp, "%s %s %d\n", out_graph[i].second.first.c_str(), out_graph[i].second.second.c_str(), out_graph[i].first);
    }
    cout << "ALL DONE" << endl;
    
    fclose (fp);
    return 0;
}
