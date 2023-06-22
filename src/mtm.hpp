//
//  mtm.hpp
//  mtm
//
//  Created by Penghang Liu on 10/27/21.
//  Copyright © 2021 Penghang Liu. All rights reserved.
//

#ifndef mtm_hpp
#define mtm_hpp

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <set>
#include <array>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <initializer_list>
#include <tuple>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <chrono>
#include <sys/stat.h>

using namespace std;

typedef chrono::duration<double> tms;
typedef string vertex;
typedef int timestamp;
typedef pair<vertex, vertex> edge;
typedef pair<timestamp, edge> event;
typedef vector<event> motif;
typedef unordered_map<string, pair<double, vector<timestamp>>> motif_count;
typedef unordered_map<string, vector<tuple<string, double, double>>> transition;
typedef pair<int, set<vertex>> counts;  
typedef unordered_map<vector<event>, pair<int, set<vertex>>> instancemap; 

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
    template<typename S, typename T> struct hash<pair<S, T>>
    {
        inline size_t operator()(const pair<S, T> & v) const
        {
            size_t seed = 0;
            ::hash_combine(seed, v.first);
            ::hash_combine(seed, v.second);
            return seed;
        }
    };
    
    template<typename S, typename T> struct hash<vector<pair<S, T>>>
    {
        inline size_t operator()(const vector<pair<S, T>> &k) const
        {
            size_t seed = 0;
            for (auto it=k.begin(); it != k.end(); ++it) {
                hash_combine(seed, it->first);
                hash_combine(seed, it->second);
            }
            return seed;
        }
    };
}

inline void print_time (FILE* fp, const string& str, tms t) {
    fprintf (fp, "%s %.6lf\n", str.c_str(), t.count());
    fflush(fp);
}

void createEvents (string filename, vector<event>& events, set<vertex>& V); //Load and sort the event list
string encodeMotif(vector<event> instance); //identify the type of motif
void countMTP (event e, motif_count& MC, set<motif>& prefixes, int N_event, int d_c, vector<event>& initial_E);
set<vertex> getNodes(vector<event> key);
vector<event> RandomizeIE(vector<event> initial_IE, set<vertex> V);
void countTransiton(motif_count& MC, transition& TR, int base, double& new_ratio);
void generateGraph(transition TR, vector<event> IE, set<vertex> V, int Max_event, vector<event>& out_graph, int K, string out_time, int graph_size, double new_ratio);
int countSize(vector<event> input);
int motifEdges(string code);
#endif /* mtm_hpp */
