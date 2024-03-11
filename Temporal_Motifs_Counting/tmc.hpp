//
//  tmc.hpp
//  tmc
//
//  Created by Penghang Liu on 10/27/21.
//  Copyright © 2021 Penghang Liu. All rights reserved.
//

#ifndef tmc_hpp
#define tmc_hpp

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <initializer_list>

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
typedef vector<event> key;  //a prefix or a motif
typedef pair<int, set<vertex>> counts;  //the count of the (prefix/motif) and the vertices in the (prefix/motif)
//typedef unordered_map<vector<event>, set<vertex>> prefix;
typedef unordered_map<vector<event>, pair<int, set<vertex>>> instancemap; //a hashtable of key and counts

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

void createEvents (string filename, vector<event>& events); //Load and sort the event list
void countInstance (event e, instancemap& imap, set<vector<event>>& keys, int N_event, int d_c, string consecutive);    //Increment the instance count and update the prefix type
string encodeMotif(vector<event> instance); //identify the type of motif
void countMotif (event e, set<key>& pre, map<string, int>& motif_count, int N_vtx, int N_event, int d_c, int d_w);
set<vertex> getNodes(vector<event> key);
void countSpecificmotif (event e, set<key>& pre, int& motif_count, string code_given, int N_vtx, int N_event, int d_c, int d_w);
char sconvert (int i);
#endif /* tmc_hpp */
