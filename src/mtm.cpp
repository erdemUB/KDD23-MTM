//
//  mtm.cpp
//  mtm
//
//  Created by Penghang Liu on 10/27/21.
//  Copyright © 2021 Penghang Liu. All rights reserved.
//

#include "mtm.hpp"

void createEvents (string filename, vector<event>& events, set<vertex>& V){
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            timestamp t;
            edge e;
            ss >> u >> v >> t;
            if (u != v) {
                V.insert(u);
                V.insert(v);
                e = make_pair(u, v);
                events.push_back(make_pair(t, e));
            }
        }
    }
    sort(events.begin(), events.end());
    events.erase(unique(events.begin(), events.end()),events.end()); //remove duplicates
    return;
}

void countMTP (event e, motif_count& MC, set<motif>& prefixes, int N_event, int d_c, vector<event>& initial_E){
    bool ie {true};
    vertex u = e.second.first;
    vertex v = e.second.second;
    vector<vector<event>> new_motifs;    //used to store the new motifs
    for (auto it = prefixes.begin(); it != prefixes.end();) {   //for each current prefix
        vector<event> pre = *it;
        if ((e.first - pre.back().first <= d_c) && (pre.size() < N_event)) {  //check memory size and the number of events
            set<vertex> nodes = getNodes(pre);
            if (nodes.find(u)!=nodes.end() || nodes.find(v)!=nodes.end()) {
                if (pre.back().first!=e.first) { //check synchronous events
                    vector<event> m_new = pre;
                    m_new.push_back(e);
                    new_motifs.push_back(m_new);
                    string code_new = encodeMotif(m_new);
                    MC[code_new].first += 1;
                    MC[code_new].second.push_back(e.first - pre.back().first);
                    ie = false;
                    it = prefixes.erase(it);
                    continue;
                }
            }
            ++it;
        } else {
            it = prefixes.erase(it);    //remove prefix if it exceeds the delta constrain
        }
    }
    if (ie) {
        initial_E.push_back(e);
        vector<event> E;
        E.push_back(e);
        prefixes.insert(E); // add the new event to the current prefix list
    }
    
    //add the new motifs to the current prefix list
    if (!new_motifs.empty()) {
        for (vector<event> const &mt: new_motifs) {
            prefixes.insert(mt);
        }
    }
    return;
}

string encodeMotif(motif instance){
    string motif_code;
    map<vertex, string> code;
    int i=0;
    unordered_map<timestamp, int> concurrent_count;
    for (auto it=instance.begin(); it!=instance.end(); ++it) {
        vertex u = it->second.first;
        vertex v = it->second.second;
        if (code.find(u)==code.end()) {
            code[u] = to_string(i);
            i++;
        }
        motif_code.append(code[u]);
        if (code.find(v)==code.end()) {
            code[v] = to_string(i);
            i++;
        }
        motif_code.append(code[v]);
    }
    return motif_code;
}

set<vertex> getNodes(vector<event> key){
    set<vertex> nodes;
    for (int i=0; i<key.size(); i++) {
        nodes.insert(key[i].second.first);
        nodes.insert(key[i].second.second);
    }
    return nodes;
}

void countTransiton(motif_count& MC, transition& TR, int base, double& new_ratio){
    double N_motifs = 0;
    double new_edge = 0;
    for (auto it=MC.begin(); it!=MC.end(); ++it) {
        string code = it->first;
        string code_0 = code.substr(0, code.size()-2);
        string code_1 = code.substr(code.size()-2);
        double c = it->second.first;
        int medge = motifEdges(code);
        new_edge += medge * c;
        N_motifs += c;
        vector<timestamp> T = it->second.second;
        double lambda = 1 / (accumulate(T.begin(), T.end(), 0.0)/ T.size());
        TR[code_0].push_back(make_tuple(code_1, c, lambda));
    }
    new_ratio = new_edge / base;
    for (auto it=TR.begin(); it!=TR.end(); ++it) {
        double c_sum=0;
        string head = it->first;
        vector<tuple<string, double, double>> temp = it->second;
        for (int i=0; i<temp.size(); i++) {
            c_sum += get<1>(temp[i]);
        }
        double s_count;
        if (head == "01") {
            s_count = base - c_sum;
        } else {
            s_count = MC[head].first - c_sum;
        }
        it->second.push_back(make_tuple("S", s_count, 0));
    }
    for (auto it=MC.begin(); it!=MC.end(); ++it) {
        string code = it->first;
        if (TR.find(code)==TR.end()) {
            TR[code] = {make_tuple("S", 1, 0)};
        }
    }
    return;
}

void generateGraph(transition TR, vector<event> IE, set<vertex> V, int Max_event, vector<event>& out_graph, int K, string out_time, int graph_size, double new_ratio){
    const auto t0 = chrono::steady_clock::now();
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
   int IE_size = countSize(IE);
//    cout << "IE size: " << IE_size << endl;
   double use_existing = (graph_size - IE_size)/(IE.size() * new_ratio);
//    cout << "use existing prob: " << use_existing << endl;
    // set<edge> current_graph;
    unordered_map<vertex, pair<set<vertex>, set<vertex>>> AJ;
    vector<vertex> node_list{V.begin(),V.end()};
    for (int i=0; i<IE.size(); i++) {
        event e = IE[i];
        out_graph.push_back(e);
        vertex u = e.second.first;
        vertex v = e.second.second;
        timestamp t = e.first;
        // current_graph.insert(make_pair(u,v));
        // current_graph.insert(make_pair(v,u));
        AJ[u].second.insert(v);
        AJ[v].first.insert(u);
        vector<vertex> nodes;
        nodes.push_back(u);
        nodes.push_back(v);
        string code_0 = "01";
        while (code_0.size() < 2*Max_event) {
            vector<tuple<string, double, double>> temp = TR[code_0];
            vector<double> PROB;
            for (int j=0; j<temp.size(); j++) {
                PROB.push_back(get<1>(TR[code_0][j]));
            }
            discrete_distribution<> D(PROB.begin(), PROB.end());
            int r = D(generator);
            string code_1 = get<0>(TR[code_0][r]);
            double lambda = get<2>(TR[code_0][r]);
            if(code_1=="S") break;
            int n1 = code_1[0] - '0';
            int n2 = code_1[1] - '0';
            vertex u_new, v_new;
            if (n1 >= nodes.size()) {
                v_new = nodes[n2];
                // int current_size = current_graph.size()/2;
                // double use_existing = (graph_size - current_size)/((IE.size()-i) * new_ratio);
                binomial_distribution<int> b_distribution(1,1-use_existing);
                int new_or_not = b_distribution(generator);
                if (new_or_not && AJ[v_new].first.size()) { // use existing
                    vector<vertex> pick_nodes{AJ[v_new].first.begin(),AJ[v_new].first.end()};
                    int pick = rand() % pick_nodes.size();
                    u_new = pick_nodes[pick];
                } else {
                    do {
                        int pick = rand() % node_list.size();
                        u_new = node_list[pick];
                    } while (find(nodes.begin(),nodes.end(),u_new)!=nodes.end());
                }
                nodes.push_back(u_new);
            } else if (n2 >= nodes.size()){
                u_new = nodes[n1];
                // int current_size = current_graph.size()/2;
                // double use_existing = (graph_size - current_size)/((IE.size()-i) * new_ratio);
                binomial_distribution<int> b_distribution(1,1-use_existing);
                int new_or_not = b_distribution(generator);
                if (new_or_not && AJ[u_new].second.size()) { // use existing
                    vector<vertex> pick_nodes {AJ[u_new].second.begin(), AJ[u_new].second.end()};
                    int pick = rand() % pick_nodes.size();
                    v_new = pick_nodes[pick];
                } else {
                    do {
                        int pick = rand() % node_list.size();
                        v_new = node_list[pick];
                    } while (find(nodes.begin(),nodes.end(),v_new)!=nodes.end());
                }
                nodes.push_back(v_new);
            } else {
                u_new = nodes[n1];
                v_new = nodes[n2];
            }
            exponential_distribution<> P(lambda);
            timestamp delta_t = P(generator);
            t+=delta_t;
            event new_e = make_pair(t, make_pair(u_new, v_new));
            out_graph.push_back(new_e);
            // current_graph.insert(make_pair(u_new,v_new));
            // current_graph.insert(make_pair(v_new,u_new));
            AJ[u_new].second.insert(v_new);
            AJ[v_new].first.insert(u_new);
            code_0.append(code_1);
        }
    }
    return;
}

vector<event> RandomizeIE(vector<event> initial_E, set<vertex> V){
    vector<event> Randomized_IE;
    unordered_map<edge, vector<timestamp>> IE_projection;
    vector<int> degrees;
    unordered_map<vertex, int> degreeCount;
    set<vertex> IE_nodes;
    vector<vertex> graph_nodes(V.begin(), V.end());

    random_device rd;
    mt19937 g(rd());
    shuffle(graph_nodes.begin(), graph_nodes.end(), g);

    for (const auto& ie : initial_E) {
        IE_projection[ie.second].push_back(ie.first);
    }
    
    for (const auto& e : IE_projection) {
        degreeCount[e.first.first]++;
        degreeCount[e.first.second]++;
    }

    for (auto it=degreeCount.begin(); it!=degreeCount.end(); ++it) {
        degrees.push_back(it->second);
    }

    vector<int> nodes;
    for (int i = 0; i < degrees.size(); i++) {
        for (int j = 0; j < degrees[i]; j++) {
            nodes.push_back(i);
        }
    }

    random_device rd1;
    mt19937 g1(rd1());
    shuffle(nodes.begin(), nodes.end(), g1);

    int index = 0;
    for (const auto& e : IE_projection) {
        for (auto t=e.second.begin(); t!=e.second.end(); ++t){
            Randomized_IE.push_back({*t, {graph_nodes[nodes[index]], graph_nodes[nodes[index + 1]]}});
        }
        index += 2;
    }

    return Randomized_IE;
}

int countSize(vector<event> input){
    set<edge> temp;
    for (auto it=input.begin(); it!=input.end(); ++it) {
        edge e = it->second;
        vertex u = e.first;
        vertex v = e.second;
        temp.insert(make_pair(u,v));
        temp.insert(make_pair(v,u));
    }
    return temp.size()/2;
}

int motifEdges(string code){
    int max_e = 0;
    for (int i=0; i<code.size(); i++) {
        int temp = code[i] - '0';
        max_e = max_e > temp ? max_e : temp;
    }
    return max_e - 1;
}
