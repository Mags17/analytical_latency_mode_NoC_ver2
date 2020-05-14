#ifndef GLOBALS_H
#define GLOBALS_H


#include <string>
#include <map>
#include <vector>

using namespace std;

enum Channel {Local,North,East,South,West};

ifstream file_mapping;
ofstream file_out;
ofstream file_feature;

int num_rows;
int num_cols;
int size_mesh;
int num_mappings;

int * mapping;      // indexes are router numbers, values are task numbers

float Hs;
float delay_routing;
float delay_link;
float serial_delay;
float pkt_size;
float buf_size;
float service_time;
float inj_rate;      // global injection rate
float booksim_latency;  // latency calculated by booksim, taken from mapping file
float *lambda_src_temp;      // injection rate for each source
float *lambda_src;      // injection rate for each source
//float ** graph_traffic;
bool bPrint_features;

struct stTraffic{
    float **bw_s2d;     // src to dest BW
    float *bw_src;      // bw for src only i.e sum for all dest
}traffic;

struct stDir{
    int router_idx;
    int ip_dir;
    int op_dir;
};

//struct channel_id{
//    int Router_no;
//    int chan_no;
//};

string str_mapFile;
string str_trafficFile;

map <int,int> task_map;  // 1st param = task, 2nd param = router
map <int,int> channel_index_map;
multimap <int,int> channel_index_multimap;

bool negative_flag;
bool flag_saturate;


#endif // GLOBALS_H

