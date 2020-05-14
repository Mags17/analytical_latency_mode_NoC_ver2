
/********************************************************************************************
 This project implements the paper "An analytical approach for NoC performance Analysis" by
Umit Y. Orgas, Radu Marculecu et. al published in IEEE/TCAD in december 2010.


This is an analytical model for latency prediction in NoCs

Author: Hitesh Manghnani
First Implementation: Nov 2019
***********************************************************************************************/


#include <iostream>
#include <armadillo>
#include <vector>
#include <fstream>
#include <string>
#include <assert.h>
#include <map>
#include <chrono>


#include "globals.h"
#include "router.h"

using namespace std;
using namespace arma;

vector <Router *> R;

void Ver_date()
{
    string date = "24th April 2020";
    string version = "3.0";

    cout << "Version: " << version << endl;
    cout << "Date: " << date << endl;

    //sleep(5);
}

void calc_checksum(string prog_name)
{
    ifstream infile;
    char c;
    int sum = 0;

    infile.open(prog_name, std::fstream::in | std::fstream::binary );

     while(infile.get(c))
        sum += c;

    cout  << "Checksum: " << std::hex << sum << std::dec << endl;
    sleep(5);
}

// Read Next value from input config file
string readNextVal(ifstream & fin)
{
    string line;
    std::size_t pos;

    getline(fin,line);
    pos = line.find("#");
    assert(pos != std::string::npos);
    return(line.substr(pos+1));
}


// read the config file
void readConfigFile()
{
    ifstream file_config;
    string temp_line;
    stringstream line_in;

    file_config.open("files/config.txt");

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> num_rows;
    line_in.clear();        // for clearing fail state of stringstream as it is empty now
    //file_config >>  num_rows;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> num_cols;
    line_in.clear();        // for clearing fail state of stringstream as it is empty now
    //file_config >> num_cols;
    size_mesh = num_cols * num_rows;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> str_mapFile;
    line_in.clear();        // for clearing fail state of stringstream as it is empty now
    //file_config >> str_mapFile;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> str_trafficFile;
    line_in.clear();        // for clearing fail state of stringstream as it is empty now
   //file_config >> str_trafficFile;

/*    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> inj_rate;
    line_in.clear();
    //file_config >> inj_rate;

/*    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> num_mappings;
    line_in.clear();
    //file_config >> num_mappings;*/

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> delay_routing;
    line_in.clear();
    //file_config >> delay_routing;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> delay_link;
    line_in.clear();
   // file_config >> delay_link;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> pkt_size;
    line_in.clear();
    //file_config >> serial_delay;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> buf_size;
    line_in.clear();
    //file_config >> buf_size;

    temp_line = readNextVal(file_config);
    line_in.str(temp_line);
    line_in >> bPrint_features;
    line_in.clear();

    file_config.close();
}


// function not used now
void readMappingFile(int mapping_num)
{    
    string line;

    mapping = new int[size_mesh];

    getline(file_mapping,line);

    string line_temp = line;
    file_out << line_temp;

    std::remove(line.begin(), line.end(), ',');
    std::remove(line.begin(), line.end(), '[');
    std::remove(line.begin(), line.end(), ']');

    stringstream line_in(line);

    line_in >> mapping_num;
    cout << "Next mapping - " << mapping_num <<":\t";

    for(int i=0;i<size_mesh;i++)
    {
        line_in>>mapping[i];
        cout<<mapping[i]<<" ";

        mapping[i]--;
        task_map.insert({mapping[i],i});
    }
    line_in >> inj_rate;
    cout << endl;
}

int readMappingLine(string line)
{
    int mapping_num;
    string line_temp = line;
    file_out << line_temp;

    std::remove(line.begin(), line.end(), ',');
    std::remove(line.begin(), line.end(), '[');
    std::remove(line.begin(), line.end(), ']');

    stringstream line_in(line);

    line_in >> mapping_num;
    cout << "Next mapping - " << mapping_num <<":\t";

    for(int i=0;i<size_mesh;i++)
    {
        line_in>>mapping[i];
        cout<<mapping[i]<<" ";

        mapping[i]--;
        task_map.insert({mapping[i],i});
    }
    line_in >> inj_rate;
    line_in >> booksim_latency;
    cout << endl;

    return(mapping_num);
}


void readTrafficFile()
{
    ifstream file_traffic;
    int num_nodes, max_traffic = 0;
    string str_in;

    file_traffic.open(str_trafficFile);
    file_traffic >> num_nodes;

    assert(num_nodes == size_mesh);

    lambda_src= new float[size_mesh];       // injection rate for each source
    lambda_src_temp= new float[size_mesh];       // injection rate for each source
    traffic.bw_src= new float[size_mesh];
    traffic.bw_s2d= new float*[size_mesh];
    for(int i = 0; i < num_nodes; i++)        
        traffic.bw_s2d[i]= new float[size_mesh];

    for(int row = 0; row < num_nodes; row++)
    {
        traffic.bw_src[row] = 0.0;
        for(int col = 0; col < num_nodes; col++)
        {
            file_traffic >> str_in;

            if("INF" == str_in)                
                traffic.bw_s2d[row][col] = 0.0;
            else
            {
                traffic.bw_s2d[row][col] = stof(str_in);
                traffic.bw_src[row] += traffic.bw_s2d[row][col];
            }
        }

        if(traffic.bw_src[row] > max_traffic)
            max_traffic = traffic.bw_src[row];
    }

    for(int i = 0; i < num_nodes; i++)
    {
        lambda_src_temp[i] = traffic.bw_src[i] / max_traffic;    // This is not final value, lambda_src is updated later
    }

    file_traffic.close();

    // debug code
/*    for(int row = 0; row < num_nodes; row++)
    {
        for(int col = 0; col < num_nodes; col++)
        {
            cout << traffic.bw_s2d[row][col] << " ";
        }
        cout << endl;
    }
*/
/*    for(int row = 0; row < num_nodes; row++)
    {
        cout << traffic.bw_src[row] << " ";
    }
    cout << endl;
*/
/*
    cout << max_traffic;

    for(int i = 0; i < num_nodes; i++)
    {
        cout <<  lambda_src[i]  << " ";
    }
*/
}

vector<Router *> createMeshNw(vector <Router *>R)
{
    R.resize(size_mesh);

    for(int i = 0;i < size_mesh; i++)
    {
        R[i] = new Router();
        R[i]->CoOrd.Y = i/num_cols; // same for a row
        R[i]->CoOrd.X = i%num_cols; // same for a column
    }

    for(int i = 0; i < size_mesh; i++)
    {
        if (0 == i)         // for bottomleft corner
        {            
            R[i]->North = R[i + num_cols];
            R[i]->East = R[i+1];
        }
        else if(num_cols - 1 == i)          // for bottomright corner
        {
            R[i]->West = R[i-1];
            R[i]->North = R[i + num_cols];
        }
        else if(size_mesh - num_cols == i)  // for topleft corner
        {
            R[i]->East = R[i+1];
            R[i]->South = R[i - num_cols];
        }
        else if (size_mesh - 1 == i)    // for topright corner
        {
            R[i]->West = R[i-1];
            R[i]->South = R[i - num_cols];
        }
        else if(0 == (i % num_cols))     // for leftmost column
        {
            R[i]->East = R[i+1];
            R[i]->North = R[i + num_cols];
            R[i]->South = R[i - num_cols];
        }
        else if(0 == ((i+1) % num_cols))    // for rightmost column
        {
            R[i]->West = R[i-1];
            R[i]->North = R[i + num_cols];
            R[i]->South = R[i - num_cols];

        }
        else if(0 == (i/num_cols))  // for bottommost row
        {
            R[i]->West = R[i-1];
            R[i]->North = R[i + num_cols];
            R[i]->East = R[i+1];
        }
        else if((num_rows-1) == (i/num_cols))      // for topmost row
        {
            R[i]->West = R[i-1];
            R[i]->South = R[i - num_cols];
            R[i]->East = R[i+1];
        }
        else            // for inner routers
        {
            R[i]->East = R[i+1];
            R[i]->West = R[i-1];
            R[i]->North = R[i + num_cols];
            R[i]->South = R[i - num_cols];
        }

    }
    return(R);
}

vector <stDir> path(const int task_src,const int task_dest)
{
    stDir stTemp;
    vector <stDir> route_path;   // 1st param = router, 2nd param = direction
    bool goSouth,goNorth, goEast, goWest;
    int n_Rsrc = task_map.at(task_src);
    int n_Rdest = task_map.at(task_dest);

    goSouth = goNorth = goEast = goWest = false;

    Router *R_src = R.at(n_Rsrc);
    Router *R_dest = R.at(n_Rdest);


    if((R_dest->CoOrd.Y - R_src->CoOrd.Y) > 0)
        goNorth = true;
    else if ((R_dest->CoOrd.Y - R_src->CoOrd.Y) < 0)
        goSouth = true;

    if((R_dest->CoOrd.X - R_src->CoOrd.X) > 0)
        goEast = true;
    else if ((R_dest->CoOrd.X - R_src->CoOrd.X) < 0)
        goWest = true;

    int yDiff = abs(R_dest->CoOrd.Y - R_src->CoOrd.Y);
    int xDiff = abs(R_dest->CoOrd.X - R_src->CoOrd.X);

    //route_path.push_back(n_Rsrc);
    Router *R_curr = R_src;
    int n_Rcurr = n_Rsrc;
    stTemp.ip_dir = Channel::Local;

    if(goEast)
    {
        for(int i =0;i < xDiff;i++)
        {            

            if(xDiff + yDiff - i > R[n_Rcurr]->out_index[Channel::East])
            {
                R[n_Rcurr]->out_index[Channel::East] = xDiff + yDiff - i;
                R[n_Rcurr]->East->in_index[Channel::West] = xDiff + yDiff - i;
            }

            if((n_Rcurr == n_Rsrc) && (R[n_Rcurr]->out_index[Channel::East] >= R[n_Rcurr]->in_index[Channel::Local]))
                R[n_Rcurr]->in_index[Channel::Local] = R[n_Rcurr]->out_index[Channel::East] + 1;


            R[n_Rcurr]->num_path++;
            stTemp.router_idx = n_Rcurr;
            stTemp.op_dir = Channel::East;
            route_path.push_back(stTemp);

            stTemp.ip_dir = Channel::West;
            n_Rcurr++;
        }
    }
    else if(goWest)
    {
        for(int i =0;i < xDiff;i++)
        {           
            if(xDiff + yDiff - i > R[n_Rcurr]->out_index[Channel::West])
            {
                R[n_Rcurr]->out_index[Channel::West] = xDiff + yDiff - i;
                R[n_Rcurr]->West->in_index[Channel::East] = xDiff + yDiff - i;
            }

           if((n_Rcurr == n_Rsrc) && (R[n_Rcurr]->out_index[Channel::West] >= R[n_Rcurr]->in_index[Channel::Local]))
                R[n_Rcurr]->in_index[Channel::Local] = R[n_Rcurr]->out_index[Channel::West]  + 1;

            R[n_Rcurr]->num_path++;
            stTemp.router_idx = n_Rcurr;
            stTemp.op_dir = Channel::West;
            route_path.push_back(stTemp);

            stTemp.ip_dir = Channel::East;
            n_Rcurr--;
        }
    }

    if(goNorth)
    {
        for(int i =0;i < yDiff;i++)
        {            
            if(yDiff - i > R[n_Rcurr]->out_index[Channel::North])
            {
                R[n_Rcurr]->out_index[Channel::North] = yDiff - i;
                R[n_Rcurr]->North->in_index[Channel::South] = yDiff - i;
            }


           if((n_Rcurr == n_Rsrc) && (R[n_Rcurr]->out_index[Channel::North] >= R[n_Rcurr]->in_index[Channel::Local]))
                R[n_Rcurr]->in_index[Channel::Local] =  R[n_Rcurr]->out_index[Channel::North] + 1;


            R[n_Rcurr]->num_path++;
            stTemp.router_idx = n_Rcurr;
            stTemp.op_dir = Channel::North;
            route_path.push_back(stTemp);

            stTemp.ip_dir = Channel::South;
            n_Rcurr+=num_cols;
        }
    }
    else if(goSouth)
    {
        for(int i =0;i < yDiff;i++)
        {           
            if(yDiff - i > R[n_Rcurr]->out_index[Channel::South])
            {
                R[n_Rcurr]->out_index[Channel::South] = yDiff - i;
                R[n_Rcurr]->South->in_index[Channel::North] = yDiff - i;
            }

           if((n_Rcurr == n_Rsrc) && (R[n_Rcurr]->out_index[Channel::South] >= R[n_Rcurr]->in_index[Channel::Local]))
                R[n_Rcurr]->in_index[Channel::Local] = R[n_Rcurr]->out_index[Channel::South] + 1;

            R[n_Rcurr]->num_path++;
            stTemp.router_idx = n_Rcurr;
            stTemp.op_dir = Channel::South;
            route_path.push_back(stTemp);

            stTemp.ip_dir = Channel::North;
            n_Rcurr-=num_cols;
        }
    }

    R[n_Rcurr]->out_index[Channel::Local] = 0;

    R[n_Rcurr]->num_eject++;
    stTemp.router_idx = n_Rcurr;
    stTemp.op_dir = Channel::Local;
    route_path.push_back(stTemp);


/*   for(int i = 0; i < route_path.size();i++)
          cout << route_path.at(i).ip_dir<< "-" << route_path.at(i).router_idx << "-" << route_path.at(i).op_dir << " ";

   cout << endl;
*/
   return(route_path);
}

void update_router_arr_rates(const vector<stDir> flow)
{
  /*  for(int i = 0; i < flow.size();i++)
        cout << flow.at(i) << " ";

    cout << endl;
    */

    int R_src,R_dest, R_curr;
    int task_src,task_dest;
    int ip_dir, op_dir;

    R_src = flow.front().router_idx;
    R_dest = flow.back().router_idx;
    task_src = mapping[R_src];
    task_dest = mapping[R_dest];

   // cout << R_src << " " << R_dest << endl;

    for(unsigned i = 0; i < flow.size();i++)
    {
        R_curr = flow.at(i).router_idx;
        ip_dir = flow.at(i).ip_dir;
        op_dir = flow.at(i).op_dir;

        assert(traffic.bw_src[task_src] != 0);
        assert(ip_dir != op_dir);

        R[R_curr]->mat_lambda[ip_dir][op_dir] += lambda_src[task_src] * traffic.bw_s2d[task_src][task_dest] / traffic.bw_src[task_src];
        R[R_curr]->arr_rate[ip_dir] +=  lambda_src[task_src] * traffic.bw_s2d[task_src][task_dest] / traffic.bw_src[task_src];
    }
}

void sortChannels()
{
    int chan_id;
    int index_no;

    for (int i = 0 ; i < size_mesh ; i++)
    {
        for(int j=0; j< MAX_CHAN;j++)
        {
            chan_id = MAX_CHAN*i + j;
            index_no = R[i]->in_index[j];
            channel_index_map.insert({chan_id,index_no});
            channel_index_multimap.insert({index_no,chan_id});
        }
    }

}

void update_service_times()
{
    std::multimap<int,int>::iterator iter;
    uint router_no;
    uint channel_no;
    uint channel_index;
    uint abs_chan_no;
    float forw_prob;
    //float buf_ratio = (float)(buf_size)/pkt_size;

    for (iter=channel_index_multimap.begin(); iter!=channel_index_multimap.end(); ++iter)
    {
        //std::cout << (*iter).first << " => " << (*iter).second << '\n';
        abs_chan_no = (*iter).second;
        router_no = abs_chan_no/ MAX_CHAN;
        channel_no = abs_chan_no % MAX_CHAN;
        channel_index = (*iter).first;

        Router *R_temp;
        float sum = 0;

        if(0 == channel_index)
        {

        }
        else if(1 == channel_index)
        {
            R[router_no]->serv_time[channel_no] = service_time-1;
        }
        else
        {
            for(int i = 0; i < MAX_CHAN; i++)
            {
                if( i != channel_no)
                {
                    forw_prob = R[router_no]->mat_forwProb[channel_no][i];
                    switch(i)
                    {
                    case 0:
                        {
                            sum += forw_prob * (service_time-1);
                            break;
                        }
                    case 1:
                        {
                            R_temp = R[router_no]->North;

                            if(NULL != R_temp)
                            {
                                sum += forw_prob * (Hs + R_temp->mat_avgWaitTime[Channel::South] + R_temp->serv_time[Channel::South]
                                                      - buf_size*1);
                            }
                            break;
                        }
                    case 2:
                        {
                            R_temp = R[router_no]->East;

                            if(NULL != R_temp)
                            {
                                sum += forw_prob * (Hs + R_temp->mat_avgWaitTime[Channel::West] + R_temp->serv_time[Channel::West]
                                                      - buf_size*1);
                            }
                            break;
                        }
                    case 3:
                        {
                            R_temp = R[router_no]->South;

                            if(NULL != R_temp)
                            {
                                sum += forw_prob * (Hs + R_temp->mat_avgWaitTime[Channel::North] + R_temp->serv_time[Channel::North]
                                                      - buf_size*1);
                            }
                            break;
                        }
                    case 4:
                        {
                            R_temp = R[router_no]->West;

                            if(NULL != R_temp)
                            {
                                sum += forw_prob * (Hs + R_temp->mat_avgWaitTime[Channel::East] + R_temp->serv_time[Channel::East]
                                                      - buf_size*1);
                            }
                            break;
                        }
                    default:
                        {
                            assert(i < 5);
                            break;
                        }

                    } // End of switch(i)

                } // End of if( i != channel_no)
            } // End of for(int i = 0; i < MAX_CHAN; i++)

           // assert (sum >= 0);
            if(sum > R[router_no]->serv_time[channel_no])
                R[router_no]->serv_time[channel_no] = sum;
        } // End of if-else of channel_indexes

    } // End of for loop for multimap of channelindexes
}

void calc_Ws(float *Ws)
{
    float serv_time_src;

    serv_time_src = delay_link + serial_delay;
    //serv_time_src = Hs + serial_delay;

    for(int i = 0; i<size_mesh;i++)
    {
        Ws[i] = (lambda_src[i] * serv_time_src * serv_time_src) / (2 *(1 - lambda_src[i] * serv_time_src)) ;
        assert (Ws[i] >= 0);
        //cout << Ws[i] << endl;
    }
}

vector <float> calc_flow_latency(const vector <vector <stDir>> flows)
{
    vector <stDir> route_path;      // represents one flow
    vector <float> latency_flow;
    float temp_latency;
    float *Ws;

    Ws = new float[size_mesh];
    calc_Ws(Ws);
    for(unsigned i = 0; i < flows.size();i++)
    {
        int R_curr, R_src, task_src;
        int temp_idx;

        route_path.clear();
        route_path = flows.at(i);
        temp_latency = 0.0;        
        R_src = route_path.front().router_idx;
        task_src = mapping[R_src];

        for(unsigned j = 0;j < route_path.size();j++)
        {
            R_curr = route_path.at(j).router_idx;
            temp_idx = route_path.at(j).ip_dir;

            //temp_latency += R[R_curr]->mat_avgWaitTime[temp_idx] + Hs;
            //temp_latency += R[R_curr]->mat_avgWaitTime[temp_idx] + R[R_curr]->T - serial_delay;
            temp_latency += R[R_curr]->mat_avgWaitTime[temp_idx] + R[R_curr]->serv_time[temp_idx] - serial_delay;
        }
        temp_latency += Ws[task_src];
        temp_latency += serial_delay;       

        latency_flow.push_back(temp_latency);
    }
    return(latency_flow);
}

float calc_avg_pkt_latency(vector <float> latency_flow, vector <vector <stDir>> flows)
{
    float avg_latency = 0.0;
    float total_bw = 0.0;
    vector <stDir> route_path;      // represents one flow

    for(int i=0;i<size_mesh;i++)
        total_bw += traffic.bw_src[i];

    for(unsigned i = 0; i < flows.size();i++)
    {
        int R_src, R_dest, task_src, task_dest;
        int temp_idx;

        route_path.clear();
        route_path = flows.at(i);

        R_src = route_path.front().router_idx;
        R_dest = route_path.back().router_idx;
        task_src = mapping[R_src];
        task_dest = mapping[R_dest];

        avg_latency += (traffic.bw_s2d[task_src][task_dest] / total_bw) * latency_flow.at(i);
    }
    return(avg_latency);

}

void print_features(float avg_pkt_latency, int map_num)
{
    stringstream str_out("");

    str_out << " " << delay_routing << " " << delay_link << " " << serial_delay << " ";

    for(int r_num = 0;r_num < size_mesh;r_num++)
    {
        str_out << "[";

        for(int channel = 0; channel < MAX_CHAN; channel++)
            str_out << R[r_num]->arr_rate[channel] << " ";

        str_out << "]" << " ";
    }

    for(int r_num = 0;r_num < size_mesh;r_num++)
    {
        str_out << "[";

        for(int i = 0; i < MAX_CHAN; i++)
            for(int j = 0; j < MAX_CHAN; j++)
                str_out <<  R[r_num]->mat_forwProb[i][j] << " ";

        str_out << "]" << " ";
    }

    for(int r_num = 0;r_num < size_mesh;r_num++)
    {
        str_out << "[";

        for(int i = 0; i < MAX_CHAN; i++)
            str_out <<  R[r_num]->mat_avgBufUtil[i] << " ";

        str_out << "]" << " ";
    }

    for(int r_num = 0;r_num < size_mesh;r_num++)
    {
          str_out <<  R[r_num]->T << " ";
    }

    file_out << str_out.str();
    file_out << endl;

    file_feature << map_num << " ";
    file_feature << inj_rate << " ";
    file_feature << str_out.str();
    file_feature << avg_pkt_latency << " ";
    file_feature << booksim_latency;
    file_feature << endl;
}

bool checkForSaturation()
{
    for(int i = 0; i < size_mesh; i++)
    {
        for (int k = 0 ; k< MAX_CHAN;k++)
        {
            if((R[i]->mat_avgBufUtil[k] < 0))// || (R[i]->mat_avgBufUtil[k] >= (float)buf_size/pkt_size))
                flag_saturate = true;
        }

        //if(R[i]->sum_avgBufUtil >=1.0)
        //    flag_saturate = true;
    }

    return flag_saturate;
}

void print_AvgBufUtil(int map_num)
{
    ofstream file[16];

    for(int i = 0; i < 16; i++)
    {
        string name = "files/bufUtil_" + to_string(i) + ".txt";
        file[i].open(name,std::ofstream::out | std::ofstream::app);

        file[i] << map_num << " " << inj_rate << " ";

        for (int k = 0 ; k< MAX_CHAN;k++)
            file[i] << R[i]->mat_avgBufUtil[k] << " ";

        file[i] << "\t" << R[i]->sum_avgBufUtil << endl;

        file[i].close();
    }


}

int main(int argc, char **argv)
{   
    vector <stDir> route_path;      // represents one flow
    vector <vector <stDir>> flows;  // represents vector of all flows
    vector <float> latency_flow;    // vector containing latencies for all flows
    float avg_pkt_latency;
    int mapping_num;
    string line;
    string program = argv[0];

    using namespace std::chrono;

    high_resolution_clock::time_point g_start_time = high_resolution_clock::now();

    Ver_date();
    calc_checksum(program);

    flag_saturate = false;

    readConfigFile();

    serial_delay = pkt_size - 1;
    Hs = delay_routing + delay_link;
    service_time = Hs + serial_delay;

    /**** derive output file name ****/
    string str_fileIn = str_mapFile;//dataset_file_1.txt";
    string str_fileOut = str_fileIn;
    string str_filefeature = str_fileIn;

    size_t len = str_fileIn.size();
    str_fileOut.resize(len-4);
    str_filefeature.resize(len-4);
    str_fileOut.append("_out.txt");
    str_filefeature.append("_features.txt");

    file_mapping.open(str_mapFile);
    file_out.open(str_fileOut);

    if(bPrint_features)
        file_feature.open(str_filefeature);
    /************************************/

    readTrafficFile();
    R = createMeshNw(R);

    mapping = new int[size_mesh];

    while(getline(file_mapping,line))
    //for(mapping_num = 1; mapping_num <= num_mappings ; mapping_num++)
    {
        flag_saturate = false;

        negative_flag = false;

        mapping_num = readMappingLine(line);

        for(int i = 0; i < size_mesh; i++)
        {
            lambda_src[i] = inj_rate *  lambda_src_temp[i];
           // cout << lambda_src[i] << " ";
        }

        for(int task_src = 0; task_src < size_mesh; task_src++)
        {
            for(int task_dest = 0; task_dest < size_mesh; task_dest++)
            {
                if((task_src != task_dest) && (traffic.bw_s2d[task_src][task_dest]))
                {
                    route_path = path(task_src,task_dest);
                    flows.push_back(route_path);
                }
            }
        }

        for(int i = 0; i < flows.size();i++)
        {
            update_router_arr_rates(flows.at(i));
        }

        for(int r_num = 0;r_num < size_mesh;r_num++)
        {
            R[r_num]->calc_serv_time();
            R[r_num]->calc_mat_fwd_prob();
            R[r_num]->calc_mat_contention();
            R[r_num]->calc_mat_delta();
            R[r_num]->calc_mat_ResServTime();
            R[r_num]->calc_mat_BuffUtil();
            R[r_num]->calc_mat_waitTime();
        }

        sortChannels();

        update_service_times();

        for(int r_num = 0;r_num < size_mesh;r_num++)
        {
            R[r_num]->calc_mat_ResServTime();
            R[r_num]->calc_mat_BuffUtil();
            R[r_num]->calc_mat_waitTime();
        }

        latency_flow = calc_flow_latency(flows);

        avg_pkt_latency = calc_avg_pkt_latency(latency_flow, flows);


        if(checkForSaturation())
            avg_pkt_latency = 500.0;

        if (avg_pkt_latency > 500.0)
            avg_pkt_latency = 500.0;

        cout << avg_pkt_latency << endl;
        file_out << "\t" << avg_pkt_latency;

        //print_AvgBufUtil(mapping_num);

        if(bPrint_features)
            print_features(avg_pkt_latency, mapping_num);
        else
            file_out << endl;

        route_path.clear();
        flows.clear();
        latency_flow.clear();
        avg_pkt_latency = 0;

        task_map.clear();
        channel_index_map.clear();
        channel_index_multimap.clear();

        //delete mapping;

        for(int r_num = 0;r_num < size_mesh;r_num++)
        {
            R[r_num]->Reset();
        }       
    }

    delete mapping;

    file_mapping.close();
    file_out.close();
    if(bPrint_features)
        file_feature.close();

    high_resolution_clock::time_point g_end_time = high_resolution_clock::now();

    duration<double> g_time_span = duration_cast<duration<double>>(g_end_time - g_start_time);

    cout<<"Total Run Time "<<g_time_span.count()<< "Seconds" << endl;

    return 0;
}
