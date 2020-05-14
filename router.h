#ifndef ROUTER_H
#define ROUTER_H

#include<iostream>
#include<array>



#define MAX_CHAN 5



class Router
{
public:
    Router();
    void calc_mat_fwd_prob();
    void calc_mat_contention();
    void calc_mat_delta();
    void calc_mat_ResServTime();
    void calc_mat_BuffUtil();
    void calc_mat_waitTime();
    void calc_serv_time();
    void Reset();

    float arr_rate[MAX_CHAN];
    float mat_lambda[MAX_CHAN][MAX_CHAN];
    float mat_delta[MAX_CHAN][MAX_CHAN];
    float mat_Contention[MAX_CHAN][MAX_CHAN];
    float mat_resServTime[MAX_CHAN];
    float mat_forwProb[MAX_CHAN][MAX_CHAN];
    float mat_avgBufUtil[MAX_CHAN];
    float mat_avgWaitTime[MAX_CHAN];
    float sum_avgBufUtil;
    float T;
    float serv_time[MAX_CHAN];

    int out_index[MAX_CHAN];
    int in_index[MAX_CHAN];
    int num_eject;
    int num_path;


    struct stCoord{
        int Y;
        int X;
    }CoOrd;

    Router *North;
    Router *East;
    Router *West;
    Router *South;
};

#endif // ROUTER_H
