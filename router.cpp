//#include "globals.h"
#include "router.h"
#include <assert.h>
#include <armadillo>


using namespace std;
using namespace arma;

//extern enum Channel ch;
extern float service_time;
extern bool negative_flag;

Router::Router() :arr_rate(), mat_Contention(), mat_forwProb(),
    mat_lambda(), mat_delta(), mat_resServTime(), mat_avgWaitTime(), serv_time(), out_index({100,0,0,0,0}), in_index()
{
/*
    North = nullptr;
    East = nullptr;
    West = nullptr;
    South = nullptr;
    CoOrd.X = 0;
    CoOrd.Y = 0;
*/

    North = NULL;
    East = NULL;
    South = NULL;
    West = NULL;

    sum_avgBufUtil = 0;
    num_eject = 0;
    num_path = 0;
    T = 0;
}


void Router::Reset()
{
    memset(arr_rate,0,MAX_CHAN*sizeof(float));
    memset(mat_lambda,0,MAX_CHAN*MAX_CHAN*sizeof(float));
    memset(mat_delta,0,MAX_CHAN*MAX_CHAN*sizeof(float));
    memset(mat_Contention,0,MAX_CHAN*MAX_CHAN*sizeof(float));
    memset(mat_resServTime,0,MAX_CHAN*sizeof(float));
    memset(mat_forwProb,0,MAX_CHAN*MAX_CHAN*sizeof(float));
    memset(mat_avgBufUtil,0,MAX_CHAN*sizeof(float));
    memset(mat_avgWaitTime,0,MAX_CHAN*sizeof(float));
    memset(serv_time,0,MAX_CHAN*sizeof(float));
    memset(out_index,0,MAX_CHAN*sizeof(int));
    memset(in_index,0,MAX_CHAN*sizeof(int));

    out_index[0] = 100;
    sum_avgBufUtil = 0;
    num_eject = 0;
    num_path = 0;
    T = 0;
}
void Router::calc_serv_time()
{
    if((num_path + num_eject) != 0)
        T = (num_path * service_time + num_eject *(service_time-1))/(num_path + num_eject);
    else
        T = 0.0;

    //temporary
    for(int i =0;i<MAX_CHAN;i++)
        serv_time[i] = T;
}

void Router::calc_mat_fwd_prob()
{
    for(int i=0;i < MAX_CHAN;i++)
    {
        for(int j=0;j<MAX_CHAN;j++)
        {
            if((arr_rate[i]) && (i != j))
                mat_forwProb[i][j] = mat_lambda[i][j]/arr_rate[i];
        }
    }
}

void Router::calc_mat_contention()
{
    for(int i=0;i < MAX_CHAN;i++)
    {
        for(int j=0;j<MAX_CHAN;j++)
        {
            for(int k=0;k<MAX_CHAN;k++)
            {
                if(i == j)
                    mat_Contention[i][j] = 1;
                else
                    mat_Contention[i][j] += mat_forwProb[i][k]*mat_forwProb[j][k];
            }
        }
    }
}

void Router::calc_mat_delta()
{
    for(int i=0;i < MAX_CHAN;i++)
    {
        for(int j=0;j<MAX_CHAN;j++)
        {
            if(i == j)
                mat_delta[i][j] = arr_rate[i];
        }
    }
}

void Router::calc_mat_ResServTime()
{
    for(int i=0;i < MAX_CHAN;i++)
    {
        //mat_resServTime[i] = arr_rate[i] * service_time * service_time / 2; //R for M/G/1
       // mat_resServTime[i] = arr_rate[i] * T * T / 2; //R for M/G/1
        mat_resServTime[i] = arr_rate[i] * serv_time[i] * serv_time[i] / 2; //R for M/G/1

        //if(mat_resServTime[i] > 1)
        //    cout << this->CoOrd.Y << " " << this->CoOrd.X << endl;
    }
}

void Router::calc_mat_BuffUtil()
{
    mat delta(MAX_CHAN,MAX_CHAN);
    mat residual(MAX_CHAN,1);
    mat C(MAX_CHAN,MAX_CHAN);
    mat I = eye( MAX_CHAN, MAX_CHAN );
    mat T_mat(MAX_CHAN,MAX_CHAN);

    for (int i = 0 ; i< MAX_CHAN;i++)
    {
        for(int j =0; j<MAX_CHAN;j++)
        {
            delta(i,j) = mat_delta[i][j];
            C(i,j) = mat_Contention[i][j];

            if(i==j)
                T_mat(i,j) = serv_time[i];
            else
                T_mat(i,j) = 0;
        }
        residual(i,0) = mat_resServTime[i];
    }
    mat alpha = delta * residual;

    mat temp1 = delta * C;
    //mat temp2 = service_time * temp1;
    //mat temp2 = T * temp1;
    mat temp2 = temp1 * T_mat;

    //bool same1 = false;
   // same1 = approx_equal(temp2, tempx, "absdiff", 0);

    mat temp3 = I - temp2;

    mat beta = inv(temp3);

    mat N = beta*alpha;

   // cout << N << endl;
    for (int i = 0 ; i< MAX_CHAN;i++)
    {
        mat_avgBufUtil[i] = N(i,0);

        sum_avgBufUtil += mat_avgBufUtil[i];
       // assert(mat_avgBufUtil[i] >= 0);
        if(mat_avgBufUtil[i] < 0)
        {
            negative_flag = true;
            //cout << temp3;
        }
      //  cout << mat_avgBufUtil[i] << endl;
    }
}

void Router::calc_mat_waitTime()
{
    for (int i = 0 ; i< MAX_CHAN;i++)
    {
        if (arr_rate[i])
            mat_avgWaitTime[i] = mat_avgBufUtil[i] / arr_rate[i];
    }
}
