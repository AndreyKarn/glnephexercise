#include<integ.h>
#include<rk.h>
#include<tuple>
using namespace std;

//void RK4(double tn, double t0, double *y0, double *yn){
tuple<double,double,double,double,double,double> RK4(double tn, double t0, double *y0){
    double del_t = tn - t0;
    double *K1 = new double[6];
    double *K2 = new double[6];
    double *K3 = new double[6];
    double *K4 = new double[6];
    double *fy1 = new double[6];
    double *fy2 = new double[6];
    double *fy3 = new double[6];

// Rschet K1
    auto diff1 = integ(y0);
    K1[0] = del_t*get<0>(diff1);
    K1[1] = del_t*get<1>(diff1);
    K1[2] = del_t*get<2>(diff1);
    K1[3] = del_t*get<3>(diff1);
    K1[4] = del_t*get<4>(diff1);
    K1[5] = del_t*get<5>(diff1);
//    integ(y0, K1);
    for (int i = 0; i < 6; i++){
        fy1[i] = y0[i] + K1[i]/2.0;
    }
    // Rschet K2
    auto diff2 = integ(fy1);
    K2[0] = del_t*get<0>(diff2);
    K2[1] = del_t*get<1>(diff2);
    K2[2] = del_t*get<2>(diff2);
    K2[3] = del_t*get<3>(diff2);
    K2[4] = del_t*get<4>(diff2);
    K2[5] = del_t*get<5>(diff2);
 //   integ(fy, K2);
    for (int i = 0; i < 6; i++){
        fy2[i] = y0[i] + K2[i]/2.0;
    }
// Rschet K3
    auto diff3 = integ(fy2);
    K3[0] = del_t*get<0>(diff3);
    K3[1] = del_t*get<1>(diff3);
    K3[2] = del_t*get<2>(diff3);
    K3[3] = del_t*get<3>(diff3);
    K3[4] = del_t*get<4>(diff3);
    K3[5] = del_t*get<5>(diff3);
//    integ(fy, K3);
    for (int i = 0; i < 6; i++){
        fy3[i] = y0[i] + K3[i];
    }
// Rschet K4
    auto diff4 = integ(fy3);
    K4[0] = del_t*get<0>(diff4);
    K4[1] = del_t*get<1>(diff4);
    K4[2] = del_t*get<2>(diff4);
    K4[3] = del_t*get<3>(diff4);
    K4[4] = del_t*get<4>(diff4);
    K4[5] = del_t*get<5>(diff4);
//    integ(fy, K4);
    double X = y0[0] + 1.0/6*(K1[0] + 2.0*K2[0] + 2.0*K3[0] + K4[0]);
    double Y = y0[1] + 1.0/6*(K1[1] + 2.0*K2[1] + 2.0*K3[1] + K4[1]);
    double Z = y0[2] + 1.0/6*(K1[2] + 2.0*K2[2] + 2.0*K3[2] + K4[2]);
    double Vx = y0[3] + 1.0/6*(K1[3] + 2.0*K2[3] + 2.0*K3[3] + K4[3]);
    double Vy = y0[4] + 1.0/6*(K1[4] + 2.0*K2[4] + 2.0*K3[4] + K4[4]);
    double Vz = y0[5] + 1.0/6*(K1[5] + 2.0*K2[5] + 2.0*K3[5] + K4[5]);
    delete[] K1;
    K1 = nullptr;
    delete[] K2;
    K2 = nullptr;
    delete[] K3;
    K3 = nullptr;
    delete[] K4;
    K4 = nullptr;
    delete[] fy1;
    fy1 = nullptr;
    delete[] fy2;
    fy2 = nullptr;
    delete[] fy3;
    fy3 = nullptr;
    return make_tuple(X,Y,Z,Vx,Vy,Vz);
}

