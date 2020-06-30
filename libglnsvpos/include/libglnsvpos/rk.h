#include<tuple>
using namespace std;
#ifndef RK_H
#define RK_H

tuple<double,double,double,double,double,double> RK4(double tn, double t0, double *y0);
//void RK4(double tn, double t0, double *y0, double *yn);

#endif // RK_H
