#include <iostream>

#include <include/libglnsvpos/func.h>
#include "include/libglnsvpos/glnsvpos.h"
#include "include/libglnsvpos/rungekutta.h"
#include "include/libglnsvpos/structures.h"

using namespace std;

int main() {

    cout << "Calculation started" << endl;

    if(glnsvpos(0, 0.1)) { // (RK_valid h): RK_valid - allows calculation RungeKutta; h - time step [s]
        cout << "Calculation finished successful" << endl;
    } else {
        cout << "Calculation finished failed " << endl;
    }


    uint64_t N = 100;
    double h = 1;

    struct Y_s *Y_data;
    Y_data = new struct Y_s[N];

    Y_data[0].X = (double)rand();
    Y_data[0].Y = (double)rand();
    Y_data[0].Z = (double)rand();
    Y_data[0].VX = (double)rand();
    Y_data[0].VY = (double)rand();
    Y_data[0].VZ = (double)rand();

    RK( N, h, Y_data);

    delete []Y_data;

    cout << "Press any key to continue..." << endl;
}
