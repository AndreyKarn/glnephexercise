#include <iostream>

#include <include/libglnsvpos/func.h>
#include "include/libglnsvpos/glnsvpos.h"
#include "include/libglnsvpos/rungekutta.h"
#include "include/libglnsvpos/structures.h"

using namespace std;

int main() {

    cout << "Calculation started" << endl;

    if(glnsvpos(1, 1)) { // (RK_valid h): RK_valid - allows calculation RungeKutta; h - time step [s]
        cout << "Calculation finished successful" << endl;
    } else {
        cout << "Calculation finished failed " << endl;
    }
    cout << "Press any key to continue..." << endl;

}
