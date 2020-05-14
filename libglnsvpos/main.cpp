#include <iostream>

#include <include/libglnsvpos/func.h>
#include "include/libglnsvpos/glnsvpos.h"
#include "include/libglnsvpos/rungekutta.h"
#include "include/libglnsvpos/structures.h"

using namespace std;

int main() {

    cout << "Calculation started" << endl;

    glnsvpos();

    cout << "Calculation finished" << endl;
}
