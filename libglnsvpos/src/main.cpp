#include <iostream>
#include <fstream>
#include <include\libglnsvpos\glnsvpos.h>
#include <include\libglnsvpos\rungekutta.h>
#include <ctime>

using namespace std;

int main()
{
    time_t start, end;
    double del_t = 1E-01;
    double max_del = 0;
    int i_max = 0;
    int n = (int) 12*3600/del_t +1;
    double **koord_raschet = new double * [n];
    for (int i = 0; i < n; i++)
    {
        koord_raschet[i] = new double [6];
    }
    double *koord_file = new double[3];
    ofstream out;
    out.open("D:\\cpp_results.txt");
    time(&start);
    ifstream in("D:\\matlab_results.txt");
    if (!in)
    {
        cout << "ERORR: File from MATLAB not open!" << endl;
    } else {
        cout << "File from MATLAB open!" << endl;
    }
    koordinate_GLONASS(koord_raschet);
    for (int i = 0; i < n; i++)
    {
        in >> koord_file[0] >> koord_file[1] >> koord_file[2];
        string koord_str1 = to_string(koord_raschet[i][0]);
        string koord_str2 = to_string(koord_raschet[i][1]);
        string koord_str3 = to_string(koord_raschet[i][2]);
        out << koord_str1 << "\t" << koord_str2 << "\t" << koord_str3 << endl;
        for (int j = 0; j < 3; j++)
        {
            if (abs(koord_raschet[i][j] - koord_file[j]) > max_del)
            {
                max_del = abs(koord_raschet[i][j] - koord_file[j]);
                i_max = i;
            }
        }
        delete [] koord_raschet[i];
        koord_raschet[i] = nullptr;
    }
    time(&end);
    in.close();
    out.close();
    delete[] koord_raschet;
    koord_raschet = nullptr;
    delete[] koord_file;
    koord_file = nullptr;
    double seconds = difftime(end, start);
    string seconds1 = to_string(seconds*1000000/n);
    cout << "Vremya vypolnenia (priblizhenno): " << seconds1 << " ms" << endl;
    string max_del1 = to_string(max_del);
    cout << "Maksimalnaia raznitca koordinat: " << max_del1 << " m" << endl;
    string imax = to_string(i_max);
    cout << "Nomer otcheta s maksimalnoy raznitcei koordinat: " << imax << endl;
}
