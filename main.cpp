#include <iostream>
#include "Matrix.h"

using namespace std;

// test addition of two matrices

int main() {

    Matrix m1(3, 3);
    Matrix m2(3, 3);

    double** data1 = new double*[9];
    double** data2 = new double*[9];

    for (int i = 0; i < 9; i++) {
        data1[i] = new double[9];
        data2[i] = new double[9];
        for (int j = 0; j < 9; j++) {
            data1[i][j] = i + j;
            data2[i][j] = i + j;
        }
    }

    m1.setValues(data1);
    m2.setValues(data2);

    Matrix m3 = m1 + m2;
    m3.print();
    m3 = m1;
    m3.print();

    m3 + m2;
    m3.print();

    return 0;


}