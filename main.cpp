#include <iostream>
#include "Matrix.h"

using namespace std;

int main() {

    Matrix m1(3, 3);
    Matrix m2(3, 3);

    Matrix m3 = m1;
    m3 = m2 + m1;
    std::cout << m2;
    return 0;


}