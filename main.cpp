#include <iostream>
#include "matrix.hpp"
#include <math.h>

int main(int argc, char* argv[]) {
    Matrix A = Matrix(1 / sqrt(3), {
        {sqrt(2), 1},
        {1, -sqrt(2)}
    });

    Matrix B = Matrix(2, 2, std::complex<double>(0, 1), {1, 0});



    std::cout << A.TransposeProduct() << '\n';

    std::cout << B;
}