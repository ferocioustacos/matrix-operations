#include <iostream>
#include "matrix.hpp"
#include "flat_matrix.hpp"
#include <math.h>

int main(int argc, char* argv[]) {
    FlatMatrix A = FlatMatrix(1 / sqrt(3), {
        {sqrt(2), 1},
        {1, -sqrt(2)}
    });

    FlatMatrix B = FlatMatrix(2, 2, std::complex<double>(0, 1), {1, 0});


    std::cout << FlatMatrix::Hadamard().Tensor(FlatMatrix::Identity(2));
    std::cout << A.TransposeProduct() << '\n';
    std::cout << B << '\n';
    std::cout << "C = \n" << FlatMatrix::IdentityLike(10, 9);
    std::cout << Matrix::Hadamard().Tensor(FlatMatrix::Identity(2));
}