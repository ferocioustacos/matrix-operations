#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <complex>

#include <cmath>
#include <random>

#include <iostream>
#include <sstream>

#include <vector>
#include <array>
#include <initializer_list>

struct Matrix {
    typedef std::complex<double> data_t;
    typedef std::vector< std::vector<data_t> > M_Type;
    typedef std::vector<data_t> Row;
    typedef std::vector<data_t> Column;

    size_t n;
    size_t m;

private:
    M_Type mat_;

    Matrix() {

    }

public:
    Matrix(size_t N, size_t M) {
        resize(N, M);
    }

    Matrix(size_t N, size_t M, double defaultValue) {
        resize(N, M, defaultValue);
    }

    Matrix(const std::initializer_list<std::initializer_list<data_t>>& mat) 
        : n(mat.size()), m( mat.begin()->size() )
    {
        resize(n, m);
        size_t i = 0, j = 0;
        for(const auto& r : mat) {
            if(r.size() != m) {
                throw std::invalid_argument("Matrix not of consistent size");
            }

            for(const auto& e : r) {
                store(i, j, e);
                j++;
            }

            i++;
            j = 0;
        }
    }

    Matrix(size_t N, size_t M, const std::initializer_list<data_t>& values) 
        : n(N), m(M)
    {
        if(values.size() != N * M) {
            throw std::invalid_argument("values param not sufficient to fill Matrix");
        }
    }

    Matrix(size_t N, size_t M, const data_t& defaultValue, const std::initializer_list<data_t>& values) 
        : n(N), m(M)
    {
        set_matrix(M_Type());
        resize_matrix(n);
        const auto& V = values.begin();
        size_t v_idx = 0;
        for(size_t i = 0; i < n; i++) {
            resize_row(i, m, defaultValue);
            for(size_t j = 0; j < m; j++) {
                if(v_idx < values.size()) {
                    store(i, j, V[v_idx]);
                    v_idx++;
                }
            }
        }
    }

    Matrix(const data_t& factor, const std::initializer_list<std::initializer_list<data_t>>& mat) 
        : n(mat.size()), m( mat.begin()->size() )
    {
        resize(n, m);
        size_t i = 0, j = 0;
        for(const auto& r : mat) {
            if(r.size() != m) {
                throw std::invalid_argument("Matrix not of consistent size");
            }

            for(const auto& e : r) {
                store(i, j, e * factor);
                j++;
            }

            i++;
            j = 0;
        }
    }


    /* Special Matrix construction */
    [[nodiscard]]
    static Matrix RandomMatrix(size_t N, size_t M) {
        std::random_device device;
        std::uniform_real_distribution<double> R(-100, 100);
        Matrix mat = Matrix(N, M);
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < M; j++) {
                mat.store(i, j, R(device));
            }
        }

        return mat;
    }

    [[nodiscard]]
    static Matrix RandomMatrix(size_t N, size_t M, double min, double max) {
        std::random_device device;
        std::uniform_real_distribution<double> R(min, max);
        
        Matrix mat = Matrix(N, M);
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < M; j++) {
                mat.store(i, j, R(device));
            }
        }

        return mat;
    }    

    static Matrix diag(size_t N, double Diag[]) {
        Matrix D = Matrix(N, N);
        for(size_t i = 0; i < N; i++) {
            D.store(i, i, Diag[i]);
        }

        return D;
    }

    template<size_t N>
    static Matrix diag(std::array<double, N> values) {
        return Matrix::diag(N, values.data());
    }

    static Matrix diag(std::vector<double> values) {
        return Matrix::diag(values.size(), values.data());
    }

    [[nodiscard]]
    static Matrix Identity(size_t N) {
        Matrix I = Matrix(N, N, 0);
        for(size_t i = 0; i < N; i++) {
            I.store(i, j, 1);
        }

        return I;
    }

    [[nodiscard]]
    static Matrix IdentityLike(size_t N, size_t M) {
        Matrix I = Matrix(N, M, 0);
        for(size_t i = 0; i < N; i++) {
            I.store(i, j, 1);
        }

        return I;
    }

    [[nodiscard]]
    static Matrix Hadamard() {
        Matrix H = Matrix(2, 2, 1);
        H.store(1, 1, -1);
        return H * (1 / sqrt(2));
    }

    [[nodiscard]]
    static Matrix PauliX() {
        Matrix X = Matrix(2, 2);
        X.at(0, 1) = 1;
        X.at(1, 0) = 1;
        return X;
    }

    [[nodiscard]]
    static Matrix PauliY() {
        Matrix Y = Matrix(2, 2);
        Y.at(0, 1) = data_t(0, -1);
        Y.at(1, 0) = data_t(0, 1);
        return Y;
    }

    [[nodiscard]]
    static Matrix PauliZ() {
        Matrix Z = Matrix::Identity(2);
        Z.at(1, 1) = -1;
        return Z;
    }

    /* Matrix access */

    [[nodiscard]]
    data_t at(size_t i, size_t j) const {
        check_bounds(i, j);
        return access(i, j);
    }

    [[nodiscard]]
    data_t& at(size_t i, size_t j) {
        check_bounds(i, j);
        return access(i, j);
    }

    [[nodiscard]]
    Row getRow(size_t i) const {
        return access(i);
    }

    void setRow(size_t i, const Row& row) {
        if(row.size() != m) {
            throw std::invalid_argument("New row size does not match Matrix row size");
        }

        store(i, row);
    }
    
    [[nodiscard]]
    Column getColumn(size_t j) const {
        Column col(n);
        for(size_t i = 0; i < n; i++) {
            col[i] = access(i, j);
        }

        return col;
    }

    [[nodiscard]]
    Matrix Transpose() const {
        Matrix T = Matrix(m, n);
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) {
                T.at(j, i) = at(i, j);
            }
        }

        return T;
    }

    [[nodiscard]]
    Matrix Conjugate() const {
        Matrix T = Matrix(m, n);
        // can reduce runtime from O(2nm) to O(nm) by doing conjugate at assignment time
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) {
                data_t z = at(i, j);
                data_t conj = data_t(z.real(), -1 * z.imag());
                T.at(j, i) = conj;
            }
        }

        return T;
    }

    Matrix ConjugateProduct() const {
        return *this * Conjugate();
    }

    void resize(size_t N, size_t M) {
        resize(N, M, 0);
    }

    void resize(size_t N, size_t M, double defaultValue) {
        n = N;
        m = M;
        if(access().size() == 0) {
            set_matrix(M_Type());
        }

        resize_matrix(n);
        for(size_t i = 0; i < n; i++) {
            resize_row(i, m);                
            for(size_t j = 0; j < m; j++) {
                store(i, j, defaultValue);
            }
        }
    }

    [[nodiscard]]
    static data_t dot(const Row& row, const Column& column) {
        data_t sum = 0;
        for(size_t i = 0; i < row.size() && i < column.size(); i++) {
            sum += row.at(i) * column.at(i);
        }

        return sum;
    }

    [[nodiscard]]
    Matrix operator*(const Matrix& other) const {
        if (m != other.n) {
            throw std::invalid_argument("Matrix Multiplication: # of Rows does not match # of Cols");
        }

        Matrix res = Matrix(m, other.n);
        for(size_t i = 0; i < res.n; i++) {
            for(size_t j = 0; j < res.m; j++) {
                res.at(i, j) = dot(getRow(i), other.getColumn(j));
            }
        }

        return res;
    }


    [[nodiscard]]
    Matrix operator*(const data_t& factor) const {
        Matrix res = *this;
        for(size_t i = 0; i < res.n; i++) {
            for(size_t j = 0; j < res.m; j++) {
                res.at(i, j) *= factor;
            }
        }

        return res;
    }

    Matrix& operator*=(const data_t& factor) {
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) {
                at(i, j) *= factor;
            }
        }

        return *this;
    }

    [[nodiscard]]
    Matrix Tensor(const Matrix& other) const {
        return Matrix::Tensor(*this, other);
    }

    [[nodiscard]]
    static Matrix Tensor(const Matrix& A, const Matrix& B) {
        Matrix result = Matrix(A.n * B.n, A.m * B.m);
        for(size_t i = 0; i < A.n; i++) {
            for(size_t j = 0; j < A.m; j++) {
                for(size_t k = 0; k < B.n; k++) {
                    for(size_t l = 0; l < B.m; l++) {
                        result.at(i * B.m + k, j * B.m + l) = A.at(i, j) * B.at(k, l);
                    }
                }                
            }
        }
        return result;
    }

    [[nodiscard]]
    bool operator==(const Matrix& other) const {
        if(n != other.n) {
            return false;
        }

        if(m != other.m) {
            return false;
        }

        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) { 
                if(at(i, j) != other.at(i, j)) {
                    return false;
                }
            }
        }

        return true;
    }

    [[nodiscard]]
    std::string asLatexMatrix() const {
        const std::string newLine = "\\\\";
        std::stringstream array;
        array << "\\begin{pmatrix} ";

        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) {
                if(j < m - 1) {
                    array << at(i, j) << " & ";
                } else {
                    array << at(i, j) << ' ';
                }
            }

            if(i < n - 1)
                array << newLine << " ";
        }

        array << "\\end{pmatrix}";
        return array.str();
    }

    static bool CheckUnitary(const Matrix& A) {
        if(A.n != A.m) {
            return false;
        }

        return (A*A) == Matrix::Identity(A.n);
    }

    [[nodiscard]]
    Matrix TransposeProduct() const { 
        return Matrix::TransposeProduct(*this);
    }

    [[nodiscard]]
    static Matrix TransposeProduct(const Matrix& A) {
        return A * A.Transpose();
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& Matrix);


private:
    constexpr void check_bounds(size_t i, size_t j) const {
        if(i >= n) {
            throw std::invalid_argument("Accessed beyond row");
        }

        if(j >= m) {
            throw std::invalid_argument("Accessed beyond column");
        }
    }

    /* Matrix access functions */
    void set_matrix(const M_Type& new_mat) {
        mat_ = new_mat;
    }

    void resize_matrix(size_t N) {
        mat_.resize(N);
    }

    void resize_row(size_t i, size_t M) {
        mat_[i].resize(M);
    }

    void resize_row(size_t i, size_t M, const data_t& defaultValue) {
        mat_[i].resize(M, defaultValue);
    }

    M_Type& access() {
        return mat_;
    }

    data_t& access(size_t i, size_t j) {
        return mat_[i][j];
    }
    
    data_t access(size_t i, size_t j) const {
        return mat_[i][j];
    }

    Row& access(size_t i) {
        return mat_[i];
    }

    Row access(size_t i) const {
        return mat_[i];
    }

    void store(size_t i, size_t j, const data_t& z) {
        mat_[i][j] = z;
    }

    void store(size_t i, const Row& r) {
        mat_[i] = r;
    }
    
};

std::ostream& operator<<(std::ostream& os, const Matrix& Matrix) {
    for(size_t i = 0; i < Matrix.n; i++) {
        for(size_t j = 0; j < Matrix.m; j++) {
            std::cout << Matrix.at(i, j) << ' ';
        }
        os << '\n';
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix::Column& col) {
    for(size_t i = 0; i < col.size(); i++) {
        std::cout << col.at(i) << ' ';
        os << '\n';
    }

    return os;
}
#endif /* MATRIX_H */