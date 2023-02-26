#include "matrix.hpp"

struct FlatMatrix : Matrix {
    typedef std::vector<data_t> M_Type;
    using Matrix::Matrix;

private:
    FlatMatrix::M_Type mat_;

    M_Type& access() {
        return mat_;
    }

    data_t& access(size_t i, size_t j) {
        return mat_[i*m + j];
    }
    
    data_t access(size_t i, size_t j) const {
        return mat_[i*m + j];
    }

    // have to manually create the row, so can't return a reference
    Row& access(size_t i) = delete;

    Row access(size_t i) const {
        Row r;
        r.resize(m);
        for(size_t j = 0; j < m; j++) {
            r[j] = access(i, j);
        }
        return r;
    }

    void store(size_t i, size_t j, const data_t& z) {
        access(i, j) = z;
    }

    void store(size_t i, const Row& r) {
        for(size_t j = 0; j < m; j++) {
            access(i, j) = r[j];
        }
    }
};