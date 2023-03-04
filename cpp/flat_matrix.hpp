#include "matrix.hpp"

struct FlatMatrix : Matrix {
    typedef std::vector<data_t> M_Type;
    size_t n;
    size_t m;


    FlatMatrix(size_t N, size_t M) {
        resize(N, M);
    }

    FlatMatrix(size_t N, size_t M, double defaultValue) {
        resize(N, M, defaultValue);
    }

    FlatMatrix(const std::initializer_list<std::initializer_list<data_t>>& mat) 
        : n(mat.size()), m( mat.begin()->size() )
    {
        resize(n, m);
        size_t i = 0, j = 0;
        for(const auto& r : mat) {
            if(r.size() != m) {
                throw std::invalid_argument("Matrix not of consistent size");
            }

            for(const auto& e : r) {
                at(i, j) = e;
                j++;
            }

            i++;
            j = 0;
        }
    }

    FlatMatrix(size_t N, size_t M, const std::initializer_list<data_t>& values) 
        : n(N), m(M)
    {
        if(values.size() != N * M) {
            throw std::invalid_argument("values param not sufficient to fill Matrix");
        }
    }

    FlatMatrix(size_t N, size_t M, const data_t& defaultValue, const std::initializer_list<data_t>& values) 
        : n(N), m(M)
    {
        set_matrix(M_Type());
        resize_matrix(n * m);
        const auto& V = values.begin();
        size_t v_idx = 0;
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) {
                if(v_idx < values.size()) {
                    store(i, j, V[v_idx]);
                    v_idx++;
                }
            }
        }
    }

    FlatMatrix(const data_t& factor, const std::initializer_list<std::initializer_list<data_t>>& mat) 
        : n(mat.size()), m( mat.begin()->size() )
    {
        resize(n, m);
        size_t i = 0, j = 0;
        for(const auto& r : mat) {
            if(r.size() != m) {
                throw std::invalid_argument("Matrix not of consistent size");
            }

            for(const auto& e : r) {
                at(i, j) = e * factor;
                j++;
            }

            i++;
            j = 0;
        }
    }


private:
    FlatMatrix::M_Type mat_;

    void set_matrix(const M_Type& new_mat) {
        mat_ = new_mat;
    }

    void resize_matrix(size_t N) {
        mat_.resize(N);
    }

    void resize_row(size_t i, size_t M) =  delete;
    void resize_row(size_t i, size_t M, const data_t& defaultValue) = delete;

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