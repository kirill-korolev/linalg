//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <exception>
#include <initializer_list>

namespace linalg {

    template <size_t N, typename T>
    class Vector;

    using shape_t = std::pair<size_t, size_t>;

    template<typename T>
    using initialized_matrix_t = std::initializer_list<std::initializer_list<T>>;

    template <size_t N, size_t M, typename T=double>
    class Matrix{
    public:

        explicit Matrix(const T& initial=0);
        Matrix(const initialized_matrix_t<T>& matrix);
        Matrix(const Matrix& cpy);
        Matrix(Matrix&& rvalue);
        Matrix& operator=(const Matrix& rhs);
        Matrix&& operator=(Matrix&& rhs);
        virtual ~Matrix();

        Matrix& operator+=(const Matrix& rhs);
        Matrix operator+(const Matrix& rhs);

        Matrix& operator-=(const Matrix& rhs);
        Matrix operator-(const Matrix& rhs);

        template <size_t K>
        Matrix<N, K, T> operator*(const Matrix<M, K, T>& rhs);
        Matrix& operator*=(const Matrix<M, M, T>& rhs);
        Matrix operator*(const T& constant);

        Matrix<M, N, T> transpose();

        size_t getRows() const;
        size_t getCols() const;
        shape_t getShape() const;

        T& operator()(size_t row, size_t col);
        const T& operator()(size_t row, size_t col) const;

        template <typename U>
        friend Matrix operator*(const U& constant, const Matrix& matrix){
            return matrix * constant;
        }

        template <size_t K>
        Matrix<N, M + K, T> hstack(const Matrix<N, K, T>& rhs) const;

        template <size_t K>
        Matrix<N + K, M, T> vstack(const Matrix<K, M, T>& rhs) const;

        void swapRows(size_t i, size_t j);
    protected:
        void assert(size_t row, size_t col) const;
        void clear();
    protected:
        T* _data;
        size_t _size;
    };

}

namespace linalg {

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>::Matrix(const T& initial){
        _size = N * M;
        _data = static_cast<T*>(operator new[](_size * sizeof(T)));

        for(auto i = 0; i != _size; ++i){
            _data[i] = initial;
        }
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>::Matrix(const initialized_matrix_t<T>& matrix){
        auto n = matrix.size();
        if(N != n) throw std::invalid_argument("row dimensions differs");
        _size = N * M;
        _data = static_cast<T*>(operator new[](_size * sizeof(T)));

        int i = 0;

        for(auto row = matrix.begin(); row != matrix.end(); ++row){
            if(M != row->size()) throw std::invalid_argument("column dimensions differs");
            for(auto col = row->begin(); col != row->end(); ++col){
                auto val = *col;
                _data[i] = val;
                ++i;
            }
        }
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>::Matrix(const Matrix& cpy){
        _size = cpy._size;
        _data = static_cast<T*>(operator new[](_size * sizeof(T)));
        for(auto i = 0; i != _size; ++i){
            _data[i] = cpy._data[i];
        }
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>::Matrix(Matrix&& rvalue){
        _size = rvalue._size;
        _data = rvalue._data;
        rvalue._size = 0;
        rvalue._data = nullptr;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator=(const Matrix& rhs){
        if(this != &rhs){
            clear();
            _size = rhs._size;
            _data = static_cast<T*>(operator new[](_size * sizeof(T)));
            for(auto i = 0; i != _size; ++i){
                _data[i] = rhs._data[i];
            }
        }
        return *this;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>&& Matrix<N, M, T>::operator=(Matrix&& rhs){
        if(this != &rhs){
            clear();
            _size = rhs._size;
            _data = rhs._data;
            rhs._size = 0;
            rhs._data = nullptr;
        }
        return *this;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>::~Matrix(){
        clear();
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator+=(const Matrix& rhs){
        Matrix<N, M, T> sum = *this + rhs;
        *this = sum;
        return *this;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T> Matrix<N, M, T>::operator+(const Matrix& rhs){
        Matrix<N, M, T> sum(0);
        for(auto i = 0; i < _size; ++i){
            sum._data[i] = this->_data[i] + rhs._data[i];
        }
        return sum;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator-=(const Matrix& rhs){
        Matrix<N, M, T> sub = *this - rhs;
        *this = sub;
        return *this;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T> Matrix<N, M, T>::operator-(const Matrix& rhs){
        Matrix<N, M, T> sub(0);
        for(auto i = 0; i < _size; ++i){
            sub._data[i] = this->_data[i] - rhs._data[i];
        }
        return sub;
    }

    template <size_t N, size_t M, typename T>
    template <size_t K>
    Matrix<N, K, T> Matrix<N, M, T>::operator*(const Matrix<M, K, T>& rhs){
        Matrix<N, K, T> product(0);

        for(auto i = 0; i < N; ++i){
            for(auto j = 0; j < M; ++j){
                for(auto k = 0; k < M; ++k){
                    product(i, j) += (*this)(i, k) * rhs(k, j);
                }
            }
        }

        return product;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator*=(const Matrix<M, M, T>& rhs){
        Matrix<N, M, T> product = (*this) * rhs;
        *this = product;
        return *this;
    }

    template <size_t N, size_t M, typename T>
    Matrix<N, M, T> Matrix<N, M, T>::operator*(const T& constant){
        Matrix<N, M, T> product;

        for(auto i = 0; i < N; ++i){
            for(auto j = 0; j < M; ++j){
                product(i, j) = (*this)(i, j) * constant;
            }
        }

        return product;
    }

    template <size_t N, size_t M, typename T>
    Matrix<M, N, T> Matrix<N, M, T>::transpose(){
        Matrix<M, N, T> matrix;

        for(auto i = 0; i != M; ++i){
            for(auto j = 0; j != N; ++j){
                matrix(i, j) = _data[j * N + i];
            }
        }

        return matrix;
    }

    template <size_t N, size_t M, typename T>
    size_t Matrix<N, M, T>::getRows() const{
        return N;
    }

    template <size_t N, size_t M, typename T>
    size_t Matrix<N, M, T>::getCols() const{
        return M;
    }

    template <size_t N, size_t M, typename T>
    shape_t Matrix<N, M, T>::getShape() const{
        return std::make_pair(N, M);
    }

    template <size_t N, size_t M, typename T>
    void Matrix<N, M, T>::clear(){
        for(auto i = 0; i != _size; ++i){
            _data[i].~T();
        }
        operator delete[](_data);
        _size = 0;
    }

    template <size_t N, size_t M, typename T>
    void Matrix<N, M, T>::assert(size_t row, size_t col) const{
        if(row >= N) throw std::invalid_argument("row value is too big");
        if(col >= M) throw std::invalid_argument("column value is too big");
    }

    template <size_t N, size_t M, typename T>
    T& Matrix<N, M, T>::operator()(size_t row, size_t col){
        assert(row, col);
        return _data[row * M + col];
    }

    template <size_t N, size_t M, typename T>
    const T& Matrix<N, M, T>::operator()(size_t row, size_t col) const{
        assert(row, col);
        return _data[row * M + col];
    }

    template <size_t N, size_t M, typename T>
    template <size_t K>
    Matrix<N, M + K, T> Matrix<N, M, T>::hstack(const Matrix<N, K, T>& rhs) const{
        Matrix<N, M + K, T> result;

        for(auto i = 0; i != N; ++i){
            for(auto j = 0; j != M; ++j){
                result(i, j) = this->operator()(i, j);
            }
            for(auto j = 0; j != K; ++j){
                result(i, M + j) = rhs(i, j);
            }
        }

        return result;
    }

    template <size_t N, size_t M, typename T>
    template <size_t K>
    Matrix<N + K, M, T> Matrix<N, M, T>::vstack(const Matrix<K, M, T>& rhs) const{
        Matrix<N + K, M, T> result;

        for(auto j = 0; j != M; ++j){
            for(auto i = 0; i != M; ++i){
                result(i, j) = this->operator()(i, j);
            }
            for(auto i = 0; i != K; ++i){
                result(N + i, j) = rhs(i, j);
            }
        }

        return result;
    }

    template <size_t N, size_t M, typename T>
    void Matrix<N, M, T>::swapRows(size_t i, size_t j){
        if(i >= N || j >= N) return;

        for(auto k = 0; k != M; ++k)
            std::swap(this->operator()(i, k), this->operator()(j, k));
    }
}

#endif //LINALG_MATRIX_H
