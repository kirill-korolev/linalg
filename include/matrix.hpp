//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <exception>
#include <initializer_list>

namespace linalg {

    template <int N, typename T>
    class Vector;

    template <int N, int M, typename T>
    class Matrix {
    public:

        Matrix();
        explicit Matrix(const T& value);
        Matrix(const std::initializer_list<std::initializer_list<T>>& list);
        Matrix(const Matrix& other);
        Matrix(Matrix&& other) noexcept;
        Matrix& operator=(const Matrix& rhs);
        Matrix& operator=(Matrix&& rhs) noexcept;
        virtual ~Matrix();

        Matrix& operator+=(const Matrix& rhs);
        Matrix operator+(const Matrix& rhs);

        Matrix& operator-=(const Matrix& rhs);
        Matrix operator-(const Matrix& rhs);

        template <int K>
        Matrix<N, K, T> operator*(const Matrix<M, K, T>& rhs);
        Matrix& operator*=(const Matrix<M, M, T>& rhs);
        Matrix operator*(const T& rhs);

        Matrix<M, N, T> transpose();

        T& operator()(int row, int col);
        const T& operator()(int row, int col) const;

        template <typename U>
        friend Matrix operator*(const U& lhs, const Matrix& rhs){
            return rhs * lhs;
        }

        template <int K>
        Matrix<N, M + K, T> hstack(const Matrix<N, K, T>& rhs) const;

        template <int K>
        Matrix<N + K, M, T> vstack(const Matrix<K, M, T>& rhs) const;

        const int rows = N;
        const int cols = M;
    protected:
        void clear();
    protected:
        T* _data;
        int _size;
    };

}

namespace linalg {

    template <int N, int M, typename T>
    Matrix<N, M, T>::Matrix(): Matrix(T()){}

    template <int N, int M, typename T>
    Matrix<N, M, T>::Matrix(const T& value){
        _size = N * M;
        _data = static_cast<T*>(operator new[](_size * sizeof(T)));

        for(auto i = 0; i != _size; ++i){
            _data[i] = value;
        }
    }

    template <int N, int M, typename T>
    Matrix<N, M, T>::Matrix(const std::initializer_list<std::initializer_list<T>>& list){
        auto n = list.size();
        if(N != n) throw std::invalid_argument("rows dimensions differs");
        _size = N * M;
        _data = static_cast<T*>(operator new[](_size * sizeof(T)));

        int i = 0;

        for(auto row = list.begin(); row != list.end(); ++row){
            if(M != row->size()) throw std::invalid_argument("columns dimensions differs");
            for(auto col = row->begin(); col != row->end(); ++col){
                auto val = *col;
                _data[i] = val;
                ++i;
            }
        }
    }

    template <int N, int M, typename T>
    Matrix<N, M, T>::Matrix(const Matrix& other){
        _size = other._size;
        _data = static_cast<T*>(operator new[](_size * sizeof(T)));
        for(auto i = 0; i != _size; ++i){
            _data[i] = other._data[i];
        }
    }

    template <int N, int M, typename T>
    Matrix<N, M, T>::Matrix(Matrix&& other) noexcept {
        _size = other._size;
        _data = other._data;
        other._size = 0;
        other._data = nullptr;
    }

    template <int N, int M, typename T>
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

    template <int N, int M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator=(Matrix&& rhs) noexcept {
        if(this != &rhs){
            clear();
            _size = rhs._size;
            _data = rhs._data;
            rhs._size = 0;
            rhs._data = nullptr;
        }
        return *this;
    }

    template <int N, int M, typename T>
    Matrix<N, M, T>::~Matrix(){
        clear();
    }

    template <int N, int M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator+=(const Matrix& rhs){
        Matrix<N, M, T> sum = *this + rhs;
        *this = sum;
        return *this;
    }

    template <int N, int M, typename T>
    Matrix<N, M, T> Matrix<N, M, T>::operator+(const Matrix& rhs){
        Matrix<N, M, T> sum(0);
        for(auto i = 0; i < _size; ++i){
            sum._data[i] = this->_data[i] + rhs._data[i];
        }
        return sum;
    }

    template <int N, int M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator-=(const Matrix& rhs){
        Matrix<N, M, T> sub = *this - rhs;
        *this = sub;
        return *this;
    }

    template <int N, int M, typename T>
    Matrix<N, M, T> Matrix<N, M, T>::operator-(const Matrix& rhs){
        Matrix<N, M, T> sub(0);
        for(auto i = 0; i < _size; ++i){
            sub._data[i] = this->_data[i] - rhs._data[i];
        }
        return sub;
    }

    template <int N, int M, typename T>
    template <int K>
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

    template <int N, int M, typename T>
    Matrix<N, M, T>& Matrix<N, M, T>::operator*=(const Matrix<M, M, T>& rhs){
        Matrix<N, M, T> product = (*this) * rhs;
        *this = product;
        return *this;
    }

    template <int N, int M, typename T>
    Matrix<N, M, T> Matrix<N, M, T>::operator*(const T& rhs){
        Matrix<N, M, T> product;

        for(auto i = 0; i < N; ++i){
            for(auto j = 0; j < M; ++j){
                product(i, j) = (*this)(i, j) * rhs;
            }
        }

        return product;
    }

    template <int N, int M, typename T>
    Matrix<M, N, T> Matrix<N, M, T>::transpose(){
        Matrix<M, N, T> matrix;

        for(auto i = 0; i != M; ++i){
            for(auto j = 0; j != N; ++j){
                matrix(i, j) = _data[j * N + i];
            }
        }

        return matrix;
    }


    template <int N, int M, typename T>
    void Matrix<N, M, T>::clear(){
        for(auto i = 0; i != _size; ++i){
            _data[i].~T();
        }
        operator delete[](_data);
        _size = 0;
    }


    template <int N, int M, typename T>
    T& Matrix<N, M, T>::operator()(int row, int col){
        return _data[row * M + col];
    }

    template <int N, int M, typename T>
    const T& Matrix<N, M, T>::operator()(int row, int col) const {
        return _data[row * M + col];
    }

    template <int N, int M, typename T>
    template <int K>
    Matrix<N, M + K, T> Matrix<N, M, T>::hstack(const Matrix<N, K, T>& rhs) const {
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

    template <int N, int M, typename T>
    template <int K>
    Matrix<N + K, M, T> Matrix<N, M, T>::vstack(const Matrix<K, M, T>& rhs) const {
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

}

#endif //LINALG_MATRIX_H
