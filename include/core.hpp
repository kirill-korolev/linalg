//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_CORE_H
#define LINALG_CORE_H

#include "matrix.hpp"
#include <cmath>
#include <exception>

namespace linalg {

    template <int N, typename T>
    Matrix<N, N, T> identity();

    template <int N, typename T>
    Matrix<N, N, T> diagonal(const T& value);

    template <int N, typename T>
    double determinant(const Matrix<N, N, T>& matrix);

    template <int N, typename T>
    Matrix<N, N, T> ludecomp(const Matrix<N, N, T>& matrix);

    template <typename T=double>
    using Matrix2d = Matrix<2, 2, T>;
    template <typename T=double>
    using Matrix3d = Matrix<3, 3, T>;
    template <typename T=double>
    using Matrix4d = Matrix<4, 4, T>;
}

namespace linalg {

    template <int N, typename T>
    Matrix<N, N, T> identity(){
        Matrix<N, N, T> id(0);
        for(auto i = 0; i < N; ++i){
            id(i, i) = 1;
        }
        return id;
    }

    template <int N, typename T>
    Matrix<N, N, T> diagonal(const T& value){
        Matrix<N, N, T> diag(0);
        for(auto i = 0; i < N; ++i){
            diag(i, i) = value;
        }
        return diag;
    }

    template <int N, typename T>
    double determinant(const Matrix<N, N, T>& matrix){
        auto lu = ludecomp(matrix);
        double det = 1.0;

        for(int i = 0; i < N; ++i){
            det *= lu(i, i);
        }

        return det;
    };


    template <int N, typename T>
    Matrix<N, N, T> ludecomp(const Matrix<N, N, T>& matrix){
        Matrix<N, N, T> m;

        /*if(matrix(0, 0) == 0) //TODO
            return m;*/

        for(int i = 0; i < N; ++i){
            m(0, i) = matrix(0, i);
            m(i, 0) = matrix(i, 0) / m(0, 0);
        }

        m(0, 0) = matrix(0, 0);

        T subu = 0;
        T subl = 0;

        for(int i = 1; i < N; ++i){
            for(int j = i; j < N; ++j){
                subu = 0;
                subl = 0;
                for(int k = 0; k < i; ++k){
                    subu += m(i, k) * m(k, j);
                    subl += m(j, k) * m(k, i);
                }
                m(i, j) = matrix(i, j) - subu;

                if(i != j) {
                    auto v = (matrix(j, i) - subl) / m(i, i);
                    if(std::isnan(v)) v = 0;
                    m(j, i) = v;
                }
            }
        }

        return m;
    }

}

#endif //LINALG_CORE_H
