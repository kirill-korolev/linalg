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

    //template <int N, typename T>
    //double determinant(const Matrix<N, N, T>& matrix);

    //template <int N, typename T>
    //Optional<lufactorized_pair_t<N, T>> lufactor(const Matrix<N, N, T>& matrix);

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

    /*template <int N, typename T>
    double determinant(const Matrix<N, N, T>& matrix){
        auto luOptional = lufactor(matrix);

        if(!luOptional.hasValue()) throw std::exception();

        double detL = 1;
        double detU = 1;

        lufactorized_pair_t<N, T> lu = luOptional;

        for(int i = 0; i < N; ++i){
            detL *= lu.first(i, i);
            detU *= lu.second(i, i);
        }

        return detL * detU;
    };


    template <int N, typename T>
    Optional<lufactorized_pair_t<N, T>> lufactor(const Matrix<N, N, T>& matrix){
        Matrix<N, N, T> l;
        Matrix<N, N, T> u;

        if(matrix(0, 0) == 0) return Optional<std::pair<Matrix<N, N, T>, Matrix<N, N, T>>>();

        for(int i = 0; i < N; ++i){
            u(0, i) = matrix(0, i);
            l(i, 0) = matrix(i, 0) / u(0, 0);
        }

        double subu = 0;
        double subl = 0;

        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                subu = 0;
                subl = 0;
                for(int k = 0; k < i; ++k){
                    subu += l(i, k) * u(k, j);
                    subl += l(j, k) * u(k, i);
                }
                u(i, j) = matrix(i, j) - subu;
                auto v = (matrix(j, i) - subl) / u(i, i);
                if(std::isnan(v)) v = 0;
                l(j, i) = v;
            }
        }

        return std::make_pair(l, u);
    }*/

}

#endif //LINALG_CORE_H
