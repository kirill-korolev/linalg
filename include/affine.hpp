//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_AFFINE_H
#define LINALG_AFFINE_H

#include "matrix.hpp"
#include "core.hpp"
#include <array>
#include <cmath>

namespace linalg {

    template <typename T=double>
    Matrix3d<T> translate(double dx, double dy);

    template <typename T=double>
    Matrix3d<T> scale(double a, double b);

    template <typename T=double>
    Matrix3d<T> rotate(double alpha);

    template <typename T=double>
    Matrix4d<T> translate(double dx, double dy, double dz);

    template <typename T=double>
    Matrix4d<T> scale(double a, double b, double c);

    template <typename T=double>
    Matrix4d<T> rotateX(double alpha);

    template <typename T=double>
    Matrix4d<T> rotateY(double alpha);

    template <typename T=double>
    Matrix4d<T> rotateZ(double alpha);
}

namespace linalg {

    template <typename T>
    Matrix3d<T> translate(double dx, double dy){
        auto matrix = identity<3, 3, T>();
        matrix(0, 2) = dx;
        matrix(1, 2) = dy;
        return matrix;
    }

    template <typename T>
    Matrix3d<T> scale(double a, double b){
        Matrix3d<T> matrix(0);
        matrix(0, 0) = a;
        matrix(1, 1) = b;
        matrix(2, 2) = 1;
        return matrix;
    }

    template <typename T>
    Matrix3d<T> rotate(double alpha){
        Matrix3d<T> matrix(0);
        matrix(0, 0) = std::cos(alpha);
        matrix(0, 1) = std::sin(alpha);
        matrix(1, 0) = -std::sin(alpha);
        matrix(1, 1) = std::cos(alpha);
        matrix(2, 2) = 1;
        return matrix;
    }

    template <typename T>
    Matrix4d<T> translate(double dx, double dy, double dz){
        auto matrix = identity<4, 4, T>();
        matrix(0, 3) = dx;
        matrix(1, 3) = dy;
        matrix(2, 3) = dz;
        return matrix;
    }

    template <typename T>
    Matrix4d<T> scale(double a, double b, double c){
        Matrix4d<T> matrix(0);
        matrix(0, 0) = a;
        matrix(1, 1) = b;
        matrix(2, 2) = c;
        matrix(3, 3) = 1;
        return matrix;
    }

    template <typename T>
    Matrix4d<T> rotateX(double alpha){
        Matrix4d<T> matrix(0);
        matrix(0, 0) = 1;
        matrix(1, 1) = std::cos(alpha);
        matrix(1, 2) = std::sin(alpha);
        matrix(2, 1) = -std::sin(alpha);
        matrix(2, 2) = std::cos(alpha);
        matrix(3, 3) = 1;
        return matrix;
    }

    template <typename T>
    Matrix4d<T> rotateY(double alpha){
        Matrix4d<T> matrix(0);
        matrix(1, 1) = 1;
        matrix(0, 0) = std::cos(alpha);
        matrix(0, 2) = std::sin(alpha);
        matrix(2, 0) = -std::sin(alpha);
        matrix(2, 2) = std::cos(alpha);
        matrix(3, 3) = 1;
        return matrix;
    }

    template <typename T>
    Matrix4d<T> rotateZ(double alpha){
        Matrix4d<T> matrix(0);
        matrix(2, 2) = 1;
        matrix(0, 0) = std::cos(alpha);
        matrix(0, 1) = std::sin(alpha);
        matrix(1, 0) = -std::sin(alpha);
        matrix(1, 1) = std::cos(alpha);
        matrix(3, 3) = 1;
        return matrix;
    }

}

#endif //LINALG_AFFINE_H
