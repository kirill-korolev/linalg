//
// Created by Kirill Korolev on 11/04/2018.
//

#ifndef LINALG_EQUATIONS_H
#define LINALG_EQUATIONS_H

#include "matrix.h"
#include "vector.h"
#include <vector>

template <size_t N, size_t M, typename T>
std::ostream& operator<<(std::ostream& os, const linalg::Matrix<N, M, T>& matrix){
    for(auto i = 0; i != matrix.getRows(); ++i){
        for(auto j = 0; j != matrix.getCols(); ++j){
            os << matrix(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

namespace linalg{
    extern const int NIL_ROW;
}

namespace linalg{

    const int NIL_ROW = -1;

    template <size_t N, size_t M, typename T>
    void gauss(Matrix<N, M, T>& matrix, std::vector<T>& where){

        size_t rows = N;
        size_t cols = M;
        size_t select = 0;

        where.assign(rows, NIL_ROW);
        const double eps = 1e-10;


        for(size_t i = 0, j = 0; j != cols && i != rows; ++j){
            select = i;

            //Select max element between rows on a j-th position
            for(size_t k = i; k != rows; ++k)
                if(abs(matrix(k, j)) > abs(matrix(select, j)))
                    select = k;

            if(abs(matrix(select, j)) < eps)
                continue;

            //Move selected row up
            matrix.swapRows(select, i);

            where[j] = i;

            //Subtract a linear combination of selected row from each row
            for(auto k = 0; k != rows; ++k){
                if(k != i){
                    double c = (double)matrix(k, j) / matrix(i, j);

                    for(auto l = j; l != cols; ++l){
                        matrix(k, l) -= c * matrix(i, l);
                    }
                }
            }

            ++i;

        }

    }

    template <size_t N, size_t M, typename T>
    size_t rang(const Matrix<N, M, T>& matrix){
        Matrix<N, M, T> tmp(matrix);
        std::vector<T> vars;
        gauss(tmp, vars);

        size_t rang = 0;

        for(auto &v: vars)
            if(v != NIL_ROW) ++rang;

        return rang;
    }
}


#endif //LINALG_EQUATIONS_H
