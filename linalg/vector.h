//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H

#include "matrix.h"

namespace linalg {

    template <size_t N, typename T=double>
    class Vector: public Matrix<1, N, T>{
    public:
        using Matrix<1, N, T>::Matrix;

        Vector(const std::initializer_list<T>& vector);
        Vector(const initialized_matrix_t<T>& matrix)=delete;
        Vector(Matrix<1, N, T>&& rvalue): Matrix<1, N, T>(rvalue){}

        T& operator()(size_t index);
        const T& operator()(size_t index) const;

    };

}


namespace linalg {

    template <size_t N, typename T>
    Vector<N, T>::Vector(const std::initializer_list<T>& vector): Matrix<1, N, T>({vector}){}

    template <size_t N, typename T>
    T& Vector<N, T>::operator()(size_t index){
        return Matrix<1, N, T>::operator()(0, index);
    }

    template <size_t N, typename T>
    const T& Vector<N, T>::operator()(size_t index) const{
        return Matrix<1, N, T>::operator()(0, index);
    }
}

#endif //LINALG_VECTOR_H
