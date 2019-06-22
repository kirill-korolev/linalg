//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H

#include "matrix.hpp"

namespace linalg {

    template <int N, typename T>
    class Vector: public Matrix<1, N, T> {
    public:
        using Matrix<1, N, T>::Matrix;

        Vector(const std::initializer_list<T>& vector);
        Vector(const std::initializer_list<std::initializer_list<T>>& list)=delete;

        T& operator()(int index);
        const T& operator()(int index) const;

    };

}


namespace linalg {

    template <int N, typename T>
    Vector<N, T>::Vector(const std::initializer_list<T>& vector): Matrix<1, N, T>({vector}){}

    template <int N, typename T>
    T& Vector<N, T>::operator()(int index){
        return Matrix<1, N, T>::operator()(0, index);
    }

    template <int N, typename T>
    const T& Vector<N, T>::operator()(int index) const{
        return Matrix<1, N, T>::operator()(0, index);
    }
}

#endif //LINALG_VECTOR_H
