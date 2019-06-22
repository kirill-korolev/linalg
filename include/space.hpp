//
// Created by Kirill Korolev on 11/04/2018.
//

#ifndef LINALG_SPACE_H
#define LINALG_SPACE_H

#include "vector.h"
#include <cmath>
#include <functional>

namespace linalg {

    template <int N, typename T>
    T product(const Vector<N, T>& a, const Vector<N, T>& b);

    template <int P, int N, typename T>
    T norm(const Vector<N, T>& vector);

    template <int N, typename T>
    T euclidian(const Vector<N, T>& vector);

    template <int N, typename T>
    T manhattan(const Vector<N, T>& vector);

    template <int N, typename T>
    using NormAlgorithm = T(*)(const linalg::Vector<N, T>&);

    template <int N, typename T>
    T distance(Vector<N, T>& a, const Vector<N, T>& b, NormAlgorithm<N, T> norm);
}

namespace linalg{

    template <int N, typename T>
    T product(const Vector<N, T>& a, const Vector<N, T>& b){
        T result = 0;
        for(auto i = 0; i != N; ++i){
            result += a(i) * b(i);
        }
        return result;
    }

    template <int P, int N, typename T>
    T norm(const Vector<N, T>& vector){
        static_assert(P >= 1, "P cannot be less than one");
        T sum = 0;
        for(auto i = 0; i != N; ++i){
            sum += pow(vector(i), P);
        }
        return pow(sum, (double)1 / P);
    }

    template <int N, typename T>
    T euclidian(const Vector<N, T>& vector){
        return norm<2>(vector);
    }

    template <int N, typename T>
    T manhattan(const Vector<N, T>& vector){
        return norm<1>(vector);
    }

    template <int N, typename T>
    T distance(Vector<N, T>& a, const Vector<N, T>& b, NormAlgorithm<N, T> norm){
        const Vector<N, T> sub = a - b;
        return norm(sub);
    }
}

#endif //LINALG_SPACE_H
