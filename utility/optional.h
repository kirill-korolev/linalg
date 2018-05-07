//
// Created by Kirill Korolev on 19/02/2018.
//

#ifndef LINALG_UTILITY_H
#define LINALG_UTILITY_H

#include <exception>

namespace linalg {

    template <typename T>
    class Optional {
    public:

        Optional();
        Optional(const T& value);
        Optional(const Optional& cpy);
        Optional(Optional&& rvalue);
        ~Optional();

        bool hasValue() const;
        void swap(Optional& optional);

        operator T&();
        Optional& operator=(const Optional& rhs);
        Optional&& operator=(Optional&& rvalue);
    private:
        T& getRef();
    private:
        bool _valid;
        alignas(T) char _data[sizeof(T)];
    };

}

namespace linalg {

    template <typename T>
    Optional<T>::Optional(): _valid(false) {}

    template <typename T>
    Optional<T>::Optional(const T& value): _valid(true) {
        new ((void*)_data) T(value);
    }

    template <typename T>
    Optional<T>::Optional(const Optional& cpy) {
        _valid = cpy._valid;

        if(_valid){
            new ((void*)_data) T();
        }
    }

    template <typename T>
    Optional<T>::Optional(Optional&& rvalue){
        this->swap(rvalue);
    }

    template <typename T>
    Optional<T>::~Optional(){
        if(_valid){
            getRef().~T();
        }
    }

    template <typename T>
    bool Optional<T>::hasValue() const{
        return _valid;
    }

    template <typename T>
    void Optional<T>::swap(Optional& optional){
        std::swap(_valid, optional._valid);
        std::swap(_data, optional._data);
    }

    template <typename T>
    Optional<T>& Optional<T>::operator=(const Optional& rhs){
        if(this != &rhs){
            Optional<T> temp(rhs);
            this->swap(temp);
        }
        return *this;
    }

    template <typename T>
    Optional<T>&& Optional<T>::operator=(Optional&& rvalue){
        if(this != &rvalue){
            this->swap(rvalue);
        }
        return *this;
    }

    template <typename T>
    Optional<T>::operator T&(){
        if(!_valid){
            throw std::bad_exception();
        }
        return getRef();
    }

    template <typename T>
    T& Optional<T>::getRef(){
        return *static_cast<T*>((void*)_data);
    }

}

#endif //LINALG_UTILITY_H
