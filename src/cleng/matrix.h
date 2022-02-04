#ifndef NAMICS_MATRIX_H
#define NAMICS_MATRIX_H

#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include "nodes/point.h"

template<class T>
class ClMatrix {
private:
    std::vector<std::vector<T> > array;
    int height;
    int width;

public:
    ClMatrix<T>(int height, int width);

    explicit ClMatrix<T>(std::vector<std::vector<T> > const &array);

    ClMatrix<T>();

    int getHeight() const;

    int getWidth() const;

    ClMatrix<T> add(const ClMatrix<T> &m) const;

    ClMatrix<T> subtract(const ClMatrix<T> &m) const;

    ClMatrix<T> multiply(const ClMatrix<T> &m) const;

    ClMatrix<T> dot(const ClMatrix<T> &m) const;

    Point dot(Point point) const;

    ClMatrix<T> negate() const;

    ClMatrix<T> transpose() const;

    ClMatrix<T> multiply(const T &value) const;

    ClMatrix<T> divide(const T &value) const;

    ClMatrix<T> applyFunction(T (*function)(T)) const;

    ClMatrix<T> subMat(int startH, int startW, int h, int w) const;

    void fill(const T &value);

    void put(int h, int w, const T &value);

    T get(int h, int w) const;

    void print(std::ostream &flux) const;

    bool operator==(const ClMatrix<T> &m);

    bool operator!=(const ClMatrix<T> &m);

    ClMatrix<T> operator+=(const ClMatrix<T> &m);

    ClMatrix<T> operator-=(const ClMatrix<T> &m);

    ClMatrix<T> operator*=(const ClMatrix<T> &m);

    ClMatrix<T> operator*=(const T &m);

    ClMatrix<T> operator/=(const T &m);

    T &operator()(int y, int x);
};

template<class T>
ClMatrix<T> operator+(const ClMatrix<T> &a, const ClMatrix<T> &b);

template<class T>
ClMatrix<T> operator-(const ClMatrix<T> &a, const ClMatrix<T> &b);

template<class T>
ClMatrix<T> operator*(const ClMatrix<T> &a, const ClMatrix<T> &b);

template<class T>
ClMatrix<T> operator*(const T &b, const ClMatrix<T> &a);

template<class T>
ClMatrix<T> operator/(const ClMatrix<T> &a, const T &b);

template<class T>
std::ostream &operator<<(std::ostream &flux, const ClMatrix<T> &m);

#endif //NAMICS_MATRIX_H


template<class T>
ClMatrix<T>::ClMatrix(int height, int width) {
    this->height = height;
    this->width = width;
    this->array = std::vector<std::vector<T> >(height, std::vector<T>(width));
}

template<class T>
ClMatrix<T>::ClMatrix(std::vector<std::vector<T> > const &array) {
    if (array.size() == 0)
        throw std::invalid_argument("Size of array must be greater than 0.");

    this->height = array.size();
    this->width = array[0].size();
    this->array = array;
}

template<class T>
ClMatrix<T>::ClMatrix() {
    height = 0;
    width = 0;
}

template<class T>
int ClMatrix<T>::getHeight() const {
    return height;
}

template<class T>
int ClMatrix<T>::getWidth() const {
    return width;
}

template<class T>
void ClMatrix<T>::fill(const T &value) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            array[i][j] = value;
        }
    }
}

template<class T>
void ClMatrix<T>::put(int h, int w, const T &value) {
    if (!(h >= 0 && h < height && w >= 0 && w < width))
        throw std::invalid_argument("Index out of bounds.");

    array[h][w] = value;
}

template<class T>
T ClMatrix<T>::get(int h, int w) const {
    if (!(h >= 0 && h < height && w >= 0 && w < width))
        throw std::invalid_argument("Index out of bounds.");

    return array[h][w];
}

template<class T>
ClMatrix<T> ClMatrix<T>::multiply(const T &value) const {
    ClMatrix result(array);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result.array[i][j] *= value;
        }
    }

    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::divide(const T &value) const {
    ClMatrix result(array);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result.array[i][j] /= value;
        }
    }

    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::add(const ClMatrix &m) const {
    if (!(height == m.height && width == m.width))
        throw std::invalid_argument("Matrix dimension must be the same.");

    ClMatrix result(height, width);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result.array[i][j] = array[i][j] + m.array[i][j];
        }
    }

    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::subtract(const ClMatrix &m) const {
    if (!(height == m.height && width == m.width))
        throw std::invalid_argument("Matrix dimension must be the same.");

    ClMatrix result(height, width);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result.array[i][j] = array[i][j] - m.array[i][j];
        }
    }
    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::multiply(const ClMatrix &m) const {
    if (!(height == m.height && width == m.width))
        throw std::invalid_argument("Matrix dimension must be the same.");

    ClMatrix result(height, width);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result.array[i][j] = array[i][j] * m.array[i][j];
        }
    }
    return result;
}

template<class T>
Point ClMatrix<T>::dot(Point point) const {
    Point result;
    for (int i = 0; i < height; i++) {
        for (int h = 0; h < width; h++) {
            result[i] += array[i][h] * point[h];
        }
    }
    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::dot(const ClMatrix &m) const {
    if (width != m.height)
        throw std::invalid_argument("Dot product not compatible.");

    T w = 0;
    int mwidth = m.width;

    ClMatrix<T> result(height, mwidth);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < mwidth; j++) {
            for (int h = 0; h < width; h++) {
                w += array[i][h] * m.array[h][j];
            }
            result.array[i][j] = w;
            w = 0;
        }
    }

    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::negate() const {
    ClMatrix<T> result(width, height);

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            result.array[i][j] = -array[i][j];
        }
    }
    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::transpose() const {
    ClMatrix<T> result(width, height);

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            result.array[i][j] = array[j][i];
        }
    }
    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::applyFunction(T (*function)(T)) const {
    ClMatrix<T> result(height, width);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result.array[i][j] = (*function)(array[i][j]);
        }
    }

    return result;
}

template<class T>
ClMatrix<T> ClMatrix<T>::subMat(int startH, int startW, int h, int w) const {
    if (!(startH >= 0 && startH + h <= height && startW >= 0 && startW + w <= width))
        throw std::invalid_argument("Index out of bounds");

    ClMatrix<T> result(h, w);
    for (int i = startH; i < startH + h; i++) {
        for (int j = startW; j < startW + w; j++) {
            result.array[i - startH][j - startW] = array[i][j];
        }
    }
    return result;
}

template<class T>
void ClMatrix<T>::print(std::ostream &flux) const {
    int maxLength[width] = {};
    std::stringstream ss;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            ss << array[i][j];
            if (maxLength[j] < (int) ss.str().size()) {
                maxLength[j] = ss.str().size();
            }
            ss.str(std::string());
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            flux << array[i][j];
            ss << array[i][j];
            for (int k = 0; k < int(maxLength[j] - ss.str().size() + 1); k++) {
                flux << " ";
            }
            ss.str(std::string());
        }
        flux << std::endl;
    }
}

template<class T>
bool ClMatrix<T>::operator==(const ClMatrix &m) {
    if (height == m.height && width == m.width) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (array[i][j] != m.array[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

template<class T>
bool ClMatrix<T>::operator!=(const ClMatrix &m) {
    return !operator==(m);
}

template<class T>
ClMatrix<T> ClMatrix<T>::operator+=(const ClMatrix &m) {
    this->array = add(m).array;
    return *this;
}

template<class T>
ClMatrix<T> ClMatrix<T>::operator-=(const ClMatrix &m) {
    this->array = subtract(m).array;
    return *this;
}

template<class T>
ClMatrix<T> ClMatrix<T>::operator*=(const ClMatrix &m) {
    this->array = multiply(m).array;
    return *this;
}

template<class T>
ClMatrix<T> ClMatrix<T>::operator*=(const T &m) {
    *this = this->multiply(m);
    return *this;
}

template<class T>
ClMatrix<T> ClMatrix<T>::operator/=(const T &m) {
    *this = this->divide(m);
    return *this;
}

template<class T>
T &ClMatrix<T>::operator()(int y, int x) {
    if (!(y >= 0 && y < height && x >= 0 && x < width))
        throw std::invalid_argument("Index out of bounds.");
    return array[y][x];
}

template<class T>
ClMatrix<T> operator+(const ClMatrix<T> &a, const ClMatrix<T> &b) {
    return a.add(b);
}

template<class T>
ClMatrix<T> operator-(const ClMatrix<T> &a, const ClMatrix<T> &b) {
    return a.subtract(b);
}

template<class T>
ClMatrix<T> operator*(const ClMatrix<T> &a, const ClMatrix<T> &b) {
    return a.multiply(b);
}

template<class T>
ClMatrix<T> operator*(const T &b, const ClMatrix<T> &a) {
    return a.multiply(b);
}

template<class T>
ClMatrix<T> operator/(const ClMatrix<T> &a, const T &b) {
    return a.divide(b);
}

template<class T>
std::ostream &operator<<(std::ostream &flux, const ClMatrix<T> &m) {
    m.print(flux);
    return flux;
}
