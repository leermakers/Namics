//
// Created by Alexander Kazakov on 7/10/18.
//

#include "point.h"

Point::Point() : Point(0, 0, 0) {

}

Point::Point(int x, int y, int z) :
    x(x),
    y(y),
    z(z)
{}

int Point::X() const {
    return x;
}

int Point::Y() const {
    return y;
}

int Point::Z() const {
    return z;
}

void Point::shift(int dx, int dy, int dz) {
    setX(X() + dx);
    setY(Y() + dy);
    setZ(Z() + dz);
}

Point Point::negate() {
    return {-X(), -Y(), -Z()};
}

void Point::setX(int x) {
    this->x = x;
}

void Point::setY(int y) {
    this->y = y;
}

void Point::setZ(int z) {
    this->z = z;
}




