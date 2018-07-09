//
// Created by Alexander Kazakov on 7/10/18.
//

#include "point.h"

Point::Point(int x, int y, int z) :
    x(x),
    y(y),
    z(z)
{}

int Point::getX() const {
    return x;
}

int Point::getY() const {
    return y;
}

int Point::getZ() const {
    return z;
}

void Point::shift(int dx, int dy, int dz) {
    setX(getX() + dx);
    setY(getY() + dy);
    setZ(getZ() + dz);
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


