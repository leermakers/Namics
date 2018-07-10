//
// Created by Alexander Kazakov on 7/10/18.
//

#include "monolit.h"

Monolit::Monolit(const std::vector<Point> &points, int box_size) {
    this->points = points;
    this->box_size = box_size;
}

int Monolit::getX() const {
    return getPrimitive(points.back()).getX();
}

int Monolit::getY() const {
    return getPrimitive(points.back()).getY();
}

int Monolit::getZ() const {
    return getPrimitive(points.back()).getZ();
}

void Monolit::shift(int dx, int dy, int dz) {
    for (auto &&point : points) {
        point.shift(dx, dy, dz);
    }
}

const std::vector<Point> Monolit::getPoints() const {
    return points;
}

Point Monolit::getPrimitive(const Point &p) const {
    return {p.getX() % box_size,
            p.getY() % box_size,
            p.getZ() % box_size};
}


