//
// Created by Alexander Kazakov on 7/10/18.
//

#include "monolit.h"

Monolit::Monolit(const std::vector<Point> &points, int box_size) {
    this->points = points;
    this->box_size = box_size;
}

int Monolit::X() const {
    return getPrimitive(points.back()).X();
}

int Monolit::Y() const {
    return getPrimitive(points.back()).Y();
}

int Monolit::Z() const {
    return getPrimitive(points.back()).Z();
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
    return {p.X() % box_size,
            p.Y() % box_size,
            p.Z() % box_size};
}


