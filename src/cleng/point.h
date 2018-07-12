//
// Created by alexander on 7/12/18.
//

#pragma once

#include <tuple>

struct Point {
public:
    int x;
    int y;
    int z;

//    Point() : x(0), y(0), z(0) {
//    }

    Point negate() const {
        return {-x, -y, -z};
    }

    Point operator +(const Point &p) const {
        return {x + p.x, y + p.y, z + p.z};
    }

    Point operator %(const Point& box) const {
        return {x % box.x, y % box.y, z % box.z};
    }

    bool operator ==(const Point& p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator <(const Point& other) const {
        return std::make_tuple(x, y, z) < std::make_tuple(other.x, other.y, other.z);
    }
};