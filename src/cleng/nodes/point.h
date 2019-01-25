#pragma once

#include <tuple>
#include <string>
#include "../../namics.h"

struct Point {

public:
    int x;
    int y;
    int z;

    Point(const Point &point) {
        x = point.x;
        y = point.y;
        z = point.z;
    }

    Point(int x, int y, int z) : x(x), y(y), z(z) {}

    Point() : x(0), y(0), z(0) {}

    Point negate() const {
        return {-x, -y, -z};
    }

    Point operator+(const Point &p) const {
        return {x + p.x, y + p.y, z + p.z};
    }

    Point operator-(const Point &p) const {
        return {x - p.x, y - p.y, z - p.z};
    }

    Point operator%(const Point &box) const {
        return {x % box.x, y % box.y, z % box.z};
    }

    bool operator==(const Point &p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator!=(const Point &p) const {
        return x != p.x || y != p.y || z != p.z;
    }

    bool less_all_elements_than(const Point &p) const {
        return x < p.x and y < p.y and z < p.z;
    }

    bool more_all_elements_than(const Point &p) const {
        return x > p.x and y > p.y and z > p.z;
    }

    bool operator<(const Point &other) const {
        return std::make_tuple(x, y, z) < std::make_tuple(other.x, other.y, other.z);
    }

    bool operator>(const Point &other) const {
        return std::make_tuple(x, y, z) > std::make_tuple(other.x, other.y, other.z);
    }


    std::string to_string() const {
        return "{ " + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + " }";
    }

    Real distance(const Point &other) const {
        Real dx = pow(x - other.x, 2);
        Real dy = pow(y - other.y, 2);
        Real dz = pow(z - other.z, 2);
        return sqrt(dx + dy + dz);

    }
};
