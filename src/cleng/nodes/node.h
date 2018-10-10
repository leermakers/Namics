#pragma once
#include "point.h"
#include <map>
#include <string>

class Node {
public:
    virtual Point point() const = 0;
    virtual void shift(const Point& shift) = 0;
    virtual void pushSystemPoints(std::map<int, Point> &pointsById) const = 0;
    virtual bool operator <(const Node& other) const {
        return point() < other.point();
    }
    virtual std::string to_string() const = 0;
    virtual bool inSubBoxRange(Point const &subBoxRange) const = 0;
//    virtual Real distance(Point const &point) const = 0;
    virtual ~Node() = default;
};
