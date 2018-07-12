//
// Created by Alexander Kazakov on 7/10/18.
//

#pragma once

#include "node.h"
#include "point.h"

class SimpleNode : public Node {
public:
    SimpleNode(const Point& p, int id, Point box_size);

    void shift(const Point& shift) override;

    std::string to_string() const override;

    Point point() const override;

    void pushSystemPoints(std::map<int, Point> &pointsById) const override;

private:
    Point system_point;
    Point box_size;
    int id;
};
