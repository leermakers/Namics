#pragma once
#include "node.h"
#include "point.h"
#include <memory>

using std::shared_ptr;

class SimpleNode : public Node {
public:
    SimpleNode(const Point &p, int id, Point box_size);

    void shift(const Point &shift) override;

    std::string to_string() const override;

    Point point() const override;

    void pushSystemPoints(std::map<int, Point> &pointsById) const override;

    void set_cnode(shared_ptr<SimpleNode> coupled_node);

    shared_ptr<SimpleNode> get_cnode();

    bool inSubBoxRange(const Point &subBoxRange, const Point &shift) const override;

    Real distance_with_shift(const Point &point, const Point &shift) const;

    Real distance(const Point &other) const;

    Point get_system_point() const;

private:
    Point system_point;
    Point box_size;
    int id;
    shared_ptr<SimpleNode> cnode;
};
