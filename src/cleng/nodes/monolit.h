#pragma once

#include "node.h"
#include "simple_node.h"
#include <vector>
#include <memory>

using std::shared_ptr;


class Monolit : public Node {
public:
    Monolit(const std::vector<shared_ptr<SimpleNode>> &nodes);

    std::string to_string() const override;

    Point point() const override;

    void shift(const Point &shift) override;

    void pushSystemPoints(std::map<int, Point> &pointsById) const override;

    bool inSubBoxRange(const Point &subBoxRange, const Point &shift) const override;

private:
    std::vector<shared_ptr<SimpleNode>> m_nodes;
};
