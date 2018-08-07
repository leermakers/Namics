#pragma once
#include "node.h"
#include "simple_node.h"
#include <vector>



class Monolit : public Node {
public:
    Monolit(const std::vector<SimpleNode> &nodes);

    std::string to_string() const override;

    Point point() const override;

    void shift(const Point &shift) override;

    void pushSystemPoints(std::map<int, Point> &pointsById) const override;

private:
    std::vector<SimpleNode> m_nodes;
};
