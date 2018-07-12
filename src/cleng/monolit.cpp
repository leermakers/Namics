//
// Created by Alexander Kazakov on 7/10/18.
//

#include "monolit.h"

using std::string;

Monolit::Monolit(const std::vector<SimpleNode> &nodes) {
    this->m_nodes = nodes;
}


Point Monolit::point() const {
    return m_nodes[0].point();
}

void Monolit::shift(const Point &shift) {
    for (auto &&node : m_nodes) {
        node.shift(shift);
    }
}

void Monolit::pushSystemPoints(std::map<int, Point> &pointsById) const {
    for (auto &&n : m_nodes) {
        n.pushSystemPoints(pointsById);
    }
}

string Monolit::to_string() const {
    string res;
    for (auto &&n : m_nodes) {
        res += n.to_string();
        res += "; ";
    }
    return res;
}


