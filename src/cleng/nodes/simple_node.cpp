#include "simple_node.h"

SimpleNode::SimpleNode(const Point &p, int id, Point box_size) :
    system_point(p),
    box_size(box_size),
    id(id) {
}

void SimpleNode::shift(const Point& shift) {
    system_point = system_point + shift;
}


Point SimpleNode::point() const {
    return system_point % box_size;
}

void SimpleNode::pushSystemPoints(std::map<int, Point> &pointsById) const {
    pointsById[id] = system_point;
}

std::string SimpleNode::to_string() const {
    return "id: " + std::to_string(id) + " { " + std::to_string(system_point.x) + ", " + std::to_string(system_point.y) + ", " + std::to_string(system_point.z) + " }";
}
