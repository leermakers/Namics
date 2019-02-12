#include <utility>
#include "simple_node.h"

SimpleNode::SimpleNode(const Point &p, int id, Point box_size) :
        system_point(p),
        box_size(box_size),
        id(id) {}

void SimpleNode::shift(const Point &shift) {
    system_point = system_point + shift;
}


Point SimpleNode::point() const {
    return system_point % box_size;
}

Point SimpleNode::get_system_point() const {
    return system_point;
}

void SimpleNode::pushSystemPoints(std::map<int, Point> &pointsById) const {
    pointsById[id] = point();
}

std::string SimpleNode::to_string() const {
    return "id: " + std::to_string(id) + " { " + std::to_string(system_point.x) + ", " +
           std::to_string(system_point.y) + ", " + std::to_string(system_point.z) + " };";
}

void SimpleNode::set_cnode(shared_ptr<SimpleNode> coupled_node) {
    cnode = std::move(coupled_node);
}

shared_ptr<SimpleNode> SimpleNode::get_cnode() {
    return this->cnode;
}

bool SimpleNode::inSubBoxRange(const Point &subBoxRange, const Point &shift) const {
    Real dist = distance_with_shift(cnode->get_system_point(), shift);
    Point distance = {(int) dist + 2, (int) dist + 2, (int) dist + 2};

    if (distance > subBoxRange) {
        cout << "Too far away from each other nodes: " << this->id << " and " << cnode->id << endl;
        return false;
    }
    if (distance == subBoxRange) {
        cout << "Too far away from each other nodes: " << this->id << " and " << cnode->id << endl;
        return false;
    }
    return true;
}

Real SimpleNode::distance_with_shift(const Point &point, const Point &shift) const {
    Real dx = pow(this->system_point.x + shift.x - point.x, 2);
    Real dy = pow(this->system_point.y + shift.y - point.y, 2);
    Real dz = pow(this->system_point.z + shift.z - point.z, 2);
    return sqrt(dx + dy + dz);
}

Real SimpleNode::distance(const Point &other) const {
    Real dx = pow(this->system_point.x - other.x, 2);
    Real dy = pow(this->system_point.y - other.y, 2);
    Real dz = pow(this->system_point.z - other.z, 2);
    return sqrt(dx + dy + dz);
}
