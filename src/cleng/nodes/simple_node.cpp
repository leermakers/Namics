#include <utility>
#include "simple_node.h"
#include "random/random.h"

SimpleNode::SimpleNode(const Point &p, int id, const Point &box_size) :
        system_point(p),
        box_size(box_size),
        id(id) {}

void SimpleNode::shift(const Point &shift) {
    system_point = system_point + shift;
}

void SimpleNode::shift(const Matrix<Real> &matrix) {
    system_point = matrix.dot(system_point);
}

Point SimpleNode::point() const {
    Point p = system_point % box_size;
    if (p.x < 0) p.x += box_size.x;
    if (p.y < 0) p.y += box_size.y;
    if (p.z < 0) p.z += box_size.z;
    return p;
}

Point SimpleNode::get_system_point() const { return system_point; }

void SimpleNode::pushSystemPoints(std::map<int, Point> &pointsById) const { pointsById[id] = system_point; }

std::string SimpleNode::to_string() const {
    return "id: " + std::to_string(id) + " { " + std::to_string(system_point.x) + ", " +
           std::to_string(system_point.y) + ", " + std::to_string(system_point.z) + " };";
}

void SimpleNode::set_cnode(shared_ptr<SimpleNode> coupled_node) { cnode = std::move(coupled_node); }

shared_ptr<SimpleNode> SimpleNode::get_cnode() { return this->cnode; }

bool SimpleNode::inSubBoxRange(const Point &subBoxRange) const {
    Real dist = distance(cnode->get_system_point());
    Point distance = {(int) dist, (int) dist, (int) dist};

//    cout << "[Simple Node]Dist: " << dist << endl;
//    cout << "[Simple Node]distance: " << distance.to_string() << endl;

// TODO: implement operator=>
    if (distance > subBoxRange) {
        cout << "Nodes are too far away from each other: " << this->id << " and " << cnode->id << endl;
        cout << "(int) Distance: " << std::to_string((int) dist) << " between of " << system_point.to_string()
             << " and " << cnode->system_point.to_string() << endl;
        return false;
    }
    if (distance == subBoxRange) {
        cout << "Nodes are too far away from each other: " << this->id << " and " << cnode->id << endl;
        cout << "(int) Distance: " << std::to_string((int) dist) << " between of " << system_point.to_string()
             << " and " << cnode->system_point.to_string() << endl;
        return false;
    }
    return true;
}

Real SimpleNode::distance(const Point &other) const {
    Real dx = pow(this->system_point.x - other.x, 2);
    Real dy = pow(this->system_point.y - other.y, 2);
    Real dz = pow(this->system_point.z - other.z, 2);
    return sqrt(dx + dy + dz);
}

Point SimpleNode::_returnSystemPoint() const {
    return get_system_point();
}

int SimpleNode::get_ID() const {
    return this->id;
}

bool SimpleNode::_isGood() const {
    bool success = true;
    int chain_length = 49;  // TODO FIX IT!
    Point p3 = cnode->get_system_point() - system_point;
    int path_length = abs(p3.x) + abs(p3.y) + abs(p3.z);

    int path_length_even = path_length % 2;
    int chain_length_even = chain_length % 2;

//    cout << "[path_length]  " << path_length << endl;
//    cout << "[chain_length] " << chain_length << endl;

    if (path_length_even == chain_length_even) success = false;
    if (path_length >= chain_length) success = false;

    cout << "[SimpleNode _isGood] success " << success << endl;
    return success;
}

bool SimpleNode::isIdInside(const int &ID) const {
    bool success = false;
    if (this->get_ID() == ID) success = true;
    return success;
}