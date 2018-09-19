#include <utility>

#include <utility>
#include "simple_node.h"

SimpleNode::SimpleNode(const Point &p, int id, Point box_size) :
    system_point(p),
    box_size(box_size),
    id(id)
    {}

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

void SimpleNode::set_cnode(shared_ptr<SimpleNode> coupled_node) {
    cnode = std::move(coupled_node);
}

bool SimpleNode::inSubBoxRange(Point const &subBoxRange) const {

//    cout << "Point: " << this->to_string() << " cnode: " << cnode->to_string() << endl;
//    Real dist = system_point.distance(cnode->system_point);

    cout << "Point: " << this->point().to_string() << " cnode: " << cnode->point().to_string() << endl;
    Real dist = point().distance(cnode->point());

    Point distance_origin = {(int)dist, (int)dist, (int)dist};
    Point distance = {(int)dist+2, (int)dist+2, (int)dist+2};
    cout << "Point distance_origin: " << distance_origin.to_string() << endl;
    cout << "Point distance: " << distance.to_string() << endl;
    cout << "subbox_range: " << subBoxRange.to_string() << endl;
    if ( distance > subBoxRange ) {
        cout << "Too far away from each other nodes: " << this->id << " and " << cnode->id << endl;
        return false;
    }
    if ( distance == subBoxRange) {
        cout << "Too far away from each other nodes: " << this->id << " and " << cnode->id << endl;
        return false;
    }
    return true;
}