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
    return true;
    //TODO need correct implementation
//    bool res;
//    bool res1;
//    Point res_point;
//    Point res_point1;
//    cout << "cnode :" << cnode->to_string() << endl;
//    res = (point() - cnode->point()).less_all_elements_than(subBoxRange);
//    res1 = (point() - cnode->point()).more_all_elements_than(subBoxRange);
//
//    cout << "res: " << res << " res1: " << res1 << endl;
//    res_point = point() - cnode->point();
//    res_point1 = point() - cnode->point() - subBoxRange;
//    cout << "In simpleNode: point - cnode " << (res_point).to_string() << endl;
//    cout << "In simpleNode: point - cnode - subBoxRange " << (res_point1).to_string() << endl;
//    return res and res1;
}

