#include "monolit.h"

using std::string;

Monolit::Monolit(const std::vector<shared_ptr<SimpleNode>> &nodes) {
    this->m_nodes = nodes;
}


Point Monolit::point() const {
    return m_nodes[0]->point();
}

void Monolit::shift(const Point &shift) {
    for (auto &&node : m_nodes) {
        node->shift(shift);
    }
}

void Monolit::shift(const Matrix<Real> &matrix) {
    for (auto &&node : m_nodes) {
        node->shift(matrix);
    }
}

void Monolit::pushSystemPoints(std::map<int, Point> &pointsById) const {
    for (auto &&n : m_nodes) {
        n->pushSystemPoints(pointsById);
    }
}

string Monolit::to_string() const {
    string res;
    for (auto &&n : m_nodes) {
        res += n->to_string();
        res += " ";
    }
    return res;
}

bool Monolit::inSubBoxRange(const Point &subBoxRange) const {
    for (auto &&n : m_nodes) if (!n->inSubBoxRange(subBoxRange)) return false;
    return true;
}

Point Monolit::_returnSystemPoint() const {
    for (auto &&n : m_nodes) return n->_returnSystemPoint();
}

bool Monolit::_isGood() const {
    bool success = false;
    vector<bool> results;
    for (auto &&n : m_nodes) results.push_back(n->_isGood());
    for (auto &&result: results) if (result) success = true;
    return success;
}

bool Monolit::isIdInside(const int &ID) const {
    bool success = false;
    vector<bool> results;
    for (auto &&n : m_nodes) results.push_back(n->isIdInside(ID));
    for (auto &&result: results) if (result) success = true;
    return success;

}