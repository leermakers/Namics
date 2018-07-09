//
// Created by Alexander Kazakov on 7/10/18.
//

#ifndef NAMICS_MONOLIT_H
#define NAMICS_MONOLIT_H


#include "Node.h"
#include "point.h"
#include <vector>



class Monolit : public Node {
public:
    Monolit(const std::vector<Point> &points, int box_size);
    int getX() const override;

    int getY() const override;

    int getZ() const override;

    void shift(int dx, int dy, int dz) override;

    const std::vector<Point> getPoints() const;

private:
    std::vector<Point> points;
    int box_size;

    void setX(int x);

    void setY(int y);

    void setZ(int z);

    Point getPrimitive(const Point &p) const;
};


#endif //NAMICS_MONOLIT_H
