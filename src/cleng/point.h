//
// Created by Alexander Kazakov on 7/10/18.
//

#ifndef NAMICS_POINT_H
#define NAMICS_POINT_H


#include "node.h"

class Point : public Node {
public:
    Point(int x, int y, int z);

    int getX() const override;

    int getY() const override;

    int getZ() const override;

    void shift(int dx, int dy, int dz) override;

private:
    int x;
    int y;
    int z;

    void setX(int x);

    void setY(int y);

    void setZ(int z);


};


#endif //NAMICS_POINT_H
