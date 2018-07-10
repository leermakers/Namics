//
// Created by Alexander Kazakov on 7/10/18.
//

#ifndef NAMICS_POINT_H
#define NAMICS_POINT_H


#include "node.h"

class Point : public Node {
public:
    Point();
    Point(int x, int y, int z);

    int X() const override;

    int Y() const override;

    int Z() const override;

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
