//
// Created by Alexander Kazakov on 7/9/18.
//

#pragma once

class Node {
public:
    virtual int getX() const = 0;
    virtual int getY() const = 0;
    virtual int getZ() const = 0;
    virtual void shift(int dx, int dy, int dz) = 0;
};