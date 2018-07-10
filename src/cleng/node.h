//
// Created by Alexander Kazakov on 7/9/18.
//

#pragma once

class Node {
public:
    virtual int X() const = 0;
    virtual int Y() const = 0;
    virtual int Z() const = 0;
    virtual void shift(int dx, int dy, int dz) = 0;
};