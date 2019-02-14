#ifndef LATTICE_ACCESSOR_H
#define LATTICE_ACCESSOR_H

#include <unistd.h> //size_t
#include "../lattice.h"
#include <map>

enum class Dimension {
        X,
        Y,
        Z,
        ALL
};

enum Dimensionality {
        one_D = 1,
        two_D = 2,
        three_D = 3
};

/* Do **NOT** inherit this publicly, accidental upcasting will cause a boatload of trouble */

class Lattice_accessor {

  public:
    Lattice_accessor(const Lattice* Lat)
    : MX{ (size_t)Lat->MX },
      MY{ (size_t)Lat->MY },
      MZ{ (size_t)Lat->MZ },
      system_size{ (size_t)Lat->M },
      dimensionality{static_cast<Dimensionality>(Lat->gradients)},
      jump_x{ (size_t)Lat->JX },
      jump_y{ (size_t)Lat->JY },
      jump_z{ (size_t)Lat->JZ }
    { }

    const size_t MX, MY, MZ;
    const size_t system_size;
    //in lattice: gradients
    const Dimensionality dimensionality;

    const std::map<Dimension, size_t> coordinate(const size_t index) {
        std::map<Dimension, size_t> coordinate;
        size_t mod = 0;
        coordinate[Dimension::X] = index / jump_x;
        mod = index % jump_x;

        if (dimensionality > 1) {
            coordinate[Dimension::Y] = mod / jump_y;
            mod = mod % jump_y;
        }

        if (dimensionality > 2) {
            coordinate[Dimension::Z] = mod / jump_z;
        }
        return coordinate;
    }

    void skip_bounds(function<void(size_t, size_t, size_t)> function) {
    size_t x{0};
    size_t y{0};

        for ( size_t z = 0 ; z < MZ+1 ; ++z ) {
            y = 1;
            do {
                x = 1;
                do {
                    function(x, y, z);
                    ++x;
                } while (x < MX+1 );
                ++y;
            } while (y < MY+1 );
        }
    }

    inline size_t index (const size_t x, const size_t y, const size_t z) {
        return (x*jump_x + y*jump_y + z*jump_z);
    }


    void x0_boundary(function<void(size_t, size_t, size_t)> function) {
    size_t x = 0, y = 0, z = 0;
    do {
        y = 0;
        do {
            function(x,y,z);
            ++y;
            } while (y < MY+2);
        ++z;
        } while (z < MZ+2);   
    }

    inline void xm_boundary(function<void(size_t, size_t, size_t)> function) {
      size_t x = MX + 1, y = 0, z = 0;
      do {
        y = 0;
        do {
            function(x,y,z);
          ++y;
        } while (y < MY+2);
        ++z;
      } while (z < MZ+2);
    }

    void y0_boundary(function<void(size_t, size_t, size_t)> function) {
        size_t x = 0, y = 0, z = 0;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+2);
        ++z;
        } while (z < MZ+2);
    }

    void ym_boundary(function<void(size_t, size_t, size_t)> function) {
        size_t x = 0, y = MY + 1, z = 0;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+2);
            ++z;
        } while (z < MZ+2);
    }

    void z0_boundary(function<void(size_t, size_t, size_t)> function) {
        size_t x = 0, y = 0, z = 0;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+2);
            ++y;
        } while (y < MY+2);
    }

    void zm_boundary(function<void(size_t, size_t, size_t)> function) {
        size_t x = 0, y = 0, z = MZ + 1;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+2);
            ++y;
        } while (y < MY+2);
    }

  private:
      //in lattice: JX
    const size_t jump_x;
    //in lattice: JY
    const size_t jump_y;
    //in lattice: JZ
    const size_t jump_z;
    //in lattice: M

};

#endif