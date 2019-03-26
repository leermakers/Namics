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
    static constexpr uint8_t SYSTEM_EDGE_OFFSET = 1;
    static constexpr uint8_t BOUNDARIES = 2;
    typedef std::map<Dimension, size_t> Coordinate;

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

    const Coordinate coordinate(const size_t index) {
        
        size_t mod = 0;
        Coordinate coordinate;

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

    inline void skip_bounds(function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x{0};
    size_t y{0};

        for ( size_t z = SYSTEM_EDGE_OFFSET ; z < MZ+SYSTEM_EDGE_OFFSET ; ++z ) {
            y = SYSTEM_EDGE_OFFSET;
            do {
                x = SYSTEM_EDGE_OFFSET;
                do {
                    function(x, y, z);
                    ++x;
                } while (x < MX+SYSTEM_EDGE_OFFSET );
                ++y;
            } while (y < MY+SYSTEM_EDGE_OFFSET );
        }
    }

    inline void full_system_plus_direction_neighborlist(function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x{0};
    size_t y{0};

        for ( size_t z = 0 ; z < MZ+SYSTEM_EDGE_OFFSET ; ++z ) {
            y = SYSTEM_EDGE_OFFSET;
            do {
                x = SYSTEM_EDGE_OFFSET;
                do {
                    function(x, y, z);
                    ++x;
                } while (x < MX+SYSTEM_EDGE_OFFSET );
                ++y;
            } while (y < MY+SYSTEM_EDGE_OFFSET );
        }
    }

    inline void full_system_minus_direction_neighborlist(function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x{0};
    size_t y{0};

        for ( size_t z = 1 ; z < MZ+SYSTEM_EDGE_OFFSET+1 ; ++z ) {
            y = SYSTEM_EDGE_OFFSET;
            do {
                x = SYSTEM_EDGE_OFFSET;
                do {
                    function(x, y, z);
                    ++x;
                } while (x < MX+SYSTEM_EDGE_OFFSET );
                ++y;
            } while (y < MY+SYSTEM_EDGE_OFFSET );
        }
    }

    inline void system_plus_bounds(function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x{0};
    size_t y{0};

        for ( size_t z = 0 ; z < MZ+BOUNDARIES ; ++z ) {
            y = 0;
            do {
                x = 0;
                do {
                    function(x, y, z);
                    ++x;
                } while (x < MX+BOUNDARIES );
                ++y;
            } while (y < MY+BOUNDARIES );
        }
    }

    inline size_t index (const size_t x, const size_t y, const size_t z) noexcept {
        return (x*jump_x + y*jump_y + z*jump_z);
    }


    inline void x0_boundary(function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = 0;
    do {
        y = 0;
        do {
            function(x,y,z);
            ++y;
            } while (y < MY+BOUNDARIES);
        ++z;
        } while (z < MZ+BOUNDARIES);   
    }

    inline void xm_boundary(function<void(size_t, size_t, size_t)> function) noexcept {
      size_t x = MX + SYSTEM_EDGE_OFFSET, y = 0, z = 0;
      do {
        y = 0;
        do {
            function(x,y,z);
          ++y;
        } while (y < MY+BOUNDARIES);
        ++z;
      } while (z < MZ+BOUNDARIES);
    }

    inline void y0_boundary(function<void(size_t, size_t, size_t)> function) noexcept {
        size_t x = 0, y = 0, z = 0;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+BOUNDARIES);
        ++z;
        } while (z < MZ+BOUNDARIES);
    }

    inline void ym_boundary(function<void(size_t, size_t, size_t)> function) noexcept {
        size_t x = 0, y = MY + SYSTEM_EDGE_OFFSET, z = 0;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+BOUNDARIES);
            ++z;
        } while (z < MZ+BOUNDARIES);
    }

    inline void z0_boundary(function<void(size_t, size_t, size_t)> function) noexcept {
        size_t x = 0, y = 0, z = 0;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+BOUNDARIES);
            ++y;
        } while (y < MY+BOUNDARIES);
    }

    inline void zm_boundary(function<void(size_t, size_t, size_t)> function) noexcept {
        size_t x = 0, y = 0, z = MZ + SYSTEM_EDGE_OFFSET;
        do {
            x = 0;
            do {
                function(x,y,z);
                ++x;
            } while (x < MX+BOUNDARIES);
            ++y;
        } while (y < MY+BOUNDARIES);
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