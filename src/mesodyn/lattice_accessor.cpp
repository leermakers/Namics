#include "lattice_accessor.h"
#include "../lattice.h"

Lattice_accessor::Lattice_accessor(const Lattice* Lat)
: MX{ (size_t)Lat->MX },
    MY{ (size_t)Lat->MY },
    MZ{ (size_t)Lat->MZ },
    system_size{ (size_t)Lat->M },
    dimensionality{static_cast<Dimensionality>(Lat->gradients)},
    jump_x{ (size_t)Lat->JX },
    jump_y{ (size_t)Lat->JY },
    jump_z{ (size_t)Lat->JZ }
{ }

const Lattice_accessor::Coordinate Lattice_accessor::coordinate(size_t index) {
    
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

void Lattice_accessor::skip_bounds(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t z{0};
size_t y{0};

    for ( size_t x = SYSTEM_EDGE_OFFSET ; x < MX+SYSTEM_EDGE_OFFSET ; ++x ) {
        y = SYSTEM_EDGE_OFFSET;
        do {
            z = SYSTEM_EDGE_OFFSET;
            do {
                function(x, y, z);
                ++z;
            } while (z < MZ+SYSTEM_EDGE_OFFSET );
            ++y;
        } while (y < MY+SYSTEM_EDGE_OFFSET );
    }
}

void Lattice_accessor::full_system_plus_direction_neighborlist(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t z{0};
size_t y{0};

    for ( size_t x = 0 ; x < MX+SYSTEM_EDGE_OFFSET ; ++x ) {
        y = SYSTEM_EDGE_OFFSET;
        do {
            z = SYSTEM_EDGE_OFFSET;
            do {
                function(x, y, z);
                ++z;
            } while (z < MZ+SYSTEM_EDGE_OFFSET );
            ++y;
        } while (y < MY+SYSTEM_EDGE_OFFSET );
    }
}

void Lattice_accessor::full_system_minus_direction_neighborlist(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t z{0};
size_t y{0};

    for ( size_t x = 1 ; x < MX+SYSTEM_EDGE_OFFSET+1 ; ++x ) {
        y = SYSTEM_EDGE_OFFSET;
        do {
            z = SYSTEM_EDGE_OFFSET;
            do {
                function(x, y, z);
                ++z;
            } while (z < MZ+SYSTEM_EDGE_OFFSET );
            ++y;
        } while (y < MY+SYSTEM_EDGE_OFFSET );
    }
}

void Lattice_accessor::system_plus_bounds(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t z{0};
size_t y{0};

    for ( size_t x = 0 ; x < MX+BOUNDARIES ; ++x ) {
        y = 0;
        do {
            z = 0;
            do {
                function(x, y, z);
                ++z;
            } while (z < MZ+BOUNDARIES );
            ++y;
        } while (y < MY+BOUNDARIES );
    }
}

size_t Lattice_accessor::index (const size_t x, const size_t y, const size_t z) noexcept {
    return (x*jump_x + y*jump_y + z*jump_z);
}

size_t Lattice_accessor::index (const size_t x, const size_t y, const size_t z) const noexcept {
    return (x*jump_x + y*jump_y + z*jump_z);
}


void Lattice_accessor::x0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t x = 0, y = 0, z = 0;
do {
    z = 0;
    do {
        function(x,y,z);
        ++z;
        } while (z < MZ+BOUNDARIES);
    ++y;
    } while (y < MY+BOUNDARIES);   
}

void Lattice_accessor::xm_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = MX + SYSTEM_EDGE_OFFSET, y = 0, z = 0;
    do {
    z = 0;
    do {
        function(x,y,z);
        ++z;
    } while (z < MZ+BOUNDARIES);
    ++y;
    } while (y < MY+BOUNDARIES);
}

void Lattice_accessor::y0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = 0;
    do {
        z = 0;
        do {
            function(x,y,z);
            ++z;
        } while (z < MZ+BOUNDARIES);
    ++x;
    } while (x < MX+BOUNDARIES);
}

void Lattice_accessor::ym_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = MY + SYSTEM_EDGE_OFFSET, z = 0;
    do {
        z = 0;
        do {
            function(x,y,z);
            ++z;
        } while (z < MZ+BOUNDARIES);
        ++x;
    } while (x < MX+BOUNDARIES);
}

void Lattice_accessor::z0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = 0;
    do {
        y = 0;
        do {
            function(x,y,z);
            ++y;
        } while (y < MY+BOUNDARIES);
        ++x;
    } while (x < MX+BOUNDARIES);
}

void Lattice_accessor::zm_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = MZ + SYSTEM_EDGE_OFFSET;
    do {
        y = 0;
        do {
            function(x,y,z);
            ++y;
        } while (y < MY+BOUNDARIES);
        ++x;
    } while (x < MX+BOUNDARIES);
}