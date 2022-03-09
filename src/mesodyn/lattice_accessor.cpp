#include "lattice_accessor.h"
#include "../lattice.h"

#include <algorithm>
#include <cstdio>

std::istringstream & operator >> (std::istringstream &in,  Coordinate &coordinate_class) 
{ 
    std::string token;
    std::vector<std::string> coordinate_tokens;

    while(std::getline(in, token, ',')) {
        coordinate_tokens.push_back(token);
    }

    coordinate_class.tokens_to_map(coordinate_tokens);

    return in; 
}

std::ostream & operator << (std::ostream& out,  const Coordinate &coordinate_class) 
{ 
    out << coordinate_class[Dimension::X] << ',' << coordinate_class[Dimension::Y] << ',' << coordinate_class[Dimension::Z];
    return out;
}

std::istringstream & operator >> (std::istringstream &in,  Range &range_class) 
{ 
    std::string token;
    std::vector<std::string> range_tokens;

    while(std::getline(in, token, ';')) {
        range_tokens.push_back(token);
    }

    if (range_tokens.size() != 2) {
        cerr << "Malformed range: " << in.str() << ". Got " << range_tokens.size() << " coordinate tokens, whereas 2 tokens, delimited by ';' were expected." << endl;
        throw 1;
    }

    std::istringstream stream_low(range_tokens[0]);
    std::istringstream stream_high(range_tokens[1]);

    stream_low >> range_class.m_low;
    stream_high >> range_class.m_high;

    return in; 
}

std::ostream & operator << (std::ostream& out,  const Range &range_class) 
{ 
    out << range_class.m_low[Dimension::X] << ',' << range_class.m_low[Dimension::Y] << ',' << range_class.m_low[Dimension::Z] << ';';
    out << range_class.m_high[Dimension::X] << ',' << range_class.m_high[Dimension::Y] << ',' << range_class.m_high[Dimension::Z];
    return out;
}

void Coordinate::tokens_to_map(std::vector<std::string>& tokens) {
    switch (tokens.size()) {
        case 3:
            sscanf(tokens.at(2).c_str(), "%zu", &m_coordinate[Dimension::Z]);
        case 2:
            sscanf(tokens.at(1).c_str(), "%zu", &m_coordinate[Dimension::Y]);
        case 1:
            sscanf(tokens.at(0).c_str(), "%zu", &m_coordinate[Dimension::X]);
    }
}

Coordinate::Coordinate(std::string in) {
    istringstream stream(in);
    stream >> *this;
}

Coordinate::Coordinate(size_t x, size_t y, size_t z) {
    m_coordinate[Dimension::X] = x;
    m_coordinate[Dimension::Y] = y;
    m_coordinate[Dimension::Z] = z;
}

size_t& Coordinate::operator [] (Dimension d) {
    return m_coordinate[d];
}

size_t Coordinate::operator [] (Dimension d) const {
    try {
        return m_coordinate.at(d);
    } catch (...) {
        return std::numeric_limits<size_t>::quiet_NaN();
    }
}

size_t Coordinate::at(Dimension d) {
    return m_coordinate.at(d);
}

Range::Range(const std::string& in) {
    istringstream stream(in);
    stream >> *this;
}

Range::Range(const Coordinate& low, const Coordinate& high)
: m_low{low}, m_high{high}
{ }

Range::~Range()
{}

void Range::loop(std::function<void(size_t, size_t, size_t)> function) noexcept
{
    size_t x = m_low[Dimension::X];
    do {
        size_t y = m_low[Dimension::Y];
        do {
            size_t z = m_low[Dimension::Z];
            do {
                function(x, y, z);
                ++z;
            } while (z <= m_high[Dimension::Z]);
            ++y;
        } while (y <= m_high[Dimension::Y] );
        ++x;
    } while (x <= m_high[Dimension::X] );
}

const stl::device_vector<size_t>& Range::get_indices(const Lattice_accessor& geometry_)
{
    if (m_indices.empty()) {
        loop( [this, geometry_](size_t x, size_t y, size_t z) {
            m_indices.push_back(geometry_.index(x,y,z));
        });
    }

    return m_indices;
}

size_t Range::size()
{
    auto box = as_box();

    return box[Dimension::X] * box[Dimension::Y] * box[Dimension::Z];
}

Coordinate Range::as_box() const {
    size_t x = (m_high[Dimension::X] - m_low[Dimension::X]);
    if (x == 0)
        x = 1;

    size_t y = (m_high[Dimension::Y] - m_low[Dimension::Y]);
    if (y == 0)
        y = 1;

    size_t z = (m_high[Dimension::Z] - m_low[Dimension::Z]);
    if (z == 0)
        z = 1;

    return Coordinate(x, y, z);
}

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

Coordinate Lattice_accessor::coordinate(size_t index) {
    
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

void Lattice_accessor::skip_bounds_x_major(std::function<void(size_t, size_t, size_t)> function) noexcept {
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

void Lattice_accessor::system_plus_bounds_x_major(std::function<void(size_t, size_t, size_t)> function) noexcept {
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