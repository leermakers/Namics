  #ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <map>
#include <cassert>
#include <functional>
#include <memory>
#include "stl_typedef.h"


// Forward declarations
class Lattice_accessor;
enum class Dimension;

template<typename T>
class Lattice_object;

//Configuration
enum Direction {
  minus = -1,
  plus = 1
};

typedef std::map<Dimension,int> Offset_map;

typedef std::function<void( std::function<void(size_t, size_t, size_t)> )> Subsystem_loop;

class INeighborlist_config {
  public:
    INeighborlist_config(Subsystem_loop);
    INeighborlist_config();
    virtual ~INeighborlist_config() {};

    virtual Offset_map get_offset();

    //optional:
    Subsystem_loop m_subsystem_loop;
  
  protected:
    Offset_map m_offset;
};

class Neighborlist_config : public INeighborlist_config {
  public:
    Neighborlist_config(Dimension, Direction, size_t, Subsystem_loop);
    Neighborlist_config(Dimension, Direction, size_t);

    Offset_map process_configuration();

  private:
    Dimension m_dimension;
    Direction m_direction;
    size_t m_parallel_offset;
};

class Neighborlist_config_extended : public INeighborlist_config {
  public:
    Neighborlist_config_extended(const Offset_map, Subsystem_loop);
    Neighborlist_config_extended(const Offset_map);
}; 

class Neighborlist {
  private:
    const Lattice_object<size_t>& m_mask;

    stl::device_vector<size_t> m_subject;
    stl::device_vector<size_t> m_neighbors;

    stl::host_vector<size_t> temp_subject;
    stl::host_vector<size_t> temp_neighbors;

    void skip_bounds(std::function<void(size_t, size_t, size_t)> function);

  protected:
    std::vector< std::shared_ptr<INeighborlist_config> > m_configurations;

  public:
    Neighborlist(const Lattice_object<size_t>& mask_) noexcept;

    virtual ~Neighborlist() { }

    // explicitly allow certain child classes
    void register_config(const Neighborlist_config& configuration_);
    void register_config(const Neighborlist_config_extended& configuration_);

    void build();
    const stl::device_vector<size_t>& get_subject() const noexcept;
    const stl::device_vector<size_t>& get_neighbors() const noexcept;
};

#endif