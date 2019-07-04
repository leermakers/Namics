  #ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <map>
#include <cassert>
#include <functional>
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

struct Neighborlist_config {
  Dimension dimension;
  int direction;
  size_t parallel_offset;
  //optional:
  std::function<void(
      std::function<void(size_t, size_t, size_t)>
    )> subsystem_loop;
};

class Neighborlist {
  private:
    std::vector<Neighborlist_config> m_configurations;
    const Lattice_object<size_t>& m_mask;

    stl::device_vector<size_t> m_subject;
    stl::device_vector<size_t> m_neighbors;

    stl::host_vector<size_t> temp_subject;
    stl::host_vector<size_t> temp_neighbors;

    void skip_bounds(std::function<void(size_t, size_t, size_t)> function);

  public:
    Neighborlist(const Lattice_object<size_t>& mask_) noexcept;

    ~Neighborlist() { }

    void register_config(const Neighborlist_config& configuration_);
    std::map<Dimension,int> process_configuration(const Neighborlist_config& config);
    void build();
    const stl::device_vector<size_t>& get_subject() const noexcept;
    const stl::device_vector<size_t>& get_neighbors() const noexcept;
};

#endif