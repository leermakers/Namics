#include "neighborlist.h"
#include "lattice_object.h"
#include "lattice_accessor.h"

Neighborlist::Neighborlist(const Lattice_object<size_t>& mask_) noexcept : m_mask{mask_}  { }

void Neighborlist::register_config(const Neighborlist_config& configuration_) {
    m_configurations.emplace_back(configuration_);

    if (!m_configurations.back().subsystem_loop)
        m_configurations.back().subsystem_loop = std::bind(&Lattice_object<size_t>::full_system_plus_direction_neighborlist, m_mask, std::placeholders::_1);

}

std::map<Dimension,int> Neighborlist::process_configuration(const Neighborlist_config& config) {
      std::map<Dimension,int> offset {
        {Dimension::X, 0},
        {Dimension::Y, 0},
        {Dimension::Z, 0}
      };

      offset[config.dimension] = config.parallel_offset*config.direction;

      return offset;
}

//Clears all registered configurations!
void::Neighborlist::build() {
  const stl::host_vector<size_t>& temp_mask = m_mask.m_data;

      for (Neighborlist_config config : m_configurations) {

        // No fluxes will ever be calculated going from the boundary size_to the system
        std::map<Dimension,int> offset = process_configuration(config);

        config.subsystem_loop([this, config, offset, temp_mask](const size_t x, const size_t y, const size_t z) mutable {  
          if ( temp_mask[m_mask.index(x, y, z)] == 1) {
            temp_subject.push_back( m_mask.index(x,y,z) );

            if ( temp_mask[m_mask.index(x + offset[Dimension::X], y + offset[Dimension::Y], z + offset[Dimension::Z])] == 1)
              temp_neighbors.push_back( m_mask.index(x + offset[Dimension::X], y + offset[Dimension::Y], z + offset[Dimension::Z]) );
          }
        });
        
      }

      m_subject.insert(m_subject.end(), temp_subject.begin(), temp_subject.end());
      m_neighbors.insert(m_neighbors.end(), temp_neighbors.begin(), temp_neighbors.end());

      m_configurations.clear();
  
}

const stl::device_vector<size_t>& Neighborlist::get_subject() const noexcept {
  return m_subject;
}

const stl::device_vector<size_t>& Neighborlist::get_neighbors() const noexcept {
  return m_neighbors;
}