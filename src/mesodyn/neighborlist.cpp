#include "neighborlist.h"
#include "lattice_object.h"
#include "lattice_accessor.h"

INeighborlist_config::INeighborlist_config(Subsystem_loop loop_) : m_subsystem_loop{loop_}, m_offset {
        {Dimension::X, 0},
        {Dimension::Y, 0},
        {Dimension::Z, 0}
      }
{ }

INeighborlist_config::INeighborlist_config() : m_offset {
        {Dimension::X, 0},
        {Dimension::Y, 0},
        {Dimension::Z, 0}
      }
{ }

Offset_map INeighborlist_config::get_offset() {
  return m_offset;
}

Neighborlist_config::Neighborlist_config(Dimension dimension_, Direction direction_, size_t offset_, Subsystem_loop loop_)
: INeighborlist_config(loop_), m_dimension{dimension_}, m_direction{direction_}, m_parallel_offset{offset_}
{ 
  m_offset = process_configuration();
}

Neighborlist_config::Neighborlist_config(Dimension dimension_, Direction direction_, size_t offset_)
: INeighborlist_config(), m_dimension{dimension_}, m_direction{direction_}, m_parallel_offset{offset_}
{ 
   m_offset = process_configuration();
}

Neighborlist_config_extended::Neighborlist_config_extended(const Offset_map offset_, Subsystem_loop loop_)
: INeighborlist_config(loop_)
{
  m_offset = offset_;
}

Neighborlist_config_extended::Neighborlist_config_extended(const Offset_map offset_)
: INeighborlist_config()
{
  m_offset = offset_;
}

Neighborlist::Neighborlist(const Lattice_object<size_t>& mask_) noexcept : m_mask{mask_}  { }

void Neighborlist::register_config(const Neighborlist_config& configuration_) {
    std::shared_ptr<INeighborlist_config> ptr = std::make_shared<Neighborlist_config>(configuration_);
    m_configurations.emplace_back(ptr);

    if (!m_configurations.back()->m_subsystem_loop)
        m_configurations.back()->m_subsystem_loop = std::bind(&Lattice_object<size_t>::full_system_plus_direction_neighborlist, m_mask, std::placeholders::_1);

}

void Neighborlist::register_config(const Neighborlist_config_extended& configuration_) {
    std::shared_ptr<INeighborlist_config> ptr = std::make_shared<Neighborlist_config_extended>(configuration_);
    m_configurations.emplace_back(ptr);

    if (!m_configurations.back()->m_subsystem_loop)
        m_configurations.back()->m_subsystem_loop = std::bind(&Lattice_object<size_t>::full_system_plus_direction_neighborlist, m_mask, std::placeholders::_1);

}

Offset_map Neighborlist_config::process_configuration() {
      m_offset[this->m_dimension] = this->m_parallel_offset*this->m_direction;

      return m_offset;
}

//Clears all registered configurations!
void Neighborlist::build() {
  const stl::host_vector<size_t>& temp_mask = m_mask.m_data;

      for (auto& config : m_configurations) {

        // No fluxes will ever be calculated going from the boundary size_to the system
        Offset_map offset = config->get_offset();

        config->m_subsystem_loop([this, offset, temp_mask](const size_t x, const size_t y, const size_t z) mutable {  
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