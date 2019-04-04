#include "collection_procedures.h"

constexpr uint8_t TOTAL_DENSITY = 1;

Norm_densities_relative::Norm_densities_relative(vector<Molecule *> mol_, vector<shared_ptr<IComponent>> components_, size_t solvent_mo)
: Norm_densities(mol_, components_, solvent_mo), m_subject_molecule{0}, m_adjustment{0}
{
  
}

void Norm_densities_relative::adjust_theta(size_t molecule_index_, Real theta_adjustment_)
{
  assert(theta_adjustment_ < 0 && "Theta adjustment > 0 can lead to densities > 1 or < 0 and would cause namics to crash.");

  m_subject_molecule = molecule_index_;
  m_adjustment = theta_adjustment_;

  m_mol[m_subject_molecule]->theta=m_mol[m_subject_molecule]->theta*(1+m_adjustment);

  m_mol[molecule_index_]->n = m_mol[molecule_index_]->theta/m_mol[molecule_index_]->chainlength;

  fetch_theta();
}

void Norm_densities_relative::execute()
{
  stl::device_vector<Real> local_adjustment(m_system_size, 0);
  stl::transform( m_components[m_subject_molecule]->rho.begin(), m_components[m_subject_molecule]->rho.end(), local_adjustment.begin(), local_adjustment.begin(), 
  [this] DEVICE_LAMBDA(const Real& x, const Real& y) { return (x / ( 1-x )) * this->m_adjustment; });

  for (size_t i = 0 ; i < m_components.size() ; ++i)
  {
    if ( i != m_subject_molecule )
      stl::transform( m_components[i]->rho.begin(), m_components[i]->rho.end(), local_adjustment.begin(), m_components[i]->rho.begin(),
        [this] DEVICE_LAMBDA(const Real& x, const Real& y) { return x * (1 - y); }
      );
  }

  stl::for_each( m_components[m_subject_molecule]->rho.begin(), m_components[m_subject_molecule]->rho.end(),
    [this] DEVICE_LAMBDA(Real& x) mutable { x *= (1 + this->m_adjustment); }
  );
}


Norm_densities::Norm_densities(vector<Molecule *> mol_, vector<shared_ptr<IComponent>> components_, size_t solvent_mol)
    : m_components{components_}, m_mol{mol_}, m_system_size{components_[0]->rho.size()}, m_solvent{solvent_mol}
{
  assert(mol_[solvent_mol]->MolMonList.size() == 1 and "Norming solvents with mutliple monomers is not supported! Please write your own script");
  fetch_theta();
}

void Norm_densities::execute()
{
  //TODO: when norming after theta in mol changed: fetch_theta before running this.
  stl::device_vector<Real> residuals(m_system_size, 0);
  for (auto component : m_components)
  {
#ifdef PAR_MESODYN
    Norm((Real *)component->rho, (theta[component] / component->theta()), m_system_size);
#else
    Real comp_theta = component->theta();
    for (Real &value : component->rho)
      value *= (theta[component] / comp_theta);
#endif
    stl::transform(component->rho.begin(), component->rho.end(), residuals.begin(), residuals.begin(), stl::plus<Real>());
  }
  stl::transform(residuals.begin(), residuals.end(), m_components[m_solvent]->rho.begin(), m_components[m_solvent]->rho.begin(),
                 [this] DEVICE_LAMBDA(const double &x, const double &y) { return (y - (x - TOTAL_DENSITY)); });
}

void Norm_densities::fetch_theta()
{
  int j = 0;
  for (auto &molecule : m_mol)
  {
    Molecule_density this_density = Molecule_density(molecule);

    vector<Real> densities = this_density.monomer_total_mass();
    for (auto density : densities)
    {

      theta[m_components[j]] = density;

      if (density == 0)
        m_solvent = j;

      ++j;
    }
  }
}

void Norm_densities::adjust_theta(size_t segment_index, Real theta_adjustment)
{
  m_mol[segment_index]->theta=m_mol[segment_index]->theta+theta_adjustment;
  
  for (auto& all_molecules : m_mol)
    all_molecules->n = all_molecules->theta/all_molecules->chainlength;

  fetch_theta();
}

Order_parameter::Order_parameter(vector<shared_ptr<IComponent>> components_, std::map<size_t, size_t> combinations_, Real boundaryless_volume_)
    : m_boundaryless_volume{boundaryless_volume_}, m_components{components_}, m_combinations{combinations_}
{
}

void Order_parameter::execute()
{
  stl::device_vector<Real> difference(m_components[0]->rho.size());
  m_order_parameter = 0;
  for (auto &index_of : m_combinations)
  {

    stl::transform(m_components[index_of.first]->rho.begin(), m_components[index_of.first]->rho.end(), m_components[index_of.second]->rho.begin(),
                   difference.begin(),
                   [this] DEVICE_LAMBDA(const Real &a, const Real &b) mutable { return pow(a - b, 2); });

#ifdef PAR_MESODYN
    m_order_parameter = thrust::reduce(difference.begin(), difference.end(), m_order_parameter);
#else
    m_order_parameter = std::accumulate(difference.begin(), difference.end(), m_order_parameter);
#endif
  }
  m_order_parameter /= m_boundaryless_volume;
}

Real& Order_parameter::get()
{
  return m_order_parameter;
}