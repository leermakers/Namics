#include "collection_procedures.h"

constexpr uint8_t TOTAL_DENSITY = 1;

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

Real Order_parameter::get()
{
  return m_order_parameter;
}