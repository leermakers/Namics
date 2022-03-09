#include <map>
#include "density_initializer.h"
#include "stl_typedef.h"
#include "../molecule.h"
#include "../tools.h"
#include "lattice_object.h"
#include "component.h"
#include <assert.h>

//Virtual base class  
class Lattice_object_collection_procedure {
  public:
    Lattice_object_collection_procedure() { }
    virtual ~Lattice_object_collection_procedure() { }
    Lattice_object_collection_procedure(const Lattice_object_collection_procedure&) { }

    virtual void execute() = 0;
};

class Order_parameter : public Lattice_object_collection_procedure {
  public:

    Order_parameter(vector< shared_ptr<IComponent> > components_, std::multimap<size_t, size_t> combinations_, Real boundaryless_volume_);
    virtual ~Order_parameter() { }

    void execute() override;
    Real& attach();

  private:
    Real m_boundaryless_volume;
    Real m_order_parameter;
    std::vector< shared_ptr<IComponent> > m_components;
    std::multimap<size_t, size_t> m_combinations;
};

class Treat_as_zero : public Lattice_object_collection_procedure {
  public:

    Treat_as_zero(std::vector< shared_ptr<IComponent> > components_, Real tolerance_);
    virtual ~Treat_as_zero() { }

    void execute() override;

  private:
    std::vector< shared_ptr<IComponent> > m_components;
    Real m_tolerance;
};

class Norm_densities : public Lattice_object_collection_procedure {
  public:

    Norm_densities(vector<Molecule*> mol_, vector< shared_ptr<IComponent> > components_, size_t solvent_mol);

    void execute() override;
    virtual void fetch_theta();
    virtual void adjust_theta(size_t, Real);

  protected:

    std::map<shared_ptr< IComponent >, Real> theta;
    std::vector< shared_ptr<IComponent> > m_components;
    std::vector<Molecule*> m_mol;
    size_t m_system_size;
    size_t m_solvent;
};

class Norm_densities_relative : public Norm_densities {
  public:

    Norm_densities_relative(vector<Molecule*> mol_, vector< shared_ptr<IComponent> > components_, size_t solvent_mol);

    void execute() override;
    void adjust_theta(size_t, Real) override;

  protected:

    size_t m_subject_molecule;
    Real m_adjustment;
};