#ifndef DENSITY_INITIALIZER_H
#define DENSITY_INITIALIZER_H

#include "../molecule.h"
#include "../system.h"
#include "lattice_object.h"
#include <memory>
#include <vector>
#include <assert.h>

class Molecule_density {
    private:
        // List of all monomers in this molecule
        std::vector<Real> m_monomer_fraction_of_molecule;
        // In namics: theta, sum of all densities.
        const Real m_molecule_total_mass;
        // Densities per monomer that the molecule is comprised of
        std::vector<Lattice_object<Real>> m_homogeneous_monomer_densities;
        // Total densities per monomer type.
        std::vector<Real> m_total_monomer_densities;
        // The lattice on which the molecules live
        const Lattice* m_lat;

    public:
        // Used for solutes
        Molecule_density(Molecule* molecule);
        // Used for solvents
        Molecule_density(Molecule* molecule, Real molecule_total_mass);

        void set_total_monomer_densities();
        //Returns a vector of homogeneously distributed densities per monomer for this molecule
        std::vector<Lattice_object<Real>> homogeneous(size_t system_volume);
        Real total_mass();
        vector<Real>& monomer_total_mass();
};

//WARNING: THIS ASSUMES THAT SYSTEM HAS CALLED PREPAREFORCALCULATIONS
class Homogeneous_system_initializer {
    private:
        const size_t SOLVENT;
        const std::vector<Molecule*> m_molecules;
        const size_t m_system_volume;
        std::vector<Lattice_object<Real>> m_densities;

    public:
        Homogeneous_system_initializer(System* system);
        void build_objects();
        void push_data_to_objects(std::vector<Lattice_object<Real>>& target);
};

#endif