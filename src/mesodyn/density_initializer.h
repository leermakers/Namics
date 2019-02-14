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
        const size_t m_molecule_total_mass;
        // Densities per monomer that the molecule is comprised of
        std::vector<Lattice_object<Real>> m_homogeneous_monomer_densities;
        // Total densities per monomer type.
        std::vector<Real> m_total_monomer_densities;
        // The lattice on which the molecules live
        const Lattice* m_lat;

    public:
        // Used for solutes
        Molecule_density(Molecule* molecule)
        : Molecule_density(molecule, molecule->theta)
        { }

        // Used for solvents
        Molecule_density(Molecule* molecule, Real molecule_total_mass)
        : m_monomer_fraction_of_molecule(0),
          m_molecule_total_mass{(size_t)molecule_total_mass},
          m_total_monomer_densities(0),
          m_lat(molecule->Lat[0])
        {
            assert(molecule->MolMonList.size() > 0);
            // The indices in MolMonList are accepted by molecule->fraction
            // and return the mass fraction of that particular monomer in the
            // molecule
            for(int& all_monomers : molecule->MolMonList)    
                m_monomer_fraction_of_molecule.push_back(
                    molecule->fraction(all_monomers)
                );

            set_total_monomer_densities();
        }

        void set_total_monomer_densities() {
            for (Real& fraction : m_monomer_fraction_of_molecule) {
                Real total_monomer_mass = m_molecule_total_mass * fraction;

                m_total_monomer_densities.push_back(total_monomer_mass);
            }
        }

        //Returns a vector of homogeneously distributed densities per monomer for this molecule
        std::vector<Lattice_object<Real>> homogeneous(size_t system_volume) {
            // argument system_volume: volume of the system minus the boundaries and solids
            for (Real& total_monomer_mass : m_total_monomer_densities) {

                // Masking objects is by in the mask class, which should probably be used after calling this function.
                Real density_per_site = total_monomer_mass / system_volume;

                Lattice_object<Real> this_monomer_density(m_lat, density_per_site);

                m_homogeneous_monomer_densities.push_back(
                    this_monomer_density
                );

            }
          
            return m_homogeneous_monomer_densities;
        }

        size_t total_mass() {
            return m_molecule_total_mass;
        }

        vector<Real>& monomer_total_mass() {
            return m_total_monomer_densities;
        }
};

//WARNING: THIS ASSUMES THAT SYSTEM HAS CALLED PREPAREFORCALCULATIONS
class Homogeneous_system_initializer {
    private:
        const size_t SOLVENT;
        const std::vector<Molecule*> m_molecules;
        const size_t m_system_volume;
        std::vector<Lattice_object<Real>> m_densities;

    public:
        Homogeneous_system_initializer(System* system)
        : SOLVENT{(size_t)system->solvent},
          m_molecules{system->Mol},
          m_system_volume{(size_t)system->boundaryless_volume}
        { 
            assert(system->boundaryless_volume > 0);
            assert(system->Mol.size() > 1);
            assert(system->Lat.size() > 0);
        }

        std::vector<Lattice_object<Real>> build_objects()
        {
            std::vector< Molecule_density* > initializer(m_molecules.size());
            Real total_solute_mass{0};

            //We need to know total solute mass to calculate total solvent mass (theta) by subtracting
            //Total solute mass from the total volume. Later, we initialize all densities at once so all
            //densities per monomer stay in the same order as the input file.
            for (size_t i = 0 ; i < m_molecules.size() ; ++i)
                if (i != SOLVENT) {
                    initializer[i] = new Molecule_density (
                        m_molecules[i]);
                    //Pool solute mass to find solvent mass
                    total_solute_mass += initializer[i]->total_mass();
                }

            Real total_solvent_mass = m_system_volume - total_solute_mass;

            //Initialize solvent with overloaded constructor
            initializer[SOLVENT] = new Molecule_density (
                    m_molecules[SOLVENT], total_solvent_mass);
            
            //Build rho in order of the input file.
            for (Molecule_density* all_initializers : initializer)
                for (auto all_densities : all_initializers->homogeneous(m_system_volume) )
                    m_densities.push_back(all_densities);

            return m_densities;
        }
};

#endif