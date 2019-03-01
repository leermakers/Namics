#include "mesodyn.h"
#include "mesodyn/density_initializer.h"
#include "mesodyn/file_reader.h"
#include "mesodyn/neighborlist.h"
#include <iterator>

/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018-2019 */
/* Most of the physics in this module is based on the work of Fraaije et al. in the 1990s  */

Mesodyn::Mesodyn(int start, vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_, string name_)
    : 
      Lattice_accessor(Lat_[0]),
      name{name_}, In{In_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sta{Sta_}, Rea{Rea_}, Sys{Sys_}, New{New_},

      //Const-correct way of initializing member variables from file, see template in header file.
      //Initialization syntax: initialize<Datatype>("option", default_value);
      D                 { initialize<Real>("diffusionconstant", 0.01) },
      dt                { initialize<Real>("delta_t", 0.1) },
      mean              { initialize<Real>("mean", 0.0) },
      stddev            { initialize<Real>("stddev", (2 * D * sqrt(dt) ) )},
      seed              { initialize<Real>("seed", -12345.6789) },
      seed_specified    { seed != -12345.6789 ? true : false },   
      timesteps         { initialize<int>("timesteps", 100) },
      timebetweensaves  { initialize<int>("timebetweensaves", 1) },
      cn_ratio          { initialize<Real>("cn_ratio", 0.5) },

      //Variables for rho initialization
      initialization_mode { INIT_HOMOGENEOUS },
      component_no        { Sys[0]->SysMolMonList.size() },

      //Output-related variables
      writes{0}
  {

  // to get the correct KSAM and volume.
  Sys[0]->PrepareForCalculations();
  boundaryless_volume = Sys[0]->boundaryless_volume;

  rho.resize(component_no*system_size);
  
    //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (size_t i = 0; i < Seg.size(); ++i) {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }
  
}

Mesodyn::~Mesodyn() {
    // We only use smart pointers here, they'll take care of deleting themselves when needed.
}

bool Mesodyn::CheckInput() {
    string empty = "";

    if ( (read_filename = initialize<std::string>("read_pro",empty)) != empty)
      initialization_mode = INIT_FROMPRO;
    else if ( (read_filename = initialize<std::string>("read_vtk",empty)) != empty) {
        if (read_filename.find(".vtk") != string::npos) {
          cerr << "Mesodyn will add the component number and extension by itself (in that order), please format the remainder of the filename in the input file as ending in ." << endl;
          exit(0);
        }
        initialization_mode = INIT_FROMVTK;
    }

  return true;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {

  //   Initialize densities
  cout << "Initializing.." << endl;
  initial_conditions();

  vector<Sanity_check*> checks;

  // Attach sanity checks
  for (size_t i = 0 ; i < component_no; ++i) {
    checks.push_back(new Check_between_zero_and_one<Real>(&components[i]->rho, i));
    checks.push_back(new Check_theta<Real>(&components[i]->rho, std::accumulate(components[i]->rho.begin(), components[i]->rho.end(), 0), i));
  }

  Check_index_unity<Real> check_rho(&components[0]->rho);
  for (size_t i = 1 ; i < component_no ; ++i)
    check_rho.register_checkable(&components[i]->rho);

  //Prepare IO
  set_filename();
  cout.precision(8);

  cout << "Mesodyn is all set, starting calculations.." << endl << endl << endl;

  // Prepare callback functions for SolveMesodyn in Newton
  function<Real*()> solver_callback = bind(&Mesodyn::solve_crank_nicolson, this);
  function<void(Real*,size_t)> loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  /**** Main MesoDyn time loop ****/
  for (int t = 1; t < timesteps+1; t++) {
    cout << "\x1b[A" << "\x1b[A" << "MESODYN: t = " << t << " / " << timesteps << endl;

    gaussian->generate(system_size);

    for (auto& all_fluxes : fluxes) all_fluxes->J.save_state();
    for (auto& component : components) component->rho.save_state();

    New[0]->SolveMesodyn(loader_callback, solver_callback);
    order_parameter->execute();

    cout << "Order parameter: " << order_parameter->get() << endl;

    sanity_check();

    if (t % timebetweensaves == 0) {
      write_output(t);
    }
  } // time loop

  std::cout << "Done." << std::endl;
  return true;
}

std::map<size_t, size_t> Mesodyn::generate_pairs(size_t N)
{
    // Returns every combination in a set of size N, e.g.
    // 0,1 - 0,2 - 0,3 - 1,2 - 1,3 - 2,3

    std::string bitmask(2, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    std::map<size_t, size_t> combinations;

    do {

      std::vector<size_t> index;

      for (size_t i = 0; i < N; ++i) // [0..N-1] integers
        if (bitmask[i]) index.push_back(i);

      combinations[ index[0] ] = index[1];

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return combinations;
}

void Mesodyn::sanity_check() {
  for (auto all_components : components)
    all_components->rho.perform_checks(); 
}

void Mesodyn::load_alpha(Real* alpha, const size_t i) {
    if (i < component_no)
     dynamic_pointer_cast<Component>( components[i] )->alpha.load_array(alpha,system_size);
}

Real* Mesodyn::solve_crank_nicolson() {
  for (auto& all_components : components) {
    all_components->rho.reinstate_previous_state();
    all_components->update_boundaries();
  }

  for (auto& all_fluxes : fluxes)
    all_fluxes->flux();

  for (auto& flux : fluxes) {
    flux->component_a->update_density(flux->J, cn_ratio, +1);
    flux->component_b->update_density(flux->J, cn_ratio, -1);
  }

  for ( size_t n = 0; n < components.size() ; ++n )
    stl::copy(components[n]->rho.begin(), components[n]->rho.end(), rho.begin()+n*system_size);

  #ifdef PAR_MESODYN
  return stl::raw_pointer_cast(rho.data());
  #else
  return rho.data();
  #endif
}

int Mesodyn::initial_conditions() {
  
  Lattice_object<int> t_mask(Lat[0]);

  #if defined(PAR_MESODYN) || ! defined(CUDA)
  stl::copy(Sys[0]->KSAM, Sys[0]->KSAM+system_size, t_mask.begin());
  #else
  TransferDataToHost(t_mask.data(), Sys[0]->KSAM, system_size);
  #endif

  Lattice_object<size_t> mask( t_mask ); 

  vector<Lattice_object<Real>> densities(component_no, Lattice_object<Real>(Lat[0]) );
  

  switch (initialization_mode) {
    case INIT_HOMOGENEOUS:
      {
        Homogeneous_system_initializer initializer(Sys[0]);
        initializer.build_objects();
        initializer.push_data_to_objects(densities);
        break;
      }
    case INIT_FROMPRO:
      {
        Readable_file file(read_filename, filetype::PRO);
        Reader pro_reader;
        pro_reader.read_objects_in(file);
        pro_reader.assert_lattice_compatible(Lat[0]);
        pro_reader.push_data_to_objects(densities);
        break;
      }
    case INIT_FROMVTK:
      {
        //Vtk files are written per component, so we need multiple files
        vector<Readable_file> vtk_files;
        for (size_t i = 0 ; i < component_no ; ++i) {
          //Because the input file doesn't specify each one, but we know how many there are.
          string filename = read_filename + to_string(i+1) + ".vtk";

          vtk_files.push_back( Readable_file( filename, filetype::VTK_STRUCTURED_GRID) );
        }
        Reader vtk_reader;
        
        for (const Readable_file& each_file : vtk_files) {
          vtk_reader.read_objects_in(each_file);
          vtk_reader.assert_lattice_compatible(Lat[0]);
          vtk_reader.push_data_to_objects(densities);
        }
        break;
      }    
  }

  Boundary::Map boundary_conditions;

  // BC0: bX0, BC1: bXm, etc.
  boundary_conditions[Dimension::X] = Boundary::Adapter[Lat[0]->BC[0]];
  boundary_conditions[Dimension::Y] = Boundary::Adapter[Lat[0]->BC[2]];
  boundary_conditions[Dimension::Z] = Boundary::Adapter[Lat[0]->BC[4]];

  Mesodyn::boundary = Boundary::Factory::Create(dimensionality, mask, boundary_conditions);

  std::map<size_t, size_t> combinations = generate_pairs(component_no);

  for (auto& density : densities) 
    components.push_back(make_shared<Component>(Lat[0], boundary, density));

  if (seed_specified == true)
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev, seed);
  else
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev);  

  for (auto& index_of : combinations) {
      fluxes.push_back(
        Flux::Factory::Create(dimensionality, Lat[0], D * dt, mask, components[index_of.first], components[index_of.second], gaussian));
    }

  norm_densities = make_unique<Norm_densities>(Mol, components, Sys[0]->solvent);
  order_parameter = make_unique<Order_parameter>(components, combinations, boundaryless_volume);

  norm_densities->execute();

  return 0; 
}

/******* Output generation *******/

void Mesodyn::set_filename() {
  filename << In[0]->output_info.getOutputPath() << "mesodyn-";
  filename << time(time_t());
}

int Mesodyn::write_output(int t) {
  //Implement below for more starts? (OutputList higher numbers?)

    //This takes a ton of time, but you apparently cannot re-use the output class without ending up with increasing duplicates per timestep in .kal files..
    Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[0], writes, timesteps / timebetweensaves));

    if (!Out[0]->CheckInput(1)) {
      cout << "input_error in output " << endl;
      return 0;
    }

    Out[0]->output_nr = writes;
    Out[0]->n_output = timesteps / timebetweensaves;

    New[0]->PushOutput();

    Out[0]->push("order_parameter", order_parameter->get());
    Out[0]->push("time",t);
    Out[0]->push("timesteps", timesteps);
    Out[0]->push("timebetweensaves", timebetweensaves);
    Out[0]->push("diffusionconstant", D);
    Out[0]->push("seed", seed);
    Out[0]->push("mean", mean);
    Out[0]->push("stddev", stddev);
    Out[0]->push("delta_t", dt);
    Out[0]->push("cn_ratio", cn_ratio);

    for (size_t i = 0 ; i < components.size() ; ++i) {
      ostringstream vtk_filename;
      vtk_filename << filename.str() << "-rho" << i << "-" << writes << ".vtk";

      stl::host_vector<Real> t_rho = components[i]->rho.m_data;
      Out[0]->vtk_structured_grid(vtk_filename.str(), t_rho.data(), i);
    }
    
    ++writes;

    Out[0]->WriteOutput(writes);
    Out.clear();

  return 0;
}

vector<string> Mesodyn::PARAMETERS;
vector<string> Mesodyn::VALUES;
vector<string> Mesodyn::KEYS
{   "read_pro",
    "read_vtk",
    "diffusionconstant",
    "delta_t",
    "mean",
    "stddev",
    "seed",
    "timesteps",
    "timebetweensaves",
    "cn_ratio"
};