#include "mesodyn.h"
#include "mesodyn/density_initializer.h"
#include "mesodyn/file_reader.h"
#include "mesodyn/neighborlist.h"
#include <iterator>

//TODO: noise: visitor pattern

/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018-2019 */
/* Most of the physics in this module is based on the work of Fraaije et al. in the 1990s  */

Mesodyn::Mesodyn(int start, vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_, string name_)
    : 
      Lattice_accessor(Lat_[0]),

      //Namics classes
      name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sta{Sta_},
      Rea{Rea_},
      Sys{Sys_},
      New{New_},

      KEYS{"read_pro",
           "read_vtk",
           "diffusionconstant",
           "delta_t",
           "mean",
           "stddev",
           "seed",
           "timesteps",
           "timebetweensaves",
           "cn_ratio"},

      input_success { In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES) },

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

      //Size helper classes.
      component(0),
      flux(0),

      //Output-related variables
      writes{0},
      //write_vtk{false}, //defunct
      order_parameter{0}
  {

  set_update_lists(component_no);

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

bool Mesodyn::CheckInput(const int start) {

  input_success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);

  if (input_success) {

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
  } 

  return input_success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {

  //   Initialize densities
  cout << "Initializing.." << endl;
  initial_conditions();

  vector<Sanity_check*> checks;

  // Attach sanity checks
 // for (size_t i = 0 ; i < component_no; ++i) {
 //   checks.push_back(new Check_between_zero_and_one<Real>(&solver_component[i]->rho, i));
 //   checks.push_back(new Check_theta<Real>(&solver_component[i]->rho, solver_component[i]->theta(), i));
  //}

 // Check_index_unity<Real> check_rho(&solver_component[0]->rho);
 // for (size_t i = 1 ; i < component_no ; ++i)
 //   check_rho.register_checkable(&solver_component[i]->rho);

  //Prepare IO
  set_filename();

  cout << "Mesodyn is all set, starting calculations.." << endl;

  // Prepare callback functions for SolveMesodyn in Newton
  function<Real*()> solver_callback = bind(&Mesodyn::solve_crank_nicolson, this);
  function<void(Real*,size_t)> loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  /**** Main MesoDyn time loop ****/
  for (int t = 0; t < timesteps; t++) {
    cout << "MESODYN: t = " << t << endl;

    gaussian->generate(system_size);

    sanity_check();

    New[0]->SolveMesodyn(loader_callback, solver_callback);

  
    //Calculate and add noise flux
    for (size_t i = 0 ; i < flux.size() ; ++i)
      flux[i]->J = solver_flux[i]->J;

    cout << "Order parameter: " << calculate_order_parameter() << endl;

    for (size_t i = 0 ; i < component.size() ; ++i)
      component[i]->rho = solver_component[i]->rho;

    if (t % timebetweensaves == 0) {
      write_output(t);
    }
  } // time loop

  std::cout << "Done." << std::endl;
  return true;
}

void Mesodyn::set_update_lists(size_t N)
{

    // Returns a flux index number coupled to every combination of components, e.g.
    // flux 0, component 0 and 1. flux 1, component 0 and 2, flux 2, component 1 and 2.
    // Components are in a pair, which is coupled to a flux index in a map.
    // Access as follows: combinations[flux].first (first component) combinations[flux].second (second component)
    // Or range based loop: for (auto& combination : combinations) combination.first <<<flux>>> combination.second <<<component as above>>>

    std::string bitmask(2, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
    size_t flux{0};

    do {

      std::vector<size_t> this_combination;

      for (size_t i = 0; i < N; ++i) // [0..N-1] integers
        if (bitmask[i]) this_combination.push_back(i);

      combinations[flux] = make_pair(this_combination[0], this_combination[1]);
      ++flux;

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

}

Real Mesodyn::calculate_order_parameter() {
  stl::device_vector<Real> difference(system_size);
  Mesodyn::order_parameter = 0;
  for (auto& map : combinations) {
    auto& index_of = map.second;

    stl::transform(solver_component[index_of.first]->rho.begin(), solver_component[index_of.first]->rho.end(), solver_component[index_of.second]->rho.begin(),
      difference.begin(),
      [this] DEVICE_LAMBDA (const double& x, const double& y) mutable {return pow(x-y,2);}
      );

    #ifdef PAR_MESODYN
    Mesodyn::order_parameter = thrust::reduce(difference.begin(), difference.end(), Mesodyn::order_parameter);
    #else
    Mesodyn::order_parameter = std::accumulate(difference.begin(), difference.end(), Mesodyn::order_parameter);
    #endif
  }
  Mesodyn::order_parameter /= boundaryless_volume;

  return Mesodyn::order_parameter;
}

void Mesodyn::sanity_check() {
  for (auto all_components : solver_component)
    all_components->rho.perform_checks(); 
}

void Mesodyn::load_alpha(Real* alpha, const size_t i) {
    if (i < component_no) {
      solver_component[i]->alpha.load_array(alpha,system_size);
    }
}

Real* Mesodyn::solve_crank_nicolson() {
  for (size_t c = 0 ; c < solver_component.size() ; ++c) {
    solver_component[c]->rho = component[c]->rho;
    solver_component[c]->update_boundaries();
  }

  for (auto& all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  for (auto& kv : combinations) {
    solver_component[kv.second.first]->update_density(flux[kv.first]->J,solver_flux[kv.first]->J, cn_ratio);
    solver_component[kv.second.second]->update_density(flux[kv.first]->J,solver_flux[kv.first]->J, cn_ratio, -1);
  }

  for ( size_t n = 0; n < solver_component.size() ; ++n )
    stl::copy(solver_component[n]->rho.begin(), solver_component[n]->rho.end(), rho.begin()+n*system_size);

  #ifdef PAR_MESODYN
  return stl::raw_pointer_cast(rho.data());
  #else
  return rho.data();
  #endif
}

int Mesodyn::initial_conditions() {
  
  stl::host_vector<stl::host_vector<Real>> rho(component_no, stl::host_vector<Real>(system_size));

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

  std::map< std::string, Boundary_type> boundary_settings;
  boundary_settings["mirror"] = Boundary_type::MIRROR;
  boundary_settings["periodic"] = Boundary_type::PERIODIC;

  Boundary_map boundary_conditions;

  // BC0: bX0, BC1: bXm, etc.
  boundary_conditions[Dimension::Z] = boundary_settings[Lat[0]->BC[4]];
  boundary_conditions[Dimension::Y] = boundary_settings[Lat[0]->BC[2]];
  boundary_conditions[Dimension::X] = boundary_settings[Lat[0]->BC[0]];

  Mesodyn::boundary = Boundary_factory::Create(dimensionality, mask, boundary_conditions);

  for (auto& density : densities) {
    component.push_back(make_shared<Component>(Lat[0], boundary, density));
    solver_component.push_back(make_shared<Component>(Lat[0], boundary, density));
  }

  if (seed_specified == true)
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev, seed);
  else
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev);  

  for (auto& index_of : combinations) {
      flux.push_back(
        Flux_factory::Create(dimensionality, Lat[0], D * dt, mask, component[index_of.second.first], component[index_of.second.second], gaussian));
      solver_flux.push_back(
        Flux_factory::Create(dimensionality, Lat[0], D * dt, mask, solver_component[index_of.second.first], solver_component[index_of.second.second], gaussian));
    }

  int j = 0;
  for (auto& molecule : Mol) {
    Molecule_density this_density = Molecule_density(molecule);

    vector<Real> densities = this_density.monomer_total_mass();
    for (auto density : densities) {

      theta[solver_component[j]] = theta[component[j]] = density;

      if (density == 0)
        m_solvent = j;

      ++j;
    }
  }

  norm_theta(component);
  norm_theta(solver_component);

  return 0;
}

int Mesodyn::norm_theta(vector< shared_ptr<Component> >& component_) {

  assert(Mol[ Sys[0]->solvent ]->MolMonList.size() == 1 and "Norming solvents with mutliple monomers is not supported! Please write your own script");
  stl::device_vector<Real> residuals(system_size,0);

  for (auto component : component_) {
      Norm( (Real*)component->rho, ( theta[component]/component->theta() ), system_size);
      stl::transform(component->rho.begin(), component->rho.end(), residuals.begin(), residuals.begin(), stl::plus<Real>());
  }

  stl::transform(residuals.begin(), residuals.end(), component_[m_solvent]->rho.begin(), component_[m_solvent]->rho.begin(),
      [this] DEVICE_LAMBDA (const double& x, const double& y) {
        return ( y-(x-1));
      }
    );

  return 0; 
}

/******* Output generation *******/

void Mesodyn::set_filename() {
  /* Open filestream and set filename to "mesodyn-datetime.csv" */
  string output_folder = "output/";
  string bin_folder = "bin";

  // Find path to Namics executable
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  string executable_path = string(result, (count > 0) ? count : 0);

  // Find the last string before the executable
  size_t found = executable_path.find_last_of("/\\");

  // Set the output folder to be one level up from the binary folder, plus the specified output folder
  output_folder = executable_path.substr(0, found - bin_folder.size()) + output_folder;

  filename << output_folder << "mesodyn-";

  time_t rawtime;
  time(&rawtime);
  filename << rawtime;
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

    Out[0]->push("order_parameter",order_parameter);
    Out[0]->push("time",t);
    Out[0]->push("timesteps", timesteps);
    Out[0]->push("timebetweensaves", timebetweensaves);
    Out[0]->push("diffusionconstant", D);
    Out[0]->push("seed", seed);
    Out[0]->push("mean", mean);
    Out[0]->push("stddev", stddev);
    Out[0]->push("delta_t", dt);
    Out[0]->push("cn_ratio", cn_ratio);

    for (size_t i = 0 ; i < solver_component.size() ; ++i) {
      ostringstream vtk_filename;
      vtk_filename << filename.str() << "-rho" << i << "-" << writes << ".vtk";

      stl::host_vector<Real> t_rho = solver_component[i]->rho.m_data;
      Out[0]->vtk_structured_grid(vtk_filename.str(), t_rho.data(), i);
    }
    
    ++writes;

    Out[0]->WriteOutput(writes);
    Out.clear();

  return 0;
}