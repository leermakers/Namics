#include "mesodyn.h"
#include "mesodyn/density_initializer.h"
#include "mesodyn/file_reader.h"
#include "mesodyn/neighborlist.h"
#include <iterator>

//TODO: noise: visitor pattern

/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018 */
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

      //Const-correct way of initializing member variables from file, see template in header file.
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

  set_update_lists();

  // to get the correct KSAM and volume.
  Sys[0]->PrepareForCalculations();

  rho.resize(component_no*system_size);

  // The right hand side of the minus sign calculates the volume of the boundaries. I know, it's hideous, but it works for all dimensions.
  boundaryless_volume = Sys[0]->volume - ((2 * dimensionality - 4) * Lat[0]->MX * Lat[0]->MY + 2 * Lat[0]->MX * Lat[0]->MZ + 2 * Lat[0]->MY * Lat[0]->MZ + (-2 + 2 * dimensionality) * (Lat[0]->MX + Lat[0]->MY + Lat[0]->MZ) + pow(2, static_cast<float>(dimensionality)));
}

Mesodyn::~Mesodyn() {
    // We only use smart pointers here, they'll take care of deleting themselves when needed.
}

bool Mesodyn::CheckInput(const int start) {

  input_success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);

  // TODO: implement this properly
  // If the user has asked for vtk output
  //if (std::find(In[0]->OutputList.begin(), In[0]->OutputList.end(), "vtk") != In[0]->OutputList.end())
  //  write_vtk = true;

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
  for (size_t i = 0 ; i < component_no; ++i) {
    checks.push_back(new Check_between_zero_and_one<Real>(&solver_component[i]->rho, i));
 //   checks.push_back(new Check_theta<Real>(&solver_component[i]->rho, solver_component[i]->theta(), i));
  }

  Check_index_unity<Real> check_rho(&solver_component[0]->rho);
  for (size_t i = 1 ; i < component_no ; ++i)
    check_rho.register_checkable(&solver_component[i]->rho);

  //Prepare IO
  set_filename();

  cout << "Mesodyn is all set, starting calculations.." << endl;

  // Do one explicit step before starting the crank nicolson scheme
  //explicit_start();

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

    //noise_flux();

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

void Mesodyn::set_update_lists() {
  int c = 0;
  vector<int> temp;
  for (size_t i = 0; i < component_no; ++i) {
    for (size_t j = 0; j < (component_no - 1) - i; ++j) {
      temp.push_back(i);
      temp.push_back(c);
      update_plus.push_back(temp);
      temp.clear();
      ++c;
    }
  }

  c = 0;
  for (size_t j = 0; j < (component_no - 1); ++j) {
    for (size_t i = 1 + j; i < component_no; ++i) {
      temp.push_back(i);
      temp.push_back(c);
      update_minus.push_back(temp);
      temp.clear();
      ++c;
    }
  }
}

Real Mesodyn::calculate_order_parameter() {
  stl::device_vector<Real> difference(system_size);
  Mesodyn::order_parameter = 0;
  for (size_t i = 0; i < component_no - 1; ++i)
      for (size_t j = i + 1; j < component_no; ++j) {
          stl::transform(solver_component[i]->rho.begin(),
                         solver_component[i]->rho.end(),
                         solver_component[j]->rho.begin(),
                         difference.begin(),
                         order_param_functor());
          #ifdef PAR_MESODYN
          Mesodyn::order_parameter = thrust::reduce(difference.begin(), difference.end(), Mesodyn::order_parameter);
          #else
          Mesodyn::order_parameter = std::accumulate(difference.begin(), difference.end(), Mesodyn::order_parameter);
          #endif
      }
  Mesodyn::order_parameter /= boundaryless_volume;

  return Mesodyn::order_parameter;
}

int Mesodyn::noise_flux() {

  for (auto all_components : solver_component) {
    gaussian->generate(system_size);
    stl::fill(all_components->alpha.begin(), all_components->alpha.end(), 0);
    gaussian->add_noise(all_components->alpha);
  }

  for (auto& all_fluxes : solver_flux) {
    all_fluxes->langevin_flux();
  }

  for (vector<int>& i : update_plus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J);

  for (vector<int>& i : update_minus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J, -1);

  for (size_t i = 0 ; i < flux.size() ; ++i)
    stl::transform(flux[i]->J.begin(), flux[i]->J.end(), solver_flux[i]->J.begin(), flux[i]->J.begin(), std::plus<Real>());

  return 0;
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

Real* Mesodyn::solve_explicit() {
  for (auto all_components : solver_component)
    all_components->update_boundaries();

  for (auto& all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  for (size_t n = 0; n < solver_component.size() ; ++n)
    stl::copy(solver_component[n]->rho.begin(), solver_component[n]->rho.end(), rho.begin()+n*system_size);

  #ifdef PAR_MESODYN
  return stl::raw_pointer_cast(&rho[0]);
  #else
  return &rho[0];
  #endif
}

Real* Mesodyn::solve_crank_nicolson() {
  for (size_t c = 0 ; c < solver_component.size() ; ++c)
    solver_component[c]->rho = component[c]->rho;

  for (auto all_components : solver_component)
    all_components->update_boundaries();

  for (auto& all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  for (vector<int>& i : update_plus)
      solver_component[i[0]]->update_density(flux[i[1]]->J, solver_flux[i[1]]->J, cn_ratio);

  for (vector<int>& i : update_minus)
      solver_component[i[0]]->update_density(flux[i[1]]->J, solver_flux[i[1]]->J, cn_ratio, -1);

  for ( size_t n = 0; n < solver_component.size() ; ++n ) {
//    norm_theta(solver_component);
    stl::copy(solver_component[n]->rho.begin(), solver_component[n]->rho.end(), rho.begin()+n*system_size);
  }

  #ifdef PAR_MESODYN
  return stl::raw_pointer_cast(&rho[0]);
  #else
  return &rho[0];
  #endif
}

int Mesodyn::initial_conditions() {
  
  //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (size_t i = 0; i < Seg.size(); ++i) {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }


  stl::host_vector<stl::host_vector<Real>> rho(component_no, stl::host_vector<Real>(system_size));

  Lattice_object<int> t_mask(Lat[0]);

  #if defined(PAR_MESODYN) || ! defined(CUDA)
  stl::copy(Sys[0]->KSAM, Sys[0]->KSAM+system_size, t_mask.begin());
  #else
  TransferDataToHost(mask.data(), Sys[0]->KSAM, system_size);
  #endif

  Lattice_object<size_t> mask( t_mask ); 

  vector<Lattice_object<Real>> densities(component_no, Lattice_object<Real>(Lat[0]) );
  

  switch (initialization_mode) {
    case INIT_HOMOGENEOUS:
      {
        densities = Homogeneous_system_initializer(Sys[0]).build_objects();
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
  

  switch (dimensionality) {
  case one_D:
    Mesodyn::boundary = make_shared<Boundary1D>(mask, boundary_conditions);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(make_shared<Component>(Lat[0], boundary, densities[i]));
      solver_component.push_back(make_shared<Component>(Lat[0], boundary, densities[i]));
    }

    if (seed_specified == true) {
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, D, system_size, mean, stddev, seed);
    } else {
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, D, system_size, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(make_unique<Flux1D>(Lat[0], D * dt, mask, component[i], component[j], gaussian));
        solver_flux.push_back(make_unique<Flux1D>(Lat[0], D * dt, mask, solver_component[i], solver_component[j], gaussian));
      }
    }
    break;
  case two_D:
    Mesodyn::boundary = make_shared<Boundary2D>(mask, boundary_conditions);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(make_shared<Component>(Lat[0], boundary, densities[i]));
      solver_component.push_back(make_shared<Component>(Lat[0], boundary, densities[i]));
    }

    if (seed_specified == true) {
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, D, system_size, mean, stddev, seed);
    } else {
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, D, system_size, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(make_unique<Flux2D>(Lat[0], D * dt, mask, component[i], component[j], gaussian));
        solver_flux.push_back(make_unique<Flux2D>(Lat[0], D * dt, mask, solver_component[i], solver_component[j], gaussian));
      }
    }
    break;
  case three_D:
    
    Mesodyn::boundary = make_shared<Boundary3D>(mask, boundary_conditions);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(make_shared<Component>(Lat[0], boundary, densities[i]));
      solver_component.push_back(make_shared<Component>(Lat[0], boundary, densities[i]));
    }

    if (seed_specified == true) {
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, D, system_size, mean, stddev, seed);
    } else {
      Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, D, system_size, mean, stddev);
    }


    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(make_unique<Flux3D>(Lat[0], D * dt, mask, component[i], component[j], gaussian));
        solver_flux.push_back(make_unique<Flux3D>(Lat[0], D * dt, mask, solver_component[i], solver_component[j], gaussian));
      }
    }

    break;
  }

  int j = 0;
  for (size_t i = 0 ; i < Mol.size() ; ++i) {
    Molecule_density this_density = Molecule_density(Mol[i]);
    vector<Real> densities = this_density.monomer_total_mass();
    for (size_t z = 0 ; z < densities.size() ; ++z) {
      theta[solver_component[j]] = densities[z];
      ++j;
    }
  }

  norm_theta(component);
  norm_theta(solver_component);

  return 0;
}

void Mesodyn::explicit_start() {
  // Prepare callbacks
  auto explicit_solver_callback = bind(&Mesodyn::solve_explicit, this);
  auto loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);
  
  // start Crank-Nicolson using one explicit step
  New[0]->SolveMesodyn(loader_callback, explicit_solver_callback);

  // Update densities after first timestep
  for (vector<int>& i : update_plus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J);

  for (vector<int>& i : update_minus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J, -1);

  // Load rho_k+1 into rho_k
  int i = 0;
  for (auto all_components : component) {
    all_components->rho = solver_component[i]->rho;
    ++i;
  }

  // Load J_k+1 into J_k
  i = 0;
  for (auto& all_fluxes : flux) {
    all_fluxes->J = solver_flux[i]->J;
    ++i;
  }
}

int Mesodyn::norm_theta(vector< shared_ptr<Component> >& component) {
  size_t solvent = Sys[0]->solvent;
  size_t solvent_component;
  stl::device_vector<Real> residuals(system_size,0);

  for (size_t i = 0 ; i < component.size() ; ++i) {
      Real sum_of_elements{0};
      if (theta[solver_component[i]] == 0)
        solvent_component = i;

      sum_of_elements = component[i]->theta();
      Norm( (Real*)component[i]->rho,(theta[solver_component[i]]/sum_of_elements),system_size);

      // We now know the total density and can adjust the solvent accodingly to add up to 1.
      // Let's first find out how much there is to adjust.
      stl::transform(component[i]->rho.begin(), component[i]->rho.end(), residuals.begin(), residuals.begin(), stl::plus<Real>());
  }

  stl::device_vector<Real> one(system_size,1);

  // Calculate excesss / defecit
  stl::transform(residuals.begin(), residuals.end(), one.begin(),residuals.begin(), stl::minus<Real>());

  // If there's only one solvent mon, this problem is easy.
  if (Mol[solvent]->MolMonList.size() == 1) {
    stl::transform(residuals.begin(), residuals.end(), component[solvent_component]->rho.begin(), component[solvent_component]->rho.begin(), saxpy_functor(-1));
  } else {
    cerr << "Norming solvents with mutliple monomers is not supported! Please write your own script" << endl;
    throw ERROR_FILE_FORMAT;
  }

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

    int component_count{1};
    for (auto all_components : solver_component) {
      stl::host_vector<Real> t_rho = all_components->rho.m_data;
      ostringstream vtk_filename;
      vtk_filename << filename.str() << "-rho" << component_count << "-" << writes << ".vtk";
      Out[0]->vtk_structured_grid(vtk_filename.str(), t_rho.data(), component_count);
      ++component_count;
    }
    ++writes;

    Out[0]->WriteOutput(writes);
    Out.clear();

  return 0;
}

/******* FLUX: TOOLS FOR CALCULATING FLUXES BETWEEN 1 PAIR OF COMPONENTS, HANDLING OF SOLIDS *********/

Flux1D::Flux1D(Lattice* Lat, const Real D, Lattice_object<size_t>& mask, shared_ptr<Component> A, shared_ptr<Component> B, shared_ptr<Gaussian_noise> gaussian)
    : J_plus(Lat), J_minus(Lat), J(Lat), A{A}, B{B}, L(Lat), mu(Lat), t_L(Lat), t_mu(Lat), D{D}, gaussian(gaussian)
  {
  
  Neighborlist_config configuration {
    Dimension::X,
    Direction::plus,
    1
  };

  shared_ptr<Neighborlist> x_neighborlist = make_shared<Neighborlist>(mask);
  x_neighborlist->register_config(configuration);
  x_neighborlist->build();
  attach_neighborlists(x_neighborlist, Dimension::X);
}

Flux2D::Flux2D(Lattice* Lat, const Real D, Lattice_object<size_t>& mask, shared_ptr<Component> A, shared_ptr<Component> B, shared_ptr<Gaussian_noise> gaussian)
    : Flux1D(Lat, D, mask, A, B, gaussian)
  {

  Neighborlist_config configuration {
    Dimension::Y,
    Direction::plus,
    1
  };

  shared_ptr<Neighborlist> y_neighborlist = make_shared<Neighborlist>(mask);
  y_neighborlist->register_config(configuration);
  y_neighborlist->build();
  attach_neighborlists(y_neighborlist, Dimension::Y);
}

Flux3D::Flux3D(Lattice* Lat, const Real D, Lattice_object<size_t>& mask, shared_ptr<Component> A, shared_ptr<Component> B, shared_ptr<Gaussian_noise> gaussian)
    : Flux2D(Lat, D, mask, A, B, gaussian)
  {
  
  Neighborlist_config configuration {
    Dimension::Z,
    Direction::plus,
    1
  };

  shared_ptr<Neighborlist> z_neighborlist = make_shared<Neighborlist>(mask);
  z_neighborlist->register_config(configuration);
  z_neighborlist->build();

  attach_neighborlists(z_neighborlist, Dimension::Z);
}

Flux1D::~Flux1D() {
}

Flux2D::~Flux2D() {
}

Flux3D::~Flux3D() {
}

void Flux1D::attach_neighborlists(shared_ptr<Neighborlist> neighborlist, Dimension dimension) {
  J.attach_neighborlist(neighborlist, dimension);
  J_plus.attach_neighborlist(neighborlist, dimension);
  J_minus.attach_neighborlist(neighborlist, dimension);
  L.attach_neighborlist(neighborlist, dimension);
  mu.attach_neighborlist(neighborlist, dimension);
  t_L.attach_neighborlist(neighborlist, dimension);
  t_mu.attach_neighborlist(neighborlist, dimension);
}

int Flux1D::langevin_flux() {

  //Zero (with bounds checking) vector J before use
  stl::fill(J.begin(), J.end(), 0);

  if (A->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  onsager_coefficient(A->rho, B->rho);
  potential_difference(A->alpha, B->alpha);

  langevin_flux(Dimension::X);

  return 0;
}

int Flux2D::langevin_flux() {
  Flux1D::langevin_flux();

  Flux1D::langevin_flux(Dimension::Y);

  return 0;
}

int Flux3D::langevin_flux() {
  Flux2D::langevin_flux();

  Flux1D::langevin_flux(Dimension::Z);

  return 0;
}

int Flux1D::onsager_coefficient(Lattice_object<Real>& A, Lattice_object<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  stl::transform(A.begin(), A.end(), B.begin(), L.begin(), stl::multiplies<Real>());

  return 0;
}

int Flux1D::potential_difference(Lattice_object<Real>& A, Lattice_object<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(A.begin(), A.end(), B.begin(), mu.begin(), stl::minus<Real>());

  gaussian->add_noise(mu);

  return 0;
}

int Flux1D::langevin_flux(Dimension dimension) {
  stl::fill(J_minus.begin(), J_minus.end(), 0);
  stl::fill(J_plus.begin(), J_plus.end(), 0);

  stl::transform(mu.available_neighbors[dimension]->begin(), mu.available_neighbors[dimension]->end(), mu.available_sites->begin(), t_mu.available_sites->begin(), stl::minus<Real>());
  stl::transform(L.available_sites->begin(), L.available_sites->end(), L.available_neighbors[dimension]->begin(), t_L.available_sites->begin(), stl::plus<Real>());
  stl::transform(t_mu.begin(), t_mu.end(), t_L.begin(), J_plus.begin(), const_multiply_functor(-D) );

  stl::transform(J_plus.available_sites->begin(), J_plus.available_sites->end(), J_minus.available_neighbors[dimension]->begin(), stl::negate<Real>());

  stl::transform(J_plus.begin(), J_plus.end(), J.begin(), J.begin(), stl::plus<Real>());
  stl::transform(J_minus.begin(), J_minus.end(), J.begin(), J.begin(), stl::plus<Real>());

  //#else
  //   for (auto itt = mask_plus.begin() ; itt < mask_plus.end(); ++itt) {
  //     auto z = *itt;
  //     J_plus[z] = -D * ((L[z] + L[z + jump]) * (mu[z + jump] - mu[z]));
  //  }
  //
  //  for (auto itt = mask_minus.begin() ; itt < mask_minus.end(); ++itt) {
  //     auto z = *itt;
  //    J_minus[z] = -J_plus[z - jump];
  //  }
  //#endif

  return 0;
}

/****************** COMPONENT: DENSITY PROFILE STORAGE AND UPDATING, BOUNDARY CONDITIONS ********************/

/******* Constructors *******/

Component::Component(Lattice* Lat, shared_ptr<Boundary1D> boundary, Lattice_object<Real> rho)
//TODO: fix alpha size.
    : rho(rho), alpha(Lat), Lat{Lat}, boundary(boundary) {
  update_boundaries();
}

Component::~Component() {
}

/******* Interface *******/

int Component::update_density(Lattice_object<Real>& J, int sign) {
  //Explicit update

  if (J.size() != rho.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(J.begin(), J.end(), rho.begin(), rho.begin(), saxpy_functor(sign) );

  return 0;
}

int Component::update_density(Lattice_object<Real>& J1, Lattice_object<Real>& J2, Real ratio, int sign) {
  //Implicit update
  if (J1.size() != rho.size() || J1.size() != J2.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  // Rho <- A * J1 + Rho
  stl::transform(J1.begin(), J1.end(), rho.begin(), rho.begin(), saxpy_functor(sign*ratio) );
  stl::transform(J2.begin(), J2.end(), rho.begin(), rho.begin(), saxpy_functor((1.0f-ratio)*sign) );

  return 0;
}

int Component::update_boundaries() {
  boundary->update_boundaries( alpha.m_data );
  boundary->update_boundaries( rho.m_data );
  return 0;
}

Real Component::theta() {
  //TODO: update once rho gets a neighborlist
  Real sum{0};
  boundary->zero_boundaries( rho.m_data );
  #ifdef PAR_MESODYN
  sum = stl::reduce(rho.begin(), rho.end(), sum);
  #else
  sum = stl::accumulate(rho.begin(), rho.end(), sum);
  #endif
  boundary->update_boundaries( rho.m_data );
  return sum;
}

/******* GAUSSIAN_NOISE: GENERATES WHITE NOISE FOR FLUXES ********/

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real D, int size, Real mean, Real stddev) : noise(size), prng{std::random_device{}()}, dist(mean, stddev), boundary{boundary} {}

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real D, int size, Real mean, Real stddev, size_t seed) : noise(size), prng(seed), dist(mean, stddev), boundary{boundary} {}

int Gaussian_noise::generate(size_t system_size) {
  stl::host_vector<Real> tmp_noise(system_size);
  for (Real& value : tmp_noise)
    value = dist(prng);

  noise = tmp_noise;

  boundary->update_boundaries(noise);

  return 0;
}

int Gaussian_noise::add_noise(stl::device_vector<Real>& target) {
  stl::transform(noise.begin(), noise.end(), target.begin(), target.begin(), stl::plus<Real>());
  return 0;
}


int Gaussian_noise::add_noise(Lattice_object<Real>& target) {
  stl::transform(noise.begin(), noise.end(), target.begin(), target.begin(), stl::plus<Real>());
  return 0;
}

/******* TOOLS: MATHEMATICS AND INTERFACE FOR THE IO CLASSES ********/

/***** Mathematics *****/


int Mesodyn::factorial(int n) {
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
  return factorial(n) / (factorial(n - k) * factorial(k));
}
