#include "mesodyn.h"

/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018-2019 */
/* Most of the physics in this module is based on the work of Fraaije et al. in the 1990s  */

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
    "save_delay",
    "cn_ratio",
    "sanity_check",
    "profile_type",
    "write_grand_potential",
    "write_free_energy",
    "write_alpha",
    "write_density",
    "grand_cannonical",
    "grand_cannonical_time_average",
    "grand_cannonical_molecule",
    "treat_lower_than_as_zero",
    "adaptive_tolerance_modifier",
    "adaptive_tolerance"
};

Mesodyn::Mesodyn(int start, vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_, string name_)
    : 
      Lattice_accessor(Lat_[0]),
      name{name_}, In{In_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sta{Sta_}, Rea{Rea_}, Sys{Sys_}, New{New_},

      //Const-correct way of initializing member variables from file, see template in header file.
      //Initialization syntax: initialize<Datatype>("option", default_value);
      D                                { initialize<Real>("diffusionconstant", 0.01) },
      dt                               { initialize<Real>("delta_t", 0.1) },
      mean                             { initialize<Real>("mean", 0.0) },
      stddev                           { initialize<Real>("stddev", sqrt(2 * D))},
      seed                             { initialize<Real>("seed", -12345.6789) },
      seed_specified                   { seed != -12345.6789 ? true : false },   
      timesteps                        { initialize<size_t>("timesteps", 100) },
      save_delay                       { initialize<size_t>("save_delay", 0) },
      timebetweensaves                 { initialize<size_t>("timebetweensaves", 1) },
      cn_ratio                         { initialize<Real>("cn_ratio", 0.5) },
      treat_lower_than_as_zero         { initialize<Real>("treat_lower_than_as_zero", New.back()->tolerance*0.1)},
      adaptive_tolerance               { initialize<bool>("adaptive_tolerance", 1)},
      adaptive_tolerance_modifier      { initialize<Real>("adaptive_tolerance_modifier", 100)},
      enable_sanity_check              { initialize<bool>("sanity_check", 0)},
      output_profile_filetype          { initialize_enum<Writable_filetype>("profile_type", Writable_filetype::VTK_STRUCTURED_GRID, Profile_writer::output_options)},
      grand_cannonical                 { initialize<bool>("grand_cannonical", 0)},
      grand_cannonical_time_average    { initialize<size_t>("grand_cannonical_time_average", timesteps > 100 ? 20 : 5 ) },
      grand_cannonical_molecule        { initialize<size_t>("grand_cannonical_molecule", Sys[0]->solvent == 0 ? 1 : 0)},

      //Variables for rho initialization
      initialization_mode              { INIT_HOMOGENEOUS },
      component_no                     { Sys.back()->SysMolMonList.size() },
      t                                { 0 }
{
  // to get the correct KSAM and volume.
  Sys.front()->PrepareForCalculations();

  callback_densities.resize(component_no*system_size);
  
    //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (size_t i = 0; i < Seg.size(); ++i)
  {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }

  CheckInput();

  cout << "Initializing.." << endl;
  initial_conditions();

  register_output();
  set_filename();

  Writable_file out_file(filename.str(), output_profile_filetype );
  profile_writers.push_back(Profile_writer::Factory::Create(output_profile_filetype, Lat[0], out_file));
  profile_writers[0]->bind_data(output_profiles);
}

Mesodyn::~Mesodyn() {
    // We only use smart pointers here, they'll take care of deleting themselves when needed.
}

bool Mesodyn::CheckInput() {
    string empty = "";

    if ( (read_filename = initialize<std::string>("read_pro",empty)) != empty)
      input_data_filetype = Readable_filetype::PRO;
    else if ( (read_filename = initialize<std::string>("read_vtk",empty)) != empty) {
      input_data_filetype = Readable_filetype::VTK_STRUCTURED_GRID;
    }

    if (input_data_filetype != Readable_filetype::NONE)
      initialization_mode = Mesodyn::INIT_FROMFILE;

    if ( find(PARAMETERS.begin(), PARAMETERS.end(), "grand_cannonical_time_average") != PARAMETERS.end() 
      or find(PARAMETERS.begin(), PARAMETERS.end(), "grand_cannonical_molecule") != PARAMETERS.end()  )
        if( grand_cannonical == false )
        {
          cout << "Please enable grand_cannonical in input or remove grand_cannonical options!" << endl;
          throw 1;
        }


  return true;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  
  if (enable_sanity_check) {
    cout << "Binding checks" << endl;
    vector<Sanity_check*> checks;

    // Attach sanity checks
    for (size_t i = 0 ; i < component_no; ++i) {
      checks.emplace_back(new Check_between_zero_and_one<Real>(&components[i]->rho, i));
      checks.emplace_back(new Check_theta<Real>(&components[i]->rho, std::accumulate(components[i]->rho.begin(), components[i]->rho.end(), 0.0), i));
    }

  //  Check_index_unity<Real> check_rho(&components[0]->rho);
  //  for (size_t i = 1 ; i < component_no ; ++i)
  //    check_rho.register_checkable(&components[i]->rho);
  }

  // Prepare IO
  
  cout.precision(8);

  cout << "Mesodyn is all set, starting calculations.." << endl;// << endl << endl;

  // Prepare callback functions for SolveMesodyn in Newton
  function<Real*()> solver_callback = bind(&Mesodyn::solve_crank_nicolson, this);
  function<void(Real*,size_t)> loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  /**** Main MesoDyn time loop ****/
  for (t = 0; t < timesteps+1; t++) {

    cout << "MESODYN: t = " << t << " / " << timesteps << endl;

    gaussian->generate(system_size);

    for (auto& all_fluxes : fluxes) all_fluxes->J.save_state();
    for (auto& all_components : components) all_components->rho.save_state();

    New[0]->SolveMesodyn(loader_callback, solver_callback);

   // norm_densities->execute();

    //Somehow if stencil_full in combination with frozen segments gives wrong densities
    // Breaks periodic boundaries
/*     if (Lat.back()->stencil_full)
      for (auto& all_components : components)
        Times((Real*)all_components->rho, (Real*)all_components->rho, Sys.back()->KSAM, system_size); */

    order_parameter->execute();

    cout << "Order parameter: " << order_parameter->attach() << endl;

    if (enable_sanity_check)
      sanity_check();

    write_parameters();

    if (t > save_delay and t % timebetweensaves == 0)
      write_profile();

    if (adaptive_tolerance) {
      adapt_tolerance();
    }

    //Zero(New.back()->xx, system_size);
    

  } // time loop

  std::cout << "Done." << std::endl;
  return true;
}

std::multimap<size_t, size_t> Mesodyn::generate_pairs(size_t N)
{
    // Returns every combination in a set of size N, e.g.
    // 0,1 - 0,2 - 0,3 - 1,2 - 1,3 - 2,3

    std::string bitmask(2, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    using combimap = std::multimap<size_t, size_t>;
    using combipair = combimap::value_type;

    combimap combinations;

    // print integers and permute bitmask
    do {
        std::deque<size_t> index;

        for (size_t i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) index.push_back(i);
        }
        combinations.insert(combipair(index.front(), index.back()));
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return combinations;
}

void Mesodyn::adapt_tolerance() {
  // Dynamically adjust tolerance so that namics doesn't give us negative values for densities
  auto minmax = stl::minmax_element(callback_densities.begin(), callback_densities.end());
  int magnitude = floor(log10(*minmax.first));
  double magnitude_2 = pow(10, magnitude) / adaptive_tolerance_modifier;
  if (magnitude_2  > treat_lower_than_as_zero and New.back()->tolerance > magnitude_2) {
    cout << "The gradients are getting steeper, lowering tolerance accordingly. New tolerance: " << magnitude_2 << endl;
    New.back()->tolerance = magnitude_2;
  }
    if (New.back()->tolerance < magnitude_2) {
    cout << "The gradients are getting smoother, increasing tolerance accordingly. New tolerance: " << magnitude_2 << endl;
    New.back()->tolerance = magnitude_2;
  } 
}

Real* Mesodyn::solve_crank_nicolson() {
  for (auto& all_components : components) {
    all_components->rho.reinstate_previous_state();
    all_components->update_boundaries();
  //  Times((Real*)all_components->rho, (Real*)all_components->rho, Sys.back()->KSAM, system_size);
  }

  for (auto& all_fluxes : fluxes)
    all_fluxes->flux();

  update_densities();

  enforce_minimum_density->execute();

  prepare_densities_for_callback();

  return device_vector_ptr_to_raw(callback_densities);
}


void Mesodyn::update_densities() {
  for (auto& all_fluxes : fluxes) {
    all_fluxes->component_a->update_density(all_fluxes->J, cn_ratio, +1);
    all_fluxes->component_b->update_density(all_fluxes->J, cn_ratio, -1);
  }

}

void Mesodyn::prepare_densities_for_callback() {
  for ( size_t n = 0; n < components.size() ; ++n )
    stl::copy(components[n]->rho.begin(), components[n]->rho.end(), callback_densities.begin()+n*system_size);
}

Real* Mesodyn::device_vector_ptr_to_raw(stl::device_vector<Real>& input_) {

  #ifdef PAR_MESODYN
    return stl::raw_pointer_cast(input_.data());
  #else
    return input_.data();
  #endif

}

void Mesodyn::sanity_check() {
  for (auto& all_components : components)
    all_components->rho.perform_checks();
}

void Mesodyn::load_alpha(Real* alpha, const size_t i) {
    if (i < component_no)
     dynamic_pointer_cast<Component>( components[i] )->alpha.load_array(alpha,system_size);
}

int Mesodyn::initial_conditions() {
  
  Lattice_object<size_t> mask = load_mask_from_sys();


  vector<Lattice_object<Real>> densities(component_no, Lattice_object<Real>(Lat[0]) );

  if (initialization_mode == INIT_FROMFILE)
    initialize_from_file(densities);
  else //if initialization_mode == INIT_HOMOGENEOUS
    initialize_homogeneous(densities);


  shared_ptr<Boundary1D> boundary = build_boundaries(mask);

  for (auto& density : densities)
    Mesodyn::components.emplace_back(make_shared<Component>(Lat[0], boundary, density));

  std::multimap<size_t, size_t> combinations = generate_pairs(components.size());

  if (seed_specified == true)
    Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev, seed);
  else
    Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev);  

  cout << "Building neighborlists.." << endl;

  for (auto& index_of : combinations)
    if (Lat.back()->stencil_full) {
      assert (three_D == dimensionality and "Full stencil is only supported in 3D when using mesodyn.");
      Mesodyn::fluxes.emplace_back(make_shared<Flux3D_extended_stencil>(Lat[0], D * dt, mask, components[index_of.first], components[index_of.second], gaussian));
    }
    else {
   // if (Seg[index_of.first]->freedom != "frozen" and Seg[index_of.second]->freedom != "frozen")
      Mesodyn::fluxes.emplace_back(
        Flux::Factory::Create(dimensionality, Lat[0], D * dt, mask, components[index_of.first], components[index_of.second], gaussian));
    }

  Mesodyn::norm_densities = make_unique<Norm_densities>(Mol, components, Sys[0]->solvent);
  Mesodyn::order_parameter = make_unique<Order_parameter>(components, combinations, Sys.front()->boundaryless_volume);
  Mesodyn::enforce_minimum_density = make_unique<Treat_as_zero>(components, treat_lower_than_as_zero);

  //norm_densities->execute();

  return 0; 
}

Lattice_object<size_t> Mesodyn::load_mask_from_sys() {
  Lattice_object<int> t_mask(Lat[0]);

  #if defined(PAR_MESODYN) || ! defined(CUDA)
  stl::copy(Sys[0]->KSAM, Sys[0]->KSAM+system_size, t_mask.begin());
  #else
  TransferDataToHost(t_mask.data(), Sys[0]->KSAM, system_size);
  #endif

  Lattice_object<size_t> mask( t_mask );

  return mask; 
}

shared_ptr<Boundary1D> Mesodyn::build_boundaries(const Lattice_object<size_t>& mask) {
  Boundary::Map boundary_conditions;

  // BC0: bX0, BC1: bXm, etc.
  boundary_conditions[Dimension::X] = Boundary::Adapter[Lat[0]->BC[0]];
  boundary_conditions[Dimension::Y] = Boundary::Adapter[Lat[0]->BC[2]];
  boundary_conditions[Dimension::Z] = Boundary::Adapter[Lat[0]->BC[4]];

  return Boundary::Factory::Create(dimensionality, mask, boundary_conditions);
}

void Mesodyn::initialize_homogeneous(vector<Lattice_object<Real>>& densities) {

  Homogeneous_system_initializer initializer(Sys[0]);
  initializer.build_objects();
  initializer.push_data_to_objects(densities);

}

void Mesodyn::initialize_from_file(vector<Lattice_object<Real>>& densities) {

  Readable_file file(read_filename, Mesodyn::input_data_filetype);
  Reader file_reader;
  file_reader.read_objects_in(file);
  file_reader.assert_lattice_compatible(Lat[0]);
  file_reader.push_data_to_objects(densities);

}

/******* Output generation *******/

void Mesodyn::set_filename() {
  filename << In[0]->output_info.getOutputPath() << "mesodyn-";
  filename << time(time_t());
}

void Mesodyn::register_output() {
    if (initialize<bool>("write_density", 1))
      for (size_t i = 0 ; i < components.size() ; ++i)
      {
        string description = "component:" + to_string(i);
        register_output_profile(description + ":density", (Real*)components[i]->rho);
      }

    if (initialize<bool>("write_alpha", 0))
      for (size_t i = 0 ; i < components.size() ; ++i)
      {
        string description = "component:" + to_string(i);
        register_output_profile(description + ":alpha", (Real*)components[i]->alpha);
      }

    if (initialize<bool>("write_free_energy", 0))
      register_output_profile("free_energy_density", Sys.back()->FreeEnergyDensity);

    if (initialize<bool>("write_grand_potential", 0))
      register_output_profile("grand_potential_density", Sys.back()->GrandPotentialDensity);
}

int Mesodyn::write_profile() {
    for (auto& parameter_writer : parameter_writers)
      parameter_writer->write();

    for (auto& pair : output_profiles) {
      dynamic_cast<Output_ptr<Real>*>(pair.second.get())->set_buffer(system_size);
    }

    for (auto& profile_writer : profile_writers)
    {
      profile_writer->prepare_for_data();
      profile_writer->write();
    }

    for (auto& pair : output_profiles) {
      dynamic_cast<Output_ptr<Real>*>(pair.second.get())->clear_buffer();
    }

    return 0;
}

void Mesodyn::write_parameters() {
    if (not In.back()->OutputList.empty()) {
     Out.emplace_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[0], (int)t, timesteps / timebetweensaves));
     Out[0]->CheckInput(1);
     Out[0]->output_nr = t;
     Out[0]->n_output = timesteps / timebetweensaves;
     New[0]->PushOutput(); 

     Out[0]->push("filename", filename.str());
     Out[0]->push("order_parameter", order_parameter->attach());
     Out[0]->push("time",(int)t);
     Out[0]->push("timesteps", (int)timesteps);
     Out[0]->push("timebetweensaves", (int)timebetweensaves);
     Out[0]->push("diffusionconstant", D);
     Out[0]->push("seed", seed);
     Out[0]->push("mean", mean);
     Out[0]->push("stddev", stddev);
     Out[0]->push("delta_t", dt);
     Out[0]->push("cn_ratio", cn_ratio);

     Out[0]->WriteOutput(t);
     Out.clear();
  }
}