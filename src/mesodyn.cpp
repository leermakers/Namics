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
    "grand_cannonical",
    "grand_cannonical_time_average",
    "grand_cannonical_molecule"
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
      stddev                           { initialize<Real>("stddev", (2 * D * sqrt(dt) ) )},
      seed                             { initialize<Real>("seed", -12345.6789) },
      seed_specified                   { seed != -12345.6789 ? true : false },   
      timesteps                        { initialize<size_t>("timesteps", 100) },
      save_delay                       { initialize<size_t>("save_delay", 0) },
      timebetweensaves                 { initialize<size_t>("timebetweensaves", 1) },
      cn_ratio                         { initialize<Real>("cn_ratio", 0.5) },
      enable_sanity_check              { initialize<bool>("sanity_check", 0)},
      output_profile_filetype          { initialize_enum<Writable_filetype>("profile_type", Writable_filetype::VTK_STRUCTURED_GRID, Profile_writer::input_options)},
      grand_cannonical                 { initialize<bool>("grand_cannonical", 0)},
      grand_cannonical_time_average    { initialize<size_t>("grand_cannonical_time_average", timesteps > 100 ? 20 : 5 ) },
      grand_cannonical_molecule        { initialize<size_t>("grand_cannonical_molecule", Sys[0]->solvent == 0 ? 1 : 0)},

      //Variables for rho initialization
      initialization_mode              { INIT_HOMOGENEOUS },
      component_no                     { Sys.front()->SysMolMonList.size() },
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
      initialization_mode = INIT_FROMFILE;

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
  
    vector<Sanity_check*> checks;

    // Attach sanity checks
    for (size_t i = 0 ; i < component_no; ++i) {
      checks.emplace_back(new Check_between_zero_and_one<Real>(&components[i]->rho, i));
      checks.emplace_back(new Check_theta<Real>(&components[i]->rho, std::accumulate(components[i]->rho.begin(), components[i]->rho.end(), 0), i));
    }

    Check_index_unity<Real> check_rho(&components[0]->rho);
    for (size_t i = 1 ; i < component_no ; ++i)
      check_rho.register_checkable(&components[i]->rho);

  //Prepare IO
  
  cout.precision(8);

  cout << "Mesodyn is all set, starting calculations.." << endl;// << endl << endl;

  // Prepare callback functions for SolveMesodyn in Newton
  function<Real*()> solver_callback = bind(&Mesodyn::solve_crank_nicolson, this);
  function<void(Real*,size_t)> loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  //Real grand_potential_average{0};

  //Norm_densities_relative relative_norm( Mol, components, Sys[0]->solvent );

  /**** Main MesoDyn time loop ****/
  //int z = 1;
  for (t = 0; t < timesteps+1; t++) {
    //cout << "\x1b[A" << "\x1b[A" << "MESODYN: t = " << t << " / " << timesteps << endl;
    cout << "MESODYN: t = " << t << " / " << timesteps << endl;

    gaussian->generate(system_size);

    for (auto& all_fluxes : fluxes) all_fluxes->J.save_state();
    for (auto& all_components : components) all_components->rho.save_state();

    New[0]->SolveMesodyn(loader_callback, solver_callback);
    order_parameter->execute();

    cout << "Order parameter: " << order_parameter->attach() << endl;

    if (enable_sanity_check)
      sanity_check();

    if (not In.back()->OutputList.empty() and t > save_delay and t % timebetweensaves == 0)
      write_output();

 /*   grand_potential_average += Sys[0]->GetGrandPotential();

    if (grand_cannonical and t > save_delay && t % grand_cannonical_time_average == 0)
    {
        grand_potential_average /= z;
        z=1;
        Real adjustment = -1;
      
        if( grand_potential_average < 0) {
          relative_norm.adjust_theta(Sys[0]->solvent, -0.001);
          relative_norm.execute();

          for (size_t j = 0 ; j < Mol.size() ; ++j)
            {
              if (j != Sys[0]->solvent)
              {
                Real sum {0};
	              for (size_t i = 0 ; i < Mol[j]->MolMonList.size() ; ++i) {
                  sum += components[Mol[j]->MolMonList[i]]->theta();
                }
                  
                Mol[j]->theta = sum;
              }
            }
        }
        else if( grand_potential_average > 0)
        {
          adjustment = -1.0*static_cast<Real>(system_size)*0.001; //(Mol[grand_cannonical_molecule]->theta*(grand_cannonical_average*10));

          for(int i = 0 ; i < (int)Mol.size() ; ++i)
            if (i != Sys[0]->solvent)
                norm_densities->adjust_theta(i, adjustment);

          norm_densities->execute();
        }

      }

      cout << grand_potential_average/z << endl;
      ++z;
        for ( auto& asdf : Mol) {
          cout << asdf->theta << " ";
        }

        grand_potential_average=0;

        cout << endl;*/
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
        if (bitmask[i]) index.emplace_back(i);

      combinations[ index[0] ] = index[1];

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return combinations;
}

Real* Mesodyn::solve_crank_nicolson() {

  for (auto& all_components : components) {
    all_components->rho.reinstate_previous_state();
    all_components->update_boundaries();
  }

  for (auto& all_fluxes : fluxes)
    all_fluxes->flux();

  update_densities();

  prepare_densities_for_callback();

 /*  for (auto& component : components) {
      stl::host_vector<Real> asdf = component->rho.m_data;
      for (size_t z = 0 ; z < 5 ; ++z)
      for (size_t y = 0 ; y < 5 ; ++y)
      for (size_t x = 0 ; x < 5 ; ++x)
        {
        cout << z << " " << y << " " << x << "  :  ";
        cout << asdf[z*component->rho.m_subject_lattice->JZ+y*component->rho.m_subject_lattice->JY+x*component->rho.m_subject_lattice->JX] << endl;
        }
      cout << endl;
    }*/

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

  std::map<size_t, size_t> combinations = generate_pairs(components.size());

  if (seed_specified == true)
    Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev, seed);
  else
    Mesodyn::gaussian = make_shared<Gaussian_noise>(boundary, mean, stddev);  

  for (auto& index_of : combinations)
    Mesodyn::fluxes.emplace_back(
      Flux::Factory::Create(dimensionality, Lat[0], D * dt, mask, components[index_of.first], components[index_of.second], gaussian));

  Mesodyn::norm_densities = make_unique<Norm_densities>(Mol, components, Sys[0]->solvent);
  Mesodyn::order_parameter = make_unique<Order_parameter>(components, combinations, Sys.front()->boundaryless_volume);

  //TODO: this norm is broken for boundaries, but that doesn't really seem to pose a problem.
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
    for (size_t i = 0 ; i < components.size() ; ++i)
    {
      string description = "component:" + to_string(i);
      register_output_profile(description + ":density", (Real*)components[i]->rho);
    }
}

int Mesodyn::write_output() {
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

    for (auto& parameter_writer : parameter_writers)
      parameter_writer->write();

    for (auto& profile_writer : profile_writers)
    {
      profile_writer->prepare_for_data();
      profile_writer->write();
    }

    return 0;
}

/*      Writable_file kal_file(filename.str(), Writable_filetype::KAL);
    parameter_writers.push_back( make_shared<Kal_writer>(kal_file) );
    for (auto& parameter_writer : parameter_writers) {
      parameter_writer->bind_data(output_params);
      parameter_writer->prepare_for_data(selected_options);
    }

    *******************

    Out.emplace_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[0], (int)t, timesteps / timebetweensaves));
    if (!Out[0]->CheckInput(1)) {
        cout << "input_error in output " << endl;
        exit(0);
    }

    for (size_t i = 0 ; i < Out[0]->OUT_key.size() ; ++i)
      if( Out[0]->OUT_key[i] == "mesodyn")
        selected_options.push_back(Out[0]->OUT_prop[i]);

    register_output_param("time", &t);
    register_output_param("timesteps", &timesteps);
    register_output_param("timebetweensaves", &timebetweensaves);
    register_output_param("diffusionconstant", &D);
    register_output_param("seed", &seed);
    register_output_param("mean", &mean);
    register_output_param("stddev", &stddev);
    register_output_param("delta_t", &dt);
    register_output_param("cn_ratio", &cn_ratio);
    register_output_param("order_parameter", &order_parameter->get());*/