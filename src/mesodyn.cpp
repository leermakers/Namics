#include "mesodyn.h"

/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018 */
/* Most of the physics in this module is based on the work of Fraaije et al. in the 1990s  */

Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_, string name_)
    : Lattice_Access(Lat_[0]),
      name{name_}, In{In_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sys{Sys_}, New{New_},
      D{0.01}, dt{0.1}, mean{0}, stddev{2 * D * sqrt(dt)}, seed{1}, seed_specified{false}, timesteps{100}, timebetweensaves{1},
      initialization_mode{INIT_HOMOGENEOUS},
      component_no{Sys[0]->SysMolMonList.size()}, edge_detection{false}, edge_detection_threshold{50}, RC{0}, cn_ratio{0.5},
      component(0), flux(0),
      writes{0}, write_vtk{false} {
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("delta_t");
  KEYS.push_back("read_pro");
  KEYS.push_back("read_vtk");
  KEYS.push_back("equilibrate");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stddev");
  KEYS.push_back("cn_ratio");
  KEYS.push_back("edge_detection");
  KEYS.push_back("edge_detection_threshold");

  // TODO: implement this properly
  // If the user has asked for vtk output
  if (std::find(In[0]->OutputList.begin(), In[0]->OutputList.end(), "vtk") != In[0]->OutputList.end()) {
    write_vtk = true;
  }

  if (debug)
    cout << "Mesodyn initialized." << endl;
}

Mesodyn::~Mesodyn() {
  for (size_t i = 0; i < flux.size(); ++i)
    delete flux[i];
  flux.clear();
  for (size_t i = 0; i < component.size(); ++i)
    delete component[i];
  component.clear();
  delete boundary;
}

bool Mesodyn::CheckInput(int start) {
  if (debug)
    cout << "CheckInput in Mesodyn" << endl;
  bool success = true;

  /* When adding new items here, please add them to the function prepareOutputFile()*/
  success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);

  if (success) {
    vector<string> options;
    if (GetValue("timesteps").size() > 0) {
      success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 100000000, "The number of timesteps should be between 1 and 10000");
    }
    if (debug)
      cout << "Timesteps is " << timesteps << endl;

    if (GetValue("timebetweensaves").size() > 0) {
      timebetweensaves = In[0]->Get_int(GetValue("timebetweensaves"), timebetweensaves);
    }
    if (debug)
      cout << "Time bewteen saves is " << timebetweensaves << endl;

    if (GetValue("delta_t").size() > 0) {
      dt = In[0]->Get_Real(GetValue("delta_t"), dt);
    }
    if (debug)
      cout << "Delta t is " << dt << endl;

    bool equilibrate = false;
    if (GetValue("equilibrate").size() > 0) {
      if (In[0]->Get_bool(GetValue("equilibrate"), equilibrate))
        initialization_mode = INIT_EQUILIBRATE;
    }

    if (GetValue("read_pro").size() > 0 && equilibrate)
      cout << "WARNING: read_pro overrides equilibrate!" << endl;

    if (GetValue("read_pro").size() > 0) {
      read_filename = In[0]->Get_string(GetValue("read_pro"), read_filename);
      initialization_mode = INIT_FROMPRO;
    }
    if (debug)
      cout << "Filename to read rho from is: " << read_filename << endl;

    if (GetValue("read_vtk").size() > 0) {
      read_filename = In[0]->Get_string(GetValue("read_vtk"), read_filename);
      if (read_filename.find(".vtk") != string::npos) {
        cerr << "Mesodyn will add the component number and extension by itself (in that order), please format the remainder of the filename accordingly." << endl;
        exit(0);
      }
      initialization_mode = INIT_FROMVTK;
    }
    if (debug)
      cout << "Filename to read rho from is: " << read_filename << endl;

    if (GetValue("diffusionconstant").size() > 0) {
      D = In[0]->Get_Real(GetValue("diffusionconstant"), D);
    }
    if (debug)
      cout << "Diffusion const is " << D << endl;

    if (GetValue("seed").size() > 0) {
      seed_specified = true;
      seed = In[0]->Get_Real(GetValue("seed"), seed);
    }
    if (debug)
      cout << "Seed is " << seed << endl;

    if (GetValue("mean").size() > 0) {
      mean = In[0]->Get_Real(GetValue("mean"), mean);
    }
    if (debug)
      cout << "Mean is " << mean << endl;

    if (GetValue("stddev").size() > 0) {
      stddev = In[0]->Get_Real(GetValue("stddev"), stddev);
    }
    if (debug)
      cout << "Stdev is " << stddev << endl;

    if (GetValue("cn_ratio").size() > 0) {
      cn_ratio = In[0]->Get_Real(GetValue("cn_ratio"), cn_ratio);
    }
    if (debug)
      cout << "Cn_ratio is " << cn_ratio << endl;

    if(GetValue("edge_detection").size() > 0) {
      vector<string> options;
			options.push_back("sobel");

			if(GetValue("edge_detection") == options[0]) {
          edge_detection = true;
			} else {
        cout << "edge_detection algorithm " << GetValue("edge_detection") << " not recognized, defaulting to no edge detection" << endl;
      }
    }

    if(GetValue("edge_detection_threshold").size() > 0) {
      edge_detection_threshold = In[0]->Get_int(GetValue("edge_detection_threshold"), edge_detection_threshold);
    }
  }

  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  if (debug)
    cout << "mesodyn in Mesodyn" << endl;

  //   Initialize densities
  initial_conditions();

  //Prepare IO
  set_filename();
  write_settings();

  // write initial conditions
  write_density(solver_component);

  cout << "Mesodyn is all set, starting calculations.." << endl;

  // Do one explicit step before starting the crank nicolson scheme

  explicit_start();

  write_density(solver_component);
  write_output();

  // Prepare callback functions for SolveMesodyn in Newton
  auto solver_callback = bind(&Mesodyn::solve_crank_nicolson, this);
  auto loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  string file = filename.str(); // TODO: remove

  /**** Main MesoDyn time loop ****/
  for (int t = 1; t < timesteps; t++) {
    cout << "MESODYN: t = " << t << endl;

    norm_theta(component);
    norm_theta(solver_component);

      for (Flux1D* all_fluxes : solver_flux) {
        all_fluxes->gaussian->generate();
      }

    New[0]->SolveMesodyn(loader_callback, solver_callback);

    //sanity_check();

    int i = 0;
    for (Component* all_components : component) {
      all_components->load_rho(solver_component[i]->rho);
      ++i;
    }

    //TODO: remove this check?
    skip_bounds([this](int x, int y, int z) mutable {
    for (Component* all_components : solver_component)
      if (val(all_components->rho, x, y, z) < 0 || val(rho, x, y, z) > 1 )
        cerr << "CRITICAL ERROR IN DENSITIES AT " << x << "," << y << "," << z  << endl;
    });

    if (t % timebetweensaves == 0) {
      write_density(component);
      write_output();

      if (edge_detection) {
        interface->detect_edges(edge_detection_threshold);
        interface->write_edges(file);
      }
    }

    i = 0;
    for (Flux1D* all_fluxes : flux) {
      all_fluxes->J = solver_flux[i]->J;
      ++i;
    }

  } // time loop



  return true;
}

int Mesodyn::sanity_check() {

  //Check mass conservation
  int i = 0;
  for (Component* all_components : component) {
    if ( fabs(all_components->theta() - solver_component[i]->theta() ) > numeric_limits<Real>::epsilon() ) {
      Real difference = all_components->theta()-solver_component[i]->theta();
        if (difference != 0) {
          cerr << "WARNING: Mass of component " << i << " is NOT conserved! ";
          if (difference > 0)
            cerr << "You're losing mass, difference: ";
          else
            cerr << "You're gaining mass, difference: ";
          cerr << difference << endl;
          }
      }
      ++i;
    }

  return 0;
}

void Mesodyn::load_alpha(vector<Real>& alpha, size_t i) {
    if (i < component_no) {
      solver_component[i]->load_alpha(alpha);
    }
}

Real* Mesodyn::solve_explicit() {
  for (Component* all_components : solver_component)
    all_components->update_boundaries();

  for (Flux1D* all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  rho.clear();
  for (Component* all_components : solver_component) {
    copy(all_components->rho.begin(), all_components->rho.end(), back_inserter(rho));
  }

  return &rho[0];
}

Real* Mesodyn::solve_crank_nicolson() {
  for (Component* all_components : solver_component)
    all_components->update_boundaries();

  for (Flux1D* all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  int c = 0;
  for (size_t i = 0; i < component_no; ++i) {
    for (size_t j = 0; j < (component_no - 1) - i; ++j) {
      solver_component[i]->update_density(component[i]->rho, flux[c]->J, solver_flux[c]->J, cn_ratio);
      ++c;
    }
  }

  c = 0;
  for (size_t j = 0; j < (component_no - 1); ++j) {
    for (size_t i = 1 + j; i < component_no; ++i) {
      solver_component[i]->update_density(component[i]->rho, flux[c]->J, solver_flux[c]->J, cn_ratio, -1.0);
      ++c;
    }
  }

  rho.clear();
  for (Component* all_components : solver_component) {
    copy(all_components->rho.begin(), all_components->rho.end(), back_inserter(rho));
  }

  return &rho[0];
}

int Mesodyn::initial_conditions() {

  //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (size_t i = 0; i < Seg.size(); ++i) {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }

  vector<vector<Real>> rho(component_no, vector<Real>(M));

  Sys[0]->PrepareForCalculations(); // to get the correct KSAM.

  //TODO: once KSAM becomes a vector, this loop won't be nessecary anymore, just pass KSAM to the flux constructor.
  vector<int> mask(M);
  for (int i = 0; i < M; ++i) {
    mask[i] = *(Sys[0]->KSAM + i);
  }

  switch (initialization_mode) {
    case INIT_HOMOGENEOUS:
      init_rho_homogeneous(rho, mask);
      break;
    case INIT_FROMPRO:
      init_rho_frompro(rho, read_filename);
      break;
    case INIT_FROMVTK:
      for (size_t i = 0 ; i < rho.size() ; ++i) {
        string filename = read_filename + to_string(i) + ".vtk";
        init_rho_fromvtk(rho[i], filename);
      }
      break;
    case INIT_EQUILIBRATE:
      init_rho_equilibrate(rho);
      break;
  }

  // BC0: bX0, BC1: bXm, etc.
  vector<Boundary1D::boundary> boundaries;
  for (string& boundary : Lat[0]->BC) {
    if (boundary == "mirror") {
      boundaries.push_back(Boundary1D::MIRROR);
    }
    if (boundary == "periodic") {
      boundaries.push_back(Boundary1D::PERIODIC);
    }
  }

  //Because noise per site messes up mirror boundaries, we need to tell mesodyn to mask them.
  if (boundaries[0] == Boundary1D::MIRROR) {
    x0_boundary( [this, &mask] (int x, int y, int z) mutable {
      *val_ptr(mask, x, y, z) = 0;
    });
  }
  if (boundaries[1] == Boundary1D::MIRROR) {
    xm_boundary( [this, &mask] (int x, int y, int z) mutable {
      *val_ptr(mask, x, y, z) = 0;
    });
  }
  if (dimensions > 1) {
    if (boundaries[2] == Boundary1D::MIRROR) {
      y0_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }

    if (boundaries[3] == Boundary1D::MIRROR) {
      ym_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }
  }
  if (dimensions > 2) {
    if (boundaries[4] == Boundary1D::MIRROR) {
      z0_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }

    if (boundaries[5] == Boundary1D::MIRROR) {
      zm_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }
  }

  Gaussian_noise* gaussian_noise;

  switch (dimensions) {
  case 1:
    boundary = new Boundary1D(Lat[0], boundaries[0], boundaries[1]);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(new Component(Lat[0], boundary, rho[i]));
      solver_component.push_back(new Component(Lat[0], boundary, rho[i]));
    }

    if (seed_specified == true) {
      gaussian_noise = new Gaussian_noise(boundary, D, M, mean, stddev, seed);
    } else {
      gaussian_noise = new Gaussian_noise(boundary, D, M, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(new Flux1D(Lat[0], gaussian_noise, D * dt, mask, component[i], component[j]));
        solver_flux.push_back(new Flux1D(Lat[0], gaussian_noise, D * dt, mask, solver_component[i], solver_component[j]));
      }
    }

    break;
  case 2:
    boundary = new Boundary2D(Lat[0], boundaries[0], boundaries[1], boundaries[2], boundaries[3]);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(new Component(Lat[0], boundary, rho[i]));
      solver_component.push_back(new Component(Lat[0], boundary, rho[i]));
    }

    if (seed_specified == true) {
      gaussian_noise = new Gaussian_noise(boundary, D, M, mean, stddev, seed);
    } else {
      gaussian_noise = new Gaussian_noise(boundary, D, M, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(new Flux2D(Lat[0], gaussian_noise, D * dt, mask, component[i], component[j]));
        solver_flux.push_back(new Flux2D(Lat[0], gaussian_noise, D * dt, mask, solver_component[i], solver_component[j]));
      }
    }

    break;
  case 3:
    boundary = new Boundary3D(Lat[0], boundaries[0], boundaries[1], boundaries[2], boundaries[3], boundaries[4], boundaries[5]);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(new Component(Lat[0], boundary, rho[i]));
      solver_component.push_back(new Component(Lat[0], boundary, rho[i]));
    }

    if (seed_specified == true) {
      gaussian_noise = new Gaussian_noise(boundary, D, M, mean, stddev, seed);
    } else {
      gaussian_noise = new Gaussian_noise(boundary, D, M, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(new Flux3D(Lat[0], gaussian_noise, D * dt, mask, component[i], component[j]));
        solver_flux.push_back(new Flux3D(Lat[0], gaussian_noise, D * dt, mask, solver_component[i], solver_component[j]));
      }
    }

    break;
  }

  norm_theta(component);
  norm_theta(solver_component);

  if (edge_detection)
    interface = new Interface(Lat[0], component);

  return 0;
}

void Mesodyn::explicit_start() {
  // Prepare callbacks
  auto explicit_solver_callback = bind(&Mesodyn::solve_explicit, this);
  auto loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  // start Crank-Nicolson using one explicit step
  New[0]->SolveMesodyn(loader_callback, explicit_solver_callback);

  // Update densities after first timestep
  int c = 0;
  for (size_t i = 0; i < component_no; ++i) {
    for (size_t j = 0; j < (component_no - 1) - i; ++j) {
      solver_component[i]->update_density(solver_flux[c]->J);
      ++c;
    }
  }

  c = 0;
  for (size_t j = 0; j < (component_no - 1); ++j) {
    for (size_t i = 1 + j; i < component_no; ++i) {
      solver_component[i]->update_density(solver_flux[c]->J, -1.0);
      ++c;
    }
  }

  // Load rho_k+1 into rho_k
  int i = 0;
  for (Component* all_components : component) {
    all_components->load_rho(solver_component[i]->rho);
    ++i;
  }

  // Load J_k+1 into J_k
  i = 0;
  for (Flux1D* all_fluxes : flux) {
    all_fluxes->J = solver_flux[i]->J;
    ++i;
  }
}

int Mesodyn::init_rho_equilibrate(vector<vector<Real>>& rho) {
  cout << "Equilibrating.." << endl;
  New[0]->Solve(true);
  for (size_t i = 0; i < component_no; ++i) {
    for (int z = 0; z < M; ++z)
      rho[i][z] = *(Seg[i]->H_phi + z);
  }
  return 0;
}

vector<string> Mesodyn::tokenize(string line, char delim) {
  istringstream stream{line};

  vector<string> tokens;
  string token;

  //Read lines as tokens (tab delimited)
  while (getline(stream, token, delim)) {
    tokens.push_back(token);
  }

  return tokens;
}

int Mesodyn::init_rho_fromvtk(vector<Real>& rho, string filename) {
  ifstream rho_input;

  rho_input.open(filename);

  if (!rho_input.is_open()) {
    cerr << "Error opening file! Is the filename correct? Is there a vtk for each component, ending in [component number].vtk, starting from 0?" << endl;
    throw ERROR_FILE_FORMAT;
  }

  string line;

  while (line.find("LOOKUP_TABLE default") == string::npos ) {
    getline(rho_input, line);
  }

  skip_bounds([this, &rho_input, &rho, &line](int x, int y, int z) mutable {
    getline(rho_input, line);
    *val_ptr(rho, x, y, z) = atof(line.c_str());
  });

  return 0;
}

int Mesodyn::init_rho_frompro(vector<vector<Real>>& rho, string filename) {

  ifstream rho_input;
  rho_input.open(filename);

  if (!rho_input.is_open()) {
    cerr << "Error opening file! Is the filename correct?" << endl;
    throw ERROR_FILE_FORMAT;
  }

  size_t i{1};

  //Discard rows that contain coordinates (will cause out of bounds read below)
  switch (dimensions) {
  case 1:
    i = 1;
    break;
  case 2:
    i = 2;
    break;
  case 3:
    i = 3;
    break;
  }

  //Find in which column the density profile starts
  //This depends on the fact that the first mon output is phi
  string line;
  getline(rho_input, line);

  if (line.find('\t') == string::npos) {
    cerr << "Wrong delimiter! Please use tabs.";
    throw ERROR_FILE_FORMAT;
  }

  vector<string> tokens = tokenize(line, '\t');
  int first_column{-1};
  int last_column{0};
  bool found = false;
  for (; i < tokens.size(); ++i) {
    vector<string> header_tokens = tokenize(tokens[i], ':');
    if (header_tokens.size() == 3) {
      if (header_tokens[0] == "mol" && header_tokens[2].size() > 3)
        if (header_tokens[2].substr(0, 4) == "phi-") {
          if (first_column == -1)
            first_column = i;
          found = true;
          last_column = i;
        }
    }
  }

  if (found != true) {
    cerr << "No headers in the format mol:[molecule]:phi-[monomer]." << endl;
    throw ERROR_FILE_FORMAT;
  }

  if (component_no != (size_t)(last_column - first_column + 1)) {
    cerr << "Not enough components detected in the headers, please adjust .pro file accordingly." << endl;
    throw ERROR_FILE_FORMAT;
  }

  int z{0};
  //Read lines one at a time
  while (getline(rho_input, line)) {

    vector<string> tokens = tokenize(line, '\t');
    //Read all densities into rho.
    for (size_t i = 0; i < component_no; ++i) {
      rho[i][z] = atof(tokens[first_column + i].c_str());
    }
    ++z;
  }

  if (z != M) { // +1 because z starts at 0 (would be M=1)
    cerr << "Input densities not of length M (" << z << "/" << M << ")" << endl;
    throw ERROR_FILE_FORMAT;
  }

  return 0;
}

int Mesodyn::norm_theta(vector<Component*>& component) {
  size_t solvent = (size_t)Sys[0]->solvent;

  Real sum_theta{0};
  int c{0};

  vector<int> solvent_mons;

  for (size_t i = 0; i < Mol.size(); ++i) {

    size_t mon_nr = Mol[i]->MolMonList.size();

    if (i != solvent) {
      Real theta = Mol[i]->theta;
      sum_theta += theta;
      for (size_t j = 0; j < mon_nr; ++j) {

        Real mon_theta{0};
        mon_theta = theta * Mol[i]->fraction(Mol[i]->MolMonList[j]);

        Real sum_of_elements{0};
        skip_bounds([this, &sum_of_elements, component, c](int x, int y, int z) mutable {
          sum_of_elements += val(component[c]->rho, x, y, z);
        });

        skip_bounds([this, &component, c, mon_theta, sum_of_elements](int x, int y, int z) mutable {
          *val_ptr(component[c]->rho, x, y, z) *= mon_theta / sum_of_elements;
        });

        ++c;
      }
    } else {
      for (size_t j = 0; j < mon_nr; ++j) {
        solvent_mons.push_back(c);
        ++c;
      }
    }
  }

  // Adjust solvent so that sum at position = 1

  // We now know the total density and can adjust the solvent accodingly to add up to 1.
  // Let's first find out how much there is to adjust.
  vector<Real> residuals(M);

  //Pool densities per position
  for (int i = 0; i < M; ++i) {
    for (size_t j = 0; j < component_no; ++j)
      residuals[i] += component[j]->rho[i];
  }

  // Calculate excesss / defecit
  for (Real& i : residuals) {
    i -= 1;
  }

  // If there's only one solvent mon, this problem is easy.
  if (solvent_mons.size() == 1) {
    skip_bounds([this, &component, residuals, solvent_mons](int x, int y, int z) mutable {
      *val_ptr(component[solvent_mons[0]]->rho, x, y, z) -= val(residuals, x, y, z);
    });
  } else {
    cerr << "Norming solvents with mutliple monomers is not supported! Please write your own script" << endl;
    throw ERROR_FILE_FORMAT;
  }

  return 0;
}

int Mesodyn::init_rho_homogeneous(vector<vector<Real>>& rho, vector<int>& mask) {
  size_t solvent = (size_t)Sys[0]->solvent; // Find which mol is the solvent

  int tMX = Lat[0]->MX;
  int tMY = Lat[0]->MY;
  int tMZ = Lat[0]->MZ;

  // The right hand side of the minus sign calculates the volume of the boundaries. I know, it's hideous.
  int volume = Sys[0]->volume - ((2 * dimensions - 4) * tMX * tMY + 2 * tMX * tMZ + 2 * tMY * tMZ + (-2 + 2 * dimensions) * (tMX + tMY + tMZ) + pow(2, dimensions));

  Real sum_theta{0};
  vector<Real> solvent_mons;
  int c{0};
  Real theta{0};
  Real solvent_mol{0};

  for (size_t i = 0; i < Mol.size(); ++i) {
    size_t mon_nr = Mol[i]->MolMonList.size();
    if (i != solvent) {
      theta = Mol[i]->theta;
      sum_theta += theta;
      for (size_t j = 0; j < mon_nr; ++j) {
        Real mon_theta{0};
        mon_theta = theta * Mol[i]->fraction(Mol[i]->MolMonList[j]);
        for (int z = 0; z < M; ++z) {
          for (size_t i = 0; i < component_no; ++i) {
            if (mask[z] == 0)
              rho[i][z] = 0;
            else
              rho[c][z] = mon_theta / volume;
          }
        }
        ++c;
      }
    } else {
      solvent_mol = i;
      for (size_t j = 0; j < mon_nr; ++j) {
        solvent_mons.push_back(c);
        ++c;
      }
    }
  }

  for (size_t i = 0; i < solvent_mons.size(); ++i) {
    Real mon_theta{0};
    mon_theta = (volume - sum_theta) * Mol[solvent_mol]->fraction(Mol[solvent_mol]->MolMonList[i]);
    for (int z = 0; z < M; ++z) {
      if (mask[z] == 0)
        rho[solvent][z] = 0;
      else
        rho[solvent_mons[i]][z] = mon_theta / volume;
    }
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

void Mesodyn::write_settings() {
  ostringstream settings_filename;

  settings_filename << filename.str() << "-settings.dat";

  mesodyn_output.open(settings_filename.str());

  mesodyn_output.precision(16);

  /* Output settings */

  ostringstream settings;

  settings << "Settings for mesodyn, read from file:,\n";
  settings << "Timesteps," << timesteps << ",\n";
  settings << "Time between saves," << timebetweensaves << ",\n";
  settings << "Diffusion constant," << D << ",\n";
  settings << "Seed for noise," << seed << ",\n";
  settings << "Mean for noise," << mean << ",\n";
  settings << "Stddev for noise," << stddev << ",\n";

  mesodyn_output << settings.str();
  mesodyn_output.close();
}

int Mesodyn::write_output() {
  //Implement below for more starts? (OutputList higher numbers?)

/* Output has changed. as you don't use output I have commented it out. New output has more arguments and we need to pass on the argument to mesodyn when you want to use it.
  for (size_t i = 0; i < In[0]->OutputList.size(); ++i) {
    Out.push_back(new Output(In, Lat, Seg, Mol, Sys, New, In[0]->OutputList[i], writes, timesteps / timebetweensaves));
    if (!Out[i]->CheckInput(1)) {
      cout << "input_error in output " << endl;
      return 0;
    }
  }*/

  PushOutput();
  New[0]->PushOutput();

  for (Output* all_output : Out) {
    all_output->WriteOutput(writes);
    delete all_output;
  }
  Out.clear();

  return 0;
}

void Mesodyn::write_density(vector<Component*>& component) {
  // Could also do this by overwriting the denisties in ... I think Mol..? right before it pushes the densities written by vtk in output.
  // Right now output does nothing with the vtk function when mesodyn is running.

  int component_count{1};
  for (Component* all_components : component) {

    ostringstream vtk_filename;
    vtk_filename << filename.str() << "-rho" << component_count << "-" << writes << ".vtk";

    mesodyn_output.open(vtk_filename.str());

    /* Generate headers for rho entries */

    ostringstream vtk;

    vtk << "# vtk DataFile Version 3.0 \n";
    vtk << "Mesodyn output \n";
    vtk << "ASCII\n";
    vtk << "DATASET STRUCTURED_GRID \n";
    vtk << "DIMENSIONS " << MX - 2 << " " << MY - 2 << " " << MZ - 2 << "\n";
    vtk << "POINTS " << (MX - 2) * (MY - 2) * (MZ - 2) << " int\n";
    skip_bounds([this, &vtk](int x, int y, int z) mutable {
      vtk << x << " " << y << " " << z << "\n";
    });

    vtk << "POINT_DATA " << (MX - 2) * (MY - 2) * (MZ - 2) << "\n";
    vtk << "SCALARS Component_" << component_count << " float\nLOOKUP_TABLE default \n";

    skip_bounds([this, &vtk, all_components](int x, int y, int z) mutable {
      vtk << val(all_components->rho, x, y, z) << "\n";
    });

    mesodyn_output << vtk.str();
    mesodyn_output.flush();

    mesodyn_output.close();
    ++component_count;
  }
  ++writes;
}

/******* INTERFACE ********/

Interface::Interface(Lattice* Lat, vector<Component*> components)
  : Lattice_Access(Lat), component(components), order_params(component[0]->rho.size()) {}

Interface::~Interface() {
  for (Component* all_components : component)
    delete all_components;
  component.clear();
}

int Interface::order_parameters(Component* A, Component* B) {
  if (order_params.size() != A->rho.size() || order_params.size() != B->rho.size() )
    throw ERROR_SIZE_INCOMPATIBLE;

  skip_bounds([this, A, B](int x, int y, int z) mutable {
    *val_ptr(order_params, x, y, z) = val(A->rho, x, y, z) * val(B->rho, x, y, z);
  });

  return 0;
}

int Interface::detect_edges(int threshold) {
  edges.clear();
  vector<Real> temp((MX-3)*(MY-3)*(MZ-3));
  temp = gaussian_blur(component[0]->rho);
  edges = sobel_edge_detector(threshold, temp);
  return 0;
}

int Interface::write_edges(string FILENAME) {
  ostringstream filename;
  filename << FILENAME;
  ofstream testfile;
  time_t rawtime;
  time(&rawtime);
  filename << "test" << rawtime << ".vtk";
  testfile.open( filename.str() );

  ostringstream vtk;

  vtk << "# vtk DataFile Version 3.0 \n";
  vtk << "Mesodyn output \n";
  vtk << "ASCII\n";
  vtk << "DATASET STRUCTURED_GRID \n";
  vtk << "DIMENSIONS " << MX-2 << " " << MY-2 << " " << MZ-2 << "\n";
  vtk << "POINTS " << (MX-2) * (MY-2) * (MZ-2) << " int\n";

  for (int x = 1; x < MX - 1; ++x)
    for (int y = 1; y < MY - 1; ++y)
      for (int z = 1 ; z < MZ - 1 ; ++z )
        vtk << x << " " << y << " " << z << "\n";

  vtk << "POINT_DATA " << (MX-2) * (MY-2) * (MZ-2) << "\n";
  vtk << "SCALARS Sobel float\nLOOKUP_TABLE default \n";

  for (Real& all_values : edges) {
    vtk << all_values << "\n";
  }

  testfile << vtk.str();

  testfile.flush();
  testfile.close();

  return 0;
}

vector<Real> Interface::sobel_edge_detector(Real tolerance, vector<Real>& rho) {
  vector<Real> result((MX - 2) * (MY - 2) * (MZ - 2));
  int threshold = tolerance;

  int i = 0;

  vector<int> Gy_minus = {-1, -3, -1, -3, -6, -3, -1, -3, -1};
  vector<int> Gy_mid = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<int> Gy_plus = {1, 3, 1, 3, 6, 3, 1, 3, 1};

  vector<int> Gx_minus = {-1, 0, 1, -3, 0, 3, -1, 0, 1};
  vector<int> Gx_mid = {-3, 0, 3, -6, 0, 6, -3, 0, 3};
  vector<int> Gx_plus = {-1, 0, 1, -3, 0, 3, -1, 0, 1};

  vector<int> Gz_minus = {1, 3, 1, 0, 0, 0, -1, -3, -1};
  vector<int> Gz_mid = {3, 6, 3, 0, 0, 0, -3, -6, -3};
  vector<int> Gz_plus = {1, 3, 1, 0, 0, 0, -1, -3, -1};

  for (int x = 0; x < MX - 2; ++x)
    for (int y = 0; y < MY - 2; ++y)
      for (int z = 0 ; z < MZ - 2 ; ++z ) {
        Real conv_x = convolution(Gx_minus, get_xy_plane(rho, x, y, z));
        conv_x += convolution(Gx_mid, get_xy_plane(rho, x, y, z+1));
        conv_x += convolution(Gx_plus, get_xy_plane(rho, x, y, z+2));

        Real conv_y = convolution(Gy_minus, get_xy_plane(rho, x, y, z));
        conv_y += convolution(Gy_mid, get_xy_plane(rho, x, y, z+1));
        conv_y += convolution(Gy_plus, get_xy_plane(rho, x, y, z+2));

        Real conv_z = convolution(Gz_minus, get_xz_plane(rho, x, y, z));
        conv_z += convolution(Gz_mid, get_xz_plane(rho, x, y+1, z));
        conv_z += convolution(Gz_plus, get_xz_plane(rho, x, y+2, z));

        result[i] = abs(conv_x) + abs(conv_y) + abs(conv_z);
        ++i;
    }

  //normalize between 0 and 255
  Real min = *min_element(result.begin(), result.end());
  Real max = *max_element(result.begin(), result.end());

  for (Real& all_elements : result)
    all_elements = (255 - 0) * ((all_elements - min) / (max - min)) + 0;

  //cut-off at threshold
  for (Real& all_elements : result)
    if (all_elements < threshold) {
      all_elements = 0;
    }

  return result;
}

vector<Real> Interface::gaussian_blur(vector<Real>& rho) {
  vector<Real> result( (MX-3)*(MY-3)*(MY-3) );
  vector<Real> G1 = {1, 2, 1, 2, 4, 2, 1, 2, 1};
  vector<Real> G2 = {1, 1, 1, 1, 2, 1, 1, 1, 1};
  vector<Real> G3 = {1, 1, 1, 1, 2, 1, 1, 1, 1};

  for (Real& all_values : G1)
    all_values = all_values * 1/16;

  for (Real& all_values : G2)
    all_values = all_values * 1/16;

  for (Real& all_values : G3)
    all_values = all_values * 1/16;

  int i{0};

  for (int x = 0; x < MX - 3; ++x)
    for (int y = 0; y < MY - 3; ++y)
      for (int z = 0 ; z < MZ - 3 ; ++z ) {
        Real conv_1 = convolution(G1, get_xy_plane(rho, x, y, z));
        Real conv_2 = convolution(G2, get_xy_plane(rho, x, y, z+1));
        Real conv_3 = convolution(G3, get_xy_plane(rho, x, y, z+2));
        result[i] = (conv_1 + conv_2 + conv_3);
        ++i;
    }
  return result;
}

Real Interface::convolution(vector<int> kernel, vector<Real> pixel) {
  if (kernel.size() != pixel.size()) {
    cerr << "Convolution: pixel and kernel not of equal size!" << endl;
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  Real accumulator = 0;

  for (size_t i = 0 ; i < kernel.size() ; ++i) {
    accumulator += kernel[i] * pixel[i];
  }

  return accumulator;
}

//TODO: template this
Real Interface::convolution(vector<Real> kernel, vector<Real> pixel) {
  if (kernel.size() != pixel.size()) {
    cerr << "Convolution: pixel and kernel not of equal size!" << endl;
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  Real accumulator = 0;

  for (size_t i = 0 ; i < kernel.size() ; ++i) {
    accumulator += kernel[i] * pixel[i];
  }

  return accumulator;
}

vector<Real> Interface::get_xy_plane(vector<Real>& rho, int x, int y, int z, int size) {
  vector<Real> pixel(size*size);

  int i = 0;

  for (int vertical = 0 ; vertical < size ; ++vertical)
    for (int horizontal = 0 ; horizontal < size ; ++horizontal) {
        pixel[i] = val(rho, x+horizontal, y+vertical, z);
        ++i;
    }

  return pixel;
}

vector<Real> Interface::get_xz_plane(vector<Real>& rho, int x, int y, int z, int size) {
  vector<Real> pixel(size*size);

  int i = 0;

  for (int horizontal = 0 ; horizontal < size ; ++horizontal)
    for (int depth = 0 ; depth < size ; ++depth) {
        pixel[i] = val(rho, x+horizontal, y, z+depth);
        ++i;
    }

  return pixel;
}

/******* FLUX: TOOLS FOR CALCULATING FLUXES BETWEEN 1 PAIR OF COMPONENTS, HANDLING OF SOLIDS *********/

Flux1D::Flux1D(Lattice* Lat, Gaussian_noise* gaussian, Real D, vector<int>& mask, Component* A, Component* B)
    : Lattice_Access(Lat), J_plus(M), J_minus(M), J(M), gaussian{gaussian}, A{A}, B{B}, L(M), mu(M), D{D}, JX{Lat->JX} {
  Flux1D::mask(mask);
}

Flux2D::Flux2D(Lattice* Lat, Gaussian_noise* gaussian, Real D, vector<int>& mask, Component* A, Component* B)
    : Flux1D(Lat, gaussian, D, mask, A, B), JY{Lat->JY} {
  Flux2D::mask(mask);
}

Flux3D::Flux3D(Lattice* Lat, Gaussian_noise* gaussian, Real D, vector<int>& mask, Component* A, Component* B)
    : Flux2D(Lat, gaussian, D, mask, A, B), JZ{Lat->JZ} {
  Flux3D::mask(mask);
}

Flux1D::~Flux1D() {
}

Flux2D::~Flux2D() {
}

Flux3D::~Flux3D() {
}

Real Flux1D::J_at(int x, int y, int z) {
  return val(J, x, y, z);
}

Real Flux1D::L_at(int x, int y, int z) {
  return val(L, x, y, z);
}

Real Flux1D::mu_at(int x, int y, int z) {
  return val(mu, x, y, z);
}

int Flux1D::mask(vector<int>& mask_in) {
  if ((int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  // No fluxes will ever be calculated going from the boundary into the system
  skip_bounds([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y, z) == 1) {
      if (val(mask_in, x + 1, y, z) == 1) {
        Mask_plus_x.push_back(index(x, y, z));
      }
      if (val(mask_in, x - 1, y, z) == 1) {
        Mask_minus_x.push_back(index(x, y, z));
      }
    }
  });

  //for the x-boundary:
  x0_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x+1, y, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_plus_x.push_back(index(x, y, z));
  });
  xm_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x-1, y, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_minus_x.push_back(index(x, y, z));
  });

  return 0;
}

int Flux2D::mask(vector<int>& mask_in) {
  if ((int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  // No fluxes will ever be calculated going from the boundary into the system
  skip_bounds([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y, z) == 1) {
      if (val(mask_in, x, y + 1, z) == 1)
        Mask_plus_y.push_back(index(x, y, z));
      if (val(mask_in, x, y - 1, z) == 1)
        Mask_minus_y.push_back(index(x, y, z));
    }
  });

  //for the y-boundary:
  y0_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y+1, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_plus_y.push_back(index(x, y, z));
  });
  ym_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y-1, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_minus_y.push_back(index(x, y, z));
  });

  return 0;
}

int Flux3D::mask(vector<int>& mask_in) {
  if ((int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  // No fluxes will ever be calculated going from the boundary into the system, messing up periodic boundaries
  skip_bounds([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y, z) == 1) {
      if (val(mask_in, x, y, z + 1) == 1) {
        Mask_plus_z.push_back(index(x, y, z));
      }
      if (val(mask_in, x, y, z - 1) == 1) {
        Mask_minus_z.push_back(index(x, y, z));
      }
    }
  });

  //for the z-boundary:
  z0_boundary([this, mask_in](int x, int y, int z) mutable {
      if (val(mask_in, x, y, z+1) == 1 && val(mask_in, x, y, z) == 1)
        Mask_plus_z.push_back(index(x, y, z));
  });
  zm_boundary([this, mask_in](int x, int y, int z) mutable {
      if (val(mask_in, x, y, z-1) == 1 && val(mask_in, x, y, z) == 1)
        Mask_minus_z.push_back(index(x, y, z));
  });

  return 0;
}

int Flux1D::langevin_flux() {

  //Zero (with bounds checking) vector J before use
  for (Real& i : J) {
    i = 0;
  }

  if (A->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  onsager_coefficient(A->rho, B->rho);
  potential_difference(A->alpha, B->alpha);
  langevin_flux(Mask_plus_x, Mask_minus_x, JX);

  return 0;
}

int Flux2D::langevin_flux() {
  Flux1D::langevin_flux();

  Flux1D::langevin_flux(Mask_plus_y, Mask_minus_y, JY);

  return 0;
}

int Flux3D::langevin_flux() {
  Flux2D::langevin_flux();

  Flux1D::langevin_flux(Mask_plus_z, Mask_minus_z, JZ);

  return 0;
}

int Flux1D::onsager_coefficient(vector<Real>& A, vector<Real>& B) {
  //TODO: maybe do this in propagator style inline / per J calculation to preserve memory

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  transform(A.begin(), A.end(), B.begin(), L.begin(), [](Real A, Real B) { return A * B; });
  return 0;
}

int Flux1D::potential_difference(vector<Real>& A, vector<Real>& B) {
  //TODO: maybe do this in propagator style inline / per J calculation to preserve memory

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  transform(A.begin(), A.end(), B.begin(), mu.begin(), [](Real A, Real B) { return A - B; });

  return 0;
}

int Flux1D::langevin_flux(vector<int>& mask_plus, vector<int>& mask_minus, int jump) {

  for (Real& i : J_minus) {
    i = 0;
  }
  for (Real& i : J_plus) {
    i = 0;
  }

  for (int& z : mask_plus) {
    J_plus[z] = -D * ((L[z] + L[z + jump]) * (mu[z + jump] - mu[z] + gaussian->noise[z]));
  }

  for (int& z : mask_minus) {
    J_minus[z] = -J_plus[z - jump]; // = -D * ((L[z - jump] + L[z]) * (mu[z - jump] - mu[z] - gaussian->noise[z]));
    // We have to do it this way because otherwise the noise will be trouble
  }

  transform(J_plus.begin(), J_plus.end(), J.begin(), J.begin(), [](Real A, Real B) { return A + B; });
  transform(J_minus.begin(), J_minus.end(), J.begin(), J.begin(), [](Real A, Real B) { return A + B; });

  return 0;
}

/****************** COMPONENT: DENSITY PROFILE STORAGE AND UPDATING, BOUNDARY CONDITIONS ********************/

/******* Constructors *******/

Component::Component(Lattice* Lat, Boundary1D* boundary, vector<Real>& rho)
    : Lattice_Access(Lat), rho{rho}, alpha(M), boundary(boundary) {
  //This check is implemented multiple times throughout mesodyn because rho and alpha are public.
  if (rho.size() != alpha.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  update_boundaries();
}

Component::~Component() {
}

/******* Interface *******/

Real Component::rho_at(int x, int y, int z) {
  return val(rho, x, y, z);
}

Real Component::alpha_at(int x, int y, int z) {
  return val(alpha, x, y, z);
}

int Component::update_density(vector<Real>& J, int sign) {
  //Explicit update

  if (J.size() != rho.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  skip_bounds([this, J, sign](int x, int y, int z) mutable {
    *val_ptr(rho, x, y, z) += val(J, x, y, z) * sign;
  });

  return 0;
}

int Component::update_density(vector<Real>& rho_old, vector<Real>& J1, vector<Real>& J2, Real ratio, int sign) {
  //Implicit update
  if (J1.size() != rho.size() || J1.size() != J2.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  skip_bounds([this, J1, ratio, sign, rho_old](int x, int y, int z) mutable {
    *val_ptr(rho, x, y, z) = (val(rho_old, x, y, z) + ratio * sign * val(J1, x, y, z));
  });

  skip_bounds([this, J2, ratio, sign](int x, int y, int z) mutable {
    *val_ptr(rho, x, y, z) = val(rho, x, y, z) + (1 - ratio) * sign * val(J2, x, y, z);
  });

  return 0;
}

int Component::load_alpha(vector<Real>& alpha) {
  Component::alpha = alpha;
  return 0;
}

int Component::load_rho(vector<Real>& rho) {
  Component::rho = rho;
  return 0;
}

int Component::update_boundaries() {
  boundary->update_boundaries(alpha);
  boundary->update_boundaries(rho);
  return 0;
}

Real Component::theta() {
  Real sum{0};
  skip_bounds([this, &sum](int x, int y, int z) mutable {
    sum += val(rho, x, y, z);
  });
  return sum;
}

/****************** Lattice_Access: AN INTERFACE FOR LATTICE ********************/

Lattice_Access::Lattice_Access(Lattice* Lat)
    : dimensions{Lat->gradients}, JX{Lat->JX}, JY{Lat->JY}, JZ{Lat->JZ}, M{Lat->M}, MX{2 + Lat->MX}, MY{setMY(Lat)}, MZ{setMZ(Lat)} {
  // If this is not true, NOTHING will work. So this check is aboslutely necessary.
  // If, at some point, the above happens to be the case every class in this module will probably need rewriting.
  assert(
      (MX > 0 && MY == 0 && MZ == 0) ||
      (MX > 0 && MY > 0 && MZ == 0) ||
      (MX > 0 && MY > 0 && MZ > 0));
}

Lattice_Access::~Lattice_Access() {
}

inline void Lattice_Access::skip_bounds(function<void(int, int, int)> function) {
  int x{1};
  int y{1};
  int z{1};
  do {
    y = 1;
    do {
      x = 1;
      do {
        function(x, y, z);
        ++x;
      } while (x < MX - 1);
      ++y;
    } while (y < MY - 1);
    ++z;
  } while (z < MZ - 1);
}

inline void Lattice_Access::bounds(function<void(int, int, int)> function) {

  x0_boundary(function);
  xm_boundary(function);

  if (dimensions > 1) {
    y0_boundary(function);
    ym_boundary(function);
  }

  if (dimensions > 2) {
    z0_boundary(function);
    zm_boundary(function);
  }
}

inline void Lattice_Access::x0_boundary(function<void(int, int, int)> function) {
  int x = 0;
  int y = 1;
  int z = 1;
  do {
    y = 1;
    do {
      function(x, y, z);
      ++y;
    } while (y < MY - 1);
    ++z;
  } while (z < MZ - 1);
}

inline void Lattice_Access::xm_boundary(function<void(int, int, int)> function) {
  int x = MX - 1, y = 1, z = 1;
  do {
    y = 1;
    do {
      function(x, y, z);
      ++y;
    } while (y < MY - 1);
    ++z;
  } while (z < MZ - 1);
}

inline void Lattice_Access::y0_boundary(function<void(int, int, int)> function) {
  int x = 1, y = 0, z = 1;
  do {
    x = 1;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX - 1);
    ++z;
  } while (z < MZ - 1);
}

inline void Lattice_Access::ym_boundary(function<void(int, int, int)> function) {
  int x = 1, y = MY - 1, z = 1;
  do {
    x = 1;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX - 1);
    ++z;
  } while (z < MZ - 1);
}

inline void Lattice_Access::z0_boundary(function<void(int, int, int)> function) {
  int x = 1, y = 1, z = 0;
  do {
    x = 1;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX - 1);
    ++y;
  } while (y < MY - 1);
}

inline void Lattice_Access::zm_boundary(function<void(int, int, int)> function) {
  int x = 1, y = 1, z = MZ - 1;
  do {
    x = 1;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX - 1);
    ++y;
  } while (y < MY - 1);
}

inline Real Lattice_Access::val(vector<Real>& v, int x, int y, int z) {
  return v[x * JX + y * JY + z * JZ];
}

inline int Lattice_Access::val(vector<int>& v, int x, int y, int z) {
  return v[x * JX + y * JY + z * JZ];
}

inline Real* Lattice_Access::val_ptr(vector<Real>& v, int x, int y, int z) {
  return &v[x * JX + y * JY + z * JZ];
}

inline int* Lattice_Access::val_ptr(vector<int>& v, int x, int y, int z) {
  return &v[x * JX + y * JY + z * JZ];
}

int Lattice_Access::setMY(Lattice* Lat) {
  //used by constructor
  if (dimensions < 2)
    return 0;
  else
    return Lat->MY + 2;
}

int Lattice_Access::setMZ(Lattice* Lat) {
  //used by constructor
  if (dimensions < 3)
    return 0;
  else
    return Lat->MZ + 2;
}

inline int Lattice_Access::index(int x, int y, int z) {
  return x * JX + y * JY + z * JZ;
}

inline vector<int> Lattice_Access::coordinate(int n) {
  int x = 0, y = 0, z = 0;
  int mod = 0;
  x = n / JX;
  mod = n % JX;

  if (dimensions > 1) {
    y = mod / JY;
    mod = mod % JY;
  }

  if (dimensions > 2) {
    z = mod / JZ;
  }
  return {x, y, z};
}

/******* Boundary conditions *******/

Boundary1D::Boundary1D(Lattice* Lat, boundary x0, boundary xm)
    : Lattice_Access(Lat) {
  set_x_boundaries(x0, xm);
}

Boundary2D::Boundary2D(Lattice* Lat, boundary x0, boundary xm, boundary y0, boundary ym)
    : Boundary1D(Lat, x0, xm) {
  set_y_boundaries(y0, ym);
}

Boundary3D::Boundary3D(Lattice* Lat, boundary x0, boundary xm, boundary y0, boundary ym, boundary z0, boundary zm)
    : Boundary2D(Lat, x0, xm, y0, ym) {
  set_z_boundaries(z0, zm);
}

Boundary1D::~Boundary1D() {}

Boundary2D::~Boundary2D() {}

Boundary3D::~Boundary3D() {}

int Boundary1D::update_boundaries(vector<Real>& target) {
  bX0(target);
  bXm(target);
  return 0;
}

int Boundary2D::update_boundaries(vector<Real>& target) {
  Boundary1D::update_boundaries(target);
  bY0(target);
  bYm(target);
  return 0;
}

int Boundary3D::update_boundaries(vector<Real>& target) {
  Boundary2D::update_boundaries(target);
  bZ0(target);
  bZm(target);
  return 0;
}

int Boundary1D::set_x_boundaries(boundary x0, boundary xm, Real bulk) {
  using namespace std::placeholders;

  switch (x0) {
  case MIRROR:
    bX0 = bind(&Boundary1D::bX0Mirror, this, _1);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bX0 = bind(&Boundary1D::bXPeriodic, this, _1);
    break;
  }

  switch (xm) {
  case MIRROR:
    bXm = bind(&Boundary1D::bXmMirror, this, _1);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    //TODO: fix this double periodic (including ones below)
    bXm = bind(&Boundary1D::bXPeriodic, this, _1);
    break;
  }

  return 0;
}

int Boundary2D::set_y_boundaries(boundary y0, boundary ym, Real bulk) {
  using namespace std::placeholders;

  switch (y0) {
  case MIRROR:
    bY0 = bind(&Boundary2D::bY0Mirror, this, _1);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bY0 = bind(&Boundary2D::bYPeriodic, this, _1);
    break;
  }

  switch (ym) {
  case MIRROR:
    bYm = bind(&Boundary2D::bYmMirror, this, _1);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bYm = bind(&Boundary2D::bYPeriodic, this, _1);
    break;
  }

  return 0;
}

int Boundary3D::set_z_boundaries(boundary z0, boundary zm, Real bulk) {
  using namespace std::placeholders;

  switch (z0) {
  case MIRROR:
    bZ0 = bind(&Boundary3D::bZ0Mirror, this, _1);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bZ0 = bind(&Boundary3D::bZPeriodic, this, _1);
    break;
  }

  switch (zm) {
  case MIRROR:
    bZm = bind(&Boundary3D::bZmMirror, this, _1);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bZm = bind(&Boundary3D::bZPeriodic, this, _1);
    break;
  }

  return 0;
}

void Boundary1D::bX0Mirror(vector<Real>& target) {
  x0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x + 1, y, z); //start
  });
}

void Boundary1D::bXmMirror(vector<Real>& target) {
  xm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x - 1, y, z); //end
  });
}

void Boundary1D::bXPeriodic(vector<Real>& target) {
  x0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, MX - 2, y, z); //start
  });
  xm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, 1, y, z); //end
  });
}

void Boundary2D::bY0Mirror(vector<Real>& target) {
  y0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y + 1, z); //start
  });
}

void Boundary2D::bYmMirror(vector<Real>& target) {
  ym_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y - 1, z); //end
  });
}

void Boundary2D::bYPeriodic(vector<Real>& target) {
  y0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, MY - 2, z); //start
  });
  ym_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, 1, z); //end
  });
}
void Boundary3D::bZ0Mirror(vector<Real>& target) {
  z0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, z + 1); //start
  });
}

void Boundary3D::bZmMirror(vector<Real>& target) {
  zm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, z - 1); //end
  });
}

void Boundary3D::bZPeriodic(vector<Real>& target) {
  z0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, MZ - 2); //start
  });
  zm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, 1); //end
  });
}

/******* GAUSSIAN_NOISE: GENERATES WHITE NOISE FOR FLUXES ********/

Gaussian_noise::Gaussian_noise(Boundary1D* boundary, Real D, int size, Real mean, Real stddev) : noise(size), prng{std::random_device{}()}, dist(mean, stddev), boundary{boundary} {}

Gaussian_noise::Gaussian_noise(Boundary1D* boundary, Real D, int size, Real mean, Real stddev, size_t seed) : noise(size), prng(seed), dist(mean, stddev), boundary{boundary} {}

int Gaussian_noise::generate() {
  for (Real& value : noise)
    value = dist(prng);

  boundary->update_boundaries(noise);

  return 0;
}

int Gaussian_noise::add_noise(vector<Real>& target) {
  transform(noise.begin(), noise.end(), target.begin(), target.begin(), [](Real A, Real B) { return A + B; });
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

/***** IO *****/

void Mesodyn::PutParameter(string new_param) { KEYS.push_back(new_param); }
string Mesodyn::GetValue(string parameter) {
  int i = 0;
  int length = PARAMETERS.size();
  while (i < length) {
    if (parameter == PARAMETERS[i]) {
      return VALUES[i];
    }
    i++;
  }
  return "";
}
void Mesodyn::push(string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Mesodyn::push(string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Mesodyn::push(string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}
void Mesodyn::push(string s, string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
void Mesodyn::PushOutput() {
  strings.clear();
  strings_value.clear();
  bools.clear();
  bools_value.clear();
  Reals.clear();
  Reals_value.clear();
  ints.clear();
  ints_value.clear();
}
int Mesodyn::GetValue(string prop, int& int_result, Real& Real_result,
                      string& string_result) {
  int i = 0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}
