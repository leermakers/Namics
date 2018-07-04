#include "mesodyn.h"
#include "newton.h"

/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018 */
/* Most of the physics in this module is based on the work of Fraaije et al. in the 1990s  */

Mesodyn::Mesodyn(vector<Input *> In_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<Molecule *> Mol_, vector<System *> Sys_, vector<Solve_scf*> New_, string name_)
  :
    Lattice_Access(Lat_[0]),
    name{name_}, In{In_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sys{Sys_}, New{New_},
    D{0.01}, mean{0}, stddev{1*D}, seed{1}, timesteps{100}, timebetweensaves{1},
    componentNo{(int)Sys[0]->SysMolMonList.size()},
    boundary(0), component(0), flux(0),
    writes{0}
{
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stddev");

  if (debug)
    cout << "Mesodyn initialized." << endl;
}

Mesodyn::~Mesodyn() {
  for (unsigned int i = 0; i < flux.size(); ++i)
    delete flux[i];
  flux.clear();
  for (unsigned int i = 0; i < component.size(); ++i)
    delete component[i];
  component.clear();
  delete boundary[0];
  boundary.clear();
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

    if (GetValue("diffusionconstant").size() > 0) {
      D = In[0]->Get_Real(GetValue("diffusionconstant"),D);
    }
    if (debug)
      cout << "Diffusion const is " << D << endl;

    if (GetValue("seed").size() > 0) {
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
  }

  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  if (debug)
    cout << "mesodyn in Mesodyn" << endl;

  initial_conditions(); // Initialize densities by running the classical method once.
  set_filename();
  write_settings();

  // write initial conditions
  write_density(solver_component);

  if (debug)
    cout << "Mesodyn is all set, starting calculations.." << endl;

  vector<Real> rho;

  for (int t = 1; t < timesteps; t++) {
    cout << "MESODYN: t = " << t << endl;

    gaussian_noise->generate();

    New[0]->SolveMesodyn(
        [this](vector<Real>& alpha, int i) {
          if (i < componentNo) solver_component[i]->load_alpha(alpha);
        },
        [this, &rho] () {
          for (Component* all_components : solver_component) all_components->update_boundaries();

          for (Flux1D* all_fluxes : solver_flux) all_fluxes->langevin_flux();

          int c = 0;
          for (int i = 0; i < componentNo; ++i) {
            for (int j = 0; j < (componentNo - 1) - i; ++j) {
              solver_component[i]->update_density(component[i]->rho, flux[c]->J, solver_flux[c]->J);
              ++c;
            }
          }

          c = 0;
          for (int j = 0; j < (componentNo - 1); ++j) {
            for (int i = 1 + j; i < componentNo; ++i) {
              solver_component[i]->update_density(component[i]->rho, flux[c]->J, solver_flux[c]->J, -1.0);
              ++c;
            }
          }

          rho.clear();
          for (Component* all_components : solver_component) {
            copy(all_components->rho.begin(), all_components->rho.end(), back_inserter(rho));
          }

          return &rho[0];
        }
      );

    int i = 0;
    for(Component* all_components : component) {
      all_components->load_rho(solver_component[i]->rho);
      ++i;
    }

    if (t % timebetweensaves == 0)
      write_density(solver_component);

    i = 0;
    for(Flux1D* all_fluxes : flux) {
      all_fluxes->J = solver_flux[i]->J;
      ++i;
    }

  } // time loop

  write_density(solver_component);
  return true;
}

int Mesodyn::initial_conditions() {

  //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (int i = 0; i < (int)Seg.size(); ++i) {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }

  vector<vector<Real>> rho(componentNo, vector<Real>(M));

  Sys[0]->PrepareForCalculations(); // to get the correct KSAM.

  //TODO: once KSAM becomes a vector, this loop won't be nessecary anymore, just pass KSAM to the flux constructor.
  vector<int> mask(M);
  for (int i = 0; i < M; ++i )
  {
    mask[i] = *(Sys[0]->KSAM+i);
  }

  init_rho(rho, mask);

  vector<Boundary1D::boundary> boundaries;
  for (string& boundary : Lat[0]->BC) {
    if (boundary == "mirror") {
      boundaries.push_back(Boundary1D::MIRROR);
    }
    if (boundary == "periodic") {
      boundaries.push_back(Boundary1D::PERIODIC);
    }
    if (boundary == "bulk") {
      boundaries.push_back(Boundary1D::BULK);
    }
  }

  boundary.push_back(new Boundary1D(Lat[0], boundaries[0], boundaries[1]));
  gaussian_noise = new Gaussian_noise(boundary[0], D, M, mean, stddev);

  switch (dimensions) {
  case 1:

    for (int i = 0; i < componentNo; ++i) {
      component.push_back(new Component(Lat[0], boundary[i], rho[i]));
      solver_component.push_back(new Component(Lat[0], boundary[0], rho[i]));
    }

    for (int i = 0; i < componentNo - 1; ++i) {
      for (int j = i + 1; j < componentNo; ++j) {
        flux.push_back(new Flux1D(Lat[0], gaussian_noise, D, mask, component[i], component[j]));
        solver_flux.push_back(new Flux1D(Lat[0], gaussian_noise, D, mask, solver_component[i], solver_component[j]));
      }
    }

    break;
  case 2:
    for (int i = 0; i < componentNo; ++i) {
      boundary.push_back(new Boundary2D(Lat[0], boundaries[0], boundaries[1], boundaries[2], boundaries[3]));

      component.push_back(new Component(Lat[0], boundary[i], rho[i]));
      solver_component.push_back(new Component(Lat[0], boundary[0], rho[i]));
    }

    for (int i = 0; i < componentNo - 1; ++i) {
      for (int j = i + 1; j < componentNo; ++j) {
        flux.push_back(new Flux2D(Lat[0], gaussian_noise, D, mask, component[i], component[j]));
        solver_flux.push_back(new Flux2D(Lat[0], gaussian_noise, D, mask, solver_component[i], solver_component[j]));
      }
    }

    break;
  case 3:
    for (int i = 0; i < componentNo; ++i) {
      boundary.push_back(new Boundary3D(Lat[0], boundaries[0], boundaries[1], boundaries[2], boundaries[3], boundaries[4], boundaries[5]));

      component.push_back(new Component(Lat[0], boundary[i], rho[i]));
      solver_component.push_back(new Component(Lat[0], boundary[0], rho[i]));
    }

    for (int i = 0; i < componentNo - 1; ++i) {
      for (int j = i + 1; j < componentNo; ++j) {
        flux.push_back(new Flux3D(Lat[0], gaussian_noise, D, mask, component[i], component[j]));
        solver_flux.push_back(new Flux3D(Lat[0], gaussian_noise, D, mask, solver_component[i], solver_component[j]));
      }
    }

    break;
  }

  return 0;
}

int Mesodyn::init_rho(vector<vector<Real>>& rho, vector<int>& mask) {
  int solvent = Sys[0]->solvent; // Find which componentNo is the solvent
//  int volume = Sys[0]->volume - (pow(M, dimensions) - pow((M - 2), dimensions));
  int tMX = Lat[0]->MX;
  int tMY = Lat[0]->MY;
  int tMZ = Lat[0]->MZ;

  // The right hand side of the minus sign calculates the volume of the boundaries. I know, it's hideous.
  int volume = Sys[0]->volume - ( (2*dimensions-4)*tMX*tMY+2*tMX*tMZ+2*tMY*tMZ+(-2+2*dimensions)*(tMX+tMY+tMZ)+pow(2,dimensions));

  Real sum_theta{0};
  Real theta{0};

  for (int i = 0; i < componentNo; ++i) {
    if (i != solvent) {
      theta = Mol[i]->theta;
      sum_theta += theta;
      for (int z = 0; z < M; ++z) {
        for (int i = 0; i < componentNo; ++i) {
          if (mask[z] == 0) rho[i][z] = 0;
          else rho[i][z] = theta / volume;
        }
      }
    }
  }

  for (int z = 0; z < M; ++z) {
    if (mask[z] == 0) rho[solvent][z] = 0;
    else rho[solvent][z] = (volume - sum_theta) / volume;
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

  settings_filename << filename.str()  << "-settings.dat";

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

void Mesodyn::write_density(vector<Component*>& component) {
  int component_count{1};
  for (Component* all_components : component) {

    ostringstream vtk_filename;
    vtk_filename << filename.str() << "-rho" << component_count << "-" << writes << ".vtk";

    mesodyn_output.open(vtk_filename.str());

    /* Generate headers for rho entries */

    ostringstream vtk;

    //int permutations = combinations(componentNo, 2);

    vtk << "# vtk DataFile Version 7.0 \nvtk output \nASCII \nDATASET STRUCTURED_POINTS \nDIMENSIONS " << MX << " " << MY << " " << MZ << "\n";
    vtk << "SPACING 1 1 1 \nORIGIN 0 0 0 \nPOINT_DATA " << MX * MY * MZ << "\n";
    vtk << "SCALARS Box_profile float\nLOOKUP_TABLE default \n";

    for (Real& value : all_components->rho)
      vtk << value << " \n";

    mesodyn_output << vtk.str();
    mesodyn_output.flush();

    mesodyn_output.close();
    ++component_count;
  }
  ++writes;
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

int Flux1D::langevin_flux() {

  //Zero (with bounds checking) vector J before use
  for (Real& i : J) {
    i = 0;
  }
  for (Real& i : J_minus) {
    i = 0;
  }
  for (Real& i : J_plus) {
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

Real Flux1D::J_at(int x, int y, int z) {
  return val(J, x, y, z);
}

Real Flux1D::L_at(int x, int y, int z) {
  return val(L, x, y, z);
}

Real Flux1D::mu_at(int x, int y, int z) {
  return val(mu, x, y, z);
}

inline int val(vector<int>&, int, int, int);
int Flux1D::mask(vector<int>& mask_in) {
  if ( (int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  int x,y,z;

  z = 1;
  do {
    y = 1;
    do{
      x = 1;
      do{
        if (val(mask_in,x,y,z) == 1) {
          if ( val(mask_in,x+1,y,z) == 1) {
            Mask_plus_x.push_back( index(x,y,z) );
          }
          if ( val(mask_in,x-1,y,z) == 1) {
            Mask_minus_x.push_back( index(x,y,z) );
          }
        }
        ++x;
      } while (x < MX - 1);
      ++y;
    } while (y < MY - 1);
    ++z;
  } while (z < MZ - 1);

  return 0;
}

int Flux2D::mask(vector<int>& mask_in) {
  if ( (int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  Flux1D::mask(mask_in);
  int x,y,z;

  z = 1;
  do {
    y = 1;
    do{
      x = 1;
      do{
        if (val(mask_in,x,y,z) == 1) {
          if ( val(mask_in,x,y+1,z) == 1) {
            Mask_plus_y.push_back( index(x,y,z) );
          }
          if (val(mask_in,x,y-1,z) == 1) {
            Mask_minus_y.push_back( index(x,y,z) );
          }
        }
        ++x;
      } while (x < MX - 1);
      ++y;
    } while (y < MY - 1);
    ++z;
  } while (z < MZ - 1);

  return 0;
}


int Flux3D::mask(vector<int>& mask_in) {
  if ( (int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  Flux2D::mask(mask_in);
  int x,y,z;

  z = 1;
  do {
    y = 1;
    do{
      x = 1;
      do{
        if (val(mask_in,x,y,z) == 1) {
          if ( val(mask_in,x,y,z+1) == 1) {
            Mask_plus_z.push_back( index(x,y,z) );
          }
          if ( val(mask_in,x,y,z-1) == 1) {
            Mask_minus_z.push_back( index(x,y,z) );
          }
        }
        ++x;
      } while (x < MX - 1);
      ++y;
    } while (y < MY - 1);
    ++z;
  } while (z < MZ - 1);

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

  transform(A.begin(), A.end(), B.begin(), mu.begin(), [](Real A, Real B) { return A - B ; });
  gaussian->add_noise(mu);

  return 0;
}

int Flux1D::langevin_flux(vector<int>& mask_plus, vector<int>& mask_minus, int jump) {

  for (int& z: mask_plus) {
    J_plus[z] += -D * ((L[z] + L[z + jump]) * (mu[z + jump] - mu[z]));
  }
  for (int& z: mask_minus) {
    J_minus[z] += -D * ((L[z - jump] + L[z]) * (mu[z - jump] - mu[z])); // = -J_plus[z-jump];
  }

  transform(J_plus.begin(), J_plus.end(), J.begin(), J.begin(), [](Real A, Real B) { return A + B; });
  transform(J_minus.begin(), J_minus.end(), J.begin(), J.begin(), [](Real A, Real B) { return A + B; });
  return 0;
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

inline Real Lattice_Access::val(vector<Real>& v, int x, int y, int z) {
  return v[x * JX + y * JY + z * JZ];
}

inline int Lattice_Access::val(vector<int>& v, int x, int y, int z) {
  return v[x * JX + y * JY + z * JZ];
}

inline Real* Lattice_Access::valPtr(vector<Real>& v, int x, int y, int z) {
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

  int i = 0;

  for (Real& Flux : J) {
    rho[i] += sign * Flux;
    ++i;
  }

  return 0;
}

int Component::update_density(vector<Real>& rho_old, vector<Real>& J1, vector<Real>& J2, int sign) {
  //Implicit update
  if (J1.size() != rho.size() || J1.size() != J2.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  int i = 0;

  for (Real& Flux : J1) {
    rho[i] = rho_old[i] + 0.5 * sign * Flux;
    ++i;
  }

  i = 0;

  for (Real& Flux : J2) {
    rho[i] += 0.5 * sign * Flux;
    ++i;
  }

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
    bX0 = bind(&Boundary1D::bX0Mirror, this, _1, MY, MZ);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bX0 = bind(&Boundary1D::bXPeriodic, this, _1, MY, MZ, MX);
    break;
  case BULK:
    bX0 = bind(&Boundary1D::bX0Bulk, this, _1, MY, MZ, bulk);
    break;
  }

  switch (xm) {
  case MIRROR:
    bXm = bind(&Boundary1D::bXmMirror, this, _1, MY, MZ, MX);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    //TODO: fix this double periodic (including ones below)
    bXm = bind(&Boundary1D::bXPeriodic, this, _1, MY, MZ, MX);
    break;
  case BULK:
    bXm = bind(&Boundary1D::bXmBulk, this, _1, MY, MZ, MX, bulk);
    break;
  }

  return 0;
}

int Boundary2D::set_y_boundaries(boundary y0, boundary ym, Real bulk) {
  using namespace std::placeholders;

  switch (y0) {
  case MIRROR:
    bY0 = bind(&Boundary2D::bY0Mirror, this, _1, MX, MZ);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bY0 = bind(&Boundary2D::bYPeriodic, this, _1, MX, MZ, MY);
    break;
  case BULK:
    bY0 = bind(&Boundary2D::bY0Bulk, this, _1, MX, MZ, bulk);
    break;
  }

  switch (ym) {
  case MIRROR:
    bYm = bind(&Boundary2D::bYmMirror, this, _1, MX, MZ, MY);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bYm = bind(&Boundary2D::bYPeriodic, this, _1, MX, MZ, MY);
    break;
  case BULK:
    bYm = bind(&Boundary2D::bYmBulk, this, _1, MX, MZ, MY, bulk);
    break;
  }

  return 0;
}

int Boundary3D::set_z_boundaries(boundary z0, boundary zm, Real bulk) {
    using namespace std::placeholders;

  switch (z0) {
  case MIRROR:
    bZ0 = bind(&Boundary3D::bZ0Mirror, this, _1, MX, MY);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bZ0 = bind(&Boundary3D::bZPeriodic, this, _1, MX, MY, MZ);
    break;
  case BULK:
    bZ0 = bind(&Boundary3D::bZ0Bulk, this, _1, MX, MY, bulk);
    break;
  }

  switch (zm) {
  case MIRROR:
    bZm = bind(&Boundary3D::bZmMirror, this, _1, MX, MY, MZ);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bZm = bind(&Boundary3D::bZPeriodic, this, _1, MX, MY, MZ);
    break;
  case BULK:
    bZm = bind(&Boundary3D::bZmBulk, this, _1, MX, MY, MZ, bulk);
    break;
  }

  return 0;
}

void Boundary1D::bX0Mirror(vector<Real>& target, int fMY, int fMZ) {
  int y = 0;
  int z = 0;
  do {
    y=0;
    do {
      *valPtr(target, 0, y, z) = val(target, 1, y, z);     //start
      ++y;
    } while (y < fMY);
    ++z;
  } while (z < fMZ);
}

void Boundary1D::bXmMirror(vector<Real>& target, int fMY, int fMZ, int fMX) {
  int y = 0;
  int z = 0;
  do {
    y=0;
    do {
      *valPtr(target, fMX - 1, y, z) = val(target, fMX - 2, y, z);     //end
      ++y;
    } while (y < fMY);
  ++z;
  } while (z < fMZ);
}

void Boundary1D::bXPeriodic(vector<Real>& target, int fMY, int fMZ, int fMX) {
  int y = 0;
  int z = 0;
  do {
    y=0;
    do {
      *valPtr(target, 0, y, z) = val(target, fMX - 2, y, z); //start
      *valPtr(target, fMX - 1, y, z) = val(target, 1, y, z); //end
      ++y;
    } while (y < fMY);
    ++z;
  } while (z < fMZ);
}

void Boundary1D::bX0Bulk(vector<Real>& target, int fMY, int fMZ, Real bulk) {
  int y = 0;
  int z = 0;
  do {
    y=0;
    do {
      *valPtr(target, 0, y, z) = bulk; //start
      ++y;
    } while (y < fMY);
    ++z;
  } while (z < fMZ);
}

void Boundary1D::bXmBulk(vector<Real>& target, int fMY, int fMZ, int fMX, Real bulk) {
  int y = 0;
  int z = 0;
  do {
    y=0;
    do {
      *valPtr(target, fMX - 1, y, z) = bulk; //end
      ++y;
    } while (y < fMY);
    ++z;
  } while (z < fMZ);
}

void Boundary2D::bY0Mirror(vector<Real>& target, int fMX, int fMZ) {
  int x = 0;
  int z = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, 0, z) = val(target, x, 1, z);     //start
      ++x;
    } while (x < fMX);
    ++z;
  } while (z < fMZ);
}

void Boundary2D::bYmMirror(vector<Real>& target, int fMX, int fMZ, int fMY) {
  int x = 0;
  int z = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, fMY - 1, z) = val(target, x, fMY - 2, z);     //end
      ++x;
    } while (x < fMX);
    ++z;
  } while (z < fMZ);
}

void Boundary2D::bYPeriodic(vector<Real>& target, int fMX, int fMZ, int fMY) {
  int x = 0;
  int z = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, 0, z) = val(target, x, fMY - 2, z); //start
      *valPtr(target, x, fMY - 1, z) = val(target, x, 1, z); //end
      ++x;
    } while (x < fMX);
    ++z;
  } while (z < fMZ);
}

void Boundary2D::bY0Bulk(vector<Real>& target, int fMX, int fMZ, Real bulk) {
  int x = 0;
  int z = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, 0, z) = bulk; //start
      ++x;
    } while (x < fMX);
    ++z;
  } while (z < fMZ);
}

void Boundary2D::bYmBulk(vector<Real>& target, int fMX, int fMZ, int fMY, Real bulk) {
  int x = 0;
  int z = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, fMY - 1, z) = bulk; //end
      ++x;
    } while (x < fMX);
    ++z;
  } while (z < fMZ);
}

void Boundary3D::bZ0Mirror(vector<Real>& target, int fMX, int fMY) {
  int x = 0;
  int y = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, y, 0) = val(target, x, y, 1);     //start
      ++x;
    } while (x < fMX);
  ++y;
  } while (y < fMY);
}

void Boundary3D::bZmMirror(vector<Real>& target, int fMX, int fMY, int fMZ) {
  int x = 0;
  int y = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, y, fMZ - 1) = val(target, x, y, fMZ - 2);     //end
      ++x;
    } while (x < fMX);
  ++y;
  } while (y < fMY);
}

void Boundary3D::bZPeriodic(vector<Real>& target, int fMX, int fMY, int fMZ) {
  int x = 0;
  int y = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, y, 0) = val(target, x, y, fMZ - 2); //start
      *valPtr(target, x, y, fMZ - 1) = val(target, x, y, 1); //end
      ++x;
    } while (x < fMX);
    ++y;
  } while (y < fMY);
}

void Boundary3D::bZ0Bulk(vector<Real>& target, int fMX, int fMY, Real bulk) {
  int x = 0;
  int y = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, y, 0) = bulk; //start
      ++x;
    } while (x < fMX);
    ++y;
  } while (y < fMY);
}

void Boundary3D::bZmBulk(vector<Real>& target, int fMX, int fMY, int fMZ, Real bulk) {
  int x = 0;
  int y = 0;
  do {
    x=0;
    do {
      *valPtr(target, x, y, fMZ - 1) = bulk; //end
      ++x;
    } while (x < fMX);
    ++y;
  } while (y < fMY);
}

/******* GAUSSIAN_NOISE: GENERATES WHITE NOISE FOR FLUXES ********/

Gaussian_noise::Gaussian_noise(Boundary1D* boundary, Real D, int size, Real mean, Real stddev) : prng { std::random_device{} () }, dist(mean, stddev), noise(size), boundary{boundary} {}

Gaussian_noise::Gaussian_noise(Boundary1D* boundary, Real D, int size, Real mean, Real stddev, size_t seed) : prng(seed), dist(mean, stddev), noise(size), boundary{boundary} {}

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
