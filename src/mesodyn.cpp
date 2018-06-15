#include "mesodyn.h"
#include "newton.h"

Mesodyn::Mesodyn(vector<Input *> In_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<Molecule *> Mol_, vector<System *> Sys_, vector<Solve_scf*> New_, string name_)
  :
    Access(Lat_[0]),
    name{name_}, In{In_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sys{Sys_}, New{New_},
    D{0.001}, mean{1}, stdev{1}, seed{1}, timesteps{100}, timebetweensaves{1},
    componentNo{(int)Sys[0]->SysMolMonList.size()},
    alpha(componentNo*M)

{
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stdev");

  initial_conditions(); // Initialize densities by running the classical method once.

  if (debug) cout << "Mesodyn initialized." << endl;
}

Mesodyn::~Mesodyn() {
  for (unsigned int i = 0 ; i < flux.size(); ++i) delete flux[i];
  flux.clear();
  for (unsigned int i = 0 ; i < component.size(); ++i) delete component[i];
  flux.clear();
}

bool Mesodyn::CheckInput(int start) {
  if (debug) cout << "CheckInput in Mesodyn" << endl;
  bool success = true;

  /* When adding new items here, please add them to the function prepareOutputFile()*/
  success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);

  if (success) {
    vector<string> options;
    if (GetValue("timesteps").size() > 0) { success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 100000000, "The number of timesteps should be between 1 and 10000");}
    if (debug) cout << "Timesteps is " << timesteps << endl;

    if (GetValue("timebetweensaves").size() > 0) {timebetweensaves = In[0]->Get_int(GetValue("timebetweensaves"), timebetweensaves);}
    if (debug) cout << "Time bewteen saves is " << timebetweensaves << endl;

    if (GetValue("diffusionconstant").size() > 0) {D = In[0]->Get_Real(GetValue("diffusionconstant"), D);}
    if (debug)  cout << "Diffusion const is " << D << endl;

    if (GetValue("seed").size() > 0) {seed = In[0]->Get_Real(GetValue("seed"), seed);}
    if (debug) cout << "Seed is " << seed << endl;

    if (GetValue("mean").size() > 0) {mean = In[0]->Get_Real(GetValue("mean"), mean);}
    if (debug) cout << "Mean is " << mean << endl;

    if (GetValue("stdev").size() > 0) { stdev = In[0]->Get_Real(GetValue("stdev"), stdev);}
    if (debug) cout << "Stdev is " << stdev << endl;
  }

  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  if (debug)
    cout << "mesodyn in Mesodyn" << endl;

  prepareOutputFile();
  writeRho(0); // write initial conditions

  if (debug) cout << "Mesodyn is all set, starting calculations.." << endl;

  for (int t = 1; t < timesteps; t++) {
    cout << "MESODYN: t = " << t << endl;

    // TODO: once phi becomes a vector, this will be much sexier.
    vector<Real> rho;
    for (Component1D* all_components : component) {
      copy( all_components->rho.begin(),  all_components->rho.end(), back_inserter(rho));
    }

    New[0]->SolveMesodyn(rho, alpha);

    for (int i = 0 ; i < componentNo ; ++i) component[i]->load_alpha(&alpha[0+i*M], M);

    for (Component1D* all_components : component) all_components->update_boundaries();

    for (Flux1D* all_fluxes : flux) all_fluxes->langevin_flux();

    if (t % timebetweensaves == 0) writeRho(t);

    int c = 0;
    for (int i = 0 ; i < componentNo ; ++i) {
      for (int j = 0 ; j < (componentNo - 1) - i ; ++j) {
        component[i]->update_density( flux[c]->J );
        ++c;
      }
    }

    c = 0;
    for (int j = 0 ; j < (componentNo - 1) ; ++j) {
      for (int i = 1+j ; i < componentNo; ++i) {
        component[i]->update_density( flux[c]->J , -1.0 );
        ++c;
      }
    }

  } // time loop
  writeRho(timesteps);
  return true;
}

int Mesodyn::initial_conditions() {

    //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
    for (Segment* seg : Seg) {
      if (seg->freedom == "pinned")
        seg->freedom = "free";
    }

    vector< vector<Real> > rho(componentNo, vector<Real>(M));

    //TODO: generalize (M-1-volume?) for 2D/3D
    Sys[0]->PrepareForCalculations();
    int solvent = Sys[0]->solvent;
    int volume = Sys[0]->volume-(pow(M,dimensions) - pow ((M-2),dimensions));

    Real sum_theta {0};
    Real theta {0};

    for (int i = 0; i < componentNo ; ++i) {
      if (i != solvent) {
        theta = Mol[i]->theta;
        sum_theta += theta;
        for (int z = M-1-volume; z < Lat[0]->M-1; z++) {
          rho[i][z] = theta/volume;
        }
      }
    }

    for (int z = M-1-volume ; z < Lat[0]->M-1 ; ++z) {
      rho[solvent][z] = (volume-sum_theta)/volume;
    }

    //TODO: wall boundary
    for (int i = 0 ; i < componentNo ; ++i) {
      rho[i][1] = 0;
    }

    vector<Component1D::boundary> boundaries;
    for (string& boundary : Lat[0]->BC) {
      if (boundary == "mirror") {
        boundaries.push_back(Component1D::MIRROR);
      }
      if (boundary == "periodic") {
        boundaries.push_back(Component1D::PERIODIC);
      }
      if (boundary == "bulk") {
        boundaries.push_back(Component1D::BULK);
      }
      if (boundary == "surface") {
        boundaries.push_back(Component1D::SURFACE);
      }
    }

    // TODO: specify boundary conditions in constructor
    switch (dimensions) {
    case 1:
      for (int i = 0 ; i < componentNo ; ++i)
          component.push_back( new Component1D( Lat[0], rho[i], boundaries[0], boundaries[1] ));

      for (int i = 0; i < componentNo - 1; ++i) {
        for (int j = i + 1; j < componentNo; ++j) {
            flux.push_back( new Flux1D( Lat[0], D, component[i], component[j] ) );
        }
      }


      break;
    case 2:
      for (int i = 0 ; i < componentNo ; ++i)
        component.push_back( new Component2D( Lat[0], rho[i], boundaries[0], boundaries[1], boundaries[2], boundaries[3] ));

      for (int i = 0; i < componentNo - 1; ++i) {
        for (int j = i + 1; j < componentNo; ++j) {
          flux.push_back( new Flux2D( Lat[0], D, component[i], component[j] ) );
        }
      }

      break;
    case 3:
      for (int i = 0 ; i < componentNo ; ++i)
        component.push_back( new Component3D( Lat[0], rho[i], boundaries[0], boundaries[1], boundaries[2], boundaries[3], boundaries[4], boundaries[5] ));

      for (int i = 0; i < componentNo - 1; ++i) {
        for (int j = i + 1; j < componentNo; ++j) {
          flux.push_back( new Flux3D( Lat[0], D, component[i], component[j] ) );
        }
      }

      break;
    }

      //  New[0]->DeAllocateMemory();
      //  New[0]->AllocateMemory(componentNo);

  return 0;
}

/******* Rho initialization *******/

/******* Output generation *******/

void Mesodyn::prepareOutputFile() {
  /* Open filestream and set filename to "mesodyn-datetime.csv" */
  ostringstream filename;

  // TODO: Get all of this stuff below from output!
  string output_folder = "output/";
  string bin_folder = "bin";

  // Find path to Namics executable
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  string executable_path = string( result, (count > 0) ? count : 0 );

  // Find the last string before the executable
  size_t found = executable_path.find_last_of("/\\");

  // Set the output folder to be one level up from the binary folder, plus the specified output folder
  output_folder = executable_path.substr(0,found - bin_folder.size() ) + output_folder;

  filename << output_folder << "mesodyn-";

  time_t rawtime;
  time(&rawtime);
  filename << rawtime << ".csv";

  mesFile.open(filename.str());

  mesFile.precision(16);

  /* Output settings */

  ostringstream settings;

  settings << "Settings for mesodyn, read from file:,\n";
  settings << "Timesteps," << timesteps << ",\n";
  settings << "Time between saves," << timebetweensaves << ",\n";
  settings << "Diffusion constant," << D << ",\n";
  settings << "Seed for noise," << seed << ",\n";
  settings << "Mean for noise," << mean << ",\n";
  settings << "Stdev for noise," << stdev << ",\n";

  mesFile << settings.str();

  /* Generate headers for rho entries */

  ostringstream headers;

  int permutations = combinations(componentNo, 2);

  headers << "x:y:z,";

  for (int i = 1; i <= componentNo; ++i) {
    headers << "rho" << i << ",";
  }

  for (int i = 1; i <= permutations; ++i) {
    headers << "L" << i << ",";
  }

  for (int i = 1; i <= componentNo; ++i) {
    headers << "alpha" << i << ",";
  }

  for (int i = 1; i <= permutations; ++i) {
    headers << "mu" << i << ",";
  }

  for (int i = 1; i <= permutations; ++i) {
    headers << "J" << i << ",";
  }

  headers << "\n";

  mesFile << headers.str();
}

void Mesodyn::writeRho(int t) {
  ostringstream headers;

  headers << t << "," << "\n";
  mesFile << headers.str();

  ostringstream rhoOutput;

  rhoOutput.precision(16);
  int x{0}, y{0}, z{0};

  do {
    do {
      do {
        rhoOutput << x << ","; //":" << y << ":" << z << ",";
        for (Component1D *all_components : component) {
          rhoOutput << all_components->rho_at(x, y, z) << ",";
        }

        for (Flux1D *all_fluxes : flux) {
          rhoOutput << all_fluxes->L_at(x, y, z) << ",";
        }

        for (Component1D *all_components : component) {
          rhoOutput << all_components->alpha_at(x, y, z) << ",";
        }

        for (Flux1D *all_fluxes : flux) {
          rhoOutput << all_fluxes->mu_at(x, y, z) << ",";
        }

        for (Flux1D *all_fluxes : flux) {
          rhoOutput << all_fluxes->J_at(x, y, z) << ",";
        }

        rhoOutput << "\n";

        mesFile << rhoOutput.str();

        // Set rhoOutput to empty
        rhoOutput.str("");
        // Clear any remaining error flags
        rhoOutput.clear();

        x++;
      } while (x < MX);
      y++;
    } while (y < MY);
    z++;
  } while (z < MZ);

  mesFile.flush();
}

/******* FLUX *********/

Flux1D::Flux1D(Lattice* Lat, Real D, Component1D* A, Component1D* B)
    : Access(Lat), J(M), A{A}, B{B}, L(M), mu(M), D{D}, JX{Lat->JX}   {
}

Flux2D::Flux2D(Lattice* Lat, Real D, Component1D* A, Component1D* B)
    : Flux1D(Lat, D, A, B), JY{Lat->JY} {
}

Flux3D::Flux3D(Lattice* Lat, Real D, Component1D* A, Component1D* B)
    : Flux2D(Lat, D, A, B), JZ{Lat->JZ} {
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

  if (A->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  onsager_coefficient(A->rho, B->rho);
  potential_difference(A->alpha, B->alpha);
  langevin_flux(JX);

  return 0;
}

int Flux2D::langevin_flux() {
  Flux1D::langevin_flux();

  Flux1D::langevin_flux(JY);

  return 0;
}

int Flux3D::langevin_flux() {
  Flux2D::langevin_flux();

  Flux1D::langevin_flux(JZ);

  return 0;
}

Real Flux1D::J_at(int x, int y, int z) {
  return val(J, x, y ,z);
}

Real Flux1D::L_at(int x, int y, int z) {
  return val(L, x, y ,z);
}

Real Flux1D::mu_at(int x, int y, int z) {
  return val(mu, x, y ,z);
}

int Flux1D::onsager_coefficient(vector<Real>& A, vector<Real>& B) {
  //TODO: maybe do this inline / per J calculation to preserve memory

  if (A.size() != B.size() ) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  transform(A.begin(), A.end(), B.begin(), L.begin(), [](Real A, Real B) { return A * B; });
  return 0;
}

int Flux1D::potential_difference(vector<Real>& A, vector<Real>& B) {
  //TODO: maybe do this inline / per J calculation to preserve memory

  if (A.size() != B.size() ) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  transform(A.begin(), A.end(), B.begin(), mu.begin(), [](Real A, Real B) { return A - B; });
  return 0;
}

int Flux1D::langevin_flux(int jump) {
  for (int z = 1; z < M-1; ++z){
  J[z] += -D * (((L[z] + L[z + jump]) * (mu[z + jump] - mu[z])) + ((L[z - jump] + L[z]) * (mu[z - jump] - mu[z])));
  }
  return 0;
}

/****************** ACCESS ********************/

Access::Access(Lattice* Lat)
  : dimensions{Lat->gradients}, JX {Lat->JX}, JY {Lat->JY}, JZ {Lat->JZ}, M {Lat->M}, MX {2+Lat->MX}, MY { setMY(Lat) }, MZ{ setMZ(Lat) }
{
  // If this is not true, NOTHING will work. So this check is aboslutely necessary.
  assert(
    (MX > 0 && MY == 0 && MZ == 0) ||
    (MX > 0 && MY > 0  && MZ == 0) ||
    (MX > 0 && MY > 0  && MZ > 0)
  );
}

Access::~Access() {

}

inline Real Access::val(vector<Real>& v, int x, int y, int z) {
  return v[x * JX + y * JY + z];
}

inline Real* Access::valPtr(vector<Real>& v, int x, int y, int z) {
  return &v[x * JX + y * JY + z];
}

inline int Access::xyz(int x, int y, int z) {
  return (x * JX + y * JY + z);
}

int Access::setMY(Lattice* Lat) {
  if (dimensions < 2) return 0;
  else return Lat->MY+2;
}

int Access::setMZ(Lattice* Lat) {
  if (dimensions < 3) return 0;
  else return Lat->MZ+2;
}

/****************** COMPONENT ********************/


          /******* Constructors *******/

Component1D::Component1D(Lattice* Lat, vector<Real>& rho, boundary x0, boundary xm)
    : Access(Lat), rho{rho}, alpha(M)
{
  //This check is implemented multiple times throughout mesodyn because rho and alpha are public.
  if ( rho.size() != alpha.size() ) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  set_x_boundaries(x0, xm);
  update_boundaries();
}

Component2D::Component2D(Lattice* Lat, vector<Real>& rho, boundary x0, boundary xm, boundary y0, boundary ym)
    : Component1D(Lat, rho, x0, xm)
{
  set_y_boundaries(y0, ym);
  update_boundaries();
}

Component3D::Component3D(Lattice* Lat, vector<Real>& rho, boundary x0, boundary xm, boundary y0, boundary ym, boundary z0, boundary zm)
    : Component2D(Lat, rho, x0, xm, y0, ym)
{
  set_z_boundaries(z0, zm);
  update_boundaries();
}

Component1D::~Component1D() {
}

Component2D::~Component2D() {
}

Component3D::~Component3D() {
}


/******* Interface *******/

Real Component1D::rho_at(int x, int y, int z) {
  return val(rho, x, y, z);
}

Real Component1D::alpha_at(int x, int y, int z) {
  return val(alpha, x, y ,z);
}

int Component1D::update_density(vector<Real>& J, int sign) {
  //Explicit update

  if (J.size() != rho.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }

  int i = 0;

  for (Real& Flux : J) {
    rho[i] += sign*Flux;
    ++i;
  }

  return 0;
}

int Component1D::update_density(vector<Real>& J1, vector<Real>& J2) {
  //Implicit update

  throw ERROR_NOT_IMPLEMENTED;
  return 1;
}

int Component1D::load_alpha(Real* alpha, int m) {
  //  TODO: C++11ify this hideous function please
  if (m != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }
  Component1D::alpha.assign(alpha,alpha+m);
  return 0;
}

int Component1D::load_rho(Real* rho, int m) {
  //  TODO: C++11ify this hideous function please
  if (m != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
    return 1;
  }
  Component1D::rho.assign(rho,rho+m);
  return 0;
}

int Component1D::update_boundaries() {
  bX0();
  bXm();
  return 0;
}

int Component2D::update_boundaries() {
  Component1D::update_boundaries();
  bY0();
  bYm();
  return 0;
}

int Component3D::update_boundaries() {
  Component2D::update_boundaries();
  bZ0();
  bZm();
  return 0;
}

/******* Boundary conditions *******/

int Component1D::set_x_boundaries(boundary x0, boundary xm) {
  switch (x0) {
  case MIRROR:
    bX0 = bind(&Component1D::bX0Mirror, this, MY, MZ);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bX0 = bind(&Component1D::bXPeriodic, this, MY, MZ, MX);
    break;
  case BULK:
    bX0 = bind(&Component1D::bX0Bulk, this, MY, MZ, rho[1]);
    break;
  case SURFACE:
    //nothing yet
    break;
  }

  switch (xm) {
  case MIRROR:
    bXm = bind(&Component1D::bXmMirror, this, MY, MZ, MX);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    break;
  case BULK:
    bXm = bind(&Component1D::bXmBulk, this, MY, MZ, MX, rho[MX-2]);
    break;
  case SURFACE:
   // nothing yet
    break;
  }

  return 0;
}

int Component2D::set_y_boundaries(boundary y0, boundary ym) {
  switch (y0) {
  case MIRROR:
    bY0 = bind(&Component2D::bY0Mirror, this, MX, MZ);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bY0 = bind(&Component2D::bYPeriodic, this, MX, MZ, MY);
    break;
  case BULK:
    //nothing yet
  case SURFACE:
   // nothing yet
    break;
  }

  switch (ym) {
  case MIRROR:
    bYm = bind(&Component2D::bYmMirror, this, MX, MZ, MY);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    break;
  case BULK:
    //nothing yet
    break;
  case SURFACE:
    // nothing yet
    break;
  }

  return 0;
}

int Component3D::set_z_boundaries(boundary z0, boundary zm) {
  switch (z0) {
  case MIRROR:
    bZ0 = bind(&Component3D::bZ0Mirror, this, MX, MY);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    bZ0 = bind(&Component3D::bZPeriodic, this, MX, MY, MZ);
    break;
  case BULK:
    //nothing yet
  case SURFACE:
   // nothing yet
    break;
  }

  switch (zm) {
  case MIRROR:
    bZm = bind(&Component3D::bZmMirror, this, MX, MY, MZ);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
      return 1;
    }
    break;
  case BULK:
    //nothing yet
    break;
  case SURFACE:
    // nothing yet
    break;
  }

  return 0;
}

void Component1D::bX0Mirror(int fMY, int fMZ) {
  int y = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, 0, y, z) = val(rho, 1, y, z);     //start
      *valPtr(alpha, 0, y, z) = val(alpha, 1, y, z); //start
      *valPtr(alpha, 1, y, z) = val(alpha, 2, y, z); //start // DELETE WALL
      ++z;
    } while (z < fMZ);

    ++y;
  } while (y < fMY);
}

void Component1D::bXmMirror(int fMY, int fMZ, int fMX) {
  int y = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, fMX - 1, y, z) = val(rho, fMX - 2, y, z);     //end
      *valPtr(alpha, fMX - 1, y, z) = val(alpha, fMX - 2, y, z); //end
      ++z;
    } while (z < fMZ);

    ++y;
  } while (y < fMY);
}

void Component1D::bXPeriodic(int fMY, int fMZ, int fMX) {
  int y = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, 0, y, z) = val(rho, fMX - 2, y, z); //start
      *valPtr(rho, fMX - 1, y, z) = val(rho, 1, y, z); //end

      *valPtr(alpha, 0, y, z) = val(alpha, fMX - 2, y, z); //start
      *valPtr(alpha, fMX - 1, y, z) = val(alpha, 1, y, z); //end
    } while (z < fMZ);
    ++y;
  } while (y < fMY);
}

void Component1D::bX0Bulk(int fMY, int fMZ, Real bulk) {
  int y = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, 0, y, z) = bulk;     //start
      ++z;
    } while (z < fMZ);

    ++y;
  } while (y < fMY);
}

void Component1D::bXmBulk(int fMY, int fMZ, int fMX, Real bulk) {
  int y = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, fMX - 1, y, z) = bulk;     //end
      ++z;
    } while (z < fMZ);
    ++y;
  } while (y < fMY);
}

void Component2D::bY0Mirror(int fMX, int fMZ) {
  int x = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, x, 0, z) = val(rho, x, 1, z);     //start
      *valPtr(alpha, x, 0, z) = val(alpha, x, 1, z); //start
    } while (z < fMZ);
    ++x;
  } while (x < fMX);
}

void Component2D::bYmMirror(int fMX, int fMZ, int fMY) {
  int x = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, x, fMY - 1, z) = val(rho, x, fMY - 2, z);     //end
      *valPtr(alpha, x, fMY - 1, z) = val(alpha, x, fMY - 2, z); //end
    } while (z < fMZ);
    ++x;
  } while (x < fMX);
}

void Component2D::bYPeriodic(int fMX, int fMZ, int fMY) {
  int x = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, x, 0, z) = val(rho, x, fMY - 2, z); //start
      *valPtr(rho, x, fMY - 1, z) = val(rho, x, 1, z); //end

      *valPtr(alpha, x, 0, z) = val(alpha, x, fMY - 2, z); //start
      *valPtr(alpha, x, fMY - 1, z) = val(alpha, x, 1, z); //end
    } while (z < fMZ);
    ++x;
  } while (x < fMX);
}

void Component2D::bY0Bulk(int fMX, int fMZ, Real bulk) {
  int x = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, x, 0, z) = bulk;     //start
    } while (z < fMZ);
    ++x;
  } while (x < fMX);
}

void Component2D::bYmBulk(int fMX, int fMZ, int fMY, Real bulk) {
  int x = 0;
  int z = 0;
  do {
    do {
      *valPtr(rho, x, fMY - 1, z) = bulk;     //end
    } while (z < fMZ);
    ++x;
  } while (x < fMX);
}

void Component3D::bZ0Mirror(int fMX, int fMY) {
  int x = 0;
  int y = 0;
  do {
    do {
      *valPtr(rho, x, y, 0) = val(rho, x, y, 1);     //start
      *valPtr(alpha, x, y, 0) = val(alpha, x, y, 1); //start
    } while (x < fMX);
    ++x;
  } while (y < fMY);
}

void Component3D::bZmMirror(int fMX, int fMY, int fMZ) {
  int x = 0;
  int y = 0;
  do {
    do {
      *valPtr(rho, x, y, fMZ - 1) = val(rho, x, y, fMZ - 2);     //end
      *valPtr(alpha, x, y, fMZ - 1) = val(alpha, x, y, fMZ - 2); //end
    } while (x < fMX);
    ++x;
  } while (y < fMY);
}

void Component3D::bZPeriodic(int fMX, int fMY, int fMZ) {
  int x = 0;
  int y = 0;
  do {
    do {
      *valPtr(rho, x, y, 0) = val(rho, x, y, fMZ - 2); //start
      *valPtr(rho, x, y, fMZ - 1) = val(rho, x, y, 1); //end

      *valPtr(alpha, x, y, 0) = val(alpha, x, y, fMZ - 2); //start
      *valPtr(alpha, x, y, fMZ - 1) = val(alpha, x, y, 1); //end
    } while (x < fMX);
    ++x;
  } while (y < fMY);
}

void Component3D::bZ0Bulk(int fMX, int fMY, Real bulk) {
  int x = 0;
  int y = 0;
  do {
    do {
      *valPtr(rho, x, y, 0) = bulk;     //start
    } while (x < fMX);
    ++x;
  } while (y < fMY);
}

void Component3D::bZmBulk(int fMX, int fMY, int fMZ, Real bulk) {
  int x = 0;
  int y = 0;
  do {
    do {
      *valPtr(rho, x, y, fMZ - 1) = bulk;     //end
    } while (x < fMX);
    ++x;
  } while (y < fMY);
}


/******* Tools ********/

int Mesodyn::factorial(int n) {
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
  return factorial(n) / (factorial(n - k) * factorial(k));
}

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
int Mesodyn::GetValue(string prop, int &int_result, Real &Real_result,
                      string &string_result) {
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
