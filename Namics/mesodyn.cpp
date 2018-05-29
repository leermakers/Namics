#include "mesodyn.h"

Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_, string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sys{Sys_},
      New{New_},
      D{0.5},
      mean{1},
      stdev{1},
      seed{1},
      timesteps{0},
      timebetweensaves{0},
      JZ{Lat[0]->JZ}, // Usage for e.g. layer z: foo[z]+foo[z+1] becomes foo[z] + foo[z+JX]
      JY{Lat[0]->JY},
      JX{Lat[0]->JX},
      M{Lat[0]->M},                                   // Neighboring component
      componentNo{(int)Sys[0]->SysMolMonList.size()}, //find how many compontents there are (e.g. head, tail, solvent)
      size{(unsigned int)componentNo * M},                          //find out how large the density vector is (needed for sizing the flux vector)
                                                      //which will be 1 flux per lattice site per component per dimension
      dimensions{findDimensions()},                    // used to decide which fluxes to calculate
      J(dimensions * size),
      L(combinations(componentNo, 2) * M),
      rho(size),
      ptrComponentStart(componentNo),
      U(componentNo*(componentNo-1)*M)
{
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stdev");
  if (debug)
    cout << "Mesodyn initialized" << endl;

}

Mesodyn::~Mesodyn() {
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
      success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 10000, "The number of timesteps should be between 1 and 10000");
    }
    if (debug)
      cout << "Timesteps is " << timesteps << endl;

    if (GetValue("timebetweensaves").size() > 0) {
      success = In[0]->Get_int(GetValue("timebetweensaves"), timebetweensaves);
    }
    if (debug)
      cout << "Time bewteen saves is " << timebetweensaves << endl;

    if (GetValue("diffusionconstant").size() > 0) {
      success = In[0]->Get_Real(GetValue("diffusionconstant"), D);
    }
    if (debug)
      cout << "Diffusion const is " << D << endl;

    if (GetValue("seed").size() > 0) {
      success = In[0]->Get_Real(GetValue("seed"), seed);
    }
    if (debug)
      cout << "Seed is " << seed << endl;

    if (GetValue("mean").size() > 0) {
      success = In[0]->Get_Real(GetValue("mean"), mean);
    }
    if (debug)
      cout << "Mean is " << mean << endl;

    if (GetValue("stdev").size() > 0) {
      success = In[0]->Get_Real(GetValue("stdev"), stdev);
    }
    if (debug)
      cout << "Stdev is " << stdev << endl;
  }

  return success;
}

//Used by langevinFlux() to calculate the correct number of fluxes
int Mesodyn::findDimensions() {
  int d = 0;
  if (Lat[0]->MX > 0)
    ++d;
  if (Lat[0]->MY > 0)
    ++d;
  if (Lat[0]->MZ > 0)
    ++d;
  return d;
}


/******** Flow control ********/

bool Mesodyn::mesodyn() {
  if (debug)
    cout << "mesodyn in Mesodyn" << endl;
  initRho(); // set initial conditions
  prepareOutputFile();
  writeRho(0); // write initial conditions
  if (debug)
    cout << "Mesodyn is all set, starting calculations.." << endl;
  for (int t = 1; t < 100; t++) { // get segment potentials by iteration so that they match the given rho.
        New[0]->SolveMesodyn(rho);
        onsagerCoefficient();
        potentialDifference();
        langevinFlux();
        updateDensity();
        cout << "And LOOP!";
        writeRho(t);
    }
  return true;
}

//defaults to homogeneous system for now
void Mesodyn::initRho() {
  if (debug) cout << "initRho in Mesodyn." << endl;

  //Not resizing will break .size() operations like in mesodyn().
  rho.resize(size);
  vector<Real>::iterator rhoIterator;

  Real homogeneous = (Real)1 / componentNo;

  for (rhoIterator = rho.begin() ; rhoIterator != rho.end() ; ++rhoIterator){
      *rhoIterator = homogeneous;
    }

  //Next, register where in rho all components are located:
  //start indices for each segment type are pointed to by ptrComponentStart.
  for (int i {0}; i < componentNo; ++i) {
    ptrComponentStart.at(i) = &rho[M * i];
  }
}

/******** Calculations ********/
void Mesodyn::onsagerCoefficient() {
  //Tested: Works as intended.

  if (debug)
    cout << "onsagerCoefficient in Mesodyn." << endl;

  //TODO: maybe do this inline / per J calculation to preserve memory
  vector<Real>::iterator lIterator;
  lIterator = L.begin();

  //all combinations of ptrComponentStart multiplications
  for (int i = 0; i < componentNo - 1; ++i) {
    for (int j = i + 1; j < componentNo; ++j) {
      //at every coordinate (pointer arithmatic)
      for (int xyz = 0; xyz < M; ++xyz) {
        *lIterator = (*(ptrComponentStart[i] + xyz) * *(ptrComponentStart[j] + xyz));

        //No bounds checking needed because it runs out of scope right after.
        lIterator++;
      }
    }
  }
}

void Mesodyn::potentialDifference() {
  //Structured in much the same way as the onsager coefficients.
  if (debug)
    cout << "potentialDifference in Mesodyn." << endl;

  //TODO: maybe do this inline / per J calculation to preserve memory
  vector<Real>::iterator uIterator;
  uIterator = U.begin();

  int jump = componentNo-1;

  //for all components
  for (int i = 0; i < componentNo; ++i) {
    //for all components-1 = number of combinations possible
    for (int j = 0; j < componentNo-1; ++j){
      //for xyz
      for (int z = 0 ; z < M ; ++z) {

    //TODO: can do this with another for statement (for ( i <= jump ; something = j+1 ))

        if (i <= jump) {
          //subtract the potential of the other components j or j+1 from the potential the current component i
          *uIterator = New[0]->xx[z+i*M] - New[0]->xx[z+j*M];
          ++uIterator;
        } else {
          *uIterator = New[0]->xx[z+i*M] - New[0]->xx[z+(j+1)*M];
          ++uIterator;
          //No bounds checking needed because it runs out of scope right after.
        }
      }
    }
  }
}

void Mesodyn::langevinFlux() {
  if (debug)
    cout << "langevinFlux in Mesodyn." << endl;

  //TODO: safer to change size calculation to xx.size()?
  //TODO: boundary condition in lattice?

  //This jIterator is going to keep track of where in the flux vector we are and keeps incrementing over all the coming loops.
  //so that fluxes are sequentially placed like Jx1, Jy1, Jz1, Jx2, Jy2, Jz2, where the number is the (indexed) lattice site.
  //(this structure is important to know when updating the densities accordingly).
  vector<Real>::iterator jIterator;

  // Zero the jIterator
  for (jIterator = J.begin(); jIterator != J.end(); ++jIterator) {
    *jIterator = 0;
  };

  jIterator = J.begin();

  /*  This next part finds a vector that selects which of the combinations of onsager coefficients
      should be used for that particular component. For example for components A, B, C, D we get:
      AB-AC-AD-BC-BD-CD. For component A (i=0) we need indices 0, 1, 2. For component B (i=1) 0,3,4. for C 1,3,5.
      It finds the correct indices, puts them in cCombinations, calculates the fluxes that correspond to that component and loops.
  */

  int j = 1;
  vector<int> cCombinations(componentNo - 1);
  vector<int>::iterator nIterator;

  // So this for loop runs over all the fluxes, each loop represents 1 component.
  for (int i = 0; i < componentNo; ++i) {

    //First, we calculate which at which indices the correct onsager coefficients can be found.
    nIterator = cCombinations.begin() + i;
    while (nIterator != cCombinations.end()) {
      *nIterator = j;
      ++j;
      ++nIterator;
    }
    if (i > 1) {
      for (int k = 0; k < i - 1; ++k) {
        cCombinations[k] = ++cCombinations[k];
      }
    }

    // Then, we cancluate the fluxes for each dimension.
    // This assumes that if we have one dimension, it's the x dimension, whereas if we have 2, it's x & y.

    //for all lattice sites
    for (int z = 0; z < Lat[0]->M; ++z) {
      //Generate noise for this flux.
      //TODO: one noise, two noise, three noise? Different per dimension?
      gaussianNoise(mean, stdev, 1);
      //for all combinations with other components
      for (int l = 0; l < componentNo - 1; ++l) {
        *jIterator += -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JX]) * (U[(i*(componentNo-1)+l)*M + JX] - U[(i*(componentNo-1)+l)*M])) - ((L[z + cCombinations[l] * M - JX] + L[z + cCombinations[l] * M]) * (U[z+(i*(componentNo-1)+l)*M] - U[(i*(componentNo-1)+l)*M - JX])) + noise[0];
      }
      ++jIterator;
    }
    if (dimensions > 1) {
      for (int z = 0; z < Lat[0]->M; ++z) {
        gaussianNoise(mean, stdev, 1);
        for (int l = 0; l < componentNo - 1; ++l) {
          *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JY]) * (U[(i*(componentNo-1)+l)*M + JY] - U[(i*(componentNo-1)+l)*M])) - ((L[z + cCombinations[l] * M - JY] + L[z + cCombinations[l] * M]) * (U[(i*(componentNo-1)+l)*M] - U[(i*(componentNo-1)+l)*M - JY])) + noise[0];
        }
        ++jIterator;
      }
    }
    if (dimensions > 2) {
      for (int z = 0; z < Lat[0]->M; ++z) {
        gaussianNoise(mean, stdev, 1);
        for (int l = 0; l < componentNo - 1; ++l) {
          *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JZ]) * (U[(i*(componentNo-1)+l)*M + JZ] - U[(i*(componentNo-1)+l)*M])) - ((L[z + cCombinations[l] * M - JZ] + L[z + cCombinations[l] * M]) * (U[(i*(componentNo-1)+l)*M] - U[(i*(componentNo-1)+l)*M - JZ])) + noise[0];
        }
        ++jIterator;
      }
    }
    // And loop!
  }
}

void Mesodyn::updateDensity() {
  if (debug)
    cout << "updateDensity in Mesodyn." << endl;

  // J looks like: Jx1 - Jy1 - Jz1 - Jx2 - Jy2 - Jz2 etc.
  // Each J has size Lat[0]->M

  for (int j = 0; j < componentNo; ++j) {
    for (int i = 0; i < M; ++i) {
      // Density = flux [lattice site * coordinate start index + component]
      rho[i + j * M] += J[i + j * M];         // x
      if (dimensions == 2) rho[i + j * M] += J[i * M + j * M];     // y
      if (dimensions == 3) rho[i + j * M] += J[i * 2 * M + j * M]; // z
    }
  }
}

/* Generates a vector of length count, contianing gaussian noise of given mean, standard deviation.
	 Noise is stored in vector<Real> Mesodyn::noise
	 Possible errors: What if count > sizeof(unsinged long)?
	 Called by langevinFlux()
*/
void Mesodyn::gaussianNoise(Real mean, Real stdev, unsigned int count) {

  random_device generator;

  seed_seq seed("something", "something else");

  //Mersenne Twister 19937 bit state size PRNG
  mt19937 prng(seed);

  normal_distribution<> dist(mean, stdev);

  this->noise.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
    this->noise[i] = dist(prng);
  }

/* Debugging code (output value in all elements):
  for (auto const &element: mesodyn.thisNoise)
				std::cout << element << ' ';
*/
}

/******* Output generation *******/

void Mesodyn::prepareOutputFile() {
  /* Open filestream and set filename to "mesodyn-datetime.csv" */
  ostringstream filename;
  filename << "mesodyn-";

  time_t rawtime;
  time (&rawtime);
  filename << rawtime << ".csv";

  mesFile.open(filename.str());

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

  headers << "x:y:z,";

  for (int i = 1; i <= componentNo; ++i) {
    headers << "rho" << i << ",";
  }

  headers << "\n";

  mesFile << headers.str();
}

void Mesodyn::writeRho(int t) {
  ostringstream timeOutput;
  timeOutput << t << ",\n";
  mesFile << timeOutput.str();

  ostringstream rhoOutput;

  for (int z = 0 ; z < Lat[0]->MZ ; ++z) {
    for (int y = 0 ; y < Lat[0]->MY ; ++y) {
      for (int x = 0 ; x < Lat[0]->MX ; ++x) {
          rhoOutput << x << ":" << y << ":" << z << ",";
          for (int c = 0; c < componentNo ; ++c) {
            rhoOutput << val(rho, c, x, y, z) << ",";
          }
          rhoOutput << "\n";
          mesFile << rhoOutput.str();

          //Set rhoOutput to empty
          rhoOutput.str("");
          //Clear any remaining error flags
          rhoOutput.clear();
      }
    }
  }

}


/******* Tools ********/

inline Real Mesodyn::val(vector<Real>& v, int c, int x, int y, int z) {
  return v[c * Lat[0]->M + x * Lat[0]->JX + y * Lat[0]->JY + z];
}

int Mesodyn::factorial(int n) {
  if (debug)
    cout << "factorial in Mesodyn." << endl;
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
  if (debug)
    cout << "combinations in Mesodyn." << endl;
  return factorial(n) / (factorial(n - k) * factorial(k));
}

void Mesodyn::PutParameter(string new_param) {
  KEYS.push_back(new_param);
}
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
int Mesodyn::GetValue(string prop, int& int_result, Real& Real_result, string& string_result) {
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
