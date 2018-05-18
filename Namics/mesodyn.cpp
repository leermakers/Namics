#include "mesodyn.h"

//Constructor
//TODO: Read D, noise seed from file.
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
      M{Lat[0]->M},                    // Neighboring component
      componentNo{(int)In[0]->MolList.size()}, //find how many compontents there are (e.g. head, tail, solvent)
      size{componentNo * Lat[0]->M},           //find out how large the density vector is (needed for sizing the flux vector)
                                               //which will be 1 flux per lattice site per component per dimension
      initRho{1 / (Real)componentNo},          // default is a homogeneous system.
      dimensions{findDimensions()}             // used to decide which fluxes to calculate
{
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stdev");
  cout << "Mesodyn initialized." << endl;
}

Mesodyn::~Mesodyn() {
}

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

void Mesodyn::AllocateMemory() {
  cout << "Allocating Memory.." << endl;

  //alocate memory for the fluxes of all components in all dimensions
  try {
    J.resize(dimensions * size);
    L.reserve(combinations(componentNo, 2) * Lat[0]->M);
    rho.reserve(size);
  } catch (...) {
    cout << "Failed to reserve enough memory. System too large for RAM?";
    abort();
  }
}

bool Mesodyn::CheckInput(int start) {
  if (debug)
    cout << "CheckInput in Mesodyn" << endl;
  bool success = true;

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

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  if (debug)
    cout << "mesodyn in Mesodyn." << endl;
  AllocateMemory(); //this HAS to be done before fillRho
  fillRho(initRho);
  if (success) {
    cout << "Mesodyn is all set, starting calculations.." << endl;
    for (int t = 0; t < timesteps; t++) { // get segment potentials by iteration so that they match the given rho.
      //debug = true;
      New[0]->Solve(rho[0]);
      cout << "We're back!" << endl;
      //debug = true;
      onsagerCoefficient();
      langevinFlux();
      updateDensity();
    }
  }
  return success;
}

//defaults to homogeneous system for now
void Mesodyn::fillRho(Real givenDensity) {
  if (debug)
    cout << "fillRho in Mesodyn." << endl;

  for (int i = 0; i < size; ++i) {
    rho[i] = givenDensity;
  }
  for (int i = 1; i <= componentNo; ++i) {
    ptrComponentStart.push_back(&rho[Lat[0]->M * i]);
  }
}

void Mesodyn::abort() {
  if (debug)
    cout << "abort in Mesodyn." << endl;
  //Once false is returned, Mesodyn automatically quits to main.
  success = false;
}

/******** Calculations ********/
void Mesodyn::onsagerCoefficient() {
  if (debug)
    cout << "onsagerCoefficient in Mesodyn." << endl;

  //TODO: maybe do this inline / per J calculation to preserve memory
  vector<Real>::iterator lIterator;
  lIterator = L.begin();

  //TODO: Untested code
  //all combinations of ptrComponentStart multiplications
  for (int i = 0; i < componentNo - 1; ++i) {
    for (int j = i + 1; j < componentNo; ++j) {
      //at every coordinate (pointer arithmatic)
      for (int xyz = 0; xyz < (Lat[0]->M); ++xyz) {
        *lIterator = (*(ptrComponentStart[i] + xyz) * *(ptrComponentStart[j] + xyz));
        ++lIterator;
      }
    }
  }
}

void Mesodyn::langevinFlux() {
  if (debug)
    cout << "langevinFlux in Mesodyn." << endl;

  //TODO: safer to change size calculation to xx.size()?
  vector<Real> u(size); //segment potential A

  //TODO: boundary condition in lattice?
  int z = 1;
  u[z] = New[0]->xx[z]; //which alphas? component a & b or one component at two sites?

  /*  This next part finds a vector that selects which of the combinations of onsager coefficients
      should be used for that particular component. For example for components A, B, C, D we get:
      AB-AC-AD-BC-BD-CD. For component A we need indices 0, 1, 2. For component B 0,3,4. for C 1,3,5.
  */
  vector<Real>::iterator jIterator;
  jIterator = J.begin();

  int j = 1;
  vector<int> cCombinations(componentNo - 1);
  vector<int>::iterator nIterator;
  for (int i = 0; i < componentNo; ++i) {
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

    for (int z = 0; z < Lat[0]->M; ++z) {
      gaussianNoise(mean, stdev, 1);
      int l = 0;
      //something like: for the number of onsager's coefficients, calculate flux according to x, y and z.
      *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JX]) * (u[z + i * M + JX] - u[z + i * M])) - ((L[z + cCombinations[l] * M - JX] + L[z + cCombinations[l] * M]) * (u[z + i * M] - u[z + i * M - JX])) + noise[0];
      if (componentNo > 2) {
        for (int l = 1; l < componentNo - 1; ++l) {
          *jIterator += -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JX]) * (u[z + i * M + JX] - u[z + i * M])) - ((L[z + cCombinations[l] * M - JX] + L[z + cCombinations[l] * M]) * (u[z + i * M] - u[z + i * M - JX])) + noise[0];
        }
      }
      ++jIterator;
    }
    if (dimensions > 1) {
      for (int z = 0; z < Lat[0]->M; ++z) {
        int l = 0;
        gaussianNoise(mean, stdev, 1);
        *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JY]) * (u[z + i * M + JY] - u[z + i * M])) - ((L[z + cCombinations[l] * M - JY] + L[z + cCombinations[l] * M]) * (u[z + i * M] - u[z + i * M - JY])) + noise[0];
        if (componentNo > 2) {
          for (int l = 1; l < componentNo - 1; ++l) {
            *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JY]) * (u[z + i * M + JY] - u[z + i * M])) - ((L[z + cCombinations[l] * M - JY] + L[z + cCombinations[l] * M]) * (u[z + i * M] - u[z + i * M - JY])) + noise[0];
          }
        }
        ++jIterator;
      }
    }
    if (dimensions > 2) {
      for (int z = 0; z < Lat[0]->M; ++z) {
        int l = 0;
        gaussianNoise(mean, stdev, 1);
        *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JZ]) * (u[z + i * M + JZ] - u[z + i * M])) - ((L[z + cCombinations[l] * M - JZ] + L[z + cCombinations[l] * M]) * (u[z + i * M] - u[z + i * M - JZ])) + noise[0];
        if (componentNo > 2) {
          for (int l = 1; l < componentNo - 1; ++l) {
            *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JZ]) * (u[z + i * M + JZ] - u[z + i * M])) - ((L[z + cCombinations[l] * M - JZ] + L[z + cCombinations[l] * M]) * (u[z + i * M] - u[z + i * M - JZ])) + noise[0];
          }
        }
        ++jIterator;
      }
    }
  }
}

inline Real Mesodyn::val(vector<Real>& v, int c, int x, int y, int z) {
  return v[c * Lat[0]->M + x * Lat[0]->JX + y * Lat[0]->JY + z];
}

void Mesodyn::updateDensity() {
  if (debug)
    cout << "updateDensity in Mesodyn." << endl;

  // J looks like: Jx1 - Jy1 - Jz1 - Jx2 - Jy2 - Jz2 etc.
  // Each J has size Lat[0]->M
/*
  for (int j = 0; j <= componentNo; ++j) {
    for (int i = 0; i < Lat[0]->M; ++i) {
      rho[i+j*M] += J[i       +j*M];
      rho[i+j*M] += J[i * M   +j*M];
      rho[i+j*M] += J[i * 2*M +j*M];
    }
  }
  */
}

/* Generates a vector of length count, contianing gaussian noise of given mean, standard deviation.
	 Noise is stored in vector<Real> Mesodyn::noise
	 Possible errors: What if count > sizeof(unsinged long)?
	 Called by langevinFlux()
*/
void Mesodyn::gaussianNoise(Real mean, Real stdev, unsigned int count) {
  if (debug)
    cout << "gaussianNoise in Mesodyn." << endl;

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

/******* Tools ********/
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
Real* Mesodyn::GetPointer(string s) {
  //vector<string> sub;
  //nothing yet
  return NULL;
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
