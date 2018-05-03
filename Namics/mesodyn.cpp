#include "mesodyn.h"

//Constructor
Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_, string name_) {
  In = In_;
  name = name_;
  Lat = Lat_;
  Mol = Mol_;
  Seg = Seg_;
  Sys = Sys_;
  New = New_;
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  xNeighbors.reserve(2);
  yNeighbors.reserve(2);
  zNeighbors.reserve(2);
  componentNo = 0;
}
Mesodyn::~Mesodyn() {
}

void Mesodyn::AllocateMemory() {
  if (debug)
    cout << "nothing to allocate in Mesodyn" << endl;
}

bool Mesodyn::CheckInput(int start) {
  if (debug)
    cout << "Check Mesodyn" << endl;
  bool success = true;
  //right now all the checks for engine are done in input.cpp. could be move back here.
  success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);
  if (success) {
    vector<string> options;
    if (GetValue("timesteps").size() > 0) {
      success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 10000, "The number of timesteps should be between 1 and 10000");
    }
    cout << "timesteps is " << timesteps << endl;
  }
  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  //TODO: How to get here from main.
  pepareForCalculations();
  if (success) {
    for (int t = 0; t < timesteps; t++) {
      New[0]->Solve(&dummyVector[0], &dummyVector[0]);
      langevinFlux(dummyVector, dummyVector, dummyVector, dummyVector);
      updateDensity();
    }
  }
  return success;
}

void Mesodyn::pepareForCalculations() {
  componentNo = findComponentNo();    //find how many compontents there are (e.g. head, tail, solvent)
  size = findDensityVectorSize();     //find out how large the density vector is (needed for sizing the flux vector)
                                      //which will be 1 flux per lattice site per component per dimension

  //find the index ranges where each component / dimension is located in the vector that
  //contains phi's and u's (which is a 3D lattice mapped onto 1D vector (Row-major))
  if (componentNo <= 1) {
      cout << "WARNING: Not enough components found for Mesodyn, aborting!";
      abort();
  } else if (componentNo > 1 && componentNo < 3) {
      setNeighborIndices(xNeighbors, yNeighbors, zNeighbors);
      setComponentStartIndices(component);
  } else {
      cout << "Unable to do Mesodyn for " << componentNo << " components, aborting!";
      abort();
  }
  //alocate memory for the fluxes of all components in all dimensions
  J.resize(size);
}

void Mesodyn::abort() {
  //Once false is returned, Mesodyn automatically quits to main.
  success = false;
}

/****** Functions for handling indices in a 1D vector that contains a 3D lattice with multiple components *******/

int Mesodyn::findComponentNo() {
  component.reserve(In[0]->MolList.size());
  return In[0]->MolList.size();
}

int Mesodyn::findDensityVectorSize() {
  //  TODO: boundary conditions
  return componentNo * Lat[0]->MX * Lat[0]->MY * Lat[0]->MZ;
}

void Mesodyn::setNeighborIndices(vector<int>& xNeighbors, vector<int>& yNeighbors, vector<int>& zNeighbors) {
  //  TODO: boundary conditions
  xNeighbors[0] = -1;
  xNeighbors[1] = 1;

  yNeighbors[0] = -Lat[0]->MX;
  yNeighbors[1] = Lat[0]->MX;

  zNeighbors[0] = -Lat[0]->MX * Lat[0]->MY;
  zNeighbors[1] = Lat[0]->MX * Lat[0]->MY;
}

void Mesodyn::setComponentStartIndices(vector<int>& component) {
  if (componentNo == 2) {
    component[0] = 0;
    component[1] = Lat[0]->MX * Lat[0]->MY * Lat[0]->MZ;
  }
  //unless, for whatever reason componentNo changed
  else {
    cout << "Component number changed, this should definitely not have happened and the programmer messed up.";
    abort();
  }
}


/******** Calculations ********/

//Two components
void Mesodyn::langevinFlux(vector<Real>& phiA, vector<Real>& phiB, vector<Real>& alphaA, vector<Real>& alphaB) {
  vector<Real> L(size); //onsager coefficient
  vector<Real> u(size); //segment potential A

  //TODO: boundary condition in lattice?

  for (int z = 1; z < size; z++) {
    L[z] = phiA[z] * phiA[z];
  }

  for (int z = 1; z < size; z++) {
    gaussianNoise(dummyMean, dummyStdev, 2);
    //segment "chemical potential" gradient
    u[z] = alphaA[z] - alphaB[z];
    J[z] = (-D * (((L[z] + L[z + 1]) * (u[z + 1] - u[z])) - ((L[z - 1] + L[z]) * (u[z] - u[z - 1]))) + noise[0]);
    J[size + z] = (-D * (((L[z] + L[z + 1]) * (-u[z + 1] - -u[z])) - ((L[z - 1] + L[z]) * (-u[z] - -u[z - 1]))) + noise[1]);
  }
}

void Mesodyn::updateDensity() {
  //old density + langevinFluxTwo
}

/* Generates a vector of length count, contianing gaussian noise of given mean, standard deviation.
	 Noise is stored in vector<Real> Mesodyn::noise
	 Possible errors: What if count > sizeof(unsinged long)?
	 Called by langevinFlux()
*/
void Mesodyn::gaussianNoise(Real mean, Real stdev, unsigned long count) {

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
