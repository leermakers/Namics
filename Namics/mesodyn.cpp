#include "mesodyn.h"

Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_, string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sys{Sys_},
      New{New_},
      D{1},
      mean{1},
      stdev{1},
      seed{1},
      timesteps{100},
      timebetweensaves{1},
      JZ{Lat[0]->JZ}, // Usage for e.g. layer z: foo[z]+foo[z+1] becomes foo[z] + foo[z+JX]
      JY{Lat[0]->JY},
      JX{Lat[0]->JX},
      MZ{Lat[0]->MZ+2},
      MY{Lat[0]->MY+2},
      MX{Lat[0]->MX+2},
      LMZ{Lat[0]->MZ},
      LMY{Lat[0]->MY},
      LMX{Lat[0]->MX},
      M{Lat[0]->M},
      componentNo{(int)Sys[0]->SysMolMonList.size()},
      cCombinations {combinations(componentNo, 2)},
      dimensions{findDimensions()},

      //These also set size == capacity. Not resizing will break all iterators.
      J(dimensions * componentNo * M),
      L(cCombinations * M),
      rho(componentNo * M),
      ptrComponentStart(componentNo),
      U(componentNo * (componentNo - 1) * M)

{
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stdev");
  if (debug)
    cout << "Mesodyn initialized" << endl;
  assert(  //Mesodyn only works if dimensions are consistently x for 1 dimension, xy for 2 dimension, xyz for 3 dimensions
       (MX > 0 && MY == 0 && MZ == 0)
    || (MX > 0 && MY > 0 && MZ == 0)
    || (MX > 0 && MY > 0 && MZ > 0)
        );
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

/* Dynamically sets the boundary condition functions to the correct type.
*/
void Mesodyn::setBoundaryPointers() {

  if (Lat[0]->BC[0] == "mirror")
    switch (dimensions) {
      case 1:
        bX0 = bind(&Mesodyn::bX0Mirror, this, 1, 1);
        break;
      case 2:
        bX0 = bind(&Mesodyn::bX0Mirror, this, MY, 1);
        break;
      case 3:
        bX0 = bind(&Mesodyn::bX0Mirror, this, MY, MZ);
        break;
    }
  else if (Lat[0]->BC[0] == "periodic")
  switch (dimensions) {
    case 1:
      bX0 = bind(&Mesodyn::bX0Periodic, this, 1, 1, MX);
      break;
    case 2:
      bX0 = bind(&Mesodyn::bX0Periodic, this, MY, 1, MX);
      break;
    case 3:
      bX0 = bind(&Mesodyn::bX0Periodic, this, MY, MZ, MX);
      break;
  }
  else
    bX0 = bind(&Mesodyn::bNothing, this);

  if (Lat[0]->BC[1] == "mirror")
  switch (dimensions) {
    case 1:
      bXm = bind(&Mesodyn::bXmMirror, this, 1, 1, MX);
      break;
    case 2:
      bXm = bind(&Mesodyn::bXmMirror, this, MY, 1, MX);
      break;
    case 3:
      bXm = bind(&Mesodyn::bXmMirror, this, MY, MZ, MX);
      break;
  }
  else if (Lat[0]->BC[1] == "periodic")
    bXm = bind(&Mesodyn::bNothing, this);
  else
    bXm = bind(&Mesodyn::bNothing, this);

  // Only true if dimensions > 1
  if (Lat[0]->BC[2] == "mirror")
  switch (dimensions) {
    case 2:
      bY0 = bind(&Mesodyn::bY0Mirror, this, MX, 1);
      break;
    case 3:
      bY0 = bind(&Mesodyn::bY0Mirror, this, MX, MZ);
      break;
  }
  else if (Lat[0]->BC[2] == "periodic")
  switch (dimensions) {
    case 2:
      bY0 = bind(&Mesodyn::bY0Periodic, this, MX, 1, MY);
      break;
    case 3:
      bY0 = bind(&Mesodyn::bY0Periodic, this, MX, MZ, MY);
      break;
  }
  else
    bY0 = bind(&Mesodyn::bNothing, this);

  // Only true if dimensions > 1
  if (Lat[0]->BC[3] == "mirror")
  switch (dimensions) {
    case 2:
      bYm = bind(&Mesodyn::bYmMirror, this, MX, 1, MY);
      break;
    case 3:
      bYm = bind(&Mesodyn::bYmMirror, this, MX, MZ, MY);
      break;
  }
  else if (Lat[0]->BC[3] == "periodic")
    bYm = bind(&Mesodyn::bNothing, this);
  else
    bYm = bind(&Mesodyn::bNothing, this);

  // Only true if dimensions > 2
  if (Lat[0]->BC[4] == "mirror")
    bZ0 = bind(&Mesodyn::bZ0Mirror, this, MX, MY);
  else if (Lat[0]->BC[4] == "periodic")
    bZ0 = bind(&Mesodyn::bZ0Periodic, this, MX, MY, MZ);
  else
    bZ0 = bind(&Mesodyn::bNothing, this);

  // Only true if dimensions > 2
  if (Lat[0]->BC[5] == "mirror")
    bZm = bind(&Mesodyn::bZmMirror, this, MX, MY, MZ);
  else if (Lat[0]->BC[5] == "periodic")
    bZm = bind(&Mesodyn::bNothing, this);
  else
    bZm = bind(&Mesodyn::bNothing, this);
}

//Used by langevinFlux() to calculate the correct number of fluxes
int Mesodyn::findDimensions() {
  int d = 0;
  if (LMX > 0)
    ++d;
  if (LMY > 0)
    ++d;
  if (LMZ > 0)
    ++d;
  return d;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  if (debug)
    cout << "mesodyn in Mesodyn" << endl;

  initRho(); // Get initial conditions (phi and potentials) by running the classical method once.
  prepareOutputFile();
//  writeRho(0); // write initial conditions
  setBoundaryPointers();

  if (debug)
    cout << "Mesodyn is all set, starting calculations.." << endl;

  int save {1};

  for (int t = 1; t < timesteps; t++) { // get segment potentials by iteration so that they match the given rho.
    cout << "MESODYN: t = " << t << endl;
    New[0]->SolveMesodyn(&rho[0]);
    onsagerCoefficient();
    potentialDifference();
    boundaryConditions();
    langevinFlux();
    if( t == timebetweensaves*save ) {
      writeRho(t);
      save+=1;
    }
    updateDensity();
  }
  return true;
}

//defaults to homogeneous system for now
void Mesodyn::initRho() {
  if (debug)
    cout << "initRho in Mesodyn." << endl;

  New[0]->Solve(true);

  //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (int i = 0; i < (int)Seg.size(); ++i) {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }

  for (int z = 0; z < Lat[0]->M; z++) {
    rho[z] = Mol[0]->phi[z];
    rho[z + M] = Mol[1]->phi[z];
  }

  //Next, register where in rho all components are located:
  //start indices for each segment type are pointed to by ptrComponentStart.
  for (int i{0}; i < componentNo; ++i) {
    ptrComponentStart.at(i) = &rho[M * i];
  }
}

void Mesodyn::updateDensity() {
  if (debug)
    cout << "updateDensity in Mesodyn." << endl;

  // J looks like: Jx1 - Jy1 - Jz1 - Jx2 - Jy2 - Jz2 etc.
  // Each J has size Lat[0]->M
  // Density = flux [lattice site * coordinate start index + component]

  for (int j = 0; j < componentNo; ++j) {
    for (int i = 0; i < M; ++i) {

      rho[i + j * M] += J[i + j * M]; // x of component j
      for (int z = 0; z < j; ++z) {
        rho[i + z * M] -= (J[i + j * M] / (componentNo-1));
      }
      for (int z = j+1; z < componentNo; ++z) {
        rho[i + z * M] -= (J[i + j * M] / (componentNo-1));
      }

      if (dimensions > 1) {
        rho[i + j * M] += J[i * M + j * M]; // y
        for (int z = 0; z < j; ++z) {
          rho[i + z * M] -= (J[i * M + j * M] / (componentNo-1));
        }
        for (int z = j+1; z < componentNo; ++z) {
          rho[i + z * M] -= (J[i * M + j * M] / (componentNo-1));
        }

        if (dimensions > 2) {
          rho[i + j * M] += J[i * 2 * M + j * M]; // z
          for (int z = 0; z < j; ++z) {
            rho[i + z * M] -= (J[i * 2 + j * M] / (componentNo-1));
          }
          for (int z = j+1; z < componentNo; ++z) {
            rho[i + z * M] -= (J[i * 2 + j * M] / (componentNo-1));
          }

        }
      }
    }
  }
}

void Mesodyn::boundaryConditions() {
  if (debug)
    cout << "boundaryConditions in Mesodyn." << endl;

  if (debug) cout << "bX0" << endl;
  bX0();
  if (debug) cout << "bXm" << endl;
  bXm();
  if (debug) cout << "bY0" << endl;
  bY0();
  if (debug) cout << "bYm" << endl;
  bYm();
  if (debug) cout << "bZ0" << endl;
  bZ0();
  if (debug) cout << "bZm" << endl;
  bZm();
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

  //for all components
  for (int i = 0; i < componentNo; ++i) {
    //for all components-1 = number of combinations possible
    for (int j = 0; j < componentNo - 1; ++j) {
      //for xyz
      for (int z = 0; z < M; ++z) {

        if (j < i) {
          //subtract the potential of the other components j or j+1 from the potential the current component i
          *uIterator = New[0]->xx[z + i * M] - New[0]->xx[z + j * M];
          ++uIterator;
        } else {
          *uIterator = New[0]->xx[z + i * M] - New[0]->xx[z + (j+1) * M];
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

  //This jIterator is going to keep track of where in the flux vector we are and keeps incrementing over all the coming loops.
  //so that fluxes are sequentially placed like Jx1, Jy1, Jz1, Jx2, Jy2, Jz2, where the number is the (indexed) lattice site.
  //(this structure is important to know when updating the densities accordingly).
  vector<Real>::iterator jIterator;

  // Zero the jIterator
  for (jIterator = J.begin(); jIterator != J.end(); ++jIterator) {
    *jIterator = 0;
  }

  jIterator = J.begin();

  /*  This next part finds a vector that selects which of the combinations of onsager coefficients
      should be used for that particular component. For example for components A, B, C, D we get:
      AB-AC-AD-BC-BD-CD. For component A (i=0) we need indices 0, 1, 2. For component B (i=1) 0,3,4. for C 1,3,5.
      It finds the correct indices, puts them in cCombinations, calculates the fluxes that correspond to that component and loops.
  */

  int j = 0;
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

    //Boundary (start)
    ++jIterator;

    // Then, we cancluate the fluxes for each dimension.
    // This assumes that if we have one dimension, it's the x dimension, whereas if we have 2, it's x & y.
    //for all lattice sites
    for (int z = 1; z < M - 1; ++z) {
      //Generate noise for this flux.
      //TODO: one noise, two noise, three noise? Different per dimension?
      gaussianNoise(mean, stdev, 1);
      //for all combinations with other components
      for (int l = 0; l < componentNo - 1; ++l) {
        *jIterator += -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JX]) * (U[z + (i*M*(componentNo-1)) + l * M + JX] - U[z + (i*M*(componentNo-1)) + l * M])) - ((L[z + cCombinations[l] * M - JX] + L[z + cCombinations[l] * M]) * (U[z + (i*M*(componentNo-1)) + l * M] - U[z + (i*M*(componentNo-1)) + l * M - JX]));
      }
      ++jIterator;
    }

    //Boundary (end)
    ++jIterator;

    if (dimensions > 1) {

      //Boundary (start)
      ++jIterator;

      for (int z = 1; z < M - 1; ++z) {
        gaussianNoise(mean, stdev, 1);
        for (int l = 0; l < componentNo - 1; ++l) {
          *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JY]) * (U[z + (i*M*(componentNo-1)) + l * M + JY] - U[z + (i*M*(componentNo-1)) + l * M])) - ((L[z + cCombinations[l] * M - JY] + L[z + cCombinations[l] * M]) * (U[z + (i*M*(componentNo-1)) + l * M] - U[z + (i*M*(componentNo-1)) + l * M - JY]));
        }
        ++jIterator;
      }

      //Boundary (end)
      ++jIterator;

      if (dimensions > 2) {

        //Boundary (start)
        ++jIterator;

        for (int z = 1; z < M - 1; ++z) {
          gaussianNoise(mean, stdev, 1);
          for (int l = 0; l < componentNo - 1; ++l) {
            *jIterator = -D * ((L[z + cCombinations[l] * M] + L[z + cCombinations[l] * M + JZ]) * (U[z + (i*M*(componentNo-1)) + l * M + JZ] - U[z + (i*M*(componentNo-1)) + l * M])) - ((L[z + cCombinations[l] * M - JZ] + L[z + cCombinations[l] * M]) * (U[z + (i*M*(componentNo-1)) + l * M] - U[z + (i*M*(componentNo-1)) + l * M - JZ]));
          }
          ++jIterator;
        }

        //Boundary (end)
        ++jIterator;
      }
    }
    // And loop!
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

/******* Boundary conditions *******/

//MX-1 and MX-2 because the last index is MX-1 (starts at 0)

void Mesodyn::bX0Mirror(int fMY, int fMZ) {
  if (debug)
    cout << "bX0Mirror in Mesodyn" << endl;

  for (int y = 0; y < fMY; ++y) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, 0, y, z) = val(L, c, 1, y, z); //start
      }

      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, 0, y, z) = val(U, c, 1, y, z); //start
      }
    }
  }
}

void Mesodyn::bXmMirror(int fMY, int fMZ, int fMX) {
  if (debug)
    cout << "bXmMirror in Mesodyn" << endl;

  for (int y = 0; y < fMY; ++y) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, fMX - 1, y, z) = val(L, c, fMX - 2, y, z); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        if (c < componentNo - 1)
          *valPtr(U, c, fMX - 1, y, z) = val(U, c, fMX - 2, y, z); //end
      }
    }
  }
}

void Mesodyn::bX0XmMirror(int fMY, int fMZ, int fMX) {
  if (debug)
    cout << "bX0XmMirror in Mesodyn" << endl;

  for (int y = 0; y < fMY; ++y) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, 0, y, z) = val(L, c, 1, y, z);           //start
        *valPtr(L, c, fMX - 1, y, z) = val(L, c, fMX - 2, y, z); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        if (c < componentNo - 1)
          *valPtr(U, c, 0, y, z) = val(U, c, 1, y, z);         //start
        *valPtr(U, c, fMX - 1, y, z) = val(U, c, fMX - 2, y, z); //end
      }
    }
  }
}

void Mesodyn::bY0Mirror(int fMX, int fMZ) {
  if (debug)
    cout << "bY0Mirror in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, 0, z) = val(L, c, x, 1, z); //start
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        if (c < componentNo - 1)
          *valPtr(U, c, x, 0, z) = val(U, c, x, 1, z); //start
      }
    }
  }
}

void Mesodyn::bYmMirror(int fMX, int fMZ, int fMY) {
  if (debug)
    cout << "bYmMirror in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, fMY - 1, z) = val(L, c, x, fMY - 2, z); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, fMY - 1, z) = val(U, c, x, fMY - 2, z); //end
      }
    }
  }
}

void Mesodyn::bY0YmMirror(int fMX, int fMZ, int fMY) {
  if (debug)
    cout << "bY0YmMirror in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, 0, z) = val(L, c, x, 1, z);           //start
        *valPtr(L, c, x, fMY - 1, z) = val(L, c, x, fMY - 2, z); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, 0, z) = val(U, c, x, 1, z);           //start
        *valPtr(U, c, x, fMY - 1, z) = val(U, c, x, fMY - 2, z); //end
      }
    }
  }
}

void Mesodyn::bZ0Mirror(int fMX, int fMY) {
  if (debug)
    cout << "bZ0Mirror in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int y = 0; y < fMY; ++y) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, y, 0) = val(L, c, x, y, 1); //start
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, y, 0) = val(U, c, x, y, 1); //start
      }
    }
  }
}

void Mesodyn::bZmMirror(int fMX, int fMY, int fMZ) {
  if (debug)
    cout << "bZmMirror in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int y = 0; y < fMY; ++y) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, y, fMZ - 1) = val(L, c, x, y, fMZ - 2); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, y, fMZ - 1) = val(U, c, x, y, fMZ - 2); //end
      }
    }
  }
}

void Mesodyn::bZ0ZmMirror(int fMX, int fMY, int fMZ) {
  if (debug)
    cout << "bZ0ZmMirror in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int y = 0; y < fMY; ++y) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, y, 0) = val(L, c, x, y, 1);           //start
        *valPtr(L, c, x, y, fMZ - 1) = val(L, c, x, y, fMZ - 2); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, y, 0) = val(U, c, x, y, 1);           //start
        *valPtr(U, c, x, y, fMZ - 1) = val(U, c, x, y, fMZ - 2); //end
      }
    }
  }
}

void Mesodyn::bX0Periodic(int fMY, int fMZ, int fMX) {
  if (debug)
    cout << "bX0Periodic in Mesodyn" << endl;

  for (int y = 0; y < fMY; ++y) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, 0, y, z) = val(L, c, fMX - 2, y, z); //start
        *valPtr(L, c, fMX - 1, y, z) = val(L, c, 1, y, z); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, 0, y, z) = val(U, c, fMX - 2, y, z); //start
        *valPtr(U, c, fMX - 1, y, z) = val(U, c, 1, y, z); //end
      }
    }
  }
}

void Mesodyn::bY0Periodic(int fMX, int fMZ, int fMY) {
  if (debug)
    cout << "bY0Periodic in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int z = 0; z < fMZ; ++z) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, 0, z) = val(L, c, x, fMY - 2, z); //start
        *valPtr(L, c, x, fMY - 1, z) = val(L, c, x, 1, z); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, 0, z) = val(U, c, x, fMY - 2, z); //start
        *valPtr(U, c, x, fMY - 1, z) = val(U, c, x, 1, z); //end
      }
    }
  }
}

void Mesodyn::bZ0Periodic(int fMX, int fMY, int fMZ) {
  if (debug)
    cout << "bZ0Periodic in Mesodyn" << endl;

  for (int x = 0; x < fMX; ++x) {
    for (int y = 0; y < fMY; ++y) {
      for (int c = 0; c < cCombinations; ++c) {
        //Onsager coefficents
        *valPtr(L, c, x, y, 0) = val(L, c, x, y, fMZ - 2); //start
        *valPtr(L, c, x, y, fMZ - 1) = val(L, c, x, y, 1); //end
      }
      for (int c = 0; c < componentNo; ++c) {
        //Potentials
        *valPtr(U, c, x, y, 0) = val(U, c, x, y, fMZ - 2); //start
        *valPtr(U, c, x, y, fMZ - 1) = val(U, c, x, y, 1); //end
      }
    }
  }
}

void Mesodyn::bNothing() {
  if (debug)
    cout << "bNothing in Mesodyn" << endl;
  //do nothing
}

/******* Output generation *******/

void Mesodyn::prepareOutputFile() {
  /* Open filestream and set filename to "mesodyn-datetime.csv" */
  ostringstream filename;
  filename << "output/mesodyn-";

  time_t rawtime;
  time(&rawtime);
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

  for (int i = 1; i <= cCombinations; ++i) {
    headers << "L" << i << ",";
  }

  for (int i = 1; i <= componentNo; ++i) {
    headers << "U" << i << ",";
  }

  for (int i = 1; i <= componentNo; ++i) {
    headers << "J" << i << ",";
  }

  headers << "\n";

  mesFile << headers.str();
}

void Mesodyn::writeRho(int t) {
  ostringstream timeOutput;
  timeOutput << t << ",\n";
  mesFile << timeOutput.str();

  ostringstream rhoOutput;

  int tLMZ {LMZ};
  int tLMY {LMY};
  int tLMX {LMX};

  switch ( dimensions ) {
    case 1:
      tLMY = 1;
      tLMZ = 1;
      break;
    case 2:
      tLMZ = 1;
      break;
    default:
      break;

  }

  for (int z = 0; z < tLMZ; ++z) {
    for (int y = 0; y < tLMY; ++y) {
      for (int x = 0; x < tLMX; ++x) {
        rhoOutput << x << ",";//":" << y << ":" << z << ",";
        for (int c = 0; c < componentNo; ++c) {
          rhoOutput << val(rho, c, x, y, z) << ",";
        }

        for (int c = 0; c < cCombinations; ++c) {
          rhoOutput << val(L, c, x, y, z) << ",";
        }

        for (int c = 0; c < componentNo; ++c) {
          rhoOutput << val(U, c, x, y, z) << ",";
        }

        for (int c = 0; c < componentNo; ++c) {
          rhoOutput << val(J, c, x, y, z) << ",";
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

inline Real* Mesodyn::valPtr(vector<Real>& v, int c, int x, int y, int z) {
  return &v[c * Lat[0]->M + x * Lat[0]->JX + y * Lat[0]->JY + z];
}

int Mesodyn::factorial(int n) {
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
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
