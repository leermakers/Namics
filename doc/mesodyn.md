# Mesodyn

## Introduction

The mesodyn classes are an implementation of dynamic mean-field theory, that was developed by Fraaije in the 1990s. It allows the user to generate the time evolution of multi-component system given the parameters used by the classic Namics SCF tools. Mesodyn uses inheritance to handle different dimensions of the system.

## Dependencies
### System
Provides a list - called a mask - (variable `KSAM`) of where immovable objects are located (such as pinned or frozen molecules). This list is used to initialize the density profile and to disable fluxes to and from lattice sites with such objects. Moreover, the number of components in the system is read from `SysMolMonList`.

### Lattice
Contains all information on the geometry of the system and passes the boundary conditions to mesodyn contained in its `BC` vector. Mesodyn implements an interface to the lattice class called Lattice_access. The interface provides easy and consistent access to the lattice using 3 coordinates instead of an indexed vector.

## Variables read from file

The following inputs can be read from file using `Mesodyn : langevin : variable : value`

| C++ Variable         | Variable in input | Description    |
| -------------------- | ----------------- | -------------- |
| Real D               | diffusionconstant | Diffusion constant for the langevin flux |
| Real mean            | mean              | Mean of gaussian noise (should be 0)     |
| Real stdev           | stdev             | Standard deviation of gaussian noise (should be 1*D) |
| Real seed            | seed              | Seed of gaussian noise, for reproducable simulations including noise |
| int timesteps        | timesteps         | Length scale of the time evolution |
| int timebetweensaves | timebetweensaves  | How many timesteps before mesodyn writes the current variables to file |

## Class definitions

### Mesodyn

###### Introduction
Mesodyn is the book-keeper class for the implementation. It handles input and initialization, output, and the overall structure of the mesodyn algorithm. Variables that are not read from file are initialized in the initialization list.

###### Functions
`CheckInput(int start)` <br />
First of two functions called from `main()`. Reads input from the .in file. Contrary to most variables that live in the mesodyn class, variables read from file are not marked `const`. If a variable is not specified in the .in file, the variable is given a default value. A developer adding a new input-readable variable in this function should also add them to `prepareOutputFile()`.
Argument `start` is passed from main via this function to the `Input` class and is not used by mesodyn.


`bool mesodyn()` <br />
Called from `main()` after `CheckInput()`. Provides the overall flow of the mesodyn algorithm and thus only calls other functions. It calls the functions that implement the following capabilities:

1. Initialize density profile and construct helper classes `Flux` and `Component`
2. Open and format output file
3. Start loop:
  1. Prepare an indexed rho vector for Newton
  2. Call Newton to find potentials that correspond to the density profile
  3. Update the boundaries given the boundary conditions
  4. Compute fluxes
  5. Update densities
  6. Loop for given number of timesteps and write data to file for every given number of timesteps
  7. Return true when everything goes according to plan.

`int initial_conditions()` <br />
Sets up the two helper classes, `Flux` and `Component` and passes the data needed for initialization. First, the mask (see system dependency) is read from system. The density profile is initialized using `init_rho(rho, mask)` (see below). It then constructs two vectors of size component number that contain pointers to the respective class instances of `Flux` and `Component`. Depending on the number of dimensions, these are instances of either the 1D, 2D or 3D version of the classes. This removes the need for lots of `if (dimension)` statements throughout the code. The boundary conditions are read from input or assigned a default value by `Lattice`. Called by `mesodyn()`.

`int init_rho(vector<vector<Real>>& rho, vector<int>& mask)` <br />
Computes the volume of the system (number of lattice sites minus solid objects) using the argument `mask` and distributes the theta (sum of all densities, read from input by `System`) over the remaining lattice sites. It also finds which component is the solvent and accounts for it in the density profile. The density profile is loaded into the (reference) argument rho, the first dimension being the component number, the second being the lattice site. Called by `initial_conditions()`.

`void prepareOutputFile()`  <br />
Opens a file in /output and prepares it with headers for density profile, langevin coefficient, potentials, potential differences, fluxes. Also writes all file-readable variables at the top of the file. A developer adding a new input-readable variable in `CheckInput` should also add them to this function. Called by `mesodyn()`.

`void writeRho(int t)` <br />
Outputs the data described in `prepareOutputFile()` for a given timestep. Called by `mesodyn()`.


### Component

###### Introduction

The component class is used to keep track of variables and properties that are specified per one component. These include the density profiles `rho` and potentials `alpha`. Furthermore it implements functions that form the interface to these variables for `Mesodyn`. Because boundary conditions are applied to the potentials and could be set per component, these also live in this class. Each derived class houses the boundary conditions for the dimension that they correspond to (X for 1D, Y for 2D, Z for 3D). Lattice_Access asserts that this is indeed the case. Boundary conditions are implemented using the `std::bind` function. This allows the class to assign the function pointers for updating the boundaries to the correct functions and fix the arguments. This design removes the need for multiple `if (dimension)` statements.

###### Inheritance

![Component inheritance](images/component_inheritance.png)

###### Functions

### Flux

###### Introduction

The flux class is used to keep track of variables and functions that are specified per pair of components. These include the fluxes `J` (and its two directions `J_plus` and `J_minus`) and the variables needed to compute them: the Onsager coefficient `L` and the chemical potential difference `mu`. Each instance of the class is assigned a pair of components for which it does its calculations. Furthermore it implements functions that form the interface to these variables for `Mesodyn`. Each derived class houses the flux (and the function to calculate it) for the dimension that they correspond to (X for 1D, Y for 2D, Z for 3D). Lattice_Access asserts that the latter is indeed the case. This design removes the need for multiple `if (dimension)` statements.

###### Inheritance

![Component inheritance](images/flux_inheritance.png)

###### Functions

## References
