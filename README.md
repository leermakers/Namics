## Namics:

Hybrid SCF-MD simulation tool. Package under development.

Developers:   
F.A.M.Leermakers - Self-consistent field caculation modules.   
R.Varadharajan - teng.cpp   
Daniel Emmery - mesodyn.cpp   
Alexander Kazakov - cleng.cpp   

### Dependencies:

GPU Accelleration:
- SCF tested to work on:
	- CUDA 9.0 & g++5
	- CUDA 9.1 & g++5
	- CUDA 9.2 & g++7
	- CUDA 10.0 & g++7
	
- NOT working:
	- CUDA 10.1 & g++8
	- CUDA 9.0 & g++ 6
	- CUDA 9.2 & g++ 5

Set the ccbin value in the NVCC flags in the makefile to the correct g++ version and replace the CUDA paths if needed. Also set the nvcc arch flag to the correct compute capability (list can be found [here](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/))

## TODO 9-4-2018
- [x] Fix cuda for 3d : likely problem is with generating arrays in branched propagator  
- [ ] In solve_scf: target function  
	- [ ] g(i) = phit-1/phit
	- [ ] and give \alpha weighting with \phi, may give better performance...  
- [x] Mesodyn is worked on by Daniel  
	- [ ] Grand canonical simulations (renormalise densities periodically)  
	- [ ] Multiple states  
- [ ] Layer analysis (first moments, second moment) (available in sfbox)  
- [ ] SOT evaluation Gibbs plane etcetera (partial available in sfbox)  
- [x] Dendrimers (symmetric) (available in sfbox)  
- [ ] Dendrimers (asymmetric) (partially available insfbox)  
- [x] Comb (available in sfbox)  
- [ ] Rings (with tether) (partial available in sfbox)  
- [x] Semi-flexibility (partially available in sfbox)  
- [ ] Force ensemble (partially available in sfbox)  
- [x] Grid refinement in 2D and 3D  
- [ ] LBFGS implementation in newton (available in sfbox)  
- [ ] MGRES implementation in newton  
- [ ] Trunctated newton  (available in sfbox)  
- [ ] Fix Picard method  
- [ ] BD-SCF hybrid  
	- [ ] Fixed mobility  
	- [ ] Quasi-Newton mobility  
- [ ] Cleng (with MC) is underway by Sacha  
	- [ ] Dynamics  
- [ ] Teng (with MC) is underway by Ram  
	- [ ] Dynamics  
- [ ] Steady states  

Funtionalities SCF module summer 2021. 
Gradients: 	
	one, (classical)
	two  (cylindrical and spherical geometry)
	treee (limited system sizes due to storage of large arrays)
Geometries: planar (in 1, 2, 3 gradients)
	    cylindrical (in 1 and 2 gradients)
	    spherical (in 1 gradient only).
Markov chains: first order (direct backfolding allowed; Freely jointed chains)
	       second order (semi flexible chains @ branch points chains are usually freely jointed)
Chain architecture
	linear chains
	branched chains (parts in [ ] are side chains)
	regular dendrimer (symmetric)
	combs
Interations
	Flory Huggins nearest neighbour interactions
	electrostatic interactions
		weak 
		strong
	fixed surface potential
	fixed surface charge
	Incompressible limit
Constraints
	Pinning constraints 
	Frozen (spectator segments with local density unity)
	Clamping (constrained at both ends)
	beta constraint (extra constraint to e.g. fix the position of an interface). 
Output	
	Density profiles
	Free energy 
		Grand potential
		chemical potentials
Lattice type 
	simple_cubic
	hexagonal 

Lattice discretisation
	FJC-choices 3: segment size equal to lattice side (classical)
	In one-gradient calculations the FJC value can be increased to 5, 7, 9, etc (lattice refinement, quasi lattice-free)
	In two-gradient calculations FJC = 5 can be set.
	When Markov =2, FJC must be FJC =3. 

Newton iterations
	Quasi Newton with storage of (large) Jacobian (Hessian) matrix
	DIIS wich is Hessian free storage
	Line search is a mix of strategies to prevent too large steps.
	
Computatial tricks
	memory saving (only part of the end-point distributions is stored, when needed others are recomuted...)
	local solutions (in 3d box to find solutions)
		
	
