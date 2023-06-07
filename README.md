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

## TODO 30-9-2021
- [x] Fix cuda for 3d : likely problem is with generating arrays in branched propagator  
- [ ] In solve_scf: target function  
	- [ ] g(i) = phit-1/phit
	- [ ] and give \alpha weighting with \phi, may give better performance...  
- [x] Mesodyn is worked on by Daniel  
	- [ ] Grand canonical simulations (renormalise densities periodically)  
	- [ ] Multiple states  
- [ ] Layer analysis (first moments, second moment) (available in sfbox)  
- [ ] SOT evaluation Gibbs plane etcetera (partial available in sfbox)  
- [x] Dendrimers (asymmetric)   
- [x] Rings (with tether)   
- [ ] Force ensemble (partially available in sfbox)  
- [x] Grid refinement in 2D and 3D  
- [x] LBFGS implementation in newton
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
- [x] Steady states in fjc=1 and 1 gradient case  

Funtionalities SCF module summer 2022. 
Gradients: 	
	one, (classical)
	two  (cylindrical and spherical geometry)
	treee (limited system sizes due to storage of large arrays)
Geometries: planar (in 1, 2, 3 gradients)
	    cylindrical (in 1 and 2 gradients)
	    spherical (in 1 gradient only).
Markov chains: first order (direct backfolding allowed; Freely jointed chains)
	       second order (semi flexible chains @ branch points chains are usually freely jointed)
	       second order in one-gradient systems works for all grit refinement values.
Chain architecture
	linear chains 'including ' ring.
	branched chains (composition rules: parts in [ ] are side chains)
	dendrimer (symmetric and aymmetric ) use @dend(?) for help
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
Boundery conditions
	In 1 and 2 gradient calculations the surface bc is implemented
	In  3 gradient systems the surface is to be placed inside the box. 
		The idea is that this will fix the 'image'-charges problem. 
		Gradients of electrostatic potential inside solid phase can be studies by using a sufficiently thick solid phase in the system.
	Freedom of molecule can be set to 'fill-range' (when it is pinned to surface layer)
	       second order in one-gradient systems works for all grit refinement values to prevent the solvent to penetrate the surface layer.
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
	In two-gradient calculations FJC <9 can be set.
	When Markov =2 and one-gradient will work for all FJC_choices, but for 2 and 3 gradients Markov =2 works only for FJC_choices==3. 

Newton iterations
	Quasi Newton with storage of (large) Jacobian (Hessian) matrix
	DIIS wich is Hessian free storage
	Line search is a mix of strategies to prevent too large steps.
	
Computatial tricks
	memory saving (only part of the end-point distributions is stored, when needed others are recomuted...)
	local solutions solutions)
		
	
