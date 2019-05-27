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
- [ ] Dendrimers (symmetric) (available in sfbox)  
- [ ] Dendrimers (asymmetric) (partially available insfbox)  
- [ ] Comb (available in sfbox)  
- [ ] Rings (with tether) (partial available in sfbox)  
- [ ] Semi-flexibility (partially available in sfbox)  
- [ ] Force ensemble (partially available in sfbox)  
- [ ] Grid refinement in 2D and 3D  
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
