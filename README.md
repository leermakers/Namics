## Namics: 

Hybrid SCF-MD simulation tool. Package under development.

Developers: F.A.M.Leermakers, R.Varadharajan.

TODO 21-8-2018
-electrostatics: 
	generalize for variable epsilon. done: testing
-fix cuda for 3d : likely problem is with generating arrays in branched propagator
-in solve_scf: target function
	-g(i) = phit-1/phit and give \alpha weighting with \phi, may give better performance...
-mesodyn is worked on by Daniel
	grand canonical simulations (renormalise densities periodically) 
	multiple states
-multiple states in namics
	thermodynamics done: testing
	restart (available in sfbox) done: testing
	variate (available in sfbox) done : testing
	output of all kinds of computables in reaction.cpp (available in sfbox) done: testing
-layer analysis (first moments, second moment) (available in sfbox)
-SOT evaluation Gibbs plane etcetera (partial available in sfbox)
-dendrimers (symmetric) (available in sfbox)
-dendrimers (asymmetric) (partially available insfbox)
-comb (available in sfbox)
-rings (with tether) (partial available in sfbox)
-semi-flexibility (partially available in sfbox)
-force ensemble (partially available in sfbox)
-grit refinement in 2D and 3D
-LBFGS implementation in newton (available in sfbox)
-MGRES implementation in newton
-Trunctated newton  (available in sfbox)
-fix Picard method
-BD-SCF hybrid
	fixed mobility
	Quasi-Newton mobility
-Cleng (with MC) is underway by Sacha
	dynamics
-Teng (with MC) is underway by Ram
	dynamics
-steady states   
	1D
	2D 
	3D

