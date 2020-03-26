#include "system.h"
#include "tools.h"


System::System(vector<Input *> In_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<State *> Sta_, vector<Reaction *> Rea_, vector<Molecule *> Mol_, string name_)
{
	Seg = Seg_;
	Mol = Mol_;
	Lat = Lat_;
	In = In_;
	name = name_;
	Sta = Sta_;
	Rea = Rea_;
	prepared = false;
	if (debug)
		cout << "Constructor for system " << endl;
  KEYS.push_back("calculation_type");
	KEYS.push_back("constraint");
	KEYS.push_back("delta_range");
	KEYS.push_back("delta_inputfile");
	KEYS.push_back("delta_molecules");
	KEYS.push_back("phi_ratio");
	KEYS.push_back("generate_guess");
	KEYS.push_back("initial_guess");
	KEYS.push_back("guess_inputfile");
	KEYS.push_back("final_guess");
	KEYS.push_back("guess_outputfile");
	KEYS.push_back("GPU");
	KEYS.push_back("X");
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++)
	  KEYS.push_back("guess-" + In[0]->MonList[i]);
	charged=false;
	constraintfields=false;
  boundaryless_volume=0;
	grad_epsilon = false;
}
System::~System()
{
	if (debug)
		cout << "Destructor for system " << endl;
	free(H_GrandPotentialDensity);
	free(H_FreeEnergyDensity);
	free(H_alpha);
	if (charged)
	{
		free(H_q);
		free(H_psi);
	}
	if (constraintfields)
	{
		free(H_beta);
		free(H_BETA);
	}
#ifdef CUDA
  cudaFree(phitot);
  cudaFree(alpha);
  cudaFree(GrandPotentialDensity);
  cudaFree(FreeEnergyDensity);
  cudaFree(TEMP);
  cudaFree(KSAM);
  if (charged) {
    cudaFree(psi);
    cudaFree(eps);
    cudaFree(q);
    cudaFree(EE);
    cudaFree(E);
    cudaFree(psiMask);
  }
  if (constraintfields) {
	cudaFree(beta);
	cudaFree(BETA);
  }
#else
  free(phitot);
  free(TEMP);
  free(KSAM);
  free(CHI);
  if (charged) {
    free(EE);
    free(E);
    free(psiMask);
    free(eps);
  }

#endif
}

void System::AllocateMemory()
{
	if (debug)
		cout << "AllocateMemory in system " << endl;
	int M = Lat[0]->M;
	//H_GN_A = new Real[n_box];
	//H_GN_B = new Real[n_box];
	H_GrandPotentialDensity = (Real *)malloc(M * sizeof(Real));
	std::fill(H_GrandPotentialDensity,H_GrandPotentialDensity+M,0);

	H_FreeEnergyDensity = (Real *)malloc(M * sizeof(Real));
	H_alpha = (Real *)malloc(M * sizeof(Real));
	if (charged)
	{
		H_q = (Real *)malloc(M * sizeof(Real));
		H_psi = (Real *)malloc(M * sizeof(Real));
	}
	if (constraintfields)
	{
		H_beta = (int *)malloc(M * sizeof(int));
		std::fill(H_beta,H_beta+M,0);
		H_BETA = (Real *)malloc(M * sizeof(Real));
		Lat[0]->FillMask(H_beta, px, py, pz, delta_inputfile);
	}

#ifdef CUDA
	phitot = (Real *)AllOnDev(M);
	alpha = (Real *)AllOnDev(M);
	GrandPotentialDensity = (Real *)AllOnDev(M);
	FreeEnergyDensity = (Real *)AllOnDev(M);
	TEMP = (Real *)AllOnDev(M);
	KSAM = (int *)AllIntOnDev(M);
	Zero(KSAM, M);
  if (charged) {
    psi = (Real*)AllOnDev(M);
    q = (Real*)AllOnDev(M);
    eps = (Real*)AllOnDev(M);
    EE = (Real*)AllOnDev(M);
     E = (Real*)AllOnDev(M);
    psiMask = (int*)AllIntOnDev(M);
  }
  if (constraintfields) {
	BETA = (Real*)AllOnDev(M);
	beta = (int*)AllIntOnDev(M);
  }
#else
  phitot = (Real*)malloc(M * sizeof(Real));
  alpha = H_alpha;
  if (charged) {
    psi = H_psi;
    q = H_q;
    eps = (Real*)malloc(M * sizeof(Real));
    EE = (Real*)malloc(M * sizeof(Real));
     E = (Real*)malloc(M * sizeof(Real));
    psiMask = (int*)malloc(M * sizeof(int));
  }
  if (constraintfields) {
	beta=H_beta;
	BETA=H_BETA;
  }
  KSAM = (int*)malloc(M * sizeof(int));
  FreeEnergyDensity = H_FreeEnergyDensity;
  GrandPotentialDensity = H_GrandPotentialDensity;
  TEMP = (Real*)malloc(M * sizeof(Real));
#endif
  Zero(KSAM, M);
  if (charged) {
    Zero(psi, M);
    Zero(EE, M);
    Zero(E,M);
  }

	n_mol = In[0]->MolList.size();
	Lat[0]->AllocateMemory();
	int n_mon = In[0]->MonList.size();
	for (int i = 0; i < n_mon; ++i)
		Seg[i]->AllocateMemory();
	for (int i = 0; i < n_mol; ++i)
		Mol[i]->AllocateMemory();
}

bool System::generate_mask()
{
	int M = Lat[0]->M;
	bool success = true;

	FrozenList.clear();
	int length = In[0]->MonList.size();
	for (int i = 0; i < length; i++)
	{
		if (Seg[i]->freedom == "frozen")
			FrozenList.push_back(i);
	}
	if (FrozenList.size() + SysMonList.size() + SysTagList.size() + SysClampList.size() != In[0]->MonList.size())
	{
		cout << " Thera are un-used monomers in system. Remove them before starting" << endl;
		success = false;
	}

	Zero(KSAM, M);

	length = FrozenList.size();
	for (int i = 0; i < length; ++i)
	{
		Add(KSAM, Seg[FrozenList[i]]->MASK, M);
	}

	length = SysTagList.size();
	for (int i = 0; i < length; ++i)
	{
		Add(KSAM, Seg[SysTagList[i]]->MASK, M);
	}

	length = SysClampList.size();
	for (int i = 0; i < length; ++i)
	{
		Add(KSAM, Seg[SysClampList[i]]->MASK, M);
	}

	Invert(KSAM, KSAM, M);

	if (Lat[0]->gradients < 3)
	{
		for (int i = 0; i < M; i++)
		{
			volume += KSAM[i] * Lat[0]->L[i];
		}
	}
	else
	{
		Sum(volume, KSAM, M);
	}

	// Used to initialize densities in mesodyn.
	// I know it's hideous, but it works for all dimensions.
	// If you really care about legibility you could do this with a switch or object.
	this->boundaryless_volume = this->volume - ((2 * Lat[0]->gradients - 4) * Lat[0]->MX * Lat[0]->MY + 2 * Lat[0]->MX * Lat[0]->MZ + 2 * Lat[0]->MY * Lat[0]->MZ + (-2 + 2 * Lat[0]->gradients) * (Lat[0]->MX + Lat[0]->MY + Lat[0]->MZ) + pow(2, Lat[0]->gradients));

	Lat[0]->Accesible_volume=volume;

	return success;
}

bool System::PrepareForCalculations(bool first_time)
{
	if (debug)
		cout << "PrepareForCalculations in System " << endl;

	bool success = true;
	int M = Lat[0]->M;

	// necessary part; essentially for cleng
	if (In.back()->MesodynList.empty() or prepared == false) {
		success = generate_mask();
		prepared = true;
	}

	if (constraintfields)
	{
		Boltzmann(BETA, BETA, M);
		//for(int i=0; i<M; i++) cout << "BETA at i: " << i << " is: " << BETA[i] << endl;
		//cin.get();
	}

	n_mol = In[0]->MolList.size();
	success = Lat[0]->PrepareForCalculations();
	int n_mon = In[0]->MonList.size();

	for (int i = 0; i < n_mon; i++)
	{
		success = Seg[i]->PrepareForCalculations(KSAM,first_time);
	}
	for (int i = 0; i < n_mol; i++)
	{
		success = Mol[i]->PrepareForCalculations(KSAM);
	}
#ifdef CUDA
	if (constraintfields)
	{
		//TransferDataToDevice(H_BETA, BETA, M);
		TransferDataToDevice(H_beta, beta, M);
	}
#endif
	if (first_time) {
  if (charged) {
    int length = FrozenList.size();
    Zero(psiMask, M);
    fixedPsi0 = false;
    for (int i = 0; i < length; ++i) {
      if (Seg[FrozenList[i]]->fixedPsi0) {
        fixedPsi0 = true;
        Add(psiMask, Seg[FrozenList[i]]->MASK, M);
      }
    }
    Real eps=Seg[0]->epsilon;
    for (int i=1; i<n_mon; i++) {
	if (Seg[i]->epsilon != eps) grad_epsilon = true;
    }
  }
  if (initial_guess == "polymer_adsorption") {
	  int length=FrozenList.size();
	  for (int i=0; i<length; i++){
		  for (int j=0; j<n_mon; j++) {if (!CHI[i+n_mon*j] == 0) Seg[j]->PutAdsorptionGuess(CHI[i+n_mon*j],Seg[i]->MASK);
		  }
	  }
  }
	}

  return success;
}

string System::GetMonName(int mon_number_I)
{
	if (debug)
		cout << "GetMonName for system " << endl;
	return Seg[mon_number_I]->name;
}

bool System::CheckInput(int start)
{
	if (debug)
		cout << "CheckInput for system " << endl;
	bool success = true;
	bool solvent_found = false;
	tag_segment = -1;
	solvent = -1; //value -1 means no solvent defined. tag_segment=-1;
	Real phibulktot = 0;
	success = In[0]->CheckParameters("sys", name, start, KEYS, PARAMETERS, VALUES);
	if (success)
	{
		success = CheckChi_values(In[0]->MonList.size());
		GPU = In[0]->Get_bool(GetValue("GPU"), false);
		if (GPU)
			if (!cuda)
			{
				cout << "You should compile the program using the CUDA=1 flag: GPU calculations are impossible; proceed with CPU computations..." << endl;
				GPU = false;
			}
		if (Lat[0]->gradients < 3)
		{
			if (GPU)
				cout << "GPU support is (for the time being) only available for three-gradient calculations " << endl;
		}
		if (cuda)
		{
			if (!GPU)
			{
				cout << "Please enable GPU in the input file (SYS : BRAND : GPU : true)." << endl;
				success = false;
			}
		}

		int length = In[0]->MonList.size();
		for (int i = 0; i < length; i++)
			if (Seg[i]->state_name.size() == 0)
				StatelessMonList.push_back(i);
		length = In[0]->MolList.size();
		int statelength = In[0]->StateList.size();
		int i = 0;
		while (i < length)
		{
			int j = 0;
			int LENGTH = Mol[i]->MolMonList.size();
			while (j < LENGTH)
			{
				SysMolMonList.push_back(Mol[i]->MolMonList[j]);
				if (!In[0]->InSet(SysMonList, Mol[i]->MolMonList[j]))
				{
					if (Seg[Mol[i]->MolMonList[j]]->freedom != "tagged" && Seg[Mol[i]->MolMonList[j]]->freedom != "clamp")
					{
						SysMonList.push_back(Mol[i]->MolMonList[j]);
						if (Seg[Mol[i]->MolMonList[j]]->state_name.size() < 1 && IsUnique(Mol[i]->MolMonList[j], -1))
						{
							ItMonList.push_back(Mol[i]->MolMonList[j]);
						}
					}
				}
				if (Seg[Mol[i]->MolMonList[j]]->freedom == "tagged")
				{
					if (In[0]->InSet(SysTagList, Mol[i]->MolMonList[j]))
					{
						//cout <<"You can not use the 'tag monomer' " + GetMonName(Mol[i]->MolMonList[j]) + " in more than one molecule." << endl; success=false;
					}
					else
						SysTagList.push_back(Mol[i]->MolMonList[j]);
				}
				if (Seg[Mol[i]->MolMonList[j]]->freedom == "clamp")
				{
					if (In[0]->InSet(SysClampList, Mol[i]->MolMonList[j]))
					{
						//cout <<"You can not use the 'clamp monomer' " + GetMonName(Mol[i]->MolMonList[j]) + " in more than one molecule." << endl; success=false;
					}
					else
						SysClampList.push_back(Mol[i]->MolMonList[j]);
				}
				j++;
			}

			if (Mol[i]->freedom == "free")
				phibulktot += Mol[i]->phibulk;
			if (Mol[i]->freedom == "solvent")
			{
				solvent_found = true;
				solvent = i;
			}
			if (Mol[i]->IsTagged())
			{
				tag_segment = Mol[i]->tag_segment;
				Mol[i]->n = 1.0 * Seg[tag_segment]->n_pos;
				Mol[i]->theta = Mol[i]->n * Mol[i]->chainlength;
			}
			i++;
		}
		for (int j = 0; j < statelength; j++)
		{
			if (IsUnique(-1, j))
			{
				ItStateList.push_back(j);
				//cout <<"itstatelist extended with " << j << endl;
			}
		}
		if (!solvent_found || (phibulktot > 0.99999999 && phibulktot < 1.0000000001))
		{
			cout << "In system '" + name + "' the 'solvent' was not found, while the volume fractions of the bulk do not add up to unity. " << endl;
			success = false;
		}
		if (IsCharged())
		{
			charged = true;
			neutralizer = -1;
			bool neutralizer_needed = false;

			int length = In[0]->MolList.size();
			for (int i = 0; i < length; i++)
				if (Mol[i]->freedom == "neutralizer")
					neutralizer = i;

			if (neutralizer < 0)
			{
				Real phibulk = 0;
				for (int i = 0; i < length; i++)
				{
					if (Mol[i]->freedom == "free")
						phibulk += Mol[i]->Charge() * Mol[i]->phibulk;
					if (Mol[i]->freedom == "restricted" && Mol[i]->IsCharged() && !Mol[i]->IsPinned())
						neutralizer_needed = true;
				}

				if (neutralizer_needed || abs(phibulk) > 1e-4)
				{
					cout << "Neutralizer needed because some un-pinned molecules have freedom 'restricted' and/or overall charge density of 'free' molecules in the bulk is not equal to zero " << endl;
					success = false;
				}
			}
		}

		if (GetValue("constraint").size() > 0)
		{
			constraintfields = true;
			px.clear();
			py.clear();
			pz.clear();
			vector<string> constraints;
			constraints.push_back("delta");
			ConstraintType = "";
			if (!In[0]->Get_string(GetValue("constraint"), ConstraintType, constraints, "Info about 'constraint' rejected"))
			{
				success = false;
			};
			if (ConstraintType == "delta")
			{
				if (GetValue("delta_range").size() > 0)
				{
					string s = GetValue("delta_range");
					vector<string> sub;
					vector<string> set;
					vector<string> coor;
					In[0]->split(s, ';', sub);
					int n_points = sub.size();
					for (int i = 0; i < n_points; i++)
					{
						set.clear();
						In[0]->split(sub[i], '(', set);
						int length = set.size();
						if (length != 2)
						{
							if (length == 1 && set[0] == "file")
							{
								if (GetValue("delta_inputfile").size() > 0)
								{
									delta_inputfile = GetValue("delta_inputfile");
								}
								else
								{
									success = false;
									cout << "When 'delta_range' is set to 'file', you should provide a 'delta_inputfile'" << endl;
								}
							}
							else
							{
								success = false;
								cout << "In 'delta_range', for position " << i << "the expected '(x,y,z)' format was not found " << endl;
							}
						}
						else
						{
							coor.clear();
							In[0]->split(set[1], ',', coor);
							int grad = Lat[0]->gradients;
							int corsize = coor.size();
							if (corsize != grad)
							{
								success = false;
								if (grad == 1)
									cout << "In 'delta_range', for position " << i << " the expected '(x)' format was not found " << endl;
								if (grad == 2)
									cout << "In 'delta_range', for position " << i << " the expected '(x,y)' format was not found " << endl;
								if (grad == 3)
									cout << "In 'delta_range', for position " << i << " the expected '(x,y,z)' format was not found " << endl;
							}
							else
							{
								px.push_back(In[0]->Get_int(coor[0], -1));
								if (grad > 1)
									py.push_back(In[0]->Get_int(coor[1], -1));
								if (grad > 2)
									pz.push_back(In[0]->Get_int(coor[2], -1));
							}
						}
					}
				}
				else
				{
					success = false;
					cout << "When 'constraint' is set to 'delta', you should specify a 'delta_range' " << endl;
				}

				if (GetValue("delta_molecules").size() > 0)
				{
					string deltamols = GetValue("delta_molecules");
					vector<string> sub;
					In[0]->split(deltamols, ';', sub);
					int length_sub = sub.size();
					if (length_sub != 2)
					{
						success = false;
						cout << " delta_molecules item should contain two 'molecule names' separated by a ';'" << endl;
					}
					else
					{
						DeltaMolList.clear();
						int length = In[0]->MolList.size();
						for (int i = 0; i < length; i++)
						{
							if (sub[0] == Mol[i]->name)
								DeltaMolList.push_back(i);
							if (sub[1] == Mol[i]->name)
								DeltaMolList.push_back(i);
						}
						if (DeltaMolList[0] == DeltaMolList[1])
						{
							success = false;
							cout << " In delta_molecules you should specify two different names " << endl;
						}
						if (DeltaMolList.size() < 2)
						{
							success = false;
							cout << " In 'delta_molecules', one or more molecule names are not recognized" << endl;
						}
					}
				}
				else
				{
					success = false;
					cout << "When 'constraint' is set to 'delta', you should specify a set of 'delta_molecules' " << endl;
				}


				phi_ratio=1.0;
				if(GetValue("phi_ratio").size()>0){phi_ratio=In[0]->Get_Real(GetValue("phi_ratio"),1);}
			}

			//if (Mol[DeltaMolList[0]]->freedom=="restricted" || Mol[DeltaMolList[1]]->freedom=="restricted" ) {success =false;  cout <<"Molecule in list of delta_molecules has not freedom 'free'"<<endl; }
		}


		vector<string> options;

		options.push_back("micro_phasesegregation");
		CalculationType = "";
		if (GetValue("calculation_type").size() > 0)
		{
			if (!In[0]->Get_string(GetValue("calculation_type"), CalculationType, options, " Info about calculation_type rejected"))
			{
				success = false;
			};
			if (CalculationType == "micro_emulsion")
			{
				if (In[0]->MolList.size() < 3)
				{
					cout << "In 'calculation_type : micro-emulsion', we expect at least three types of molecules " << endl;
					success = false;
				}
				int length = In[0]->MolList.size();
				int i = 0;
				int number_of_solvents = 0;
				int number_of_surfactants = 0;
				while (i < length)
				{
					if (Mol[i]->IsTagged())
					{
						number_of_surfactants++;
					}
					else
					{
						number_of_solvents++;
					}
					i++;
				}
				if (number_of_solvents < 2 || number_of_surfactants < 1)
				{
					cout << "In 'calculation_type : micro-emulsion', we expect at least two solvents and one surfactant which is tagged. " << endl;
					success = false;
				}
			}
			if (CalculationType == "micro_phasesegregation")
			{
				int length = In[0]->MolList.size();
				bool found = false;
				int i = 0;
				while (i < length && !found)
				{
					if (Mol[i]->MolMonList.size() > 1)
						found = true;
				}
				if (!found)
				{
					cout << "In 'calculation_type : micro_phasesegregation', we expect at least copolymer in the system " << endl;
					success = false;
				}
			}
		}
		options.clear();
		options.push_back("lamellae");
		options.push_back("Im3m");
		GuessType = "";
		MonA = 0;
		MonB = 0;
		if (GetValue("generate_guess").size() > 0)
		{
			if (!In[0]->Get_string(GetValue("generate_guess"), GuessType, options, " Info about 'generate_guess' rejected"))
			{
				success = false;
			};
			if (GuessType == "lamellae" || GuessType == "Im3m")
			{
				int length = In[0]->MonList.size();
				int n_guess = 0;
				for (int i = 0; i < length; i++)
				{
					string s = "guess-" + In[0]->MonList[i];
					if (GetValue(s).size() > 0)
					{
						Seg[i]->guess_u = In[0]->Get_Real(GetValue(s), 0);
						if (Seg[i]->guess_u < -2 || Seg[i]->guess_u > 2)
						{
							cout << "Suggested 'guess value' for 'u' for mon '" + In[0]->MonList[i] + "' out of range: current range is -2, .., 2. Value ignored. " << endl;
							Seg[i]->guess_u = 0;
						}
						else
						{
							n_guess++;
							if (i == 0 && n_guess == 1)
								MonA = i;
							if (i == 1 && n_guess == 2)
								MonB = i;
						}
					}
				}
				if (n_guess < 2)
				{
					cout << "generate_guess needs two valid non-zero (real) values for 'guess_xxx' quantities. Here xxx is a monomomer name; 'generate_guess : " + GuessType + "' is most likely ineffective.  " << endl;
					cout << "The idea is that the first 'mon' quantity will be the 'oil' and the second 'mon' quantity the 'water' component. Each filling one half-space. Adjust input accordingly. " << endl;
				}
			}
		}
		initial_guess = "previous_result";
		if (GetValue("initial_guess").size() > 0)
		{
			options.clear();
			options.push_back("previous_result");
			options.push_back("file");
			options.push_back("polymer_adsorption");
			In[0]->Get_string(GetValue("initial_guess"), initial_guess, options, " Info about 'initial_guess' rejected;");
			if (initial_guess == "file")
			{
				if (GetValue("guess_inputfile").size() > 0)
				{
					guess_inputfile = GetValue("guess_inputfile");
				}
				else
				{
					success = false;
					cout << " When 'initial_guess' is set to 'file', you need to supply 'guess_inputfile', but this entry is missing. Problem terminated " << endl;
				}
			}
			if (initial_guess=="polymer_adsorption") {
				//cout <<"guess for polymer adsorption not implemented" << endl;
			}
		}
		final_guess = "next_problem";
		if (GetValue("final_guess").size() > 0)
		{
			options.clear();
			options.push_back("next_problem");
			options.push_back("file");
			if (!In[0]->Get_string(GetValue("final_guess"), final_guess, options, " Info about 'final_guess' rejected; default: 'next_problem' used."))
			{
				final_guess = "next_problem";
			}
			if (final_guess == "file")
			{
				if (GetValue("guess_outputfile").size() > 0)
				{
					guess_outputfile = GetValue("guess_outputfile");
				}
				else
				{
					guess_outputfile = "";
					cout << "Filename not found for 'output_guess'. Default with inputfilename and extention '.outiv' is used. " << endl;
				}
			}
		}
	}

	internal_states = false;

	if (In[0]->StateList.size() > 1)
	{
		internal_states = true;
		int num_of_Seg_with_states = 0;
		int num_of_Eqns = In[0]->ReactionList.size();
		int num_of_alphabulk_fixed = 0;
		int num_of_states = In[0]->StateList.size();
		for (int k = 0; k < num_of_states; k++)
			if (Sta[k]->fixed)
				num_of_alphabulk_fixed++;
		int length = In[0]->MonList.size();
		for (int k = 0; k < length; k++)
			if (Seg[k]->state_name.size() > 1)
				num_of_Seg_with_states++;
		if (num_of_Seg_with_states + num_of_Eqns + num_of_alphabulk_fixed != num_of_states)
		{
			cout << " num_of_Seg_with_states+num_of_Eqns+num_of_alphabulk_fixed !=num_of_states" << endl;
			if (num_of_alphabulk_fixed == 0)
			{
				cout << " Consider to define for one of the states an alphabulk value " << endl;
			}
			else
			{
				if (num_of_alphabulk_fixed > 1)
				{
					cout << " possibly you have specified too many alphabulk values for multiple states " << endl;
				}
				else
				{
					if ((num_of_Seg_with_states + num_of_Eqns + num_of_alphabulk_fixed > num_of_states))
					{
						cout << " Possibly you have defined too many equations ... " << endl;
					}
					else
					{
						cout << " Possibly you have defined too few equations ... " << endl;
					}
				}
			}
			success = false;
		}
	}

	if (GetValue("X").size() > 0)
	{
		XmolList.clear();
		XstateList_1.clear();
		XstateList_2.clear();
		Xn_1.clear();
		string s = GetValue("X");
		vector<string> sub;
		vector<string> SUB;
		In[0]->split(s, '-', sub);
		if (sub[0] != "F")
		{
			cout << "X is the characteristic function specified by user." << endl;
			cout << "Example of how a characteristic function is defined:" << endl;
			cout << "Case 1: no internal states: 'F - molname_1  -molname_2 - ...' " << endl;
			cout << "Here molname_1 etc are names of molecules in the system. " << endl;
			cout << "Case 2: internal states : F - molname_1 -(statename_1,statename_2,#number) - ... " << endl;
			cout << "Note that in ... you can add as many molname's and (...,...,...) combinations as you wish" << endl;
			cout << "Here F = Helmholtz energy. " << endl;
			cout << "and the statement '-molname_i' implies that 'n_i times mu_i' is subtracted from F. " << endl;
			cout << "The (statename_i,statename_j,n) implies that 'n times theta_i times mu_j' is subtracted from F " << endl;
			cout << endl;
			cout << "Error found: The first item is not the expected 'F' " << endl;
			success = false;
		}
		else
		{
			int length_sub = sub.size();
			int length_mol = In[0]->MolList.size();
			int length_state = In[0]->StateList.size();
			for (int i = 1; i < length_sub; i++)
			{
				SUB.clear();
				In[0]->split(sub[i], ',', SUB);
				if (SUB.size() == 1)
				{ //want to see mol name
					bool found = false;
					for (int k = 0; k < length_mol; k++)
					{
						if (SUB[0] == In[0]->MolList[k])
						{
							XmolList.push_back(k);
							found = true;
						}
					}
					if (!found)
					{
						success = false;
						cout << "In characteristic function X, the entry '" + SUB[0] + "' is not a molecule name " << endl;
					}
				}
				else
				{ //want to see (statename1,statename2)
					if (SUB.size() != 3)
					{
						cout << "In characteristic function X, the entry '" + sub[i] + "' is not recognised as '(statename_1, statename_2,n_1 )' " << endl;
						success = false;
					}
					else
					{
						bool found_1 = false, found_2 = false;
						for (int k = 0; k < length_state; k++)
						{
							if (SUB[0].substr(1, SUB[0].length() - 1) == In[0]->StateList[k])
							{
								found_1 = true;
								XstateList_1.push_back(k);
							}
							if (SUB[1] == In[0]->StateList[k])
							{
								found_2 = true;
								XstateList_2.push_back(k);
							}
						}
						int sto = In[0]->Get_int(SUB[2].substr(0, SUB[2].length() - 1), -1);
						if (sto < 0)
						{
							success = false;
							cout << "In characteristic function X, the entry '" + sub[i] + "' does not include a positive integer at the third argument. " << endl;
						}
						else
							Xn_1.push_back(sto);
						if (!(found_1 && found_2))
						{
							cout << "In characteristic function X, the entry '" + sub[i] + "' is not coding for '(statename_1,statename_2,n_1)" << endl;
							if (!found_1)
								cout << "first state name " + SUB[0].substr(1, SUB[0].length() - 1) + " does not exist" << endl;
							if (!found_2)
								cout << "second state name " + SUB[1] + " does not exist" << endl;
							success = false;
						}
					}
				}
			}
			//cout <<"this is what I have found" << endl;
			//int length_1= XmolList.size();
			//int length_2= Xn_1.size();
			//cout <<"the molecules that are subtracted are " << endl;
			//for (int i=0; i<length_1; i++) cout <<XmolList[i] << " "; cout <<endl;
			//cout <<"the states that are subtracted are " << endl;
			//for (int i=0; i<length_2; i++) cout <<XstateList_1[i] <<" " <<XstateList_2[i]<< " " << Xn_1[i] << endl;
		}
	}
	return success;
}

bool System::IsUnique(int Segnr_, int Statenr_)
{
	if (debug)
		cout << "System::IsUnique: Segnr = " << Segnr_ << " Statenr = " << Statenr_ << endl;
	bool is_unique = true;
	bool is_equal = true;
	if (In[0]->MesodynList.size() > 0)
		return is_unique;
	int Segnr = Segnr_, Statenr = Statenr_;
	int length = 0;
	if (Segnr > -1)
		length = Seg[Segnr]->chi.size();
	else
		length = Sta[Statenr]->chi.size();
	int itmonlength = ItMonList.size();
	int itstatelength = ItStateList.size();
	if (Statenr < 0)
		for (int i = 0; i < itmonlength; i++)
		{
			if (is_unique)
			{
				for (int k = 0; k < length; k++)
					if (Seg[Segnr]->chi[k] != Seg[ItMonList[i]]->chi[k])
						is_equal = false;
				if (is_equal)
				{
					is_unique = false;
					Seg[Segnr]->unique = false;
					Seg[Segnr]->seg_nr_of_copy = ItMonList[i];
				}
				else
					is_equal = true;
			}
		}
	else
	{
		for (int i = 0; i < itmonlength; i++)
		{
			if (is_unique)
			{
				for (int k = 0; k < length; k++)
					if (Sta[Statenr]->chi[k] != Seg[ItMonList[i]]->chi[k])
						is_equal = false;
				if (is_equal)
				{
					is_unique = false;
					Sta[Statenr]->unique = false;
					Sta[Statenr]->seg_nr_of_copy = ItMonList[i];
				}
				else
					is_equal = true;
			}
		}
		for (int i = 0; i < itstatelength; i++)
		{
			if (is_unique)
			{
				for (int k = 0; k < length; k++)
					if (Sta[Statenr]->chi[k] != Sta[ItStateList[i]]->chi[k])
						is_equal = false;
				if (is_equal)
				{
					is_unique = false;
					Sta[Statenr]->unique = false;
					Sta[Statenr]->state_nr_of_copy = ItStateList[i];
				}
				else
					is_equal = true;
			}
		}
	}
	return is_unique;
}

bool System::PutVarInfo(string Var_type_, string Var_target_, Real Var_target_value_)
{
	if (debug)
		cout << "System::PutVarInfo " << endl;
	bool success = true;
	Var_target = -1;
	if (Var_type_ != "target")
		success = false;
	if (Var_target_ == "free_energy")
		Var_target = 0;
	if (Var_target_ == "grand_potential")
	{
		Var_target = 1;
	}
	if (Var_target < 0 || Var_target > 1)
	{
		success = false;
		cout << "Var target " + Var_target_ + " rejected in PutVarInfo in System " << endl;
	}
	Var_target_value = Var_target_value_;
	if (Var_target_value < -1e4 || Var_target_value > 1e4)
		success = false;
	return success;
}

Real System::GetError()
{
	if (debug)
		cout << "System::GetError " << endl;
	Real Error = 0;
	switch (Var_target)
	{
	case 0:
		Error = FreeEnergy - Var_target_value;
		break;
	case 1:
		Error = -1.0 * (GrandPotential - Var_target_value);
		break;
	default:
		cout << "Program error in GetVarError" << endl;
		break;
	}
	return Error;
}

bool System::IsCharged()
{
	if (debug)
		cout << "System::IsCharged " << endl;
	bool success = false;
	int length = In[0]->MolList.size();
	for (int i = 0; i < length; i++)
	{
		if (Mol[i]->IsCharged())
			success = true;
	}
	return success;
}

void System::PutParameter(string new_param)
{
	if (debug)
		cout << "PutParameter for system " << endl;
	KEYS.push_back(new_param);
}

string System::GetValue(string parameter)
{
	if (debug)
		cout << "GetValue " + parameter + " for system " << endl;
	int length = PARAMETERS.size();
	for (int i = 0; i < length; ++i)
	{
		if (parameter == PARAMETERS[i])
		{
			return VALUES[i];
		}
	}
	return "";
}

void System::push(string s, Real X)
{
	if (debug)
		cout << "push (Real) for system " << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void System::push(string s, int X)
{
	if (debug)
		cout << "push (int) for system " << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void System::push(string s, bool X)
{
	if (debug)
		cout << "push (bool) for system " << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void System::push(string s, string X)
{
	if (debug)
		cout << "push (string) for system " << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}
void System::PushOutput()
{
	if (debug)
		cout << "PushOutput for system " << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("e", e);
	push("k_B", k_B);
	push("eps0", eps0);
	push("temperature", T);
	push("free_energy", FreeEnergy);
	push("grand_potential", GrandPotential);
	if (Lat[0]->gradients == 1 && Lat[0]->geometry == "planar")
	{
		push("KJ0", -Lat[0]->Moment(GrandPotentialDensity, 1));
		push("Kbar", Lat[0]->Moment(GrandPotentialDensity, 2));
	}
	Real X = 0;
	if (Xn_1.size() > 0 || XmolList.size() > 0)
	{
		X = FreeEnergy;
		int length_mol = XmolList.size();
		for (int i = 0; i < length_mol; i++)
		{
			X -= Mol[XmolList[i]]->n * Mol[XmolList[i]]->Mu;
		}
		int length_state = Xn_1.size();
		Real mu = -999;
		for (int i = 0; i < length_state; i++)
		{
			length_mol = In[0]->MolList.size();
			for (int k = 0; k < length_mol; k++)
			{
				if (Mol[k]->chainlength == 1)
				{
					int seg = Mol[k]->MolMonList[0];
					if (Seg[seg]->ns > 1)
					{
						for (int j = 0; j < Seg[seg]->ns; j++)
						{
							if (Seg[seg]->state_name[j] == In[0]->StateList[XstateList_2[i]])
							{
								mu = Mol[k]->mu_state[j];
							}
						}
					}
				}
			}
			if (mu == -999)
			{
				cout << "Failed to find chemical potential for state " + In[0]->StateList[XstateList_2[i]] + ": (not a monomer?) In characteristic function X, mu is set to zero." << endl;
				mu = 0;
			}
			X -= Seg[Sta[XstateList_1[i]]->mon_nr]->state_theta[Sta[XstateList_1[i]]->state_nr] * Xn_1[i] * mu;
		}
		push("X", X);
		cout << " X  = " << X << endl;
	}
	push("calculation_type", CalculationType);
	push("guess_type", GuessType);
	push("cuda", cuda);
	push("solvent", Mol[solvent]->name);

	string s = "profile;0";
	push("alpha", s);
	s = "profile;1";
	push("GrandPotentialDensity", s);
	s = "profile;2";
	push("FreeEnergyDensity", s);
	if (charged)
	{
		s = "profile;3";
		push("psi", s);
		s = "profile;4";
		push("q", s);
		s = "profile;5";
		push("eps", s);
	}
#ifdef CUDA
	int M = Lat[0]->M;
	TransferDataToHost(H_alpha, alpha, M);
	TransferDataToHost(H_GrandPotentialDensity, GrandPotentialDensity, M);
	TransferDataToHost(H_FreeEnergyDensity, FreeEnergyDensity, M);
#endif
}

Real *System::GetPointer(string s, int &SIZE)
{
	if (debug)
		cout << "GetPointer for system " << endl;
	vector<string> sub;
	SIZE = Lat[0]->M;
	In[0]->split(s, ';', sub);
	if (sub[1] == "0")
		return H_alpha;
	if (sub[1] == "1")
		return H_GrandPotentialDensity;
	if (sub[1] == "2")
		return H_FreeEnergyDensity;
	if (sub[1] == "3")
		return psi;
	if (sub[1] == "4")
		return q;
	if (sub[1] == "5")
		return eps;
	return NULL;
}
int *System::GetPointerInt(string s, int &SIZE)
{
	if (debug)
		cout << "GetPointerInt for system " << endl;
	vector<string> sub;
	In[0]->split(s, ';', sub);
	if (sub[0] == "array")
	{ //set SIZE and return pointer of int array
	}
	return NULL;
}

int System::GetValue(string prop, int &int_result, Real &Real_result, string &string_result)
{
	if (debug)
		cout << "GetValue (long) for system " << endl;
	int length = ints.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == ints[i])
		{
			int_result = ints_value[i];
			return 1;
		}
	}
	length = Reals.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == Reals[i])
		{
			Real_result = Reals_value[i];
			return 2;
		}
	}
	length = bools.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == bools[i])
		{
			if (bools_value[i])
				string_result = "true";
			else
				string_result = "false";
			return 3;
		}
	}
	length = strings.size();
	for (int i = 0; i < length; ++i)
	{
		if (prop == strings[i])
		{
			string_result = strings_value[i];
			return 3;
		}
	}
	return 0;
}

bool System::CheckChi_values(int n_seg)
{
	if (debug)
		cout << "CheckChi_values for system " << endl;
	bool success = true;
	CHI = (Real *)malloc(n_seg * n_seg * sizeof(Real));
	for (int i = 0; i < n_seg; i++)
		for (int k = 0; k < n_seg; k++)
		{
			CHI[i * n_seg + k] = In[0]->Get_Real(Seg[i]->GetValue("chi-" + Seg[k]->name), 123);
		}
	for (int i = 0; i < n_seg; i++)
		for (int k = 0; k < n_seg; k++)
			if (CHI[i * n_seg + k] == 123)
				CHI[i * n_seg + k] = CHI[k * n_seg + i];
	for (int i = 0; i < n_seg; i++)
		for (int k = 0; k < n_seg; k++)
			if (CHI[i * n_seg + k] == 123)
				CHI[i * n_seg + k] = 0;
	for (int i = 0; i < n_seg; i++)
		for (int k = i + 1; k < n_seg; k++)
			if (CHI[i * n_seg + k] != CHI[k * n_seg + i])
			{
				cout << "CHI-value symmetry violated: chi(" << Seg[i]->name << "," << Seg[k]->name << ") is not equal to chi(" << Seg[k]->name << "," << Seg[i]->name << ")" << endl;
				success = false;
			}
	for (int i = 0; i < n_seg; i++)
	{
		if (CHI[i * n_seg + i] != 0)
		{
			cout << "CHI-values for 'same'-segments (e.g. CHI(x,x) ) should be zero. " << endl;
			success = false;
		}
	}

	int n_segments = In[0]->MonList.size();
	int n_states = In[0]->StateList.size();
	if (n_states == 1)
		n_states = 0;
	int n_chi = n_segments + n_states;

	//for (int i=0; i<n_segments; i++)
	//{for (int k=0; k<n_chi; k++) cout <<Seg[i]-> chi[k] << " "; cout << endl; }
	//for (int i=0; i<n_states; i++)
	//{for (int k=0; k<n_chi; k++) cout <<Sta[i]-> chi[k] << " "; cout << endl; }

	for (int i = 0; i < n_segments; i++)
		for (int j = 0; j < n_segments; j++)
		{
			if (Seg[i]->chi[j] == -999 && Seg[j]->chi[i] == -999)
			{
				Seg[i]->chi[j] = Seg[j]->chi[i] = 0;
			}
			else
			{
				if (Seg[i]->chi[j] != -999 && Seg[j]->chi[i] == -999)
				{
					Seg[j]->chi[i] = Seg[i]->chi[j];
				}
				else
				{
					if (Seg[i]->chi[j] == -999 && Seg[j]->chi[i] != -999)
					{
						Seg[i]->chi[j] = Seg[j]->chi[i];
					}
					else
					{
						if (Seg[i]->chi[j] != Seg[j]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Seg[i]->name << "," << Seg[j]->name << ") != chi (" << Seg[j]->name << "," << Seg[i]->name << ")" << endl;
						}
					}
				}
			}
		}

	for (int i = n_segments; i < n_chi; i++)
		for (int j = n_segments; j < n_chi; j++)
		{
			if (Sta[i - n_segments]->chi[j] == -999 && Sta[j - n_segments]->chi[i] == -999)
			{
				Sta[i - n_segments]->chi[j] = Sta[j - n_segments]->chi[i] = Seg[Sta[i - n_segments]->mon_nr]->chi[Sta[j - n_segments]->mon_nr];
			}
			else
			{
				if (Sta[i - n_segments]->chi[j] != -999 && Sta[j - n_segments]->chi[i] == -999)
				{
					Sta[j - n_segments]->chi[i] = Sta[i - n_segments]->chi[j];
				}
				else
				{
					if (Sta[i - n_segments]->chi[j] == -999 && Sta[j - n_segments]->chi[i] != -999)
					{
						Sta[i - n_segments]->chi[j] = Sta[j - n_segments]->chi[i];
					}
					else
					{
						if (Sta[i - n_segments]->chi[j] != Sta[j - n_segments]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Sta[i - n_segments]->name << "," << Sta[j - n_segments]->name << ") != chi (" << Sta[j - n_segments]->name << "," << Sta[i - n_segments]->name << ")" << endl;
						}
					}
				}
			}
		}

	for (int i = 0; i < n_segments; i++)
		for (int j = n_segments; j < n_chi; j++)
		{
			if (Seg[i]->chi[j] == -999 && Sta[j - n_segments]->chi[i] == -999)
			{
				Seg[i]->chi[j] = Sta[j - n_segments]->chi[i] = Seg[i]->chi[Sta[j - n_segments]->mon_nr];
			}
			else
			{
				if (Seg[i]->chi[j] != -999 && Sta[j - n_segments]->chi[i] == -999)
				{
					Sta[j - n_segments]->chi[i] = Seg[i]->chi[j];
				}
				else
				{
					if (Seg[i]->chi[j] == -999 && Sta[j - n_segments]->chi[i] != -999)
					{
						Seg[i]->chi[j] = Sta[j - n_segments]->chi[i];
					}
					else
					{
						if (Seg[i]->chi[j] != Sta[j - n_segments]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Seg[i]->name << "," << Sta[j - n_segments]->name << ") != chi (" << Sta[j - n_segments]->name << "," << Seg[i]->name << ")" << endl;
						}
					}
				}
			}
		}

	for (int i = n_segments; i < n_chi; i++)
		for (int j = 0; j < n_segments; j++)
		{
			if (Sta[i - n_segments]->chi[j] == -999 && Seg[j]->chi[i] == -999)
			{
				Sta[i - n_segments]->chi[j] = Seg[j]->chi[i] = Seg[j]->chi[Sta[i - n_segments]->mon_nr];
			}
			else
			{
				if (Sta[i - n_segments]->chi[j] != -999 && Seg[j]->chi[i] == -999)
				{
					Seg[j]->chi[i] = Sta[i - n_segments]->chi[j];
				}
				else
				{
					if (Sta[i - n_segments]->chi[j] == -999 && Seg[j]->chi[i] != -999)
					{
						Sta[i - n_segments]->chi[j] = Seg[j]->chi[i];
					}
					else
					{
						if (Sta[i - n_segments]->chi[j] != Seg[j]->chi[i])
						{
							success = false;
							cout << " conflict in chi values! chi(" << Sta[i - n_segments]->name << "," << Seg[j]->name << ") != chi (" << Seg[j]->name << "," << Sta[i - n_segments]->name << ")" << endl;
						}
					}
				}
			}
			if (Sta[i - n_segments]->mon_nr == j && Seg[j]->chi[i] != 0 && Seg[j]->chi[i] != -999)
			{
				success = false;
				cout << " chi between mon-type and one of its states is not allowed for chi(" << Sta[i - n_segments]->name << "," << Seg[j]->name << ")" << endl;
			}
		}

	return success;
}

void System::DoElectrostatics(Real *g, Real *x)
{
	int M = Lat[0]->M;
	int n_seg = In[0]->MonList.size();
	Zero(q, M);
	Zero(eps, M);
	for (int i = 0; i < n_seg; i++)
	{
		if (Seg[i]->ns < 2)
		{
			YplusisCtimesX(q, Seg[i]->phi, Seg[i]->valence, M);
		}
		YplusisCtimesX(eps, Seg[i]->phi, Seg[i]->epsilon, M);
	}
	int statelistlength = In[0]->StateList.size();
	for (int i = 0; i < statelistlength; i++)
	{
		YplusisCtimesX(q, Seg[Sta[i]->mon_nr]->phi_state + Sta[i]->state_nr * M, Sta[i]->valence, M);
	}

	Div(q, phitot, M);
	Div(eps, phitot, M);
	Lat[0]->set_bounds(eps);
	Cp(psi,x,M); Cp(g,psi,M);
	Lat[0]->set_bounds(psi);
	if (fixedPsi0) {
		int length=FrozenList.size();
		for (int i=0; i<length; i++) {
			Seg[FrozenList[i]]->UpdateValence(g,psi,q,eps,grad_epsilon);
		}
	}
}

bool System::ComputePhis(){
if(debug) cout <<"ComputePhis in system" << endl;
	int M= Lat[0]->M;
	Real A=0, B=0; //A should contain sum_phi*charge; B should contain sum_phi
	bool success=true;
	Zero(phitot,M);
	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		Real *phi_frozen=Seg[FrozenList[i]]->phi;
		Add(phitot,phi_frozen,M);
	}

	for (int i=0; i<n_mol; i++) {
		if (constraintfields) {
			if (i==DeltaMolList[0]) {
				success=Mol[i]->ComputePhi(BETA,1);
			} else {
				if (i==DeltaMolList[1]) {
					success=Mol[i]->ComputePhi(BETA,-1);
				} else {
					success=Mol[i]->ComputePhi(BETA,0);
				}
			}
		}
		else
			success = Mol[i]->ComputePhi(BETA, 0);
	}

	for (int i = 0; i < n_mol; i++)
	{
		Real norm = 0;
		if (Mol[i]->freedom == "free")
		{
			norm = Mol[i]->phibulk / Mol[i]->chainlength;
			Mol[i]->n = norm * Mol[i]->GN;

			A += Mol[i]->phibulk * Mol[i]->Charge();
			B += Mol[i]->phibulk;
		}

		if (Mol[i]->freedom == "restricted")
		{
			if (Mol[i]->GN > 0)
			{
				norm = Mol[i]->n / Mol[i]->GN;
				if (Mol[i]->IsPinned())
				{
					Mol[i]->phibulk = 0;
				}
				else
				{
					Mol[i]->phibulk = Mol[i]->chainlength * norm;
					A += Mol[i]->phibulk * Mol[i]->Charge();
					B += Mol[i]->phibulk;
				}
			}
			else
			{
				norm = 0;
				cout << "GN for molecule " << i << " is not larger than zero..." << endl;
				throw - 1;
			}
		}

		if (Mol[i]->IsTagged() || Mol[i]->IsPinned())
		{
			if (Mol[i]->GN > 0)
				norm = Mol[i]->n / Mol[i]->GN;
			else
			{
				norm = 0;
				cout << "GN for molecule " << i << " is not larger than zero..." << endl;
			}
			Mol[i]->phibulk = 0;
		}
		if (Mol[i]->IsClamped())
		{
			norm = 1;
			Mol[i]->phibulk = 0;
		}

		int k = 0;
		Mol[i]->norm = norm;
		int length = Mol[i]->MolMonList.size();

		while (k < length)
		{
			if (!(Seg[Mol[i]->MolMonList[k]]->freedom == "clamp" || Mol[i]->freedom == "frozen"))
			{
				Real *phi = Mol[i]->phi + k * M;
				//Real *G1=Mol[i]->G1+k*M;
				Real *G1 = Seg[Mol[i]->MolMonList[k]]->G1;
				Div(phi, G1, M);
				if (norm > 0)
					Norm(phi, norm, M);

				if (debug)
				{
					Real sum;
					Sum(sum, phi, M);
					cout << "Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl;
				}
			}
			k++;
		}

		if (Mol[i]->freedom == "range_restricted")
		{
			Real *phit = Mol[i]->phitot;
			Zero(phit, M);
			int k = 0;
			while (k < length)
			{
				Real *phi = Mol[i]->phi + k * M;
				Add(phit, phi, M);
				k++;
			}
			OverwriteA(phit, Mol[i]->R_mask, phit, M);
			//Lat[0]->remove_bounds(phit);
			Real theta = Lat[0]->ComputeTheta(phit);
			norm = Mol[i]->theta_range / theta;
			Mol[i]->norm = norm;
			Mol[i]->phibulk = Mol[i]->chainlength * norm;

			A += Mol[i]->phibulk * Mol[i]->Charge();
			B += Mol[i]->phibulk;

			Mol[i]->n = norm * Mol[i]->GN;
			Mol[i]->theta = Mol[i]->n * Mol[i]->chainlength;
			Zero(phit, M);
			k = 0;
			while (k < length)
			{
				Real *phi = Mol[i]->phi + k * M;
				Norm(phi, norm, M);
				if (debug)
				{
					Real sum = Lat[0]->ComputeTheta(phi);
					cout << "Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl;
				}
				k++;
			}
		}
	}
	if (charged && neutralizer > -1)
	{
		//Real phib=0;
		//for (int i=0; i<n_mol; i++) {
			//if (i!=neutralizer && i!=solvent) phib+=Mol[i]->phibulk*Mol[i]->Charge();
		//}

		//Mol[neutralizer]->phibulk = -A/Mol[neutralizer]->Charge();
		if (Mol[neutralizer]->Charge()==Mol[solvent]->Charge()) {
			cout << "WARNING: solvent charge equals neutralizer charge; outcome problematic...." << endl;
		} else Mol[neutralizer]->phibulk= ((B-1.0)*Mol[solvent]->Charge() -A)/(Mol[neutralizer]->Charge()-Mol[solvent]->Charge());

		if (Mol[neutralizer]->phibulk<0) {
			cout << "WARNING: neutralizer has negative phibulk. Consider changing neutralizer...: outcome problematic...." << endl;
cout <<"A is " << A << endl;
for (int j=0; j<n_mol; j++) {
	cout << " mol : " << Mol[j]->name << " phibulk " << Mol[j]->phibulk << endl;

}


		}
		B += Mol[neutralizer]->phibulk;
		Real norm = Mol[neutralizer]->phibulk / Mol[neutralizer]->chainlength;
		Mol[neutralizer]->n = norm * Mol[neutralizer]->GN;
		//cout <<"n neutralizer: " << Mol[neutralizer]->n << endl;
		//cout <<"phibulk neutr: " << Mol[neutralizer]->phibulk << endl;
		Mol[neutralizer]->theta = Mol[neutralizer]->n * Mol[neutralizer]->chainlength;
		Mol[neutralizer]->norm = norm;
	}

	Mol[solvent]->phibulk = 1.0 - B;
	//cout <<"phibulk solv: " << Mol[solvent]->phibulk << endl;
	if (Mol[solvent]->phibulk < 0)
	{
		cout << "WARNING: solvent has negative phibulk. outcome problematic " << endl;
		throw - 4;
	}
	Real norm = Mol[solvent]->phibulk / Mol[solvent]->chainlength;
	Mol[solvent]->n = norm * Mol[solvent]->GN;
	Mol[solvent]->theta = Mol[solvent]->n * Mol[solvent]->chainlength;
	Mol[solvent]->norm = norm;

	int k = 0;
	length = Mol[solvent]->MolMonList.size();
	while (k < length)
	{
		Real *phi = Mol[solvent]->phi + k * M;
		if (norm > 0)
			Norm(phi, norm, M);
		if (debug)
		{
			Real sum;
			Sum(sum, phi, M);
			cout << "Sumphi in mol " << solvent << "for mon " << k << ":" << sum << endl;
		}
		k++;
	}
	if (charged && neutralizer > -1)
	{
		k = 0;
		length = Mol[neutralizer]->MolMonList.size();
		while (k < length)
		{
			Real *phi = Mol[neutralizer]->phi + k * M;
			if (Mol[neutralizer]->norm > 0)
				Norm(phi, Mol[neutralizer]->norm, M);
			if (debug)
			{
				Real sum;
				Sum(sum, phi, M);
				cout << "Sumphi in mol " << neutralizer << "for mon " << k << ":" << sum << endl;
			}
			k++;
		}
	}

	for (int i = 0; i < n_mol; i++)
	{

		int length = Mol[i]->MolMonList.size();
		int k = 0;
		while (k < length)
		{
			Real *phi_mon = Seg[Mol[i]->MolMonList[k]]->phi;
			Real *mol_phitot = Mol[i]->phitot;
			Real *phi_molmon = Mol[i]->phi + k * M;
			Add(phi_mon, phi_molmon, M);
			//if (!(Seg[Mol[i]->MolMonList[k]]->freedom == "tagged"))
			Add(phitot, phi_molmon, M);
			Add(mol_phitot, phi_molmon, M);
			Seg[Mol[i]->MolMonList[k]]->phibulk += Mol[i]->fraction(Mol[i]->MolMonList[k]) * Mol[i]->phibulk;
			k++;
		}
		length = SysTagList.size();
		k = 0;
		while (k < length)
		{
			Cp(Seg[SysTagList[k]]->phi, Seg[SysTagList[k]]->MASK, M);
			k++;
		}
		length = SysClampList.size();
		k = 0;
		while (k < length)
		{
			Cp(Seg[SysClampList[k]]->phi, Seg[SysClampList[k]]->MASK, M);
			k++;
		}
	}
	int n_seg = In[0]->MonList.size();
	for (int i = 0; i < n_seg; i++)
		Seg[i]->SetPhiSide();
	return success;
}

bool System::CheckResults(bool e_info_)
{
	if (debug)
		cout << "CheckResults for system " << endl;

	bool e_info = e_info_;
	bool success = true;

	FreeEnergy = GetFreeEnergy();
	GrandPotential = GetGrandPotential();
	CreateMu();

	if (e_info)
		cout << endl;
	if (e_info)
	{
		cout << "free energy                 = " << FreeEnergy << endl;
		cout << "grand potential             = " << GrandPotential << endl;
	}
	Real n_times_mu = 0;
	for (int i = 0; i < n_mol; i++)
	{
		Real Mu = Mol[i]->Mu;
		Real n = Mol[i]->n;
		if (Mol[i]->IsClamped())
			n = Mol[i]->n_box;
		n_times_mu += n * Mu;
	}
	if (e_info)
	{
		cout << "free energy     (GP + n*mu) = " << GrandPotential + n_times_mu << endl;
		cout << "Grand potential (F - n*mu)  = " << FreeEnergy - n_times_mu << endl;
		//	}

		//for (int i=0; i<n_mol; i++) { //NEED FIX . densities are not yet computed correctly that is why it is turned off.....!!!
		//if (Mol[i]->MolAlList.size()>0) {
		//	Mol[i]->compute_phi_alias=true;
		//Mol[i]->ComputePhi();
		//}
		//ComputePhis();
		//}
		cout << endl;
		int M = Lat[0]->M;
		for (int i = 0; i < n_mol; i++)
		{
			int n_molmon = Mol[i]->MolMonList.size();
			Real theta_tot = Mol[i]->n * Mol[i]->chainlength;
			for (int j = 0; j < n_molmon; j++)
			{
				Real FRACTION = Mol[i]->fraction(Mol[i]->MolMonList[j]);
				if (Seg[Mol[i]->MolMonList[j]]->freedom != "clamp")
				{
					Real THETA = Lat[0]->WeightedSum(Mol[i]->phi + j * M);
					cout << "MOL " << Mol[i]->name << " Fraction " << Seg[Mol[i]->MolMonList[j]]->name << ": " << FRACTION << "=?=" << THETA / theta_tot << " or " << THETA << " of " << theta_tot << endl;
				}
			}
		}
		cout << endl;
	}
	return success;
}

Real System::GetFreeEnergy(void)
{ //eqn 2.91 of thesis of J.v.Male;
	if (debug)
		cout << "GetFreeEnergy for system " << endl;
	int M = Lat[0]->M;
	Real FreeEnergy = 0;
	Real *F = FreeEnergyDensity;
	Real constant = 0;
	int n_seg = In[0]->MonList.size();
	int n_mol = In[0]->MolList.size();
	int n_states = In[0]->StateList.size();
	//for (int i=0; i<n_mol; i++) Lat[0]->remove_bounds(Mol[i]->phitot);
	int n_mon = In[0]->MonList.size();
	for (int i = 0; i < n_mon; i++)
	{
		if (Seg[i]->ns < 2)
			Lat[0]->remove_bounds(Seg[i]->phi_side);
		else
		{
			for (int j = 0; j < Seg[i]->ns; j++)
			{
				Lat[0]->remove_bounds(Seg[i]->phi_side + j * M);
			}
		}
	}

	Zero(F, M);

	for (int i = 0; i < n_mol; i++)
	{
		if (Mol[i]->freedom == "clamped")
		{
			int n_box = Mol[i]->n_box;
			#ifdef CUDA
			Real* hgn = NULL;
			TransferDataToHost(hgn, Mol[i]->gn, n_box);
			for (int p = 0; p < n_box; p++)
			{
				FreeEnergy -= log(hgn[p]);
			}
			#else
			for (int p=0; p<n_box; p++)
				FreeEnergy -= log(Mol[i]->gn[p]);
			#endif
		}
		else
		{
			Real n = Mol[i]->n;
			Real GN = Mol[i]->GN;
			int N = Mol[i]->chainlength;
			Real *phi = Mol[i]->phitot; //contains also the tagged segment
			if (Mol[i]->IsTagged())
				N--; //assuming there is just one tagged segment per molecule
			constant = log(N * n / GN) / N;
			Cp(TEMP, phi, M);
			Norm(TEMP, constant, M);
			Add(F, TEMP, M);
		}
	}
	Real *phi;
	Real *phi_side;
	Real *g;
	Real chi;
	int n_sysmon = SysMonList.size();
	for (int j = 0; j < n_sysmon; j++)
	{
		phi = Seg[SysMonList[j]]->phi;
		g = Seg[SysMonList[j]]->G1;
		MinLog(TEMP, g, M);
		Times(TEMP, phi, TEMP, M);
		Norm(TEMP, -1, M);
		Add(F, TEMP, M);
	}

	for (int j = 0; j < n_seg; j++)
	{
		if (Seg[j]->ns < 2)
		{
			phi = Seg[j]->phi;
			for (int k = 0; k < n_seg; k++)
			{
				if (Seg[k]->ns < 2)
				{
					if (Seg[k]->freedom == "frozen" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "tagged")
						chi = Seg[j]->chi[k];
					else
						chi = Seg[j]->chi[k] / 2; //double counted.
					phi_side = Seg[k]->phi_side;
					if (!(Seg[j]->freedom == "frozen" || Seg[j]->freedom == "clamp" || Seg[j]->freedom == "tagged" || chi == 0))
					{
						Times(TEMP, phi, phi_side, M);
						Norm(TEMP, chi, M);
						Add(F, TEMP, M);
					}
				}
			}
			for (int i = 0; i < n_states; i++)
			{
				chi = Seg[j]->chi[n_seg + i]; //not double counted;
				phi_side = Seg[Sta[i]->mon_nr]->phi_side + Sta[i]->state_nr * M;
				if (!(Seg[j]->freedom == "frozen" || Seg[j]->freedom == "clamp" || Seg[j]->freedom == "tagged" || chi == 0))
				{
					Times(TEMP, phi, phi_side, M);
					Norm(TEMP, chi, M);
					Add(F, TEMP, M);
				}
			}
		}
	}

	for (int j = 0; j < n_states; j++)
	{
		phi = Seg[Sta[j]->mon_nr]->phi_state + Sta[j]->state_nr * M;
		for (int k = 0; k < n_seg; k++)
			if (Seg[k]->ns < 2)
			{
				chi = Sta[j]->chi[k];
				if (!(Seg[k]->freedom == "frozen" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "tagged"))
					chi = chi / 2;
				phi_side = Seg[k]->phi_side;
				if (chi != 0)
				{
					Times(TEMP, phi, phi_side, M);
					Norm(TEMP, chi, M);
					Add(F, TEMP, M);
				}
			}
		for (int l = 0; l < n_states; l++)
		{
			chi = Sta[j]->chi[n_seg + l] / 2;
			phi_side = Seg[Sta[j]->mon_nr]->phi_side + Sta[j]->state_nr * M;
			if (chi != 0)
			{
				Times(TEMP, phi, phi_side, M);
				Norm(TEMP, chi, M);
				Add(F, TEMP, M);
			}
		}
		//Cp(TEMP,phi,M); Norm(TEMP,log(Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr]),M); Add(F,TEMP,M);
	}

	for (int i = 0; i < n_mol; i++)
	{
		constant = 0;
		int n_molmon = Mol[i]->MolMonList.size();
		for (int j = 0; j < n_molmon; j++)
			for (int k = 0; k < n_molmon; k++)
			{
				Real fA = Mol[i]->fraction(Mol[i]->MolMonList[j]);
				Real fB = Mol[i]->fraction(Mol[i]->MolMonList[k]);
				if (Mol[i]->IsTagged())
				{
					int N = Mol[i]->chainlength;
					if (N > 1)
					{
						fA = fA * N / (N - 1);
						fB = fB * N / (N - 1);
					}
					else
					{
						fA = 0;
						fB = 0;
					}
				}
				Real chi = CHI[Mol[i]->MolMonList[j] * n_mon + Mol[i]->MolMonList[k]] / 2;
				constant -= fA * fB * chi;
			}
		Real *phi = Mol[i]->phitot;
		Cp(TEMP, phi, M);
		Norm(TEMP, constant, M);
		Add(F, TEMP, M);
	}
	//Lat[0]->remove_bounds(F);
	Times(F, F, KSAM, M); //clean up contributions in frozen and tagged sites.

Zero(TEMP,M);
	if (charged) {
		//Times(TEMP,EE,eps,M);
cout <<"Sum EE*eps = " << Lat[0]->WeightedSum(TEMP) << endl;
		//Norm(TEMP,1.0,M);
		AddTimes(TEMP,q,psi,M);
		Norm(TEMP,0.5,M);

cout <<"El to F = " << Lat[0]->WeightedSum(TEMP)  << endl;
		Add(F,TEMP,M);
	}
	return FreeEnergy + Lat[0]->WeightedSum(F);
}

Real System::GetSpontaneousCurvature()
{
	Real *GP = GrandPotentialDensity;
	//cout <<"Get spontaneous Curvature not yet inplemented in system " << endl;
	return Lat[0]->Moment(GP, 1);
};

Real System::GetGrandPotential(void)
{ //Eqn 293
	if (debug)
		cout << "GetGrandPotential for system " << endl;
	int M = Lat[0]->M;
	Real *GP = GrandPotentialDensity;
	int n_mol = In[0]->MolList.size();
	//int n_mon=In[0]->MonList.size();
	Zero(GP, M);

	for (int i = 0; i < n_mol; i++)
	{
		Real *phi = Mol[i]->phitot;
		Real phibulk = Mol[i]->phibulk;
		int N = Mol[i]->chainlength;
		if (Mol[i]->IsTagged())
		{
			N--;
			phibulk = 0;
		} //One segment of the tagged molecule is tagged and then removed from GP through KSAM
		if (Mol[i]->IsClamped())
		{
			N = N - 2;
			phibulk = 0;
		}
		Cp(TEMP, phi, M);
		YisAplusC(TEMP, TEMP, -phibulk, M);
		Norm(TEMP, 1.0 / N, M); //GP has wrong sign. will be corrected at end of this routine;
		Add(GP, TEMP, M);
	}

	Add(GP,alpha,M);
	Real phibulkA;
	Real phibulkB;
	Real chi;
	Real *phi;
	Real *phi_side;
	int n_seg = In[0]->MonList.size();
	int n_states = In[0]->StateList.size();

	for (int j = 0; j < n_seg; j++)
		if (!(Seg[j]->freedom == "tagged" || Seg[j]->freedom == "clamp" || Seg[j]->freedom == "frozen"))
		{
			if (Seg[j]->ns < 2)
			{
				phi = Seg[j]->phi;
				phibulkA = Seg[j]->phibulk;
				for (int k = 0; k < n_seg; k++)
					if (!(Seg[k]->freedom == "tagged" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "frozen"))
					{
						if (Seg[k]->ns < 2)
						{
							chi = Seg[j]->chi[k] / 2;
							phi_side = Seg[k]->phi_side;
							phibulkB = Seg[k]->phibulk;
							if (chi != 0)
							{
								Times(TEMP, phi, phi_side, M);
								YisAplusC(TEMP, TEMP, -phibulkA * phibulkB, M);
								Norm(TEMP, chi, M);
								Add(GP, TEMP, M);
							}
						}
						else
						{
							for (int k = 0; k < n_states; k++)
							{
								chi = Seg[j]->chi[n_seg + k] / 2; //not double counte;
								phi_side = Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr;
								phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
								if (chi != 0)
								{
									Times(TEMP, phi, phi_side, M);
									YisAplusC(TEMP, TEMP, -phibulkA * phibulkB, M);
									Norm(TEMP, chi, M);
									Add(GP, TEMP, M);
								}
							}
						}
					}
			}
		}

	for (int j = 0; j < n_states; j++)
	{
		phi = Seg[Sta[j]->mon_nr]->phi_state + Sta[j]->state_nr * M;
		phibulkA = Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
		for (int k = 0; k < n_seg; k++)
			if (!(Seg[k]->freedom == "frozen" || Seg[k]->freedom == "clamp" || Seg[k]->freedom == "tagged"))
			{
				if (Seg[k]->ns < 2)
				{
					chi = Sta[j]->chi[k] / 2;
					phibulkB = Seg[k]->phibulk;
					phi_side = Seg[Sta[j]->mon_nr]->phi_side + Sta[k]->state_nr * M;
					if (chi != 0)
					{
						Times(TEMP, phi, phi_side, M);
						YisAplusC(TEMP, TEMP, -phibulkA * phibulkB, M);
						Norm(TEMP, chi, M);
						Add(GP, TEMP, M);
					}
				}
			}
		for (int k = 0; k < n_states; k++)
		{
			chi = Sta[j]->chi[n_seg + k] / 2; //double counted
			phi_side = Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr * M;
			phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
			if (chi != 0)
			{
				Times(TEMP, phi, phi_side, M);
				YisAplusC(TEMP, TEMP, -phibulkA * phibulkB, M);
				Norm(TEMP, chi, M);
				Add(GP, TEMP, M);
			}
		}
	}

	/*
	int n_sysmon=SysMonList.size();
	for (int j=0; j<n_sysmon; j++)for (int k=0; k<n_sysmon; k++){
			if (!(Seg[SysMonList[j]]->freedom=="tagged" || Seg[SysMonList[j]]->freedom=="clamp"  )){
			phibulkA=Seg[SysMonList[j]]->phibulk;
			phibulkB=Seg[SysMonList[k]]->phibulk;
			chi = CHI[SysMonList[j]*n_mon+SysMonList[k]]/2;
			phi=Seg[SysMonList[j]]->phi;
			phi_side=Seg[SysMonList[k]]->phi_side;
			Times(TEMP,phi,phi_side,M); YisAplusC(TEMP,TEMP,-phibulkA*phibulkB,M); Norm(TEMP,chi,M); Add(GP,TEMP,M);
		}
	}*/
	Norm(GP,-1.0,M); //correct the sign.

	Zero(TEMP,M);
if (charged) {
	Times(TEMP,EE,eps,M);
//cout <<"eps EE/2 " << Lat[0]->WeightedSum(TEMP) << endl;
	Norm(TEMP,-2.0,M);

	AddTimes(TEMP,q,psi,M);
	Norm(TEMP,-1.0/2.0,M);

cout <<"el to G " << Lat[0]->WeightedSum(TEMP) << endl;

	Add(GP,TEMP,M); Times(GP,GP,KSAM,M);
	Times(TEMP,q,KSAM,M); YisAminB(TEMP,q,TEMP,M); Times(TEMP,TEMP,psi,M); Norm(TEMP,0.5,M);
	Add(GP,TEMP,M);
}

	if (!charged) Times(GP,GP,KSAM,M); //necessary to make sure that there are no contribution from solid, tagged or clamped sites in GP.

	return  Lat[0]->WeightedSum(GP);

}

bool System::CreateMu()
{
	if (debug)
		cout << "CreateMu for system " << endl;
	//int M=Lat[0]->M;
	bool success = true;
	Real constant;
	Real n;
	Real GN;
	int n_mol = In[0]->MolList.size();
	int n_mon = In[0]->MonList.size();
	for (int i = 0; i < n_mol; i++)
	{
		Real Mu = 0;
		Real NA = Mol[i]->chainlength;
		if (Mol[i]->IsTagged())
			NA = NA - 1;
		if (Mol[i]->IsClamped())
			NA = NA - 2;
		if (Mol[i]->IsClamped())
		{
			n = Mol[i]->n_box;
			int n_box = Mol[i]->n_box;
			GN=0;
			#ifdef CUDA
			Real* hgn = NULL;
			TransferDataToHost(hgn, Mol[i]->gn, n_box);
			for (int p=0; p<n_box; p++) {
				GN += log(hgn[p]);
			}
			#else
			for (int p=0; p<n_box; p++)
				GN += log(Mol[i]->gn[p]);
			#endif
			Mu = -GN / n +1;
		}
		else
		{
			n = Mol[i]->n;
			GN = Mol[i]->GN;
			Mu = log(NA * n / GN) + 1;
		}
		constant = 0;
		for (int k = 0; k < n_mol; k++)
		{
			Real NB = Mol[k]->chainlength;
			if (Mol[k]->IsTagged())
				NB = NB - 1;
			if (Mol[k]->IsClamped())
				NB = NB - 2;
			Real phibulkB = Mol[k]->phibulk;
			constant += phibulkB / NB;
		}
		Mu = Mu - NA * constant;

		//Real theta;
		Real phibulkA;
		Real phibulkB;
		Real FA;
		Real FB;
		Real chi;
		//Real *phi;
		//Real *alpha;
		//int n_states;
		//int n_seg = Mol[i]->MolMonList.size();
		int statelistlength = In[0]->StateList.size();

		for (int j = 0; j < n_mon; j++)
		{
			if (Seg[j]->ns < 2)
			{
				phibulkA = Seg[j]->phibulk;
				FA = Mol[i]->fraction(j);
				if (Mol[i]->IsTagged())
					FA = (NA + 1) / (NA); //works only in case of homopolymers?
				if (Mol[i]->IsClamped())
					FA = (NA + 2) / (NA); //works only when homopolymers are clamped....needs probably a fix.
				for (int k = 0; k < n_mon; k++)
				{
					if (Seg[k]->ns < 2)
					{
						phibulkB = Seg[k]->phibulk;
						FB = Mol[i]->fraction(k);
						if (Mol[i]->IsTagged())
							FB *= (NA + 1) / (NA);
						if (Mol[i]->IsClamped())
							FB *= (NA + 2) / (NA);
						chi = Seg[j]->chi[k] / 2;
						Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
					}
					else
					{
						for (int l = 0; l < statelistlength; l++)
						{
							phibulkB = Seg[Sta[l]->mon_nr]->state_phibulk[Sta[l]->state_nr];
							FB = Mol[i]->fraction(Sta[l]->mon_nr) * Seg[Sta[l]->mon_nr]->state_alphabulk[Sta[l]->state_nr];
							if (Mol[i]->IsTagged())
								FB *= (NA + 1) / (NA);
							if (Mol[i]->IsClamped())
								FB *= (NA + 2) / (NA);
							chi = Seg[j]->chi[n_mon + l] / 2;
							Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
						}
					}
				}
			}
		}
		for (int j = 0; j < statelistlength; j++)
		{
			phibulkA = Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
			FA = Mol[i]->fraction(Sta[j]->mon_nr) * Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr];
			if (Mol[i]->IsTagged())
				FA *= (NA + 1) / (NA);
			if (Mol[i]->IsClamped())
				FA *= (NA + 2) / (NA);
			for (int k = 0; k < n_mon; k++)
				if (Seg[k]->ns < 2)
				{
					phibulkB = Seg[k]->phibulk;
					FB = Mol[i]->fraction(k);
					if (Mol[i]->IsTagged())
						FB *= (NA + 1) / (NA);
					if (Mol[i]->IsClamped())
						FB *= (NA + 2) / (NA);
					chi = Sta[j]->chi[k] / 2;
					Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
				}
			for (int k = 0; k < statelistlength; k++)
			{
				phibulkB = Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
				FB = Mol[i]->fraction(Sta[j]->mon_nr) * Seg[Sta[j]->mon_nr]->state_alphabulk[Sta[j]->state_nr];
				if (Mol[i]->IsTagged())
					FB *= (NA + 1) / (NA);
				if (Mol[i]->IsClamped())
					FB *= (NA + 2) / (NA);
				chi = Sta[j]->chi[n_mon + k] / 2;
				Mu = Mu - NA * chi * (phibulkA - FA) * (phibulkB - FB);
			}
		}
		/*
		for (int j=0; j<n_seg; j++) {
			phi=Mol[i]->phi+j*M;
			n_states=Seg[Mol[i]->MolMonList[j]]->ns;
			if (n_states >1) {
				for (int k=0; k<n_states; k++) {
					alpha=Seg[Mol[i]->MolMonList[j]]->alpha+k*M;
					Times(TEMP,phi,alpha,M); Times(TEMP,TEMP,KSAM,M);
					theta=Lat[0]->WeightedSum(TEMP);
					Mu+=theta*log(Seg[Mol[i]->MolMonList[j]]->state_alphabulk[k])/n;
				}
			}
		}*/

		Mol[i]->Mu = Mu;
	}
	return success;
}

//TODO //make sure that the tagged positions are not in frozen range.

/* ------------------------------------------------------keep for a while---was in mu'
		for (int j=0; j<n_mon; j++) for (int k=0; k<n_mon; k++) {
			Real chi= CHI[j*n_mon+k]/2;
			Real phibulkA=Seg[j]->phibulk;
			Real phibulkB=Seg[k]->phibulk;
			Real Fa=Mol[i]->fraction(j);
			Real Fb=Mol[i]->fraction(k);
			if (Mol[i]->IsTagged()) {Fa=Fa*(NA+1)/(NA); Fb=Fb*(NA+1)/NA;}
			if (Mol[i]->IsClamped()) {Fa=Fa*(NA+2)/(NA); Fb=Fb*(NA+2)/NA;}
			Mu = Mu-NA*chi*(phibulkA-Fa)*(phibulkB-Fb);
			//Mu=Mu-NA*chi*phibulkA*phibulkB;
		}

		for (int j=0; j<n_seg; j++) {
			phi=Mol[i]->phi+j*M;
			n_states=Seg[j]->ns;
			if (n_states <2) {
				theta=Lat[0]->WeightedSum(phi)/n;
				for (int i=0; i<n_mon; i++) {
					if (Seg[i]->ns <2) {
						chi=Seg[j]->chi[i];
						phibulkB=Seg[i]->phibulk;
						Mu+=theta*chi*phibulkB;
					}
				}
				for (int l=0; l<statelistlength; l++) {
					chi = Seg[j]->chi[n_mon+l];  //not double counting
					phibulkB = Seg[Sta[l]->mon_nr]->state_phibulk[Sta[l]->state_nr];
					Mu+=theta*chi*phibulkB;
				}
			} else {
				for (int k=0;k<n_states; k++) {
					alpha=Seg[j]->alpha+k*M;
					Times(TEMP,phi,alpha,M);
					theta=Lat[0]->WeightedSum(TEMP)/n;
					Mu+=theta*log(Seg[j]->state_alphabulk[k]);
				}
			}
		}

		theta = NA;
		for (int j=0; j<n_mon; j++) {
			n_states=Seg[j]->ns;
			if (n_states<2) {
				phibulkA=Seg[j]->phibulk;
				for (int k=0; k<n_mon; k++) {
					if (Seg[k]->ns<2) {
						chi=Seg[j]->chi[k]/2;
						phibulkB=Seg[k]->phibulk;
						Mu-=theta*phibulkA*chi*phibulkB;
					}
				}
				for (int k=0; k<statelistlength; k++) {
					chi=Seg[j]->chi[n_mon+k];
					phibulkB=Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr];
					Mu-=theta*phibulkA*chi*phibulkB;
				}
			}
		}
		for (int j=0; j<statelistlength; j++) {
			phibulkA=Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
			for (int k=0; k<statelistlength; k++) {
				chi=Sta[j]->chi[n_mon+k]/2;
				phibulkB=Seg[Sta[j]->mon_nr]->state_phibulk[Sta[j]->state_nr];
				Mu-=theta*phibulkA*chi*phibulkB;
			}
		}*/

/*         ---------------------------in free energy -----------------------

	for (int j=0; j<n_mon; j++) for (int k=0; k<n_mon; k++) {
		Real chi;
		if (Seg[k]->freedom=="frozen") chi=CHI[j*n_mon+k]; else  chi = CHI[j*n_mon+k]/2;
		phi_side = Seg[k]->phi_side;
		phi = Seg[j]->phi;
		if (Seg[j]->freedom !="frozen") {Times(TEMP,phi,phi_side,M); Norm(TEMP,chi,M); Add(F,TEMP,M);}
	}

	for (int i=0; i<n_mol; i++){
		constant=0;
		int n_molmon=Mol[i]->MolMonList.size();
		for (int j=0; j<n_molmon; j++) for (int k=0; k<n_molmon; k++) {
			Real fA=Mol[i]->fraction(Mol[i]->MolMonList[j]);
			Real fB=Mol[i]->fraction(Mol[i]->MolMonList[k]);
			if (Mol[i]->IsTagged()) {
				int N=Mol[i]->chainlength;
				if (N>1) {fA=fA*N/(N-1); fB=fB*N/(N-1); } else {fA =0; fB=0;}
			}
			Real chi = CHI[Mol[i]->MolMonList[j]*n_mon+Mol[i]->MolMonList[k]]/2;
			constant -=fA*fB*chi;
		}
		Real* phi=Mol[i]->phitot;
		Cp(TEMP,phi,M); Norm(TEMP,constant,M); Add(F,TEMP,M);
	}*/
/*Real System::GetFreeEnergyOld(void) {
  if (debug)
    cout << "GetFreeEnergy for system " << endl;
  int M = Lat[0]->M;
  Real FreeEnergy = 0;
  Real* F = FreeEnergyDensity;
  Real constant = 0;
  int n_mol = In[0]->MolList.size();
  //for (int i=0; i<n_mol; i++) Lat[0]->remove_bounds(Mol[i]->phitot);
  int n_mon = In[0]->MonList.size();
  for (int i = 0; i < n_mon; i++) {
    Lat[0]->remove_bounds(Seg[i]->phi_side);
  }

  Zero(F, M);
  if (charged) {
    Times(F, EE, eps, M);
    Norm(F, -2.0, M); //EE bevat een deling door 2 en de norm met -2 is two pre-correct de latere deling door 2.
    AddTimes(F, q, psi, M);
    Norm(F, 0.5, M);
  }


  for (int i = 0; i < n_mol; i++) {
    if (Mol[i]->freedom == "clamped") {
      int n_box = Mol[i]->n_box;
      for (int p = 0; p < n_box; p++) {
        FreeEnergy -= log(Mol[i]->gn[p]);
      }
    } else {
      Real n = Mol[i]->n;
      Real GN = Mol[i]->GN;
      int N = Mol[i]->chainlength;

      Real* phi = Mol[i]->phitot; //contains also the tagged segment
   	  if (Mol[i]->IsTagged()) N--;                      //assuming there is just one tagged segment per molecule
	  constant = log(N * n / GN)/N;
      Cp(TEMP, phi, M);
      Norm(TEMP, constant, M);
      Add(F, TEMP, M);
    }
  }


	int n_sysmon=SysMonList.size();
	for (int j=0; j<n_sysmon; j++) {
		Real* phi=Seg[SysMonList[j]]->phi;
		Real* u=Seg[SysMonList[j]]->u;
		Times(TEMP,phi,u,M); Norm(TEMP,-1,M); Add(F,TEMP,M);
	}



  for (int i = 0; i < n_mol; i++) {
    int length = Mol[i]->MolMonList.size();
    for (int j = 0; j < length; j++) {
      Real* phi = Mol[i]->phi + M * j;
      Real* u = Mol[i]->u + M * j;
      Times(TEMP, phi, u, M);
      Norm(TEMP, -1, M);
      Add(F, TEMP, M);
    }
  }

 for (int j = 0; j < n_mon; j++)
    for (int k = 0; k < n_mon; k++) {
      Real chi;
      if (Seg[k]->freedom == "frozen")
        chi = CHI[j * n_mon + k];
      else
        chi = CHI[j * n_mon + k] / 2;
      Real* phi_side = Seg[k]->phi_side;
      Real* phi = Seg[j]->phi;
      if (Seg[j]->freedom != "frozen") {
        Times(TEMP, phi, phi_side, M);
        Norm(TEMP, chi, M);
        Add(F, TEMP, M);
      }
    }


  for (int i = 0; i < n_mol; i++) {
    constant = 0;
    int n_molmon = Mol[i]->MolMonList.size();
    for (int j = 0; j < n_molmon; j++)
      for (int k = 0; k < n_molmon; k++) {
        Real fA = Mol[i]->fraction(Mol[i]->MolMonList[j]);
        Real fB = Mol[i]->fraction(Mol[i]->MolMonList[k]);
        if (Mol[i]->IsTagged()) {
          int N = Mol[i]->chainlength;
          if (N > 1) {
            fA = fA * N / (N - 1);
            fB = fB * N / (N - 1);
          } else {
            fA = 0;
            fB = 0;
          }
        }
        Real chi = CHI[Mol[i]->MolMonList[j] * n_mon + Mol[i]->MolMonList[k]] / 2;
        constant -= fA * fB * chi;
		//constant =0;
      }
    Real* phi = Mol[i]->phitot;
    Cp(TEMP, phi, M);
    Norm(TEMP, constant, M);
    Add(F, TEMP, M);
  }

  Lat[0]->remove_bounds(F);
  Times(F, F, KSAM, M);
  return FreeEnergy + Lat[0]->WeightedSum(F);
}


bool System::CreateMuOld() {
  if (debug)
    cout << "CreateMu for system " << endl;
  bool success = true;
  Real constant;
  Real n;
  Real GN;
  int n_mol = In[0]->MolList.size();
  int n_mon = In[0]->MonList.size();
  for (int i = 0; i < n_mol; i++) {
    Real Mu = 0;
    Real NA = Mol[i]->chainlength;
    if (Mol[i]->IsTagged())
      NA = NA - 1;
    if (Mol[i]->IsClamped())
      NA = NA - 2;
    if (Mol[i]->IsClamped()) {
      n = Mol[i]->n_box;
      int n_box = Mol[i]->n_box;
      GN = 0;
      for (int p = 0; p < n_box; p++) {
        GN += log(Mol[i]->gn[p]);
      }
      Mu = -GN / n + 1;
    } else {
      n = Mol[i]->n;
      GN = Mol[i]->GN;
      Mu = log(NA * n / GN) + 1;
    }
    constant = 0;
    for (int k = 0; k < n_mol; k++) {
      Real NB = Mol[k]->chainlength;
      if (Mol[k]->IsTagged())
        NB = NB - 1;
      if (Mol[k]->IsClamped())
        NB = NB - 2;
      Real phibulkB = Mol[k]->phibulk;
      constant += phibulkB / NB;
    }
    Mu = Mu - NA * constant;

    for (int j = 0; j < n_mon; j++)
      for (int k = 0; k < n_mon; k++) {
        Real chi = CHI[j * n_mon + k] / 2;
        Real phibulkA = Seg[j]->phibulk;
        Real phibulkB = Seg[k]->phibulk;
        Real Fa = Mol[i]->fraction(j);
        Real Fb = Mol[i]->fraction(k);
        if (Mol[i]->IsTagged()) {
          Fa = Fa * (NA + 1) / (NA);
          Fb = Fb * (NA + 1) / NA;
        }
        if (Mol[i]->IsClamped()) {
          Fa = Fa * (NA + 2) / (NA);
          Fb = Fb * (NA + 2) / NA;
        }
        Mu = Mu - NA * chi * (phibulkA - Fa) * (phibulkB - Fb);
      }

    Mol[i]->Mu = Mu;
    //cout <<"mol" << i << " = " << Mu << endl;
  }
  return success;
}*/
