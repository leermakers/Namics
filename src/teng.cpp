// TODO: Code review required. (Look at Checkinput for language reference) 
// Preprocessor Dependencies.
#include "teng.h"
#include "output.h"
#include <string>
// Constructor
Teng::Teng(vector<Input *> In_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<State *> Sta_, vector<Reaction *> Rea_, vector<Molecule *> Mol_, vector<System *> Sys_, vector<Solve_scf *> New_, string name_)
		: name{name_},
			In{In_},
			Lat{Lat_},
			Mol{Mol_},
			Seg{Seg_},
			Sta{Sta_},
			Rea{Rea_},
			Sys{Sys_},
			New{New_}

{
	if (debug)
		cout << "Teng initialized" << endl;
	KEYS.push_back("MCS");
	KEYS.push_back("save_interval");
	KEYS.push_back("move");
	KEYS.push_back("engine");
	//KEYS.push_back("save_filename"); // TODO: Implement options to userdefine filename. 
	//KEYS.push_back("seed");				// seed is not used now and redundant. Can be removed if needed.
}
// Destructor
Teng::~Teng()
{
}

// MonteCarlo Engine to drive the tagged molecules.
// TODO: Needs generalization. 
bool Teng::MonteCarlo()
{
	if (debug)
		cout << "Monte Carlo in Teng" << endl;
	bool success;
	bool solved = false;
	int time = 0;
	Real *copyitvar;
	copyitvar = (Real *)malloc(New[0]->iv * sizeof(Real));
	New[0]->i_info = 100;
	solved = New[0]->Solve(true);
	success = CP(to_teng);
	for (int i = 0; i < New[0]->iv; i++)
	{
		copyitvar[i] = New[0]->xx[i];
	}
	Real F_bm = Sys[0]->FreeEnergy; //Real G_bm = Sys[0]->GrandPotential;
	Real F_am;
	WriteOutput(time);
	Real accepted = 0.0;
	Real rejected = 0.0;
	Real Percentageaccepted = 0.0;
	Real Percentagerejected = 0.0;
	Real acceptance;

	for (time = 1; time <= MCS; time++)
	{
		success = CP(to_bm);
		ChangeMode();
		solved = New[0]->Solve(true);
		F_am = Sys[0]->FreeEnergy;

		acceptance = GetRandom(1.0);

		if ((F_am - F_bm <= 0 || acceptance < exp(-1.0 * (F_am - F_bm))) && solved)
		{
			F_bm = Sys[0]->FreeEnergy; // G_bm=Sys[0]->GrandPotential;
			for (int i = 0; i < New[0]->iv; i++)
			{
				copyitvar[i] = New[0]->xx[i];
			}
			accepted += 1.0;
			cout << "Accepted. Number of MC moves accepted so far: " << accepted << endl;
		}
		else
		{
			success = CP(reset);
			success = CP(to_segment);
			for (int j = 0; j < New[0]->iv; j++)
			{
				New[0]->xx[j] = copyitvar[j];
			}
			solved = New[0]->Solve(true);
			F_bm = Sys[0]->FreeEnergy; //	G_bm=Sys[0]->GrandPotential;
			rejected += 1.0;
			cout << "Rejected. Number of MC moves rejected so far: " << rejected << endl;
		}

		Percentageaccepted = accepted / time;
		Percentagerejected = rejected / time;

		if (time % save_interval == 0)
			WriteOutput(time);
		cout << "MonteCarlo step: " << time << "		Percentage moves accepted: " << Percentageaccepted * 100 << " 	Percentage of rejected moves: " << Percentagerejected * 100 << endl;
	}
	return success;
}
// TODO: Function overloaded to obtain random numbers. Could go to end. 
int Teng::GetRandom(int maxvalue)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distance(1, maxvalue);
	int randomnumber = distance(gen);
	return randomnumber;
}

Real Teng::GetRandom(Real maxvalue)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> distance(0, maxvalue);
	Real randomnumber = distance(gen);
	return randomnumber;
}
// TODO: Changes the constraints as modes. 
// more moves required. 
bool Teng::ChangeMode()
{
	bool success = false;
	bool copy = false;
	while (!success)
	{
		Real Amplitude;
		Real Wavenumber;
		//Real pi = 4.0 * atan(1.0);
		Amplitude = GetRandom(1.0);
		Wavenumber = round(GetRandom(Lat[0]->MZ / 2.0)) * 2; //creates even wavenumbers in uniform space

		//TODO : Select random number of particles and translate them randomly

		if (debug)
			cout << "amplitude: " << Amplitude << " \t wavenumber: " << Wavenumber << endl;
		for (int i = 0; i < n_particles; i++)
		{
			X[i] = X[i];// + round(0.5 - round(GetRandom(1.0)));
			Y[i] = Y[i];// + round(0.5 - round(GetRandom(1.0)));
			Z[i] = Z[i] + round(0.5-round(GetRandom(1.0))); //round(Amplitude * (sin(Wavenumber * pi * X[i] / Lat[0]->MX) * sin(Wavenumber * pi * Y[i] / Lat[0]->MY)));
		}
		success = IsLegal();
		if (success)
			copy = CP(to_segment);
	}
	return copy;
}

// Checks for illegal montecarlo moves. i.e., particle collisions alone. Boundary checky is still not implemented.
// However, boundary checks might be redundant in mode based montecarlo as the modes can never be moving particles out of boundary.
bool Teng::IsLegal()
{
	bool success = true;
	int i, j;
	int xbox = Lat[0]->MX;
	int ybox = Lat[0]->MY;
	int zbox = Lat[0]->MZ;
	// Put molecules back in periodic box or reflect them back based on boundaries.
	for (i = 0; i < n_particles; i++)
	{
		if (X[i] < 1 && Lat[0]->BC[0] == "periodic")
			X[i] += xbox;
		if (Y[i] < 1 && Lat[0]->BC[2] == "periodic")
			Y[i] += ybox;
		if (Z[i] < 1 && Lat[0]->BC[4] == "periodic")
			Z[i] += zbox;
		if (X[i] > xbox && Lat[0]->BC[0] == "periodic")
			X[i] -= xbox;
		if (Y[i] > ybox && Lat[0]->BC[2] == "periodic")
			Y[i] -= ybox;
		if (Z[i] > zbox && Lat[0]->BC[4] == "periodic")
			Z[i] -= zbox;
		if (X[i] < 1 && Lat[0]->BC[0] == "mirror")
			X[i] = 1;
		if (Y[i] < 1 && Lat[0]->BC[2] == "mirror")
			Y[i] = 1;
		if (Z[i] < 1 && Lat[0]->BC[4] == "mirror")
			Z[i] = 1;
		if (X[i] > xbox && Lat[0]->BC[0] == "mirror")
			X[i] = xbox;
		if (Y[i] > ybox && Lat[0]->BC[2] == "mirror")
			Y[i] = ybox;
		if (Z[i] > zbox && Lat[0]->BC[4] == "mirror")
			Z[i] = zbox;
	}

	// Checking for particle out of bounds
	for (i = 0; i < n_particles; i++)
	{
		if (X[i] > Lat[0]->MX || X[i] < 1)
		{
			success = false;
			cout << "This particle with particle id: " << i << "wanted to leave the box in x-direction." << endl;
		}
		if (Y[i] > Lat[0]->MY || Y[i] < 1)
		{
			success = false;
			cout << "This particle with particle id: " << i << "wanted to leave the box in y-direction." << endl;
		}
		if (Z[i] > Lat[0]->MZ || Z[i] < 1)
		{
			success = false;
			cout << "This particle with particle id: " << i << "wanted to leave the box in z-direction." << endl;
		}
	}
	// Checking for particle collisions again
	for (i = 0; i < n_particles; i++)
	{
		for (j = 0; j < i; j++)
		{
			if (i != j)
			{
				if (X[i] == X[j] && Y[i] == Y[j] && Z[i] == Z[j])
				{
					success = CP(reset);
					success = false;
					cout << "Particle collided." << endl;
				}
			}
		}
	}
	return success;
}

// Transfer the particle locations from segment to teng, and vice versa.,
// TODO: delta locations should be now transfered insted. 
// Can it just use Mask file (is Mask file updated everytime?)
bool Teng::CP(transfer tofrom)
{
	int JX = Lat[0]->JX;
	int JY = Lat[0]->JY;
	int M = Lat[0]->M;

	bool success = true;
	int i;
	switch (tofrom)
	{
	case to_teng:
		if(MoveType=="tags"){
			for (i = 0; i < n_particles; i++){
				PX = Seg[tag_seg]->H_P[i] / JX;
				PY = (Seg[tag_seg]->H_P[i] - PX * JX) / JY;
				PZ = (Seg[tag_seg]->H_P[i] - PX * JX - PY * JY);
				X.push_back(PX);
				Y.push_back(PY);
				Z.push_back(PZ);
				X_bm.push_back(PX);
				Y_bm.push_back(PY);
				Z_bm.push_back(PZ);
			}
		}
		else{
			success=TrackInterface();
		}
		break;
	case to_segment:
		if(MoveType=="tags"){
			Zero(Seg[tag_seg]->H_MASK, M); //TODO: (Ram) Not sure If I have to remove the mask, all calculations use H_P after this.
			Zero(Seg[tag_seg]->H_P, n_particles);
			for (i = 0; i < n_particles; i++){
				Seg[tag_seg]->H_P[i] = X[i] * JX + Y[i] * JY + Z[i];
				Seg[tag_seg]->H_MASK[X[i] * JX + Y[i] * JY + Z[i]] = 1;
			}
		}
		else {
			Sys[0]->constraintfields=true;
			Zero(Sys[0]->H_beta,M);
			for (i = 0; i < n_particles; i++){
				Sys[0]->H_beta[X[i] * JX + Y[i] * JY + Z[i]] = 1;
			}

		}
		break;
	case to_bm:
		for (i = 0; i < n_particles; i++)
		{
			X_bm[i] = X[i];
			Y_bm[i] = Y[i];
			Z_bm[i] = Z[i];
		}
		break;
	case reset:
		for (i = 0; i < n_particles; i++)
		{
			X[i] = X_bm[i];
			Y[i] = Y_bm[i];
			Z[i] = Z_bm[i];
		}
		break;
	default:
		success = false;
		cout << "error in transfer" << endl;
		break;
	}
	return success;
}

bool Teng::TrackInterface(){
	bool success=true;
	n_particles=0;
	int JX=Lat[0]->JX; int JY=Lat[0]->JY; int MX=Lat[0]->MX; int MY=Lat[0]->MY; int MZ=Lat[0]->MZ;
	Real sum;
	Real radius;
	for(int i=1; i<=MX; i++){
		for(int j=1; j<=MY; j++){
			sum=0;
			// Now tracks at every x and y and makes a node there. It should rather be limited to a certain number of points.
			// Also assumes that there is only one interface at given (x,y). 
			for(int k=1; k<=MZ; k++){
				sum+=Mol[0]->phi[i*JX+j*JY+k]-Mol[0]->phi[i*JX+j*JY+MZ];
				if(k==MZ) {
					n_particles+=1;
					radius=sum/(Mol[0]->phi[i*JX+j*JY+1]-Mol[0]->phi[i*JX+j*JY+MZ]);
					cout << radius << endl;
					PX = i;
					PY = j;
					PZ = round(radius);
					X.push_back(PX);
					Y.push_back(PY);
					Z.push_back(PZ);
					X_bm.push_back(PX);
					Y_bm.push_back(PY);
					Z_bm.push_back(PZ);

				}
			}
		}
	}
	return success;
}




// TODO: Output modules should be rearranged. 
// Push outputs from this and other classes to Output class
void Teng::WriteOutput(int subloop_)
{
	int subloop = subloop_;
	PushOutput();
	WritePdb(subloop);
	New[0]->PushOutput();
	for (int i = 0; i < n_out; i++)
	{
		Out[i]->WriteOutput(subloop);
	}
}

// TODO: Engine assumes any constraint as particles. This
// part writes the evolution of such points 'delta' or 'tags'
// as the engine runs. 
void Teng::WritePdb(int time)
{
	FILE *fp;
	string filename;
	filename = name;
	filename.append(".pdb");
	fp = fopen(filename.c_str(), "a");
	fprintf(fp, "%s %9d\n", "MODEL", time);
	for (int i = 0; i < n_particles; i++)
		fprintf(fp, "%s%7d%s%12d%7d%7d%7d\n", "ATOM", (time * n_particles) + (i + 1), "  H", (i + 1) + (time * n_particles), X[i], Y[i], Z[i]);
	fprintf(fp, "%s\n", "ENDML");
	fclose(fp);
}

// Additional information to pass outputs from Teng class.
void Teng::PushOutput()
{
	int *point;
	for (int i = 0; i < n_out; i++)
	{
		Out[i]->PointerVectorInt.clear();
		Out[i]->PointerVectorReal.clear();
		Out[i]->SizeVectorInt.clear();
		Out[i]->SizeVectorReal.clear();
		Out[i]->strings.clear();
		Out[i]->strings_value.clear();
		Out[i]->bools.clear();
		Out[i]->bools_value.clear();
		Out[i]->Reals.clear();
		Out[i]->Reals_value.clear();
		Out[i]->ints.clear();
		Out[i]->ints_value.clear();
		if (Out[i]->name == "ana" || Out[i]->name == "kal")
			Out[i]->push("t", t);
		if (Out[i]->name == "ana" || Out[i]->name == "kal")
			Out[i]->push("MCS", MCS);
		if (Out[i]->name == "ana" || Out[i]->name == "vec")
		{																															//example for putting an array of Reals of arbitrary length to output
			string s = "vector;0";																			//the keyword 'vector' is used for Reals; the value 0 is the first vector, use 1 for the next etc,
			Out[i]->push("gn", s);																			//"gn" is the name that will appear in output file
			if(MoveType=="tags"){
				Out[i]->PointerVectorReal.push_back(Mol[tag_mol]->gn);			//this is the pointer to the start of the 'vector' that is reported to output.
				Out[i]->SizeVectorReal.push_back(sizeof(Mol[tag_mol]->gn)); //this is the size of the 'vector' that is reported to output
			}
		}
		if (Out[i]->name == "ana" || Out[i]->name == "pos")
		{ //example for putting 3 array's of Integers of arbitrary length to output
			string s = "array;0";
			Out[i]->push("X", s);
			point = X.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(X.size());
			s = "array;1";
			Out[i]->push("Y", s);
			point = Y.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Y.size());
			s = "array;2";
			Out[i]->push("Z", s);
			point = Z.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Z.size());
		}
	}
}

// Standard procedures.
void Teng::PutParameter(string new_param)
{
	KEYS.push_back(new_param);
}

// Gets value. See if this object could be made a const object. Most parameters obtained are probably not changed in runtime. 
string Teng::GetValue(string parameter)
{
	int i = 0;
	int length = PARAMETERS.size();
	while (i < length)
	{
		if (parameter == PARAMETERS[i])
		{
			return VALUES[i];
		}
		i++;
	}
	return "";
}

// Procedure that acquires the given inputs and checks
// Primary call to all engines (montecarlo or langevin dynamics) are also located here.
// teng : engine_name : MC or MD : true/false
// teng : engine_name : MCS : INTEGER
// teng : engine_name : save_interval : INTEGER
bool Teng::CheckInput(int start)
{
	if (debug)
		cout << "CheckInput in Teng" << endl;
	bool success = true;
	success = In[0]->CheckParameters("teng", name, start, KEYS, PARAMETERS, VALUES);
	// Checks what kind of engine is selected
	if (success)
	{
		vector<string> options;
		// TODO: This check is relevant only if we have more than two engines. Should this now be removed?
		if (GetValue("engine").size() > 0)
		{
			vector<string> engines;
			engines.push_back("MC");
			engines.push_back("MD");
			EngineType = "";
			if (!In[0]->Get_string(GetValue("engine"), EngineType, engines, "At present TransientEngine (Teng) module can only perform 'MC' or 'MD.'"))
			{
				success = false;
			};
		}
		else
		{
			success = false;
			cout << "While specifying Teng please specify 'engine'. Here is the format for doing so ---> Teng: name : engine : move_type " << endl;
		}
		if (success && EngineType == "MC")
		{
			if (GetValue("MCS").size() > 0)
				success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of timesteps should be between 1 and 10000");
			if (debug)
				cout << "MCS is " << MCS << endl;
			if (GetValue("save_interval").size() > 0)
				success = In[0]->Get_int(GetValue("save_interval"), save_interval, 1, MCS, "The save interval nr should be between 1 and MCS (specified)");
			if (debug)
				cout << "Save_interval " << save_interval << endl;
		}
		else
		{
			success = false;
			cerr << "Only MonteCarlo engine has been implemented." << endl;
		}
	}
	// Checks which kind of particle should be moved
	if (success)
	{
		vector<string> options;
		if (GetValue("move").size() > 0)
		{
			vector<string> moves;
			moves.push_back("tags");
			moves.push_back("interface");
			MoveType = "";
			if (!In[0]->Get_string(GetValue("move"), MoveType, moves, "At present TransientEngine (Teng) module can only perform moves on 'tags' or 'interface.'"))
			{
				success = false;
			};
		}
		else
		{
			success = false;
			cout << "While specifying Teng please specify 'move'. Here is the format for doing so ---> Teng: name : move : move_type " << endl;
		}
			
		if (MoveType=="tags"){
			if (Sys[0]->SysTagList.size() < 1){
				cout << "Teng needs to have tagged molecules in the system" << endl;
				success = false;
			}
			else {
				tag_seg = Sys[0]->SysTagList[0];
				if (Sys[0]->SysTagList.size() > 1){
					success = false;
					cout << "Currently the Tagging is limited to one molecule per system. " << endl;
				}
			}
			if (success){
				n_particles = Seg[tag_seg]->n_pos;
				tag_mol = -1;
				int length = In[0]->MolList.size();
				for (int i = 0; i < length; i++){
					if (Mol[i]->freedom == "tagged") tag_mol = i;}
			}
		}
		else if (MoveType=="interface"){	
			//TODO: Track the interface at every x,y
			//Choose the number of points to be made as constraint
			//define n_particles and mc_mol. Change tag_mol to mc_mol to generalize rest of program. 


		}
		else{
			success=false; cerr << "Problem in selecting Movetype" << endl;
		}
	}
	// Creates output class and calls Engine (Currently only Montecarlo is implemented.)
	if (success)
	{
		n_out = In[0]->OutputList.size();
		if (n_out == 0)
			cout << "Warning: no output defined!" << endl;
		for (int i = 0; i < n_out; i++)
		{
			Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
			if (!Out[i]->CheckInput(start))
			{
				cout << "input_error in output " << endl;
				success = false;
			}
		}

		MonteCarlo();
	}
	return success;
}
