#include "mol_preview.h"


mol_preview::mol_preview(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) {
	In=In_; Seg=Seg_; name=name_;  Lat=Lat_;
if (debug) cout <<"Constructor for Mol " + name << endl;
	KEYS.push_back("freedom");
	KEYS.push_back("composition");
	KEYS.push_back("ring");
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
	KEYS.push_back("save_memory");
	KEYS.push_back("restricted_range");
	KEYS.push_back("compute_width_interface");
	KEYS.push_back("Kw");
	KEYS.push_back("Markov");
	KEYS.push_back("k_stiff");
	save_memory=false;
}

mol_preview::~mol_preview() {
}


bool mol_preview::CheckInput(int start_) {
if (debug) cout <<"In Preview Molecule:: CheckInput  for Mol " + name << endl;
start=start_;
	bool success=true;
	if (!In[0]->CheckParameters("mol",name,start,KEYS,PARAMETERS,VALUES)) {
		success=false;
	} else {
		if (GetValue("composition").size()==0) {cout << "For mol '" + name + "' the definition of 'composition' is required" << endl; success = false;
		} else {
			if (!Decomposition(GetValue("composition")))
				{cout << "For mol '" + name + "' the composition is rejected. " << endl; success=false; return success;}
		}
 	}
	if ( Seg[mon_nr[0]]->freedom == "clamp") freedom="clamped";
	Markov=1;
	if (GetValue("Markov").size()>0) Markov=In[0]->Get_int(GetValue("Markov"),1);
	if (Markov<0 || Markov>2) {
			success=false;
			cout <<"Markov input in mol " +name + " is out of bounds. Integer value should be either '1' or '2'. " << endl;
	}

	return success;
}


int mol_preview::GetAlNr(string s){
if (debug) cout <<"GetAlNr for Mol " + name << endl;
	int n_als=MolAlList.size();
	int found=-1;
	int i=0;
	while(i<n_als) {
		if (Al[i]->name ==s) found=i;
		i++;
	}
	return found;
}

int mol_preview::GetMonNr(string s){
if (debug) cout <<"GetMonNr for Mon " + name << endl;
	int n_segments=In[0]->MonList.size();
	int found=-1;
	int i=0;
	while(i<n_segments) {
		if (Seg[i]->name ==s) found=i;
		i++;
	}
	return found;
}

bool mol_preview::ExpandAlias(vector<string> sub, string &s) {
if (debug) cout <<"Molecule:: ExpandAlias" << endl;
	bool success=true;
	vector<int> open;
	vector<int> close;
	int length_al=sub.size();
	int i=0;
	while (i<length_al-1) {
		string sA;
		sA=sub[i+1];
		if (!In[0]->InSet(In[0]->AliasList,sA)) {
			cout <<"In composition of mol '" + name + "' Alias '" + sA + "' was not found"<<endl; success=false;
		} else {
			int Alnr =GetAlNr(sA);
			if (Alnr<0) {
				Al.push_back(new Alias(In,Lat,sA));
				Alnr=Al.size();
				if (!Al[Alnr-1]->CheckInput(start)) {cout <<"Alias '" + sA + "' in composition not recognised " << endl; return false;}
				MolAlList.push_back(Alnr);
			}
			Alnr =GetAlNr(sA);
			int iv = Al[Alnr]->value;
			string al_comp=Al[Alnr]->composition;
			if (iv < 0) {
				string si;
				stringstream sstm;
				sstm << Alnr;
				si = sstm.str();
				string ssub="";
				sub[i+1]=":"+si+":"+al_comp+":"+si+":";

			} else {
				string sss;
				stringstream sstm;
				sstm << iv;
				sss = sstm.str();
				sub[i+1]=sss;
			}
		}
		i+=2;
	}

	string ss;
	for (int i=0; i<length_al; i++) {
		ss=ss.append(sub[i]);
	}
	s=ss;
//cout <<"expand alias results in " << s << endl;
	return success;
}

bool mol_preview::ExpandBrackets(string &s) {
if (debug) cout <<"Molecule:: ExpandBrackets" << endl;
	bool success=true;
	if (s[0] != '(') {cout <<"illegal composition: " << s << endl; return false;}
	vector<int> open;
	vector<int> close;
	bool done=false; //now interpreted the (expanded) composition
	while (!done) { done = true;
		open.clear(); close.clear();
		if (!In[0]->EvenBrackets(s,open,close)) {
			cout << "In composition of mol '" + name + "' the backets are not balanced."<<endl; success=false;
		 }
		int length=open.size();
		int pos_open;
		int pos_close;
		int pos_low=0;
		int i_open=0; {pos_open=open[0]; pos_low=open[0];}
		int i_close=0; pos_close=close[0];
		if (pos_open > pos_close) {cout << "Brackets open in composition not correct" << endl; return false;}
		while (i_open < length-1 && done) {
			i_open++;
			pos_open=open[i_open];
			if (pos_open < pos_close && done) {
				i_close++; if (i_close<length) pos_close=close[i_close];
			} else {
				if (pos_low==open[i_open-1] ) {
					pos_low=pos_open;
					i_close++; if (i_close<length) pos_close=close[i_close];
					if (pos_open > pos_close) {cout << "Brackets open in composition not correct" << endl; return false;}
				} else {
					done=false;
					int x=In[0]->Get_int(s.substr(pos_close+1),1);
					string sA,sB,sC;
					if (s.substr(pos_open-1,1)=="]") {pos_open --;  }
					sA=s.substr(0,pos_low);
					sB=s.substr(pos_low+1,pos_close-pos_low-1);
					sC=s.substr(pos_open,s.size()-pos_open+1);
					s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
				}
			}
		}
		if (pos_low < open[length-1]&& done) {
			done=false;
			pos_close=close[length-1];
			int x=In[0]->Get_int(s.substr(pos_close+1),0);
			string sA,sB,sC;
			sA=s.substr(0,pos_low);
			sB=s.substr(pos_low+1,pos_close-pos_low-1);
			sC="";
			s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
		}
	}
	return success;
}

bool mol_preview::Interpret(string s,int generation){
if (debug) cout <<"Molecule:: Interpret" << endl;
	bool success=true;
	vector<string>sub;
	vector<int>open;
	vector<int>close;
	In[0]->split(s,':',sub);
	int length_sub =sub.size();
	int AlListLength=MolAlList.size();
	int i=0;

	while (i<length_sub) {
		open.clear(); close.clear();
		In[0]->EvenBrackets(sub[i],open,close);
		if (open.size()==0) {
			int a=In[0]->Get_int(sub[i],0);
			if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;
		} else {
			int k=0;
			int length=open.size();

			while (k<length) {
				string segname=sub[i].substr(open[k]+1,close[k]-open[k]-1);
				int mnr=GetMonNr(segname);
				if (mnr <0)  {cerr <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised"  << endl; success = false;
				// The entire chain of successes seems to be broken, as they are not transferred to or handled correctly by functions down the stack.
				// This function is called multiple times with different results, causing the success handling to break by only transferring the last success result.
				// Throwing to prevent segfaults and other undefined behavior. Caught by CheckInput.
				throw "Composition Error";
				} else {

					int length=Gnr.size();
					if (length>0) {//fragments at branchpoint need to be just 1 segment long.
						if (Gnr[length-1]<generation) {
							if (n_mon[length-1]>1) {
								n_mon[length-1]--;
								n_mon.push_back(1);
								mon_nr.push_back(mon_nr[length-1]);
								Gnr.push_back(Gnr[length-1]);
							 	last_b[Gnr[length-1]]++;
							}
						}
					}
					mon_nr.push_back(mnr);
					Gnr.push_back(generation);
					if (first_s[generation] < 0) first_s[generation]=chainlength;
					if (first_b[generation] < 0) first_b[generation]=mon_nr.size()-1;
					last_b[generation]=mon_nr.size()-1;
					for (int i=0; i<AlListLength; i++) {if (Al[i]->active) Al[i]->frag.push_back(1); else Al[i]->frag.push_back(0);}
				}
				int nn = In[0]->Get_int(sub[i].substr(close[k]+1,s.size()-close[k]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity " << endl; success=false;
				throw "Composition error";
				} else {
					n_mon.push_back(nn);
				}
				chainlength +=nn; last_s[generation]=chainlength;
				k++;
			}
		}
		i++;
	}
	return success;
}

bool mol_preview::GenerateTree(string s,int generation,int &pos, vector<int> open,vector<int> close) {
if (debug) cout <<"Molecule:: GenerateTree" << endl;
	bool success=true;
	string ss;
	int i=0;
	int newgeneration=0;
	int new_generation=0;
	int length=open.size();
	int pos_open=0;
	int pos_close=s.length();
	bool openfound,closedfound;
	while  (pos_open<pos_close && success) {
		pos_open =s.length()+1;
		pos_close=s.length();
		openfound=closedfound=false;
		i=0;
		while (i<length && !(openfound && closedfound) ){
			if (close[i]>pos && !closedfound) {closedfound=true; pos_close=close[i]+1; new_generation=i+1;}
			if (open[i]>pos && !openfound) {openfound=true; pos_open=open[i]+1; newgeneration=i+1;}
			i++;
		}

		if (pos_close<pos_open) {
			ss=s.substr(pos,pos_close-pos);
			if (ss.substr(0,1)=="[") {
				pos=pos+1;
				first_s.push_back(-1);
				last_s.push_back(-1);
				first_b.push_back(-1);
				last_b.push_back(-1);
				GenerateTree(s,new_generation,pos,open,close);
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
				success = Interpret(ss,generation);
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
			success = Interpret(ss,generation);
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			GenerateTree(s,newgeneration,pos,open,close);
		}
	}
	return success;
}

bool mol_preview::Decomposition(string s){
if (debug) cout <<"Decomposition for Mol " + name << endl;
	bool success = true;
	bool aliases = true;
	MolType=linear;//default;
	chainlength=0;
	int loopnr=0;
	vector<int> open;
	vector<int> close;
	vector<string>sub;

	sub.clear();
	In[0]->split(s,'@',sub);

	if (s!=sub[0]) {
		bool keyfound=false;
		vector<string>SUB;
		In[0]->split(sub[1],'(',SUB);
		if (SUB[0]=="water") {
			MolType=water; keyfound=true;
			s=s.substr(7,s.length()-8);
			int mnr=GetMonNr(s);
			if (mnr<0) {
				cout <<"Language for MolType 'water' is as follows: @water(mon_name) wherein mon_name is a valid monomer name with freedom free" << endl;
				success=false; return false;
			} else {
				mon_nr.push_back(mnr);
				n_mon.push_back(1);
			}
		}
		if (SUB[0]=="dend") {
			MolType=dendrimer; keyfound=true;
			s=s.substr(6,s.length()-7);
			if (save_memory) {success=false; cout <<"In dendrimer no save_memory implemented yet. Can be done relatively easily....." << endl; return success;}
			//if (aliases) {success=false; cout <<"In dendrimer no aliases are not implemented yet. Can be done relatively easily....." << endl; return success;}
		}
		if (SUB[0]=="comb") {
			MolType=comb; keyfound=true;
			s=s.substr(6,s.length()-7);
		}
		if (!keyfound) {
				success =false ;
				cout << "Keyword specifying Moltype not recognised: select from @dend, @comb.@water Problem terminated "<< endl ;
				cout << "Perhaps you intent to use an 'alias' in the 'composition'. Do not use @ but bracket the alias as #alias_name# instead. "<< endl ;
			return false;
		 }
	}

	while (aliases && success) {
		loopnr++;
		sub.clear();
		In[0]->split(s,'#',sub);
		if (sub.size()%2!=1) {
			if (s.back()!='#') {
				cout << " Alias in composition should be bracketed on both sides by '#'. For example, (A)10#alias_name#(B)3, when e.g., 'alias : alias_name : value : (X)5' is defined, so that you oubtain (A)10(X)5(B)3 " << endl; success=false;
			}
		}
		aliases=(s!=sub[0]);
		if (aliases) if (!ExpandAlias(sub,s)) {cout <<"alias not properly processed" <<endl; success=false; return false;};
		if (loopnr == 20) {
			cout << "Nesting nr 20 reached in aliases for mol " + name + " -composition. It is decided that this is too deep to continue; Possible, you have defined an alias-A inside alias-B which itself refers to alias-A. This is not allowed. Problem terminated. " << endl;
			success=false; return 0;
		}
	}
	sub.clear();
	In[0]->split(s,',',sub);
	int length=sub.size();
	int length_open;
	int i,j,k,a,f,dd;
	string ss;
	switch(MolType) {
		case water:
			break;
		case dendrimer:
			for (i=0; i<length; i++) {
				open.clear(); close.clear();
				In[0]->EvenBrackets(sub[i],open,close);
				length_open=open.size();
				//if (sub[i].substr(0,1)=="(") {
				if (length_open>0) {
					if (!ExpandBrackets(sub[i])) cout <<"brackets '(' ')' not well positioned in "+sub[i] << endl;
          				success=false;
				}
			}
			for (i=0; i<length; i++) {
				ss.append(sub[i]);
				if (i<length-1) ss.append(",");
			}
			s=ss;
		break;
		case comb:
			//success=false;
		break;;
		default:
			if (!ExpandBrackets(s)) success=false;
		break;
	}
	if (!In[0]->EvenSquareBrackets(s,open,close)) {
		// Another success bool that doesn't travel down the stack.
		// The error message drowns in the rest of the output really quickly, throwing for safety instead.
		// Caught by CheckInput
		cout << "Error in composition of mol '" + name + "'; the square brackets are not balanced in: " << s << endl;
		success=false;
		throw "Composition error";
	}
	if (open.size()>0) {
		if (MolType==linear) { //overwrite default.
			MolType=branched;
		} else {
			success=false;
			cout <<" In 'composition' you can not combine special keywords such as 'dend, or comb with branched compositions. This implies that square brackets are not allowed. " <<endl ;
			return success;
		}
	}
	int generation=0;
	int pos=0;

	MolMonList.clear();
	int AlListLength=MolAlList.size();
	for (int i=0;i<AlListLength; i++) Al[i]->active=false;
	vector<string>sub_gen;
	vector<string>sub_dd;
	int length_g,length_dd,mnr,nn,arm=-1,degeneracy=1,arms=0;
	string segname;
	int N=0;
	int chainlength_backbone,chainlength_arm;
	switch(MolType) {
		case water:
			break;
		case linear:
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			success = GenerateTree(s,generation,pos,open,close);
			break;
		case branched:
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			success = GenerateTree(s,generation,pos,open,close);
			break;
		case dendrimer:
			first_a.clear();
			last_a.clear();
			first_b.clear();
			last_b.clear();
			first_s.clear();
			last_s.clear();
			n_arm.clear();
			mon_nr.clear();
			n_mon.clear();
			In[0]->split(s,';',sub_gen);
			n_generations=sub_gen.size();
			sym_dend=true; //default.
			//cout <<"n_generations " << n_generations << endl;
			chainlength=0; N=-1;
			for (i=0; i<n_generations; i++) {
				sub.clear();	arms=0;
				In[0]->split(sub_gen[i],',',sub);
				int sublength=sub.size();
				if ((sublength-1)%2 >0) {success=false; cout << "In composition of dend for generation " << i << " that is, in " + sub_gen[i] + " the number of arguments is not 3; use @dend(?) for help. " << endl; }
				if (sublength>3) MolType=asym_dendrimer;
				if (sublength<2) { success=false;

					cout<<" ------Dendrimer language:------ " << endl;
					cout<<" example 1: consider segment name 'A' for the branch points and spacers (B)2 with functionality 3 and 4 generations: " << endl;
					cout<<" @dend(A,(B)2,3; A,(B)2,3; A,(B)2,3; A,(B)2,3)  " << endl;
					cout<<" hence generations are seperated by ';', branch points are 'monomer_names' without '(', ')'. Each generation can have unique numbers and structures. " << endl;
					cout<<" expample 2: asymmetric dendrimers, with asymmetry in branches of first generation:" << endl;
					cout <<" @dend(A,(B)2,3,(B)4,1; A,(B)2,3,(B)2,2; A,(B)2,3; A,(B)2,3)  " << endl;
					cout <<" example 2: star with 10 arms " << endl;
					cout <<" @dend(A,(B)100,10) " << endl;
					cout <<" Note that in the 'standard' dendrimer the branches are symmetric -all are equal- " << endl;
					cout <<" example 4: 10 arm star end-grafted by one arm has end-segment 'X' which, e.g., might be pinned to a surface" << endl;
					cout <<" @dend(A,(B)100,9,(B)99(X)1,1)" << endl;
					cout<<" ------Dendrimer language:------ " << endl;
				}

				length_g = sub.size();
				sub_dd.clear();
				In[0]->split(sub[0],':',sub_dd);
				length_dd =sub_dd.size();
				if (sub_dd.size()>1) {
					cout <<"In dendrimer alias not yet implemented" << endl; success=false; return success;
				}

				if (length_dd==4) {
					mnr=GetMonNr(sub_dd[2]);
					a=In[0]->Get_int(sub_dd[2],0);
					if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;

				} else mnr=GetMonNr(sub[0]);
				if (mnr <0)  {success=false; cout <<"In composition of mol '" + name + "', segment name '" + sub_dd[0] + "' is not recognised"  << endl; }
				for (a=0; a<AlListLength; a++) {if (Al[a]->active) Al[a]->frag.push_back(1); else Al[a]->frag.push_back(0);}

				n_mon.push_back(1);
				mon_nr.push_back(mnr);
				d_mon.push_back(degeneracy);
				first_a.push_back(-1);
				last_a.push_back(-1);
				k=1; chainlength+=degeneracy; N++;

				while (k<length_g-1) {
					arm++;
					first_s.push_back(N+1);
          				last_s.push_back(-1);
        				first_b.push_back(-1);
         				last_b.push_back(-1);
					if (first_a[first_a.size()-1]==-1) first_a[first_a.size()-1]=arm;
					last_a[last_a.size()-1]=arm;

					f=In[0]->Get_int(sub[k+1],0); //should not contain double dots....
					if (f<1) {
						success=false; cout <<"In dendrimer-composition, in generation "<<i << " an integer number is expected at argument " << k+1 << " problem terminated" << endl;
					}
					n_arm.push_back(f); arms+=f;

					sub_dd.clear();
					In[0]->split(sub[k],':',sub_dd);
					dd=0; length_dd=sub_dd.size();
					while (dd<length_dd) {
						open.clear(); close.clear();
            					In[0]->EvenBrackets(sub_dd[dd],open,close);
						if (open.size()==0) {
							a=In[0]->Get_int(sub_dd[dd],-1);
							if (a==-1) {
								cout <<"No integer found. Possibly you have a segment name in composition that is not surrounded by brackets " << endl; success=false;
							} else {if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;}
						} else {
							j=0; length=open.size();
							while (j<length) {
								segname=sub_dd[dd].substr(open[j]+1,close[j]-open[j]-1);
								mnr=GetMonNr(segname);
								if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; this occurs at generation " <<i << " arm " << k << "."  << endl; success=false;}
								mon_nr.push_back(mnr);
								nn=In[0]->Get_int(sub_dd[dd].substr(close[j]+1,s.size()-close[j]-1),0);
								if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; this occurs at generation " <<i << " arm " <<k <<"."<< endl; success=false;}
								n_mon.push_back(nn); N+=nn;
								d_mon.push_back(degeneracy*f);
								chainlength +=degeneracy*nn*f;
                if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
                last_b[first_b.size()-1]=mon_nr.size()-1;
								last_s[last_s.size()-1]=N;
								for (a=0; a<AlListLength; a++) {if (Al[a]->active) Al[a]->frag.push_back(1); else Al[a]->frag.push_back(0);}
								j++;
							}
						}
						dd++;
					}
					k=k+2;
				}
				degeneracy*=arms;
			}
   //anticipated that asymmetric dendrimers will be of interest in the future
//for (int i=0; i<n_generations; i++) cout<<"generation " << i << " first_a " << first_a[i] << " last_a " << last_a[i] << endl;
//length=n_arm.size();
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_b " << first_b[i] << " last_b " << last_b[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_s " << first_s[i] << " last_s " << last_s[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " n_arm " << n_arm[i] << endl;
//length=n_mon.size();
//for (int i=0; i<length; i++) cout <<"block " << i << " n_mon " << n_mon[i] << " mon_nr " << mon_nr[i] << " d_mon " << d_mon[i] << endl;
//cout <<"Chain length =" << chainlength << endl;
//cout <<" N = " << N << endl;
//cout <<"Composition = " <<  s << endl;

			break;
		case comb:
			first_a.clear();
			last_a.clear();
			first_b.clear();
			last_b.clear();
			first_s.clear();
			last_s.clear();
			n_arm.clear();
			mon_nr.clear();
			n_mon.clear();
			In[0]->split(s,';',sub_gen);
			n_generations=sub_gen.size();
			int j;
			int mnr;
			int nn;

			sub_dd.clear();
			In[0]->split(s,':',sub_dd);
			length_dd=sub_dd.size();
			if (length_dd>1) {
				success = false;
				cout <<"In comb the use of 'aliases' is not allowed (yet). " << endl; return success;
			}
			if (save_memory) {
				success = false;
				cout <<"In comb the use of 'save_memory' is not allowed (yet). " << endl; return success;
			}

			if (n_generations != 3) {
				success=false;
				cout<<" ------Comb language:------ " << endl;
				cout<<" generic example: @comb((A)5; A,(B)3,(A)6,4; (A)7 )" << endl;
				cout<<" Backbnone has a first part of 5 A-segments, then spacers with 6 A units and a trailing of 7 A-segments. " << endl;
				cout<<" Teeth are composed of 3 B-segments " << endl;
				cout<<" Number of repeats (inbetween the two ';') is eqaul to 4. " << endl;
				cout<<" Number of segments in molelecule is 5+(1+3+6)x4+7 = 52 segments" << endl;
				cout<<" Backbone length is 5+(1+6)x4+7 = 40 segments " << endl;
				cout<<" ------Comb language:------ " << endl;
				return success;
			}
			chainlength=0; N=-1;
			i=0;
			//success=Interpret(sub_gen[0],i);

			first_s.push_back(N+1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			open.clear(); close.clear();
			In[0]->EvenBrackets(sub_gen[0],open,close);
			j=0; length=open.size();
			while (j<length) {
				segname=sub_gen[0].substr(open[j]+1,close[j]-open[j]-1);
				mnr=GetMonNr(segname);
				if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; this occurs at generation " << endl; success=false;}
				mon_nr.push_back(mnr); d_mon.push_back(1);
				nn=In[0]->Get_int(sub_gen[0].substr(close[j]+1,s.size()-close[j]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; this occurs at generation " << endl; success=false;}
				n_mon.push_back(nn); N+=nn;
				chainlength +=nn;
				if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
				last_b[first_b.size()-1]=mon_nr.size()-1;
				last_s[last_s.size()-1]=N;
				j++;
			}

			chainlength_backbone=chainlength;
			i=1;
			sub.clear();
			In[0]->split(sub_gen[1],',',sub);
			//cout <<"sub_gen[1] " << sub_gen[1] << "size of sub : " << sub.size() << endl;
			if (sub.size() != 4) {
				cout << "Central part in comb definition should contain four arguments and three ',' to separare them. Use comb(? ) for details."  << endl;
				success=false; return success;
			}

			n_arm.push_back(In[0]->Get_int(sub[3],0));
			first_s.push_back(N+1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			open.clear(); close.clear();
			In[0]->EvenBrackets(sub[1],open,close);
			j=0; length=open.size();
			while (j<length) {
				segname=sub[1].substr(open[j]+1,close[j]-open[j]-1);
				mnr=GetMonNr(segname);
				if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; this occurs at generation " << endl; success=false;}
				mon_nr.push_back(mnr); d_mon.push_back(n_arm[0]);
				nn=In[0]->Get_int(sub[1].substr(close[j]+1,s.size()-close[j]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; this occurs at generation " << endl; success=false;}
				n_mon.push_back(nn); N+=nn;
				chainlength +=nn;
				if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
				last_b[first_b.size()-1]=mon_nr.size()-1;
				last_s[last_s.size()-1]=N;
				j++;
			}

			chainlength_arm=chainlength-chainlength_backbone;
			chainlength=chainlength_backbone;

			if (n_arm[0] <1) {
				success=false; cout <<" Error in composition of mol "+ name + " number of arms is less than unity. Use comb(? ) for details. " << endl; return success;
			}
			n_generations=3;
			for (int a=1; a<=n_arm[0]; a++) {
				mnr=GetMonNr(sub[0]);
				if (mnr <0)  {success=false; cout <<"In composition of mol '" + name + "', segment name '" + sub_dd[0] + "' is not recognised. For the branching point the ( ) are not needed."  << endl; return success;}


				n_mon.push_back(1); d_mon.push_back(1);
				mon_nr.push_back(mnr); N++;
				chainlength +=chainlength_arm+1;
				n_generations++;
				//success=Interpret(sub[2],n_generations-1);

				first_s.push_back(N+1);
				last_s.push_back(-1);
				first_b.push_back(-1);
				last_b.push_back(-1);
				open.clear(); close.clear();
				In[0]->EvenBrackets(sub[2],open,close);
				j=0; length=open.size();
				while (j<length) {
					segname=sub[2].substr(open[j]+1,close[j]-open[j]-1);
					mnr=GetMonNr(segname);
					if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; Use ring(?) for details." << endl; success=false;}
					mon_nr.push_back(mnr); d_mon.push_back(1);
					nn=In[0]->Get_int(sub[2].substr(close[j]+1,s.size()-close[j]-1),0);
					if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; Use ring(? ) for details. "<< endl; success=false;}
					n_mon.push_back(nn); N+=nn;
					chainlength +=nn;
					if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
					last_b[first_b.size()-1]=mon_nr.size()-1;
					last_s[last_s.size()-1]=N;
					j++;
				}

			}
			n_generations++;

			//success=Interpret(sub_gen[2],n_generations-1);
			first_s.push_back(N+1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			open.clear(); close.clear();
			In[0]->EvenBrackets(sub_gen[2],open,close);
			j=0; length=open.size();
			while (j<length) {
				segname=sub_gen[2].substr(open[j]+1,close[j]-open[j]-1);
				mnr=GetMonNr(segname);
				if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; " << endl; success=false;}
				mon_nr.push_back(mnr); d_mon.push_back(1);
				nn=In[0]->Get_int(sub_gen[2].substr(close[j]+1,s.size()-close[j]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity;  " << endl; success=false;}
				n_mon.push_back(nn); N+=nn;
				chainlength +=nn;
				if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
				last_b[first_b.size()-1]=mon_nr.size()-1;
				last_s[last_s.size()-1]=N;
				j++;
			}

			//length=first_b.size();
			//for (int i=0; i<length; i++) cout <<"first_b " <<  first_b[i] << " last_b " << last_b[i] << endl;
			//for (int i=0; i<length; i++) cout <<"first_s " <<  first_s[i] << " last_s " << last_s[i] << endl;
			break;
					default:
			break;
	}

	if (MolType==branched) { //invert numbers;
		int g_length=first_s.size();
		int length = n_mon.size();
		int xxx;
		int ChainLength=last_s[0];

		for (int i=0; i<g_length; i++) {
			first_s[i] = ChainLength-first_s[i]-1;
			last_s[i] = ChainLength-last_s[i]-1;
			xxx=first_s[i]; first_s[i]=last_s[i]+1; last_s[i]=xxx;
			first_b[i]=length-first_b[i]-1;
			last_b[i]=length-last_b[i]-1;
			xxx=first_b[i]; first_b[i]=last_b[i]; last_b[i]=xxx;
		}

		for (int i=0; i<length/2; i++) {
			xxx=Gnr[i]; Gnr[i]=Gnr[length-1-i]; Gnr[length-1-i]=xxx;
			xxx=n_mon[i]; n_mon[i]=n_mon[length-1-i]; n_mon[length-1-i]=xxx;
			xxx=mon_nr[i]; mon_nr[i]=mon_nr[length-1-i]; mon_nr[length-1-i]=xxx;
			for (int j=0; j<AlListLength; j++) {
				xxx=Al[j]->frag[i]; Al[j]->frag[i]=Al[j]->frag[length-1-i]; Al[j]->frag[length-1-i]=xxx;
			}
		}
	}

	success=MakeMonList();
	if (chainlength==1 && MolType!=water) MolType=monomer;
	return success;
}

int mol_preview::GetChainlength(void){
if (debug) cout <<"GetChainlength for Mol " + name << endl;
	return chainlength;
}

bool mol_preview:: MakeMonList(void) {
if (debug) cout <<"Molecule:: MakeMonList" << endl;
	bool success=true;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (!In[0]->InSet(MolMonList,mon_nr[i])) {
			if (Seg[mon_nr[i]]->GetFreedom()=="frozen") {
				success = false;
				cout << "In 'composition of mol " + name + ", a segment was found with freedom 'frozen'. This is not permitted. " << endl;

			}
			MolMonList.push_back(mon_nr[i]);
		}
		i++;
	}
	i=0;
	int pos;
	while (i<length) {
		if (In[0]->InSet(MolMonList,pos,mon_nr[i])) {molmon_nr.push_back(pos);
		//cout << "in frag i " << i << " there is segment nr " <<  mon_nr[i] << " and it is on molmon  pos " << pos << endl;
		} else {cout <<"program error in mol PrepareForCalcualations" << endl; }
		i++;
	}
	return success;
}



string mol_preview::GetValue(string parameter) {
if (debug) cout <<"GetValue " + parameter + " for Mol " + name << endl;
	int length = PARAMETERS.size();
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}
		i++;
	}
	return "";
}

