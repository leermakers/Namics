mon : C : epsilon : 2
mon : N : epsilon : 10
mon : O : epsilon : 10
mon : Np1 : epsilon : 10
//mon : Np1 : valence : 1
//pKa 12.1
mon : Np2 : epsilon : 10
//mon : Np2 : valence : 0.5
//pKA 6.04
mon : Np3 : epsilon : 10
//mon : Np3 : valence : 1
//pKA 10.87
mon : On : epsilon : 10
//mon : On : valence : -1
//pKa 3.71 ...4.15
mon : S : epsilon : 10

mon : W: epsilon : 80
mol : water : composition : @water(W)
mol : water : freedom : solvent
mol : water : Kw : 50

mon : W : chi_C : 0.8
//estimated from lipid solubilities (cmc vs tail length)
mon : W : chi_N : 0.1 
//estimate by solubility of oligoglycine data

mon : W : chi_O : 0.1
mon : W : chi_On : 0
mon : W : chi_NP1 : 0
mon : W : chi_NP2 : 0
mon : W : chi_NP3 : 0
mon : W : chi_S : 0

mon : C : chi_N : 0.5
mon : C : chi_O : 0.5
mon : C : chi_On : 1
mon : C : chi_NP1 : 1
mon : C : chi_NP2 : 1
mon : C : chi_NP3 : 1
mon : C : chi_S : 0.5

mon : N : chi_O : -1
mon : N : chi_On : 0
mon : N : chi_NP1 : 0
mon : N : chi_NP2 : 0
mon : N : chi_NP3 : 0
mon : N : chi_S : 0

mon : O : chi_On : 0
mon : O : chi_NP1 : 0
mon : O : chi_NP2 : 0
mon : O : chi_NP3 : 0
mon : O : chi_S : 0

mon : On : chi_NP1 : 0
mon : On : chi_NP2 : 0
mon : On : chi_NP3 : 0
mon : On : chi_S : 0

mon : NP1 : chi_NP2 : 0
mon : NP1 : chi_NP3 : 0
mon : NP1 : chi_S : 0

mon : NP2 : chi_NP3 : 0
mon : NP2 : chi_S : 0

mon : NP3 : chi_S : 0

alias : R : value : (N)1(C)1[(C)3(N)1(C)1[(N)1][(Np1)1]](C)1[(O)1]
alias : R_ : value :    (C)1[(C)3(N)1(C)1[(N)1][(Np1)1]](C)1[(O)1]
//R = arginine
alias : H : value : (N)1(C)1[(C)2[(C)2(N)1(C)](Np2)1](C)1[(O)1]
alias : H_ : value :    (C)1[(C)2[(C)2(N)1(C)](Np2)1](C)1[(O)1]
//H = histidine
alias : K : value : (N)1(C)1[(C)4(Np3)1](C)1[(O)1]
alias : K_ : value :    (C)1[(C)4(Np3)1](C)1[(O)1]
//K = lysine
alias : D : value : (N)1(C)1[(C)2[(0)1](On)1](C)1[(O)1]
alias : D_ : value :    (C)1[(C)2[(0)1](On)1](C)1[(O)1]
//D = aspartic Acid
alias : E : value : (N)1(C)1[(C)3[(0)1](On)1](C)1[(O)1]
alias : E_ : value :    (C)1[(C)3[(0)1](On)1](C)1[(O)1]
//E = glutamic Acid
alias : S : value : (N)1(C)1[(C)1(O)1](C)1[(O)1]
alias : S_ : value :    (C)1[(C)1(O)1](C)1[(O)1]
//S = serine
alias : T : value : (N)1(C)1[(C)1[(O)1](C)1](C)1[(O)1]
alias : T_ : value :    (C)1[(C)1[(O)1](C)1](C)1[(O)1]
//T = threonine
alias : N : value : (N)1(C)1[(C)2[(O)1](N)1](C)1[(O)1]
alias : N_ : value :    (C)1[(C)2[(O)1](N)1](C)1[(O)1]
//N = asparagine
alias : Q : value : (N)1(C)1[(C)3[(O)1](N)1](C)1[(O)1]
alias : Q_ : value :    (C)1[(C)3[(O)1](N)1](C)1[(O)1]
//Q = glutamine
alias : C : value : (N)1(C)1[(C)1(S)1](C)1[(O)1]
alias : C_ : value :    (C)1[(C)1(S)1](C)1[(O)1]
//C = cysteine
//C and U (selenocysteine) is the same.
alias : G : value : (N)1(C)1(C)1[(O)1]
alias : G_ : value :    (C)1(C)1[(O)1]
//G = glycine
alias : P : value : (N)1(C)1[(C)2][(C)1](C)1[(O)1]
alias : P_ : value :    (C)1[(C)2][(C)1](C)1[(O)1]
//P = proline
alias : A : value : (N)1(C)1[(C)1](C)1[(O)1]
alias : A_ : value :    (C)1[(C)1](C)1[(O)1]
//A = alanine
alias : I : value : (N)1(C)1[(C)1[(C)1](C)2](C)1[(O)1]
alias : I_ : value :    (C)1[(C)1[(C)1](C)2](C)1[(O)1]
//I = isoleucine
alias : L : value : (N)1(C)1[(C)2[(C)1](C)1](C)1[(O)1]
alias : L_ : value :    (C)1[(C)2[(C)1](C)1](C)1[(O)1]
//L = leucine
alias : M : value : (N)1(C)1[(C)2(S)1(C)1](C)1[(O)1]
alias : M_ : value :    (C)1[(C)2(S)1(C)1](C)1[(O)1]
//M = methionine
alias : F : value : (N)1(C)1[(C)2[(C)2](C)2](C)1[(O)1]
alias : F_ : value :    (C)1[(C)2[(C)2](C)2](C)1[(O)1]
//F = phenylalanine
alias : W : value : (N)1(C)1[(C)2[(C)1(N)1](C)1[(C)2](C)3](C)1[(O)1]
alias : W_ : value :    (C)1[(C)2[(C)1(N)1](C)1[(C)2](C)3](C)1[(O)1]
//W = tryptophan
alias : Y : value : (N)1(C)1[(C)2[(C)2](C)3(O)1](C)1[(O)1]
alias : Y_ : value :    (C)1[(C)2[(C)2](C)3(O)1](C)1[(O)1]
//Y = tyrosine
alias : V : value : (N)1(C)1[(C)1[(C)1](C)1](C)1[(O)1]
alias : V_ : value :    (C)1[(C)1[(C)1](C)1](C)1[(O)1]
//V = valine





