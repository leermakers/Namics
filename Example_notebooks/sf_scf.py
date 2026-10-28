import numpy as np

class lattice():
    def __init__(self,input,name):
        self.name=name; self.input=input; self.gradients=-1; self.geometry =""; self.n_layers=[];
        self.bc_lb=[]; self.bc_ub=[] 
        
    def check_input(self):
        success=True
        In=self.input()
        self.n_layers=In.lats[-1].n_layers
        self.gradients=len(In.lats[-1].n_layers)
        self.bc_lb=In.lats[-1].bc_lb
        self.bc_ub=In.lats[-1].bc_ub    
        self.geometry=In.lats[-1].geometry
        success = len(self.n_layers)==len(self.bc_lb) and len(self.n_layers)==len(self.bc_ub)
        if (not success): print("In lat, dimensions of 'n_layers' not consitent with 'bc_lb' or 'bc_ub'")
        return success
    
class segment():
    def __init__(self,input,name,lats):
        self.name=name; self.input=input; self.lats=lats;
        self.size=1; self.epsilon=80; self.valence=0; self.seg_nr=-1;
        self.freedom="free"; self.pinned_range=""; self.frozen_range="";
        self.phi=np.zeros(lats[-1].n_layers[0])
            
    def check_input(self):
        success=True
        In=self.input()
        i=0
        for seg in In.segs:
            if seg.name==self.name :
                self.freedom=seg.freedom
                self.size=seg.size
                self.epsilon=seg.epsilon
                self.valence=seg.valence
                self.seg_nr=i
                i=i+1
        return success
 
class molecule():
    def __init__(self,input,name,segs,lats):
        self.name=name; self.input=input; self.segs=segs; self.lats=lats;
        self.composition =""; self.mol_nr=-1;
        self.freedom ="free"; self.theta = 0; self.n =0; self.phibulk=0; self.chainlength=0;
        self.phi=np.random.rand(lats[-1].n_layers[0])
            
    def check_input(self):
        success=True
        In=self.input()
        i=0;
        for mol in In.mols :
            if mol.name==self.name :
                self.composition=mol.composition
                self.freedom = mol.freedom
                self.mol_nr =i;
                if mol.freedom=="free":
                    self.phibulk=mol.phibulk
                elif mol.freedom=="restricted":
                    self.theta=mol.theta
                elif mol.freedom=="solvent":
                    self.theta=0; self.phibulk=0;
                elif mol.freedom=="neutralizer": 
                       self.phbulk=0
                else :
                    print("freedom error in mol "+ name)
                    success=False
            i = i+1 
            return success
            
class system():
    def __init__(self,input,name,mols,segs,lats):
        self.name=name; self.mols=mols; self.input=input; self.mols=mols; self.segs=segs; self.lats=lats;
        self.charged=False; self.neutralizer=-1; self.solvent=-1; 
    def check_input(self):
        success=True
        In=self.input()
        self.charged=In.syss[-1].charged
        for mol in self.mols:
            if mol.freedom == "neutralizer":
                self.neutralizer=mol.mol_nr
            if mol.freedom == "solvent":
                self.solvent=mol.mol_nr
            if self.charged and self.neutralizer<0 : 
                print("system is charged but there is no neutralizer.");
                success=False
            if self.solvent<0:
                print("system has no solvent defined")
                success=False
        return success
    
class sfbox(): 
    def __init__(self,input):
        self.input=input; self.lats=[]; self.segs=[]; self.mols=[]; self.syss=[];
        self.charged=False;
    
    def doit(self,mode):
        if mode == "verbose" : 
            info_level=10
        elif mode == "silent" :
            info_level=2
        else: info_level=1;

        if info_level>2: print("sfbox is notified")
        success=True
        In=self.input()
        for lat in In.lats:
            name=lat.name
            if success : 
                self.lats.extend([lattice(self.input,name)])
                success=self.lats[-1].check_input()
        
        for seg in In.segs:
            if seg.valence != 0: 
                self.charged = True;
            name=seg.name
            if success: 
                self.segs.extend([segment(self.input,name,self.lats)])
                success=self.segs[-1].check_input()
        In.syss[-1].charged=self.charged
         
        for mol in In.mols:
            name=mol.name
            if success:
                self.mols.extend([molecule(self.input,name,self.segs,self.lats)])
                success=self.mols[-1].check_input()

        for sys in In.syss:
            name=sys.name
            if success:
                self.syss.extend([system(self.input,name,self.mols,self.segs,self.lats)])
                success=self.syss[-1].check_input()

        if success:
            #self.syss[-1].doit()
            if info_level>2 : print("sfbox done. You can start data mining now")
        else:
            if info_level>0 : print("sfbox failed to start up.")  