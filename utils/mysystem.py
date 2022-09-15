class MySystem:
    """ Main class for structure"""
    def __init__(self,filename,ref=None):
        import mdtraj as md
        import MDAnalysis as mda
        from MDAnalysis.topology.guessers import guess_types
        _fil=filename
        _res='name CA and (residue 10 to 32 or residue 61 to 183 or residue 190 to 483 or residue 510 to 561)'   
        self.t=md.load(filename)
        self.top=self.t.topology
        if not ref == None:
            print('superposed to reference')
            _iref=ref.topology.select(_res)
            _itar =self.top.select(_res)
            self.t.superpose(ref,frame=0,atom_indices=_itar,
                             ref_atom_indices=_iref,parallel=True)
            self.t.save('system.gro')
            _fil='system.gro'
        self.filename=filename
        self._fil=_fil
        self.u=mda.Universe(_fil,guess_bonds=True)
        _u = mda.Universe('equilibration.tpr', 'equilibration.gro')
        _charges=_u.select_atoms("protein or resname SXD or resname IOD or resname REO").charges
        self.u.add_TopologyAttr('charges',_charges)
        _resnames=self.u.residues.resnames
        _resnames[-3:]=['NA','NA','IOD']
        self.u.add_TopologyAttr('resnames',_resnames)
        _names=self.u.atoms.names
        _names[-3:]=['NA','NA','I']
        self.u.add_TopologyAttr('names',_names)
        _elements = guess_types(self.u.atoms.names)
        self.u.add_TopologyAttr('elements',_elements)
        _sel="protein or resname NA or resname IOD or resname REO"
        self.protein=self.t.atom_slice(self.top.select("protein"))
        self.protein.save('protein.pdb')
        self.protein.top=self.protein.topology
        self.protandIons=self.u.select_atoms(_sel)
        self.protein.res=self.u.select_atoms("protein").residues
        self.protandIons.write('protandions.pdb')
    def __repr__(self):
        return self.t.__repr__().replace("mdtraj.Trajectory with",self._fil+" with")
