from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the AutoModel class

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        # Add some restraints from a file:
        # rsr.append(file=’my_rsrs1.rsr’)
        # Residues 20 through 30 should be an alpha helix:
        rsr.add(secondary_structure.Alpha(self.residue_range('1:A', '8:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('44:A', '49:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('55:A', '61:A')))
        #alpha fold predicted
        rsr.add(secondary_structure.Alpha(self.residue_range('30:A', '38:A')))
        #
        r1=RigiBody(self.residue_range('10:A','32:A') or 
                    self.residue_range('61:A','183:A') or 
                    self.residue_range('190:A','483:A') or 
                    self.residue_range('510:A','561:A'))
#        r2=RigiBody(self.residue_range('61:A','183:A'))
#        r3=RigiBody(self.residue_range('190:A','483:A'))
#        r4=RigiBody(self.residue_range('510:A','561:A'))
        rsr.rigid_bodies.append(r1)
#        rsr.rigid_bodies.append(r2)
#        rsr.rigid_bodies.append(r3)
#        rsr.rigid_bodies.append(r4)
        # Two beta-strands:
        # rsr.add(secondary_structure.strand(self.residue_range('1:', '6:'))
        # An anti-parallel sheet composed of the two strands:
        # rsr.add(secondary_structure.sheet(at[’N:1’], at[’O:14’],
        # sheet_h_bonds=-5))
        # Use the following instead for a *parallel* sheet:
        # rsr.add(secondary_structure.sheet(at[’N:1’], at[’O:9’],
        # sheet_h_bonds=5))
        # Restrain the specified CA-CA distance to 10 angstroms (st. dev.=0.1)
        # Use a harmonic potential and X-Y distance group.
        #TOFIX Better edscription of  the binding constraints for the peptide, more?
        #
        #rsr.add(forms.gaussian(group=physical.xy_distance,
        #                       feature=features.distance(at['CB:724:A'],at['CB:125:A']),
        #                       mean=5.0, stdev=0.1))
        #rsr.add(forms.gaussian(group=physical.xy_distance,
        #                       feature=features.distance(at['CB:722:A'],at['CB:126:A']),
        #                       mean=5.4, stdev=0.1))
        rsr.write(file='NIS_modeller.rsr')
