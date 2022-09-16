Molecular simulation of the  Na<sup>+</sup>/I<sup>-</sup> symporter (NIS) in its apo and ions bound forms
Molecular dynamics simulations were carried out using the NIS apo- and holo-structures (2 Na<sup>+</sup>/I<sup>-</sup>
) embedded in a pre-equilibrated 1,2-Dioleoyl-sn-glycero-3-phosphocholine
(DOPC) bilayer. We obtained the parameters and scripts for the minimization, equilibration, and
production in GROMACS using CHARMM-GUI14.  All simulations were carriedout with GROMACS, version 2020, 
in conjunction with the CharMM36 force field, using a TIP3P water model, Berger-derived DOPC lipids, 
and ion parameters by Joung and Cheatham. Van der Waals interactions were cut off at 1 nm, and electrostatics 
were treated by PME beyond 1 nm. Temperature and pressure were kept at 310.5 K and 1 bar using the V-Rescale 
thermostat and Parrinello-Rahman barostat, respectively. All bonds were restrainedusing LINCS, and an integration 
time step of 2 fs was used. We energy-minimized (steepest descent) using a double-precision version of GROMACS, 
and six steps (125, 125, 125, 500, 500, and 500 ns) of position restraint with a gradually lifted harmonic force 
constant (Fc) with different values for backbone atoms (100000, 2000, 1000, 500, 200, and 50 kJ/mol/nm<sup>2</sup>), 
side-chain atoms (2000, 2000, 1000, 500, 200, 50, and 0 kJ/mol/nm<sup>2</sup>), residue dihedrals (1000, 200, 200, 
100, and kJ/mol/nm<sup>2</sup>), and lipids (1000, 400, 400, 200, 40, and 0 kJ/mol/nm<sup>2</sup>). Van der Waals parameters
for I<sup>-</sup> were taken from Li et al. Trajectories (only including protein atoms and ions), the python 3.8 Jupyter notebook 
and its auxiliary files used in the calculations reported in "Structural Insigths into the mechanism of the 
sodium/iodide symporter (NIS)" by Ravera <ital>et al.</ital> 2022.
<img src="https://github.com/mabianchet/Images/NISIO.png">
