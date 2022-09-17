Molecular dynamics simulation of the Na<sup>+</sup>/I<sup>-</sup> symporter (NIS)  ions-bound form
was carried out using the cryoEM structure (2 Na<sup>+</sup>/I<sup>-</sup>
) embedded in a pre-equilibrated 1,2-Dioleoyl-sn-glycero-3-phosphocholine
(DOPC) bilayer. We obtained the parameters and scripts for the minimization, equilibration, and
production using the molecular dynamics program <a href=https://www.gromacs.org/>GROMACS</a> with the script generator web-page <a href=https://charmm-gui.org> CHARMM-GUI </a>. 
The simulation was carried out with GROMACS, version 2020, 
in conjunction with the <a href=https://onlinelibrary.wiley.com/doi/10.1002/jcc.23354>CharMM36 </a>force field, using a TIP3P water model, Berger-derived DOPC lipids, 
and ion parameters by <a hrf=https://doi.org/10.1021/jp8001614> Joung and Cheatham 2008 </a>. Van der Waals interactions were cut off at 1 nm, and electrostatics 
were treated by PME beyond 1 nm. Temperature and pressure was kept at 310.5 K and 1 bar using the V-Rescale 
thermostat and Parrinello-Rahman barostat, respectively. All bonds were restrained using LINCS, and an integration 
time step of 2 fs was used. We energy-minimized (steepest descent)the experimental structure using a double-precision version of GROMACS, 
and performed a six steps (125, 125, 125, 500, 500, and 500 ns) of position restraint equilibration dynamics with a gradually lifted harmonic force 
constant (Fc) with different values of Fc for backbone atoms (100000, 2000, 1000, 500, 200, and 50 kJ/mol/nm<sup>2</sup>), 
side-chain atoms (2000, 2000, 1000, 500, 200, 50, and 0 kJ/mol/nm<sup>2</sup>), residue dihedrals (1000, 200, 200, 
100, and kJ/mol/nm<sup>2</sup>), and lipids (1000, 400, 400, 200, 40, and 0 kJ/mol/nm<sup>2</sup>). Van der Waals parameters
for I<sup>-</sup> were taken from <a href=https://pubs.acs.org/doi/pdf/10.1021/ct500918t>Li <it>et al.</it> 2015</a> This repository include the molecular dynamics trajectories
sampled every 100 ps (only including protein atoms and ions), the python 3.8 Jupyter notebook 
and its auxiliary files used in the calculations reported in "Structural Insigths into the mechanism of the 
sodium/iodide symporter (NIS)"  by <a href=doi: https://doi.org/10.1101/2022.04.07.487502 > Ravera <ital>et al.</ital> 2022</a>.
