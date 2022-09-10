# Selection strings

"""Transmembrane helices definition
Helix TMS1    14-37
Break 37-52
TMS2    54-65  |
Loop C  66-73  |--- TMS2
TMS2    74-83  |
TMS3    86-90 (3-10) 93-113
Helix 4  117-125
TMS4     126-161
TMS5     163-183
TMS6     186-214
Helix 8   216-227
TMS7     243-260
Helix 10  262-272
TMS8      275-308
Helix 12  323-332
TMS9     338-374
TMS10  378-408
TMS11  411-425-kink-426-437
TMS12  441-452-kink-453-466
Break   512-518 
TMS13     525-547
"""
Ca     = "name CA"
TMS1   = "residue 14 to 37"
TMS2   = "residue 61 to 82"
TMS3   = "residue 86 to 113"
Helix4 = "residue 117 to 125"
TMS4   = "residue 126 to 161"
TMS5   = "residue 163 to 183"
TMS6   = "residue 186 to 214"
Helix8 = "residue 216 to 227"
TMS7   = "residue 243 to 260"
Helix10= "residue 262 to 272"
TMS8   = "residue 275 to 308"
Helix12= "residue 323 to 332"
TMS9   = "residue 338 to 374"
TMS10  = "residue 378 to 408"
TMS11  = "residue 411 to 437"
TMS12  = "residue 441 to 466"
TMS13  = "residue 525 to 547"

Phe67  = "residue 67"
Repeat1= "residue 4 to 213"
Repeat2= "residue 243 to 436"

indices_TMS2   =NIS.top.select(TMS2)
indices_Repeat1=NIS.top.select(Repeat1)
indices_Repeat2=NIS.top.select(Repeat2)

interfaceA     =NIS.u.select_atoms("((not name H* and resid 4:213) and around 3.5 (not name H* and resid 243:436) )")
interfaceB     =NIS.u.select_atoms("((not name H* and resid 243:436) and around 3.5 (not name H* and resid 4:213))")
I_sphere       =NIS.u.select_atoms("( not name H* and around 5.0 (not name H* and resname I))")
