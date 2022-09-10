def look_for_pairs(u,order_pairs=None):
    import numpy as np
    _pairs=[]
    ifaceA=u.select_atoms("(not name H* and resid   4:213) and around 3.3 (not name H* and resid 243:436)")
    ifaceB=u.select_atoms("(not name H* and resid 243:436) and around 3.3 (not name H* and resid   4:213)")
    for ix in ifaceA.ix:
        sel='not name H* and resid 243:436 and around 3.3 bynum {:d}'.format(ix+1)
        neib=u.select_atoms(sel)
        for inx in neib.ix:
            aa=[ix,inx]
            if aa not in _pairs:
                _pairs.append(aa)
   # Unique elements way 2:
    _pairs=np.vstack({tuple(row) for row in _pairs}).tolist()
    _pairs=sorted(_pairs, key = lambda x: x[0] )
    if order_pairs == None :
        order_pairs=list(range(len(_pairs)))
    _pairs=np.array(_pairs)[order_pairs].tolist()
    print(f"Number of pairs: {len(_pairs)}")
    return _pairs   

# Unique elements way 1: pairs=np.vstack(set(map(tuple,pairs))).tolist()

# Python code to sort the tuples using second element  
# of sublist Function to sort using sorted() 
def Sort3(sub_li): 
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using the first element or the second element when is an interaction of monomer B  
    # sublist lambda has been used 
    return(sorted(sub_li, key = lambda x: x[0] if x[0] < isep else x[1]))  
