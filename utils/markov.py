import numpy as np
def cluster_centroids(traj,cluster,toppdb):
        centroids=[]
        _c=np.concatenate(cluster.dtrajs)
        atom_indices = [a.index for a in toppdb.topology.atoms if a.element.symbol != 'H']
        for i in range(cluster.n_clusters):
            idx=np.where( _c == i )
            _ac=traj[idx[0]]
            for it in idx[1:]:
                _ac=_ac+traj[it]
            distances = np.empty((_ac.n_frames, _ac.n_frames))
            for ix in range(_ac.n_frames):
                distances[ix] = md.rmsd(_ac,_ac, ix, atom_indices=atom_indices)   
            beta = 1
            index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()
            centroids.append(_ac[index])
        return centroids

def prune_pairs(a,b,description=None):
    try:
        b.remove(a)
        if description != None:
            print("removing pair:",description[[i for i in range(0,len(description)) if description[i,0] == a],1])
    except:
        print(a," no in the pairs list")
    return b
#pairs=prune_pairs([4581, 6842],pairs,description=pair_description)
#pairs=prune_pairs([61, 9112],pairs,description=pair_description)
#pairs


def remove_constant(X, threshold=0.001):
    if isinstance(X, np.ndarray):
        X = [X]
    Ds = [np.max(x, axis=0) - np.min(x, axis=0) for x in X]
    D = np.min(np.array(Ds), axis=0)
    Ivar = np.where(D > 0.001)[0]
    Y = [x[:, Ivar] for x in X]
    if len(Y) == 1:
        Y = Y[0]
    return Y

def project_and_cluster(trajfiles, featurizer, sparsify=False, tica=True, lag=100, 
                        scale=True, var_cutoff=0.95, nclusters=100):
    """
    Returns
    -------
    trans_obj, Y, clustering

    """
    X = coor.load(trajfiles, featurizer)
    if sparsify:
        X = remove_constant(X)
    if tica:
        trans_obj = coor.tica(X, lag=lag, var_cutoff=var_cutoff)
    else:
        trans_obj = coor.pca(X, dim=-1, var_cutoff=var_cutoff)
    Y = trans_obj.get_output()
    if scale:
        for y in Y:
            y *= trans_obj.eigenvalues[:trans_obj.dimension()]
    cl_obj = coor.cluster_kmeans(Y, k=ncluster, max_iter=5, fixed_seed=True)
    return trans_obj, Y, cl_obj

def eval_transformer(trans_obj):
    import numpy as np
    # Effective dimension (Really? If we just underestimate the Eigenvalues this value also shrinks...))
    print('Evaluating transformer: ', str(trans_obj.__class__))
    print('effective dimension', np.sum(1.0 - trans_obj.cumvar))
    print('eigenvalues', trans_obj.eigenvalues[:5])
    print('partial eigensum', np.sum(trans_obj.eigenvalues[:10]))
    print('total variance', np.sum(trans_obj.eigenvalues ** 2))
    print()

def search_frame(i1,traj):
    ir=0
    ilen=len(traj[0])
    if i1 > ilen-1:
        while i1 > ilen-1 :
            inext=i1-ilen
            ilen=ilen+len(traj[ir])
            ir=ir+1
        print('tr')
        return traj[ir][inext]
    else:
        return traj[0][i1]
    return traj[ir][inext]
def cluster2trajs(i,trajs,ctrajs):
    idx=np.array(np.where(ctrajs == i ))[0]
    newtraj=search_frame(trajs,idx[0])
    for ir in range(1,len(idx)):
        newtraj=newtraj+search_frame(trajs,idx[ir])
    return newtraj
def rmstraj(irow):
    i,t=irow
    return md.rmsd(t,t, i, atom_indices=atom_indices.ix,parallel=False)
def centroid(t,atom_indices):
    from multiprocessing import Pool
    distances = np.empty((t.n_frames, t.n_frames))
    with Pool(processes=16) as pool:
                distances=np.array(pool.map(rmstraj,[(i,t) for i in range(t.n_frames)]))
    beta = 1
    index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()
    return t[index]
def cluster_centroid(icluster,ctrajs):
    _t=cluster2trajs(icluster,trajs,ctrajs)
    return icluster,centroid(_t,atom_indices) 
    

def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
    """Compute a cross-validated VAMP2 score.

    We randomly split the list of independent trajectories into
    a training and a validation set, compute the VAMP2 score,
    and repeat this process several times.

    Parameters
    ----------
    data : list of numpy.ndarrays
        The input data.
    dim : int
        Number of processes to score; equivalent to the dimension
        after projecting the data with VAMP2.
    lag : int
        Lag time for the VAMP2 scoring.
    number_of_splits : int, optional, default=10
        How often do we repeat the splitting and score calculation.
    validation_fraction : int, optional, default=0.5
        Fraction of trajectories which should go into the validation
        set during a split.
    """
    # we temporarily suppress very short-lived progress bars
    with pyemma.util.contexts.settings(show_progress_bars=False):
        nval = int(len(data) * validation_fraction)
        scores = np.zeros(number_of_splits)
        for n in range(number_of_splits):
            ival = np.random.choice(len(data), size=nval, replace=False)
            vamp = pyemma.coordinates.vamp(
                [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
            scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
    return scores
