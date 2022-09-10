def shiftedColorMap(cmap, min_val, max_val, name):
    '''Function to offset the "center" of a colormap. Useful for data with a negative min and positive max and you want the middle of the colormap's dynamic range to be at zero. Adapted from https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    Input
    -----
      cmap : The matplotlib colormap to be altered.
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.'''
    epsilon = 0.001
    start, stop = 0.0, 1.0
    min_val, max_val = min(0.0, min_val), max(0.0, max_val) # Edit #2
    midpoint = 1.0 - max_val/(max_val + abs(min_val))
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
    # shifted index to match the data
    shift_index = np.hstack([np.linspace(0.0, midpoint, 128, endpoint=False), np.linspace(midpoint, 1.0, 129, endpoint=True)])
    for ri, si in zip(reg_index, shift_index):
        if abs(si - midpoint) < epsilon:
            r, g, b, a = cmap(0.5) # 0.5 = original midpoint.
        else:
            r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap

from matplotlib import cbook
from matplotlib.colors import Normalize
from numpy import ma

class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint

class FigManager(object):
    #Camille Scott Context manager
    def __init__(self, fn='', exts=['svg', 'pdf', 'png', 'eps'], 
                 show=False, nrows=1, ncols=1, 
                 figsize=(18,12), tight_layout=True,
                 **fig_kwds):
        if plt.gcf():
            print('leaving context of', repr(plt.gcf()))
        self.fig, self.ax = plt.subplots(nrows=nrows,
                                         ncols=ncols,
                                         figsize=figsize,
                                         tight_layout=tight_layout,
                                         **fig_kwds)

        self.fn = fn
        self.exts = exts
        self.show = show

        assert self.fig == plt.gcf()

    def __enter__(self):

        return self.fig, self.ax

    def __exit__(self, exc_t, exc_v, trace_b):

        if exc_t is not None:
            print ('ERROR', exc_t, exc_v, trace_b)

        if self.fn:
            print ('saving figure', repr(self.fig))
            for ext in self.exts:
                self.fig.savefig('{}.{}'.format(self.fn, ext))

        if self.show:
            assert self.fig == plt.gcf()
            print ('showing figure', repr(self.fig))
            plt.show(self.fig)

        print ('closing figure', repr(self.fig))
        self.fig.delaxes(self.ax)
        plt.close(self.fig)
        del self.ax
        del self.fig
        print('returning context to', repr(plt.gcf()))

def my_plot_feature_histograms(xyzall,
                            feature_labels=None,
                            ax=None,
                            ylog=False,
                            outfile=None,
                            n_bins=50,
                            ignore_dim_warning=False,
                            **kwargs):
    import numpy as _np
    r"""Feature histogram plot

    Parameters
    ----------
    xyzall : np.ndarray(T, d)
        (Concatenated list of) input features; containing time series data to be plotted.
        Array of T data points in d dimensions (features).
    feature_labels : iterable of str or pyemma.Featurizer, optional, default=None
        Labels of histogramed features, defaults to feature index.
    ax : matplotlib.Axes object, optional, default=None.
        The ax to plot to; if ax=None, a new ax (and fig) is created.
    ylog : boolean, default=False
        If True, plot logarithm of histogram values.
    n_bins : int, default=50
        Number of bins the histogram uses.
    outfile : str, default=None
        If not None, saves plot to this file.
    ignore_dim_warning : boolean, default=False
        Enable plotting for more than 50 dimensions (on your own risk).
    **kwargs: kwargs passed to pyplot.fill_between. See the doc of pyplot for options.

    Returns
    -------
    fig : matplotlib.Figure object
        The figure in which the used ax resides.
    ax : matplotlib.Axes object
        The ax in which the historams were plotted.

    """

    if not isinstance(xyzall, _np.ndarray):
        raise ValueError('Input data hast to be a numpy array. Did you concatenate your data?')

    if xyzall.shape[1] > 50 and not ignore_dim_warning:
        raise RuntimeError('This function is only useful for less than 50 dimensions. Turn-off this warning '
                           'at your own risk with ignore_dim_warning=True.')

    if feature_labels is not None:
        if not isinstance(feature_labels, list):
            from pyemma.coordinates.data.featurization.featurizer import MDFeaturizer as _MDFeaturizer
            if isinstance(feature_labels, _MDFeaturizer):
                feature_labels = feature_labels.describe()
            else:
                raise ValueError('feature_labels must be a list of feature labels, '
                                 'a pyemma featurizer object or None.')
        if not xyzall.shape[1] == len(feature_labels):
            raise ValueError('feature_labels must have the same dimension as the input data xyzall.')

    # make nice plots if user does not decide on color and transparency
#    if 'color' not in kwargs.keys():
#        kwargs['color'] = 'b'
    import matplotlib.colors as _colors
    import matplotlib.pyplot as _plt
    if 'color' not in kwargs.keys():
        lcolor=list(_colors.TABLEAU_COLORS.keys())
    else:
        lcolor=kwargs['color']
    lcolor.reverse()
    if 'alpha' not in kwargs.keys():
        kwargs['alpha'] = .25
    # check input
    if ax is None:
        fig, ax = _plt.subplots()
    else:
        fig = ax.get_figure()

    hist_offset = -.2
    for h, coordinate in enumerate(reversed(xyzall.T)):
        hist, edges = _np.histogram(coordinate, bins=n_bins)
        kwargs['color']=lcolor[h]
        if not ylog:
            y = hist / hist.max()
        else:
            y = _np.zeros_like(hist) + _np.NaN
            pos_idx = hist > 0
            y[pos_idx] = _np.log(hist[pos_idx]) / _np.log(hist[pos_idx]).max()
        ax.fill_between(edges[:-1], y + h + hist_offset, y2=h + hist_offset, **kwargs)
        ax.axhline(y=h + hist_offset, xmin=0, xmax=1, color='k', linewidth=.2)
    ax.set_ylim(hist_offset, h + hist_offset + 1)

    # formatting
    if feature_labels is None:
        feature_labels = [str(n) for n in range(xyzall.shape[1])]
        ax.set_ylabel('Feature histograms')

    ax.set_yticks(_np.array(range(len(feature_labels))) + .3)
    ax.set_yticklabels(feature_labels[::-1])
    ax.set_xlabel('Feature values')

    # save
    if outfile is not None:
        fig.savefig(outfile)
    return fig, ax

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
