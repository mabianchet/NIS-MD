from seaborn import xkcd_rgb
from matplotlib import colors
from matplotlib.colors import cnames as mpl_colors

msme_rgb = {
            'rawdenim': '#547098',
            'beryl': '#68bcbe',
            'dijon': '#c4a15d',
            'carbon': '#525055',
            'pomegranate': '#d64d67',
            'oxblood': '#321613',
            'marush': '#706f2c',
            'mar': '#6a8485',
            'tarragon': '#ffb24d',
            'cochineal': '#ff3b44'
            }


__all__ = ['all_colors', 'msme_rgb', 'xkcd_rgb']

all_colors = {}
for d in (msme_rgb, xkcd_rgb, mpl_colors):
    all_colors.update(d)

# Add the single letter colors.
for name, rgb in colors.ColorConverter.colors.items():
    hex_ = colors.rgb2hex(rgb)
    all_colors.update(((name, hex_), ))
