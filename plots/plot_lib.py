import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np

def savefig(plt, filename):
    plt.savefig(filename + '.png', dpi=300)
    plt.savefig(filename + '.pdf', dpi=300)

def set_plot_style(fs=(9.15, 8.7)):
    matplotlib.rcParams.update({
                               #'axes.grid': True,
                               #'axes.titlesize': 'medium',
                               'font.family': 'serif',
                               'font.serif': 'Palatino', #'Helvetica Neue',
                               'font.size': 10,
                               #'grid.color': 'w',
                               #'grid.linestyle': '-',
                               #'grid.alpha': 0.5,
                               #'grid.linewidth': 1,
                               'legend.frameon': False,
                               'legend.fancybox': False,
                               'legend.fontsize': 11,
                               'legend.numpoints': 1,
                               'legend.loc': 'best',
                               #'legend.framealpha': 0.7,
                               #'legend.handletextpad': 0.1,
                               #'legend.labelspacing': 0.2,
                               'lines.linewidth': 1.4,
                               'savefig.bbox': 'tight',
                               #'savefig.pad_inches': 0.02,
                               'text.usetex': True,
                               #'text.latex.preamble': r'\usepackage{txfonts}',
                               'xtick.labelsize': 10,
                               'ytick.labelsize': 10,
                               'xtick.direction': 'in',
                               'ytick.direction': 'in',
                               'axes.labelpad': 2,
                               'figure.autolayout': True,
                               })
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(1, 1, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.8)
    ax.minorticks_on()
    ax.tick_params('both', length=6, width=0.8, which='major', pad=5, bottom=True, top=True, left=True, right=True)
    ax.tick_params('both', length=0, width=1, which='minor', pad=10, bottom=True, top=True, left=True, right=True)
    return fig, ax

