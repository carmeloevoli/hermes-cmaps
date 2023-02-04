import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./hermes.mpltstyle')
import numpy as np


def set_axes(ax):
    ax.set_yscale('log')
    ax.set_ylim([3e-3, 0.1])
    ax.set_xscale('log')
    ax.set_xlim([1e1, 1e6])
    ax.set_xlabel(r'E$_\gamma$ [GeV]')
    ax.set_ylabel(r'E$^{2.5}$ $\Phi$ []')

def plot_flux(ax, filename, color, label, linestyle='-'):
    E, flux = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    E25 = np.power(E, 2.5)
    ax.plot(E, E25 * flux, color=color, linestyle=linestyle, label=label)

def gamma_plot():
    fig = plt.figure(figsize=(12.5, 9.0))
    ax = fig.add_subplot(111)
    set_axes(ax)
    
    plot_flux(ax, 'spectrum-AAfrag-Pi0-HI-128_inner.txt', 'tab:red', 'AAfrag')
    plot_flux(ax, 'spectrum-KA-Pi0-HI-128_inner.txt', 'tab:green', 'KelnerAharonian')

    ax.legend()
    plt.savefig('gamma_plot.pdf')

if __name__== "__main__":
    gamma_plot()
