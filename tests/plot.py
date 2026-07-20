import numpy as np
import matplotlib.pyplot as plt
import mrock_centralized_scripts.path_appender as __ap
__ap.append()
from get_data import load_panda, hubbard_params
import continued_fraction_pandas as cf
import plot_settings as ps
import sys
import os


def load_data(name):
    """ Returns the data for the Hubbard test with the name 'name'.
    sc_cdw: SC-CDW phase
    sc: SC phase
    cdw: CDW phase
    afm: AFM phase
    """
    
    __CUBE_DIR__ = os.path.join("hubbard", "cube")
    __SQUARE_DIR__ = os.path.join("hubbard", "square")
    if name == "sc_cdw":
        return load_panda(__SQUARE_DIR__, "test", "resolvents.json.gz", **hubbard_params(0.0, -2.5, 0.))
    elif name == "sc":
        return load_panda(__CUBE_DIR__, "test", "resolvents.json.gz", **hubbard_params(0.0, -2.5, -0.1))
    elif name == "cdw":
        return load_panda(__SQUARE_DIR__, "test", "resolvents.json.gz", **hubbard_params(0.0, -2.5, 0.1))
    elif name == "afm":
        return load_panda(__SQUARE_DIR__, "test", "resolvents.json.gz", **hubbard_params(0.0, 2.5, 0.1))
    else:
        raise ValueError("Hubbard test: Invalid index")

def create_plot(name):
    pd_data = load_data(name)
    resolvents = cf.ContinuedFraction(pd_data)

    fig, ax = plt.subplots()
    ax.set_ylim(-0.05, 1.)
    ax.set_xlabel(r"$\omega / t$")
    ax.set_ylabel(r"$\mathcal{A} (\omega)  / t^{-1}$")

    colors = np.array(( "blue", "orange", "black", "limegreen", "deepskyblue", "magenta" ))
    linestyles = ["-", "-.", "--", "-", "--", ":"]

    w_lin = np.linspace(-0.01, pd_data["continuum_boundaries"][1] + 0.3, 5000, dtype=complex)
    w_lin += 1e-6j

    ax.plot(w_lin.real, resolvents.spectral_density(w_lin, "phase_SC"),                  label="Phase", ls=linestyles[0], c=colors[0])
    ax.plot(w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC"),              label="Higgs", ls=linestyles[1], c=colors[1])
    ax.plot(w_lin.real, resolvents.spectral_density(w_lin, "amplitude_CDW"),             label="CDW"  , ls=linestyles[2], c=colors[2])
    ax.plot(w_lin.real, resolvents.spectral_density(w_lin, "amplitude_AFM"),             label="l.AFM", ls=linestyles[3], c=colors[3])
    ax.plot(w_lin.real, resolvents.spectral_density(w_lin, "amplitude_AFM_transversal"), label="t.AFM", ls=linestyles[4], c=colors[4])

    resolvents.mark_continuum(ax)
    ax.set_title(f"Hubbard: {name}")
    ax.legend()
    fig.tight_layout()
    plt.show()
    
if len(sys.argv) > 1:
    create_plot(sys.argv[1])
else:
    print("Please provide the name of the test you would like to plot.")