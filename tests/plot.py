import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Import the continued-fraction module.
# This module is used to evaluate resolvents and spectral densities from the
# continued-fraction coefficients stored in the loaded data.
import mrock.continued_fraction as cf

# Import data-loading utilities and parameter helper functions.
# This includes DataLoader and continuum_params.
from mrock.get_data import *

# Create a global DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()


def load_data(name):
    """Returns the data for the Hubbard test case with the given name.

    Available names
    ---------------
    sc_cdw
        Superconducting-charge-density-wave phase.

    sc
        Superconducting phase.

    cdw
        Charge-density-wave phase.

    afm
        Antiferromagnetic phase.
    """

    # Define the two lattice/data subdirectories used by the test cases.
    #
    # os.path.join is used so that the path works on different operating
    # systems, e.g. Linux/macOS and Windows.
    CUBE_DIR = os.path.join("cube", "test")
    SQUARE_DIR = os.path.join("square", "test")

    # Each branch below loads one specific Hubbard-model simulation.
    #
    # The data are stored under:
    #
    #     data / hubbard / <subdir> / parameter-folders / resolvents.json.gz
    #
    # where <subdir> is either:
    #
    #     cube/test
    #
    # or:
    #
    #     square/test
    #
    # The parameter-folder path is generated from hubbard_params(T, U, V).
    #
    # Here:
    #
    # - T is the temperature,
    # - U is the on-site interaction,
    # - V is the non-local interaction.
    #
    # The returned file contains continued-fraction coefficients for several
    # resolvents, such as phase, Higgs, CDW, and AFM channels.

    if name == "sc_cdw":
        # Load the superconducting-charge-density-wave test case.
        #
        # Parameters:
        # T = 0.0
        # U = -2.5, attractive on-site interaction
        # V = 0.0
        return data_loader.load_panda(
            model="hubbard",
            subdir=SQUARE_DIR,
            file="resolvents.json.gz",
            **hubbard_params(0.0, -2.5, 0.)
        )

    elif name == "sc":
        # Load the superconducting test case.
        #
        # This uses the cube geometry/data set and a small attractive
        # non-local interaction V.
        return data_loader.load_panda(
            model="hubbard",
            subdir=CUBE_DIR,
            file="resolvents.json.gz",
            **hubbard_params(0.0, -2.5, -0.1)
        )

    elif name == "cdw":
        # Load the charge-density-wave test case.
        #
        # Compared to the SC case, this uses a positive V, which favors
        # charge-density-wave ordering.
        return data_loader.load_panda(
            model="hubbard",
            subdir=SQUARE_DIR,
            file="resolvents.json.gz",
            **hubbard_params(0.0, -2.5, 0.1)
        )

    elif name == "afm":
        # Load the antiferromagnetic test case.
        #
        # Here U is positive, corresponding to a repulsive on-site interaction,
        # which favors antiferromagnetic correlations.
        return data_loader.load_panda(
            model="hubbard",
            subdir=SQUARE_DIR,
            file="resolvents.json.gz",
            **hubbard_params(0.0, 2.5, 0.1)
        )

    else:
        # Stop with a clear error if the user requests an unknown test case.
        raise ValueError("Hubbard test: Invalid index")


def create_plot(name):
    """Load one Hubbard test case and plot several spectral densities."""

    # Load the selected Hubbard data set.
    #
    # pd_data contains metadata, continuum boundaries, and continued-fraction
    # coefficients for the available resolvents.
    pd_data = load_data(name)

    # Create a ContinuedFraction evaluator.
    #
    # The object reads the continued-fraction coefficients from pd_data and
    # uses them to evaluate the corresponding resolvents.
    resolvents = cf.ContinuedFraction(pd_data)

    fig, ax = plt.subplots()
    ax.set_ylim(-0.05, 1.)
    ax.set_xlabel(r"$\omega / t$")
    ax.set_ylabel(r"$\mathcal{A} (\omega)  / t^{-1}$")
    colors = np.array((
        "blue",
        "orange",
        "black",
        "limegreen",
        "deepskyblue",
        "magenta"
    ))
    linestyles = ["-", "-.", "--", "-", "--", ":"]

    # Construct the frequency grid at which the spectral functions are
    # evaluated.
    #
    # The grid starts slightly below zero and extends somewhat beyond the upper
    # continuum boundary:
    #
    #     upper limit = continuum upper edge + 0.3
    #
    # dtype=complex is used because the continued fraction is evaluated at
    # complex frequencies.
    w_lin = np.linspace(
        -0.01,
        pd_data["continuum_boundaries"][1] + 0.3,
        5000,
        dtype=complex
    )

    # Add a small positive imaginary part.
    #
    # This corresponds to evaluating the retarded response slightly above the
    # real frequency axis:
    #
    #     omega -> omega + i eta
    #
    # with eta = 1e-6.
    #
    # Numerically, this also gives narrow broadening to sharp peaks.
    w_lin += 1e-6j

    # Plot the superconducting phase spectral density.
    #
    # The resolvent name "phase_SC" refers to the superconducting phase channel.
    ax.plot(
        w_lin.real,
        resolvents.spectral_density(w_lin, "phase_SC"),
        label="Phase",
        ls=linestyles[0],
        c=colors[0]
    )

    # Plot the superconducting amplitude, or Higgs, spectral density.
    #
    # The resolvent name "amplitude_SC" refers to the amplitude channel of the
    # superconducting order parameter.
    ax.plot(
        w_lin.real,
        resolvents.spectral_density(w_lin, "amplitude_SC"),
        label="Higgs",
        ls=linestyles[1],
        c=colors[1]
    )

    # Plot the charge-density-wave amplitude spectral density.
    #
    # This channel describes CDW amplitude fluctuations.
    ax.plot(
        w_lin.real,
        resolvents.spectral_density(w_lin, "amplitude_CDW"),
        label="CDW",
        ls=linestyles[2],
        c=colors[2]
    )

    # Plot the longitudinal antiferromagnetic amplitude spectral density.
    #
    # The label "l.AFM" stands for longitudinal AFM.
    ax.plot(
        w_lin.real,
        resolvents.spectral_density(w_lin, "amplitude_AFM"),
        label="l.AFM",
        ls=linestyles[3],
        c=colors[3]
    )

    # Plot the transverse antiferromagnetic amplitude spectral density.
    #
    # The label "t.AFM" stands for transverse AFM.
    ax.plot(
        w_lin.real,
        resolvents.spectral_density(w_lin, "amplitude_AFM_transversal"),
        label="t.AFM",
        ls=linestyles[4],
        c=colors[4]
    )

    # Shade the continuum region.
    #
    # The continuum boundaries are taken from the loaded data and visualized as
    # a grey vertical span. Since no scale factor is passed here, the continuum
    # is shown in the same units as the frequency grid, i.e. omega / t.
    resolvents.mark_continuum(ax)

    ax.set_title(f"Hubbard: {name}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"test_result_{name}.pdf")


# Command-line interface:
#
# If the script is run with an additional argument, that argument is interpreted
# as the Hubbard test case to plot.
#
# Example:
#
#     python plot_hubbard.py sc
#
# or:
#
#     python plot_hubbard.py afm
#
# If no argument is provided, print a short instruction instead.
if len(sys.argv) > 1:
    create_plot(sys.argv[1])
else:
    print("Please provide the name of the test you would like to plot.")