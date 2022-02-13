#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

PYTOV:
======

    PYTOV is a simple Python Tolman-Oppenheimer-Volkoff (TOV) equation 
    integrator. 
    
    Summary:
    ========    
    -   A high density equation of state file and a low density equation of 
        state file will be read and then combined in a simple manner. 
    -   Logarithmic interpolation of the combined equation of state is used. 
    -   A simple fixed step fourth order Runge-Kutta method is used to 
        integrate the TOV equations. A central density is specified and the 
        TOV equations are integrated out to the surface of the compact star. 
        This is repeated for a range of densities. 
    -   The required data for a mass vs radius curve is output to a file along
        with a file containing the details of the maximum mass compact star.
        
@author: D. L. Whittenbury
"""

from re import L
import numpy as np
import matplotlib.pyplot as plt
import os

# Useful constants
from constants import MeVfm_3Tokm_2, km_2ToMeVfm_3, kmToModot, ModotTokm


def read_eos(
    filename: str,
    order: dict[str, int] = {"density": 0, "pressure": 1, "energy_density": 2},
):
    """
    Simple EoS file reader. The eos is output in the default order density, pressure, energy density.

    filename: Filename of EoS to be read.
    order: Dictionary with keys density, pressure, energy_density and the
        values indicate which column number, e.g.,
        {"density":0, "pressure":1, "energy_density":2}.

    """

    # Paths to files
    current_dir = os.getcwd()
    INPUT_PATH = os.path.join(current_dir, "input")

    in_eos = np.loadtxt(os.path.join(INPUT_PATH, filename))
    eos = np.zeros(in_eos.shape)

    # Density
    eos[:, 0] = in_eos[:, order["density"]]
    eos[:, 1] = in_eos[:, order["pressure"]]
    eos[:, 2] = in_eos[:, order["energy_density"]]

    return eos


def join_eos(low_density_eos: np.ndarray, high_density_eos: np.ndarray):
    """
    Simple EoS joiner.

    """

    i1 = np.argmax(low_density_eos[:, 0] > high_density_eos[0, 0])
    i2 = np.argmax(low_density_eos[:, 1] > high_density_eos[0, 1])
    i3 = np.argmax(low_density_eos[:, 2] > high_density_eos[0, 2])

    index = min(i1, i2, i3)

    selected_low_density_eos = low_density_eos[:index, :]
    combined_eos = np.row_stack((selected_low_density_eos, high_density_eos))

    return combined_eos


def density(P):
    """Calculate the baryonic density as a function of pressure

    INPUT: Pressure [km^-2]

    OUTPUT: Baryonic density [fm^-3 ]

    NOTES: Only interpolation and no extrapolation.

    """
    if P < 0:
        return 0

    else:

        i = np.argmax(peos > P)
        # Finds the index of the first p in peos > P

        if i == 0:
            return neos[0]  # No extrapolation to lower density/pressure used

        else:

            last = len(neos) - 1
            # The last index

            if i <= last:  # Then interpolate

                den = neos[i - 1] * np.exp(
                    np.log(P / peos[i - 1])
                    * np.log(neos[i - 1] / neos[i])
                    / np.log(peos[i - 1] / peos[i])
                )

                return den

            else:
                return neos[last]  # No extrapolation to higher density/pressure


def energydensity(P):

    """Calculate the energy density as a function of pressure

    INPUT: Pressure [km^-2]

    OUTPUT: Energy density [km^-2]

    NOTES: Only logarithmic interpolation and no extrapolation.

    """
    if P < 0:
        return 0
    else:

        i = np.argmax(peos > P)
        # Finds the index of the first p in peos > P

        if i == 0:
            return eeos[0]  # No extrapolation to lower en. density/pressure used

        else:

            last = len(neos) - 1
            # The last index

            if i <= last:  # Then interpolate

                eden = eeos[i - 1] * np.exp(
                    np.log(P / peos[i - 1])
                    * np.log(eeos[i - 1] / eeos[i])
                    / np.log(peos[i - 1] / peos[i])
                )

                return eden

            else:
                return eeos[
                    last
                ]  # No extrapolation to higher en. density/pressure used


def pressure(n):
    """Calculate the pressure as a function of baryonic density

    INPUT: Baryonic density [fm^-3 ]

    OUTPUT: Pressure [km^-2]

    NOTES: Only logarithmic interpolation and no extrapolation.

    """
    if n < 0:
        return 0

    else:

        i = np.argmax(neos > n)
        # Finds the index of the first n in neos > n

        if i == 0:
            return peos[0]  # No extrapolation to lower density/pressure used

        else:

            last = len(neos) - 1
            # The last index

            if i <= last:  # Then interpolate

                press = peos[i - 1] * np.exp(
                    np.log(n / neos[i - 1])
                    * np.log(peos[i - 1] / peos[i])
                    / np.log(neos[i - 1] / neos[i])
                )

                return press

            else:
                return peos[last]  # No extrapolation to higher density/pressure used


def TOVderivs(y, t):
    """Calculate the derivatives in the TOV equations

    INPUT:  Dependent variables y = y(t), where y is gravitational mass M,
            graviational field Phi, pressure P and baryon number of baryons A.
            The independent variable is r.

    OUTPUT: Derivatives of dependent variables mass M, graviational field Phi,
            pressure P and baryon number of baryons A.

    NOTES:  Uses the functions density() and energydensity() which utilise
            logarithmic interpolation.

    """

    # Radius in km
    r = t
    r2 = r**2
    r3 = r**3

    # The dependent variables y
    M = y[0]  # Mass in km , note 1 solar mass = 1.4766 km

    # Phi = y[1]; # Gravitational potential, later match to suface

    P = y[2]  # Pressure km^-2

    # A = y[3] # Number of baryons

    # Lookup energy density and baryonic density
    nB = density(P)  # fm^-3
    E = energydensity(P)  # km^-2

    # Array of derivative equations
    dydx = np.zeros(4)

    # Eq. for graviational mass
    dydx[0] = 4.0 * np.pi * r2 * E

    # Gravitational potential - to be continued ... Needs to be matched!!
    dydx[1] = (M + 4.0 * np.pi * r3 * P) / (r2 * (1.0 - 2.0 * M / r))

    # Hydrostatic equilibrium ... equation for pressure
    dydx[2] = -(E + P) * dydx[1]

    # Eq. for baryon number A
    dydx[3] = (4.0 * np.pi * r2 * nB) / np.sqrt(1.0 - 2.0 * M / r)

    return dydx


def RK4Step(u, t, dt, derivs):
    """Fourth order Runge Kutta method

    INPUT:  Dependent variables u, independent variable t, step size dt,
            function to calculate derivatives of dependent variables
            with respect to t, i.e., expect derivs(u,t).

    OUTPUT: Advanced independent variables unew.

    NOTES:  See Applied numerical methods for engineers and scientists,
            S. S. Rao, page 659.

    """

    n = len(u)
    unew = np.zeros(n)
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)

    k1 = derivs(u, t)
    k2 = derivs(u + 0.5 * dt * k1, t + 0.5 * dt)
    k3 = derivs(u + 0.5 * dt * k2, t + 0.5 * dt)
    k4 = derivs(u + dt * k3, t + dt)

    unew = u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    return unew


def TOV_initial(rs: float, nB: float):
    """
    Constructs the initial vector y for the dependent variables, gravitational
    mass M, graviational field Phi, pressure P and baryon number of baryons A.

    rs: Small radial offset [km].
    nB: Total baryonic density [fm^-3].

    """

    # Dependent variables
    y = np.zeros(4)

    # We are not at the centre of the star, but rather slightly off centre,
    # starting at rs.
    # Gravitational mass M[r=0] = 0, but we are not at the centre so we have a
    # small contribution
    y[0] = (4.0 * np.pi / 3.0) * energydensity(pressure(nB)) * rs**3

    # Graviational field Phi
    # Phi[r=0] is not known. Actually phi will be matched at the surface, so we
    # could set the y[1] = 0 and there would be no difference! BUT we can also
    # choose phi0 = 0 at the centre and add a second term which is the
    # correction because we are not starting at the centre.
    y[1] = (
        0.0
        + (2.0 * np.pi / 3.0)
        * (energydensity(pressure(nB)) + 3.0 * pressure(nB))
        * rs**2
    )

    # Pressrue P[r=0] = P0, but we are not at the centre, so we have a small
    # additional contribution
    y[2] = pressure(nB)
    -(2.0 * np.pi / 3.0) * (energydensity(pressure(nB)) + pressure(nB)) * (
        energydensity(pressure(nB)) + 3.0 * pressure(nB)
    ) * rs**2

    # Number of baryons A[r=0] = A0 = 0 , but we are not at the centre so we
    # have a small number of baryons
    y[3] = (4.0 * np.pi / 3.0) * nB * rs**3

    return y


def IntegrateTOV(uStart, tStart, dt, derivs):
    """Integrate the TOV equations

    INPUT:  Starting values for the dependent variables uStart, starting
            value for the dependent variable tStart, step size dt, the
            function to evaluate derivatives of the dependent variables u
            with respect to the independent variable t, i.e., expect
            derivs(u,t).

    OUTPUT: The variable out which contains (R, MG, Phi, P, Ab),
            radius R, gravitational mass MG, graviational field Phi (but
            not matched yet), pressure P (when integration stopped)),
            number of baryons Ab.

    NOTES:

    """

    # Maximum no. of steps
    step_limit = 100000

    n = len(uStart)

    out = np.zeros(n + 1)  # r and n variables

    Var = np.zeros(n)
    Var = uStart  # Starting variable values

    r = tStart  # Starting r
    dr = dt  # Step size

    # Solutions to be stored here
    sol = np.zeros((step_limit, n))
    sol[0] = Var

    counter = 0
    while sol[counter, 2] > peos[0] and counter < step_limit - 1:
        # Stopping criteria: Continue while pressure is greater than the smallest
        #                    value in the EoS table and while the no. of steps is
        #                    less than the chosen limit. No extrapolation to lower
        #                    values of pressure, we take the smallest value as
        #                    approx. defining the location of the surface.
        Var = RK4Step(Var, r, dr, derivs)  # RK4 step forward
        counter = counter + 1
        sol[counter] = Var
        r = r + dr

    if sol[counter, 2] < 0:  # Do not want a negative pressure solution so subtract 1
        counter = counter - 1

    # out = (r, sol[counter,0], sol[counter,1], sol[counter,2], sol[counter,3]) =(R, MG,Phi, P, Ab)
    out[0] = r
    out[1 : n + 1] = sol[counter]

    return out


def determine_max_mass(
    masses: np.ndarray, radii: np.ndarray, rhocs: np.ndarray, phis: np.ndarray
):
    """
    A simple determination of the maximum mass.

    One could of course be more accurate and interpolate and find a maximum.

    """

    # Maximum mass
    mass_max = max(masses)

    # Index of maximum mass neutron star
    idx = [i for i, v in enumerate(masses) if v == mass_max]

    # The corresponding radius
    radius_max = (radii[idx])[0]

    # The corresponding central density
    rho_c_max = rhocs[idx][0]

    # The gravitational field on the surface
    phi_max = (phis[idx])[0]

    return mass_max, radius_max, rho_c_max, phi_max


if __name__ == "__main__":

    # PREPARE THE EQUATION OF STATE(S)
    # ================================

    # Reads in EoS. If order is not specified, then it is assumed to be
    # density, pressure and energy density. read_eos outputs the eos in the
    # default order density, pressure, energy density.
    low_density_eos = read_eos(filename="bps.dat")
    high_density_eos = read_eos(filename="beta_eos.dat")

    # When you combine EoS they should be in same order and have the same units. Unit conversion is on you, the user.
    combined_eos = join_eos(
        low_density_eos=low_density_eos, high_density_eos=high_density_eos
    )

    # Separate EoS and convert units if necessary for solving the TOV equations. Unit conversion is on you, the user.
    # Density should be fm^-3
    # Pressure and energy density should be km^-2
    neos = combined_eos[:, 0]
    peos = combined_eos[:, 1] * MeVfm_3Tokm_2
    eeos = combined_eos[:, 2] * MeVfm_3Tokm_2

    # Select a central density and integrate out to the surface. Repeat for a
    # range of densities. This will give you the required data for a mass vs
    # radius curve.

    # SOLVE THE TOV EQUATIONS
    # =======================

    # Saturation density
    n0 = 0.16

    # Starting density
    startingDensity = 0.8 * n0

    # Stepsize in density
    stepsize = 0.02 * n0

    # central density loop
    # Number of iterations (0.8*0.16 + 0.02*0.16*335 = 1.2 [fm^-3])
    iterations = 335

    # Start off centre
    rs = 0.00001  # [km]
    # Stepsize (APPROX. TIME: dr = 0.01 about ~1 min, dr = 0.001 about ~10 min )
    # Best to use 0.001 or smaller when you require accurate results. It all
    # depends on how long you are willing to wait. Although, you could improve this
    # integrator by upgrading to a variable step size integrator.
    dr = 0.001  # 0.01; 0.001; 0.0001;

    radii = np.zeros(iterations)
    masses = np.zeros(iterations)
    phis = np.zeros(iterations)
    rhocs = np.zeros(iterations)

    # Loop over densities
    for i in range(iterations):

        # Central density
        nB = startingDensity + stepsize * float(i)

        # Central density
        rhocs[i] = nB

        y = TOV_initial(rs, nB)

        # Solutions (r, M, Phi (not matched yet), P, A)
        solution = IntegrateTOV(y, rs, dr, TOVderivs)

        # Convert to correct units before output to file
        radii[i] = solution[0]
        masses[i] = solution[1] * kmToModot
        phis[i] = solution[2]

        # Not yet matched at the surface
        # Rescale phi to match at surface with Schwarzchild solution
        const = 0.5 * np.log(1.0 - 2.0 * masses[i] * ModotTokm / radii[i]) - phis[i]
        phis[i] = phis[i] + const

        Pre = solution[3] * km_2ToMeVfm_3
        Ab = solution[4] * 1.0e54  # fm^-3 <-> km^-3

        # Print to screen as calculation progresses
        print(i, nB, radii[i], masses[i])

    # Determine maximum mass configuration
    mass_max, radius_max, rho_c_max, phi_max = determine_max_mass(
        masses=masses, radii=radii, rhocs=rhocs, phis=phis
    )

    # PROCESS THE OUTPUT: CREATE DATA FILES AND FIGURES
    # =================================================

    # Save calculations to file
    # =========================

    # Output to file: compact_stars.dat
    OUT = np.c_[rhocs, radii, masses]
    np.savetxt(os.path.join("output", "compact_stars.dat"), OUT, fmt="%1.6e")

    # Output to file: max_mass_star.dat
    NAMES = np.array([r"# \rho_c [fm^{-3}]", "# R [km]", r"# M [M_{odot}]"])
    DAT = np.array([rho_c_max, radius_max, mass_max])
    OUT = np.zeros(NAMES.size, dtype=[("var1", float), ("var2", "U20")])
    OUT["var1"] = DAT
    OUT["var2"] = NAMES
    np.savetxt(
        os.path.join("output", "max_mass_star.dat"),
        OUT,
        delimiter="    ",
        fmt="%-1.6e %-20s",
    )
    # the dash `-` in the format specifer tells python to align to the left

    # Print to screen
    # ===============

    print("\n\nMax mass determined without interpolation!")
    print("Mmax: ", mass_max, " [M_sol]")
    print("Rmax: ", radius_max, " [km]")

    # Figures
    # =======

    fs = 15  # fontsize

    # Mass vs radius
    fig1 = plt.figure(figsize=(8, 6))
    plt.plot(radii, masses, "b*-")
    plt.plot(radius_max, mass_max, "ro")
    plt.xlabel("Radius [km]", fontsize=fs)
    plt.ylabel(r"Mass [$M_{\odot}$]", fontsize=fs)
    plt.show()
    fig1.savefig(os.path.join("output", "mass_vs_radius.png"))

    # Pressue vs density
    fig2 = plt.figure(figsize=(8, 6))
    plt.plot(neos, peos * km_2ToMeVfm_3, "b-")
    plt.plot(rho_c_max, pressure(rho_c_max) * km_2ToMeVfm_3, "ro")
    plt.xlabel(r"Density $\rho$ [fm$^{-3}$]", fontsize=fs)
    plt.ylabel("Pressure [MeVfm$^{-3}$]", fontsize=fs)
    plt.show()
    fig2.savefig(os.path.join("output", "pressure_vs_density.png"))

    # Gravitational potential vs radius

    # Phi from the surface out to 20~km say
    rkm = np.linspace(radius_max, 20, 200)

    # Phi external
    phiExt = 0.5 * np.log(1.0 - np.divide(2.0 * mass_max * ModotTokm, rkm))

    fig3 = plt.figure(figsize=(8, 6))
    plt.plot(radius_max, phi_max, "r*")  # Graviational field a surface
    plt.plot(rkm, phiExt)  # External gravitational field
    plt.xlabel("Radius [km]", fontsize=fs)
    plt.ylabel(r"$\Phi$", fontsize=fs)
    plt.show()
    fig3.savefig(os.path.join("output", "grav_vs_radius.png"))
