# Tolman-Oppenheimer-Volkoff equation solver

## PYTOV:
PYTOV is a simple Python Tolman-Oppenheimer-Volkoff (TOV) equation integrator.

The TOV equations

<p align="center"><img alt="$$&#10;\frac{dP}{dr} = - \frac{G}{(cr)^{2}} \frac{(\epsilon (r) + P(r))(M(r) + 4\pi r^{3} \frac{P(r)}{c^{2}})}{(1 - \frac{2G M(r)}{c^{2}r})}&#10;$$" src="https://github.com/dlwhittenbury/PYTOV/blob/master/svgs/df0261d65270bc97c0b09d7835b14a41.svg" align="middle" width="328.35165pt" height="49.145085pt"/></p>


### Summary:
-   A high density equation of state file and a low density equation of state file will be read and then combined in a simple manner.

-   Logarithmic interpolation of the combined equation of state is used.

-   A simple fixed step fourth order Runge-Kutta method is used to integrate the TOV equations. A central density is specified and the TOV equations are integrated out to the surface of the compact star. This is repeated for a range of densities.

-   The required data for a mass vs radius curve is output to a file along with a file containing the details of the maximum mass compact star.



# References
1. Reference 1
2. Reference 2
3. Reference 3
