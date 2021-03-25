# bubbles_multi/3D_finite

`bubbles_multi/3D_finite` is a package using the capacitance matrix approximation and the multipole expansion method to study Helmholtz resonance problems in three dimensions. Specifically, it handles scattering by a finite number of spherical material inclusions.

Please cite the following reference when using this source:

[1] Ammari H, Davies B & Hiltunen E O. 2021 *Functional analytic methods for discrete approximations of subwavelength resonator systems*. To appear.

## Multipole expansion method

In the case of spherical resonators, the solutions to scattering problems can be expanded in terms of spherical harmonics. This approach is explained in detail in the appendix of [2].

## Capacitance matrix

A significant improvement in computational efficiency can be achieved by using the eigenvalues of the generalized capacitance matrix to approximate the resonant frequencies. The function `capacitance.m` use the multipole expansion method (with just one multipole) to approximate the capacitance coefficients associated to the resonators' geometry. We can then compute the eigenvalues of the generalized capacitance matrix and use the fact that <img src="https://latex.codecogs.com/svg.latex?\omega_j=\sqrt{\lambda_j}+O(\delta)"> to find the resonant frequencies.

## Dilute approximation of the capacitance matrix

In the case that the resonators are small compared to the distance between them, we have an explicit formula for the capacitance matrix at leading order. In particular, if the resonators are given by

<img src="https://latex.codecogs.com/svg.latex?\large&space;D=\bigcup_{j=1}^N \left(\epsilon B + z_j\right),">

then we have that

<img src="https://latex.codecogs.com/svg.latex?\large&space;C_{ij}%20=%20\begin{cases}\epsilon\mathrm{Cap}_B%20+%20O(\epsilon^3),%20&\quad%20i=j,\\-\epsilon^2%20\frac{(\mathrm{Cap}_B)^2}{4\pi%20|z_i-z_j|}%20+%20O(\epsilon^3),%20&\quad%20i\neq%20j.\end{cases}">

# Notes on the current version

The current version of the multipole method assumes that the resonators are collinear, in the sense that their centres all lie on a straight line.


# References

[2] Ammari H, Davies B, Hiltunen E O and Yu S. 2020 *Topologically protected edge modes in one-dimensional chains of subwavelength resonators.* Journal de Mathématiques Pures et Appliquées 144: 17-49. (https://doi.org/10.1016/j.matpur.2020.08.007)
