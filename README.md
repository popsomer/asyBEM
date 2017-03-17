# asyBEM
Exploit high-frequency scattering asymptotics in a Boundary Element Method. Tested with Matlab 2016b (www.matlab.com), Chebfun 5.6.0 (www.chebfun.org), Sage 6.9 (www.sagemath.org) and Julia 0.4.7-pre (www.julialang.org).

This implementation accompanies the article "Coupling modes in high-frequency multiple scattering problems: the case of two circles", D. Huybrechs and P. Opsomer, (in preparation). For the case of two circular scatterers, we compute a Taylor approximation of the equilibrium phase of the density in ray tracing in the Sagemath worksheet LimitCycleTwoCircles.sws. The Taylor coefficients can be computed independently of the wavenumber and incident wave.

This phase is also the eigenvector of a matrix representing a full cycle of reflections, computed in the script makeVb1.m in the folder `matlab'. After adding the path to chebfun, the script phase2Circles.m validates our symbolic results and the geometric interpretation of the stationary point and phase: the latter is the distance to the periodic orbit via an infinite number of reflections at equal angles, where one substracts the distance between the circles at each reflection.

There is also some code to apply asymptotic compression in 2D in the matlab folder (github.com/popsomer/bempp for validation), and in 3D in the julia folder.


