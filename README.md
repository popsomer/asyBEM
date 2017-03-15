# asyBEM
Exploit high-frequency scattering asymptotics in a Boundary Element Method.

The script makeVb1.m in the folder `matlab' constructs the eigenvector of a full cycle of reflections between two circles. This is used by the script phase2Circles.m to validate our symbolic results. The Sagemath worksheet LimitCycleTwoCircles.sws computes these symbolic coefficients.

There is also some code to apply asymptotic compression in 2D in the matlab folder (see https://github.com/popsomer/bempp for validation), and in 3D in the julia folder.

Tested with Matlab 2016b, Sage 6.9 and Julia 0.4.7-pre.
