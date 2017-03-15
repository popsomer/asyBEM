# Test the simple 3D BEM implementation and asymptotic compression
# See popsomer/bempp/meshes for additional meshes
using PyCall 
@pyimport matplotlib.pyplot as plt
nrMesh = 0
k = 15

include("bea.jl")
include("window.jl")

(xyz, elmt) = readMesh(string("sphere", nrMesh, ".msh") )
# Prints to ignore first enty which is origin:
xyz = xyz[2:end,:]

A = @time(getSystemMatrix(xyz,elmt,k))
b = -exp(1im*k*xyz[:,1])

c = A\b
f = open(string("ck", k, "sph", nrMesh), "w") # save results
writedlm(f,[real(c) imag(c)])
close(f)
resid = A*c-b

print("Potential at -1 0 0 is ", potential(xyz,elmt,k,[-1 0 0],c), k, " =k, nrMesh = ", nrMesh, "\n\n")
checkSolution(xyz,elmt,k,c)


# Compare with compressed system matrix
Acom = @time(getComprSystemMatrix(xyz,elmt,k,cuto))
ccom = Acom\b
f = open(string("windowck", k, "sph", nrMesh), "w")
writedlm(f,[real(ccom) imag(ccom)])
close(f)
resid = A*ccom-b
residcom = Acom*ccom-b

print(resid[1:10], "=resid, residcom = ", residcom[1:10], '\n')
print("Potential at -1 0 0 is ", potential(xyz,elmt,k,[-1 0 0],ccom), "\n\n")
checkSolution(xyz,elmt,k,ccom)
