# bea.jl: contains numerical computations using triangular elements
# About
#   Author       	- Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      	- Created December 2014
# References
#   Shape functions	- Wu, T.W., "Boundary element acoustics", WIT Press, 2000,ISBN 1-85312-570-9, 2005.
#   Meshes (gmsh)	- Geuzaine, C. and Remacle, J.-F., "Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities", International Journal for Numerical Methods in Engineering, volume 79, number 11, pages 1309-1331, 2009
#   Gauss-Legendre      - http://pomax.github.io/bezierinfo/legendre-gauss.html

# Using 6-point tensor Gauss-Legendre
qno = [0.6612093864662645, 0.2386191860831969, 0.9324695142031521]
wei = [0.3607615730481386, 0.4679139345726910, 0.1713244923791704]

# Interpolate on a triangle over the first dimension
# Input
#	xi	- Coordinates in the parameter domain
#	c	- Coefficients
# Output
#		- Triangular interpolation
function triang(xi,c,result)
	for i = 1:length(result)
		result[i] = c[1,i]*xi[1]+c[2,i]*xi[2]+c[3,i]*(1-xi[1]-xi[2])
	end # Does not seem possible: c[1,:]*xi[1]+c[2,:]*xi[2]+c[3,:]*(1-xi[1]-xi[2])
end

# Calculate the Jacobian of a triangular element to xi
# Input
#	x,y,z	- Triangle is between (x[1],y[1],z[1]), (x[2],y[2],z[2]) and (x[3],y[3],z[3])
# Output
#		- The Jacobian
function jacob(x,y,z)
	return sqrt( ( (y[1]-y[3])*(z[2]-z[3]) - (z[1]-z[3])*(y[2]-y[3]))^2 + ( (z[1]-z[3])*(x[2]-x[3]) - (x[1]-x[3])*(z[2]-z[3]))^2 + ( (x[1]-x[3])*(y[2]-y[3])- (y[1]-y[3])*(x[2]-x[3]) )^2 )/4
end

# Read in an input file
# Input
#	name	- Name of the file
# Output
#	xyz	- x,y and z coordinates of nodes
# 	elmts	- Element connectivity; for example, elmt[45,:] = [3 89 462] means triangle 45 has node 3 as first node etc
function readMesh(name)
	f = open(name)
	if (readline(f) != "\$MeshFormat\n") | (readline(f) != "2.2 0 8\n") | (readline(f) != "\$EndMeshFormat\n") | (readline(f) != "\$Nodes\n")
		print(STDERR, " Wrong mesh type, check whether file ", name, " starts with $MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n")
		exit()
	end
	nnode = int(readline(f))
	xyz = zeros(nnode,3)
	for ni = 1:nnode
		i = 1
		q = readline(f)
		while(q[i] != ' ') i = i+1 end
		if(int(q[1:i]) != ni) # Actually q[1:(i-1)] but same result
			print(STDOUT, "Wrong mesh format in ", name, ": expected index ", ni, " for node but got ", int(q[1:(i-1)]), ".\n")
			quit()
		end
		for dim = 1:3
			i = i+1
			start = i
			while( (q[i] != ' ') & (q[i] != '\n') ) i = i+1 end
			xyz[ni,dim] = float(q[start:i])
		end
		if(abs(norm(xyz[ni,:])-1) > 1e-6)
			print(STDERR, "Ignore entry ", ni, xyz[ni,:], '\n')
		end
	end
	if (readline(f) != "\$EndNodes\n") | (readline(f) != "\$Elements\n")
		print(STDERR, "Wrong mesh format in ", name, " after nodes.\n")
		exit()
	end
	nelem = int(readline(f)) # Not final, ignore junk elements when returning
	elmts = zeros(Int64,nelem,3)
	ei = 0
	for ni = 1:nelem
		i = 1
		q = readline(f)
		while(q[i] != ' ') i = i+1 end
		i = i+1
		start = i
		while(q[i] != ' ') i = i+1 end
		if(q[start:(i-1)] != "2")
			continue
		end
		# Now only type 2 which is 3D triangular elements
		ei = ei+1
		i = i+1
		start = i
		while(q[i] != ' ') i = i+1 end
		nbtag = int(q[start:i])
		for ti = 1:nbtag
			i = i+1
			while(q[i] != ' ') i = i+1 end
		end
		for ti = 1:3
			i = i+1
			start = i
			while( (q[i] != ' ') & (q[i] != '\n') ) i = i+1 end
			elmts[ei,ti] = int(q[start:i])-1 # because we delete origin
		end
	end
	close(f)
	return (xyz,elmts[1:ei,:])
end

# Make the system matrix for triangular elements.
# Input
#   xyz		- Mesh
#   elmts 	- Element connectivity
#   k 		- Wave number
# Output
#   A		- The system matrix
function getSystemMatrix(xyz, elmts, k)
	nnode = size(xyz,1)
	A = zeros(Complex{Float64}, nnode, nnode)
	for ni = 1:nnode
		if (ni%10 == 0)
			print(ni/nnode, " ")
		end
		updateA(xyz,elmts,k,ni,A)
	end
	return A
end

# Update row of system matrix
# Input
#   xyz		- Mesh
#   elmts 	- Element connectivity
#   k 		- Wave number
#   ni		- Row to update
#   A		- System matrix
function updateA(xyz, elmts,k, ni, A)
	nelem = size(elmts,1)
	xk = [xyz[ni,:]; xyz[ni,:]; xyz[ni,:]]
	J = 0.0
	common = 0.0+0.0im
	xm = zeros(size(xk))
	indir = [1, 2, 3]
	tri = [0.0 0.0 0.0]
	for ei = 1:nelem
		xm[:,:] = xyz[vec(elmts[ei,:]),:]
		# Sort so that third node of the element (xi1=0=xi2) is closest to or on singularity (xm = xyz[ni,:])
		(val,idx) = findmin(sqrt( sum((xm-xk).^2,2) ) )
		indir[:] = [1, 2, 3]
		indir[idx] = 3
		indir[3] = idx
		J = jacob(xm[indir,1], xm[indir,2], xm[indir,3])/4
		for qi = 1:3, ri=1:3, si = -1:2:1, ti = -1:2:1
			triang([(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2), (1+si*qno[qi])*(1+ti*qno[ri])/4],xm[indir,:],tri)
			common = wei[qi]*wei[ri]*exp(1im*k*norm(tri -xk[1,:]))/(4*pi)/sqrt(sum( (xm[indir[1],:]+(1+ti*qno[ri])/2*(xm[indir[2],:]-xm[indir[1],:]) -xm[indir[3],:]+(xm[indir[3],:]-xk[1,:])/(1+si*qno[qi])*2).^2,2)[1] )*J

			A[ni,elmts[ei,indir[1]] ] += common*((1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2))
			A[ni,elmts[ei,indir[2]] ] += common*(1+si*qno[qi])*(1+ti*qno[ri])/4
			A[ni,elmts[ei,indir[3]] ] += common*(1-(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2)- (1+si*qno[qi])*(1+ti*qno[ri])/4)
		end
	end
end

# Compute potential from density. 
# Input
#   xyz		- Mesh
#   elmts 	- Element connectivity
#   k 		- Wave number
#   sing	- Point at which to evaluate the potential
#   v		- Known density
# Output
#   pot		- Potential at point
function potential(xyz, elmts,k, sing, v)
	nelem = size(elmts,1)
	pot = 0+0im
	J = 0.0
	xm = zeros(size([sing;sing;sing]))
	indir = [1, 2, 3]
	tri = [0.0 0.0 0.0]
	for ei = 1:nelem
		xm[:,:] = xyz[vec(elmts[ei,:]),:]
		# Sort so that third node of the element (xi1=0=xi2) is closest to or on singularity (xm = xyz[ni,:])
		(val,idx) = findmin(sqrt( sum((xm-[sing;sing;sing]).^2,2) ) )
		indir[:] = [1, 2, 3]
		indir[idx] = 3
		indir[3] = idx

		J = jacob(xm[indir,1], xm[indir,2], xm[indir,3])/4
		for qi = 1:3, ri=1:3, si = -1:2:1, ti = -1:2:1
			triang([(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2), (1+si*qno[qi])*(1+ti*qno[ri])/4],xm[indir,:],tri)
			pot += wei[qi]*wei[ri]*exp(1im*k*norm(tri -sing))/(4*pi)/sqrt(sum( (xm[indir[1],:]+(1+ti*qno[ri])/2*(xm[indir[2],:]-xm[indir[1],:]) -xm[indir[3],:]+(xm[indir[3],:]-sing)/(1+si*qno[qi])*2).^2,2)[1] )*J*( v[elmts[ei,indir[1]]]*(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2)     + v[elmts[ei,indir[2]]]*(1+si*qno[qi])*(1+ti*qno[ri])/4    +   v[elmts[ei,indir[3]]]*(1-(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2)- (1+si*qno[qi])*(1+ti*qno[ri])/4)   )
		end
	end
	return pot
end

# Interpolate the grid somewhere by taking the closest mesh point
# Input
#   where		- Coordinate where to interpolate
#   c			- Coefficient vector of density solution
#   xyz			- Mesh
# Output
#   			- Interpolated value
function interpolateGrid(where,c,xyz)
	(x, y, z) = where
	res = 0.0im*x + 0.0*y + 0.0im*z
	for pt = 1:length(x)
		prevR = 2
		for ni =1:length(c)
			r = sqrt( (xyz[ni,1] -x[pt])^2 + (xyz[ni,2] -y[pt])^2 + (xyz[ni,3] -z[pt])^2)	
			if(r < prevR) # just take value of closest by
				prevR = r
				res[pt] = c[ni]
			end
		end
	end
	return res
end

