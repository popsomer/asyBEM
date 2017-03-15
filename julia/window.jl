# window.jl: contains numerical test and windowed version of our simple 2D BEM
# About
#   Author       	- Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      	- Created December 2014
# References
#   Smooth cutoff function	- Huybrechs, Daan, "Phase extraction without phase extraction: asymptotic compression of dense BEM matrices", talk on high frequency scattering workshop, Reading, 16/12/2014.
include("bea.jl")

# Test the solution on BC, Helmholtz and Sommerfeld and save some density values
# Input
#	xyz	- Mesh
#	elmt	- Element connectivity
#	k	- Wave number
#   	v	- Known density
function checkSolution(xyz,elmt,k,c)
	rests = zeros(6,7)
	mrt = zeros(6,7)
	hh = 1e-5
	cdf = zeros(Complex,3,2)
	sommer = zeros(6,7)
	for ti = 1:size(sommer,1), ppi = 1:size(sommer,2)
		theta = ti*pi/(size(sommer,2)-1)
		phi = ppi*pi/(size(sommer,1)/2)
		x = sin(theta)*cos(phi)
		y = sin(theta)*sin(phi)
		z = cos(theta)
		rests[ti,ppi] = abs( (exp(1im*k*x)+potential(xyz,elmt,k,[x y z],c))/exp(1im*k*x))
		# Dont take pt on boundary when checking central differences
		pt = [x y z]*1.2
		potpt = exp(1im*k*pt[1])+potential(xyz,elmt,k,pt,c)
		for dim = 1:3, si = -1:2:1
			pt[dim] = pt[dim] + si*hh
			cdf[dim,int((si+3)/2)] = exp(1im*k*pt[1])+potential	(xyz,elmt,k,pt,c)
			pt[dim] = pt[dim] - si*hh
		end
		mrt[ti,ppi] = abs( (( (sum(sum(cdf,1),2)-6*potpt)[1])/hh^2+k^2*potpt)/(k^2*potpt))
		if mrt[ti,ppi] > 1
			print("Err ", ti, ppi, ", potpt and cdf are ", potpt, cdf, ", su = ", ( (sum(sum(cdf,1),2)-6*potpt)[1]), '\n')
		end
		pt = [sin(theta) sin(theta) cos(theta)]
		sommer[ti,ppi] = abs(potential(xyz,elmt,k,2^ppi*pt,c) )*2^ppi
	end
	print(norm(rests), "is Rel. error on boundary conditions and Rel. errors on central differences and going away are ", norm(mrt), sommer[3,:], '\n')

	nPointsTh = 4
	nPointsPh = 4
	theta = linspace(0.1,pi-0.1,nPointsTh)
	phi = linspace(pi/2,3*pi/2-0.1,nPointsPh)
	
	results= zeros(Complex,3+nPointsTh^2)
	results[1:3] = interpolateGrid( ([1,-1,-1/sqrt(3)],[0,0,-1/sqrt(3)], [0,0,1/sqrt(3)]), c, xyz)

	f = open(string("resultsBEJ", nrMesh, "k", k), "w")
	writedlm(f,k)
	write(f,"was k, ansatz is x, x&y&z are \n")
	writedlm(f,[1 0 0; -1 0 0; -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)])
	for nt = 1:nPointsTh
		x = sin(theta[nt])*cos(phi)
		y = sin(theta[nt])*sin(phi)
		z = cos(theta[nt])*ones(size(phi))
		results[(4+(nt-1)*nPointsTh):(3+(nt)*nPointsTh)] = interpolateGrid((x, y, z ), c, xyz)
		writedlm(f,[x y z])
	end
	write(f,"were x&y&z, results are \n")
	writedlm(f,[real(results) imag(results)])
	close(f)
end

# Apply window function
# Input
#	r	- Distance
#	m	- Cutoff
# Output
#		- 0 if r>m, 1 if r<0.7*m and else C^\infty-smooth
function window(r,m) 
	if (r <= 0.7*m) return 1; 
	elseif (r < m) return exp(2*exp(-(0.3*m)/(r-0.7*m) )/( (r-0.7*m)/(0.3*m)-1) ); 
	else return 0; 
	end
end

# Make the compressed system matrix for triangular elements.
# Input
#   xyz		- Mesh
#   elmts 	- Element connectivity
#   k 		- Wave number
#   m		- Cutoff
# Output
#   A		- The compressed system matrix
function getComprSystemMatrix(xyz, elmts, k, m)
	nnode = size(xyz,1)
	A = zeros(Complex{Float64}, nnode, nnode)
	for ni = 1:nnode
		if (ni%10 == 0)
			print(ni/nnode, " ")
		end
		if(xyz[ni,1] < -0.2) # Illuminated region
			updateComprA(xyz,elmts,k,ni,A,m);
		else
			updateA(xyz,elmts,k,ni,A);
		end
	end
	return A
end

# Update row of the compressed system matrix
# Input
#   xyz		- Mesh
#   elmts 	- Element connectivity
#   k 		- Wave number
#   ni		- Row to update
#   A		- Compressed system matrix
#   m		- Cutoff
function updateComprA(xyz, elmts,k, ni, A, m)
	nelem = size(elmts,1)
	xk = [xyz[ni,:]; xyz[ni,:]; xyz[ni,:]]
	for ei = 1:nelem
		xm = xyz[vec(elmts[ei,:]),:]
		# Sort so that third node of the element (xi1=0=xi2) is closest to or on singularity (xm = xyz[ni,:])
		(val,idx) = findmin(sqrt( sum((xm-xk).^2,2) ) )
		indir = [1, 2, 3]
		indir[idx] = 3
		indir[3] = idx
		J = jacob(xm[indir,1], xm[indir,2], xm[indir,3])/4
		for qi = 1:3, ri=1:3, si = -1:2:1, ti = -1:2:1
			dist = norm(triang([(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2), (1+si*qno[qi])*(1+ti*qno[ri])/4],xm[indir,:]) -xk[1,:])
			common = wei[qi]*wei[ri]*exp(1im*k*dist)/(4*pi)/sqrt(sum( (xm[indir[1],:]+(1+ti*qno[ri])/2*(xm[indir[2],:]-xm[indir[1],:]) -xm[indir[3],:]+(xm[indir[3],:]-xk[1,:])/(1+si*qno[qi])*2).^2,2)[1] )*J*window(dist,m)

			A[ni,elmts[ei,indir[1]] ] += common*((1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2))
			A[ni,elmts[ei,indir[2]] ] += common*(1+si*qno[qi])*(1+ti*qno[ri])/4
			A[ni,elmts[ei,indir[3]] ] += common*(1-(1+si*qno[qi])/2*(1-(1+ti*qno[ri])/2)- (1+si*qno[qi])*(1+ti*qno[ri])/4)
		end
	end
end

