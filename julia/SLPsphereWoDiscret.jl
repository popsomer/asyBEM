# Compute windowed versions of the Single Layer potential applied to the known density function of the solution of an incident plane wave on a Dirichlet sphere.
using Cubature;
using GSL;

function field(k,the,r)
	res = exp(1im*k*cos(the)*r); # Incident field
	prevTerm = 20*k;
	term = 0;
	for n =0:round(1.4*k+40)
		n = Int64(n);
		term = (2*n+1)*1im^n*sf_bessel_jl(n,k)*(sf_bessel_jl(n,k*r)+1im*sf_bessel_yl(n,k*r))/(sf_bessel_jl(n,k)+1im*sf_bessel_yl(n,k))*sf_legendre_Plm(n,0,cos(the))
		if((abs(term) > abs(prevTerm) ) & (n > 10))
			return res;
		end
		res = res - term;
		prevTerm = term;
	end
	return res
end

function qDens(k,the)
	res = 1im*k*cos(the)*exp(1im*k*cos(the)); # Normal derivative of the incident field
	prevTerm = 20*k;
	term = 0;
	for n =0:round(1.4*k+40)
		n = Int64(n);
		if(n==0)
			derivHank = exp(1im*k)/k*(1+1im/k)
		else
			derivHank = (sf_bessel_jl(n-1,k)+1im*sf_bessel_yl(n-1,k) -(n+1)/k*(sf_bessel_jl(n,k)+1im*sf_bessel_yl(n,k)))
		end
		term = (2*n+1)*1im^n*sf_bessel_jl(n,k)*k*derivHank/(sf_bessel_jl(n,k)+1im*sf_bessel_yl(n,k))*sf_legendre_Plm(n,0,cos(the))
		if((abs(term) > abs(prevTerm) ) & (n > 10))
			return res;
		end
		res = res - term;
		prevTerm = term;
	end
	return res
end

function dis(x,y) 
	phis = x[1];
	thetas=x[2];
	rstar = x[3];

	phi = y[1];
	theta=y[2];
	sqrt((cos(theta)-rstar*cos(thetas))^2 + (sin(theta)*sin(phi)-rstar*sin(thetas)*sin(phis))^2 + (sin(theta)*cos(phi)-rstar*sin(thetas)*cos(phis))^2);
end

function inte(x,y,k)
	theta=y[2]; # Jacobian is sin(theta)
	exp(k*1im*dis(x,y) )/4/pi/dis(x,y)*qDens(k,theta)*sin(theta)
end

function windInte(x,y,k,a,b)
	phis = x[1];
	thetas=x[2]; 
	phi = y[1];
	theta=y[2];
	r = sqrt( (sin(theta)*sin(phi)-sin(thetas)*sin(phis))^2 + (sin(theta)*cos(phi)-sin(thetas)*cos(phis))^2)
	if ((r <= a) || (abs(r-1) < 1e-8))
		return	exp(k*1im*dis(x,y) )/4/pi/dis(x,y)*qDens(k,theta)*sin(theta);
	elseif r < b
		return 	exp(k*1im*dis(x,y) )/4/pi/dis(x,y)*qDens(k,theta)*sin(theta)*exp(2*exp((a-b)/(r-a))/( (r-a)/(b-a)-1));
	else
		return 0;
	end
end

function pCmplx(rea,ima)
	if ima > 0
		return string(rea, " + ", ima, "im ");
	else
		return string(rea, " ",   ima, "im ");
	end
end

function intePol(x,y,k)
	rPol = y[1]
	thetaPol = y[2];
	if((rPol > 2) || (rPol < 0) || (x[3]^2*sin(thetaPol)^2/x[1]^2 + cos(thetaPol)^2 + (sin(thetaPol)-x[2]*cos(thetaPol)/x[3])^2 < 0) )
		print("Error: r=", rPol, ", theta=", thetaPol, ", under sqrt= " , x[3]^2*sin(thetaPol)^2/x[1]^2 + cos(thetaPol)^2 + (sin(thetaPol)-x[2]*cos(thetaPol)/x[3])^2, "\n")
	end
	xCo = -(1/2)*rPol*sqrt(-rPol^2+4)*x[3]*sin(thetaPol)/(sqrt(x[3]^2*sin(thetaPol)^2/x[1]^2 + cos(thetaPol)^2 + (sin(thetaPol)-x[2]*cos(thetaPol)/x[3])^2)*x[1])+(1-(1/2)*rPol^2)*x[1];
	if(abs(xCo) > 1)
		print("Error: r=", rPol, ", theta=", thetaPol, ", xCo= " , xCo, ", x = ", x, "\n")
	end
	exp(k*1im*rPol)/4/pi*qDens(k,acos(xCo))
end

xt = [0.8, pi/4, 1.0];
xtCar = xt[3]*[cos(xt[2]), sin(xt[2])*sin(xt[1]), sin(xt[2])*cos(xt[1])];
xo = [0.8, pi/4, 1.4];
xi = [2*pi/3, -0.1, 0.5];
percDecay = 0.5;

ks = 2.^(3:5);
Ts = linspace(0.9,1.6,2);
tol = 1e-8;

print(field(ks[1],xo[2],xo[3] ), "= total field ext, BC outside = ", exp(1im*ks[1]*cos(xo[2])*xo[3]), "\n")


# Test to check whether the integral without and with window satisfies the boundary condition:
(valr,err) = hcubature(y -> begin real(inte(xt,y,kt)); end, [0,0],[2*pi,pi])
(vali,eri) = hcubature(y -> begin imag(inte(xt,y,kt)); end, [0,0],[2*pi,pi])
print(kt , " =k, integral at x=", xt, " is ", pCmplx(valr,vali), ", b(x) = ", exp(1im*kt*cos(xt[2])), "\n")

bt = 1.9;
at = 1.2;
(valr,err) = hcubature(y -> begin real(windInte(xt,y,kt,at,bt)); end, [0,0],[2*pi,pi])
(vali,eri) = hcubature(y -> begin imag(windInte(xt,y,kt,at,bt)); end, [0,0],[2*pi,pi])
print(kt , " =k, wind integral at x=", xt, " is ", pCmplx(valr,vali), ", b(x) = ", exp(1im*kt*cos(xt[2])), ", a=", at, ", b=", bt, "\n")
# Alternatives to hcubature: hcubature_v, hquadrature, hquadrature_v: the _v should execute in parallel

# Test the integral using polar coordinates
if false
	# This integral should be small as there is no critical point in the range
	(valr,err) = hcubature(y -> begin real(intePol(xtCar,y,ks[ki])); end, [0.67,1.77],[0.7,1.8])
	(valr,err) = hcubature(y -> begin imag(intePol(xtCar,y,ks[ki])); end, [0.67,1.77],[0.7,1.8])
	print(valr, " + (", vali, ") polar coords \n");
elseif false
	(valr,err) = hcubature(y -> begin real(intePol(xtCar,y,ks[ki])); end, [0,0],[2,2*pi])
	(vali,eri) = hcubature(y -> begin imag(intePol(xtCar,y,ks[ki])); end, [0,0],[2,2*pi])
	print(valr, " + (", vali, ") polar coords \n");
end

if false
	# Integrate along some path of steepest descent
	function integr(p) exp(50*im*(-0.72764585274627917e-1*p^2+(.59216178402453730*1im)*p+(.40824829046386303*(-0.72764585274627917e-1*p^2+(.59216178402453730*1im)*p))*sqrt(-(-0.72764585274627917e-1*p^2+(.59216178402453730*1im)*p)^2+4)+(1/3*(1-(1/2)*(-0.72764585274627917e-1*p^2+(.59216178402453730*1im)*p)^2))*sqrt(3)))*(0.45572809000084120e-1-0.72764585274627917e-1*p^2+(.59216178402453730*1im)*p+(-0.72764585274627917e-1*p^2+(.59216178402453730*1im)*p)^2)*(-.14552917054925583*p+.59216178402453730*1im) end # thi=1,sdi=1
	maxP = (17/50)*log(10);
	(valr,err) = hcubature(x -> begin real(integr(x[1])); end, 0, maxP)
	(vali,eri) = hcubature(x -> begin imag(integr(x[1])); end, 0, maxP)
	print(valr+vali*1im, " =val, err = ", err+eri*1im, " , maxP =",maxP ," , integr at maxP = ", integr(maxP), '\n')
end



intVals = (0+0im)*zeros(length(ks),2)
intWindVals = (0+0im)*zeros(length(ks),length(Ts));
extVals = (0+0im)*zeros(length(ks),2)
extWindVals = (0+0im)*zeros(length(ks),length(Ts));
realVals = (0+0im)*zeros(length(ks),2)
windVals = (0+0im)*zeros(length(ks),length(Ts));
realErrs = (0+0im)*zeros(length(ks),1)
windErrs = (0+0im)*zeros(length(ks),length(Ts));


for ki = 1:length(ks)
	print("\n k = ", ks[ki], " ", now());
	valr = 0.0; err = 0.0; vali = 0.0; eri = 0.0;
	(valr,err) = hcubature(y -> begin real(inte(xt,y,ks[ki])); end, [0,0],[2*pi,pi], reltol = tol)
	(vali,eri) = hcubature(y -> begin imag(inte(xt,y,ks[ki])); end, [0,0],[2*pi,pi], reltol = tol)
	realVals[ki,2] = exp(1im*ks[ki]*cos(xt[2]));
	realVals[ki,1] = valr+1im*vali - realVals[ki,2];
	realErrs[ki,1] = err+1im*eri;

	(valr,err) = hcubature(y -> begin real(inte(xo,y,ks[ki])); end, [0,0],[2*pi,pi], reltol = tol)
	(vali,eri) = hcubature(y -> begin imag(inte(xo,y,ks[ki])); end, [0,0],[2*pi,pi], reltol = tol)
	extVals[ki,2] = field(ks[ki],xo[2], xo[3]);
	extVals[ki,1] = valr+1im*vali + extVals[ki,2]; # Maybe also try valr+1im*vali + exp(1im*ks[ki]*cos(xo[2])*xo[3]);
	print("Did BC and exterior, starting interior at ", now());
	
	(valr,err) = hcubature(y -> begin real(inte(xi,y,ks[ki])); end, [0,0],[2*pi,pi], reltol = tol)
	(vali,eri) = hcubature(y -> begin imag(inte(xi,y,ks[ki])); end, [0,0],[2*pi,pi], reltol = tol)
	intVals[ki,2] = exp(1im*ks[ki]*cos(xi[2])*xi[3]);
	intVals[ki,1] = valr+1im*vali - intVals[ki,2];

	for ti = 1:length(Ts)
		print("Starting windowed version with T = ", Ts[ti], " at ", now());
		valr = 0.0; err = 0.0; vali = 0.0; eri = 0.0;
		(valr,err) = hcubature(y -> begin real(windInte(xt,y,ks[ki],(1-percDecay)*Ts[ti],Ts[ti])); end, [0,0],[2*pi,pi], reltol = tol)
		(vali,eri) = hcubature(y -> begin imag(windInte(xt,y,ks[ki],(1-percDecay)*Ts[ti],Ts[ti])); end, [0,0],[2*pi,pi], reltol = tol)
		windVals[ki,ti] = valr+1im*vali - realVals[ki,2];
		windErrs[ki,ti] = err+1im*eri;

		(valr,err) = hcubature(y -> begin real(windInte(xo,y,ks[ki],(1-percDecay)*Ts[ti],Ts[ti])); end, [0,0],[2*pi,pi], reltol = tol)
		(vali,eri) = hcubature(y -> begin imag(windInte(xo,y,ks[ki],(1-percDecay)*Ts[ti],Ts[ti])); end, [0,0],[2*pi,pi], reltol = tol)
		extWindVals[ki,ti] = valr+1im*vali + extVals[ki,2];
		print(", T ext/int, ", now());

		(valr,err) = hcubature(y -> begin real(windInte(xi,y,ks[ki],(1-percDecay)*Ts[ti],Ts[ti])); end, [0,0],[2*pi,pi], reltol = tol)
		(vali,eri) = hcubature(y -> begin imag(windInte(xi,y,ks[ki],(1-percDecay)*Ts[ti],Ts[ti])); end, [0,0],[2*pi,pi], reltol = tol)
		intWindVals[ki,ti] = valr+1im*vali - intVals[ki,2];
	end
	f = open("simRes.txt","w");
	write(f,"xboundary=", string(xt), "xouter=", string(xo), "xinterior=", string(xi), "\n");
	write(f,"percDecay=", string(percDecay), ", tol=", tol, "\n")
	write(f,"ks=", string(ks), "\n Ts=", string(Ts), "\n")
	write(f,"intVals [scat-inc, inc] ", string(intVals), "\n")
	write(f,"intWindVals [scat-inc] ", string(intWindVals), "\n")
	write(f,"extVals [scat+exact field, exact field] ", string(extVals), "\n")
	write(f,"extWindVals [scat+exact field] ", string(extWindVals), "\n")
	write(f,"realVals [scat-inc, inc] ", string(realVals), "\n")
	write(f,"windVals (only [scattered - incident]) ", string(windVals), "\n")
	write(f,"realErrs", string(realErrs), "\n")
	write(f,"windErrs", string(windErrs), "\n")
	close(f);

	writedlm("intVals", intVals)
	writedlm("intWindVals", intWindVals)
	writedlm("extVals", extVals)
	writedlm("extWindVals", extWindVals)
	writedlm("realVals", realVals)
	writedlm("realErrs", realErrs)
	writedlm("windVals", windVals)
	writedlm("windErrs", windErrs)
end
print("\n\n")
run(`cat simRes.txt`)
print("\n\n")

