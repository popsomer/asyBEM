% Return the points and weights of a composite quadrature rule that has a hp-type grading towards the point a.
% Input
%   a		- The singularity on the collocation point
%   b		- The next collocation point
%   n		- Number of points, choose 7 for example
%   sigma	- The grading parameter, choose 0.15 for example
%   mu		- Parameter controlling degree p, choose 1 for example
%   minsize	- Minimum size of the interval, choose 1e-10 for example
% Output
%   x,w		- Nodes and weights
% References:
%   W. Gautschi, ``Orthogonal Polynomials: Computation and Approximation'', Clarendon Press, Oxford, 2004.
function [x,w] = quadHp(a, b, n, sigma, mu, minsize)
if nargin == 5
    minsize = 0;
end
p = max(2, floor(mu*n) + 1);
x = []; % We cannot easily preallocate because of the condition abs(u-v) < minsize.
w = [];
t2 = 1;
t1 = 1;
i = n;
while i >= 0
	if i == 0
		u = a;
		v = t1*(b-a)+a;
		p = 2;
	else
		t1 = t2 * sigma;
		u = t1*(b-a)+a;
		v = t2*(b-a)+a;
		if abs(u-v) < minsize
			% undo the last change and go to the last iteration
			t1 = t1/sigma;
			i = 0;
			continue
		end
	end
	% The following two lines are r_jacobi(N=p,0,0)
	nab = 2*(2:p-1)';
	ab = [zeros(p,1) [2; 1/3; nab.^2./(4*(nab+1).*(nab-1))]];
	
	% The following lines are gauss(N=p,ab)
	J = zeros(p);
	for n=1:p, J(n,n) = ab(n,1); end
	for n=2:p % The input-n is not used in the main loop, so we can re-use it
	        J(n,n-1) = sqrt(ab(n,2));
        	J(n-1,n) = J(n,n-1);
	end
	[V,D]=eig(J);
	[D,I]=sort(diag(D));
	xw=[D ab(1,2)*V(1,I)'.^2];
	
	% Scale [-1,1] to [u,v]
	x = [x; (v-u)*(xw(:,1)+1)/2 + u];
	w = [w; xw(:,2)*(v-u)/2];
	t2 = t1;
	i = i-1;
end
