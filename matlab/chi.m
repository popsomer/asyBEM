% The C-infinity window function.
% Input
%   t                     - Abcis
%   lambda < l < r < rho  - chiw is 0 for t not in (lambda, rho) and 1 for t in [l, r]
%   firstCall		  - nonzero if windows might have to be added periodically
% Output
%   chiw                  - Window value
function chiw = chi(t, lambda, l, r, rho, firstCall)
if firstCall
    chiw = chi(t, lambda, l, r, rho, 0) + chi(t-1, lambda, l, r, rho, 0) + chi(t+1, lambda, l, r, rho, 0);
    return
end
chiw = zeros(size(t));
for ti = 1:length(t)
    if (t(ti) <= lambda) || (t(ti) >= rho)
        chiw(ti) = 0;
    elseif (t(ti) >= l) && (t(ti) <= r)
        chiw(ti) = 1;
    else
        if t(ti) < l
            u = (t(ti) - l)/(lambda - l);
        else
            u = (t(ti) - r)/(rho - r);
        end
        if (u < 0) || (u > 1)
            error('t invalidy parsed');
        end
        % We don't need a series expansion for u near 0 nor 1, nor equal to 0 but u could be 1 and still t(ti) > r.
        if u == 1
            chiw(ti) = 0;
        else
            chiw(ti) = exp( 2*exp(-1/u)/(u-1) );
            if chiw(ti) < eps % Try to increase compression, but this could go wrong when multiplying by large matrix entries.
                chiw(ti) = 0;
            end
        end
    end
end
