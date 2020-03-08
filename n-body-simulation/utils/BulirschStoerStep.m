function [zF, info] = BulirschStoerStep(dynFun, tSpan, z0, tol)
% [zF, info] = BulirschStoerStep(dynFun, tSpan, z0, tol)
%
% Computes a single step using the Bulirsch-Stoer method
%
% INPUTS:
%   dynFun = function handle for the system dynamics
%       dz = dynFun(t,z)
%           t = scalar time
%           z = [nz,1] = state as column vector
%           dz = [nz,1] = derivative of state as column vector
%   tSpan = [1,2] = [t0, tF] = time span for the step
%   z0 = [nz,1] = initial state vector
%   tol = [nz,1] = error tolerance along each dimension. If tol is a
%       scalar, then all dimensions will satisfy that error tolerance.
%
% OUTPUTS:     (nt = n+1)
%   zF = [nz,1] = final state
%   info = struct with solver information:
%       .exit = exit condition
%           'converge' = successful convergence
%           'maxRefine' = reached max refinement; did not converge
%       .error = [nz,1] = error estimate along each dimension
%       .nFunEval = scalar int = count calls to dynFun
%       .nRefine = how many refinement steps were required?
%
% NOTES:
%   Implementation details:
%   http://web.mit.edu/ehliu/Public/Spring2006/18.304/implementation_bulirsch_stoer.pdf
%

% Set an upper limit on the number of mesh refinements in the sequence
nRefineMax = 8; 

% Simple logistics and memory allocation
n = 2*(1:nRefineMax);
nz = size(z0,1);
if length(tol)==1
    tol = tol*ones(size(z0));
end
T = zeros(nz,nRefineMax,nRefineMax);   %Extrapolation table
E = zeros(nz,nRefineMax);   %Error estimate table

info.exit = 'maxRefine';  %Assume that we fail to meet tolerance
for j=1:nRefineMax  %Loop over the sequence of improving meshes
    
    % Compute the estimate of the solution on the current mesh
    [~,z] = modifiedMidpointRule(dynFun, tSpan, z0, n(j));
    T(:,j,1) = z(:,end);
    
    if j>1
        
        % Compute the extrapolation table entries:
        for k=2:j
            num = T(:,j,k-1) - T(:,j-1,k-1);
            den = (n(j)/(n(j-k+1)))^2 - 1;
            T(:,j,k) = T(:,j,k-1) + num/den;
        end
        
        % Compute the error estimates:
        E(:,j) = abs(T(:,j,j-1) - T(:,j,j));
        
        % Check convergence:
        if all(E(:,j)<tol)
            info.exit = 'converged';
            break;
        end
    end
    
end

% Other useful things:
info.error = E(:,j);     %Error estimate
info.nFunEval = sum(n(1:j));    %number of function evaluations
info.nRefine = j;

% Return the estimate of the solution:
zF = T(:,j,j);

end
