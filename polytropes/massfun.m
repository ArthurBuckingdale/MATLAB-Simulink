function mass=massfun(n,dtdx,xi_1,density)
%the purpose of this function is to compute the mass of the polytrope given
%a central density and a polytropic index. 
G=6.67384.*(10.^(-11));
K=1;% set this to 1 for now
mass_const=(4.*3.14159)*((K./G).*((n+1)./(4.*3.14159))).^(3./2);
density=density.^((3-n)./(2*n));
mass=mass_const.*density.*(-((xi_1).^2).*(dtdx));
