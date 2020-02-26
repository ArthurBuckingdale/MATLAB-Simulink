function radius=radfun(n,xi_one,density)
G=6.67384*(10.^(-11));
K=1;% set this to 1 for now
radius_const=sqrt((K/G)*((n+1)./(4*3.14159)));
central_density=density;
radius=radius_const*central_density.^(1-n/(2*n))*xi_one;
end