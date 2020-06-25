function t_cross_sec=thomson_cross_section(particle_charge,particle_mass)
%the purpose of this function is to compute the thompson cross section for
%a given particle. The cross section is the only output from this function.
%   Inputs: particle_charge:double, charge of particle [C]
%           particle_mass:double, mass of particle [Kg]
%   Output:t_cross_sec:double, 

const=8*pi/3;
denom_const=4*pi*(8.854187*10.^-12);
speed_of_light=299792458;
t_cross_sec=const*(particle_charge.^2/(denom_const*particle_mass*(speed_of_light.^2))).^2;