function energy=compton_scatter_energy(lambda,theta)
%the purpose of this function is to compute the energy of the compton
%scattered particle.
%   Inputs: lambda:double, incident wavelength [m]
%           theta:double, scatter angle [deg]
%   Outputs:energy:double, energy of output particle [J]

planck_constant=6.626070*(10.^-34);
speed_of_light=299792458;
electron_rest_mass=9.109383*(10.^-31);
incident_energy=(planck_constant*speed_of_light)./lambda;
energy=incident_energy./(1+(incident_energy/electron_rest_mass*speed_of_light.^2)*(1-cosd(theta)));