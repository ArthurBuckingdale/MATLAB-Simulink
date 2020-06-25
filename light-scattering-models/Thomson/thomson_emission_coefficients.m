function [radial_em,tangential_em]=thomson_emission_coefficients(energy,time,...
    area,wavelength,particle_density,particle_mass,particle_charge,chi)
%the purpose of this function is to compute the emission coefficients as a
%result of Thomson scattering. We will obtain a radial and tangential
%component here. 
%   Inputs: energy:double, energy of the incident light [J]
%           time:double, time over which this energy passes [s]
%           area:double, area over which this energy passes [m^2]
%           wavelength:double, wavelengths of light in this energy. [m]
%           particle_density:double, density of particles at scattering
%               location
%           particle_mass:double, mass of particle [Kg]
%           particle_charge:double, charge of particle [C]
%           chi:double, angle between incidence and observer [deg]
%   Outputs:radial_em:double, radial emission coefficient
%           tangential_em:double, tangential emission coefficient

incident_flux=energy/(time*area*wavelength);
tangential_em=(pi*thomson_cross_section(particle_charge,particle_mass)*incident_flux ...
    *particle_density)/2;
radial_em=tangential_em*(cosd(chi).^2);
