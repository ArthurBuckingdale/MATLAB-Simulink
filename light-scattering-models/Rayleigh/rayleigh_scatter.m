function scattered_intensity=rayleigh_scatter(intensity_nought,theta,lambda,...
    refractive_index,distance_from_particle,sphere_diameter)
%the purpose of this function is to compute the rayleigh scatter of light
%as the result of interation with particles. 
%   Inputs: intensity_nought:double, initial light intensity [W/m^2]
%           theta:double, scattering angle [rad]
%           lambda:double, wavelength of light [m]
%           refractive_index:double, refractive index
%           distance_from_particle:double, distance away from particle [m]
%           sphere_diameter:double diameter of scattering sphere. 
%   Outputs:scattered_intensity:double,scattered light intensity [w/m^2]


%% firstly, we want to validate the particle/wavelength size
%remember that Rayleigh scattering is a special case of Mie scattering. We
%want to ensure that these parameters are met. 
if sphere_diameter > (lambda/10)
    disp('Rayleigh approximations fail here')
    disp('Particle diameter must be <1/10 of lambda')
    warning('Change input parameter, paritcle size not correct for Rayleigh approximation')
end

%% computing the intensity of scattered light.

refraction_term=(((refractive_index.^2)-1)/((refractive_index.^2)+2)).^2;
diameter_term=(sphere_diameter/2).^6;
lambda_term=(2*pi/lambda).^4;
spherical_term=(1+cos(theta).^2)/(2*distance_from_particle.^2);
scattered_intensity=intensity_nought.*spherical_term.*lambda_term.*diameter_term.*refraction_term;

