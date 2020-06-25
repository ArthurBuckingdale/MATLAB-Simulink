function lambda_prime=compton_scatter(lambda,theta)
%the purpose of this function is to compute the compton scattering angle
%for a given wavelength of light. It is a very basic model for how light is
%scattered, which makes it the perfect location to begin. We will model the
%new wavelength as a function of input wavelength and scatter angle. If the
%user desires the scattered wavelength as a function of angle with a
%constant input wavelength, or vice versa, they can use this function. We
%can also input both parameters to create a heat map for this.
%   Inputs: lambda:double, incident wavelength [m]
%           theta:double, scatter angle [deg]
%   Outputs:lambda_prime: double, sacttered wavelength of the light [m]

compton_wavelength=2.43*(10.^-12);
lambda_prime=compton_wavelength*(1-cosd(theta))+lambda;





