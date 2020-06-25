function theta=inverse_compton_scatter(lambda,lambda_prime)
%the purpose of this function is to take as inputs the inital and final
%wavelengths to determine what the scatter angle is. 
%   Inputs: lambda:double, incident wavelength [m]
%           lambda_prime:double, output wavelength [m]
%   Outputs:theta:double, angle between light rays [deg]

compton_wavelength=2.43*(10.^-12);
theta=acosd(1-((1/compton_wavelength)*(lambda_prime-lambda)));

