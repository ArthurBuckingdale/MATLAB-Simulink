%The purpose of this script is to come up with a Rayleigh scattering model
%for light propagating through an atmosphere. The idea is we want to
%determine what the light will look like after passing through the
%atmosphere of a planet which is transiting its star. We will require a
%more complicated model based on Mie scattering, but a simpler start to
%identify a methodology is presented here. Once we begin with Mie
%scattering we factor in absorbtion as well and things start becoming very
%complicated. This will not be an accurate model necessarily, but once we
%switch out the physics for more complicated analysis, the accuracy will
%increase a lot. 

%please not we will be using m,s,kg for units here.


%% defining the planet
%we need to define how the planet is shaped along with its atmosphere and
%distance from host star.


core_to_surface_radius = 1000000; %[m]
atmospheric_thickness = 20000; %[m]
atmospheric_index_of_refraction = 1.1; 
distance_from_host_star = 20000000; %[m]
eccentricity_of_planet = 0.01;
tilt_of_planet = 13; %[deg]

%% defining the star
%just as important as the planet we must define the star. We must tell the
%program the spectral itnesity of each star, or we can use the model at one
%particular wavelength, or a couple wavelengths spanning a certain range.
%This will allow us to obtain a certain 
