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


close all
clear all


%% defining the planet (morphology)
%we need to define how the planet is shaped along with its atmosphere and
%distance from host star. 

core_to_surface_radius = 1000000; %[m]
atmospheric_thickness = 20000; %[m]
distance_from_host_star = 20000000; %[m]
eccentricity_of_planet = 0.01;
tilt_of_planet = 13; %[deg]
dV=1.0; %infintessimal volume segment, stepper in the for loop

%% defining the planet's atmosphere
%we then need to define the type of atmosphere that we want to model. For
%rayleigh scattering this is just the simple index of refraction and size
%of molecules which are doing the scattering. For Mie scattering regime, this 
%will be elaborate for physics and for coding. 

atmospheric_index_of_refraction = 1.1;
diameter_of_scattering_spheres=300*10^-11; %must be at least 10 times smaller

%% defining the star
%just as important as the planet we must define the star. We must tell the
%program the spectral itnesity of each star, or we can use the model at one
%particular wavelength, or a couple wavelengths spanning a certain range.
%This will allow us to obtain a certain 

intensity_nought=4*10^30; %lol this won't be the actual number.
lambda=400*10^-9; %this will be an upper and lower limit in future.  

%% important constants for this version of code
%so for this first version, we assume theta is always 0 i.e. light is only
%scattered perfectly forwards.
theta=0;

%% performing the calculations
%for this first version, let's assume that the planet is perfectly
%spehrical, so we must only compute the columns we want for one row. We can
%then just sum this up over the whole 2D circle that we would be seeing
%form earth. 

%step 1 using a cube of 1m compute the column length(the total number of
%times light is scattered). 
nn=1;
for i=core_to_surface_radius:dV:(core_to_surface_radius+atmospheric_thickness)
    %calculate the scatter column length
    total_radius=core_to_surface_radius+atmospheric_thickness;
    theta_one=acos(i/total_radius);
    column_length=sin(theta_one)*total_radius;
    mm=1;
    for j=1:dV:column_length
        scattered_intensity(nn).val(mm)=rayleigh_scatter(intensity_nought,theta,lambda,...
         atmospheric_index_of_refraction,dV,diameter_of_scattering_spheres);
        intensity_nought=scattered_intensity(nn).val(mm);
        mm=mm+1;
    end
    nn=nn+1;
    pause
end


















