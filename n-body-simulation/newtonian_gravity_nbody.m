function [transit_timings,integration_error]=newtonian_gravity_nbody(body_name,body_mass,init_state_vector,DEBUG)
%The purpose of this script is to perform an n-body simulation for
%planetary motion. It is going to be a re-write of my undergrad thesis. My
%programming skills have increased a lot since and I think I can make a
%better performing model in MATLAB. On top of this the data visualisation
%will be much better in MATLAB than in c++. I am going to work with the
%Bulirsh-Stoer integrator to solve this problem. It is a very accurate
%integrator and converges in a reasonable time frame. In the future this
%will be used by a RL agent to solve for the mass of bodies in an exoplanet
%system. It will also be able to play around with the number of bodies and
%the orbital parameters such as SMA.
%   inputs: body_name:cell array N length with N names of bodies
%           body_mass:vector, N length with masses of bodies in body_name
%           init_state_vector: 6xN matrix of (P,V) state vectors
%   outputs:transit_timings......
%           integration_error:double, contains the error from solving the
%           ode's.


% %let's now go ahead and solve a couple more cases and plot them
% [t1,y1] = ode45(@(t,y) odefcn(t,y,n), tspan, ic);
% %here is the odefun that is listed above.
%     function dydt = odefcn(t,y,n)
%         dydt = zeros(2,1);
%         dydt(1) = y(2);
%         dydt(2) = -y(1).^n-(2/t)*y(2);
%     end

%% initialise the system
%first, we are going to specify a struct for the planet variables
%we want to do this to keep code legible and store all info related to a
%body in a database style format. We are then going to compute all relative
%distances and the

for i=1:length(body_name)
    body_information(i).name=body_name{i};
    body_information(i).mass=body_mass(i);
    body_information(i).position=init_state_vector(1:3,i);
    body_information(i).velocity=init_state_vector(4:6,i);
end
kinetic_energy=calc_kinetic_energy(body_information);
[relative_distance,grav_potential]=calc_relative_distances(body_information);
initial_system_energy=sum(unique(grav_potential))+sum([kinetic_energy(:).value]);
if DEBUG == 1
    disp('information about all bodies input to the system')
    unfold(body_information)
    disp('Kinetic energy of all bodies in the system')
    unfold(kinetic_energy)
    disp('Relative distance and gravitational potential energy matrix of the system')
    unfold(relative_distance)
    disp(grav_potential)
    fprintf('Initial enrgy contained in the system %d \n',initial_system_energy)
end

%% prep the differential equation calculations
%This section will contain the function that calculates the force which is
%exerted on each object. This calculation is quite awful because the force
%is a vector. Each component acting on each body needs to be calculated. We
%are going to borrow part of the architecture from the
%calc_relative_distances routine. We are not going to try and access the
%information in here because we would be comparing a load of strings to
%make sure we have the correct denominator for the force we are looking to
%calculate.

accel_of_gravity=calc_grav_accel(body_information);
disp('acceleration due to gravity for the various components of the bodies')
unfold(accel_of_gravity)



   










































end




