function [transit_timings,integration_error]=newtonian_gravity_nbody(body_name,body_mass,init_state_vector,DEBUG)
%The purpose of this function is to perform an n-body simulation for
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
%%%%%%%%%%%% %let's now go ahead and solve a couple more cases and plot them
%below is the required format that the ode's need to be in if you feel the
%need to make some changes to the code.
% [t1,y1] = ode45(@(t,y) odefcn(t,y,n), tspan, ic);
% %here is the odefun that is listed above.
%     function dydt = odefcn(t,y,n)
%         dydt = zeros(2,1);
%         dydt(1) = y(2);
%         dydt(2) = -y(1).^n-(2/t)*y(2);
%     end
%%%%%%%%%%%%

%this are temporary for development
close all

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
kinetic_energy_start=calc_kinetic_energy(body_information);
[relative_distance,grav_potential_start]=calc_relative_distances(body_information);
initial_system_energy=sum(unique(grav_potential_start))+sum([kinetic_energy_start(:).value]);
if DEBUG == 1
    disp('information about all bodies input to the system')
    unfold(body_information)
    disp('Kinetic energy of all bodies in the system')
    unfold(kinetic_energy_start)
    disp('Relative distance and gravitational potential energy matrix of the system')
    unfold(relative_distance)
    disp(grav_potential_start)
    fprintf('Initial enrgy contained in the system %d \n',initial_system_energy)
end

%% prep the differential equation values
%This section will contain the function that calculates the force which is
%exerted on each object. This calculation is quite awful because the force
%is a vector. Each component acting on each body needs to be calculated. We
%are going to borrow part of the architecture from the
%calc_relative_distances routine. We are not going to try and access the
%information in here because we would be comparing a load of strings to
%make sure we have the correct denominator for the force we are looking to
%calculate.

accel_of_gravity=calc_grav_accel(body_information);
if DEBUG == 1
    disp('acceleration due to gravity for the various components of the bodies')
    unfold(accel_of_gravity)
end

%% fill the format that the BS method requires as inputs
%this is going to create the structure that is the required input. We could
%go ahead and fill it with the previous function, but it' nice to have
%things that are discreet. And the structure to solve the ode's requires
%the velocity, which is the integral of the acceleration. In order to solve
%ODE's of degree two or higher, we require a system of equations. For
%example if we have an Nth degree ODE, we require N first order ode's to
%solve for this. This will get a bit messy in terms of code, simply because
%it is a relatively convoluted calculation, we need the effect of each
%component on every other body in the system.

%%%%%%%%%%%%%%%%
x=struct_to_vector(body_information);
%this is the equivalent of the odefun listed above.
dxdt=distance_derivatives_func(1,x,body_information);
%%%%%%%%%%%%%%%%

%% set the initial conditions and the time duration of integration.
%we need to set the duration of time that we wish to integrate the fucntion
%for. The initial conditions also need to be set,fortunately, we have
%created the struct_to_vector function which will do just that. We are
%going to test out the ode45 method first as it is a matlab suuported one.
%We will move on to the BS method after testing. Or it will be a callable
%option. We can observe the performance deltas here as well.

tspan = [0.1 10];
t=1;
ic=[struct_to_vector(body_information)]';
[t1,y1] = ode45(@(t,x) distance_derivatives_func(t,x,body_information), tspan, ic);

hold on 
plot(y1(:,7),y1(:,8),'r*')
plot(y1(:,1),y1(:,2),'b*')
plot(y1(:,13),y1(:,14),'y*')
hold off


%% Error from the ode solving
%we decided to use energy to monitor the error here. We must therfore grab
%our energy calculation function to determine if energy has been conserved.
%We must first append the state vectors of each body to the data format
%that our calculate energy function is expecting. The most recent state
%vectors for our system will be the last line of y1. remember that our
%calc_relative_distances function yields the gravitational potential.

body_information=update_body_position(y1(end,:),body_information);
unfold(body_information)
kinetic_energy_final=calc_kinetic_energy(body_information);
[~,grav_potential_final]=calc_relative_distances(body_information);
final_system_energy=sum(unique(grav_potential_final))+sum([kinetic_energy_final(:).value]);
energy_conservation=(final_system_energy./initial_system_energy);
fprintf('The final energy divided by the initial is %d \n',energy_conservation)
fprintf('The inital energy was %d, the final is %d \n',initial_system_energy,final_system_energy)




















































end




