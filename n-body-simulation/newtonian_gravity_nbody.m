function [transit_timings,integration_error]=newtonian_gravity_nbody(varargin)
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

%declaring some constants
newtonsg = 6.647e-11;

%let's now go ahead and solve a couple more cases and plot them
[t1,y1] = ode45(@(t,y) odefcn(t,y,n), tspan, ic);


%here is the odefun that is listed above.
    function dydt = odefcn(t,y,n)
        dydt = zeros(2,1);
        dydt(1) = y(2);
        dydt(2) = -y(1).^n-(2/t)*y(2);
    end

%first, we are going to specify a struct for the planet variables
%we want to do this to keep code legible and store all info related to a
%body in a database style format.

for i=1:length(input_variable)
    body_information(i).name=body_name(i);
    body_information(i).mass=body_mass(i);
    body_information(i).position=initial_state_vector(1:3,i);
    body_information(i).velocity=initial_state_vector(4:6,i);
end

%compute kinetic and potential energies for error monitoring
kinetic_energy=calc_kinetic_energy(body_information);
%create a function to compute relative distances between bodies 
relative_distance=calc_relative_distances(body_information);
%create a function to compute gravitational potential energy



    function scalar_gravitational_potential=calc_grav_potential_energy(body_information,relative_distance)
        newtonsg = 6.647e-11;
        for i=1:length(body_information)
            for j=1:length(body_information)
                scalar_gravitational_potential(i).value=newtonsg.*body_information
            end
        end
        
        
        
    end



end




