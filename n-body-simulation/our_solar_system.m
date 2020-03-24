%the purpose of the script is a user friendly way to interact with the
%n-body simulator. It will give an example where we set up the solar
%system. This way, users of the n-body code will have a good framework to
%move from our solar system to others.This is going to serve as a test
%platform for now. 

%%%%%%%%%%%%%%%%IMPORTANT NOTE
%this is not the correct format of our solar system, im initialising the
%bodies with the correct orbital speed and SMA, but not the correct
%locations relative to each other. This causes the simulation to fail since
%our solar system has formed important resonances in its orbits when the
%bodies were born. 

%% a cell array containing the names of the bodies

names{1}='sun';
names{2}='mercury';
names{3}='venus';
names{4}='earth';
names{5}='mars';
names{6}='jupiter';
%names{7}='saturn';
%names{8}='uranus';
%names{9}='neptune';
%going to exclude pluto for the sim RIP

%% matrix containing the masses of these bodies
mass(1)=1.989e30;
mass(2)=3.285e23;
mass(3)=4.867e24;
mass(4)=5.9274e24;
mass(5)=6.417e23;
mass(6)=1.898e27;
%mass(7)=5.683e26;
%mass(8)=8.681e25;
%mass(9)=1.024e26;

%% matrix containing the initial state vectors for the system
%so, we need to be careful here. We cannot simply just line up all of the
%planets on one axis and start the simulation. We need to initialise from a
%state that the planets currently occupy because of the different resonant
%orbits that we can have. There will be another document which illustrates
%how we can determine the initial state vectors for a system of planets.
state_vectors(:,1)=[0.1;0.1;0.1;0.0000001;0.0000001;0.00];
state_vectors(:,4)=[149000000000.0;0;0;0;29770.0;0];
state_vectors(:,6)=[-778000000000.0;0;0;0;-13100.0;0];
state_vectors(:,2)=[0;69499000000.0;0;-48000.0;0;0];
state_vectors(:,3)=[0;-107480000000.0;0;35000.0;0;0];
state_vectors(:,5)=[-223780000000.0;0;0;0;-24000.0;0];

%% the orbital periods 
% in order to map the transit timimgs, we must subtract the number of
% periods which have elapsed, since we are observing a time which is
% increasing. we need to subtract the period in order to see the transit 
%shift irrelevant of how much time has passed. We will append these times
%to the body information struct which is in the newtonian_gravity_nbody
%routine below. Remember if there are N-bodies, the star will be at the
%center, so we will have only N-1 periods to deal with here(with some assumptions)
%stars to wiggle a bit obviously. The orbital periods here will need to be
%in seconds, as this is the unit in which time is specified for this
%routine. 
orbital_period(1)=7600521.6; %mercury
orbital_period(2)=19409760.0; %venus
orbital_period(3)=31450000.0; %earth 
orbital_period(4)=59354294.4; %mars
orbital_period(5)=374080032.0; %jupiter



%% choose the debugging boolean
%the debugging option will display all of the calculated outputs from
%various functions. If something is not making sense enable this option to
%see how the data is being passed into the function and used. Check that
%the values make sense. This script is partially a sanity check. We know
%our solar system well and can use it for a validation. 
DEBUG = 1;

%% options
%there will be various options to configure the output of the function or
%tune it to the particular use case that we are having. (no options yet
%code is still being developped.One is the time span. This will yield the
%duration which we integrate the system for. 
plot_xy = 0; %makes a overhead plot of the system
plot_trans_time = 0;

%% set the timespan for integration 
tspan = [0.001 92000000];

%% executing the function
%this section will call and execute the fucntion
[transit_timings,integration_error]=newtonian_gravity_nbody(names,mass,state_vectors,DEBUG,tspan,plot_xy,plot_trans_time,orbital_period);

%% various post processing and such things
%show all the data that comes from this function...

























