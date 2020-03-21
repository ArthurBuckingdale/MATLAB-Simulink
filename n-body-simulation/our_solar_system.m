%the purpose of the script is a user friendly way to interact with the
%n-body simulator. It will give an example where we set up the solar
%system. This way, users of the n-body code will have a good framework to
%move from our solar system to others.

%% a cell array containing the names of the bodies

names{1}='sun';
names{2}='mercury';
names{3}='venus';
names{4}='earth';
names{5}='mars';
names{6}='jupiter';
names{7}='saturn';
names{8}='uranus';
names{9}='neptune';
%going to exclude pluto for the sim RIP

%% matrix containing the masses of these bodies
mass(1)=1.989e30;
mass(2)=3.285e23;
mass(3)=4.867e24;
mass(4)=5.9274e24;
mass(5)=6.417e23;
mass(6)=1.898e27;
mass(7)=5.683e26;
mass(8)=8.681e25;
mass(9)=1.024e26;

%% matrix containing the initial state vectors for the system
%so, we need to be careful here. We cannot simply just line up all of the
%planets on one axis and start the simulation. We need to initialise from a
%state that the planets currently occupy because of the different resonant
%orbits that we can have. There will be another document which illustrates
%how we can determine the initial state vectors for a system of planets.



%% choose the debugging boolean
%the debugging option will display all of the calculated outputs from
%various functions. If something is not making sense enable this option to
%see how the data is being passed into the function and used. Check that
%the values make sense. This script is partially a sanity check. We know
%our solar system well and can use it for a validation. 
DEBUG =0;

%% options
%there will be various options to configure the output of the function or
%tune it to the particular use case that we are having. (no options yet
%code is still being developped.

%% executing the function
%this section will call and execute the fucntion
[transit_timings,integration_error]=newtonian_gravity_nbody(body_name,body_mass,init_state_vector,DEBUG);

%% various post processing and such things
%show all the data that comes from this function...

























